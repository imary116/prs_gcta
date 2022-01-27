import hailtop.batch as hb
import hail as hl

# Step 1: Perform clumping so that we can approximately capture the right level of causal signal in which only weakly correlated SNPs are retained but preferentially the ones that are most associated with the phenotype under study
def clumping(b, input_files, sumstats, label):
    j = b.new_job(name=f'{label}-clumping')  # define job and label it
    j.image('hailgenetics/genetics:0.2.67')  # use a Docker image that contains PLINK
    j.memory('10Gi')  # set memory - value obtain from the Hail batch cookbook
    j.storage('40Gi')  # increase storage
    select_col = "'NR!=1{print $3}'"  # remove header and which column to keep - index snp id

    # run PLINK clumping command
    j.command(f'''
plink --bfile {input_files} \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump {sumstats} \
    --clump-snp-field SNP \
    --clump-field P 

# From the .clumped file produced by PLINK, select the 3rd column (index SNP ID - snp is in the 3rd column in the clumping output) and save that to a temporary Batch file so it can be used in the PRS calculation 
# save the log file too 
awk {select_col} plink.clumped > {j.clumped}
mv plink.log {j.log} 
    ''')
    return j

# Step 2: Calculate PRS with different p-value thresholds
def calculate_prs(b, input_files, sumstats, valid_snp, label):
    j = b.new_job(name=f'{label}-calculating PRS')  # define job and label it
    j.image('hailgenetics/genetics:0.2.67')  # use a Docker image that contains PLINK
    j.memory('10Gi')  # set memory
    j.storage('40Gi')  # increase storage
    select_col = "'{print $2,$7}'"  # which columns to keep - snp and p

    # specify desired outputs
    j.declare_resource_group(ofile={
        '0.5.profile': '{root}.0.5.profile',
        '0.4.profile': '{root}.0.4.profile',
        '0.3.profile': '{root}.0.3.profile',
        '0.2.profile': '{root}.0.2.profile',
        '0.1.profile': '{root}.0.1.profile',
        '0.05.profile': '{root}.0.05.profile',
        '0.001.profile': '{root}.0.001.profile',
        'log': '{root}.log'
    })

    # run commands
    j.command(f'''

    # obtain SNP IDs and corresponding p-values from the transformed sumstats file 
    awk {select_col} {sumstats} > snp_pvalue

    # create a file containing different p-value thresholds for inclusion of SNPs in the PRS calculation 
    # the threshold boundaries are inclusive
    echo "0.001 0 0.001" > range_list 
	echo "0.05 0 0.05" >> range_list 
	echo "0.1 0 0.1" >> range_list 
	echo "0.2 0 0.2" >> range_list 
	echo "0.3 0 0.3" >> range_list 
	echo "0.4 0 0.4" >> range_list 
	echo "0.5 0 0.5" >> range_list

    # run PLINK PRS command - --score values: 2nd column is the SNP ID, 4th column is the effective allele information, the 6th column is the effect size estimate
    plink --bfile {input_files} \
        --score {sumstats} 2 4 6 header \
        --q-score-range range_list snp_pvalue \
        --extract {valid_snp} \
        --out {j.ofile}
    ''')
    j.command(f'echo {j.ofile.log}')
    return j


if __name__ == '__main__':
    backend = hb.ServiceBackend(billing_project='daly-neale-sczmeta', bucket='imary116')  # set up backend

    b = hb.Batch(backend=backend, name='PRS-PLINK-multi-ds(BIP-updated)')  # define batch

    ### import needed datasets ###

    # read in QCed sum stats file
    sumstats = b.read_input('gs://imary116/prs/PLINK/input/sumstats/bip_updated.QC.Transformed') # bip
    #sumstats = b.read_input('gs://imary116/prs/PLINK/input/sumstats/mdd_updated.QC.Transformed') # mdd

    # get the GCS paths of the plink files from a file and pick out the identifying name for each file
    with hl.hadoop_open('gs://imary116/prs/PLINK/input/file_paths.txt') as file_paths:
        names = []
        for path in file_paths:
            n = path.split('.')[-3]
            names.append(n)

    # recreate the paths using the identifying names so that correct QC files are imported for each file set
    for i in set(names):
        default_path = f'gs://imary116/prs/PLINK/input/{i}/scz.'
        # read in GWAS data set that's in PLINK binary PED format
        input_files = b.read_input_group(
            bed=default_path + i + ".QC.bed",
            bim=default_path + i + ".QC.bim",
            fam=default_path + i + ".QC.fam")


        # Step 1: Clumping
        run_clumping = clumping(b, input_files, sumstats, i)
        b.write_output(run_clumping.log, f'gs://imary116/prs/PLINK/output/log_files/bip_clump_{i}')  # write out log file

        # Step 2: PRS calculation
        run_prs = calculate_prs(b, input_files, sumstats, run_clumping.clumped, i)
        b.write_output(run_prs.ofile.log, f'gs://imary116/prs/PLINK/output/log_files/bip_prs_{i}')  # write out log file

        ### write out final results ###
        b.write_output(run_prs.ofile, f'gs://imary116/prs/PLINK/output/bip_updated/{i}/{i}')

    b.run(open=True, wait=False)  # run batch

    backend.close()








