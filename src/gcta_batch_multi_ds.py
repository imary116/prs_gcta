import hailtop.batch as hb
import hail as hl

# Step 1: Create the phenotype file - from the .fam file - needed to run greml
def create_phen(b, fam_file, label):
    j = b.new_job(name=f'{label}-create_phenotype_file')  # define job and label it
    j.storage('5Gi')  # increase storage
    columns = "'{print $1, $2, $6}'" # which columns to keep - family ID, individual ID, and phenotypes (in the data we have rn cases = 2, controls = 1, missing = "-9" or "NA")
    j.command(f''' awk {columns} {fam_file}  > {j.ofile}''') # subset .fam file to only include the specified columns and no header line
    return j

# Step 2: GCTA-GRM - calculating the genetic relationship matrix (GRM) from all the autosomal SNPs (assuming the maf is 0.01 by default)
def grm(b, input_files, label, maf_val: float = 0.01, thread_val: int = 10):
    j = b.new_job(name=f'{label}-GCTA-GRM')  # define job and label it
    # use a Docker image that contains gcta64
    j.image('hailgenetics/genetics:0.2.67')
    j.cpu(4)  # set cpu
    j.storage('10Gi')  # increase storage
    # the genetic relationship matrix will be saved in the files bg_files.grm.bin, bg_files.grm.N.bin and bg_files.grm.id
    j.declare_resource_group(ofile={
        'grm.bin': '{root}.grm.bin',
        'grm.N.bin': '{root}.grm.N.bin',
        'grm.id': '{root}.grm.id',
        'log': '{root}.log'
    })
    # run GCTA-GRM analysis
    j.command(
        f'''/gcta/gcta_1.93.1beta/gcta64 --bfile {input_files} --autosome --maf {maf_val} --make-grm --out {j.ofile} --thread-num {thread_val}''')
    j.command(f'echo {j.ofile.log}')
    return j

# Step 3: GCTA-GREML - estimate the variance explained by all the autosomal SNPs on the observed 0-1 scale and transform the estimate to that on the underlying liability scale (assuming the disease prevalence is 0.01 by default)
def greml(b, input_files, phenotype_file, label, preval: float = 0.01, thread_val: int = 10):
    j = b.new_job(name=f'{label}-GCTA-GREML')  # define job and label it
    # use a Docker image that contains gcta64
    j.image('hailgenetics/genetics:0.2.67')
    j.cpu(4)  # set cpu
    j.storage('10G')  # increase storage
    # run GCTA-GREM analysis
    j.declare_resource_group(ofile={
        'hsq': '{root}.hsq',
        'log': '{root}.log'
    })
    j.command(
        f'''/gcta/gcta_1.93.1beta/gcta64 --grm {input_files} --pheno {phenotype_file} --reml --prevalence {preval} --out {j.ofile} --thread-num {thread_val}''')
    j.command(f'echo {j.ofile.log}')
    return j

# Step 4: Parse the greml output into a table - basically remove the two rows with sentences
def parse_output(b, file, label):
    j = b.new_job(name=f'{label}-parse_output')  # define job and label it
    # Docker image containing python3
    j.image('hailgenetics/genetics:0.2.67')
    # parse table and remove unwanted rows using pandas
    j.command(f'''
python3 -c "
import pandas as pd

data = pd.read_table('{file}.hsq')
data = data[data.Variance.notnull()]
data.to_csv('{j.ofile}', index=False, sep='\t')"

    ''')
    return j


if __name__ == '__main__':
    backend = hb.ServiceBackend(billing_project='daly-neale-sczmeta', bucket='imary116')  # set up backend

    b = hb.Batch(backend=backend, name='GCTA-multi-ds')  # define batch

    # MAF, thread num, prevalence as options in functions - might be able to set them up as arguments later that can be indicated from the command line
    maf_val = 0.01
    thread_val = 10
    preval = 0.01

    ### import needed data files ###
    # get the GCS paths of the files from a file and pick out the identifying name for each file
    with hl.hadoop_open('gs://imary116/prs/PLINK/input/file_paths.txt') as file_paths:
        names = []
        for path in file_paths:
            n = path.split('.')[-3]
            names.append(n)

    # recreate the paths using the identifying names so that correct files are imported for each file set
    for i in set(names):
        default_path = f'gs://imary116/prs/PLINK/input/{i}/scz.'
        # read in GWAS data set that's in PLINK binary PED format
        input_files = b.read_input_group(
            bed=default_path + i + ".QC.bed",
            bim=default_path + i + ".QC.bim",
            fam=default_path + i + ".QC.fam")

        ### run the following functions for each file set###

        # Step 1: .fam -> phenotype file
        run_create = create_phen(b, input_files.fam, i)
        #b.write_output(run_create.ofile, 'gs://imary116/gcta/one_ds/input_files/bg.phen')

        # Step 2: GCTA-GRM
        run_grm = grm(b, input_files, i, maf_val, thread_val)
        b.write_output(run_grm.ofile.log, f'gs://imary116/gcta/multi_ds/output/log_files/grm_{i}')  # write out log file

        # Step 3: GCTA-GREML
        run_greml = greml(b, run_grm.ofile, run_create.ofile, i, preval, thread_val)
        b.write_output(run_greml.ofile.log, f'gs://imary116/gcta/multi_ds/output/log_files/greml_{i}')  # write out log file

        # Step 4: output -> desired table format
        run_parse = parse_output(b, run_greml.ofile, i)

        ### write out final results ###
        b.write_output(run_parse.ofile, f'gs://imary116/gcta/multi_ds/output/{i}.txt')

    b.run(open=True, wait=False)  # run batch

    backend.close()








