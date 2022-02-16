import hailtop.batch as hb
import hail as hl


# Step 1: merge files
def merge(b, file1, file2_bed, file2_bim, file2_fam):
    j = b.new_job(name="Merge_files")  # define job and label it
    j.cpu(8)  # set cpu
    j.memory('highmem')  # increase memory
    j.storage('50Gi')  # increase storage
    # use a Docker image that contains plink
    j.image('hailgenetics/genetics:0.2.67')
    # run has multiple outputs
    j.declare_resource_group(ofile={
        'bed': '{root}.bed',
        'bim': '{root}.bim',
        'fam': '{root}.fam',
        'log': '{root}.log'
    })
    j.command(f'''
plink --bfile {file1} --bmerge {file2_bed} {file2_bim} {file2_fam} --make-bed --out {j.ofile}''')
    j.command(f'echo {j.ofile.log}')
    return j


# Step 2: run standard GWAS QC on the merged dataset
def GWAS_qc(b, merged_file):
    j = b.new_job(name="Run_GWAS_QC")  # define job and label it
    j.cpu(4)  # set cpu
    j.memory('highmem')  # increase memory
    j.storage('30Gi')  # increase storage
    # use a Docker image that contains plink
    j.image('hailgenetics/genetics:0.2.67')
    # run has multiple outputs
    j.declare_resource_group(ofile={
        'bed': '{root}.bed',
        'bim': '{root}.bim',
        'fam': '{root}.fam',
        'log': '{root}.log'
    })
    j.command(f''' 
plink --bfile {merged_file} --maf 0.01 --hwe 1e-6 --geno 0.01 --make-bed --out {j.ofile}
''')
    j.command(f'echo {j.ofile.log}')
    return j


# Step 3: Create pheno file - cases = 2 and controls = 1
def create_pheno(b, fam1, fam2, fam_merged):
    j = b.new_job(name="Create_phenotype_file")  # define job and label it
    j.cpu(8)  # set cpu - don't need more than 1 cpu
    j.memory('highmem')  # increase memory - how much memory each cpu has
    j.storage('50Gi')  # increase storage - how big are your files?
    # use a Docker image that contains python
    j.image('hailgenetics/genetics:0.2.67')
    j.command(f''' 
python3 -c "
import pandas as pd
import numpy as np

# read-in fam files
file1 = pd.read_csv('{fam1}', sep=' ', header=None) 
file2 = pd.read_csv('{fam2}', sep=' ', header=None)
mfile = pd.read_csv('{fam_merged}', sep=' ', header=None)

# 4 columns: FID, IID, and phenotypes from two data separated into two columns (1, 2 & missing value "-9") 
mfile['trait1'] = np.where((mfile[1].isin(file1[1])), mfile[5].astype('Int64'), -9) 
mfile['trait2'] = np.where((mfile[1].isin(file2[1])), mfile[5].astype('Int64'), -9)
mfile_pheno = mfile[[0, 1, 'trait1', 'trait2']]

# write it out 
mfile_pheno.to_csv(r'{j.ofile}', sep=' ', index=False, header=False)" 
''')
    return j


# Step 4: GCTA-GRM - calculating the genetic relationship matrix (GRM) from all the autosomal SNPs (assuming the maf is 0.01 by default)
def grm(b, input, snp_list1, snp_list2, maf_val: float = 0.01, thread_val: int = 10):
    j = b.new_job(name="GCTA-GRM")  # define job and label it
    # use a Docker image that contains gcta64
    j.image('hailgenetics/genetics:0.2.67')
    j.cpu(4)  # set cpu - j.cpu({thread value}) - 4
    j.storage('20Gi')  # increase storage
    # the genetic relationship matrix will be saved in the files input.grm.bin, input.grm.N.bin and input.grm.id
    j.declare_resource_group(ofile={
        'grm.bin': '{root}.grm.bin',
        'grm.N.bin': '{root}.grm.N.bin',
        'grm.id': '{root}.grm.id',
        'log': '{root}.log'
    })
    # run GCTA-GRM analysis - use subset of autosomal SNPs with MAF < 0.01
    j.command(f'''

# merge snp lists first (obtained from QC of individual data - prune.in)
cat {snp_list1} {snp_list2} > snps_to_keep 

/gcta/gcta_1.93.1beta/gcta64 --bfile {input} --extract snps_to_keep --autosome --maf {maf_val} --make-grm --out {j.ofile} --thread-num {thread_val}

''')
    j.command(f'echo {j.ofile.log}')
    return j


# Step 5: GCTA-GREML - bivariate GREML analysis to estimate the genetic correlation between two sample sets with the same binary disease trait from case control studies
# assuming the disease prevalence is 0.01 by default
def greml(b, input, phenotype, preval1: float = 0.01, preval2: float = 0.01, thread_val: int = 10):
    j = b.new_job(name="GCTA-GREML")  # define job and label it
    # use a Docker image that contains gcta64
    j.image('hailgenetics/genetics:0.2.67')
    j.cpu(4)  # set cpu - j.cpu({thread value}) - 4
    j.storage('20Gi')  # increase storage
    # run GCTA-GREM analysis
    j.declare_resource_group(ofile={
        'hsq': '{root}.hsq',
        'log': '{root}.log'
    })
    # for bivariant use "--reml-bivar" and "--reml-bivar-prevalence" instead of "--reml" and "--prevalence" - --reml-maxit
    j.command(
        f'''/gcta/gcta_1.93.1beta/gcta64 --grm {input} --pheno {phenotype} --reml-bivar --reml-bivar-prevalence {preval1} {preval2} --out {j.ofile} --thread-num {thread_val}''')
    j.command(f'echo {j.ofile.log}')
    return j


# Step 6: Parse greml output into a table (remove unwanted lines)
def parse_output(b, file):
    j = b.new_job(name="Parse_output")  # define job and label it
    # Docker image containing python3
    j.image('hailgenetics/genetics:0.2.67')
    # j.cpu(4)  # set cpu
    # j.storage('10G')  # increase storage
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

    # MAF, thread num, prevalence as options in functions - might be able to set them up as arguments later that can be indicated from the command line?
    maf_val = 0.01
    thread_val = 10
    preval1 = 0.01
    preval2 = 0.01

    # pair up samples for the bivariate run
    sample_names = ['celso', 'cgs1c', 'gawli', 'xclm2', 'xclo3', 'xswe5', 'xswe6']
    for first_sample in sample_names:
        for second_sample in sample_names:
            if second_sample > first_sample:
                b = hb.Batch(backend=backend, name=f'GCTA-bivariate: {first_sample}_{second_sample}')  # define batch

                # read individual fam files for phenotype file
                file1 = b.read_input_group(
                    bed=f'gs://imary116/prs/PLINK/input/{first_sample}/scz.{first_sample}.QC.bed',
                    bim=f'gs://imary116/prs/PLINK/input/{first_sample}/scz.{first_sample}.QC.bim',
                    fam=f'gs://imary116/prs/PLINK/input/{first_sample}/scz.{first_sample}.QC.fam')

                file2 = b.read_input_group(
                    bed=f'gs://imary116/prs/PLINK/input/{second_sample}/scz.{second_sample}.QC.bed',
                    bim=f'gs://imary116/prs/PLINK/input/{second_sample}/scz.{second_sample}.QC.bim',
                    fam=f'gs://imary116/prs/PLINK/input/{second_sample}/scz.{second_sample}.QC.fam')

                # read in files with snps to keep
                snp_list1 = b.read_input(f'gs://imary116/prs/PLINK/input/{first_sample}/scz.{first_sample}.QC.prune.in')
                snp_list2 = b.read_input(
                    f'gs://imary116/prs/PLINK/input/{second_sample}/scz.{second_sample}.QC.prune.in')

                # file1_size = hl.utils.hadoop_stat(file1.bed)['size'] # format is number + the letter 'G'

                # Step 1: merging
                run_merge = merge(b, file1, file2.bed, file2.bim,
                                  file2.fam)  # outputs from this run are merged bed, bim and fam - run_merge.ofile
                b.write_output(run_merge.ofile.log,
                               f'gs://imary116/gcta/bivariate_run/output/log_files/merge_{first_sample}_{second_sample}')  # write out log file

                # Step 2: standard GWAS QC
                run_GWAS_qc = GWAS_qc(b,
                                      run_merge.ofile)  # outputs from this run are QCed [merged] bed, bim and fam - run_GWAS_qc.ofile
                b.write_output(run_GWAS_qc.ofile.log,
                               f'gs://imary116/gcta/bivariate_run/output/log_files/GWASqc_{first_sample}_{second_sample}')  # write out log file

                # Step 3: phenotype file
                run_pheno = create_pheno(b, file1.fam, file2.fam,
                                         run_GWAS_qc.ofile.fam)  # output is a phenotype file - run_pheno.ofile

                ## save intermediate files?
                # b.write_output(run_GWAS_qc.ofile, 'gs://imary116/gcta/multi_ds/output_files/bg.phen')
                # b.write_output(run_create.ofile, 'gs://imary116/gcta/one_ds/input_files/bg.phen')

                # Step 4: GCTA-GRM
                run_grm = grm(b, run_GWAS_qc.ofile, snp_list1, snp_list2, maf_val, thread_val)
                b.write_output(run_grm.ofile.log,
                               f'gs://imary116/gcta/bivariate_run/output/log_files/grm_{first_sample}_{second_sample}')  # write out log file

                # Step 5: GCTA-GREML
                run_greml = greml(b, run_grm.ofile, run_pheno.ofile, preval1, preval2, thread_val)
                b.write_output(run_greml.ofile.log,
                               f'gs://imary116/gcta/bivariate_run/output/log_files/greml_{first_sample}_{second_sample}')  # write out log file

                # Step 6: Parse output
                run_parse = parse_output(b, run_greml.ofile)

                b.write_output(run_parse.ofile,
                               f'gs://imary116/gcta/bivariate_run/output/{first_sample}_{second_sample}.txt')  # write out final output

                b.run(open=False, wait=False)  # run batch and don't open the browser

    backend.close()