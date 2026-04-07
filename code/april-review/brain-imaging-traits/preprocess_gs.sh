###########################################################################################
######                          Download summary statistics                          ######
###########################################################################################
# download the brain related traits:
cd /u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0124.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0125.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0126.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0127.txt.gz

curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0013.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0014.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0015.txt.gz
curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0016.txt.gz

###########################################################################################
######                                    Create step1                               ######
###########################################################################################
# specify common file paths:
magma_dir="/u/home/l/lixinzhe/project-geschwind/software/magma"

# specify input file:
${magma_dir}/magma \
    --annotate window=10,10 \
    --snp-loc ${magma_dir}/aux/g1000_eur.bim \
    --gene-loc ${magma_dir}/aux/NCBI37.3.gene.loc \
    --out ${magma_dir}/aux/step1

out_magma_dir="/u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas/magma_output/"

###########################################################################################
######           preprocess the summary statistics and creates magma file        ######
###########################################################################################
for item in /u/home/l/lixinzhe/project-cluo/data/brain_volume_gwas/*.txt.gz; do
    qsub /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/april-review/brain-imaging-traits/submitter.sh "$item"
done