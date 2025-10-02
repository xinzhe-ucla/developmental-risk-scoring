### covaraite-process.py ##########################################################################
# load in the meta data:

meta = read.csv('/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz')
rownames(meta) = meta$ID

scDRS.directory = '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/'
mdd = data.frame(
    data.table::fread(
        file =  paste0(scDRS.directory, "/PASS_MDD_Howard2019.score.gz"),
        sep = '\t',
        header = TRUE,
        data.table = FALSE
        ),
        row.names = 1
    )

# subset to the same rows :
# format the meta data:
new_name = gsub('-0-0-0', '', rownames(mdd))
new_name = gsub('-1-0$', '', new_name)
new_name = gsub('-1$', '', new_name)
new_name = gsub('-1-0-0$', '', new_name)
stopifnot(all(new_name %in% rownames(meta)))

# get the common cells:
meta = meta[new_name, ]
rownames(meta) = rownames(mdd)

covariate = data.frame(
    cell = rownames(meta),
    const= 1,
    global= meta[,'mCG.CG_global']
)

# output:
write.table(
    covariate,
    row.names = FALSE,
    col.names = TRUE,
    sep = '\t',
    quote = FALSE,
    file = '/u/project/cluo/lixinzhe/data/BICAN3/unnormalized_genebody_blacklist_allgenes_merged_cov.cov'
    )