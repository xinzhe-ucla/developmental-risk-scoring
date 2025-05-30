### cell_type_by_developmental_risk.R #############################################################
# purpose: visualize how the associated risk change with developmental time as well as in a cell type specific manner

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(dplyr)
library(data.table)

# load in the score files:
meta_path <- "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv"

# load in the data:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta_path
    );

# label the '' as 'unknown'
meta$L1[meta$L1 == ''] = 'unknown'
meta$L2[meta$L2 == ''] = 'unknown'
meta$L3[meta$L3 == ''] = 'unknown'

# recode the fine_age:
meta$fine_age[meta$fine_age == '1m'] = '1mo'
meta$fine_age[meta$fine_age == '4m'] = '4mo'
meta$fine_age[meta$fine_age == '21'] = '21 years old'
meta$fine_age[meta$fine_age == '29'] = '29 years old'
meta$fine_age[meta$fine_age == '31'] = '31 years old'
meta$fine_age[meta$fine_age == '37'] = '37 years old'


# read data:
scDRS_directory = "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/"
score_files <- list.files(scDRS_directory, pattern = '\\.score.gz', full.names = TRUE);
risk_score <- vector('list', length = length(score_files));
names(risk_score) <- score_files;

# read into the empty list:
for (result in score_files) {
    met_scdrs_result <- data.frame(
        fread(
        file = result,
        sep = '\t',
        header = TRUE,
        data.table = FALSE
        ),
        row.names = 1)
    rownames(met_scdrs_result) = gsub('\\.allc.tsv.gz', '', rownames(met_scdrs_result))
    met_scdrs_result = met_scdrs_result[rownames(meta), ]
    met_scdrs_result$fdr = p.adjust(met_scdrs_result$pval, method = 'fdr')
    risk_score[[result]] = met_scdrs_result
    }

# simplify list names:
list_names <- gsub(scDRS_directory, '', score_files);
list_names <- gsub('/', '', list_names);
list_names <- gsub('\\.score.gz', '', list_names);

# rename the list names:
names(risk_score) <- list_names;

### PROCESS #######################################################################################
# subset by developmental age and summarize the mean and variance:
for (disease in names(risk_score)){
    # obtain the z score for different developmental time points:
    qc_score = cbind(risk_score[[disease]][rownames(meta), ], meta)
    zscore_by_age = qc_score %>% group_by(fine_age) %>% summarize(zscore_mean = mean(zscore), zscore_sd = sd(zscore))

    # obtain the z score grouped by both time and cell type:
    zscore_by_age_l1 = qc_score %>% group_by(fine_age, L1) %>% summarize(zscore_mean = mean(zscore), zscore_sd = sd(zscore))


    

    }


    # obtain check if the z socre vary by cell type:
    interactive_model = lm(zscore ~ L1*fine_age, data = qc_score)
    summary_table = summary(interactive_model)$coefficients
