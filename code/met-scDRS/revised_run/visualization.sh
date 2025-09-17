current_date=$(date +"%Y-%m-%d")

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'newL1' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-l1-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'newL2' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-l2-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'newL3' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-l3-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/revised_run/proportion-heatmap.R \
    --dir '/u/home/l/lixinzhe/project-geschwind/result/met-scDRS/dev-revised/' \
    --meta_data "/u/project/cluo/heffel/BICAN3/DATA/metadata_09122025.csv.gz" \
    --field 'fine_age_groups' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-revised-BICAN-mcg-fine-age-proportion.png" \
    --plot_type "proportion"

