current_date=$(date +"%Y-%m-%d")

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/proportion-heatmap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv" \
    --field 'L1' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-initial-BICAN-mcg-l1-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/proportion-heatmap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv" \
    --field 'L2' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-initial-BICAN-mcg-l2-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/proportion-heatmap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv" \
    --field 'L3' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-initial-BICAN-mcg-L3-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/proportion-heatmap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv" \
    --field 'fine_age' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-initial-BICAN-mcg-fine-age-proportion.png" \
    --plot_type "proportion"

Rscript /u/home/l/lixinzhe/project-github/developmental-risk-scoring/code/met-scDRS/initial-run/proportion-heatmap.R \
    --dir "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/single_cell_baseline/" \
    --meta_data "/u/home/l/lixinzhe/project-cluo/data/2025-05-02-combined-meta-QCed.csv" \
    --field 'Region' \
    --cutoff 0.1 \
    --out "/u/home/l/lixinzhe/project-geschwind/plot/${current_date}-initial-BICAN-mcg-region-proportion.png" \
    --plot_type "proportion"
