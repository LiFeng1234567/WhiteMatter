# WhiteMatter Analysis Pipeline
# Spatial transcriptomics pipeline for primate white matter heterogeneity analysis
# Workflow: Raw Stereo-seq data → raw2h5ad → cell2location → bin200 aggregation → Random Forest models
# Scripts

| Script | Purpose |
|--------|---------|
| `raw2h5ad` | Convert raw count tables to AnnData (.h5ad) |
| `cell2location` | Map snRNA-seq cell types to spatial bins |
| `process_bin200_intra.R` | Aggregate intra-cellular data into 200×200 pixel bins → Seurat object |
| `process_bin200_extra.R` | Aggregate extra-cellular data (egRNAs) into bins |
| `RandomForest_OpticTract.R` | Classify optic tract vs other tracts |
| `RandomForest_VisualFibers.R` | Classify visual fiber types (SS/OVF/Finfer/Tap) |
| `RandomForest_CallosalProj.R` | Predict callosal projection clusters (eg1/eg2/eg3) |

# Quick Run
# 1. Convert raw data
python raw2h5ad counts.txt output.h5ad positions.csv

# 2. Run cell2location (see cell2location docs)

# 3. Generate bin200 Seurat objects
Rscript process_bin200_intra.R sample_id
Rscript process_bin200_extra.R sample_id

# 4. Train random forest models
Rscript RandomForest_OpticTract.R
Rscript RandomForest_VisualFibers.R
Rscript RandomForest_CallosalProj.R

# Requirements

- R ≥ 4.2 (Seurat, randomForest, caret, pROC, MLmetrics)
- Python ≥ 3.9 (scanpy, cell2location)

# Data Availability

Raw data: [OEP00003224, private]. Access available upon request.


## License

CC BY-NC 4.0
