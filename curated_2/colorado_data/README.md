# Colorado raw-data placeholder

This directory is reserved for restricted/raw Colorado genotype inputs if the processed genomic matrix is rebuilt in the future. Raw genotype data should not be committed unless its data-use agreement explicitly permits that.

The manuscript analyses in `curated_2/` do **not** read from this directory. They start from the processed `../ILD_TOP_DRIVERS_DATA.csv`.

A future raw-data rebuild would typically place the source genotype matrix and its header/index material here (for example `genomicdata.csv` and/or `genomicdataheader.csv`) and use an explicit, auditable preprocessing script to recreate `ILD_TOP_DRIVERS_DATA.csv`.
