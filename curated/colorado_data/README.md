# Colorado genomic data staging area

This directory is reserved for the Colorado genomic inputs used by the curated manuscript workflow.

For a full regeneration from raw genotypes, stage the files expected by `01_GATHER_GENDRIVER_DATA.ipynb` here and set that notebook's `GENDIR` variable to `./colorado_data`.

Expected/useful files include:

- `genomicdata.csv` — full genotype matrix used by the original v1 data-assembly notebook.
- `genomicdataheader.csv` — one genomic column/header name per row; used by `../scripts/find_phenotype_genomic_headers.py`.
- `biological_drivers.csv` — baseline biological-driver list used by the original v1 assembly workflow, if reconstructing the historical 171-locus matrix from source data.

The full SNP/genotype data file can be added here later. No patient-level genomic data are committed in this placeholder directory at present.
