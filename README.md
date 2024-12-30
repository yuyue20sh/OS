### Environment

```shell
conda create -n sc python=3.12 -y
conda activate sc
pip install 'scanpy[leiden]'
pip install scrublet
pip install harmonypy
pip install notebook

conda create -n r433 python=3.10 -y
conda activate r433
conda install conda-forge::r-base=4.3.3 -y
conda install bioconda::bioconductor-biomart -y
conda install conda-forge::r-seurat -y
conda install conda-forge::r-devtools -y
conda install conda-forge::r-biocmanager -y
conda install conda-forge::r-rjags -y
R
> devtools::install_github("navinlabcode/copykat")
> BiocManager::install("infercnv")
> q()
```

### Run

```shell
conda activate sc
python 00_create_adata.py
python 01_qc.py
python 02_concatenate.py
python 03a_h5ad2mtx.py

conda activate r433
Rscript 03b_copykat.R > ./logs/copykat.log

conda activate sc
jupyter execute 04_integrate.ipynb --inplace
jupyter execute 05_annotate_lvl1_small.ipynb --inplace

python 06a_prepare_data.py
wget https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt -P ./outs/infercnv/inputs/
conda activate r433
Rscript 06b_infercnv_small.R > ./logs/infercnv.log
```

### Results

#### Integration

1. Harmony was good.
2. Batch effect was removed using harmony while some clusters were predominently found in only one sample.

#### Copykat

1. BC3, DLJM44, WBXM16, WBXM16_2, WYQM12 and XZHM13 had 'low confidence in classification' (as well as 'WARNING! NOT CONVERGENT!')
2. ZCLM12 had warning 'WARNING! NOT CONVERGENT!'
