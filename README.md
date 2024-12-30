## Methods

### Environment

```shell
conda create -n sc python=3.12 -y
conda activate sc
pip install 'scanpy[leiden]'
pip install scrublet harmonypy notebook
```

### Run

```shell
conda activate sc
python 00_create_adata.py
python 01_qc.py
python 02_concatenate.py
python 03_harmony.py
jupyter execute 04_annotate_lvl1.ipynb --inplace
```



## Results

### Quality control

**QC Strategy:**

1. Cells with < 500 genes were excluded.
2. Genes expressed in < 0 cells were exclude (not filtering genes before merging).
3. Percentage of counts in mitochondrial genes was calculated in each cell.
4. Doublet scores were calculated using Scrublet.
5. Leiden clustering was done.
6. Cells with percentage of counts in mitochondrial genes greater than 2
   standard deviations from the mean percentage of counts in mitochondrial
   genes of their leiden clusters were removed.
7. The predicted doublets were removed.
8. Leiden clusters with mean doublet score greater than 2 standard deviations
   from the mean doublet score of all cells the were excluded.

**After concatenation,** genes with ≥ 10 UMI in at least 10 cells were included for downstream analysis.

**Finally,** 128637 cells and 32982 genes passed QC.

### Integration

1. Harmony was good.
2. Batch effect was removed using Harmony while some clusters were predominently found in one sample. **[TODO: exclude?]**

### Annotation

Add cell type annotation using Leiden.
