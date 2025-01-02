"""Concatenate adata.
"""


import scanpy as sc
from tqdm import tqdm

from pathlib import Path


if __name__ == '__main__':
    
    data_dir = './outs/adata/qc/'
    out_file = './outs/adata/concatenated.h5ad'
    
    # load adatas
    print('loading data...')
    adatas = {}
    for f in tqdm(list(Path(data_dir).glob('./*.h5ad'))):
        sample = f.stem
        adatas[sample] = sc.read_h5ad(f)
        adatas[sample].var_names = list(adatas[sample].var['ens_id'])  # set var_names to gene ids

    # concatenate
    print('concatenating...')
    adata = sc.concat(adatas, label='sample', merge='same')
    adata.obs_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=10)  # filter genes
    sc.pp.filter_genes(adata, min_counts=5)  # filter genes
    adata.var_names = list(adata.var['symbol'])  # set var_names back to gene symbols
    adata.var_names_make_unique()
    print(adata)

    # save
    adata.write(out_file)

    print('Done!')
