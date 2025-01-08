"""Prepare data for copykat.
"""


import scanpy as sc
from tqdm import tqdm

from utils import io


if __name__ == '__main__':

    adata_file = './outs/adata/annotated_lvl1.h5ad'
    out_dir = './outs/mtx/annotated_lvl1_split_samples/'

    print('loading data...')
    adata = sc.read_h5ad(adata_file)
    adata.X = adata.layers['counts']
    adatas = {sample: adata[adata.obs['sample'] == sample] for sample in adata.obs['sample'].unique()}

    print('saving mtx...')
    for sample, sample_adata in tqdm(adatas.items()):
        # counts
        io.h5ad2mtx(sample_adata, '%s/%s' % (out_dir, sample))

    print('Done!')
