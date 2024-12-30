"""Integrate samples using Harmony.
"""


import scanpy as sc


if __name__ == '__main__':
    
    adata_file = './outs/adata/concatenated.h5ad'
    out_file = './outs/adata/harmony.h5ad'
    
    # read data
    print('loading data...')
    adata = sc.read_h5ad(adata_file)
    print(adata)

    # preprocess
    print('preprocessing...')
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, exclude_highly_expressed=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='sample')
    sc.pp.pca(adata)

    # harmony
    print('running harmony...')
    sc.external.pp.harmony_integrate(adata, key='sample')

    adata.write(out_file)

    print('Done!')
