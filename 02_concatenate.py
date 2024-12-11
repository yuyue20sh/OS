import scanpy as sc
from pathlib import Path
from tqdm import tqdm


if __name__ == '__main__':
    
    data_dir = './outs/adata/qc/'
    out_dir = './outs/adata/'
    
    # load adatas
    adatas = {}
    for f in tqdm(Path(data_dir).glob('./*.h5ad'), total=len(list(Path(data_dir).glob('./*.h5ad')))):
        id = f.stem
        adatas[id] = sc.read_h5ad(f)
        adatas[id].var_names = list(adatas[id].var['gene_ids'])

    # concatenate
    adata = sc.concat(adatas, label='sample', merge='first')
    adata.obs_names_make_unique()
    adata.obs['dataset'] = 'new'
    adata.obs.loc[adata.obs['sample'].str.contains('BC'), 'dataset'] = 'nc'
    print(adata)

    # save
    adata.write('%s/concatenated.h5ad' % out_dir)

    print('Done!')
