import scanpy as sc
from tqdm import tqdm

from pathlib import Path


if __name__ == '__main__':
    
    data_dir = './outs/adata/qc/'
    out_file = './outs/adata/concatenated.h5ad'
    
    # load adatas
    adatas = {}
    for f in tqdm(Path(data_dir).glob('./*.h5ad'), total=len(list(Path(data_dir).glob('./*.h5ad')))):
        sample = f.stem
        adatas[sample] = sc.read_h5ad(f)
        adatas[sample].var_names = list(adatas[sample].var['ens_id'])  # set var_names to gene ids

    # concatenate
    adata = sc.concat(adatas, label='sample')
    adata.obs_names_make_unique()
    adata.var_names = list(adata.var['symbol'])  # set var_names back to gene symbols
    adata.var_names_make_unique()
    print(adata)

    # save
    adata.write(out_file)

    print('Done!')
