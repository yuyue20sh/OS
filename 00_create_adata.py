import scanpy as sc
import pandas as pd
from pathlib import Path
from tqdm import tqdm


def mtx2adata(mtx_dir, make_unique=True):
    '''
    read matrix.mtx, barcodes.tsv, and features.tsv files in mtx_dir, and return an AnnData object

    Args:
        mtx_dir: str, path to the directory with matrix.mtx, barcodes.tsv, and features.tsv files
        make_unique: bool, whether to make obs_names and var_names unique
    
    Returns:
        matrix: AnnData object

    '''
    matrix_file = '%s/matrix.mtx' % mtx_dir
    barcodes_file = '%s/barcodes.tsv' % mtx_dir
    features_file = '%s/features.tsv' % mtx_dir if Path('%s/features.tsv' % mtx_dir).exists() else '%s/genes.tsv' % mtx_dir

    matrix = sc.read_mtx(matrix_file).T
    barcodes = pd.read_csv(barcodes_file, sep='\t', header=None, index_col=None)
    features = pd.read_csv(features_file, sep='\t', header=None, index_col=None)

    obs_names = barcodes[0].values
    var_names = features[1].values
    matrix.obs_names = obs_names
    matrix.var_names = var_names

    matrix.var['gene_ids'] = features[0].values
    matrix.var['gene_names'] = features[1].values

    if make_unique:
        matrix.obs_names_make_unique()
        matrix.var_names_make_unique()

    return matrix


if __name__ == "__main__":

    data_dir = './data/sc/'
    out_dir = './outs/adata/raw/'

    # load data and create adata
    print('loading data...')
    adatas = {}
    for d in Path(data_dir).glob('./*'):
        id = d.stem.split('_')
        id = '%s_%s' %(id[0], id[1]) if len(id) == 6 else id[0]   # get id from the directory name
        print(id)
        adatas[id] = mtx2adata(d)
        print(adatas[id])

    # save adata
    print('saving adata...')
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    for id, adata in tqdm(adatas.items()):
        adata.write('%s/%s.h5ad' % (out_dir, id))

    print('Done!')
