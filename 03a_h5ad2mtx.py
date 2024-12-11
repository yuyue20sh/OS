import scanpy as sc
from pathlib import Path
import scipy.io as sio
from tqdm import tqdm


def h5ad2mtx(adata, out_dir):
    '''
    convert anndata object to mtx format
    output files: matrix.mtx, barcodes.tsv, features.tsv, meta.tsv

    Args:
        adata: AnnData object, adata to convert
        out_dir: str, path to the output directory

    Returns:
        None

    '''
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    matrix = adata.X.T
    barcodes = adata.obs.index.to_frame()
    features = adata.var[['gene_ids', 'gene_names']]
    meta = adata.obs

    sio.mmwrite('%s/matrix.mtx' % out_dir, matrix)
    barcodes.to_csv('%s/barcodes.tsv' % out_dir, sep='\t', header=False, index=False)
    features.to_csv('%s/features.tsv' % out_dir, sep='\t', header=False, index=False)
    meta.to_csv('%s/meta.tsv' % out_dir, sep='\t', header=True, index=True)

    return


if __name__ == '__main__':

    data_file = './outs/adata/concatenated.h5ad'
    out_dir = './outs/mtx/concatenated/'

    adata = sc.read_h5ad(data_file)
    adatas = {id: adata[adata.obs['sample'] == id] for id in adata.obs['sample'].unique()}

    for id, sample_adata in tqdm(adatas.items()):
        h5ad2mtx(sample_adata, '%s/%s' % (out_dir, id))

    print('Done!')
