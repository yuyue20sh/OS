import scanpy as sc
import scipy.io as sio

from pathlib import Path


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

    adata_file = './outs/adata/small_annotated_lvl1.h5ad'
    mtx_dir = './outs/mtx/small_annotated_lvl1/'
    annotation_file = './outs/infercnv/inputs/annotation.tsv'

    adata = sc.read_h5ad(adata_file)

    # counts
    adata.X = adata.layers['counts']
    h5ad2mtx(adata, mtx_dir)

    # annotation
    annotation = adata.obs
    annotation['annotation'] = annotation['cell_type_lvl1'].astype(str)

    msc_cells = annotation['annotation'] == 'MSC Derived'
    annotation.loc[msc_cells, 'annotation'] = ['_'.join(i) for i in annotation.loc[msc_cells, ['annotation', 'sample']].values]

    # save
    annotation['annotation'].to_csv(annotation_file, sep='\t', index=True, header=False)

    print('Done!')
