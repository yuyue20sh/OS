"""Create adata from mtx files.
"""


from tqdm import tqdm

from pathlib import Path

from utils import io


if __name__ == "__main__":

    data_dir = './data/sc/'
    out_dir = './outs/adata/raw/'

    # read data
    adatas = {}
    for d in Path(data_dir).glob('./*'):
        info = d.stem.split('_')
        sample = '%s_%s' % (info[0], info[1]) if len(info) == 6 else info[0]
        print(sample)
        adatas[sample] = io.read_mtx('%s/matrix.mtx' % d)

    # save adata
    print('saving adata...')
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    for id, adata in tqdm(adatas.items()):
        adata.write('%s/%s.h5ad' % (out_dir, id))

    print('Done!')
