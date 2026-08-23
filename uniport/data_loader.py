#!/usr/bin/env
"""
# Author: Kai Cao
# Modified from SCALEX
"""

import numpy as np
import torch
from torch.utils.data import BatchSampler, DataLoader, Dataset
from torch.utils.data import RandomSampler, SequentialSampler
from scipy.sparse import issparse

__all__ = ['SingleCellDataset', 'SingleCellDataset_vertical', 'load_data']


def _to_dense_f32(x):
    """Densify (if needed) and cast to float32.

    The training loop casts every mini batch with ``.float()`` anyway, so doing
    it once up front is numerically identical while halving the memory held by
    the dataset and removing a per-batch conversion.
    """
    if issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def _pad_genes(x, max_gene):
    """Right-pad a cell x gene block with zero columns up to ``max_gene``."""
    if max_gene is None or x.shape[1] >= max_gene:
        return x
    pad = np.zeros((x.shape[0], max_gene - x.shape[1]), dtype=np.float32)
    return np.hstack((x, pad))


class SingleCellDataset(Dataset):

    def __init__(self, data, batch):

        self.data = torch.as_tensor(np.ascontiguousarray(_to_dense_f32(data)))
        self.batch = torch.as_tensor(np.asarray(batch), dtype=torch.int64)
        self.shape = tuple(self.data.shape)

    def __len__(self):
        return self.data.shape[0]

    def __getitem__(self, idx):
        # ``idx`` is a whole batch of indices when the loader is built by
        # ``load_data`` (see the BatchSampler there); the single-index form is
        # kept working so the dataset stays usable on its own.
        index = torch.as_tensor(idx) if isinstance(idx, (list, np.ndarray)) else idx

        return self.data[index], self.batch[index], index


class SingleCellDataset_vertical(Dataset):

    def __init__(self, adatas):

        # One tensor per modality rather than one pre-concatenated copy. For a
        # dense float32 ``.X`` this shares memory with the AnnData array, so
        # nothing is duplicated (concatenating up front would double the
        # footprint of an already large vertical dataset), and the per-batch
        # torch.cat is still one vectorised op instead of one np.concatenate
        # per sample. Sparse matrices are densified here, which slicing
        # ``adata.X`` directly could not handle at all.
        self.data = [torch.as_tensor(_to_dense_f32(adata.X)) for adata in adatas]
        self.shape = (self.data[0].shape[0], sum(int(d.shape[1]) for d in self.data))

    def __len__(self):
        return self.data[0].shape[0]

    def __getitem__(self, idx):
        index = torch.as_tensor(idx) if isinstance(idx, (list, np.ndarray)) else idx

        return torch.cat([d[index] for d in self.data], dim=-1), index


def _make_loader(scdata, batch_size, drop_last, shuffle, num_workers):
    """DataLoader that gathers a whole mini batch in one indexing operation.

    ``batch_size=None`` turns off the automatic collation, so the dataset gets
    the full list of indices and can slice the underlying tensor once instead of
    materialising ``batch_size`` python objects and stacking them. The sampler
    is the same ``BatchSampler(RandomSampler(...))`` that ``DataLoader`` would
    build internally, so the mini batches are identical.
    """
    sampler = RandomSampler(scdata) if shuffle else SequentialSampler(scdata)

    # Note: no persistent_workers. It would avoid respawning workers each epoch,
    # but it also changes how many values the loader draws from the global torch
    # RNG per epoch, which would silently change the shuffling (and therefore
    # every downstream result) relative to previous releases.
    return DataLoader(
        scdata,
        batch_size=None,
        sampler=BatchSampler(sampler, batch_size=batch_size, drop_last=drop_last),
        num_workers=num_workers,
        pin_memory=torch.cuda.is_available(),
    )


def load_data(adatas, mode='h', use_rep=['X', 'X'], max_gene=None, adata_cm=None, use_specific=False, domain_name='domain_id', batch_size=256, \
    drop_last=True, shuffle=True, num_workers=4):

    '''
    Load data for training.

    Parameters
    ----------
    adatas
        A list of AnnData matrice.
    mode
        training mode. Choose between ['h', 'd', 'v'].
    use_rep
        use '.X' or '.obsm'.
    max_gene
        maximum number of genes of each adata in adatas.
    adata_cm
        adata with common genes of adatas.
    use_specific
        use dataset-specific genes.
    domain_name
        domain name of each adata in adatas.
    batch_size
        size of each mini batch for training.
    drop_last
        drop the last samples that not up to one batch.
    shuffle
        shuffle the data
    num_workers
        number parallel load processes according to cpu cores.

    Returns
    -------
    trainloader
        data loader for training
    testloader
        data loader for testing
    '''

    if mode == 'd':
        blocks = []
        batches = []
        for i, adata in enumerate(adatas):
            rep = adata.X if use_rep[i] == 'X' else adata.obsm[use_rep[i]]
            blocks.append(_pad_genes(_to_dense_f32(rep), max_gene))
            batches.append(np.asarray(adata.obs[domain_name].astype(int)))

        x = np.vstack(blocks) if len(blocks) > 1 else blocks[0]
        batches = np.concatenate(batches)

        scdata = SingleCellDataset(x, batches)

    elif mode == 'h':
        domains = adata_cm.obs[domain_name].cat.categories.tolist()

        if use_specific:
            blocks = []
            for i, adata in enumerate(adatas):
                adata_tmp = adata_cm[adata_cm.obs[domain_name] == domains[i]]

                x_c = _to_dense_f32(adata_tmp.X)
                x_s = _pad_genes(_to_dense_f32(adata.X), max_gene)

                blocks.append(np.hstack((x_c, x_s)))
            x = np.vstack(blocks) if len(blocks) > 1 else blocks[0]
        else:
            x = _to_dense_f32(adata_cm.X)

        scdata = SingleCellDataset(x, np.asarray(adata_cm.obs[domain_name].astype(int)))

    else:
        scdata = SingleCellDataset_vertical(adatas)

    # DataLoader for train and test
    trainloader = _make_loader(scdata, batch_size, drop_last, shuffle, num_workers)
    testloader = _make_loader(scdata, batch_size, False, False, 0)

    return trainloader, testloader
