#!/usr/bin/env 
"""
# Author: Kai Cao
# Modified from SCALEX
"""

import numpy as np
import scipy
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score

__all__ = ['batch_entropy_mixing_score', 'silhouette', 'label_transfer']


def batch_entropy_mixing_score(data, batches, n_neighbors=100, n_pools=100, n_samples_per_pool=100, n_jobs=-1):
    """
    Calculate batch entropy mixing score
    
    Algorithm
    ---------
        * 1. Calculate the regional mixing entropies at the location of 100 randomly chosen cells from all batches
        * 2. Define 100 nearest neighbors for each randomly chosen cell
        * 3. Calculate the mean mixing entropy as the mean of the regional entropies
        * 4. Repeat above procedure for 100 iterations with different randomly chosen cells.
    
    Parameters
    ----------
    data
        np.array of shape nsamples x nfeatures.
    batches
        batch labels of nsamples.
    n_neighbors
        The number of nearest neighbors for each randomly chosen cell. By default, n_neighbors=100.
    n_samples_per_pool
        The number of randomly chosen cells from all batches per iteration. By default, n_samples_per_pool=100.
    n_pools
        The number of iterations with different randomly chosen cells. By default, n_pools=100.
    n_jobs
        Number of parallel jobs for the nearest-neighbor search. By default, n_jobs=-1 (all cores).

    Returns
    -------
    Batch entropy mixing score
    """
#     print("Start calculating Entropy mixing score")
    def entropy(counts):
        # `counts[i]` is how many of the cell's neighbours belong to batch i, so
        # the batch frequencies come straight out of a bincount instead of one
        # full pass over the neighbour labels per batch.
        p = counts / counts.sum()
        a = 0
        for i in range(N_batches):
            a = a + p[i]/P[i]
        entropy = 0
        for i in range(N_batches):
            adapt_p = (p[i]/P[i])/a
            entropy = entropy - adapt_p*np.log(adapt_p+10**-8)
        return entropy

    n_neighbors = min(n_neighbors, len(data) - 1)
    nne = NearestNeighbors(n_neighbors=1 + n_neighbors, n_jobs=n_jobs)
    nne.fit(data)
    kmatrix = nne.kneighbors_graph(data) - scipy.sparse.identity(data.shape[0])

    score = 0
    batches = np.asarray(batches)
    batches_, batch_codes = np.unique(batches, return_inverse=True)
    N_batches = len(batches_)
    if N_batches < 2:
        raise ValueError("Should be more than one cluster for batch mixing")
    P = np.bincount(batch_codes, minlength=N_batches) / len(batch_codes)

    for t in range(n_pools):
        indices = np.random.choice(np.arange(data.shape[0]), size=n_samples_per_pool)
        # One nonzero() pass per pool. The original called it twice for each of
        # the `n_samples_per_pool` cells, i.e. 2 x n_samples_per_pool full scans
        # of the same sliced matrix on every iteration.
        rows, cols = kmatrix[indices].nonzero()
        # nonzero() emits entries in row-major order, so each row is one slice.
        row_start = np.searchsorted(rows, np.arange(n_samples_per_pool), side='left')
        row_end = np.searchsorted(rows, np.arange(n_samples_per_pool), side='right')
        neighbor_codes = batch_codes[cols]
        score += np.mean([
            entropy(np.bincount(neighbor_codes[row_start[i]:row_end[i]], minlength=N_batches))
            for i in range(n_samples_per_pool)
        ])
    Score = score / float(n_pools)
    return Score / float(np.log2(N_batches))


def silhouette(
        X,
        cell_type,
        metric='euclidean',
        scale=True
):
    """
    Wrapper for sklearn silhouette function values range from [-1, 1] with
        1 being an ideal fit
        0 indicating overlapping clusters and
        -1 indicating misclassified cells
    By default, the score is scaled between 0 and 1. This is controlled `scale=True`

    :param group_key: key in adata.obs of cell labels
    :param embed: embedding key in adata.obsm, default: 'X_pca'
    :param scale: default True, scale between 0 (worst) and 1 (best)
    """
    asw = silhouette_score(
        X,
        cell_type,
        metric=metric
    )
    if scale:
        asw = (asw + 1) / 2
    return asw

def label_transfer(ref, query, rep='latent', label='celltype'):
    """
    Label transfer
    
    Parameters
    -----------
    ref
        reference containing the projected representations and labels
    query
        query data to transfer label
    rep
        representations to train the classifier. Default is `latent`
    label
        label name. Defautl is `celltype` stored in ref.obs
    
    Returns
    --------
    transfered label
    """

    from sklearn.neighbors import KNeighborsClassifier
    
    X_train = ref.obsm[rep]
    y_train = ref.obs[label]
    X_test = query.obsm[rep]
    
    knn = KNeighborsClassifier(n_jobs=-1).fit(X_train, y_train)
    y_test = knn.predict(X_test)
    
    return y_test