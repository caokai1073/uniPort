#!/usr/bin/env 
"""
# Author: Kai Cao
"""

import torch
import numpy as np

import os
import scanpy as sc
from anndata import AnnData
import scipy
import sklearn
import pandas as pd
from scipy.sparse import issparse

from .model.vae import VAE
from .model.utils import EarlyStopping
from .logger import create_logger
from .data_loader import load_data
from .metrics import *

from anndata import AnnData
from sklearn.preprocessing import MaxAbsScaler

from glob import glob

np.warnings.filterwarnings('ignore')
DATA_PATH = os.path.expanduser("~")+'/.uniport/'
CHUNK_SIZE = 20000

def read_mtx(path):
    """
    Read mtx format data folder including: 
    
        * matrix file: e.g. count.mtx or matrix.mtx or their gz format
        * barcode file: e.g. barcode.txt
        * feature file: e.g. feature.txt
        
    Parameters
    ----------
    path
        the path store the mtx files  
        
    Return
    ------
    AnnData
    """
    for filename in glob(path+'/*'):
        if ('count' in filename or 'matrix' in filename or 'data' in filename) and ('mtx' in filename):
            adata = sc.read_mtx(filename).T
    for filename in glob(path+'/*'):
        if 'barcode' in filename:
            barcode = pd.read_csv(filename, sep='\t', header=None).iloc[:, -1].values
            adata.obs = pd.DataFrame(index=barcode)
        if 'gene' in filename or 'peaks' in filename:
            gene = pd.read_csv(filename, sep='\t', header=None).iloc[:, -1].values
            adata.var = pd.DataFrame(index=gene)
        elif 'feature' in filename:
            gene = pd.read_csv(filename, sep='\t', header=None).iloc[:, 1].values
            adata.var = pd.DataFrame(index=gene)
             
    return adata


def load_file(path):  
    """
    Load single cell dataset from file
    
    Parameters
    ----------
    path
        the path store the file
        
    Return
    ------
    AnnData
    """
    if os.path.exists(DATA_PATH+path+'.h5ad'):
        adata = sc.read_h5ad(DATA_PATH+path+'.h5ad')
    elif os.path.isdir(path): # mtx format
        adata = read_mtx(path)
    elif os.path.isfile(path):
        if path.endswith(('.csv', '.csv.gz')):
            adata = sc.read_csv(path).T
        elif path.endswith(('.txt', '.txt.gz', '.tsv', '.tsv.gz')):
            df = pd.read_csv(path, sep='\t', index_col=0).T
            adata = AnnData(df.values, dict(obs_names=df.index.values), dict(var_names=df.columns.values))
        elif path.endswith('.h5ad'):
            adata = sc.read_h5ad(path)
    else:
        raise ValueError("File {} not exists".format(path))
        
    if not issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)

    return adata

def tfidf(X, n_components, binarize=True, random_state=0):
    from sklearn.feature_extraction.text import TfidfTransformer
    
    sc_count = np.copy(X)
    if binarize:
        sc_count = np.where(sc_count < 1, sc_count, 1)
    
    tfidf = TfidfTransformer(norm='l2', sublinear_tf=True)
    normed_count = tfidf.fit_transform(sc_count)

    lsi = sklearn.decomposition.TruncatedSVD(n_components=n_components, random_state=random_state)
    lsi_r = lsi.fit_transform(normed_count)
    
    X_lsi = lsi_r[:,1:]

    return X_lsi
    
def TFIDF_LSI(adata, n_comps=50, binarize=True, random_state=0):
    '''
    Computes LSI based on a TF-IDF transformation of the data from MultiMap. Putative dimensionality 
    reduction for scATAC-seq data. Adds an ``.obsm['X_lsi']`` field to the object it was ran on. 
    
    Input
    -----
    adata : ``AnnData``
        The object to run TFIDF + LSI on. Will use ``.X`` as the input data.
    n_comps : ``int``
        The number of components to generate. Default: 50
    binarize : ``bool``
        Whether to binarize the data prior to the computation. Often done during scATAC-seq 
        processing. Default: True
    random_state : ``int``
        The seed to use for randon number generation. Default: 0
    '''
    
    #this is just a very basic wrapper for the non-adata function
    if scipy.sparse.issparse(adata.X):
        adata.obsm['X_lsi'] = tfidf(adata.X.todense(), n_components=n_comps, binarize=binarize, random_state=random_state)
    else:
        adata.obsm['X_lsi'] = tfidf(adata.X, n_components=n_comps, binarize=binarize, random_state=random_state)

def filter_data(
        adata: AnnData,
        min_features: int = 0, 
        min_cells: int = 0,     
        log=None
    ):
    """
    Filter cells and genes
    
    Parameters
    ----------
    adata
        An AnnData matrice of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    min_features
        Filtered out cells that are detected in less than n genes. Default: 0.
    min_cells
        Filtered out genes that are detected in less than n cells. Default: 0.
        
    """
    

    if log: log.info('Filtering cells')
    sc.pp.filter_cells(adata, min_genes=min_features)
    
    if log: log.info('Filtering features')
    sc.pp.filter_genes(adata, min_cells=min_cells)

def batch_scale(adata, use_rep='X', chunk_size=CHUNK_SIZE):
    """
    Batch-specific scale data
    
    Parameters
    ----------
    adata
        AnnData
    use_rep
        use '.X' or '.obsm'
    chunk_size
        chunk large data into small chunks
    
    """
    for b in adata.obs['source'].unique():
        idx = np.where(adata.obs['source']==b)[0]
        if use_rep == 'X':
            scaler = MaxAbsScaler(copy=False).fit(adata.X[idx])
            for i in range(len(idx)//chunk_size+1):
                adata.X[idx[i*chunk_size:(i+1)*chunk_size]] = scaler.transform(adata.X[idx[i*chunk_size:(i+1)*chunk_size]])
        else:
            scaler = MaxAbsScaler(copy=False).fit(adata.obsm[use_rep][idx])
            for i in range(len(idx)//chunk_size+1):
                adata.obsm[use_rep][idx[i*chunk_size:(i+1)*chunk_size]] = scaler.transform(adata.obsm[use_rep][idx[i*chunk_size:(i+1)*chunk_size]])

def get_prior(celltype1, celltype2, alpha=2):

    """
    Create a prior correspondence matrix according to cell labels
    
    Parameters
    ----------
    celltype1
        cell labels of dataset X
    celltype2
        cell labels of dataset Y
    alpha
        the confidence of label, ranges from (1, inf). Higher alpha means better confidence. Default: 2.0

    Return
    ------
    torch.tensor
        a prior correspondence matrix between cells
    """

    Couple = alpha*torch.ones(len(celltype1), len(celltype2))
    
    for i in set(celltype1):
        index1 = np.where(celltype1==i)
        if i in set(celltype2):
            index2 = np.where(celltype2==i)
            for j in index1[0]:
                Couple[j, index2[0]]=1/alpha

    return Couple

def label_reweight(celltype):

    """
    Reweight labels to make all cell types share the same total weight 
    
    Parameters
    ----------
    celltype
        cell labels

    Return
    ------
    torch.tensor
        a vector of weights of cells 
    """

    n = len(celltype)
    unique, count = np.unique(celltype, return_counts=True)
    p = torch.zeros(n,1)

    for i in range(n):
        idx = np.where(unique==celltype[i])[0]
        tmp = 1/(len(unique)*count[idx])
        p[i] = torch.from_numpy(tmp)

    weights = p * len(celltype)

    return weights

# @profile
def Run(
        adatas=None,     
        adata_cm=None,   
        mode='h',
        lambda_s=0.5,
        lambda_recon=1.0,
        lambda_kl=0.5,
        lambda_ot=1.0,
        iteration=30000,
        ref_id=None,    
        save_OT=False,
        use_rep=['X', 'X'],
        out='latent',
        label_weight=None,
        reg=0.1,
        reg_m=1.0,
        batch_size=256, 
        lr=2e-4, 
        enc=None,
        gpu=0, 
        prior=None,
        loss_type='BCE',
        outdir='output/', 
        input_id=0,
        pred_id=1,
        seed=124, 
        num_workers=4,
        patience=30,
        batch_key='domain_id',
        source_name='source',
        model_info=False,
        verbose=False,
    ):

    """
    Run data integration
    
    Parameters
    ----------
    adatas
        List of AnnData matrices, e.g. [adata1, adata2].
    adata_cm
        AnnData matrices containing common genes.
    mode
        Choose from ['h', 'v', 'd']
        If 'h', integrate data with common genes (Horizontal integration)
        If 'v', integrate data profiled from the same cells (Vertical integration)
        If 'd', inetrgate data without common genes (Diagonal integration)
        Default: 'h'.
    lambda_s
        Balanced parameter for common and specific genes. Default: 0.5
    lambda_recon: 
        Balanced parameter for reconstruct term. Default: 1.0
    lambda_kl: 
        Balanced parameter for KL divergence. Default: 0.5
    lambda_ot:
        Balanced parameter for OT. Default: 1.0
    iteration
        Max iterations for training. Training one batch_size samples is one iteration. Default: 30000
    ref_id
        Id of reference dataset. Default: None
    save_OT
        If True, output a global OT plan. Need more memory. Default: False
    use_rep
        Use '.X' or '.obsm'. For mode='d' only.
        If use_rep=['X','X'], use 'adatas[0].X' and 'adatas[1].X' for integration.
        If use_rep=['X','X_lsi'],  use 'adatas[0].X' and 'adatas[1].obsm['X_lsi']' for integration.
        If use_rep=['X_pca', 'X_lsi'], use 'adatas[0].obsm['X_pca']' and 'adatas[1].obsm['X_lsi']' for integration.
        Default: ['X','X']
    out
        Output of uniPort. Choose from ['latent', 'project', 'predict'].
        If out=='latent', train the network and output cell embeddings.
        If out=='project', project data into the latent space and output cell embeddings. 
        If out=='predict', project data into the latent space and output cell embeddings through a specified decoder.
        Default: 'latent'. 
    label_weight
        Prior-guided weighted vectors. Default: None
    reg:
        Entropy regularization parameter in OT. Default: 0.1
    reg_m:
        Unbalanced OT parameter. Larger values means more balanced OT. Default: 1.0
    batch_size
        Number of samples per batch to load. Default: 256
    lr
        Learning rate. Default: 2e-4
    enc
        Structure of encoder
    gpu
        Index of GPU to use if GPU is available. Default: 0
    prior
        Prior correspondence matrix. Default: None
    loss_type
        type of loss. 'BCE', 'MSE' or 'L1'. Default: 'BCE'
    outdir
        Output directory. Default: 'output/'
    input_id
        Only used when mode=='d' and out=='predict' to choose a encoder to project data. Default: 0
    pred_id
        Only used when out=='predict' to choose a decoder to predict data. Default: 1
    seed
        Random seed for torch and numpy. Default: 124
    patience
        early stopping patience. Default: 10
    batch_key
        Name of batch in AnnData. Default: domain_id
    source_name
        Name of source in AnnData. Default: source
    rep_celltype
        Names of cell-type annotation in AnnData. Default: 'cell_type'   
    umap
        If True, perform UMAP for visualization. Default: False
    model_info
        If True, show structures of encoder and decoders.
    verbose
        Verbosity, True or False. Default: False
    assess
        If True, calculate the entropy_batch_mixing score and silhouette score to evaluate integration results. Default: False
    show
        If True, show the UMAP visualization of latent space. Default: False

    Returns
    -------
    adata.h5ad
        The AnnData matrice after integration. The representation of the data is stored at adata.obsm['latent'], adata.obsm['project'] or adata.obsm['predict'].
    checkpoint
        model.pt contains the variables of the model and config.pt contains the parameters of the model.
    log.txt
        Records model parameters.
    umap.pdf 
        UMAP plot for visualization if umap=True.
    """

    if mode == 'h' and adata_cm is None:
        raise AssertionError('adata_cm is needed when mode == "h"!')

    if mode not in ['h', 'd', 'v']:
        raise AssertionError('mode must be "h", "v" or "d" ')

    if adatas is None and adata_cm is None:
        raise AssertionError('at least one of adatas and adata_cm should be given!')

    np.random.seed(seed) # seed
    torch.manual_seed(seed)

    if torch.cuda.is_available(): # cuda device
        device='cuda'
        torch.cuda.set_device(gpu)
    else:
        device='cpu'
    
    print('Device:', device)

    outdir = outdir+'/'
    os.makedirs(outdir+'/checkpoint', exist_ok=True)
    log = create_logger('', fh=outdir+'log.txt')

    use_specific=True

    # split adata_cm to adatas
    if adatas is None:  
        use_specific = False
        _, idx = np.unique(adata_cm.obs[source_name], return_index=True)
        batches = adata_cm.obs[source_name][np.sort(idx)]
        flagged = []
        for batch in batches:
            flagged.append(adata_cm[adata_cm.obs[source_name]==batch].copy())
        adatas = flagged

    n_domain = len(adatas)

    # give reference datasets
    if ref_id is None:  
        ref_id = n_domain-1

    tran = {}
    num_cell = []
    num_gene = []

    for i, adata in enumerate(adatas):
        if use_rep[i]=='X':
            num_cell.append(adata.X.shape[0])
            num_gene.append(adata.X.shape[1])
        else:
            num_cell.append(adata.obsm[use_rep[i]].shape[0])
            num_gene.append(adata.obsm[use_rep[i]].shape[1])


    # training
    if out == 'latent':

        for i, adata in enumerate(adatas):
            print('Dataset {}:'.format(i), adata.obs[source_name][0])
            print(adata)

        print('Reference dataset is dataset {}'.format(ref_id))
        print('\n')

        if adata_cm is not None:
            print('Data with common HVG')
            print(adata_cm)
            print('\n')

        if save_OT:
            for i in range(n_domain):
                if i != ref_id:
                    ns = num_cell[i]
                    nt = num_cell[ref_id]
                    tran_tmp = np.ones((ns, nt)) / (ns * nt)
                    tran[i] = tran_tmp.astype(np.float32)

                    print('Size of transport plan between datasets {} and {}:'.format(i, ref_id), np.shape(tran[i]))

        trainloader, testloader = load_data(
            adatas=adatas, 
            mode=mode,
            use_rep=use_rep,
            num_cell=num_cell,
            max_gene=max(num_gene), 
            adata_cm=adata_cm,
            use_specific=use_specific, 
            domain_name=batch_key,
            batch_size=batch_size,
            num_workers=num_workers
        )

        early_stopping = EarlyStopping(patience=patience, checkpoint_file=outdir+'/checkpoint/model.pt')
        
        # encoder structure
        if enc is None:
            enc = [['fc', 1024, 1, 'relu'], ['fc', 16, '', '']]      

        # decoder structure
        dec = {} 
        if mode == 'd':     
            for i in range(n_domain):          
                dec[i] = [['fc', num_gene[i], 1, 'sigmoid']]

        elif mode == 'h':      
            num_gene.append(adata_cm.X.shape[1]) 
            dec[0] = [['fc', num_gene[n_domain], n_domain, 'sigmoid']]  # common decoder
            if use_specific: 
                for i in range(1, n_domain+1):
                    dec[i] = [['fc', num_gene[i-1], 1, 'sigmoid']]   # dataset-specific decoder

        else:
            for i in range(n_domain):
                dec[i] = [['fc', num_gene[i], 1, 'sigmoid']]    # dataset-specific decoder

        # init model
        model = VAE(enc, dec, ref_id=ref_id, n_domain=n_domain, mode=mode)

        if model_info:
            log.info('model\n'+model.__repr__())

        model.fit(
            trainloader, 
            tran,
            num_cell,
            num_gene,
            mode=mode,
            label_weight=label_weight,
            Prior=prior,
            save_OT=save_OT,
            use_specific=use_specific,
            lambda_s=lambda_s,
            lambda_recon=lambda_recon,
            lambda_kl=lambda_kl,
            lambda_ot=lambda_ot,
            reg=reg,
            reg_m=reg_m,
            lr=lr, 
            max_iteration=iteration, 
            device=device, 
            early_stopping=early_stopping, 
            verbose=verbose,
            loss_type=loss_type,
        )
        torch.save({'enc':enc, 'dec':dec, 'n_domain':n_domain, 'ref_id':ref_id, 'num_gene':num_gene}, outdir+'/checkpoint/config.pt')     


    # project or predict
    else:
        state = torch.load(outdir+'/checkpoint/config.pt')
        enc, dec, n_domain, ref_id, num_gene = state['enc'], state['dec'], state['n_domain'], state['ref_id'], state['num_gene']
        model = VAE(enc, dec, ref_id=ref_id, n_domain=n_domain, mode=mode)
        model.load_model(outdir+'/checkpoint/model.pt')
        model.to(device)
        
        _, testloader = load_data(
            adatas=adatas, 
            max_gene=max(num_gene), 
            num_cell=num_cell,
            adata_cm=adata_cm, 
            domain_name=batch_key,
            batch_size=batch_size, 
            mode=mode
        )

    if mode == 'v':
        adatas[0].obsm[out] = model.encodeBatch(testloader, num_gene, pred_id=pred_id, device=device, mode=mode, out=out)
        return adatas[0]

    elif mode == 'd':
        if out == 'latent' or out == 'project':
            for i in range(n_domain):
                adatas[i].obsm[out] = model.encodeBatch(testloader, num_gene, batch_id=i, device=device, mode=mode, out=out)
            for i in range(n_domain-1):
                adata_concat = adatas[i].concatenate(adatas[i+1])
        elif out == 'predict':
            adatas[0].obsm[out] = model.encodeBatch(testloader, num_gene, batch_id=input_id, pred_id=pred_id, device=device, mode=mode, out=out)

    elif mode == 'h':
        if out == 'latent' or out == 'project':
            adata_cm.obsm[out] = model.encodeBatch(testloader, num_gene, device=device, mode=mode, out=out) # save latent rep
        elif out == 'predict':
            adata_cm.obsm[out] = model.encodeBatch(testloader, num_gene, pred_id=pred_id, device=device, mode=mode, out=out)

    if mode == 'h':
        if save_OT:
            return adata_cm, tran
        return adata_cm

    else:
        if save_OT:
            return adata_concat, tran
        return adata_concat 

