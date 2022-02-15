#!/usr/bin/env 
"""
# Author: Kai Cao
# Modified from SCALEX
"""

import torch
import numpy as np

import os
import scanpy as sc
from anndata import AnnData
import scipy

from .modal.vae import VAE
from .modal.utils import EarlyStopping
from .metrics import batch_entropy_mixing_score, silhouette_score
from .logger import create_logger
from .data_loader import load_data

def Get_label_Prior(alpha, celltype1, celltype2):

    """
    Create a prior correspondence matrix according to cell labels
    
    Parameters
    ----------
    alpha
        the confidence of label, ranges from (1, inf). Higher alpha means better confidence.
    celltype1
        cell labels of dataset X
    celltype2
        cell labels of dataset Y

    Return
    ------
    Couple
        a prior correspondence matrix 
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
    Reweight labels if all cell types share the same weight 
    
    Parameters
    ----------
    celltype1
        cell labels 

    Return
    ------
    Weights
        a vector of weight 
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
        adata_cm = None,   
        mode='h',
        Prior = None,
        ref_id=None,    
        save_OT=False,
        rep_celltype='cell_type',
        use_specific = True,
        lambda_1=0.5,
        Lambda=0.5,
        gamma=1.0,
        batch_size=256, 
        lr=2e-4, 
        max_iteration=60000,
        seed=124, 
        gpu=0, 
        outdir='output/', 
        out='latent',
        input_id=0,
        pred_id=1,
        ignore_umap=False,
        verbose=False,
        assess=False,
        show=False,
        source_name='source',
        batch_key='domain_id',
        label_weight=None,
    ):

    """
    main function
    
    Parameters
    ----------
    adatas
        List of AnnData matrices for each dataset.
    adata_cm
        AnnData containing common genes.
    mode
        Choose from ['h', 'v', 'd']
        If 'h', integrate data with common genes
        If 'v', integrate data profiled from the same cells
        If 'd', inetrgate data without common genes
        Default: 'h'.
    Prior
        Prior correspondence matrix.
    ref_id
        Id of reference dataset.
    save_OT
        If True, output a global OT plan. Default: False.
    rep_celltype
        Names of cell-type annotation in AnnData. Default: 'cell_type'.
    use_specific
        If True, specific genes in each dataset will be considered. Default: True.
    lambda_1
        Balanced parameter for specific genes. Default: 0.5.
    Lambda: 
        Balanced parameter for KL divergence. Default: 0.5.
    gamma:
        Balanced parameter for OT. Default: 1.0.
    batch_size
        Number of samples per batch to load. Default: 256.
    lr
        Learning rate. Default: 2e-4.
    max_iteration
        Max iterations for training. Training one batch_size samples is one iteration. Default: 60000.
    seed
        Random seed for torch and numpy. Default: 124.
    gpu
        Index of GPU to use if GPU is available. Default: 0.
    outdir
        Output directory. Default: 'output/'.
    out
        Output of uniPort. Choose from ['latent', 'project', 'predict']. 
        If out=='latent', train the network and output cell embeddings. 
        If out=='project', project data into the latent space and output cell embeddings. 
        If out=='predict', project data into the latent space and output cell embeddings through a specified decoder. 
        Default: 'latent'.
    input_id
        Only used when mode=='d' and out=='predict' to choose a encoder to project data. Default: 0.
    pred_id
        Only used when out=='predict' to choose a decoder to predict data. Default: 1.
    ignore_umap
        If True, do not perform UMAP for visualization and leiden for clustering. Default: False.
    verbose
        Verbosity, True or False. Default: False.
    assess
        If True, calculate the entropy_batch_mixing score and silhouette score to evaluate integration results. Default: False.
    show
        If True, show the UMAP visualization of latent space. Default: False.
    source_name
        Name of source in AnnData. Default: source.
    batch_key
        Name of batch in AnnData. Default: domain_id.
    label_weight
        Prior-guided weighted vectors. Default: None.
    
    Return
    ------
    Weights
        a vector of weight 
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
    
    outdir = outdir+'/'
    os.makedirs(outdir+'/checkpoint', exist_ok=True)
    log = create_logger('', fh=outdir+'log.txt')

    if adatas is None:
        use_specific = False
        _, idx = np.unique(adata_cm.obs[source_name], return_index=True)
        batches = adata_cm.obs[source_name][np.sort(idx)]
        print(batches)
        flagged = []
        for batch in batches:
            flagged.append(adata_cm[adata_cm.obs[source_name]==batch].copy())

        adatas = flagged

    n_domain = len(adatas)
    if ref_id is None:
        ref_id = n_domain-1

    tran = {}
    num_cell = []
    num_gene = []

    for i, adata in enumerate(adatas):
        print('dataset {}:'.format(i))
        print(adata)

        num_cell.append(adata.X.shape[0])
        num_gene.append(adata.X.shape[1])
        print('reference dataset {}'.format(ref_id))

    if out == 'latent':

        if save_OT:           
            memory = 0
            for i in range(n_domain):
                if i != ref_id:
                    memory += num_cell[i]*num_cell[ref_id]*32/(8*1e9) * 4

            print('Warning! Saving Optimal Transport plan needs extra {:.2f} GB memory, please set save_OT=False if no enough memory!'.format(memory)) 
            for i in range(n_domain):
                if i != ref_id:
                    ns = num_cell[i]
                    nt = num_cell[ref_id]
                    tran_tmp = np.ones((ns, nt)) / (ns * nt)
                    tran[i] = tran_tmp.astype(np.float32)

                    print(tran[i].dtype)

                    print('Size of transport plan between datasets {} and {}:'.format(i, ref_id), np.shape(tran[i]))

        if adata_cm is not None:
            print(adata_cm)

        trainloader, testloader = load_data(
            adatas=adatas, 
            num_cell=num_cell,
            max_gene=max(num_gene), 
            adata_cm=adata_cm,
            use_specific=use_specific, 
            domain_name=batch_key,
            batch_size=batch_size, 
            mode=mode
        )

        early_stopping = EarlyStopping(patience=10, checkpoint_file=outdir+'/checkpoint/model.pt')
       
        dec = {}
        enc = [['fc', 1024, 1, 'relu'],['fc', 16, '', '']]
        if mode == 'd':     
            for i in range(n_domain):          
                dec[i] = [['fc', num_gene[i], 1, 'sigmoid']]

        elif mode == 'h':      
            num_gene.append(adata_cm.X.shape[1]) 
            dec[0] = [['fc', num_gene[n_domain], n_domain, 'sigmoid']]
            if use_specific:
                for i in range(1, n_domain+1):
                    dec[i] = [['fc', num_gene[i-1], 1, 'sigmoid']]

        else:
            for i in range(n_domain):
                dec[i] = [['fc', num_gene[i], 1, 'sigmoid']]

        model = VAE(enc, dec, ref_id=ref_id, n_domain=n_domain, mode=mode)

        log.info('model\n'+model.__repr__())

        model.fit(
            trainloader, 
            tran,
            num_cell,
            num_gene,
            label_weight=label_weight,
            rep_celltype=rep_celltype,
            Prior=Prior,
            save_OT=save_OT,
            use_specific=use_specific,
            lambda_1=lambda_1,
            Lambda=Lambda,
            gamma=gamma,
            lr=lr, 
            max_iteration=max_iteration, 
            device=device, 
            early_stopping=early_stopping, 
            verbose=verbose,
            mode=mode,
        )
        torch.save({'enc':enc, 'dec':dec, 'n_domain':n_domain, 'ref_id':ref_id, 'num_gene':num_gene}, outdir+'/checkpoint/config.pt')     

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


    if not ignore_umap: #and adata.shape[0]<1e6:
        log.info('Plot umap')
        sc.settings.figdir = outdir
        sc.set_figure_params(dpi=200, fontsize=10)
 
        if mode == 'h':
            sc.pp.neighbors(adata_cm, n_neighbors=30, use_rep=out)
            sc.tl.umap(adata_cm, min_dist=0.1)
            sc.tl.leiden(adata_cm)
            
            cols = [source_name, rep_celltype]
            color = [c for c in cols if c in adata_cm.obs]

            if len(color) > 0:
                sc.pl.umap(adata_cm, color=color, save='result2.pdf', title=['',''], legend_fontsize=10, s=2, show=show, \
                        wspace=0.4)
 
            if assess:
                if len(adata_cm.obs[batch_key].cat.categories) > 1:
                    entropy_score = batch_entropy_mixing_score(adata_cm.obsm['X_umap'], adata_cm.obs[batch_key])
                    log.info('batch_entropy_mixing_score: {:.3f}'.format(entropy_score))

                if rep_celltype in adata_cm.obs:
                    sil_score = silhouette_score(adata_cm.obsm['X_umap'], adata_cm.obs[rep_celltype].cat.codes)
                    log.info("silhouette_score: {:.3f}".format(sil_score))

        else:
            log.info('Plot umap')
            for i in range(n_domain-1):
                adata_concat = adatas[i].concatenate(adatas[i+1])
            sc.pp.neighbors(adata_concat, n_neighbors=30, use_rep=out)
            sc.tl.umap(adata_concat, min_dist=0.1)
            sc.pl.umap(adata_concat, color=[source_name, rep_celltype], save='specific.pdf', wspace=0.3, legend_fontsize=14, \
                title=['',''], s=5, show=show)

            entropy_score = batch_entropy_mixing_score(adata_concat.obsm['X_umap'], adata_concat.obs[source_name])
            log.info('batch_entropy_mixing_score: {:.3f}'.format(entropy_score))
            sil_score = silhouette_score(adata_concat.obsm['X_umap'], adata_concat.obs[rep_celltype].cat.codes)
            log.info("silhouette_score: {:.3f}".format(sil_score))
            

    if mode == 'h':
        if save_OT:
            return adata_cm, tran
        return adata_cm

    else:
        if save_OT:
            return adata_concat, tran
        return adata_concat 
    

def label_transfer(ref, query, rep='latent', label='celltype'):
    """
    From SCALEX
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
    
    knn = knn = KNeighborsClassifier().fit(X_train, y_train)
    y_test = knn.predict(X_test)
    
    return y_test

