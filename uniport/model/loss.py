#!/usr/bin/env 
"""
# Author: Kai Cao
# Modified from RAE
"""

import torch
from torch.distributions import Normal, kl_divergence

__all__ = ['kl_div', 'distance_matrix', 'distance_gmm', 'unbalanced_ot']

def kl_div(mu, var, weight=None):
    loss = kl_divergence(Normal(mu, var.sqrt()), Normal(torch.zeros_like(mu),torch.ones_like(var))).sum(dim=1)
    
    # if weight is not None:
    #     loss = loss * weight.squeeze(dim=1)
    return loss.mean()
 
# def balanced_binary_cross_entropy(recon_x, x):

#     return -torch.sum(x * torch.log(recon_x + 1e-8) + (1 - x) * torch.log(1 - recon_x + 1e-8), dim=-1)

def distance_matrix(pts_src: torch.Tensor, pts_dst: torch.Tensor, p: int = 2):
    """
    Returns the matrix of ||x_i-y_j||_p^p.

    Parameters
    ----------
    pts_src
        [R, D] matrix
    pts_dst
        C, D] matrix
    p
        p-norm
    
    Return
    ------
    [R, C] matrix
        distance matrix
    """
    diff = pts_src.unsqueeze(1) - pts_dst.unsqueeze(0)
    # |d|**2 is bit-identical to d*d but avoids an abs pass and the pow kernel.
    distance = torch.sum(diff * diff if p == 2 else torch.abs(diff) ** p, 2)
    return distance

def distance_gmm(mu_src: torch.Tensor, mu_dst: torch.Tensor, var_src: torch.Tensor, var_dst: torch.Tensor):
    """
    Calculate a Wasserstein distance matrix between the gmm distributions with diagonal variances

    Parameters
    ----------
    mu_src
        [R, D] matrix, the means of R Gaussian distributions
    mu_dst
        [C, D] matrix, the means of C Gaussian distributions
    logvar_src
        [R, D] matrix, the log(variance) of R Gaussian distributions
    logvar_dst
        [C, D] matrix, the log(variance) of C Gaussian distributions
    
    Return
    ------
    [R, C] matrix 
        distance matrix
    """
    std_src = var_src.sqrt()
    std_dst = var_dst.sqrt()
    distance_mean = distance_matrix(mu_src, mu_dst, p=2)
    distance_var = distance_matrix(std_src, std_dst, p=2)

    # distance_var = torch.sum(sum_matrix(std_src, std_dst) - 2 * (prod_matrix(std_src, std_dst) ** 0.5), 2)
    
    return distance_mean + distance_var + 1e-6

def unbalanced_ot(tran, mu1, var1, mu2, var2, reg=0.1, reg_m=1.0, Couple=None, device='cpu', \
    idx_q=None, idx_r=None, query_weight=None, ref_weight=None):
    '''
    Calculate a unbalanced optimal transport matrix between mini batches.

    Parameters
    ----------
    tran
        transport matrix between the two batches sampling from the global OT matrix. 
    mu1
        mean vector of batch 1 from the encoder
    var1
        standard deviation vector of batch 1 from the encoder
    mu2
        mean vector of batch 2 from the encoder
    var2
        standard deviation vector of batch 2 from the encoder
    reg:
        Entropy regularization parameter in OT. Default: 0.1
    reg_m:
        Unbalanced OT parameter. Larger values means more balanced OT. Default: 1.0
    Couple
        prior information about weights between cell correspondence. Default: None
    device
        training device
    idx_q
        domain_id of query batch
    idx_r
        domain_id of reference batch
    query_weight
        reweighted vectors of query batch
    ref_weight
        reweighted vectors of reference batch

    Returns
    -------
    float
        minibatch unbalanced optimal transport loss
    matrix
        minibatch unbalanced optimal transport matrix
    '''

    ns = mu1.size(0)
    nt = mu2.size(0)

    cost_pp = distance_gmm(mu1, mu2, var1, var2)

    if query_weight is None:
        p_s = torch.ones(ns, 1, device=device) / ns
    else:
        query_batch_weight = query_weight[idx_q]
        p_s = (query_batch_weight/torch.sum(query_batch_weight)).to(device)

    if ref_weight is None:
        p_t = torch.ones(nt, 1, device=device) / nt
    else:
        ref_batch_weight = ref_weight[idx_r]
        p_t = (ref_batch_weight/torch.sum(ref_batch_weight)).to(device)

    if tran is None:
        tran = torch.ones(ns, nt, device=device) / (ns * nt)

    dual = torch.ones(ns, 1, device=device) / ns
    f = reg_m / (reg_m + reg)

    # The Sinkhorn iterations below only ever reach the caller through
    # `tran.detach()`, so no autograd graph has to be recorded for them; the
    # gradient path of the returned loss runs through `cost_pp` alone.
    with torch.no_grad():
        cost = cost_pp.detach()
        if Couple is not None:
            cost = cost * Couple

        # `cost` is loop-invariant, so the exponentiated kernel is too: hoist it
        # out instead of recomputing exp() and max() on every outer iteration.
        kernel_base = torch.exp(-cost / (reg * torch.max(torch.abs(cost))))

        for m in range(10):
            kernel = kernel_base * tran
            kernel_t = kernel.t()
            b = p_t / (kernel_t @ dual)
            for i in range(10):
                dual = (p_s / (kernel @ b)) ** f
                b = (p_t / (kernel_t @ dual)) ** f
            tran = (dual @ b.t()) * kernel

        if torch.isnan(tran).any():
            tran = torch.ones(ns, nt, device=device) / (ns * nt)

    # pho = tran.mean()
    # h_func = 1 - 0.5 * ( 1 + torch.sign(pho - tran) )
    # hat_tran = tran * h_func
    # d_fgw1 = (cost_pp * hat_tran.detach().data).sum()
    # d_fgw2 = ((tran.detach().data - hat_tran.detach().data) * torch.log(1 + torch.exp(-cost_pp))).sum()
    # d_fgw = d_fgw1 + d_fgw2

    d_fgw = (cost_pp * tran).sum()

    return d_fgw, tran
















