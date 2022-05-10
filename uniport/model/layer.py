#!/usr/bin/env 
"""
# Author: Kai Cao
# Modified from SCALEX
"""

import math
import numpy as np

import torch
from torch import nn as nn
import torch.nn.functional as F
from torch.distributions import Normal
from torch.nn.parameter import Parameter
from torch.nn import init
from torch.autograd import Function


activation = {
    'relu':nn.ReLU(),
    'rrelu':nn.RReLU(),
    'sigmoid':nn.Sigmoid(),
    'leaky_relu':nn.LeakyReLU(),
    'tanh':nn.Tanh(),
    '':None
}


class DSBatchNorm(nn.Module):
    """
    Domain-specific Batch Normalization
    """
    def __init__(self, num_features, n_domain, eps=1e-5, momentum=0.1):
        """
        Parameters
        ----------
        num_features
            dimension of the features
        n_domain
            domain number
        """
        super().__init__()
        self.n_domain = n_domain
        self.num_features = num_features
        self.bns = nn.ModuleList([nn.BatchNorm1d(num_features, eps=eps, momentum=momentum) for i in range(n_domain)])
        
    def reset_running_stats(self):
        for bn in self.bns:
            bn.reset_running_stats()
            
    def reset_parameters(self):
        for bn in self.bns:
            bn.reset_parameters()
            
    def _check_input_dim(self, input):
        raise NotImplementedError
            
    def forward(self, x, y):
        out = torch.zeros(x.size(0), self.num_features, device=x.device) #, requires_grad=False)
        for i in range(self.n_domain):
            indices = np.where(y.cpu().numpy()==i)[0]

            if len(indices) > 1:
                out[indices] = self.bns[i](x[indices])
            elif len(indices) == 1:
                # out[indices] = x[indices]
                self.bns[i].training = False
                out[indices] = self.bns[i](x[indices])
                self.bns[i].training = True
        return out
        

class Block(nn.Module):
    """
    Basic block consist of:
        fc -> bn -> act -> dropout
    """
    def __init__(
            self,
            input_dim, 
            output_dim, 
            norm='', 
            act='', 
            dropout=0
        ):
        """
        Parameters
        ----------
        input_dim
            dimension of input
        output_dim
            dimension of output
        norm
            batch normalization, 
                * '' represent no batch normalization
                * 1 represent regular batch normalization
                * int>1 represent domain-specific batch normalization of n domain
        act
            activation function,
                * relu -> nn.ReLU
                * rrelu -> nn.RReLU
                * sigmoid -> nn.Sigmoid()
                * leaky_relu -> nn.LeakyReLU()
                * tanh -> nn.Tanh()
                * '' -> None
        dropout
            dropout rate
        """
        super().__init__()
        self.fc = nn.Linear(input_dim, output_dim)
        
        if type(norm) == int:
            if norm==1: # TO DO
                self.norm = nn.BatchNorm1d(output_dim)
            else:
                self.norm = DSBatchNorm(output_dim, norm)
        else:
            self.norm = None
            
        self.act = activation[act]
            
        if dropout >0:
            self.dropout = nn.Dropout(dropout)
        else:
            self.dropout = None
            
    def forward(self, x, y=None):
        h = self.fc(x)
        if self.norm:
            if len(x) == 1:
                pass
            elif self.norm.__class__.__name__ == 'DSBatchNorm':
                h = self.norm(h, y)
            else:
                h = self.norm(h)
        if self.act:
            h = self.act(h)
        if self.dropout:
            h = self.dropout(h)
        return h
        
    
class NN(nn.Module):
    """
    Neural network consist of multi Blocks
    """
    def __init__(self, input_dim, cfg):
        """
        Parameters
        ----------
        input_dim
            input dimension
        cfg
            model structure configuration, 'fc' -> fully connected layer
            
        Example
        -------
        >>> latent_dim = 10
        >>> dec_cfg = [['fc', x_dim, n_domain, 'sigmoid']]
        >>> decoder = NN(latent_dim, dec_cfg)
        """
        super().__init__()
        net = []
        for i, layer in enumerate(cfg):
            if i==0:
                d_in = input_dim
            if layer[0] == 'fc':
                net.append(Block(d_in, *layer[1:]))
            d_in = layer[1]
        self.net = nn.ModuleList(net)
    
    def forward(self, x, y=None):
        for layer in self.net:
            x = layer(x, y)
        return x
    
        
class Encoder(nn.Module):
    """
    VAE Encoder
    """
    def __init__(self, input_dim, cfg, mode):
        """
        Parameters
        ----------
        input_dim
            input dimension
        cfg
            encoder configuration, e.g. enc_cfg = [['fc', 1024, 1, 'relu'],['fc', 10, '', '']]
        mode
            training mode. ['h', 'd', 'v']
        """
        super().__init__()

        enc = []
        mu_enc = []
        var_enc = []

        h_dim = cfg[-2][1]
        
        if mode == 'd':
            for i in input_dim.keys():             
                enc.append(NN(input_dim[i], cfg[:-1]))
                mu_enc.append(NN(h_dim, cfg[-1:]))
                var_enc.append(NN(h_dim, cfg[-1:]))

        else:
            enc.append(NN(input_dim[0], cfg[:-1]))
            mu_enc.append(NN(h_dim, cfg[-1:]))
            var_enc.append(NN(h_dim, cfg[-1:]))

        self.enc = nn.ModuleList(enc)
        self.mu_enc = nn.ModuleList(mu_enc)
        self.var_enc = nn.ModuleList(var_enc)
        
    def reparameterize(self, mu, var):
        return Normal(mu, var.sqrt()).rsample()

    def forward(self, x, domain, y=None):
        """
        """
        q = self.enc[domain](x, y)
        mu = self.mu_enc[domain](q, y)
        var = torch.exp(self.var_enc[domain](q, y))
        z = self.reparameterize(mu, var)
        return z, mu, var


class Decoder(nn.Module):
    """
    VAE Decoder
    """
    def __init__(self, z_dim, cfg):
        """
        Parameters
        ----------
        z_dim
            latent dimension
        cfg
            decoder configuration, e.g. dec_cfg = [['fc', adatas[i].obsm[obsm[i]].shape[1], 1, 'sigmoid']]
        """
        super().__init__()

        dec = []       
        for i in cfg.keys():
            dec.append(NN(z_dim, cfg[i]))

        self.dec = nn.ModuleList(dec)


    def forward(self, z, domain, y=None):
        """
        """
        reconx_x = self.dec[domain](z, y)

        return reconx_x
    