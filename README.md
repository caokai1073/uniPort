[![license-badge](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/uniport.svg)](https://badge.fury.io/py/POT)
[![docs-badge](https://readthedocs.org/projects/uniport/badge/?version=latest)](https://uniport.readthedocs.io/en/latest/?badge=latest)
[![Downloads](https://pepy.tech/badge/uniport)](https://pepy.tech/project/uniport)

The original paper: 
[a unified single-cell data integration framework with optimal transport](https://www.biorxiv.org/content/10.1101/2022.02.14.480323v1)

![Overview](docs/source/_static/net.png)

Website and documentation: [https://uniport.readthedocs.io](https://uniport.readthedocs.io)

Source Code (MIT): [https://github.com/caokai1073/uniport](https://github.com/caokai1073/uniport)

Author's Homepage: [www.caokai.site](https://www.caokai.site)

## Installation

The `uniport` package can be installed via pip:

```sh
pip3 install uniport
```

## Usage


### Please checkout the documentations and tutorials at
**[uniport.readthedocs.io](https://uniport.readthedocs.io)**.

```Python
# integration function
adata = up.Run(adatas=None, adata_cm=None, mode='h', lambda_s=0.5, labmda_recon=1.0, lambda_kl=0.5, lambda_ot=1.0, reg=0.1, reg_m=1.0, batch_size=256, lr=2e-4, max_iteration=30000, seed=124, gpu=0, Prior=None, label_weight=None, ref_id=None, save_OT=False, use_specific=True, loss_type='BCE', outdir='output/', out='latent', input_id=0, pred_id=1, source_name='source', rep_celltype='cell_type', batch_key='domain_id', enc=None, dec=None, umap=False, verbose=False, assess=False, show=False)
```

<font color='Dodgerblue'>**Key parameters includes:**</font>

+ **adatas**: List of AnnData matrices for each dataset.
+ **adata_cm**: AnnData matrix containing common genes from different datasets.
+ **mode**: Choose from ['h', 'v', 'd'] If 'mode=h', integrate data with common genes (Horizontal integration). If 'mode=v', integrate data profiled from the same cells (Vertical integration). If 'mode=d', inetrgate data without common genes (Diagonal integration). Default: 'h'.
+ **lambda_s**: balanced parameter for common and specific genes. Default: 0.5
+ **lambda_recon**: balanced parameter for reconstruct term. Default: 1.0
+ **lambda_kl**: balanced parameter for KL divergence. Default: 0.5
+ **lambda_ot**: balanced parameter for OT. Default: 1.0
+ **max_iteration**: max iterations for training. Training one batch_size samples is one iteration. Default: 30000
+ **ref_id**: id of reference dataset. Default: The domain_id of last dataset
+ **save_OT**: if True, output a global OT plan. Need more memory. Default: False
+ **out**: output of uniPort. Choose from ['latent', 'project', 'predict']. If out=='latent', train the network and output cell embeddings. If out=='project', project data into the latent space and output cell embeddings. If out=='predict', project data into the latent space and output cell embeddings through a specified decoder. Default: 'latent'



<font color='Dodgerblue'>**Other parameters include:**</font>

+ **label_weight**: prior-guided weighted vectors. Default: None
+ **reg**: entropy regularization parameter in OT. Default: 0.1
+ **reg_m**: unbalanced OT parameter. Larger values means more balanced OT. Default: 1.0
+ **batch_size**: number of samples per batch to load. Default: 256
+ **lr**: learning rate. Default: 2e-4
+ **enc**: structure of encoder. For example: enc=[['fc', '1024', 1, 'relu'], ['fc', 16, '', '']] means that the encoder contains two layers. The first layer is fully connected with 1024 neurons, a [DSBN](https://openaccess.thecvf.com/content_CVPR_2019/papers/Chang_Domain-Specific_Batch_Normalization_for_Unsupervised_Domain_Adaptation_CVPR_2019_paper.pdf) and activative function `relu`. The second layer is fully connected with 16 neurons without DSBN or activative function.
+ **gpu**: index of GPU to use if GPU is available. Default: 0
+ **Prior**: prior correspondence matrix. Default: None
+ **loss_type**: type of loss function. Choose from ['BCE', 'MSE', 'L1']. Default: 'BCE'
+ **outdir**: output directory. Default: 'output/'
+ **input_id**: only used when mode=='d' and out=='predict' to choose a encoder to project data. Default: 0
+ **pred_id**: only used when out=='predict' to choose a decoder to predict data. Default: 1
+ **seed**: random seed for torch and numpy. Default: 124
+ **batch_key**: name of batch in AnnData. Default: domain_id
+ **source_name**: name of source in AnnData. Default: source
+ **rep_celltype**: names of cell-type annotation in AnnData. Default: 'cell_type'
+ **umap**: if True, perform UMAP for visualization. Default: False
+ **assess**: if True, calculate the entropy_batch_mixing score and silhouette score to evaluate integration results. Default: False
+ **show**: if True, show the UMAP visualization of latent space. Default: False

## Example
```Python
import uniport as up
import scanpy as sc

# HVG: highly variable genes
adata1 = sc.read_h5ad('adata1.h5ad') # preprocessed data with data1 specific HVG
adata2 = sc.read_h5ad('adata2.h5ad') # preprocessed data with data2 specific HVG, as reference data
adata_cm = sc.read_h5ad('adata_cm.h5ad') # preprocesssed data with common HVG

# integration with both common and dataset-specific genes
adata = up.Run(adatas=[adata1, adata2], adata_cm=adata_cm)
# save global optimal transport matrix
adata, OT = up.Run(adatas=[adata1, adata2], adata_cm=adata_cm, save_OT=True)
# integration with only common genes
adata = up.Run(adata_cm=adata_cm)

# integration without common genes
adata = up.Run(adatas=[adata1, adata2], mode='d')

# integration with paired datasets
adata = up.Run(adatas=[adata1, adata2], mode='v')
```

## Citation
```
@article {Cao2022.02.14.480323,
	author = {Cao, Kai and Gong, Qiyu and Hong, Yiguang and Wan, Lin},
	title = {uniPort: a unified computational framework for single-cell data integration with optimal transport},
	year = {2022},
	doi = {10.1101/2022.02.14.480323},
	publisher = {Cold Spring Harbor Laboratory},
	journal = {bioRxiv}
}
```
