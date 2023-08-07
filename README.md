<img src="docs/source/_static/uniPort.jpg" width="40%" height="40%">

[![license-badge](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/uniport.svg)](https://badge.fury.io/py/POT)
[![docs-badge](https://readthedocs.org/projects/uniport/badge/?version=latest)](https://uniport.readthedocs.io/en/latest/?badge=latest)
[![Downloads](https://pepy.tech/badge/uniport)](https://pepy.tech/project/uniport)
[![Downloads](https://pepy.tech/badge/uniport/month)](https://pepy.tech/project/uniport)

The original paper: 
[A unified single-cell data integration framework with optimal transport](https://www.nature.com/articles/s41467-022-35094-8)

Website and documentation: [https://uniport.readthedocs.io](https://uniport.readthedocs.io)

Source Code (MIT): [https://github.com/caokai1073/uniport](https://github.com/caokai1073/uniport)

Author's Homepage: [www.caokai.site](https://www.caokai.site)

![Overview](docs/source/_static/net.png)

## Installation

The **uniport** package can be installed via pip3:

```sh
pip3 install uniport
```

## Tutorials

### Please checkout the documentations and tutorials for more information at **[uniport.readthedocs.io](https://uniport.readthedocs.io)**.

### Main function: **uniport.Run**()


**Key parameters includes:**

+ **adatas**: List of AnnData matrices for each dataset.
+ **adata_cm**: AnnData matrix containing common genes from different datasets.
+ **mode**: Choose from ['h', 'v', 'd'] If 'mode=h', integrate data with common genes (Horizontal integration). If 'mode=v', integrate data profiled from the same cells (Vertical integration). If 'mode=d', inetrgate data without common genes (Diagonal integration). Default: 'h'.
+ **lambda_s**: balanced parameter for common and specific genes. Default: 0.5
+ **lambda_recon**: balanced parameter for reconstruct term. Default: 1.0
+ **lambda_kl**: balanced parameter for KL divergence. Default: 0.5
+ **lambda_ot**: balanced parameter for OT. Default: 1.0
+ **iteration**: max iterations for training. Training one batch_size samples is one iteration. Default: 30000
+ **ref_id**: id of reference dataset. Default: The domain_id of last dataset
+ **save_OT**: if True, output a global OT plan. Need more memory. Default: False
+ **out**: output of uniPort. Choose from ['latent', 'project', 'predict']. If out=='latent', train the network and output cell embeddings. If out=='project', project data into the latent space and output cell embeddings. If out=='predict', project data into the latent space and output cell embeddings through a specified decoder. Default: 'latent'

## Data
+ [**Google Drive**](https://drive.google.com/drive/folders/19WAzlfBAtI_EUIv0rgkZW1zwzUopsHAv?usp=share_link)
+ [**Baidu Drive**](https://pan.baidu.com/s/1CikO0K-ZWcKlxAgK7fkMcQ) Code: 1122

## Example
```Python
import uniport as up
import scanpy as sc

# HVG: highly variable genes
adata1 = sc.read_h5ad('adata1.h5ad') # preprocessed data with data1 specific HVG
adata2 = sc.read_h5ad('adata2.h5ad') # preprocessed data with data2 specific HVG, as reference data
adata_cm = sc.read_h5ad('adata_cm.h5ad') # preprocesssed data with common HVG

# integration with both common and dataset-specific genes
# latent representation are stored in adata.obs['latent']
adata = up.Run(adatas=[adata1, adata2], adata_cm=adata_cm)
# save global optimal transport matrix: adata, OT = up.Run(adatas=[adata1, adata2], adata_cm=adata_cm, save_OT=True)
# integration with only common genes: adata = up.Run(adata_cm=adata_cm)

```

## Citation
	@Article{Cao2022,
	author={Cao, Kai and Gong, Qiyu and Hong, Yiguang and Wan, Lin},
	title={A unified computational framework for single-cell data integration with optimal transport},
	journal={Nature Communications},
	year={2022},
	month={Dec},
	day={01},
	volume={13},
	number={1},
	pages={7419},
	issn={2041-1723},
	doi={10.1038/s41467-022-35094-8}}

Contact via caokai1073@gmail.com


