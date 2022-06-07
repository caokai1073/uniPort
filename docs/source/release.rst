Release
=======
v1.1.1
------
- add model_log parameter to *Run()* function. 
    + If model_log=True, show structures of encoder and decoders. Default: False.
- fix bugs

v1.1.0
------
 
- *Get_label_Prior()* function changes: 
    + **Get_label_Prior()** -> **get_prior()**
    + set **alpha=2** as default
- *Run()* function parameter changes: 
    + **labmda_recon** -> **lambda_recon** 
    + **Prior** -> **prior**
    + **max_iteration** -> **iteration**
    + add **use_rep** for **mode=d**
- add *TFIDF_LSI()* function for scATAC preprecess.
- remove AnnData returns in *filter_data()* and *batch_scale()* functions.
    + e.g., change **adata=up.filter_data(adata)** to **up.filter_data(adata)**.

v1.0.5
------

Parameter fixes. 
