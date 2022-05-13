Release
=======

v1.1.0
------
 
- function **Get_label_Prior** changes: 
    + **Get_label_Prior()** -> **get_prior()**
    + set **alpha=2** as default
- function **Run** parameter changes: 
    + **labmda_recon** -> **lambda_recon** 
    + **Prior** -> **prior**
    + **max_iteration** -> **iteration**
    + add **use_rep** for **mode=d**
- add function TFIDF_LSI() for scATAC preprecess.
- remove AnnData return of **filter_data** and **batch_scale**.
    + e.g., change **adata=up.filter_data(adata)** to **up.filter_data(adata)**.

v1.0.5
------

Parameter fixes. 
