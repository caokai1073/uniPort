Release
=======

v1.1.0
------
 
- function name changes: `Get_label_Prior()` -> `get_prior()`, and set `alpha=2` as default.
- parameter name changes:  change `labmda_recon` -> `lambda_recon`, `Prior` -> `prior`, `max_iteration` -> `iteration` and add use_rep for `mode=d` in up.Run().
- add function TFIDF_LSI() for scATAC preprecess.
- remove AnnData return of `filter_data` and `batch_scale`, e.g., change `adata=up.filter_data(adata)` to `up.filter_data(adata)`.

v1.0.5
------

Parameter fixes. 
