.. module:: uniport
.. automodule:: uniport
   :noindex:

API and modules
===============

Function
--------
.. module:: uniport
.. currentmodule:: uniport

.. autosummary::
    :toctree: .
    
    Run
    Get_label_Prior
    label_reweight
    load_file
    filter_data
    batch_scale

DataLoader
-----------
.. module:: uniport.data_loader
.. currentmodule:: uniport

.. autosummary::
    :toctree: .

    data_loader.BatchSampler
    data_loader.BatchSampler_balance
    data_loader.SingleCellDataset
    data_loader.SingleCellDataset_vertical
    data_loader.load_data

Model
-----
.. module:: uniport.Model
.. currentmodule:: uniport

.. autosummary::
    :toctree: .

    model.vae.VAE
    model.layer.DSBatchNorm
    model.layer.Block
    model.layer.NN
    model.layer.Encoder
    model.layer.Decoder
    model.loss.kl_div
    model.loss.distance_matrix
    model.loss.distance_gmm
    model.loss.unbalanced_ot
    model.utils.onehot
    model.utils.EarlyStopping

Evaluation
----------
.. module:: uniport.metrics
.. currentmodule:: uniport

.. autosummary::
    :toctree: .

    metrics.batch_entropy_mixing_score
    metrics.silhouette
    metrics.label_transfer

Logger
------
.. module:: uniport.logger
.. currentmodule:: uniport

.. autosummary::
    :toctree: .
    
    logger.create_logger
   