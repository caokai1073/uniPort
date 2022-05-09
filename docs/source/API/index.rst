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

    Model.vae.VAE
    Model.layer.DSBatchNorm
    Model.layer.Block
    Model.layer.NN
    Model.layer.Encoder
    Model.layer.Decoder
    Model.loss.kl_div
    Model.loss.distance_matrix
    Model.loss.distance_gmm
    Model.utils.onehot
    Model.utils.EarlyStopping

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
   