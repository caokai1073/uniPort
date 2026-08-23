"""uniPort: a unified single-cell data integration framework with optimal transport."""

__version__ = '1.2.0'
__author__ = 'Kai Cao'
__email__ = 'caokai@amss.ac.cn'

from .function import Run, get_prior, label_reweight, load_file, filter_data, batch_scale, TFIDF_LSI
from .metrics import batch_entropy_mixing_score, silhouette, label_transfer

__all__ = [
    'Run',
    'get_prior',
    'label_reweight',
    'load_file',
    'filter_data',
    'batch_scale',
    'TFIDF_LSI',
    'batch_entropy_mixing_score',
    'silhouette',
    'label_transfer',
]
