# from pkg_resources import get_distribution

# __version__ = get_distribution('scalex').version
__version__ = '1.0'
__author__ = 'Kai Cao'
__email__ = 'caokai@amss.ac.cn'

from .function import Run, label_transfer, Get_label_Prior, label_reweight
from .data_process import preprocess