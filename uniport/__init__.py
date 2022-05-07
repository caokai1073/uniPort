# from pkg_resources import get_distribution

# __version__ = get_distribution('scalex').version
__version__ = '1.0.1'
__author__ = 'Kai Cao'
__email__ = 'caokai@amss.ac.cn'

from .function import Run, Get_label_Prior, label_reweight, load_file, filter_data, batch_scale