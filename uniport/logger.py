#!/usr/bin/env 
"""
# Author: Kai Cao
# Modified from SCALEX
"""

import logging

def create_logger(name='', ch=True, fh='', levelname=logging.INFO):
    logger = logging.getLogger(name)
    logger.setLevel(levelname)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # handler 
    if ch:
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    if fh:
        fh = logging.FileHandler(fh, mode='w')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    return logger

