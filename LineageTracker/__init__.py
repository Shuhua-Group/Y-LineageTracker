#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import sys

'''
information about version and module
'''

__VERSION__ = '1.3.0'
__MODULE__ = 'Y-LineageTracker'

current_dir = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(current_dir)
sys.path.append('..')
