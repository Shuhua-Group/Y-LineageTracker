import os
import sys


current_dir = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(current_dir)
sys.path.append(current_dir+'/..')
