# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 09:13:49 2024

@author: jacob
"""

import os

def correct_path_change_dir(uncorrected_path):
    corrected_path = os.path.normpath(uncorrected_path)
    os.chdir(corrected_path)
