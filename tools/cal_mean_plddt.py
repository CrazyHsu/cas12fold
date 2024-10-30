#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @File   : cal_mean_plddt.py
# @Author : Feng Xu @ crazyhsu9527@gmail.com
# @Created: 2024-10-24

import biotite.structure.io as bsio
import sys

pdb_file = sys.argv[1]
struct = bsio.load_structure(pdb_file, extra_fields=["b_factor"])
print("{}\t{}".format(pdb_file, struct.b_factor.mean()))