#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @File   : configure.py
# @Author : Feng Xu @ crazyhsu9527@gmail.com
# @Created: 2024-09-13

import sys
import os
import re
import glob
import subprocess, argparse
import shutil

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--template_option_file', type=str, required=True)
    parser.add_argument('--conda_env_dir', type=str, required=True)
    # parser.add_argument('--install_dir', type=str, required=True)
    # parser.add_argument('--multicom3_db_dir', type=str, required=True)
    parser.add_argument('--tools_dir', type=str, required=True)
    parser.add_argument('--afdb_dir', type=str, required=True)
    parser.add_argument('--af2_dir', type=str, required=True)
    parser.add_argument('--out_option_file', default="option_file", type=str, required=False)
    args = parser.parse_args()

    # configure db_option file
    template_option_file = args.template_option_file
    newlines = []
    keywords_dict = {'YOUR_ENV': os.path.abspath(args.conda_env_dir.rstrip('/')),
                     'INSTALLDIR_TOOLS': args.tools_dir.rstrip('/'),
                     'AFDB_DIR': args.afdb_dir.rstrip('/'),
                     'AF_DIR': args.af2_dir.rstrip('/')}

    for line in open(template_option_file):
        newline = line
        for keyword in keywords_dict:
            newline = newline.replace(keyword, keywords_dict[keyword])
        newlines += [newline]

    with open(args.out_option_file, 'w') as fw:
        fw.writelines(''.join(newlines))

    print("\nConfiguration....Done")
