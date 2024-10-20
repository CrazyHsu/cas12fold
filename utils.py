#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @File   : utils.py
# @Author : Feng Xu @ crazyhsu9527@gmail.com
# @Created: 2024-09-12

"""Common utilities for data pipeline tools."""

import os, sys, argparse
import contextlib
import shutil
import tempfile
import time
import numpy as np
from typing import Optional
from absl import logging
from typing import Any, Mapping, MutableMapping, Optional, Sequence, Union
from alphafold.data import msa_identifiers
from alphafold.common import residue_constants


@contextlib.contextmanager
def tmpdir_manager(base_dir: Optional[str] = None):
  """Context manager that deletes a temporary directory on exit."""
  tmpdir = tempfile.mkdtemp(dir=base_dir)
  try:
    yield tmpdir
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


@contextlib.contextmanager
def timing(msg: str):
  logging.info('Started %s', msg)
  tic = time.time()
  yield
  toc = time.time()
  logging.info('Finished %s in %.3f seconds', msg, toc - tic)


def die(msg):
    print(msg)
    sys.exit(1)


def check_dirs(params, keys, isdir=True):
    errmsg = ''
    for key in keys:
        dirpath = params[key]
        # print(f"{key}:{params[key]}")
        if isdir and not os.path.isdir(dirpath):
            errmsg = errmsg + '{}({})\n'.format(key, dirpath)

        if not isdir and not os.path.exists(dirpath):
            errmsg = errmsg + '{}({})\n'.format(key, dirpath)

    if len(errmsg) > 0:
        errmsg = 'Directories or files are not exist:\n' + errmsg
        raise argparse.ArgumentTypeError(errmsg)


def check_contents(params, keys):
    errmsg = ''
    for key in keys:
        name = params[key]
        if len(name) == 0:
            errmsg = errmsg + '{}\n'.format(key)

    if len(errmsg) > 0:
        errmsg = 'These contents are emply:\n' + errmsg
        raise argparse.ArgumentTypeError(errmsg)


def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def check_dir(dirname):
    return is_dir(dirname)


def is_file(filename):
    """Checks if a file is an invalid file"""
    if not os.path.exists(filename):
        msg = "{0} doesn't exist".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def check_file(dirname):
    return is_file(dirname)


def makedir_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory = os.path.abspath(directory)
    return directory


def read_option_file(option_file):
    if not os.path.exists(option_file):
        die("Option file %s not exists." % option_file)
    params = {}
    for line in open(option_file):
        line = line.rstrip()
        if line.startswith('#'):
            continue
        tmp = line.split('=')
        if len(tmp) != 2:
            continue
        key = tmp[0].lstrip().rstrip()
        value = tmp[1].lstrip().rstrip()
        params[key] = value
    return params


def clean_dir(dir):
    if os.path.exists(dir):
        os.system(f'rm -rf {dir}')
    os.makedirs(dir)


def create_file(file):
    f = open(file, 'w')
    f.close()


# def tmpdir_manager(base_dir):
#     """Context manager that deletes a temporary directory on exit."""
#     tmpdir = tempfile.mkdtemp(dir=base_dir)
#     try:
#         yield tmpdir
#     finally:
#         shutil.rmtree(tmpdir, ignore_errors=True)
#
#
# def timing(msg: str):
#     logging.info('Started %s', msg)
#     tic = time.time()
#     yield
#     toc = time.time()
#     logging.info('Finished %s in %.3f seconds', msg, toc - tic)


def build_alignment_indices(sequence, start_index):
    indices_list = []
    counter = start_index
    for symbol in sequence:
        if symbol == '-':
            indices_list.append(-1)
        else:
            indices_list.append(counter)
            counter += 1
    return indices_list


def convert_taln_seq_to_a3m(query_non_gaps, aln):
    for is_query_res_non_gap, sequence_res in zip(query_non_gaps, aln):
        if is_query_res_non_gap:
            yield sequence_res


def combine_a3ms(infiles, outfile):
    descriptions = []
    seqs = []
    for infile in infiles:
        for line in open(infile):
            line = line.rstrip('\n')
            if line.startswith('>'):
                descriptions += [line]
            else:
                seqs += [line]

    with open(outfile, 'w') as fw:
        for (desc, seq) in zip(descriptions, seqs):
            fw.write(f"{desc}\n{seq}\n")

