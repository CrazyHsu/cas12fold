#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @File   : run_cas12fold_refiner_controller.py
# @Author : Feng Xu @ crazyhsu9527@gmail.com
# @Created: 2024-10-16
import os

# from pipeline import *
# from pipeline_custom import *
import utils
from cas12fold_refiner import cas12fold_refiner_pipeline
# from alphafold.data import parsers
from pathlib import Path

from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('pdb_name', None, 'Name of the pdb. Default is the prefix of the fasta file')
flags.DEFINE_string('fasta_path', None, 'Path to cas12 fasta')
flags.DEFINE_string('pdb_path', None, 'Path to af2 predicted pdb used to refine')
flags.DEFINE_string('pkl_path', None, 'Path to af2 predicted pkl used to refine')
flags.DEFINE_string('msas_dir', None, 'Path to directory of af2 predicted msas')
flags.DEFINE_string('output_dir', default="cas12fold_refine", help='Output directory')
flags.DEFINE_string('gpu_device', default="0", help='The GPU devices to predict')
flags.DEFINE_boolean('result_overwrite', default=False, help='Whether to overwrite the output. Default: False')
flags.DEFINE_boolean('template_select_foldseek_global', default=False,
                     help="Whether to select foldseek global results in structure_templates.csv.")
FLAGS = flags.FLAGS


class RefinementInput:
    def __init__(self, fasta_path, pdb_path, pkl_path, msa_path):
        self.fasta_path = fasta_path
        self.pdb_path = pdb_path
        self.pkl_path = pkl_path
        self.msa_path = msa_path


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    os.environ['TF_FORCE_UNIFIED_MEMORY'] = '1'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '4.0'
    os.environ['CUDA_VISIBLE_DEVICES'] = FLAGS.gpu_device

    FLAGS.fasta_path = os.path.abspath(FLAGS.fasta_path)
    FLAGS.output_dir = os.path.abspath(FLAGS.output_dir)

    utils.check_file(FLAGS.option_file)

    params = utils.read_option_file(FLAGS.option_file)

    # makedir_if_not_exists(FLAGS.output_dir)

    # check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    utils.check_file(FLAGS.fasta_path)
    utils.check_file(FLAGS.pdb_path)
    utils.check_file(FLAGS.pkl_path)
    utils.check_file(FLAGS.msas_dir)

    # targetname = pathlib.Path(FLAGS.fasta_path).stem
    # sequence = open(FLAGS.fasta_path).readlines()[1].rstrip('\n').strip()

    outdir = FLAGS.output_dir
    outdir_final = outdir + "_final"
    tmp_input = os.path.join(outdir, "tmp_input")

    utils.makedir_if_not_exists(outdir)
    utils.makedir_if_not_exists(tmp_input)

    # refine_inputs = []

    pdb_name = FLAGS.pdb_name if not FLAGS.pdb_name else Path(FLAGS.fasta_path).stem
    renamed_fasta_file = os.path.join(tmp_input, f"{pdb_name}.fa")
    renamed_pdb_file = os.path.join(tmp_input, f"{pdb_name}.pdb")
    renamed_pkl_file = os.path.join(tmp_input, f"{pdb_name}.pkl")
    renamed_msa_file = os.path.join(tmp_input, f"{pdb_name}.a3m")

    os.system(f"cp {FLAGS.fasta_path} {renamed_fasta_file}")
    os.system(f"cp {FLAGS.pdb_path} {renamed_pdb_file}")
    os.system(f"cp {FLAGS.pkl_path} {renamed_pkl_file}")
    cas12fold_refiner_pipeline.merge_msa(FLAGS.msas_dir, tmp_input, f"{pdb_name}.a3m")
    refine_input = RefinementInput(fasta_path=renamed_fasta_file,
                                   pdb_path=renamed_pdb_file,
                                   pkl_path=renamed_pkl_file,
                                   msa_path=renamed_msa_file)

    # refine_inputs += [refine_input]
    cas12fold_refiner_pipeline.cas12fold_refiner(params=params, refinement_input=refine_input, outdir=outdir,
                                                 finaldir=outdir_final, prefix="refine",
                                                 result_overwrite=FLAGS.result_overwrite, gpu_device=FLAGS.gpu_device,
                                                 template_select_foldseek_global=FLAGS.template_select_foldseek_global)


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'pdb_path',
        'pkl_path',
        'msas_dir'
    ])
    app.run(main)
