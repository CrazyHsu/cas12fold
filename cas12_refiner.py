#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @File   : cas12_refiner.py
# @Author : Feng Xu @ crazyhsu9527@gmail.com
# @Created: 2024-09-12
import os

from pipeline import *
from pipeline_custom import *
from pathlib import Path

from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('pdb_name', None, 'Name of the pdb. Default is the prefix of the fasta file')
flags.DEFINE_string('fasta_path', None, 'Path to cas12 fasta')
flags.DEFINE_string('pdb_path', None, 'Path to af2 predicted pdb used to refine')
flags.DEFINE_string('pkl_path', None, 'Path to af2 predicted pkl used to refine')
flags.DEFINE_string('msas_dir', None, 'Path to directory of af2 predicted msas')
flags.DEFINE_string('output_dir', default="cas12_refine", help='Output directory')
flags.DEFINE_string('gpu_device', default="0", help='Output directory')
flags.DEFINE_boolean('run_refinement', True, 'Whether to run model refinement')
FLAGS = flags.FLAGS

def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    os.environ['TF_FORCE_UNIFIED_MEMORY'] = '1'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '4.0'
    os.environ['CUDA_VISIBLE_DEVICES'] = FLAGS.gpu_device

    FLAGS.fasta_path = os.path.abspath(FLAGS.fasta_path)
    FLAGS.output_dir = os.path.abspath(FLAGS.output_dir)

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(FLAGS.fasta_path)

    targetname = pathlib.Path(FLAGS.fasta_path).stem
    sequence = open(FLAGS.fasta_path).readlines()[1].rstrip('\n').strip()

    outdir = FLAGS.output_dir
    outdir_final = FLAGS.output_dir + "_final"
    tmp_input = os.path.join(outdir, "tmp_input")

    makedir_if_not_exists(outdir)
    makedir_if_not_exists(tmp_input)

    refine_inputs = []

    pdb_name = FLAGS.pdb_name if not FLAGS.pdb_name else Path(FLAGS.fasta_path).stem
    renamed_fasta_file = os.path.join(tmp_input, f"{pdb_name}.fa")
    renamed_pdb_file = os.path.join(tmp_input, f"{pdb_name}.pdb")
    renamed_pkl_file = os.path.join(tmp_input, f"{pdb_name}.pkl")
    renamed_msa_file = os.path.join(tmp_input, f"{pdb_name}.a3m")

    os.system(f"cp {FLAGS.fasta_path} {renamed_fasta_file}")
    os.system(f"cp {FLAGS.pdb_path} {renamed_pdb_file}")
    os.system(f"cp {FLAGS.pkl_path} {renamed_pkl_file}")
    merge_msa(FLAGS.msas_dir, tmp_input, f"{pdb_name}.a3m")
    refine_input = refinement_input(fasta_path=renamed_fasta_file,
                                    pdb_path=renamed_pdb_file,
                                    pkl_path=renamed_pkl_file,
                                    msa_path=renamed_msa_file)

    refine_inputs += [refine_input]

    # for i in range(5):
    #     pdb_name = ref_ranking_avg.loc[i, 'model']
    #     refine_input = refinement_input(fasta_path=FLAGS.fasta_path,
    #                                     pdb_path=FLAGS.pdb_path,
    #                                     pkl_path=os.path.join(N4_outdir, 'pkl',
    #                                                           pdb_name.replace('.pdb', '.pkl')),
    #                                     msa_path=os.path.join(N4_outdir, 'msa', pdb_name.replace('.pdb', '.a3m')))
    #     refine_inputs += [refine_input]

    run_cas12_refinement_pipeline(params=params, refinement_inputs=refine_inputs, outdir=outdir, finaldir=outdir_final,
                                  prefix="refine")



if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'pdb_path',
        'pkl_path',
        'msas_dir'
    ])
    app.run(main)

# flags.DEFINE_string('option_file', None, 'option file')
# flags.DEFINE_string('pdb_name', None, 'Name of the pdb. Default is the prefix of the fasta file')
# flags.DEFINE_string('fasta_path', None, 'Path to cas12 fasta')
# flags.DEFINE_string('pdb_path', None, 'Path to af2 predicted pdb used to refine')
# flags.DEFINE_string('pkl_path', None, 'Path to af2 predicted pkl used to refine')
# flags.DEFINE_string('msas', None, 'Path to directory of af2 predicted msas')
# flags.DEFINE_string('output_dir', default="cas12_refine", help='Output directory')
# flags.DEFINE_boolean('run_refinement', True, 'Whether to run model refinement')
# python cas12_refiner.py \
#     --option_file=db_option \
#     --pdb_path=/path/to/pdb \
#     --feat_pkl=/path/to/feature_pkl \
#     --msas_dir=/path/to/msas
#     --output_dir=$OUTDIR