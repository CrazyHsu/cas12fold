#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @File   : run_cas12fold_controller.py
# @Author : Feng Xu @ crazyhsu9527@gmail.com
# @Created: 2024-09-18

import enum
from utils import *
# from pipeline import *
# from pipeline_custom import *
# from pathlib import Path

from absl import flags
from absl import app

@enum.unique
class Cas12folddbMergeMode(enum.Enum):
    single = 0
    merge = 1
    both = 2


flags.DEFINE_string('option_file', None, 'Path to option file')
flags.DEFINE_string('fasta_path', None, 'Path to fasta file')
flags.DEFINE_string('model_preset', default="monomer", help='Choose preset model configuration - the monomer model, '
                                                            'the monomer model with extra ensembling, monomer model with pTM head, '
                                                            'or multimer model. Default: monomer')
flags.DEFINE_string('max_template_date', default="2099-01-01", help='Maximum template release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). '
                                                                    'Important if folding historical test sets. Default: 2099-01-01')
flags.DEFINE_string('output_dir', default="cas12fold", help='Output directory')
flags.DEFINE_string('gpu_device', default="0", help='The GPU devices to predict')
flags.DEFINE_boolean('run_af2_wt', default=False, help='Whether to run raw alphafold2')
flags.DEFINE_boolean('af2_use_precomputed_msas', default=False, help='Whether to read existing MSAs')
flags.DEFINE_boolean('run_cas12fold', default=False, help='Whether to run cas12fold')
flags.DEFINE_boolean('use_cas12fold_template_db', default=False, help='Whether to use cas12fold template database')
flags.DEFINE_enum_class('cas12folddb_merge_mode', default=Cas12folddbMergeMode.both, enum_class=Cas12folddbMergeMode, help="Modes to merge cas12folddb")
flags.DEFINE_float('max_subsequence_ratio', default=1, help="The max subsequence ratio used to filter pdb hits")
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

    # makedir_if_not_exists(FLAGS.output_dir)

    # check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(FLAGS.fasta_path)

    # targetname = pathlib.Path(FLAGS.fasta_path).stem
    # sequence = open(FLAGS.fasta_path).readlines()[1].rstrip('\n').strip()

    outdir = FLAGS.output_dir
    makedir_if_not_exists(outdir)

    if FLAGS.run_af2_wt:
        run_af2_cmd = f"python {params['run_alphafold_program']}"
        run_af2_cmd += f" --fasta_paths={FLAGS.fasta_path}"
        run_af2_cmd += f" --output_dir={outdir}"
        run_af2_cmd += f" --hhblits_binary_path={params['hhblits_program']}"
        run_af2_cmd += f" --hhsearch_binary_path={params['hhsearch_program']}"
        run_af2_cmd += f" --jackhmmer_binary_path={params['jackhmmer_program']}"
        run_af2_cmd += f" --kalign_binary_path={params['kalign_program']}"
        run_af2_cmd += f" --bfd_database_path={params['bfd_database']}"
        run_af2_cmd += f" --mgnify_database_path={params['mgnify_database']}"
        run_af2_cmd += f" --template_mmcif_dir={params['template_mmcif_dir']}"
        run_af2_cmd += f" --obsolete_pdbs_path={params['obsolete_pdbs_path']}"
        run_af2_cmd += f" --pdb70_database_path={params['pdb70_hhsuite_database']}"
        run_af2_cmd += f" --uniref30_database_path={params['uniref_db']}"
        run_af2_cmd += f" --uniref90_database_path={params['uniref90_fasta']}"
        run_af2_cmd += f" --data_dir={params['alphafold_database_dir']}"
        run_af2_cmd += f" --model_preset={FLAGS.model_preset}"
        run_af2_cmd += f" --max_template_date={FLAGS.max_template_date}"
        run_af2_cmd += f" --use_precomputed_msas={str(FLAGS.af2_use_precomputed_msas).lower()}"
        run_af2_cmd += f" --num_multimer_predictions_per_model=1"
        run_af2_cmd += f" --db_preset=full_dbs --benchmark=false --use_gpu_relax=true --logtostderr"
        try:
            # print(run_af2_cmd)
            os.system(run_af2_cmd)
        except Exception as e:
            print(e)

    if FLAGS.run_cas12fold:
        if FLAGS.cas12folddb_merge_mode == Cas12folddbMergeMode.both:
            cas12folddb_merge_mode = "both"
        elif FLAGS.cas12folddb_merge_mode == Cas12folddbMergeMode.single:
            cas12folddb_merge_mode = "single"
        else:
            cas12folddb_merge_mode = "merge"
        run_cas12fold_cmd = f"python {params['run_cas12fold_program']}"
        run_cas12fold_cmd += f" --fasta_paths={FLAGS.fasta_path}"
        run_cas12fold_cmd += f" --output_dir={outdir}"
        run_cas12fold_cmd += f" --hhblits_binary_path={params['hhblits_program']}"
        run_cas12fold_cmd += f" --hhsearch_binary_path={params['hhsearch_program']}"
        run_cas12fold_cmd += f" --jackhmmer_binary_path={params['jackhmmer_program']}"
        run_cas12fold_cmd += f" --kalign_binary_path={params['kalign_program']}"
        run_cas12fold_cmd += f" --bfd_database_path={params['bfd_database']}"
        run_cas12fold_cmd += f" --mgnify_database_path={params['mgnify_database']}"
        run_cas12fold_cmd += f" --template_mmcif_dir={params['template_mmcif_dir']}"
        run_cas12fold_cmd += f" --obsolete_pdbs_path={params['obsolete_pdbs_path']}"
        run_cas12fold_cmd += f" --pdb70_database_path={params['pdb70_hhsuite_database']}"
        run_cas12fold_cmd += f" --uniref30_database_path={params['uniref_db']}"
        run_cas12fold_cmd += f" --uniref90_database_path={params['uniref90_fasta']}"
        run_cas12fold_cmd += f" --data_dir={params['alphafold_database_dir']}"
        run_cas12fold_cmd += f" --model_preset={FLAGS.model_preset}"
        run_cas12fold_cmd += f" --max_template_date={FLAGS.max_template_date}"
        run_cas12fold_cmd += f" --use_precomputed_msas={str(FLAGS.af2_use_precomputed_msas).lower()}"
        run_cas12fold_cmd += f" --cas12folddb_merge_mode={cas12folddb_merge_mode}"
        run_cas12fold_cmd += f" --cas12folddb_database_path={params['cas12folddb_database_path']}"
        run_cas12fold_cmd += f" --max_subsequence_ratio={FLAGS.max_subsequence_ratio}"
        if flags.use_cas12fold_template_db:
            run_cas12fold_cmd += f" --cas12folddb_template_database_path={params['cas12folddb_template_database_path']}"
        run_cas12fold_cmd += f" --num_multimer_predictions_per_model=1"
        run_cas12fold_cmd += f" --db_preset=full_dbs --benchmark=false --use_gpu_relax=true --logtostderr"
        try:
            print(run_cas12fold_cmd)
            os.system(run_cas12fold_cmd)
        except Exception as e:
            print(e)


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path'
    ])
    app.run(main)