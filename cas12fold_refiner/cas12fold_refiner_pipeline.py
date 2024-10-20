#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @File   : cas12fold_refiner_pipeline.py
# @Author : Feng Xu @ crazyhsu9527@gmail.com
# @Created: 2024-10-16

import os, pathlib, pickle, json
import datetime, copy
import pandas as pd
# import numpy as np
from cas12fold_refiner import cas12fold_refiner_parsers

from foldseek import Foldseek
from utils import *
from protein import complete_result
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO
from Bio.PDB import PDBList

# import os
from typing import Any, Mapping, MutableMapping, Optional, Sequence, Union
from absl import logging
from alphafold.common import residue_constants
# from alphafold.data_custom import custom_params
from alphafold.data import msa_identifiers
from alphafold.data import parsers as af2_parsers
from alphafold.data import templates as af2_template
from alphafold.data.tools import hhblits
from alphafold.data.tools import hhsearch
from alphafold.data.tools import hmmsearch
from alphafold.data.tools import jackhmmer
from cas12fold_refiner import cas12fold_refiner_templates
import numpy as np

# Internal import (7716).


class PrefilterError(Exception):
    """A base class for template prefilter exceptions."""


class DateError(PrefilterError):
    """An error indicating that the hit date was after the max allowed date."""


class AlignRatioError(PrefilterError):
    """An error indicating that the hit align ratio to the query was too small."""


class DuplicateError(PrefilterError):
    """An error indicating that the hit was an exact subsequence of the query."""


class LengthError(PrefilterError):
    """An error indicating that the hit was too short."""


# class refinement_input:
#     def __init__(self, fasta_path, pdb_path, pkl_path, msa_path):
#         self.fasta_path = fasta_path
#         self.pdb_path = pdb_path
#         self.pkl_path = pkl_path
#         self.msa_path = msa_path


FeatureDict = MutableMapping[str, np.ndarray]
TemplateSearcher = Union[hhsearch.HHSearch, hmmsearch.Hmmsearch]

def read_fasta(fileobj):
    current_sequence = ""
    current_id = None
    for line in fileobj:
        if line.startswith(">"):
            if current_id is not None:
                yield current_id, current_sequence
            current_id = line.rstrip()[1:]
            current_sequence = ""
        elif not line.startswith(";"):
            current_sequence += line.rstrip()
    yield current_id, current_sequence


def combine_a3ms(a3ms, outa3m):
    with open(outa3m, 'w') as fw:
        query_name, query_seq = None, None
        for a3m in a3ms:
            with open(a3m, 'r') as fileobj:
                for i, (seq_id, seq) in enumerate(read_fasta(fileobj)):
                    if i == 0:
                        if query_name is None and query_seq is None:
                            query_name = seq_id
                            query_seq = seq
                            fw.write(f">{seq_id}\n{seq}\n")
                        elif query_name != seq_id or query_seq != seq:
                            raise ValueError("The input a3ms don't have the same query name or query sequences")
                    else:
                        fw.write(f">{seq_id}\n{seq}\n")


def make_sequence_features(
        sequence: str, description: str, num_res: int) -> FeatureDict:
    """Constructs a feature dict of sequence features."""
    features = {}
    features['aatype'] = residue_constants.sequence_to_onehot(
        sequence=sequence,
        mapping=residue_constants.restype_order_with_x,
        map_unknown_to_x=True)
    features['between_segment_residues'] = np.zeros((num_res,), dtype=np.int32)
    features['domain_name'] = np.array([description.encode('utf-8')],
                                       dtype=np.object_)
    features['residue_index'] = np.array(range(num_res), dtype=np.int32)
    features['seq_length'] = np.array([num_res] * num_res, dtype=np.int32)
    features['sequence'] = np.array([sequence.encode('utf-8')], dtype=np.object_)
    return features


def make_msa_features(msas: Sequence[af2_parsers.Msa], msa_output_dir: str, msa_save_path: str, filter=True) -> FeatureDict:
    """Constructs a feature dict of MSA features."""
    if not msas:
        raise ValueError('At least one MSA must be provided.')

    int_msa = []
    deletion_matrix = []
    species_ids = []
    seen_desc = []
    seen_sequences = []
    for msa_index, msa in enumerate(msas):
        if not msa:
            raise ValueError(f'MSA {msa_index} must contain at least one sequence.')
        for sequence_index, sequence in enumerate(msa.sequences):
            if filter and sequence in seen_sequences:
                continue
            seen_sequences += [sequence]
            seen_desc += [msa.descriptions[sequence_index]]
            int_msa.append([residue_constants.HHBLITS_AA_TO_ID[res] for res in sequence])
            deletion_matrix.append(msa.deletion_matrix[sequence_index])
            identifiers = msa_identifiers.get_identifiers(msa.descriptions[sequence_index])
            species_ids.append(identifiers.species_id.encode('utf-8'))

    num_res = len(msas[0].sequences[0])
    num_alignments = len(int_msa)
    features = {}
    features['deletion_matrix_int'] = np.array(deletion_matrix, dtype=np.int32)
    features['msa'] = np.array(int_msa, dtype=np.int32)
    features['num_alignments'] = np.array([num_alignments] * num_res, dtype=np.int32)
    features['msa_species_identifiers'] = np.array(species_ids, dtype=np.object_)

    with open(os.path.join(msa_output_dir, msa_save_path), 'w') as fw:
        for (desc, seq) in zip(seen_desc, seen_sequences):
            fw.write(f'>{desc}\n{seq}\n')

    return features


def merge_msa(af2_msas_dir, msa_output_dir, msa_name):
    bfd_msa_file = os.path.join(af2_msas_dir, "bfd_uniref_hits.a3m")
    mgnify_msa_file = os.path.join(af2_msas_dir, "mgnify_hits.sto")
    uniref90_msa_file = os.path.join(af2_msas_dir, "uniref90_hits.sto")
    cas12folddb_msa_file = os.path.join(af2_msas_dir, "cas12folddb_hits.sto")
    with open(bfd_msa_file, 'r') as f:
        bfd_msa = af2_parsers.parse_a3m(f.read())
    with open(mgnify_msa_file, 'r') as f:
        mgnify_msa = af2_parsers.parse_stockholm(f.read())
    with open(uniref90_msa_file, 'r') as f:
        uniref90_msa = af2_parsers.parse_stockholm(f.read())
    if os.path.exists(cas12folddb_msa_file):
        with open(cas12folddb_msa_file, 'r') as f:
            cas12folddb_msa = af2_parsers.parse_stockholm(f.read())
        make_msa_features((uniref90_msa, bfd_msa, mgnify_msa, cas12folddb_msa), msa_output_dir, msa_save_path=msa_name)
    else:
        make_msa_features((uniref90_msa, bfd_msa, mgnify_msa), msa_output_dir, msa_save_path=msa_name)


def run_msa_tool(msa_runner, input_fasta_path: str, msa_out_path: str,
                 msa_format: str, use_precomputed_msas: bool,
                 max_sto_sequences: Optional[int] = None
                 ) -> Mapping[str, Any]:
    """Runs an MSA tool, checking if output already exists first."""
    if not use_precomputed_msas or not os.path.exists(msa_out_path):
        if msa_format == 'sto' and max_sto_sequences is not None:
            result = msa_runner.query(input_fasta_path, max_sto_sequences)[0]  # pytype: disable=wrong-arg-count
        else:
            result = msa_runner.query(input_fasta_path)[0]
        with open(msa_out_path, 'w') as f:
            f.write(result[msa_format])
    else:
        logging.warning('Reading MSA from file %s', msa_out_path)
        if msa_format == 'sto' and max_sto_sequences is not None:
            precomputed_msa = af2_parsers.truncate_stockholm_msa(
                msa_out_path, max_sto_sequences)
            result = {'sto': precomputed_msa}
        else:
            with open(msa_out_path, 'r') as f:
                result = {msa_format: f.read()}
    return result


class DataPipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self,
                 jackhmmer_binary_path: str,
                 hhblits_binary_path: str,
                 uniref90_database_path: str,
                 mgnify_database_path: str,
                 bfd_database_path: Optional[str],
                 uniref30_database_path: Optional[str],
                 small_bfd_database_path: Optional[str],
                 template_searcher: TemplateSearcher,
                 template_featurizer: cas12fold_refiner_templates.CustomizedMonomerHitFeaturizer,
                 use_small_bfd: bool,
                 mgnify_max_hits: int = 501,
                 uniref_max_hits: int = 10000,
                 use_precomputed_msas: bool = False):
        """Initializes the data pipeline."""
        self._use_small_bfd = use_small_bfd
        self.jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
            binary_path=jackhmmer_binary_path,
            database_path=uniref90_database_path)
        if use_small_bfd:
            self.jackhmmer_small_bfd_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=small_bfd_database_path)
        else:
            self.hhblits_bfd_uniref_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[bfd_database_path, uniref30_database_path])
        self.jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
            binary_path=jackhmmer_binary_path,
            database_path=mgnify_database_path)
        self.template_searcher = template_searcher
        self.template_featurizer = template_featurizer
        self.mgnify_max_hits = mgnify_max_hits
        self.uniref_max_hits = uniref_max_hits
        self.use_precomputed_msas = use_precomputed_msas

    def process(self, input_fasta_path: str, msa_output_dir: str, template_output_dir: str, custom_inputs, template_select_foldseek_global: bool = False) -> FeatureDict:
        """Runs alignment tools on the input sequence and creates features."""
        with open(input_fasta_path) as f:
            input_fasta_str = f.read()
        input_seqs, input_descs = af2_parsers.parse_fasta(input_fasta_str)
        if len(input_seqs) != 1:
            raise ValueError(
                f'More than one input sequence found in {input_fasta_path}.')
        input_sequence = input_seqs[0]
        input_description = input_descs[0]
        num_res = len(input_sequence)

        sequence_features = make_sequence_features(
            sequence=input_sequence,
            description=input_description,
            num_res=num_res)

        custom_result = None
        if custom_inputs.custom_msa is not None:
            if custom_inputs.custom_msa.find('.a3m') > 0:
                custom_msa_out_path = os.path.join(msa_output_dir, 'custom.a3m')
                os.system(f"cp {custom_inputs.custom_msa} {custom_msa_out_path}")
                with open(custom_msa_out_path, 'r') as f:
                    custom_result = af2_parsers.parse_a3m(f.read())
            else:
                custom_msa_out_path = os.path.join(msa_output_dir, 'custom.sto')
                os.system(f"cp {custom_inputs.custom_msa} {custom_msa_out_path}")
                with open(custom_msa_out_path, 'r') as f:
                    custom_result = af2_parsers.parse_stockholm(f.read())

        msa_features = None
        jackhmmer_uniref90_result = None
        if custom_result is None:
            uniref90_out_path = os.path.join(msa_output_dir, 'uniref90_hits.sto')
            if custom_inputs.uniref90_sto is not None:
                os.system(f"cp {custom_inputs.uniref90_sto} {uniref90_out_path}")
            jackhmmer_uniref90_result = run_msa_tool(
                msa_runner=self.jackhmmer_uniref90_runner,
                input_fasta_path=input_fasta_path,
                msa_out_path=uniref90_out_path,
                msa_format='sto',
                use_precomputed_msas=self.use_precomputed_msas,
                max_sto_sequences=self.uniref_max_hits)

            uniref90_msa = af2_parsers.parse_stockholm(jackhmmer_uniref90_result['sto'])
            uniref90_msa = uniref90_msa.truncate(max_seqs=self.uniref_max_hits)

            mgnify_out_path = os.path.join(msa_output_dir, 'mgnify_hits.sto')
            if custom_inputs.mgnify_sto is not None:
                os.system(f"cp {custom_inputs.mgnify_sto} {mgnify_out_path}")
            jackhmmer_mgnify_result = run_msa_tool(
                msa_runner=self.jackhmmer_mgnify_runner,
                input_fasta_path=input_fasta_path,
                msa_out_path=mgnify_out_path,
                msa_format='sto',
                use_precomputed_msas=self.use_precomputed_msas,
                max_sto_sequences=self.mgnify_max_hits)

            mgnify_msa = af2_parsers.parse_stockholm(jackhmmer_mgnify_result['sto'])
            mgnify_msa = mgnify_msa.truncate(max_seqs=self.mgnify_max_hits)

            bfd_out_path = os.path.join(msa_output_dir, 'bfd_uniref_hits.a3m')
            if custom_inputs.bfd_uniref_a3m is not None:
                os.system(f"cp {custom_inputs.bfd_uniref_a3m} {bfd_out_path}")
            elif custom_inputs.bfd_a3m is not None and custom_inputs.uniref_a3m is not None:
                combine_a3ms([custom_inputs.bfd_a3m, custom_inputs.uniref_a3m], bfd_out_path)

            hhblits_bfd_uniref_result = run_msa_tool(
                msa_runner=self.hhblits_bfd_uniref_runner,
                input_fasta_path=input_fasta_path,
                msa_out_path=bfd_out_path,
                msa_format='a3m',
                use_precomputed_msas=self.use_precomputed_msas)
            bfd_msa = af2_parsers.parse_a3m(hhblits_bfd_uniref_result['a3m'])
            msa_features = make_msa_features((uniref90_msa, bfd_msa, mgnify_msa), msa_output_dir,
                                             msa_save_path='monomer_final.a3m')
            logging.info('Uniref90 MSA size: %d sequences.', len(uniref90_msa))
            logging.info('BFD MSA size: %d sequences.', len(bfd_msa))
            logging.info('MGnify MSA size: %d sequences.', len(mgnify_msa))
        else:
            logging.info('Custom MSA size: %d sequences.', len(custom_result))
            msa_features = make_msa_features([custom_result], msa_output_dir, msa_save_path='monomer_final.a3m')

        logging.info('Final (deduplicated) MSA size: %d sequences.',
                     msa_features['num_alignments'][0])

        templates_result_features = None
        if custom_inputs.notemplate:
            templates_result_features = mk_mock_template(input_sequence)
        elif custom_inputs.temp_struct_csv is not None:
            templates_result_features = self.template_featurizer.get_templates(query_sequence=input_sequence,
                                                                               template_pdb_dir=template_output_dir,
                                                                               hits_file=custom_inputs.temp_struct_csv,
                                                                               template_select_foldseek_global=template_select_foldseek_global).features
        else:
            if jackhmmer_uniref90_result is None:
                uniref90_out_path = os.path.join(msa_output_dir, 'uniref90_hits.sto')
                if custom_inputs.uniref90_sto is not None:
                    os.system(f"cp {custom_inputs.uniref90_sto} {uniref90_out_path}")
                jackhmmer_uniref90_result = run_msa_tool(
                    msa_runner=self.jackhmmer_uniref90_runner,
                    input_fasta_path=input_fasta_path,
                    msa_out_path=uniref90_out_path,
                    msa_format='sto',
                    use_precomputed_msas=self.use_precomputed_msas,
                    max_sto_sequences=self.uniref_max_hits)

            msa_for_templates = jackhmmer_uniref90_result['sto']
            msa_for_templates = af2_parsers.deduplicate_stockholm_msa(msa_for_templates)
            msa_for_templates = af2_parsers.remove_empty_columns_from_stockholm_msa(msa_for_templates)

            if self.template_searcher.input_format == 'sto':
                pdb_templates_result = self.template_searcher.query(msa_for_templates)
            elif self.template_searcher.input_format == 'a3m':
                uniref90_msa_as_a3m = af2_parsers.convert_stockholm_to_a3m(msa_for_templates)
                pdb_templates_result = self.template_searcher.query(uniref90_msa_as_a3m)
            else:
                raise ValueError('Unrecognized template input format: '
                                 f'{self.template_searcher.input_format}')

            pdb_hits_out_path = os.path.join(
                msa_output_dir, f'pdb_hits.{self.template_searcher.output_format}')
            with open(pdb_hits_out_path, 'w') as f:
                f.write(pdb_templates_result)

            pdb_template_hits = self.template_searcher.get_template_hits(
                output_string=pdb_templates_result, input_sequence=input_sequence)

            templates_result = self.template_featurizer.get_templates(
                query_sequence=input_sequence,
                hits=pdb_template_hits)
            templates_result_features = templates_result.features

        if custom_inputs.notemplate:
            logging.info('Total number of templates (NB: this can include bad '
                         'templates and is later filtered to top 4): %d.',
                         len(templates_result_features['template_domain_names']))
        else:
            logging.info('Total number of templates (NB: this can include bad '
                         'templates and is later filtered to top 4): %d.',
                         templates_result_features['template_domain_names'].shape[0])

        return {**sequence_features, **msa_features, **templates_result_features}


def assess_foldseek_hit(
        hit: cas12fold_refiner_parsers.TemplateHit,
        query_sequence: str,
        max_subsequence_ratio: float = 0.95,
        min_align_ratio: float = 0.1) -> bool:

    aligned_cols = hit.aligned_cols
    align_ratio = aligned_cols / len(query_sequence)

    template_sequence = hit.hit_sequence.replace('-', '')
    length_ratio = float(len(template_sequence)) / len(query_sequence)

    duplicate = (template_sequence in query_sequence and
                 length_ratio > max_subsequence_ratio)

    if align_ratio <= min_align_ratio:
        raise AlignRatioError('Proportion of residues aligned to query too small. '
                              f'Align ratio: {align_ratio}.')

    if duplicate:
        raise DuplicateError('Template is an exact subsequence of query with large '
                             f'coverage. Length ratio: {length_ratio}.')

    if len(template_sequence) < 10:
        raise LengthError(f'Template too short. Length: {len(template_sequence)}.')

    return True


class cas12fold_refinement_iterative_pipeline:

    def __init__(self, params, max_template_count=50):

        self.params = params

        self.max_iteration = 5

        self.max_template_count = max_template_count

        # release_date_df = pd.read_csv(params['pdb_release_date_file'])
        # self._release_dates = dict(zip(release_date_df['pdbcode'], release_date_df['release_date']))
        self._release_dates = {}
        self._max_template_date = datetime.datetime.strptime(params['max_template_date'], '%Y-%m-%d')

    def search_templates_old(self, inpdb, outdir):
        makedir_if_not_exists(outdir)
        foldseek_program = self.params['foldseek_program']
        foldseek_pdb_database = self.params['foldseek_pdb_database']
        foldseek_af_database = self.params['foldseek_af_database']
        foldseek_runner = Foldseek(binary_path=foldseek_program, pdb_database=foldseek_pdb_database,
                                   max_template_date=self._max_template_date, release_dates=self._release_dates,
                                   other_databases=[foldseek_af_database])
        return foldseek_runner.query(pdb=inpdb, outdir=outdir, progressive_threshold=2000)

    def search_templates(self, inpdb, outdir, result_overwrite):
        makedir_if_not_exists(outdir)
        foldseek_program = self.params['foldseek_program']
        rcsb_foldseek_database = self.params['rcsb_foldseek_database']
        cas12_foldseek_database = self.params['cas12_foldseek_database']
        foldseek_runner = Foldseek(binary_path=foldseek_program, pdb_database=rcsb_foldseek_database,
                                   max_template_date=self._max_template_date, release_dates=self._release_dates,
                                   other_databases=[cas12_foldseek_database])
        return foldseek_runner.query(pdb=inpdb, outdir=outdir, progressive_threshold=2000,
                                     result_overwrite=result_overwrite)

    def check_and_rank_templates(self, template_result, outfile, query_sequence):

        evalue_keep_indices = []
        for i in range(len(template_result['local_alignment'])):
            hit = cas12fold_refiner_parsers.TemplateHit(index=i,
                              name=template_result['local_alignment'].loc[i, 'target'].split('.')[0],
                              aligned_cols=int(template_result['local_alignment'].loc[i, 'alnlen']),
                              query=template_result['local_alignment'].loc[i, 'qaln'],
                              hit_sequence=template_result['local_alignment'].loc[i, 'taln'],
                              indices_query=build_alignment_indices(template_result['local_alignment'].loc[i, 'qaln'],
                                                                    template_result['local_alignment'].loc[i, 'qstart']),
                              indices_hit=build_alignment_indices(template_result['local_alignment'].loc[i, 'taln'],
                                                                  template_result['local_alignment'].loc[i, 'tstart']),
                              sum_probs=0.0)
            try:
                assess_foldseek_hit(hit=hit, query_sequence=query_sequence)
            except PrefilterError as e:
                msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                print(msg)
                continue
            evalue_keep_indices += [i]

        tmscore_keep_indices = []
        for i in range(len(template_result['global_alignment'])):
            hit = cas12fold_refiner_parsers.TemplateHit(index=i,
                              name=template_result['global_alignment'].loc[i, 'target'].split('.')[0],
                              aligned_cols=int(template_result['global_alignment'].loc[i, 'alnlen']),
                              query=template_result['global_alignment'].loc[i, 'qaln'],
                              hit_sequence=template_result['global_alignment'].loc[i, 'taln'],
                              indices_query=build_alignment_indices(template_result['global_alignment'].loc[i, 'qaln'],
                                                                    template_result['global_alignment'].loc[
                                                                        i, 'qstart']),
                              indices_hit=build_alignment_indices(template_result['global_alignment'].loc[i, 'taln'],
                                                                  template_result['global_alignment'].loc[i, 'tstart']),
                              sum_probs=0.0)
            try:
                assess_foldseek_hit(hit=hit, query_sequence=query_sequence)
            except PrefilterError as e:
                msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                print(msg)
                continue
            tmscore_keep_indices += [i]

        if len(evalue_keep_indices) == 0 and len(tmscore_keep_indices) == 0:
            return False

        evalue_thresholds = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
        tmscore_thresholds = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3]

        templates_sorted = pd.DataFrame(
            columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])

        evalue_af_indices = []
        evalue_pdb_indices = []
        tmscore_af_indices = []
        tmscore_pdb_indices = []
        for evalue_threshold, tmscore_threshold in zip(evalue_thresholds, tmscore_thresholds):
            evalue_af_indices = []
            evalue_pdb_indices = []
            for i in evalue_keep_indices:
                target = template_result['local_alignment'].loc[i, 'target']
                evalue = float(template_result['local_alignment'].loc[i, 'evalue'])
                if evalue < evalue_threshold:
                    if target.find('.atom.gz') > 0:
                        evalue_pdb_indices += [i]
                    else:
                        evalue_af_indices += [i]

            tmscore_af_indices = []
            tmscore_pdb_indices = []
            for i in tmscore_keep_indices:
                target = template_result['global_alignment'].loc[i, 'target']
                evalue = float(template_result['global_alignment'].loc[i, 'evalue'])
                if evalue > tmscore_threshold:
                    if target.find('.atom.gz') > 0:
                        tmscore_pdb_indices += [i]
                    else:
                        tmscore_af_indices += [i]

            if len(evalue_af_indices) + len(evalue_pdb_indices) \
                    + len(tmscore_af_indices) + len(tmscore_pdb_indices) >= self.max_template_count:
                break

        templates_sorted = copy.deepcopy(template_result['local_alignment'].iloc[evalue_pdb_indices])
        templates_sorted = templates_sorted.append(
            copy.deepcopy(template_result['global_alignment'].iloc[tmscore_pdb_indices]))
        templates_sorted = templates_sorted.append(
            copy.deepcopy(template_result['local_alignment'].iloc[evalue_af_indices]))
        templates_sorted = templates_sorted.append(
            copy.deepcopy(template_result['global_alignment'].iloc[tmscore_af_indices]))

        templates_sorted.drop(templates_sorted.filter(regex="Unnamed"), axis=1, inplace=True)
        templates_sorted.reset_index(inplace=True, drop=True)
        templates_sorted.to_csv(outfile, sep='\t')
        return True

    def generate_msa_from_templates(self, fasta_file, start_msa, template_file, outfile):
        # print(fasta_file, start_msa, template_file, outfile)
        targetname = None
        seq = ''
        for line in open(fasta_file):
            line = line.rstrip('\n')
            if line.startswith('>'):
                targetname = line[1:]
            else:
                seq += line

        templates = pd.read_csv(template_file, sep='\t')

        alignments = {targetname: seq}
        seen_seq = []
        for i in range(len(templates)):
            target = templates.loc[i, 'target']
            qaln = templates.loc[i, 'qaln']
            qstart = int(templates.loc[i, 'qstart'])
            qend = int(templates.loc[i, 'qend'])
            taln = templates.loc[i, 'taln']
            tstart = templates.loc[i, 'tstart']
            tend = templates.loc[i, 'tend']

            query_non_gaps = [res != '-' for res in qaln]
            out_sequence = ''.join(convert_taln_seq_to_a3m(query_non_gaps, taln))

            aln_full = ['-'] * len(seq)
            aln_full[qstart - 1:qend] = out_sequence
            taln_full_seq = ''.join(aln_full)
            if taln_full_seq in seen_seq:
                continue
            alignments[target] = taln_full_seq
            seen_seq += [taln_full_seq]

        fasta_chunks = (f">{k}\n{alignments[k]}" for k in alignments)

        with open(outfile + '.temp', 'w') as fw:
            fw.write('\n'.join(fasta_chunks) + '\n')

        combine_a3ms([start_msa, f"{outfile}.temp"], f"{outfile}.comb")

        cmd = f"{self.params['hhfilter_program']} -diff 50000 -i {outfile}.comb -o {outfile} -id 90"

        os.system(cmd)

    def copy_atoms_and_unzip(self, template_csv, outdir):
        os.chdir(outdir)
        templates = pd.read_csv(template_csv, sep='\t')
        for i in range(len(templates)):
            template_pdb = templates.loc[i, 'target']
            if template_pdb.find('.pdb.gz') > 0:
                template_path = os.path.join(self.params['foldseek_af_database_dir'], template_pdb)
                os.system(f"cp {template_path} {outdir}")
            else:
                template_path = os.path.join(self.params['foldseek_pdb_database_dir'], template_pdb)
                os.system(f"cp {template_path} {outdir}")
            os.system(f"gunzip -f {template_pdb}")

    def extract_target_pdbs(self, template_csv, outdir):
        os.chdir(outdir)
        templates = pd.read_csv(template_csv, sep='\t')
        templates[["target_id", "target_group"]] = templates.target.str.split('.', n=1, expand=True)

        afdb50_target_list = os.path.join(outdir, "afdb50_target.lst")
        afdb50_templates = templates[templates['target_id'].str.contains('AF')]
        if not afdb50_templates.empty:
            afdb50_templates.target_id.to_csv(afdb50_target_list, header=None, index=False)
            foldcomp_cmd = f"{self.params['foldcomp_program']} decompress -t 10 -y --id-list {afdb50_target_list} {self.params['cas12_afdb50_foldcomp_database']} {outdir} 1>/dev/null 2>&1"
            os.system(foldcomp_cmd)

        # afdb50_target_list = os.path.join(outdir, "afdb50_target.lst")
        rcsb_templates = templates[~templates['target_id'].str.contains('AF')]
        if not rcsb_templates.empty:
            rcsb_templates_copy = copy.deepcopy(rcsb_templates)
            rcsb_templates_copy[["rcsb_id", "rcsb_chain"]] = rcsb_templates_copy.target_id.str.split('_', n=1, expand=True)
            for template_pdb in rcsb_templates_copy.rcsb_id.unique():
                # template_cif_path = os.path.join(self.params['template_mmcif_dir'], f"{template_pdb}.cif")
                # template_pdb_path = f"{outdir}/{template_pdb}.pdb"
                chain_id = rcsb_templates_copy.loc[rcsb_templates_copy.rcsb_id == template_pdb].rcsb_chain.unique()[0]

                pdbl = PDBList(server='http://files.wwpdb.org')
                if not os.path.exists(f"{outdir}/pdb{template_pdb}.ent"):
                    pdbl.retrieve_pdb_file(template_pdb, pdir=f"{outdir}", file_format='pdb', overwrite=True)
                os.system(f"{self.params['pdb_selchain_program']} -{chain_id} {outdir}/pdb{template_pdb}.ent > {outdir}/{template_pdb}.pdb")

            # cif2pdb(template_cif_path, template_pdb_path)
            # os.system(f"{self.params['cif2pdb_program']} -mol 1 {template_path} {outdir}/{template_pdb}.pdb")
            # os.system(f"cp {template_path} {outdir}")
        #     pass
        # for i in range(len(rcsb_templates_copy)):
        #     template_pdb = rcsb_templates.loc[i, 'target']
        #     if template_pdb.find('.pdb.gz') > 0:
        #         template_path = os.path.join(self.params['foldseek_af_database_dir'], template_pdb)
        #         os.system(f"cp {template_path} {outdir}")
        #     else:
        #         template_path = os.path.join(self.params['foldseek_pdb_database_dir'], template_pdb)
        #         os.system(f"cp {template_path} {outdir}")
        #     os.system(f"gunzip -f {template_pdb}")


    def search_single(self, fasta_path, pdb_path, pkl_path, msa_path, outdir, result_overwrite, gpu_device=0,
                      template_select_foldseek_global=False):

        query_sequence = ""
        for line in open(fasta_path):
            line = line.rstrip('\n')
            if line.startswith('>'):
                continue
            else:
                query_sequence += line

        makedir_if_not_exists(outdir)

        cwd = os.getcwd()

        ref_start_pdb = pdb_path
        ref_start_pkl = pkl_path
        ref_start_msa = msa_path

        model_iteration_scores = []

        print(f"Start to refine {pdb_path}")

        for num_iteration in range(self.max_iteration)[0:1]:
            os.chdir(cwd)
            current_work_dir = os.path.join(outdir, f"iteration{num_iteration + 1}")
            makedir_if_not_exists(current_work_dir)

            start_pdb = os.path.join(current_work_dir, "start.pdb")
            start_msa = os.path.join(current_work_dir, "start.a3m")
            start_pkl = os.path.join(current_work_dir, "start.pkl")

            os.system(f"cp {ref_start_pdb} {start_pdb}")
            os.system(f"cp {ref_start_msa} {start_msa}")
            os.system(f"cp {ref_start_pkl} {start_pkl}")

            with open(ref_start_pkl, 'rb') as f:
                ref_avg_lddt = np.mean(pickle.load(f)['plddt'])

            model_iteration_scores += [ref_avg_lddt]

            out_model_dir = os.path.join(current_work_dir, "alphafold")
            if not complete_result(out_model_dir, 5 * int(self.params['num_monomer_predictions_per_model'])) or result_overwrite:

                foldseek_res = self.search_templates(inpdb=start_pdb, outdir=os.path.join(current_work_dir, 'foldseek'),
                                                     result_overwrite=result_overwrite)

                if not self.check_and_rank_templates(foldseek_res, os.path.join(current_work_dir, "structure_templates.csv"),
                                                     query_sequence):
                    print(f"Cannot find any templates in iteration {num_iteration + 1}")
                    break

                self.generate_msa_from_templates(fasta_file=fasta_path,
                                                 template_file=os.path.join(current_work_dir, "structure_templates.csv"),
                                                 start_msa=start_msa,
                                                 outfile=os.path.join(current_work_dir, f"iteration{num_iteration + 1}.a3m"))

                out_template_dir = os.path.join(current_work_dir, "template_pdbs")
                makedir_if_not_exists(out_template_dir)
                self.extract_target_pdbs(template_csv=os.path.join(current_work_dir, "structure_templates.csv"),
                                         outdir=out_template_dir)

                makedir_if_not_exists(out_model_dir)
                custom_msa = os.path.join(current_work_dir, f"iteration{num_iteration + 1}.a3m")
                temp_struct_csv = os.path.join(current_work_dir, "structure_templates.csv")
                os.chdir(cwd)
                cmd = f"python {self.params['run_cas12fold_refiner_program']} " \
                      f"--fasta_path={fasta_path} " \
                      f"--env_dir={self.params['alphafold_env_dir']} " \
                      f"--database_dir={self.params['alphafold_database_dir']} " \
                      f"--monomer_num_ensemble={self.params['monomer_num_ensemble']} " \
                      f"--monomer_num_recycle={self.params['monomer_num_recycle']} " \
                      f"--num_monomer_predictions_per_model={self.params['num_monomer_predictions_per_model']} " \
                      f"--model_preset={self.params['monomer_model_preset']} " \
                      f"--benchmark={self.params['alphafold_benchmark']} " \
                      f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                      f"--models_to_relax={self.params['models_to_relax']} " \
                      f"--max_template_date={self.params['max_template_date']} " \
                      f"--custom_msa={custom_msa} " \
                      f"--temp_struct_csv={temp_struct_csv} " \
                      f"--struct_atom_dir={out_template_dir} " \
                      f"--output_dir={out_model_dir} --gpu_device={gpu_device}"
                if template_select_foldseek_global:
                    cmd = cmd + " --template_select_foldseek_global"
                if result_overwrite:
                    cmd = cmd + " --result_overwrite"

                try:
                    # os.chdir(self.params['alphafold_program_dir'])
                    # print(cmd)
                    os.system(cmd)
                except Exception as e:
                    print(e)
        #
        #     new_ranking_json_file = os.path.join(out_model_dir, "ranking_debug.json")
        #     new_ranking_json = json.loads(open(new_ranking_json_file).read())
        #     max_lddt_score = new_ranking_json["plddts"][list(new_ranking_json["order"])[0]]
        #
        #     print(f'#########Iteration: {num_iteration + 1}#############')
        #     print(f"plddt before: {ref_avg_lddt}")
        #     print(f"plddt after: {max_lddt_score}")
        #     if max_lddt_score > ref_avg_lddt:
        #         print("Continue to refine")
        #         ref_start_pdb = os.path.join(out_model_dir, "ranked_0.pdb")
        #         model_name = list(new_ranking_json["order"])[0]
        #         ref_start_pkl = os.path.join(out_model_dir, f"result_{model_name}.pkl")
        #         ref_start_msa = os.path.join(out_model_dir, 'msas', "monomer_final.a3m")
        #         print('##################################################')
        #         if num_iteration + 1 >= self.max_iteration:
        #             print("Reach maximum iteration")
        #             model_iteration_scores += [max_lddt_score]
        #     else:
        #         # keep the models in iteration 1 even through the plddt score decreases
        #         if num_iteration == 0:
        #             ref_start_pdb = os.path.join(out_model_dir, "ranked_0.pdb")
        #             model_name = list(new_ranking_json["order"])[0]
        #             ref_start_pkl = os.path.join(out_model_dir, f"result_{model_name}.pkl")
        #             ref_start_msa = os.path.join(out_model_dir, 'msas', "monomer_final.a3m")
        #             model_iteration_scores += [max_lddt_score]
        #         break
        #
        # while len(model_iteration_scores) <= self.max_iteration:
        #     model_iteration_scores += [0]
        #
        # print(model_iteration_scores)
        # df = pd.DataFrame(model_iteration_scores)
        # df.to_csv(os.path.join(outdir, 'summary.csv'))
        #
        # final_model_dir = os.path.join(outdir, 'final')
        #
        # makedir_if_not_exists(final_model_dir)
        #
        # os.system("cp " + ref_start_pdb + " " + os.path.join(final_model_dir, "final.pdb"))
        # os.system("cp " + ref_start_pkl + " " + os.path.join(final_model_dir, "final.pkl"))
        # os.system("cp " + ref_start_msa + " " + os.path.join(final_model_dir, "final.a3m"))
        # os.chdir(cwd)
        #
        # return final_model_dir


# class cas12fold_refinement_iteration:
#     def __init__(self, params):
#         self.params = params
#
#     def search(self, refinement_inputs, outdir):
#         result_dirs = []
#
#         pipeline = cas12_iterative_refinement_pipeline(self.params)
#
#         for refine_param in refinement_inputs:
#             result_dir = pipeline.search_single(fasta_path=refine_param.fasta_path, pdb_path=refine_param.pdb_path,
#                                                 pkl_path=refine_param.pkl_path, msa_path=refine_param.msa_path,
#                                                 outdir=os.path.join(outdir, pathlib.Path(refine_param.pdb_path).stem))
#             result_dirs += [result_dir]
#
#         return result_dirs


# cas12_iterative_refinement_pipeline_controller
def cas12fold_refiner(params, refinement_input, outdir, finaldir, prefix, result_overwrite, gpu_device, template_select_foldseek_global):
    pipeline = cas12fold_refinement_iterative_pipeline(params)
    pipeline.search_single(fasta_path=refinement_input.fasta_path, pdb_path=refinement_input.pdb_path,
                           pkl_path=refinement_input.pkl_path, msa_path=refinement_input.msa_path,
                           outdir=os.path.join(outdir, pathlib.Path(refinement_input.pdb_path).stem),
                           result_overwrite=result_overwrite, gpu_device=gpu_device,
                           template_select_foldseek_global=template_select_foldseek_global)
    # for refine_param in refinement_inputs:
    #     pipeline.search_single(fasta_path=refine_param.fasta_path, pdb_path=refine_param.pdb_path,
    #                            pkl_path=refine_param.pkl_path, msa_path=refine_param.msa_path,
    #                            outdir=os.path.join(outdir, pathlib.Path(refine_param.pdb_path).stem))


    # pipeline = cas12fold_refinement_iteration(params=params)
    # pipeline.search(refinement_inputs=refinement_inputs, outdir=outdir)

    # makedir_if_not_exists(finaldir)
    #
    # pipeline = iterative_refine_pipeline.Monomer_refinement_model_selection(params)
    # pipeline.select_v1(indir=outdir, outdir=finaldir, prefix=prefix)