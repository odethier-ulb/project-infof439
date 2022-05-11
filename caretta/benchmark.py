"""
Provide tools for caretta benchmarking
"""

import random
import shutil
import numpy as np
import requests
import json
import xml.etree.ElementTree as ET
from time import time
from Bio.PDB import PDBParser, PDBIO
from multiple_alignment import multiple_alignment
from utils.files import *
from utils.pdb import Protein
from utils.helper import get_ratio_gapless, get_equivalence_score
from utils.kabsch import get_final_rmsd

# change me if needed
SEED_DIRECTORY = '../test_data/seeds'

HOMSTRAD_DATASET_PATH = 'E:/Cours_ULB/INFO-F439/project/Caretta/homstrad/homstrad'
HOMSTRAD_ALIGN_MTM_URL = 'https://yanglab.nankai.edu.cn/mTM-align/benchmark/homstrad/'
HOMSTRAD_ALIGN_MATT_URL = 'http://cb.csail.mit.edu/cb/matt/homstrad/'
HOMSTRAD_ALIGN_DIR = '../results/homstrad/alignments'

SABMARK_DATASET_PATH = 'E:/Cours_ULB/INFO-F439/project/Caretta/SABmark-sup/sup'
SABMARK_ALIGN_DIR = '../results/sabmark-sup/alignments'
SABMARK_ALIGN_MTM_URL = 'https://yanglab.nankai.edu.cn/mTM-align/benchmark/SABmark-sup/'
SABMARK_ALIGN_MATT_URL = 'http://cb.csail.mit.edu/cb/matt/sabmark/'

N = 25
# GROUP_SIZES = [13, 33, 63, 93]
GROUP_SIZES = [3, 6, 12, 24]


def reformat_files(dir_path):
    """ reformat files to be able to use them with dssp """
    io = PDBIO()
    for group in listdir(dir_path):
        print(group)
        for file in listdir(path.join(dir_path, group, 'pdb')):
            if file.endswith('.ent'):
                prot_path = path.join(dir_path, group, 'pdb', file)
                prot_name = file[:-4]
                out_path = path.join(dir_path, group, 'pdb', prot_name + '.atm')
                structure = PDBParser().get_structure(prot_name, prot_path)[0]
                io.set_structure(structure)
                io.save(out_path)


def generate_seeds(n, dir_path, out_path):
    """ Randomly select n proteins from the given dir to be used as seed """
    file_list = []
    for file in listdir(dir_path):
        f_path = path.join(dir_path, file)
        if path.isdir(f_path):
            file_list += get_files(f_path, '.atm')
    random.shuffle(file_list)
    for i in range(n):
        p = Protein(file_list[i][0], file_list[i][1], None)
        shutil.copy(file_list[i][1], path.join(out_path, 'seed_{}.atm'.format(len(p))))


def generate_group_from_seed(group_size, seed):
    """ Generate a group of size protein from the given seed """
    group, p = [], Protein('seed', seed, None)
    for i in range(group_size):
        coordinates = [coord * np.around(np.random.uniform(0.9, 1.1), 3) for coord in p.get_ca_coordinates()]
        group.append(Protein.from_intermediate_node(str(i), coordinates, p.get_dssp(), p.get_consensus()))
    return group


def measure_runtime(seed_dir, result_file, group_sizes):
    """ Measure the runtime for different configurations """
    with open(result_file, 'wt') as f:
        f.write('group_size;protein_length;runtime\n')
    for seed in sorted(get_files(seed_dir, '.atm'), key=lambda x: int(x[1][:-4].split('_')[-1])):
        for group_size in group_sizes:
            group = generate_group_from_seed(group_size, seed[1])
            protein_length = seed[1][:-4].split('_')[-1]
            print('Aligning {} proteins of size {}...'.format(group_size, protein_length))
            t0 = time()
            multiple_alignment(group)
            t1 = time()
            with open(result_file, 'at') as f:
                f.write('{};{};{}\n'.format(group_size, protein_length, np.around((t1 - t0) / 60., 2)))


def fetch_alignment_results_from_web(base_url, home_dir, output_dir):
    """ Retrieve homstrad results from online """
    for family in [file for file in listdir(home_dir)]:
        print('Downloading {}...'.format(family))
        # # homstrad
        # homstrad_url = '{}/{}/{}.ali'.format(base_url, family, family)
        # response = requests.get(homstrad_url)
        # open(path.join(output_dir, family + '.ali'), 'wb').write(response.content)
        # # mtm
        # mtm_url = '{}/{}/{}_result.fasta'.format(base_url, family, family)
        # response = requests.get(mtm_url)
        # open(path.join(output_dir, family + '_mtm.fasta'), 'wb').write(response.content)
        # # matt
        # matt_url = '{}/{}.fasta'.format(base_url, family)
        # response = requests.get(matt_url)
        # open(path.join(output_dir, family + '_matt.fasta'), 'wb').write(response.content)
        pass


def compute_homstrad_mtm_metrics(homstrad_dir, alignment_dir, out_file):
    """ Create the homstrad-mtm metrics """
    with open(out_file, 'wt') as f:
        f.write('family;ratio_gapless;ratio_homstrad\n')
        for family in [file for file in listdir(homstrad_dir)]:
            try:
                mtm_aln = __parse_fasta_file(path.join(alignment_dir, family + '_mtm.fasta'),
                                             lambda x: x[1:].split('.')[0])
                if len(mtm_aln) <= 2:
                    continue
                homstrad_aln = __parse_ali_file(path.join(alignment_dir, family + '.ali'))
                f.write('{};{};{}\n'.format(family, get_ratio_gapless(mtm_aln), get_equivalence_score(homstrad_aln, mtm_aln)))
            except:
                pass


def compute_homstrad_matt_metrics(homstrad_dir, alignment_dir, out_file):
    """ Create the homstrad-mtm metrics """
    with open(out_file, 'wt') as f:
        f.write('family;ratio_gapless;ratio_homstrad\n')
        for family in [file for file in listdir(homstrad_dir)]:
            try:
                mtm_aln = __parse_fasta_file(path.join(alignment_dir, family + '_matt.fasta'), None, True)
                if len(mtm_aln) <= 2:
                    continue
                homstrad_aln = __parse_ali_file(path.join(alignment_dir, family + '.ali'), True)
                f.write('{};{};{}\n'.format(family, get_ratio_gapless(mtm_aln),
                                            get_equivalence_score(homstrad_aln, mtm_aln)))
            except:
                pass


def compute_homstrad_caretta_metrics(homstrad_dir, alignment_dir, out_file):
    """ Create the homstrad-mtm metrics """
    with open(out_file, 'wt') as f:
        f.write('family;ratio_gapless;ratio_homstrad\n')
        for family in [file for file in listdir(homstrad_dir)]:
            try:
                caretta_aln = __parse_fasta_file(path.join(alignment_dir, family + '_caretta.fasta'), lambda x: x[1:])
                if len(caretta_aln) <= 2:
                    continue
                homstrad_aln = __parse_ali_file(path.join(alignment_dir, family + '.ali'))
                f.write('{};{};{}\n'.format(family, get_ratio_gapless(caretta_aln),
                                            get_equivalence_score(homstrad_aln, caretta_aln)))
            except:
                pass


def compute_homstrad_caretta_alignments(homstrad_dir, output_dir):
    """ Compute the alignments for homstrad with caretta"""
    for family in [file for file in listdir(homstrad_dir)]:
        try:
            proteins = [Protein(p[0], p[1], None) for p in get_files(path.join(homstrad_dir, family), '.atm')]
            print('Aligning {}...'.format(family))
            alignments, _ = multiple_alignment(proteins)
            write_alignments(alignments, path.join(output_dir, '{}_caretta.fasta'.format(family)))
        except Exception as e:
            with open(path.join(output_dir, '01_caretta.logs'), 'at') as f:
                f.write('ERROR with {} : {}'.format(family, e))


def compute_sabmark_caretta_alignments(sabmark_dir, output_dir):
    for group in [group for group in listdir(sabmark_dir)]:
        try:
            proteins = [Protein(p[0], p[1], None) for p in get_files(path.join(sabmark_dir, group, 'pdb'), '.atm')]
            print('Aligning {}...'.format(group))
            alignments, _ = multiple_alignment(proteins)
            rmsd = get_final_rmsd(alignments, proteins)
            write_alignments(alignments, path.join(output_dir, '{}_caretta.fasta'.format(group)))
            with open(path.join(output_dir, 'rmsd_caretta.csv'), 'at') as f:
                f.write('{};{}\n'.format(group, str(rmsd)))
        except Exception as e:
            with open(path.join(output_dir, '01_caretta.logs'), 'at') as f:
                f.write('ERROR with {} : {}\n'.format(group, e))


def extract_mtm_rmsd(output_dir):
    """ Extract rmsd from mtm """
    with open('../results/sabmark-sup/mtm.json', 'r') as f:
        data = json.load(f)
    with open(path.join(output_dir, 'rmsd_mtm.csv'), 'wt') as f:
        for group in data:
            group_name = group['Group information']['family_name']
            rmsd = group['mTM-align Result']['pc_RMSD']
            f.write('{};{}\n'.format(group_name, rmsd))


def extract_matt_rmsd(output_dir):
    """ Extract rmsd from matt """
    tree = ET.parse('../results/sabmark-sup/matt.xml')
    root = tree.getroot()
    with open(path.join(output_dir, 'rmsd_matt.csv'), 'wt') as f:
        for i in range(1, len(root)):
            group = ''.join(root[i][0].text.split())
            rmsd = root[i][2].text
            f.write('{};{}\n'.format(group, rmsd))


def compute_sabmark_mtm_metrics(sabmark_dir, alignment_dir, out_file):
    """ Create the sabmark-mtm metrics """
    rmsd = {}
    with open(path.join(alignment_dir, 'rmsd_mtm.csv'), 'rt') as f:
        for line in f:
            l_splt = line.strip().split(';')
            if len(l_splt) > 1:
                rmsd[l_splt[0]] = np.around(float(l_splt[1]), 3)
    with open(out_file, 'wt') as f:
        f.write('group;ratio_gapless;rmsd\n')
        for group in [file for file in listdir(sabmark_dir)]:
            try:
                mtm_aln = __parse_fasta_file(path.join(alignment_dir, group + '_mtm.fasta'),
                                             lambda x: x[1:].split('.')[0])
                if len(mtm_aln) <= 2:
                    continue
                f.write('{};{};{}\n'.format(group, get_ratio_gapless(mtm_aln), rmsd[group]))
            except:
                pass


def compute_sabmark_matt_metrics(sabmark_dir, alignment_dir, out_file):
    """ Create the sabmark-mtm metrics """
    rmsd = {}
    with open(path.join(alignment_dir, 'rmsd_matt.csv'), 'rt') as f:
        for line in f:
            l_splt = line.strip().split(';')
            if len(l_splt) > 1:
                rmsd[l_splt[0]] = np.around(float(l_splt[1]), 3)
    with open(out_file, 'wt') as f:
        f.write('group;ratio_gapless;rmsd\n')
        for group in [file for file in listdir(sabmark_dir)]:
            try:
                mtm_aln = __parse_fasta_file(path.join(alignment_dir, group + '_matt.fasta'), None, True)
                if len(mtm_aln) <= 2:
                    continue
                f.write('{};{};{}\n'.format(group, get_ratio_gapless(mtm_aln), rmsd[group], '.'))
            except:
                pass


def compute_sabmark_caretta_metrics(sabmark_dir, alignment_dir, out_file):
    """ Create the homstrad-mtm metrics """
    rmsd = {}
    with open(path.join(alignment_dir, 'rmsd_caretta.csv'), 'rt') as f:
        for line in f:
            l_splt = line.strip().split(';')
            if len(l_splt) > 1:
                rmsd[l_splt[0]] = np.around(float(l_splt[1]), 3)
    with open(out_file, 'wt') as f:
        f.write('group;ratio_gapless;rmsd\n')
        for group in [file for file in listdir(sabmark_dir)]:
            try:
                caretta_aln = __parse_fasta_file(path.join(alignment_dir, group + '_caretta.fasta'), lambda x: x[1:])
                if len(caretta_aln) <= 2:
                    continue
                f.write('{};{};{}\n'.format(group, get_ratio_gapless(caretta_aln), rmsd[group], '.'))
            except:
                pass


def __parse_ali_file(file_path, num_id=False):
    alignment = {}
    cur_alignment, prot_id = '', ''
    id = 1
    with open(file_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                prot_id = str(id) if num_id else line[1:].split(';')[1]
                next(f, '')
                id += 1
            elif len(prot_id) > 0 and not line.endswith('*'):
                cur_alignment += line
            elif len(prot_id) > 0 and line.endswith('*'):
                cur_alignment += line[:-1]
                alignment[prot_id] = cur_alignment
                cur_alignment, prot_id = '', ''
    return alignment


def __parse_fasta_file(file_path, extract_id, num_id=False):
    alignment = {}
    cur_alignment, prot_id = '', ''
    id = 1
    with open(file_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if len(prot_id) > 0:
                    alignment[prot_id] = cur_alignment
                    cur_alignment, prot_id = '', ''
                    id += 1
                prot_id = str(id) if num_id else extract_id(line)
            elif len(prot_id) > 0:
                cur_alignment += line
        if len(prot_id) > 0:
            alignment[prot_id] = cur_alignment
    return alignment


if __name__ == '__main__':
    # generate_seeds(N, HOMSTRAD_DATASET_PATH, SEED_DIRECTORY)
    # measure_runtime(SEED_DIRECTORY, 'runtime_results.txt', GROUP_SIZES)
    # fetch_alignment_results_from_web(HOMSTRAD_ALIGN_MATT_URL, HOMSTRAD_DATASET_PATH, HOMSTRAD_ALIGN_DIR)
    # compute_homstrad_caretta_alignments(HOMSTRAD_DATASET_PATH, HOMSTRAD_ALIGN_DIR)
    # compute_homstrad_mtm_metrics(HOMSTRAD_DATASET_PATH, HOMSTRAD_ALIGN_DIR, '../results/homstrad/homstrad_mtm_metrics.csv')
    # compute_homstrad_matt_metrics(HOMSTRAD_DATASET_PATH, HOMSTRAD_ALIGN_DIR, '../results/homstrad/homstrad_matt_metrics.csv')
    # compute_homstrad_caretta_metrics(HOMSTRAD_DATASET_PATH, HOMSTRAD_ALIGN_DIR, '../results/homstrad/homstrad_caretta_metrics.csv')
    # reformat_files(SABMARK_DATASET_PATH)
    # compute_sabmark_caretta_alignments(SABMARK_DATASET_PATH, SABMARK_ALIGN_DIR)
    # fetch_alignment_results_from_web(SABMARK_ALIGN_MATT_URL, SABMARK_DATASET_PATH, SABMARK_ALIGN_DIR)
    # extract_mtm_rmsd(SABMARK_ALIGN_DIR)
    # extract_matt_rmsd(SABMARK_ALIGN_DIR)
    # compute_sabmark_mtm_metrics(SABMARK_DATASET_PATH, SABMARK_ALIGN_DIR, '../results/sabmark-sup/sabmark_mtm_metrics.csv')
    # compute_sabmark_matt_metrics(SABMARK_DATASET_PATH, SABMARK_ALIGN_DIR, '../results/sabmark-sup/sabmark_matt_metrics.csv')
    compute_sabmark_caretta_metrics(SABMARK_DATASET_PATH, SABMARK_ALIGN_DIR, '../results/sabmark-sup/sabmark_caretta_metrics.csv')
    pass
