"""
Util. functions to manipulate files
"""

import os.path
from os import listdir, path, makedirs
from config import ALN_GAP
from utils.pdb import replace_coordinates_and_save
from utils.helper import get_aln_idx


def create_output_directory(out_dir, clean=True):
    """
    Create the output directory if it doesn't already exist

    :param out_dir: the output directory name
    :param clean: if True, will remove existing *.pdb and *.fasta
    """
    if not path.exists(out_dir):
        makedirs(out_dir)
    else:
        for f in listdir(out_dir):
            if f.lower().endswith('.pdb') or f.lower().endswith('.fasta') or f.lower().endswith('.csv'):
                os.remove(path.join(out_dir, f))


def get_files(dir_path, extension='.pdb'):
    """
    Get all *{extension} files in the given directory

    :param dir_path: path to a directory
    :param extension
    :return: a list of tupple (file_id, file_path) of each file contained in the given directory
    """
    return [(fn[:-4], path.join(dir_path, fn))
            for fn in list(filter(lambda x: x.lower().endswith(extension), listdir(dir_path)))]


def write_alignments(alignments, out_file):
    """
    Write the result of a multiple alignments

    :param alignments: a dictionary containing the alignments
    :param out_file: the output filename
    """
    with open(out_file, 'wt') as f:
        for k in sorted(alignments.keys()):
            f.write('>{}\n{}\n'.format(k, alignments[k]))


def write_features(alignments, proteins, out_file):
    """
    Write the resulting features from a multiple alignments

    :param alignments: a dictionary containing the alignments
    :param proteins: a list of proteins
    :param out_file: the output filename
    """
    with open(out_file, 'wt') as f:
        for k in sorted(alignments.keys()):
            p = list(filter(lambda x: x.name == k, proteins))[0]
            aln_idx = get_aln_idx(alignments[k])
            aln_dssp = [p.get_dssp()[aln_idx[i]] if aln_idx[i] != ALN_GAP else 'NA' for i in range(len(aln_idx))]
            f.write('{}\n'.format(';'.join(aln_dssp)))


def write_aligned_pdbs(proteins, aligned_coordinates, out_dir):
    """
    Write aligned proteins to new *.pdb files

    :param proteins: a list of pdb.Protein
    :param aligned_coordinates: the new coordinates
    :param out_dir: the output directory
    """
    for p in proteins:
        replace_coordinates_and_save(p, path.join(out_dir, p.name + '.pdb'), aligned_coordinates[p.name])
