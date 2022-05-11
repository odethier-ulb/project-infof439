"""
Caretta - a multiple protein structure alignment and feature extraction suite

Implementation for the course INFO-F439
"""

import argparse
import sys
import logging
from pairwise_alignment import pairwise_alignment
from multiple_alignment import multiple_alignment
from utils.files import *
from utils.pdb import Protein
from utils.helper import get_alignment
from utils.kabsch import *


def main():
    logging.basicConfig(format='%(asctime)s | %(levelname)s: %(message)s', level=logging.NOTSET)
    args = parse_arguments()

    if args.pdb_dir is None:
        sys.exit('Error: no input directory')

    create_output_directory(args.out_dir)

    fasta_out = path.join(args.out_dir, 'alignments.fasta')
    features_out = path.join(args.out_dir, 'features.csv')
    chain = None if args.chain == 'all' else args.chain

    proteins = [Protein(p[0], p[1], chain) for p in get_files(args.pdb_dir, args.extension)]

    if len(proteins) == 0:
        sys.exit('Error: {} does not contain at least two *.{} files'.format(args.pdb_dir, args.extension))
    else:
        logging.info('Found {} proteins to align'.format(len(proteins)))
        if len(proteins) == 2:
            _, aln_a, aln_b = pairwise_alignment(proteins[0], proteins[1])
            alignments = {proteins[0].name: get_alignment(proteins[0].get_sequence(), aln_a),
                          proteins[1].name: get_alignment(proteins[1].get_sequence(), aln_b)}
            write_alignments(alignments, fasta_out)
            rmsd = 'NA'
        else:
            alignments, final_protein = multiple_alignment(proteins)
            write_alignments(alignments, fasta_out)
            aligned_coordinates = get_aligned_coordinates(final_protein, proteins, alignments)
            write_aligned_pdbs(proteins, aligned_coordinates, args.out_dir)
            write_features(alignments, proteins, features_out)
            rmsd = get_final_rmsd(alignments, proteins)

    logging.info('Done, the alignment results have been written to {}\nComputed RMSD : {} A'.format(args.out_dir, rmsd))


def parse_arguments():
    parser = argparse.ArgumentParser(description='Caretta: a multiple alignment suite meant for homologous but '
                                                 'sequentially divergent protein families',
                                     epilog='output: a multiple alignments fasta file')
    parser.add_argument('pdb_dir', type=str, help='A path to the directory containing *.pdb files to align')
    parser.add_argument('--out_dir', type=str, help='Name of the output directory', default='alignment_results')
    parser.add_argument('--extension', type=str, help='File extension to consider', default='.pdb')
    parser.add_argument('--chain', type=str, help='Protein chain to consider', default='all')
    return parser.parse_args()


if __name__ == '__main__':
    main()
