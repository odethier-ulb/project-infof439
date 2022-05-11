"""
Util. functions to manipulate PDB files
"""

import numpy as np
from Bio.PDB import PDBParser, Selection, DSSP, PDBIO, Select
from Bio.SeqUtils import IUPACData


class Residue:
    """
    Represent interesting information of a residue
    """

    def __init__(self, coord, dssp_code, res_code='X', seq_id=None):
        """
        Create a Residue object

        :param coord: (x, y, z) coordinate of the alpha carbon
        :param dssp_code: secondary structure dssp code
        :param res_code: amino acid single letter code
        """
        self.coord = coord
        self.dssp_code = dssp_code
        self.res_code = res_code
        self.seq_id = seq_id


class Protein:
    """
    Read and extract information from a PDB file
    """

    def __init__(self, name, file_path=None, chain='A', extract=True, residues=None, consensus=None):
        """
        Create a Protein

        :param name: name of the protein
        :param file_path: path to the PDB file
        :param chain: the chain to extract
        :param extract: True if we want to extract Residues
        """
        self.name = name
        self.file_path = file_path
        self.chain = chain
        self.atoms_coordinates = []
        if extract:
            self.residues = []
            self.__extract_residues()
            self.consensus = [1 for _ in range(len(self.residues))]
        else:
            self.residues = residues
            self.consensus = consensus

    def __len__(self):
        return len(self.residues)

    @classmethod
    def from_intermediate_node(cls, name, coord, sec, con):
        """
        Create a Protein from an intermediate node
        :param name: the protein name
        :param coord: a list if alpha carbon coordinates
        :param sec: a list of secondary structure dssp codes
        :param con: a consensus row
        :return: a new pdb.Protein
        """
        residues = [Residue(coord[i], sec[i]) for i in range(len(coord))]
        return cls(name, extract=False, residues=residues, consensus=con)

    def get_len_aln(self):
        """
        :return: the number of proteins used to construct the current one
        """
        return len(self.name.split(";"))

    def get_ca_coordinates(self):
        """
        :return: the (x, y, z) coordinates of all the residues (alpha-carbon)
        """
        return np.array([r.coord for r in self.residues])

    def get_coordinates(self):
        """
        :return: the (x, y, z) coordinates of all the atoms
        """
        return self.atoms_coordinates

    def get_dssp(self):
        """
        :return: the secondary structure as the DSSP code for all the residues
        """
        return [r.dssp_code for r in self.residues]

    def set_consensus(self, cw):
        """
        Set the consensus row

        :param cw: the consensus
        """
        self.consensus = cw

    def get_consensus(self):
        """
        :return: the consensus row
        """
        return self.consensus

    def get_sequence(self):
        """
        :return: the sequence of the Protein
        """
        return ''.join(r.res_code for r in self.residues)

    def __extract_residues(self):
        """
        Extract a list of Residue from the PDB file
        """
        structure = PDBParser(QUIET=True).get_structure(self.name, self.file_path)[0]
        dssp = DSSP(structure, self.file_path)

        for chain in Selection.unfold_entities(structure, 'C'):
            if self.chain is None or chain.id == self.chain:
                for residue in Selection.unfold_entities(chain, 'R'):
                    for atom in Selection.unfold_entities(residue, 'A'):
                        self.atoms_coordinates.append(atom.coord)
                        if atom.name == 'CA':
                            self.residues.append(Residue(
                                atom.coord,
                                dssp[(chain.id, residue.get_id())][2],
                                IUPACData.protein_letters_3to1[residue.resname.upper().title()],
                                residue.id[1]))


class ChainSelect(Select):

    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        return 1 if self.chain is None or chain.id == self.chain else 0


def replace_coordinates_and_save(protein, out_path, coordinates):
    """
     Write a new pdb file with replaced coordinates

    :param protein: a pdb.Protein
    :param out_path: the file path for the new *.pdb
    :param coordinates: the new coordinates
    """
    idx, structure = 0, PDBParser(QUIET=True).get_structure(protein.name, protein.file_path)[0]
    for chain in Selection.unfold_entities(structure, 'C'):
        for residue in Selection.unfold_entities(chain, 'R'):
            for atom in Selection.unfold_entities(residue, 'A'):
                atom.coord = coordinates[idx]
                idx += 1
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, ChainSelect(protein.chain))
