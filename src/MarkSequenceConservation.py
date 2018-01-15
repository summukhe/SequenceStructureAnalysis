import os
import argparse
import numpy as np
import PDBProcessor
import matplotlib.pyplot as plt


def is_valid_file( parser, f ):
    if os.path.isfile(f):
        return f
    parser.error('Error: invalid file location [%s]' % f)


def is_new_file(parser, f):
    topDir = os.path.dirname(os.path.realpath(f))
    fileName = os.path.realpath(f)
    if os.path.isdir(topDir) and (os.path.exists(fileName) is False):
        return f
    parser.error('Error: invalid out file name [%s]' % f)


if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description= "Utility reads a sequence alignment, " +
                                                  "for corresponding sequence tag, " +
                                                  "PDB structure, calculate " +
                                                  "the conservation using alignment " +
                                                  "mark them on the PDB bfactor." )

    parser.add_argument('--alignment', dest="alignment", metavar="FILE", required=True,
                        help="sequence alignment file input in fasta format",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument('--seq-tag' , dest='target', metavar='ID', required=True,
                        help="ID of the sequence in the alignment",
                        type=str)

    parser.add_argument('--pdb', dest='pdb', metavar='FILE', required=True,
                        help="ID full / relative path to the pdb structure",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument('--out', dest='out', metavar='FILE', required=True,
                        help="Output pdb file name",
                        type=lambda x: is_new_file(parser,x) )


    args = parser.parse_args()

    fname = args.alignment
    fasta_id = args.target
    pdb_file = args.pdb
    out_file = args.out

    data = PDBProcessor.read_fasta_alignment( fname )
    if fasta_id not in data:
        raise ValueError('Error: missing key in the file!')
    keys = sorted( list(data.keys()) )

    bit_score = PDBProcessor.calculate_amino_alignment_bitscore( [x  for x in data.values()], conservation=True, percentage=True )
    proteins = PDBProcessor.read_pdb(pdb_file)
    if len(proteins) > 1:
        raise ValueError('Error: multiple protein in the file!')

    if str(data[fasta_id]).replace(PDBProcessor.AminoAcid.gap_character(), "") != proteins[0].amino_sequence():
        raise ValueError("Error: target alignement sequence differs from pdb amino sequence.")

    gapped_sequence = data[fasta_id]
    all_protein_residues = PDBProcessor.get_all_residue_pos_in_alignment(gapped_sequence)
    residue_ids = proteins[0].residue_ids()
    assert len(residue_ids) == len(all_protein_residues)
    for i in range(len(residue_ids)):
        resId = residue_ids[i]
        alnPos = all_protein_residues[i]
        score  = bit_score[alnPos]
        proteins[0].set_residue_bfactor(resId, score)

    with open(out_file, 'w+') as f:
        f.write( str(proteins[0]) )

    exit(0)
