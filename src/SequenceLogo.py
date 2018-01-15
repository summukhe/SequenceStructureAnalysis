import os
import argparse
import PDBProcessor


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Utility reads a sequence alignment, " +
                                                 "for corresponding sequence tag, " +
                                                 "for a given PDB file and specific residues " +
                                                 "it reads the sequence conseveration from the alignment file" )

    parser.add_argument('--alignment', dest="alignment", metavar="FILE", required=True,
                        help="sequence alignment file input in fasta format",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument('--seq-tag', dest='target', metavar='ID', required=True,
                        help="ID of the sequence in the alignment",
                        type=str)

    parser.add_argument('--pdb', dest='pdb', metavar='FILE', required=True,
                        help="ID full / relative path to the pdb structure",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument('--residues', dest='residueFile', metavar='FILE', required=True,
                        help="Output pdb file name",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument('--out', dest='out', metavar='FILE', required=True,
                        help="Output fasta file name, where all the sequence logo will be written.",
                        type=lambda x: is_new_file(parser, x))

    args = parser.parse_args()

    alnFile = args.alignment
    targetId = args.target
    pdbFile = args.pdb
    residueFile = args.residueFile
    outFasta = args.out

    data = PDBProcessor.read_fasta_alignment(alnFile)

    if targetId not in data:
        raise ValueError('Error: missing sequence id in the alignment file [%s]' % targetId )

    proteins = PDBProcessor.read_pdb(pdbFile)
    if len(proteins) > 1:
        raise ValueError('Error: PDB contains more than one chain [%d]' % len(proteins))

    if str(data[targetId]).replace(PDBProcessor.AminoAcid.gap_character(), '') != proteins[0].amino_sequence():
        raise ValueError('Error: alignment sequence does not match PDB sequence' )

    target_residues = PDBProcessor.read_residue_list( residueFile )

    gapped_sequence = data[targetId]
    all_protein_residues = PDBProcessor.get_all_residue_pos_in_alignment(gapped_sequence)
    residue_ids = proteins[0].residue_ids()
    target_alnpos = [all_protein_residues[residue_ids.index(r)] for r,t in target_residues]

    sequence_signature = dict()
    for pdb_id , aln_seq in data.items():
        sequence_signature[pdb_id] = "".join( [ aln_seq[i] for i in target_alnpos] )

    #outFasta = "/home/sumanta/Project/structural_analysis/sunaina/pdbs/5MFX_rna_contacts/rna_contacts.fasta"

    with open(outFasta, "w+") as f:
        for k, s in sequence_signature.items():
            f.write(">%s\n%s\n" % (k,s))


    #nfile = "/home/sumanta/Project/structural_analysis/sunaina/pdbs/5MFX_consv.pdb"
    #with open(nfile, "w+") as f:
    #    f.write( str(proteins[0]) )

    """
    residues = proteins[0].residue_ids()
    for i in range(nsize):
        idx = PDBProcessor.find_ith_residue_from_alignment(gapped_sequence, i+1)
        proteins[0].set_residue_bfactor( residues[i], bit_score[idx] )
    with open('/home/sumanta/Project/structural_analysis/sunaina/pdbs/res_185.pdb', 'w') as f:
        f.write(str(proteins[0]))
    """

