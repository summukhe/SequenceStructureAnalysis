import os
import re
import PDBProcessor


def read_fasta_file(fname, sequence_pattern = '^[A-Z]+$'):
    if not os.path.exists(fname):
        raise FileExistsError('Error: can not find file %s' % fname )

    header, sequence = None, None
    sequences = dict()
    with open(fname, 'r') as infile:
        for line in infile:
            line = line.strip()
            if re.match('^>(\w+)', line):
                smatch = re.match('^>(\w+)', line)
                if header is not None and sequence is not None:
                    sequences[header] = sequence
                header = smatch.group(1)
                sequence = None
            elif header is not None and re.match(sequence_pattern, line.upper()):
                if sequence is None:
                    sequence = line.upper()
                else:
                    sequence = sequence + line.upper()

    if sequence is not None:
        sequences[header] = sequence
    return sequences


def read_fasta_protein( fname ) :
    sequence_pattern = '^[ACDEFGHIKLMNPQRSTVWY]+$'
    data = read_fasta_file( fname, sequence_pattern )
    res = dict()
    for hdr, seq in data.items():
        res[hdr] = PDBProcessor.AminoSequence(seq)
    return res


def read_fasta_alignment( fname ):
    sequence_pattern = '^[ACDEFGHIKLMNPQRSTVWY-]+$'
    data = read_fasta_file(fname, sequence_pattern)
    res = dict()
    for hdr, seq in data.items():
        res[hdr] = PDBProcessor.GappedAminoSequence(seq)
    return res


def read_fasta_dna( fname ):
    sequence_pattern = '^[ACGT]+$'
    data = read_fasta_file(fname, sequence_pattern)
    res = dict()
    for hdr, seq in data.items():
        res[hdr] = PDBProcessor.DNASequence(seq)
    return res