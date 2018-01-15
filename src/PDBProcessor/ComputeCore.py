import numpy as np
import PDBProcessor


def calculate_amino_alignment_bitscore( alignments , conservation = False, percentage = False ):
    assert isinstance(alignments, list)
    nseq = len(alignments)
    nlen = None
    for aln_seq in alignments:
        assert isinstance(aln_seq, PDBProcessor.GappedAminoSequence)
        if nlen is None:
            nlen = aln_seq.size()
        assert nlen == aln_seq.size()
    bit_score = [0.0 for i in range(nlen)]
    for i in range(nlen):
        amino_freq = {aa: 0. for aa in PDBProcessor.AminoAcid.all_aminos()}
        for s in alignments:
            aa = s[i+1]
            if aa == '-':
                for a in amino_freq.keys():
                     amino_freq[a] += 1./len(PDBProcessor.AminoAcid.all_aminos())
            else:
                amino_freq[aa] += 1.
        prob = np.array([f for f in amino_freq.values() if f > 0])/nseq
        bit_score[i] = -1. * np.sum( prob * np.log2(prob) )
    if conservation is True:
        max_score =  np.log2( len(PDBProcessor.AminoAcid.all_aminos()) )
        bit_score =  max_score - bit_score
        if percentage is True:
            bit_score = (100.0 / max_score) * bit_score
    return bit_score.tolist()


def find_exact_protein_in_alignment( aligments, target_sequence):
    assert isinstance(aligments, list)
    assert isinstance(target_sequence, PDBProcessor.AminoSequence)
    aln_length = None
    seq = str(target_sequence)
    matches = list()
    for i, aln in enumerate(aligments):
        assert isinstance(aln, PDBProcessor.GappedAminoSequence)
        if aln_length is None:
            aln_length = aln.size()
        assert aln_length == aln.size()
        if str(aln).replace("-", "") == seq:
            matches.append(i)
    return matches


def find_ith_residue_from_alignment( alignment, position):
    assert isinstance( alignment, PDBProcessor.GappedAminoSequence )
    sequence = str(alignment)
    pos = 0
    for i, s in enumerate(sequence):
        if s != '-':
            pos += 1
        if pos == position:
            return i
    return None


def get_all_residue_pos_in_alignment( gapped_sequence ):
    assert isinstance( gapped_sequence, PDBProcessor.GappedAminoSequence )
    sequence = str(gapped_sequence)
    pos = list()
    for i, s in enumerate(sequence):
        if s != '-':
            pos.append(i)
    return pos