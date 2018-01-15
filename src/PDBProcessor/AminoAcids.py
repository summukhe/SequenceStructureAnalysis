

class AminoAcid:

    @staticmethod
    def all_aminos():
        return 'ACDEFGHIKLMNPQRSTVWY'

    @staticmethod
    def  to_three_letter_code( chr ):
        lookup = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP',
                  'E': 'GLU', 'F': 'PHE', 'G': 'GLY',
                  'H': 'HIS', 'I': 'ILE', 'K': 'LYS',
                  'L': 'LEU', 'M': 'MET', 'N': 'ASN',
                  'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
                  'S': 'SER', 'T': 'THR', 'V': 'VAL',
                  'W': 'TRP', 'Y': 'TYR'}
        if chr in lookup.keys():
            return lookup[chr]
        else:
            raise ValueError('Error: invalid amino code %s' % chr)


    @staticmethod
    def to_one_letter_code( chr ):
        lookup = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP',
                  'E': 'GLU', 'F': 'PHE', 'G': 'GLY',
                  'H': 'HIS', 'I': 'ILE', 'K': 'LYS',
                  'L': 'LEU', 'M': 'MET', 'N': 'ASN',
                  'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
                  'S': 'SER', 'T': 'THR', 'V': 'VAL',
                  'W': 'TRP', 'Y': 'TYR'}
        rlookup = dict()
        for k,v in lookup.items():
            rlookup[v] = k
        if chr in rlookup.keys():
            return rlookup[chr]
        else:
            raise ValueError('Error: invalid amino name %s' % chr)

    @staticmethod
    def gap_character():
        return '-'

    @staticmethod
    def check_amino(chr):
        if isinstance(chr, str) is False:
            return False
        if len(chr) == 1:
            try:
                AminoAcid.to_three_letter_code(chr.upper())
                return True
            except:
                pass
        elif len(chr) == 3:
            try:
                AminoAcid.to_one_letter_code(chr.upper())
                return True
            except:
                pass
        return False