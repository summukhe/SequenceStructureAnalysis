import re


class Sequence:
    def __init__(self, seq , pattern = '^[A-Z]+$'):
        assert isinstance(seq, str)
        if re.match( pattern , seq.upper() ) :
            self._sequence = seq.upper()
            self._pattern = pattern
        else:
            raise ValueError('Error: invalid sequence')

    def size(self):
        return len( self._sequence )

    def __getitem__(self, i):
        if i > 0 and i <= self.size():
            return self._sequence[i-1]
        else:
            raise IndexError('Error: invalid sequence position %d' % i )

    def __str__(self):
        return self._sequence



class MutableSequence(Sequence):

    def __setitem__(self, i, c):
        if re.match(self._pattern, c.upper()) and len(c) == 1:
            if i > 0 and i <= self.size():
                self._sequence = self._sequence[:(i-1)] + c.upper() + self._sequence[i:]
            else:
                raise IndexError('Error: invalid sequence position %d' % i)
        else:
            raise ValueError('Error: invalid insert character %s' % c)

    def mutate(self, i , c):
        return self.__setitem__(i, c)

    def insert(self, i, c):
        if re.match(self._pattern, c.upper()):
            if i > 0 and i <= self.size():
                self._sequence = self._sequence[:(i-1)] + c.upper() + self._sequence[(i-1):]
            else:
                raise IndexError('Error: invalid sequence position %d' % i)
        else:
            raise ValueError('Error: invalid insert character %s' % c)

    def delete(self, i):
        if i > 0 and i <= self.size():
            self._sequence = self._sequence[:(i-1)] + self._sequence[i:]



class AminoSequence(Sequence):
    def __init__(self, seq):
        super(AminoSequence, self).__init__(seq, pattern='^[ACDEFGHIKLMNPQRSTVWY]+$')


class DNASequence(Sequence):
    def __init__(self, seq):
        super(DNASequence, self).__init__(seq, pattern='^[ACGT]+$')


class ProteinSequence(MutableSequence):
    def __init__(self, seq):
        super(ProteinSequence, self).__init__(seq, pattern='^[ACDEFGHIKLMNPQRSTVWY]+$')

class GeneSequence(MutableSequence):
    def __init__(self, seq):
        super(GeneSequence, self).__init__(seq, pattern='^[ACGT]+$')


class GappedAminoSequence(Sequence):
    def __init__(self, seq):
        super(GappedAminoSequence, self).__init__(seq, pattern='^[ACDEFGHIKLMNPQRSTVWY-]+$')


class GappedDNASequence(Sequence):
    def __init__(self, seq):
        super(GappedDNASequence, self).__init__(seq, pattern='^[ACGT-]+$')