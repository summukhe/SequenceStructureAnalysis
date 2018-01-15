import os
import PDBProcessor


class ProteinStructure:
    def __init__(self, name, coordinates, resIds, resNames, atomIds, atomNames, bfactors = None):
        assert isinstance(coordinates, list)
        assert isinstance(resIds, list)
        assert isinstance(resNames, list)
        assert isinstance(atomIds, list)
        assert isinstance(atomNames, list)
        n = len(coordinates)
        if bfactors is None:
            bfactors = [ 0.0 for i in range(n) ]
        assert isinstance(bfactors, list)
        assert n == len(resIds)
        assert n == len(resNames)
        assert n == len(atomIds)
        assert n == len(atomNames)
        assert n == len(bfactors)

        self._name = name
        self._coordinates = coordinates
        self._atom_ids = atomIds
        self._atom_names = atomNames
        self._bfactors = bfactors
        self.__build_protein( resIds, resNames )

    def __build_protein( self, resIds, resNames ):
        self.__residues = dict()
        n = len(resIds)

        for i, resId in enumerate(resIds):
            if resId in self.__residues:
                self.__residues[resId]['atoms'].append(i)
                assert self.__residues[resId]['resname'] == resNames[i]
            else:
                self.__residues[resId] = {'resname' : resNames[i], 'atoms' : [i] }

    def chain(self):
        return self._name

    def residue_ids(self):
        return sorted( [int(i) for i in list(self.__residues.keys())])

    def amino_sequence(self):
        ss = ''
        for r in sorted( [ int(i) for i in list(self.__residues.keys()) ] ):
            ss += PDBProcessor.AminoAcid.to_one_letter_code(self.__residues[str(r)]['resname'] )
        return ss

    def amino_frequency(self):
        freq = dict()
        for r in list(self.__residues.keys()):
            aa = self.__residues[r]['resname']
            if aa not in freq:
                freq[aa] = 0
            freq[aa] += 1
        return freq

    def set_atom_bfactor(self, atomId, bfactor):
        assert isinstance(bfactor, float) and bfactor <= 100.
        atomId = str(atomId)
        try:
            idx = self._atom_ids.index(atomId)
            self._bfactors[idx] = bfactor
        except:
            pass

    def set_residue_bfactor(self, resId, bfactor):
        assert isinstance(bfactor, float) and bfactor <= 100.
        resId = str(resId)
        if resId in self.__residues:
            for i in self.__residues[resId]['atoms']:
                self._bfactors[i] = bfactor


    def get_residue(self, resId):
        resId = str(resId)
        if resId not in self.__residues:
            raise ValueError('Error: invalid residue id requested <%s>' % resId)
        seq = self.__residues[resId]['atoms']
        resname = self.__residues[resId]['resname']
        return ProteinStructure( self._name,
                                 [self._coordinates[i] for i in seq],
                                 [resId for i in seq],
                                 [resname for i in seq],
                                 [self._atom_ids[i] for i in seq],
                                 [self._atom_names[i] for i in seq],
                                 [self._bfactors[i] for i in seq])

    def get_residues(self, resIdList):
        assert isinstance(resIdList, list)
        resIdList = [ str(resId) for resId in resIdList ]
        seq, resname = list(), list()
        for resId in resIdList:
            if resId in self.__residues:
                seq = seq + self.__residues[resId]['atoms'][:]
                resname = resname + [ self.__residues[resId]['resname'] for i in range( len(self.__residues[resId]['atoms']) )]
        if len(seq) == 0:
            raise ValueError('Error: invalid residue ids requested')
        return ProteinStructure(self._name,
                                [self._coordinates[i] for i in seq],
                                [resId for i in seq],
                                resname,
                                [self._atom_ids[i] for i in seq],
                                [self._atom_names[i] for i in seq],
                                [self._bfactors[i] for i in seq])


    def __len__(self):
        return len(self.__residues)

    def __str__(self):
        output = ''
        chain = self._name
        resIds = sorted([ int(i) for i in list( self.__residues.keys())] )
        for r in resIds:
            resname = self.__residues[str(r)]['resname']
            for i in self.__residues[str(r)]['atoms']:
                output = output + 'ATOM  %5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f \n' % (int(self._atom_ids[i]),self._atom_names[i],resname,chain,int(r),float(self._coordinates[i][0]),float(self._coordinates[i][1]), float(self._coordinates[i][2]),1.0, float(self._bfactors[i]) )
        return output




def pdb_atomline_parser( pdbline ):
    assert pdbline.startswith('ATOM  ')
    chainId = pdbline[21].strip()
    resName = pdbline[17:20].strip()
    resId   = pdbline[22:26].strip()
    atomName = pdbline[12:16].strip()
    atomId   = pdbline[6:11].strip()
    x,y,z    = pdbline[30:38].strip(), pdbline[38:46].strip(), pdbline[46:54].strip()
    tempF    = pdbline[60:66].strip()
    return { 'chain' : chainId,
             'resid' : resId,
             'resname' : resName,
             'atomid'  : atomId,
             'atomname' : atomName,
             'coord' : (x, y, z),
             'bfactor' : tempF }


def read_pdb(pdbfile):
    assert os.path.exists(pdbfile)
    chains = dict()
    with open(pdbfile, 'r') as infile:
        for line in infile:
            if line.startswith('ATOM  '):
                fields = pdb_atomline_parser(line)
                if fields['chain'] not in chains:
                    chains[fields['chain']] = { 'atomnames': list() , 'atomids' : list(),
                                                'coordinates': list(), 'resnames': list(),
                                                'resids' : list(), 'bfactors' : list() }

                chains[fields['chain']]['atomnames'].append(fields['atomname'])
                chains[fields['chain']]['atomids'].append(fields['atomid'])
                chains[fields['chain']]['coordinates'].append( fields['coord'] )
                chains[fields['chain']]['resnames'].append(fields['resname'])
                chains[fields['chain']]['resids'].append(fields['resid'])
                chains[fields['chain']]['bfactors'].append(fields['bfactor'])

    proteins = list()
    for c in chains:
        proteins.append( ProteinStructure(c,
                                          chains[c]['coordinates'],
                                          chains[c]['resids'],
                                          chains[c]['resnames'],
                                          chains[c]['atomids'],
                                          chains[c]['atomnames'],
                                          chains[c]['bfactors']))
    return proteins


def read_residue_list( residueFile ):
    residues = list()
    with open(residueFile, "r") as fInput:
        for line in fInput:
            if line.strip() == "" or line.strip().startswith('#'):
                continue
            flds = line.strip().split()
            if len(flds) == 2 and PDBProcessor.AminoAcid.check_amino(flds[0]):
                residues.append( (int(flds[1]), flds[0]) )
    return sorted(residues)

