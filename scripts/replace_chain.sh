#! /bin/bash

if [ $# -lt 3 ] || [ ! -f "$1" ] || [ $(echo -n "$2" | wc -c) -ne 3 ]
then
  echo "Usage: $0 <PDB file> <Replacing ligand name> <chain id list>" 1>&2;
  exit 1;
fi

pdb=$(readlink -f "$1")

shift

ligTag=$1

shift

chains=""

while [ $# -ne 0 ]
do
  chains=${chains}$1
  shift
done

awk -vchainId="[${chains}]" -vligTag="${ligTag}"  '/^ATOM  / || /^TER /{ 
                              if(substr($0,22,1) ~ chainId ){ 
                                    printf( "%6s%s%s%s%4d%s\n", "HETATM",substr($0,7,11),ligTag,substr($0,21,2),1,substr($0,27) ); 
                              }else{ 
                                    print($0); 
                              } } ' ${pdb}
