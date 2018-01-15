#! /bin/bash

exedir=$(dirname $(readlink -f $0))
find_nbr=${exedir}/find_neighbor_pdb


if [ ! -x ${find_nbr} ]
then
  echo "Error: missing executable ${find_nbr}" 1>&2
  exit 1
fi


if [ $# -ne 1 ] || [ ! -d $1 ]
then
  echo "Usage: $0 <target directory>" 1>&2
  exit 1
fi


target_dir=$(readlink -f "$1")

ls ${target_dir}/*.pdb | grep -P "\.pdb$" | while read pdb; 
   do 
      echo $pdb $(find_neighbor_pdb -f $pdb -C 2>/dev/null) ; 
   done
