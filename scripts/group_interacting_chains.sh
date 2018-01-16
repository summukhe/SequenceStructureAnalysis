#! /bin/bash

exedir=$(dirname $(readlink -f $0))
exeName="${exedir}/replace_chain.sh"

if [ ! -f "${exeName}" ]
then
  echo "Error: can not find the dependent executable [${exeName} ]" 1>&2;
  exit 1;
fi


findNbr=find_neighbor_pdb
findNbrExe=$(which $findNbr)


if [ "${findNbrExe}" == "" ] || [ ! -x "${findNbrExe}" ] 
then
  echo "Error: can not find dependent executable [$findNbr]" 1>&2;
  exit 1;
fi


if [ $# -lt 1 ] || [ ! -f "$1" ] || [ "$(echo $1 | grep -P '.pdb$')" != "$1" ]
then
  echo "
=========================================================================================
Usage: $0 <pdb file name> [chain length cutoff] [three letter ligand] [output directory]
=========================================================================================
default chain length cutoff is 15
default ligand code is RNA 
default output directory is same as the base pdb directory
" 1>&2;
  exit 1;
fi

pdb="$1"
cutoff_len=15
ligName="RNA"

pdb=$(readlink -f "${pdb}")
dirName=$(dirname "${pdb}")
baseName=$(basename "${pdb}")
outdir="${dirName}"

if [ $# -gt 1 ] && [ "$(echo -n $2 | grep -P '^\d+$')" -eq "$2" ]
then
  cutoff_len=$2
fi

if [ $# -gt 2 ] && [ "$(echo -n $3 | wc -c)" -eq 3 ]
then
  ligName=$3
fi

if [ $# -gt 3 ] 
then
  outdir=$(readlink -f "$4")
fi

if [ ! -f "${pdb}" ]
then
  echo "Error: can not find the PDB input [${pdb}] " 1>&2;
  exit 1;
fi

if [ ! -e "${outdir}" ]
then
  mkdir -p "${outdir}"
elif [ ! -d "${outdir}" ]
then
  echo "Error: not a valid directory [${outdir}], execution halted" 1>&2;
  exit 1;
fi

targetPDB="/tmp/temporaryPDB_$$.pdb"

ligChains=$(${findNbrExe} -f "${pdb}" -C 2>/dev/null | awk -F: -vcutoff=${cutoff_len} '$2 < cutoff { printf("%s\n",$1)}')

${SHELL} "${exeName}" "$pdb" "${ligName}" ${ligChains} > $targetPDB

if [ ! -f $targetPDB ] || [ "$(echo ${targetPDB} | grep -P '\.pdb$')" != "${targetPDB}" ]
then
   echo "Error: can not find the target pdb, or the file name is not ending with .pdb extension" 1>&2;
   exit 1;
fi



${findNbrExe} -f "${targetPDB}" -l | while read lig;
do
  ligChain=$(echo $lig | awk -F: '{print $2}')
  protChain=$(${findNbrExe} -f "${targetPDB}" -n "${lig}" -r 4.5 -t | 
              awk '{print substr($0,22,1)}' | uniq -c | 
              sort -k1,1nr | head -1 | awk '{print $2}');
  
  chains="${ligChain}${protChain}"
  awk -vchains="[${chains}]" 'substr($0,22,1) ~ chains' "${targetPDB}" > "${outdir}/${baseName%%.pdb}_${protChain}_${ligChain}.pdb"
done

rm -f "${targetPDB}"
