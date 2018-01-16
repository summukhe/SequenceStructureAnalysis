#! /bin/bash

exedir=$(dirname $(readlink -f $0))
intExe=${exedir}/group_interacting_chains.sh
msaExe=${exedir}/perform_msa.sh

for exe in "${intExe}" "${msaExe}" 
do
  if [ ! -f "${exe}" ]
  then
    echo "Error: dependent script is missing" 1>&2;
    exit 1;
  fi
done

findNbr="find_neighbor_pdb"
findNbrExe=$(which ${findNbr})

if [ "${findNbrExe}" == "" ]  || [ ! -x "${findNbrExe}" ]
then
  echo "Error: missing installation of ${findNbr}" 1>&2;
  exit 1;
fi

targetDir="ns3_rna"
cutoff_length=15
ligName="RNA"
minContacts=2

targetDir="$(readlink -f "${targetDir}")"

if [ ! -d "${targetDir}" ]
then
  echo "Error: directory not found [${targetDir}]" 1>&2;
  exit 1;
fi

if [ $(ls ${targetDir}/*.pdb | wc -l) -lt 1 ]
then
  echo "Error: directory [${targetDir}] has no pdb file" 1>&2;
  exit 1;
fi

out_pdbdir="${targetDir}/interacting_chains_$$"
out_contactdir="${targetDir}/interacting_residues_$$"
out_pdbChain="${targetDir}/interacting_protein_$$"
out_alignment="${targetDir}/interacting_alignment_$$"

for pdb in ${targetDir}/*.pdb;
do
  ${SHELL} ${intExe} "${pdb}" ${cutoff_length} ${ligName} "${out_pdbdir}"
done

if [ ! -d "${out_contactdir}" ]
then
  mkdir -p "${out_contactdir}"
fi

if [ ! -d "${out_pdbChain}" ]
then
  mkdir -p "${out_pdbChain}"
fi

if [ ! -d "${out_alignment}" ]
then
  mkdir -p "${out_alignment}"
fi

for pdb in ${out_pdbdir}/*.pdb; 
do 
  pdbName=$(basename "${pdb}")
  chainId=$(echo ${pdbName} | awk -F_ '{print $(NF-1)}')
  ${findNbrExe} -f "${pdb}" -r 4.5 | awk '{print substr($0,18,3),substr($0,23,4)}' | 
                                   sort | uniq -c | 
                                   awk -vmin_contacts=${minContacts} '$1 > min_contacts {
                                                                        printf("%s %d\n",$2,$3);
                                                                      }' > ${out_contactdir}/${pdbName%%.pdb}.contacts;
  ${findNbrExe} -f "${pdb}" -c ${chainId} -e > ${out_pdbChain}/${pdbName}
done

${SHELL} ${msaExe} ${out_pdbChain} ${out_alignment}
