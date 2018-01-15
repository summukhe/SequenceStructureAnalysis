#! /bin/bash

PYTHON=/usr/bin/python3.5
consv_code=/home/sumanta/PycharmProjects/RLECode/MarkSequenceConservation.py

exe_dir=$(dirname $(readlink -f $0))
msa_exe=${exe_dir}/perform_msa.sh

if [ ! -f ${msa_exe} ] || [ ! -f ${consv_code} ]
then
  echo "Error: missing dependent executable " 1>&2
  exit 1;
fi

if [ $# -ne 2 ]
then
  echo "Usage: $0 <pdb directory> <target pdb>" 1>&2;
  exit 1;
fi


if [ ! -d "$1" ] || [ $(ls "$1"/*.pdb | wc -l) -lt 2 ]
then
  echo "Error: not sufficient pdb files in the input directory" 1>&2;
  exit 1;
fi

input_dir=$(readlink -f "$1")
pdb_tag=$2

if [ ! -e "${input_dir}/${pdb_tag}.pdb" ]
then
  echo "Error: the target directory is not present in the input folder" 1>&2;
  exit 1;
fi

out_dir="test_$$"

bash ${msa_exe} "${input_dir}" "${out_dir}"

status_code=$?

if [ ${status_code} -ne 0 ]
then
  echo "Error: improper status code ${status_code}" 1>&2;
  if [ ! -d "${out_dir}" ]
  then
     rm -fr "${out_dir}"
  fi
  exit 1;
fi

alignment_file=$(find ${out_dir} -iname "mustang_*.afasta" -type f)

if [ ! -f "${alignment_file}" ]
then
  echo "Error: can not find sequence alignment file" 1>&2;
  exit 1;
fi

$PYTHON ${consv_code} --alignment ${alignment_file} --seq-tag ${pdb_tag} --pdb "${input_dir}/${pdb_tag}.pdb" --out ${out_dir}/${pdb_tag}_consv.pdb
