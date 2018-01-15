#! /bin/bash


struct_aln_exe=mustang

if [ "$(which $struct_aln_exe)" != "" ];
then
     struct_aln_exe="$(which $struct_aln_exe)"
else
     echo "Fatal Error: Executable ($struct_aln_exe) not found !" 1>&2;
     exit 1;
fi

if [ $# -ne 2 ]
then
  echo "Usage: $0 <Input PDB directory> <output directory>" 1>&2;
  exit 1;
fi


input_dir=$(readlink -f "$1")
out_dir=$(readlink -f "$2")


if [ ! -d "${input_dir}" ] || [ $(ls ${input_dir}/*.pdb | wc -l) -lt 0 ] || [ $(bash check_chains.sh "${input_dir}" | awk 'NF != 2' | wc -l) -ne 0 ] 
then
   echo "Error: error in the input directory ${input_dir}" 1>&2;
   exit 1;
fi

if [ ! -e "${out_dir}" ]
then
  mkdir -p ${out_dir}
fi

uniq_id=mustang_$$

curr_dir=${PWD}

cd "${out_dir}"

${struct_aln_exe} -p "${input_dir}/"  -i $(ls ${input_dir}/*.pdb | xargs -i basename {}) -F fasta -o  ${uniq_id} -s OFF -r ON  2>&1

cd "${curr_dir}"

exit 0;
