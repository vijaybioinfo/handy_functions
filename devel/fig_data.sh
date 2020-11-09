#!/usr/bin/bash

#########################
# Figure data gathering #
#########################

# This script will place all the necessary files for a figures folder

# OUTDIR
# FIGDATA
UNLINK1ST=${UNLINK1ST:-FALSE}
INFNAMES=${INFNAMES:-NONE} # array of scripts

echo
echo -e "\033[0;36m*** Linking  files - Vijay Lab\033[0m"
echo "Destination folder: ${OUTDIR}"
echo "Files links and names: ${FIGDATA}"
echo "Unlink first: ${UNLINK1ST}"

if [[ ! -d ${OUTDIR} ]]; then mkdir ${OUTDIR}; fi
cd ${OUTDIR}; echo "Jumped to `pwd`"
echo -e "\033[0;31m=========== Existing files ============\033[0m"
echo "Directories"; ls
echo -e "\033[0;34m===========\033[0m data folder"
if [[ ! -d data ]]; then mkdir data; fi
ls -loh data

## Create unique names for your files
FILE0=(`cut -f1 -d, ${FIGDATA}`)
FILE1=(`cut -f2 -d, ${FIGDATA}`)
if [[ $((`wc -l ${FIGDATA} | sed 's/ .*//g'`)) -gt 1 ]]; then FN=${#FILE0[@]}; else FN=1; fi

# If you need to unlink
if [[ "${UNLINK1ST}" == "TRUE" ]]; then
  echo "Unlinking files"
  for FITER in `seq 0 $((${FN} - 1))`; do
    unlink data/${FILE1[FITER]}
  done
fi

echo "Creating ${FN} link(s)"
for FITER in `seq 0 $((${FN} - 1))`; do
  ln -s ${FILE0[FITER]} data/${FILE1[FITER]}
done
echo -e "\033[0;32m=========== Newly  created ============\033[0m"
ls -lho data

# Modify files with these links
if [[ "${INFNAMES[0]}" != "NONE" ]]; then
  echo "Modifying ${#INFNAMES[@]} file(s)"
  for FNAME in ${INFNAMES[@]}; do
    echo "In ${FNAME}"
    for FITER in `seq 0 $((${FN} - 1))`; do
      printf "."
      sed -i 's/'${FILE0[FITER]//\//\\/}'/..\/data\/'${FILE1[FITER]}'/g' ${FNAME}
    done
    echo
  done
fi

cd -
