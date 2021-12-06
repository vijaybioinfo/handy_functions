#!/bin/bash

# set -euo pipefail

##################################
# Create new project directories #
##################################

# This script creates the necessary folders for a new project
# using the 1st template from https://github.com/cramirezs/project_structure

PROJECT_NAME="${1}"
HOST=${USER}
SMALL_DISC=${HOME} # ${3}
if [[ "${4}" == "" ]]; then
  if [[ -d "/mnt/BioAdHoc/Groups/vd-vijay" ]]; then
    BIG_DISC=/mnt/BioAdHoc/Groups/vd-vijay
  else
    BIG_DISC=${HOME}
  fi
fi

echo ' '
echo -e "\033[0;32m--- Configuration ---\033[0m"
echo "Name: ${PROJECT_NAME}"
echo "User: ${HOST}"
echo "Small disc: ${SMALL_DISC}"
echo "Big disc: ${BIG_DISC}"

if [[ "${PROJECT_NAME}" == "" ]]; then echo "No project name provided"; exit; fi
if [[ -d "${PROJECT_NAME}" ]]; then
  echo -e "\033[0;31mProject name exists as a folder in current location\n `pwd`\033[0m";
  echo ' '
  exit;
fi

while true; do
    read -p "Do you wish to continue with this structure? " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

echo -e "\033[0;32m--- Creating directories ---\033[0m"

mkdir ${SMALL_DISC}/${PROJECT_NAME}
mkdir ${BIG_DISC}/${HOST}/${PROJECT_NAME}
mkdir ${SMALL_DISC}/${PROJECT_NAME}/{info,code}
mkdir ${BIG_DISC}/${HOST}/${PROJECT_NAME}/{data,results,figures}

if [[ "${SMALL_DISC}" != "${BIG_DISC}" ]]; then
  echo -e "\033[0;36mCreating symbolic links\033[0m"
  for DNAME in $(echo {data,results,figures}); do
    ln -s ${BIG_DISC}/${HOST}/${PROJECT_NAME}/${DNAME} ${SMALL_DISC}/${PROJECT_NAME}/${DNAME}
  done
  for DNAME in $(echo {info,code}); do
    ln -s ${SMALL_DISC}/${PROJECT_NAME}/${DNAME} ${BIG_DISC}/${HOST}/${PROJECT_NAME}/${DNAME}
  done
fi

echo ' '
echo -e "\033[0;32m--- Directories contain: ---\033[0m"
echo ${SMALL_DISC}/${PROJECT_NAME}/
ls -loh --color=auto ${SMALL_DISC}/${PROJECT_NAME}/
echo ' '
echo ${BIG_DISC}/${HOST}/${PROJECT_NAME}/
ls -loh --color=auto ${BIG_DISC}/${HOST}/${PROJECT_NAME}/

echo -e "\033[0;32m--- ######## ########### ---\033[0m"
echo ' '
