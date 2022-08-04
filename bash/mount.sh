#!/usr/bin/bash

# ------------------------------------------------------------------------------
# title: Mounting.
# purpose: This script helps you mount/unmount.
# author:
#   - name: Ciro Ramírez-Suástegui
#     email: cramsuig@gmail.com
# date: 2022-01-12
# ------------------------------------------------------------------------------

# If Mac: http://osxfuse.github.io install FUSE and SSHFS
if [[ ! -s $(which sshfs) ]]; then
  echo "sudo apt-get install sshfs"
else
  echo "SSHFS present: $(which sshfs)"
fi

echo "Might need to write your local password"
USERNAME=${1:-ciro.suastegui}
HNAME=${2:-hpc-submit01}
DIR_REMOTE=${3:-/lustre/home/icb/${USERNAME}}
DIR_LOCAL=${4:-/Volumes/drive_hzm}
REMNAME=${USERNAME}@${HNAME}:${DIR_REMOTE}

echo "User: ${USERNAME}"
echo "Host: ${HNAME}"
echo "Remote: ${DIR_REMOTE}"
echo "Local: ${DIR_LOCAL}"

ls -hola /Volumes
# if [[ ! -z "$(ls -A ${DIR_LOCAL})" ]]; then
if [[ -d ${DIR_LOCAL} ]]; then
  echo "\033[0;31mMounted\033[0m, or dir exists, would you like to unmount (y/n)?"
  read -p "" choice
  case "$choice" in
    y|Y ) sudo umount ${DIR_LOCAL};;
    n|N ) echo "\033[0;34mCool\033[0m";;
    * ) echo "invalid";;
  esac
else
  sudo mkdir -p ${DIR_LOCAL}
  echo "\033[0;32mWrite your remote password\033[0m"
  sudo sshfs -o allow_other,default_permissions,umask=0022 ${REMNAME} ${DIR_LOCAL}
fi
if [[ -d ${DIR_LOCAL} ]]; then
  echo "\033[0;32mContent:\033[0m"
  ls -hol ${DIR_LOCAL}/;
fi

# If the mounted unit becomes non-responsive:
# `ps -ef` and `pkill <sshfs process id>`
