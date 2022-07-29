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
USERNAME=ciro.suastegui
DIRNAME=/Volumes/mdrive
REMNAME=${USERNAME}@vicb-submit-02:/home/icb/${USERNAME}

ls -hola /Volumes
# if [[ ! -z "$(ls -A ${DIRNAME})" ]]; then
if [[ -d ${DIRNAME} ]]; then
  echo "\033[0;31mMounted\033[0m, or dir exists, would you like to unmount (y/n)?"
  read -p "" choice
  case "$choice" in
    y|Y ) sudo umount ${DIRNAME};;
    n|N ) echo "\033[0;34mCool\033[0m";;
    * ) echo "invalid";;
  esac
else
  sudo mkdir -p ${DIRNAME}
  echo "\033[0;32mWrite your remote password\033[0m"
  sudo sshfs -o allow_other,default_permissions ${REMNAME} ${DIRNAME}
fi
if [[ -d ${DIRNAME} ]]; then
  echo "\033[0;32mContent:\033[0m"
  ls -hol ${DIRNAME}/;
fi

# If the mounted unit becomes non-responsive:
# `ps -ef` and `pkill <sshfs process id>`
