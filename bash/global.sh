#!/bin/bash

# This script contains functions used for bash-based programmes

check_dir_content () {
  clea=`ls $1 | wc -l`
  if [ ${clea} -gt 0 ]; then
    DELIM="=========== $1: ==========="
    echo -e "\033[0;31m${DELIM}\033[0m"
    echo -e "Contains:\n`ls $1`"
    printf "\033[0;31m%0.s=\033[0m" $(seq 1 ${#DELIM}); echo
    while true; do
      read -p "Do you wish to erase previous information and continue with the analysis? (y/n/s): " yn
      case $yn in
        [Yy]* ) read -p "Are you completely sure? " yn
          case $yn in
            [Yy]* ) echo "Cleaning..."; rm -r ${1}/*; break;;
            [Nn]* ) exit;;
            [Ss]* ) break;;
            * ) echo "Please answer yes (Y/y) or no (N/n) or skip (S/s).";;
          esac
          ;;
        [Nn]* ) echo "Aborting..."; exit;;
        [Ss]* ) break;;
        * ) echo "Please answer yes (Y/y) or no (N/n) or skip (S/s).";;
      esac
    done
  fi
}

# ## Enclose words in delimiters
# function delim_words () {
#   DELIM="=========== $1: ==========="
#   echo -e "\033[0;31m${DELIM}\033[0m"
#   printf "\033[0;31m%0.s=\033[0m" $(seq 1 ${#DELIM}); echo
# }

# Seconds converter
function convertsecs () {
	((h=${1}/3600))
	((m=(${1}%3600)/60))
	((s=${1}%60))
	printf "%02d:%02d:%02d\n" $h $m $s
}

# This script will monitor your job's utilized number of CPUs, actual RAM used, and actual VM (Virtual Memory aka swap).

function resourceSampler(){
    MONITOR_PID="$1"
    while true; do
        if [[ -f "$RSAMPLER_TERM" ]]; then
            rm -f "$RSAMPLER_TERM"
            break
        fi
        DATE=$(date +'%Y-%m-%d %H:%M:%S')
        ps --forest -o pcpu=,rss=,vsz= -g $MONITOR_PID \
              | awk -v date="$DATE" '{pcpu+=$1; rss+=$2; vsz+=$3} END {print date "," pcpu "," rss "," vsz}' >> "${RSAMPLER_OUTPUT}"
        sleep ${RSAMPLER_FREQ_SEC}
    done

    # cpu usage divided by 100 to convert percent usage into number of cpus in use
    CPU_USAGE=$(cut -f2 -d',' "$RSAMPLER_OUTPUT" \
                | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1}
                       END {printf "%15.2f%15.2f%15.2f", min/100, total/count/100, max/100}')
    RAM_USAGE=$(cut -f3 -d',' "$RSAMPLER_OUTPUT" \
                | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1}
                       END {printf "%15.2f%15.2f%15.2f", min/1024, total/count/1024, max/1024}')
    VM_USAGE=$(cut -f4 -d',' "$RSAMPLER_OUTPUT" \
                | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1}
                       END {printf "%15.2f%15.2f%15.2f", min/1024, total/count/1024, max/1024}')
    echo "CPUs Utilized min/avg/max: $CPU_USAGE"
    echo "RAM Used (MB) min/avg/max: $RAM_USAGE"
    echo "VM Used (MB)  min/avg/max: $VM_USAGE"
}


# Join an array?
function join_by { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

# Create interactive job
cjob () {
	NEWJ=`screen -ls | grep Attached | grep -o '\.[a-z]*' | grep -o '[a-z]*'`
        NEWJ=${NEWJ:-ajob}
        if [ "${1}" == "b" ]; then intjob ${NEWJ}; fi
        if [ "${1}" == "m" ]; then mjob ${NEWJ}; fi
        if [ "${1}" == "s" ]; then sjob ${NEWJ}; fi
        if [ "${1}" == "ms" ]; then msjob ${NEWJ}; fi
	if [ "${1}" == "day" ]; then dayjob ${NEWJ}; fi
}

# Get job id from a pattern [in the job's name]
jobid () {
	qstat -f | grep -e "Job Id" -e exec_host | grep -C 1 ${1} | tail -n 2
}

function backup () {
	MYPATH=`dirname ${1}`
	NEWNAME=${MYPATH}/.`basename ${1}`_`date +'%Y_%m_%d_%H_%M_%S'`
	echo "Backing to ${NEWNAME}"
	if [ -d ${1} ]; then
		echo "It is a folder: changing name to ${MYPATH}/.`basename ${1}`_`date +'%Y_%m_%d'`"
		cp -R ${1} "${MYPATH}/.`basename ${1}`_`date +'%Y_%m_%d'`"
	else
		cp ${1} ${NEWNAME}
	fi
}

function fastq_reads_n () {
	if echo "${1}" | grep -q "gz$" ; then
		echo "Using zcat: ${1}"
		NREADS=$(($(zcat $1 | wc -l) / 4 | bc))
	else
		echo "Using cat: ${1}"
		NREADS=$(($(cat $1 | wc -l) / 4 | bc))
	fi
	echo ${NREADS}
}

killsc () {
	screen -ls | grep .
        echo 'Enter screen to be killed:'; read session
        screen -S ${session} -X quit
}

scrf () {
	screen -ls | grep \(
        echo 'Enter screen to be attached:'; read session
        screen -r ${session}
}

qsum () {
  echo -e "User\t\tQ\tR\tC\tH\tQueue\t\tNodes"
  for u in `qstat -a | awk '/herman/ { print $2 }' | sort -u`; do
    nr="`qstat -u $u | grep herman | grep " R " | wc -l`"
    nq="`qstat -u $u | grep herman | grep " Q " | wc -l`"
    nc="`qstat -u $u | grep herman | grep " C " | wc -l`"
    nh="`qstat -u $u | grep herman | grep " H " | wc -l`"
    queue=(`qstat -u $u | awk '/herman/ { print $3 }' | sort -u`)
    nody=(`qstat -n1u $u | awk '/herman/ { print $12 }' | sort -u | grep compute | sed 's/compute.//'`)
    if [ ${#u} -gt 7 ]; then
      echo -e "${u}\t${nq}\t${nr}\t${nc}\t${nh}\t${queue[@]}\t\t${nody[@]}"
    else
      echo -e "${u}\t\t${nq}\t${nr}\t${nc}\t${nh}\t${queue[@]}\t\t${nody[@]}"
    fi
  done
}

findmapf () {
	echo 'Pattern:'; read IDT
        PATHIES=(`ls /mnt/NGSAnalyses/RNA-Seq/Mapping/ | grep $IDT`)
        for PATHY in ${PATHIES[@]}; do
                RFILE=/mnt/NGSAnalyses/RNA-Seq/Mapping/${PATHY}/report.html; if [ ! -s ${RFILE} ]; then continue; fi
                echo "/mnt/NGSAnalyses/RNA-Seq/Mapping/${PATHY}/"
                echo "`grep "Genome:" /mnt/NGSAnalyses/RNA-Seq/Mapping/${PATHY}/report.html`"
        done
}

function join_by {
  local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"
}

function find_job(){
  # qstat -fu ${USER} | grep -E 'Job_Name|Job Id|job_state'
  MYIDS=$(
    qstat -f | grep -E 'Job_Name|Job Id|job_state' |
      grep -EB 2 " = R| = Q| = H" |
      sed 's/Id: /Id_/g; s/ = /: /g; s/.herman.*/:/g' |
      grep -C 3 ${1} | grep "Id_" | sed -E 's/.*Id_|://g'
  )
  join_by ":" ${MYIDS}
}
# qsub -W depend=afterok:$(find_job "jobs_pattern") pbs_file

summrep () {
	qstat -q
        showq | grep "active" | more
}

function archive () {
  R -e "file.archive('${1}')" --slave
}
