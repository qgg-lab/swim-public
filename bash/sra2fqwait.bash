#! /bin/bash
# ============================================================

# bash to call sbatch sra2fq and work with parallel
# ============================================================

# get command line options
# ============================================================

args=`getopt -o "r:,d:,s:,i" -l "resource:,dir:,sample:,srr:" -- "$@"`
>&2 echo "command arguments given: $args"
eval set -- "$args"

# parse arguments
# ============================================================

while true;
do
  case $1 in

    -r|--resource)
      resource=$2
      shift 2;;

    -d|--dir)
			dir=$2
			shift 2;; 

    -s|--sample)
		  sample=$2
			shift 2;;
    
    -i|--srr)
      srr=$2
      shift 2;;

    --)
      shift
      break;;
      
  esac
done

# check for file availability
# ============================================================

if [[ -e $resource/resource.env && -e $resource/sbatch/sra2fq.sbatch ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: all slurm scripts found."
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: cannot find some slurm scripts."
  exit 1
fi

# run sbatch
# ============================================================

iter=`echo $srr | sed 's/,/+/g'`

mkdir "$dir"/"$sample"

jobID=$(sbatch --export=env=/mnt/research/qgg/resource/swim/resource.env,dir=$dir/$sample,srr=$iter --output="$dir"/"$sample"/sra.prepSeq.out --error="$dir"/"$sample"/sra.prepSeq.err $resource/sbatch/sra2fq.sbatch | cut -d " " -f 4)

sleep 1m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
  sleep 5m
  jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

reRunCounter=0
until [ `wc -l "$dir"/"$sample"/sra.prepSeq.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$dir"/"$sample"/sra.prepSeq.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting sra job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $dir/$sample/sra.fail.out
	sraTime=$(($reRunCounter * 4 + 28))
	sraMem=$(($reRunCounter * 5 + 10))G
  if [[ -e "$dir"/"$sample"/file.info ]]
  then
    rm "$dir"/"$sample"/file.info
    rm -f "$dir"/"$sample"/*.fastq.gz
  fi
  jobID=$(sbatch --export=env=/mnt/research/qgg/resource/swim/resource.env,dir=$dir/$sample,srr=$iter --time=$sraTime:00:00 --mem=$sraMem --output="$dir"/"$sample"/sra.prepSeq.out --error="$dir"/"$sample"/sra.prepSeq.err $resource/sbatch/sra2fq.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l "$dir"/"$sample"/sra.prepSeq.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$dir"/"$sample"/sra.prepSeq.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed sra download."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> "$dir"/"$sample"/sra.prepSeq.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):info: error for sra download."
	exit 1
fi
