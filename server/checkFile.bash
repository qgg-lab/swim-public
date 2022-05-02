#! /bin/bash
# ============================================================

# get command line options
# ============================================================

args=`getopt -o "r:,p:,m:" -l "resource:,ped:,map:" -- "$@"`
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

    -p|--ped)
			ped=$2
			shift 2;; 

  	-m|--map)
			map=$2
			shift 2;;

    --)
      shift
      break;;

  esac
done

# check for file availability
# ============================================================

if [[ -e $resource ]]
then
  echo "$(date +"%m-%d-%Y-%T"):info: found resource dir."
else
  echo "$(date +"%m-%d-%Y-%T"):error: cannot find resource dir."
  exit 1
fi

source $resource/swim.resource.env

# check if ped map is expected 
# ============================================================

pednCol=`head -n 1 $ped | awk '{print NF}'`
mapnCol=`head -n 1 $map | awk '{print NF}'`

if [ $mapnCol -ne 4 ]
then
  echo "$(date +"%m-%d-%Y-%T"):error: map file has number of columns not equal to 4 (did you upload a ped file by mistake?)."
  exit 1
fi

if [ $pednCol -lt 10 ]
then
  echo "$(date +"%m-%d-%Y-%T"):error: ped file has number of columns smaller than 10 (did you upload a map file by mistake?)."
  exit 1
fi

mapnRow=`wc -l $map | awk '{print $1}'`

if [ $mapnRow -lt 10 ]
then
  echo "$(date +"%m-%d-%Y-%T"):error: map file has number of rows smaller than 10."
  exit 1
fi

# check if there are two many individuals
# ============================================================

pednRow=`wc -l $ped | awk '{print $1}'`

if [ $pednRow -gt 20000 ]
then
  echo "$(date +"%m-%d-%Y-%T"):error: two many (> 20,000) individuals. The SWIM server has a limit of imputing less than 20,000 individuals at a time."
  exit 1
fi

# check if ped contains non atcg character
# ============================================================

nonATCG=`head -n 5 $ped | sed 's/\t/ /g' | cut -d " " -f 7- | awk '$1 !~ /[ATCG0atcg-]/' | wc -l`

if [ $nonATCG -gt 0 ]
then
  echo "$(date +"%m-%d-%Y-%T"):error: ped file contains non A/T/C/G/0/- character as alleles."
  exit 1
fi

# check if map file contains only snps not supported by the platform
# ============================================================

chipSNP=`awk '{print $2}' $map | sort | join -t $'\t' - $resource/chip.info | wc -l`

if [ $chipSNP -eq 0 ]
then
  echo "$(date +"%m-%d-%Y-%T"):error: none of the SNPs matches our supported SNP IDs."
  exit 1
fi

# output success
# ============================================================

echo "$(date +"%m-%d-%Y-%T"):progress: ped/map pass tests."
