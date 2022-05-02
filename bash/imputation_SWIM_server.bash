###### imputation for SWIM server online ######
######  check for file availability    ######

if [[ -e $resource/PorcineSNPillumina_final.info && -e $resource/PorcineSNPillumina_final_SNPchr && \
	  -e $resource/PorcineSNPillumina_final_SNPpos && -e $resource/PorcineSNPillumina_final_SNP && \
	  -e $resource/sbatch/swim_QC.sbatch && -e $resource/sbatch/swim_snpchip_allele.sbatch && \
	  -e $resource/sbatch/swim_UpdateAllele0_2VCF.sbatch && -e $resource/sbatch/swim_flip_swap.sbatch  &&\
	  -e $resource/sbatch/swim_phasing_target.sbatch && -e $resource/sbatch/imputation_impute5_final.sbatch ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: all slurm scripts found."
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: cannot find some slurm scripts."
  exit 1
fi

######  Only supports input vcf.gz file format  ######
######  check for input file availability    ######
inputfile=test.vcf.gz
inputname=$(echo $inputfile |sed 's%.vcf.gz%%g')
tmpfile=/mnt/home/dingrong/test

case "$inputfile" in
*.vcf.gz) echo "input file format is "vcf.gz"."
;;
*) 
echo "Sorry, please import the "vcf.gz" format."
exit 1 ;;
esac

if [ `gunzip -c $inputfile |grep -v ^# |awk '{print $1 }' |sort |uniq |wc -l` -eq 1]
then 
chr=`gunzip -c $inputfile |grep -v ^# |awk '{print $1 }' |sort |uniq`
else
echo "Sorry, the input file can only contain one autosome."
exit 1
fi

######  QC  ######
######  Software PLINK   ######

cd $tmpfile

if [ ! -d $tmpfile/log ]
then mkdir -p $tmpfile/log
fi

# run sbatch
jobID=$(sbatch --export=tmpfile=$tmpfile,inputname=$inputname --output=log/swim_QC.out --error=log/swim_QC.err $resource/sbatch/swim_QC.sbatch | cut -d " " -f 4)
sleep 2m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 2m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success

reRunCounter=0

until [ `wc -l $tmpfile/log/swim_QC1.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $tmpfile/log/swim_QC1.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 1 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting swim_QC job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> log/swim_QC.out
	QCTime=$(($reRunCounter))
	jobID=$(sbatch --export=tmpfile=$tmpfile,inputname=$inputname --time=$QCTime:30:00 --output=log/swim_QC.out --error=log/swim_QC.err $resource/sbatch/swim_QC.sbatch | cut -d " " -f 4)
	sleep 2m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 2m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $tmpfile/log/swim_QC1.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $tmpfile/log/swim_QC1.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed swim_QC."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> log/swim_QC.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for swim_QC."
	exit 1
fi

######  2.remove duplicate-vars, update Alleles and convert to vcf  ######

#Select1.filter duplicate-vars 

if [ ` wc -l  "$inputname".dupvar | awk '{print $1}'` -ne 0 ]
then

awk '{if ($1>0 && $1<19) print $0}' "$inputname".dupvar \
| awk -F " " '{for (i=4;i<=NF;i++)printf("%s ", $i);print ""}' \
| sed 's% %\n%g' |sed '/^$/d' >"$inputname"_duplicate.txt

else

cat /dev/null > "$inputname"_duplicate.txt

fi

#Select2.Allele is 0 converted to ref

if [ `awk '{if($5==0 || $6==0) print $0}' "$inputname"_snpPosChr.bim |wc -l |  awk '{print $1}' `  -ne 0 ]
then
# run sbatch

jobID=$(sbatch --export=infofile=$infofile,bimfile="$inputname"_snpPosChr.bim,outname="$tmpfile"/"$inputname" --output=log/"$inputname"_allele0_renew.out --error=log/"$inputname"_allele0_renew.err $resource/sbatch/swim_snpchip_allele.sbatch | cut -d " " -f 4)
sleep 2m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 2m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success

reRunCounter=0

until [ `wc -l $tmpfile/log/"$inputname"_allele0_renew.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $tmpfile/log/"$inputname"_allele0_renew.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 1 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting swim_snpchip_allele job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $tmpfile/log/"$inputname"_allele0_renew.out
	allele0_renewMem=$(($reRunCounter * 16 ))G
	allele0_renewTime=$($reRunCounter))
	jobID=$(sbatch --export=infofile=$infofile,bimfile="$inputname"_snpPosChr.bim,outname="$tmpfile"/"$inputname" --mem=$allele0_renewMem --time=$allele0_renewTime:20:00 --output=log/"$inputname"_allele0_renew.out --error=log/"$inputname"_allele0_renew.err $resource/sbatch/swim_snpchip_allele.sbatch | cut -d " " -f 4)
	sleep 2m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 2m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $tmpfile/log/"$inputname"_allele0_renew.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $tmpfile/log/"$inputname"_allele0_renew.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed allele0_renew."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $tmpfile/log/"$inputname"_allele0_renew.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for allele0_renew."
	exit 1
fi

else

cat /dev/null > "$inputname"_UpdateAlleles.txt

fi

#convert to vcf
# run sbatch

jobID=$(sbatch --export=inputname=$inputname --output=log/"$inputname"_UpdateAllele0_2VCF.out --error=log/"$inputname"_UpdateAllele0_2VCF.err $resource/sbatch/swim_UpdateAllele0_2VCF.sbatch | cut -d " " -f 4)
sleep 2m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 2m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success

reRunCounter=0

until [ `wc -l $tmpfile/log/"$inputname"_UpdateAllele0_2VCF.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $tmpfile/log/"$inputname"_UpdateAllele0_2VCF.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 1 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting UpdateAllele0_2VCF job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $tmpfile/log/"$inputname"_UpdateAllele0_2VCF.out
	UpdateAllele0_2VCFMem=$(($reRunCounter * 32 ))G
	UpdateAllele0_2VCFTime=$($reRunCounter))
	jobID=$(sbatch --export=inputname=$inputname --output=log/"$inputname"_UpdateAllele0_2VCF.out --time=$UpdateAllele0_2VCFTime:30:00 --mem=$UpdateAllele0_2VCFMem --error=log/"$inputname"_UpdateAllele0_2VCF.err $resource/sbatch/swim_UpdateAllele0_2VCF.sbatch | cut -d " " -f 4)
	sleep 10m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 10m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $tmpfile/log/"$inputname"_UpdateAllele0_2VCF.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $tmpfile/log/"$inputname"_UpdateAllele0_2VCF.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed UpdateAllele0_2VCF."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> log/"$inputname"_flip_swap.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for UpdateAllele0_2VCF."
	exit 1
fi

######  3.flip and swap  ######
# run sbatch

jobID=$(sbatch --export=inputname=$inputname,tmpfile=$tmpfile --output=log/"$inputname"_flip_swap.out --error=log/"$inputname"_flip_swap.err $resource/sbatch/swim_flip_swap.sbatch | cut -d " " -f 4)
sleep 2m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 2m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success

reRunCounter=0

until [ `wc -l log/"$inputname"_flip_swap.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" log/"$inputname"_flip_swap.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 1 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting flip and swap job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> log/"$inputname"_flip_swap.out
	flip_swapMem=$(($reRunCounter * 32 ))G
	flip_swapTime=$($reRunCounter))
	jobID=$(sbatch --export=inputname=$inputname,tmpfile=$tmpfile --output=log/"$inputname"_flip_swap.out --time=$flip_swapTime:30:00 --mem=$flip_swapMem --error=log/"$inputname"_flip_swap.err $resource/sbatch/swim_flip_swap.sbatch | cut -d " " -f 4)
	sleep 10m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 10m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $tmpfile/log/"$inputname"_flip_swap.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $tmpfile/log/"$inputname"_flip_swap.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed flip and swap."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> log/"$inputname"_flip_swap.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for flip and swap."
	exit 1
fi

######  4.phasing  ######

# run sbatch

jobID=$(sbatch --export=tmpfile=$tmpfile,inputname=$inputname,chr="$chr" --output=log/"$inputname".chr"$chr".shapeit4.phasing.out --error=log/"$inputname".chr"$chr".shapeit4.phasing.err $resource/sbatch/swim_phasing_target.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success

reRunCounter=0

until [ `wc -l log/"$inputname".chr"$chr".shapeit4.phasing.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" log/"$inputname".chr"$chr".shapeit4.phasing.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 1 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting phase job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> log/"$inputname".chr"$chr".shapeit4.phasing.out
	phaseMem=$(($reRunCounter * 32 ))G
	phaseTime=$($reRunCounter))
	jobID=$(sbatch --export=tmpfile=$tmpfile,inputname=$inputname,chr="$chr" --output=log/"$inputname".chr"$chr".shapeit4.phasing.out --error=log/"$inputname".chr"$chr".shapeit4.phasing.err $resource/sbatch/swim_phasing_target.sbatch | cut -d " " -f 4)
	sleep 10m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 10m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l log/"$inputname".chr"$chr".shapeit4.phasing.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" log/"$inputname".chr"$chr".shapeit4.phasing.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed phase."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> log/"$inputname".chr"$chr".shapeit4.phasing.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for phase."
	exit 1
fi


######  5.imputation  ######
# run sbatch

jobID=$(sbatch --export=tmpfile=$tmpfile,inputname=$inputname,chr="$chr" --output=log/"$inputname".chr"$chr".imputation.out --error=log/"$inputname".chr"$chr".imputation.err $resource/sbatch/swim_imputation.sbatch | cut -d " " -f 4)
sleep 10m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success

reRunCounter=0

until [ `wc -l log/"$inputname".chr"$chr".shapeit4.phasing.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" log/"$inputname".chr"$chr".shapeit4.phasing.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 1 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting imputation job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> "$inputname".chr"$chr".imputation.err
	imputationMem=$(($reRunCounter * 64 ))G
	imputationTime=$($reRunCounter + 1 ))
	jobID=$(sbatch --export=tmpfile=$tmpfile,inputname=$inputname,chr="$chr" --time=$imputationTime:50:00 --mem=$imputationMem  --output=log/"$inputname".chr"$chr".imputation.out --error=log/"$inputname".chr"$chr".imputation.err $resource/sbatch/swim_imputation.sbatch | cut -d " " -f 4)
	sleep 20m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 10m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l log/"$inputname".chr"$chr".shapeit4.phasing.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" log/"$inputname".chr"$chr".shapeit4.phasing.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed imputation."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> log/"$inputname".chr"$chr".imputation.err
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for imputation."
	exit 1
fi


