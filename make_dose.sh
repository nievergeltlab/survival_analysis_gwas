#Specify where probability format genotypes are stored
probdir=/home/caroline/Desktop/pgbd/qc/imputation/dasuqc1_bip_pgbd_mix_am-qc.hg19.ch.fl/qc1

#Specify where output will be stored 
outdir=/home/caroline/Desktop/pgbd/qc/imputation/dasuqc1_bip_pgbd_mix_am-qc.hg19.ch.fl/dosecnt

#Specify working directory (error logs will go here)
workingdir=/home/caroline/Desktop/pgbd/qc/imputation

#Get the number of subjects (assuming that it is the same across all files)
nsub=$(wc -l "$probdir"/dos_bip_pgbd_mix_am-qc.hg19.ch.fl.chr5_117_120.out.dosage.fam  | awk '{print $1}')


ls $probdir | grep .gz$ > files_To_make_dose.txt

#Number of commands to run is a function of the number of files
ncommands=$(wc -l files_To_make_dose.txt | awk '{print $1}' )

#Make a job code, where 'nodesize' processes will run on each node simultaneously
 nodesize=8
 nodeuse=$(($nodesize - 1))

#Total number of jobs = Number of commands / number of commands used per job, rounded up 
totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

qsub make_dose.pbs  -t 1-$totjobs -d $workingdir -F "-s $nsub -d $outdir -p $probdir -n $nodeuse -o files_To_make_dose.txt"

