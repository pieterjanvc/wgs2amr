#!/bin/bash
set -e #Stop script if any of the commands fails

while getopts ":d:r:f:s:o:" opt; do
  case $opt in
    d) diamondPackage="$OPTARG"
    ;;
	r) rPackage="$OPTARG"
    ;;
	f) readFile1="$OPTARG"
    ;;
	s) readFile2="$OPTARG"
    ;;
	o) outputFolder="$OPTARG"
    ;;	  
    \?) echo "Unknown argument provided"
	;;
  esac  
done

projectFolder=`dirname $0`/

if [ -z ${diamondPackage+x} ]; then echo "Diamond script location missing. Set -d /pathToDiamond/diamond"; exit 1; fi;
if [ -z ${rPackage+x} ]; then echo "R location missing. Set -r to /pathToR/bin/Rscript or the name of the R module e.g. R/3.5.0"; exit 1; fi;
if [ -z ${readFile1+x} ]; then echo "No read file specified. Use -f /path/readFile1.fastq.gz"; exit 1; fi;
if [ -z ${readFile2+x} ]; then echo "One read file provided, assuming single read file used"; fi;
if [ -z ${outputFolder+x} ]; then outputFolder=$projectFolder\RESULTS/; fi;

#--- START SCRIPT ---
touch $projectFolder\log
echo -e "\n--" `date` "-- Start AMR prediction pipeline" >> $projectFolder\log
mkdir -p $outputFolder

#Make sure the temp folder is clean
rm -f $projectFolder\temp/*

#Count the number of reads in the file
echo `date +"%T"` "Counting reads..." >> $projectFolder\log
readCounts=`zcat $readFile1 | wc -l`
echo `date +"%T"` " done" >> $projectFolder\log

#Run DIAMOND aligner on the files
echo `date +"%T"` "Start DIAMOND alignment File 1 ..." >> $projectFolder\log
$diamondPackage blastx \
-d $projectFolder\scriptsAndData/ncbi_ab_resistance_genes.dmnd \
-q $readFile1 \
-o $projectFolder\temp/diamondFile1.diamondOutput \
-f 6 qseqid qframe qcovhsp qlen qstart qend sseqid slen sstart send evalue length pident nident gapopen
echo `date +"%T"` " done" >> $projectFolder\log

if [ ! -z ${readFile2+x} ]
then
	echo `date +"%T"` "Start DIAMOND alignment File 2 ..." >> $projectFolder\log
	$diamondPackage blastx \
	-d $projectFolder\scriptsAndData/ncbi_ab_resistance_genes.dmnd \
	-q $readFile2 \
	-o $projectFolder\temp/diamondFile2.diamondOutput \
	-f 6 qseqid qframe qcovhsp qlen qstart qend sseqid slen sstart send evalue length pident nident gapopen
	echo `date +"%T"` " done" >> $projectFolder\log
fi

#Run the R script
echo `date +"%T"` "Start processing in R ..." >> $projectFolder\log

if grep -q Rscript <<< $rPackage
then 
	#Pointed to RScript location
	$rPackage $projectFolder\scriptsAndData/amrPrediction.R $projectFolder $outputFolder $readFile1 $readCounts
else
	#Loading R module
	module load $rPackage
	Rscript $projectFolder\scriptsAndData/amrPrediction.R $projectFolder $outputFolder $readFile1 $readCounts
fi

echo `date +"%T"`  " done" >> $projectFolder\log

#Finishing script
echo `date +"%T"`  "Pipeline finished -- results are in " $outputFolder >> $projectFolder\log
