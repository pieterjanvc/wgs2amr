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

if [ -z ${diamondPackage+x} ]; then echo -e "\nDiamond script location missing. Set -d /pathToDiamond/diamond \n"; exit 1; fi;
if [ -z ${rPackage+x} ]; then echo -e "\nR location missing. Set -r to /pathToR/bin/Rscript or the name of the R module e.g. R/3.5.0"; exit 1; fi;
if [ -z ${readFile1+x} ]; then echo -e "\nNo read file specified. Use -f /path/readFile1.fastq.gz"; exit 1; fi;
if [ -z ${readFile2+x} ]; then echo -e "\nOne read file provided, assuming single read file used \n"; fi;
if [ -z ${outputFolder+x} ]; then outputFolder=$projectFolder\RESULTS/; fi;

#--- START SCRIPT ---
sessionID=`date +%s`

touch $projectFolder\log
echo -e "\n--" `date` "-- Start AMR prediction pipeline \n" \
	"Session ID:" $sessionID "\n Files:\n\t" $readFile1 \
	`if [ ! -z ${readFile2+x} ]; then echo "\n\t" $readFile2; fi;`\
	"\n" >> $projectFolder\log

echo -e "\n--" `date` "-- Start AMR prediction pipeline \n" \
	"Session ID:" $sessionID "\n Files:\n\t" $readFile1 \
	`if [ ! -z ${readFile2+x} ]; then echo "\n\t" $readFile2; fi;`
	
mkdir -p $outputFolder

#Count the number of reads in the file
echo `date +"%T"` "Counting reads..."
readCounts=`zcat $readFile1 | wc -l`
echo `date +"%T"` " done"

#Run DIAMOND aligner on the files
echo `date +"%T"` "Start DIAMOND alignment File 1 ..."
$diamondPackage blastx \
-d $projectFolder\scriptsAndData/ncbi_ab_resistance_genes.dmnd \
-q $readFile1 \
-o $projectFolder\temp/$sessionID\_1.diamondOutput \
-f 6 qseqid qframe qcovhsp qlen qstart qend sseqid slen sstart send evalue length pident nident gapopen
echo `date +"%T"` "DIAMOND alignment File 1 finished"

if [ ! -z ${readFile2+x} ]
then
	echo `date +"%T"` "Start DIAMOND alignment File 2 ..."
	$diamondPackage blastx \
	-d $projectFolder\scriptsAndData/ncbi_ab_resistance_genes.dmnd \
	-q $readFile2 \
	-o $projectFolder\temp/$sessionID\_2.diamondOutput \
	-f 6 qseqid qframe qcovhsp qlen qstart qend sseqid slen sstart send evalue length pident nident gapopen
	echo `date +"%T"` "DIAMOND alignment File 2 finished"
fi

#Run the R script
echo `date +"%T"` "Start processing in R ..."

if grep -q Rscript <<< $rPackage
then 
	#Pointed to RScript location
	$rPackage $projectFolder\scriptsAndData/amrPrediction.R $projectFolder $outputFolder $readFile1 $readCounts $sessionID
else
	#Loading R module
	module load $rPackage
	Rscript $projectFolder\scriptsAndData/amrPrediction.R $projectFolder $outputFolder $readFile1 $readCounts $sessionID
fi

echo `date +"%T"`  " done"

#Finishing script
echo -e `date +"%T"`  "Pipeline with sessionID" $sessionID "finished\n\tresults saved in" $outputFolder >> $projectFolder\log
echo -e "\n\n---------------\n" `date`  \
	"Pipeline successfully finished\n\t results saved in" $outputFolder "\n---------------\n"
