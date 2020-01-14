# WGS2AMR
This pipeline was built to predict antimicrobial resistance (AMR) to 8 common antibiotics in 5 different Gram-negative species based on whole-genome sequencing data.

___

## Installation of dependencies
This project was tested with installation of all dependencies on a linux based machine.

### DIAMOND
This tool is used to align the bacterial sequencing files (fastq format) to a precompiled database of antimicrobial resistance genes. To install DIAMOND, follow the instruction on their [GitHub](https://github.com/bbuchfink/diamond) page. 

The database to align to is provided in this project and does not need to be compiled.

### R
Processing of the output as created by the DIAMOND aligner and predicing the AMR is done in R. For this project, R 3.5.0 was used, so this or any later version should work (though older versions are expected to work as well). After [installing R](https://www.r-project.org/) on your machine, the following packages need to be installed before running the pipeline:
```R
install.packages(c("tibble", "dplyr", "seqinr", "stringr", "purrr", "DescTools", "xgboost"))
```
___

## Input
The models were created to make predictions for the following bacterial species: Acinetobacter baumannii, Enterobacter cloacae, Escherichia Coli, Klebsiella aerogenes and Klebsiella pneumoniae. They were trained to be species independent and could theoretically make predictions on other Gram-negative species, but have not been validated for those thus far.

The files used as input should be next-generation-sequencing files in either fastq or compressed fastq.gz format. This can be in one single file or two separate files in case of pair-end read files. The `testFiles` folder in the `wgs2amr` folder contains dummy files that can be used to test the pipeline when all dependencies are installed.

___

## Running the pipeline
### Download the project folder
Download the `wgs2amr` folder from this GitHub page and put it on the same machine where you installed the dependencies.

### Call the wgs2amr.sh script
The `wgs2amr` folder contains the master script `wgs2amr.sh` that can be called with the following arguments

* -d : Path to the diamond script in the diamond folder
* -r : Either path to the Rscript in the R bin folder, or the name of the R module to load on a linux machine (e.g. R/3.5.0)
* -f : The first sequence file. If there is only one sequence file, this one should be set. Both fastq and fastq.gz are supported
* -s : The second sequence file in case of two pair-end reads files. Again, both fastq and fastq.gz are supported. In case of only one file, this argument can be omitted
* -o : The location of the output folder where the prediction results in csv format will be stored. If not set, this will default to the `RESULTS` folder within the `wgs2amr` folder

Here are some examples of a typical wgs2amr.sh call using the test data provided in the wgs2amr 
```Bash
#Using two pair-end reads input files and the link to the Rscript
/pathToScript/wgs2amr.sh \
-d '/pathToDiamond/diamond' \
-r '/pathToR/bin/Rscript' \
-f '/pathTo/wgs2amr/testFiles/readsFile1.fastq.gz' \
-s '/pathTo/wgs2amr/testFiles/readsFile2.fastq.gz'
```
```Bash
#Using one combined reads file and the name of the R module
/pathToScript/wgs2amr.sh \
-d '/pathToDiamond/diamond' \
-r 'R/3.5.0' \
-f '/pathTo/wgs2amr/testFiles/testFile.fastq.gz'
```
*Note:If an output folder is specified, make sure the path ends with a backslash*

## Output
The output file will be named after the first sequence file used in the format `<filename>_AMRpredictions.csv` and found in the default or specified output folder. 
The first column contains the names of the 8 antibiotics for which predictions were made, the second column the predicted susceptibility (susceptible or resistant), the last column the reliability, offering some idea of the degree the model thinks the output is correct (0 = very uncertain to 1 = near certain)

Example of the test file output:

|antibiotic    |prediction | reliability|
|:-------------|:----------|-----------:|
|cefepime      |resistant  |        0.70|
|cefotaxime    |resistant  |        0.78|
|ceftriaxone   |resistant  |        0.98|
|ciprofloxacin |resistant  |        0.96|
|gentamicin    |suseptible |        0.88|
|levofloxacin  |resistant  |        0.96|
|meropenem     |suseptible |        0.40|
|tobramycin    |resistant  |        0.66|
___

*NOTE: Temporary files like the diamond alignment files are stored in the `wgs2amr/temp` folder and can safely be removed after finishing the pipeline.*
