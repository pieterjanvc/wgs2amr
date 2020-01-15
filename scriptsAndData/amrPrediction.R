#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE) #arguments specified in the bash file

suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(xgboost))

options(stringsAsFactors = F)

#######################################################################
####### Create gene coverage table for metagenomic BLAST output #######
#######################################################################

##### INPUT needed ###########
##############################

# ASSIGN GLOBAL VARIABLES
#-----------------------
projectFolder = args[[1]]
resultsFolder = args[[2]]
fileName = str_extract(args[[3]], "\\w+(?=\\.fastq)")
readCounts = as.integer(args[[4]])
sessionID = args[[5]]

setwd(projectFolder)
refGenes = read.csv("scriptsAndData/refGenes.csv")
seqIdentityPercent = 90 #filter 1: minimum sequence identity of aligned parts
lengthMatchPercent = 90 #filter 2: minimum length percentage of read used in alignment


##### FUNCTIONS ###########
###########################

#FILTER FUNCTION 1 - Seqence identity
#------------------------------------

seqIdentityFilter = function(myData, filterPercentage){
  myData = myData[myData$pident >= filterPercentage,] #Percentage sequence identity is already calculated by BLAST (pident)
  
  return(myData)
}

#FILTER FUNCTION 2 - Length match
#--------------------------------

lengthMatchFilter = function(myData, filterPercentage, allowGaps = F){
  
  if(!allowGaps){
    #We can't use gaps given this messes up the length calculations
    #Gaps are not counted as length and we don't know the size
    myData = myData[myData$gapopen == 0, ] #Filter out reads with gaps
  }
  
  
  #Convert read length to AA (is in bases)
  myData$qlen = floor(myData$qlen / 3) #Some are not exact multiple of 3
  myData$lengthMatch = (myData$length / myData$qlen) * 100 #Calculate percentage of read used in alignment
  myData = myData[myData$lengthMatch >= filterPercentage, ] #Filter
  
  return(myData)
}


#Create a coverage table
#---------------------------------------------
buildCoverageTable = function(diamondFolder, sampleName, seqIdentityPercent, lengthMatchPercent, refGenes, readCounts, sessionID){
  
  myFiles = list.files(diamondFolder, pattern = paste0(sessionID, "_[1,2]\\.diamondOutput"), full.names = T)

  diamondTable = map_df(myFiles, function(x){
    if(file.info(x)$size == 0){
      return(data.frame())
    } else {
      read.table(x, sep = "\t", stringsAsFactors = F,
                 blank.lines.skip = TRUE, skipNul = TRUE, fill=T, quote="")
    }
  })

  colnames(diamondTable) = c("qseqid", "qframe", "qcovhsp", "qlen", "qstart",
                             "qend", "sseqid", "slen", "sstart", "send", "evalue",
                             "length", "pident", "nident", "gapopen")
  

  ##### PART 1 - Filtering #####
  ##############################
  myResult = seqIdentityFilter(diamondTable, seqIdentityPercent)
  myResult = lengthMatchFilter(myResult, lengthMatchPercent)

  if(nrow(myResult) > 0){
    ##### PART 2 - Calculate the coverage of the genes #####
    ########################################################
    
    #Create 3 new columns:
    #Gene coverage without taking matching percentage into account
    #Gene coverage taking highest similatrity percentage into account (pIdent)
    #Total number of alignments to the gene
    
    coverageTable = refGenes %>% add_column(mixedSampleName = sampleName, .before = T)
    
    #For every gene present, calculate its coverage
    genesCoverage = map_df(unique(myResult$sseqid), function(seqId){
      
      #Get all reads that align to a gene
      allReads = myResult %>% filter(sseqid == seqId)
      
      cover = rep(0, allReads$slen[1])
      cover = t(apply(allReads %>% select(sstart, send, pident), 1, function(x) {
        cover[x[1]:x[2]] = x[3]
        cover
      }))
      
      cDepth = apply(cover, 2, function(x) sum(x > 0))
      adaptiveC = apply(cover, 2, function(x) {
        if(sum(x) == 0){
          0
        } else {
          max(x[x > 0])
        }
      })
      simpleCoverage = mean(ifelse(colSums(cover) > 0, 1, 0))
      adaptiveCoverage = mean(adaptiveC)/100
      
      nAlignments = nrow(allReads)
      
      return(data.frame(seqId = seqId, simpleCoverage = simpleCoverage, 
                        adaptiveCoverage= adaptiveCoverage, nAlignments = nAlignments))
      })
    
    coverageTable = coverageTable %>% left_join(genesCoverage, by = c("geneIdNCBI" = "seqId")) %>% 
      filter(adaptiveCoverage > 0)
    
    if(nrow(coverageTable) == 0){
      
      return(data.frame())
      
    } else {
      #Adjust number of aligments to seq depth
      coverageTable$adjustedReadCount = coverageTable$nAlignments / (coverageTable$geneLength * ((readCounts * length(myFiles)) / 10000000))
      
      return(coverageTable)
    }
    
    
  } else {
    return(data.frame())
  }
  
}


##### Generate data for ML models ###########
#############################################

#Create coverage table from diamondoutput using the ref genes
cat(paste0(format(Sys.time(), "%H:%M:%S"), "  creating coverage table\n"))

coverageTable = buildCoverageTable("temp/", fileName, seqIdentityPercent, 
                              lengthMatchPercent, refGenes, readCounts, sessionID) %>% 
  filter(adaptiveCoverage > 0.9)

write.csv(coverageTable, paste0("temp/", sessionID, ".csv"), row.names = F)


#Generate inputvector for ML
if(nrow(coverageTable) > 0){
  inputVector = refGenes %>% distinct(simplifiedId) %>% 
    left_join(coverageTable %>% 
                group_by(simplifiedId) %>% 
                summarise(arc = mean(adjustedReadCount)), by = c("simplifiedId" = "simplifiedId")) %>% 
    arrange(simplifiedId)
  
  inputVector = ifelse(is.na(inputVector$arc), 0, inputVector$arc)  
} else {
  inputVector = rep(0, 1027)
}

##### Run models ###########
############################
cat(paste0(format(Sys.time(), "%H:%M:%S"), "  running predictions\n"))

models = readRDS("scriptsAndData/medianModels.rds")
toKeep = readRDS("scriptsAndData/colsToKeep.rds")
inputVector = inputVector[toKeep]
ABnames = c( "cefepime", "cefotaxime", "ceftriaxone", "ciprofloxacin", 
             "gentamicin","levofloxacin", "meropenem", "tobramycin")
names(models) = ABnames
  
newPredictions = map_df(1:length(ABnames), function(i){
  
  #Get all inputs for that AB model
  myInputs = t(matrix(inputVector[models[[ABnames[i]]]$colsToKeep]))
  colnames(myInputs) = models[[ABnames[i]]]$model$feature_names
  myMisClass = models[[i]]$myConfidence
  
  #Convert output to range between 0.5 - 1 (close to 0 or to 1 is supposed to be better prediction)
  myMisClass$adjustedOutput = ifelse(myMisClass$confidence < 0.5, 1 - myMisClass$confidence, myMisClass$confidence)
  myMisClass = myMisClass[order(myMisClass$adjustedOutput),]
  
  #Get misclass rate
  nFalse = sum(myMisClass$prediction == F)
  myMisClass$misclassRate  = sapply(myMisClass$adjustedOutput, function(x){
    sum(myMisClass$prediction == F & myMisClass$adjustedOutput >= x) / nFalse
  })
  
  #Get 3dr polynomial approximation
  # myPoly = unlist(lm(misclassRate ~ poly(adjustedOutput, poly = 3), data = myMisClass)$coefficients)
  myPoly = lm(misclassRate ~ poly(adjustedOutput, poly = 3), data = myMisClass)
  
  #Predict new input, but clip to min / max AO of test set 
  myPrediction = predict(models[[ABnames[i]]]$model, myInputs) 
  
  newInput = ifelse(myPrediction < 0.5, 1 - myPrediction, myPrediction)
  myMax = max(myMisClass$confidence)
  myMin = min(myMisClass$confidence)
  newInput = sapply(newInput, function(x){
    if(x > myMax){
      myMax
    } else if(x < myMin){
      myMin
    } else {
      x
    }
  })

  myReliability = 1- predict(myPoly, 
                             data.frame(adjustedOutput = newInput), 
                             type = "response")
  
  #Clip results > 1 or < 0
  myReliability = sapply(myReliability, function(x){
    if(x > 1){1} else if(x < 0){0} else{RoundTo(x, 0.02)}
  })
  
  #Save to data frame
  data.frame(antibiotic = ABnames[i], prediction = ifelse(myPrediction < 0.5, "suseptible", "resistant"),
                         reliability = myReliability)
})

#Write the output to file
write.csv(newPredictions, paste0(resultsFolder, fileName, "_AMRpredictions.csv"), row.names = F)
