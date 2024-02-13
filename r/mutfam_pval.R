# MutFam p-val calcs
# 17/02/2016  Ash
# Done on cluster for accurate random trials

# 22/02/2016
# Allow passing of task/batchID.  This will be passed via $SGE_TASK_ID from
# the CS cluster run job in array mode.

library(data.table);
library(plyr)
source('mutfam_permut.R')

#BATCH_SIZE <- 50;
# 21/03/16 - non-completetion issues in some GMFP01 at size 50.
BATCH_SIZE <- 25;
PVAL_RANDS <- 10000;

# Check arguments including if numeric TASKID passed
args <- commandArgs(trailingOnly=TRUE);
if (length(args) == 0) {
  stop("Args required - 1: RData file, [2: task/batch no] ")
}

if ((length(args) == 1) | (length(args) == 2)) {
  datFile <- args[1];
  if (file.exists(datFile) == 0 ){
    stop(paste("File ", datFile, " not found.", sep=""))
  }
}

if ((length(args)==2) & !is.na(as.numeric(args[2]))){
  taskID <- as.integer(args[2]);
}else{
  taskID <- 0;
}

load(file=datFile);

# Deal with just a subset of the MutFams using whatever TASKID we're called with...
if (taskID==0){
  range <- 1:nrow(dt_combo);
}else{
  startBatch <- ((taskID - 1) * BATCH_SIZE) + 1; 

  # TaskID puts us out of range!
  if (startBatch > nrow(dt_combo)){
    stop(paste("Start index of taskID ", taskID, " is ", startBatch," > nrow(dt_combo) [",nrow(dt_combo),"]", sep=""));
  }else{
    endBatch <- (taskID * BATCH_SIZE);
    if (endBatch > nrow(dt_combo)){
      endBatch <- nrow(dt_combo);
    }
  } 
  range <- startBatch:endBatch;
}

#test
writeLines(datFile);
print(taskID);
print(startBatch);
print(endBatch);
print(range)
#stop("stop");

# Calc p-values
# "Presentation" (for convenience!)
# 19/02/2016 v2.1 - have changed to report fm_ratio not FF_MUTS twice!
dt_Combo_p_val <- arrange(dt_combo
                          [range,.(SUPERFAMILY_ID,
                              FUNFAM_NUMBER,
                              FF_NAME = NAME,
                              Superfamily_Name,
                              FF_MUTS,
                              ef,
                              GENE_MUTS_IN_FFS,
                              FF_SUM_GENE_LEN,
                              fm_ratio,
                              SUM_FF,
                              pval=1
                          )],
                          desc(ef));
rownames(dt_Combo_p_val) <- range;
setkey(dt_Combo_p_val,"SUPERFAMILY_ID","FUNFAM_NUMBER");

dt_Combo_p_val[, testCol:=paste("mfp.mut_permut(",PVAL_RANDS,",", GENE_MUTS_IN_FFS,",", FF_SUM_GENE_LEN,",", FF_MUTS,",", 1,",",SUM_FF,")",sep="")]
# Eval fn & update...
evs <- function(x){return (eval(x))};
dt_Combo_p_val[, pval:=evs(parse(text=testCol)), by=1:nrow(dt_Combo_p_val)];
dt_Combo_p_val <- arrange(dt_Combo_p_val, desc(ef));

# 25/02/16 - Added this (needs testing!)
# return Benjami-Hochberg adjusted p-vals too and sort by them...
# 23/03/16 - ** Can't do this here - don't have full table during batch jobs!!! - see mutfam_pval_adjust.R/adjust_pval **
####dt_Combo_pvc <- arrange(dt_Combo_p_val[,pval_corr:=p.adjust(dt_Combo_p_val[,pval], method="BH")], pval_corr);

out_file <- paste(sub(".RData","",datFile),"_PVAL_", PVAL_RANDS, "_BATCH_", BATCH_SIZE, "_TASKID_", taskID ,".csv", sep="");
# 25/02/16 write the BH version...
#write.csv(dt_Combo_p_val, file=out_file, row.names=range);
##write.csv(dt_Combo_pvc, file=out_file, row.names=range);

write.csv(dt_Combo_p_val, file=out_file, row.names=range);


# p<0.05 only
#dt_Combo_pval005 <- dt_Combo_p_val[pval<0.05,];
#data_file <- paste(datFile, ".csv",sep="");
#write.csv(dt_Combo_pval005, file=data_file);
