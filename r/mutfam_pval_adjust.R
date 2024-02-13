# Adjust p-vals
# 23/03/2016
# Apply Benjami-Hochberg to overall MutFam_p tables AFTER concatenating CS cluster batch jobs
# into one big csv

adjust_pval_bh <- function(comboCSV){
	library(data.table);
	library(plyr);
  mutFam_p_ALL <- data.table(read.csv(file=comboCSV, 
  										sep=",", 
  										stringsAsFactors=FALSE, 
  										as.is=TRUE, 
  										header=TRUE
  									 )
  							);

	pvc <- arrange(mutFam_p_ALL[,pval_corr:=p.adjust(mutFam_p_ALL[,pval], method="BH")], pval_corr);

	file_out <- paste(sub(".csv", "", comboCSV),"_BH_corr.csv", sep="");
	write.csv(pvc, file=file_out);

}




