# Run p-val corrections
# Benjami-Hochberg correction of MutFam P-values

source('mutfam_pval_adjust.R');

#Example
# GMFP05
group_dir <- "data/GMFP05";
# MFP028 CHOL /I55
runid <- "MFP028";
tcga <- "CHOL";
task_max <-1;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep="")
);
adjust_pval_bh(mfdat);

