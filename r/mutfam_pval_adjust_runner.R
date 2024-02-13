# Run p-val corrections

setwd("~/woofgit/fgfr-net-mutations/genome_wide/mutfams/R");
source('mutfam_pval_adjust.R');

#17/05/17
# GMFP05
group_dir <- "~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP05";
# MFP028 CHOL /I55
runid <- "MFP028";
tcga <- "CHOL";
task_max <-1;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep="")
);
adjust_pval_bh(mfdat);

# LIHC / I56 
runid <- "MFP029";
tcga <- "LIHC";
task_max <- 5;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep="")
);
adjust_pval_bh(mfdat);

# ESCA / I57
#tail -q -n +2 mutfam_enrichment_v2.2*.csv > mutfam_enrichment_v2.2_g3d_link_ESCA_PVAL_10000_BATCH_25_TASKID_1-2.csv
#cat ../mutfam_v2.2_HDR.txt mutfam_enrichment_v2.2_g3d_link_ESCA_PVAL_10000_BATCH_25_TASKID_1-2.csv > mutfam_enrichment_v2.2_g3d_link_ESCA_PVAL_10000_BATCH_25_TASKID_1-2_HDR.csv
group_dir <- "~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP05";
runid <- "MFP030";
tcga <- "ESCA";
task_max <- 2;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep=""));
adjust_pval_bh(mfdat);

# CESC / I58
#tail -q -n +2 mutfam_enrichment_v2.2*.csv > mutfam_enrichment_v2.2_g3d_link_CESC_PVAL_10000_BATCH_25_TASKID_1-1.csv
#cat ../mutfam_v2.2_HDR.txt mutfam_enrichment_v2.2_g3d_link_CESC_PVAL_10000_BATCH_25_TASKID_1-1.csv > mutfam_enrichment_v2.2_g3d_link_CESC_PVAL_10000_BATCH_25_TASKID_1-1_HDR.csv
group_dir <- "~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP05";
runid <- "MFP031";
tcga <- "CESC";
task_max <- 1;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep=""));
adjust_pval_bh(mfdat);

# PAAD / I59
#tail -q -n +2 mutfam_enrichment_v2.2*.csv > mutfam_enrichment_v2.2_g3d_link_PAAD_PVAL_10000_BATCH_25_TASKID_1-2.csv
#cat ../mutfam_v2.2_HDR.txt mutfam_enrichment_v2.2_g3d_link_PAAD_PVAL_10000_BATCH_25_TASKID_1-2.csv > mutfam_enrichment_v2.2_g3d_link_PAAD_PVAL_10000_BATCH_25_TASKID_1-2_HDR.csv
group_dir <- "~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP05";
runid <- "MFP032";
tcga <- "PAAD";
task_max <- 2;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep=""));
adjust_pval_bh(mfdat);

# UCEC / I60
#tail -q -n +2 mutfam_enrichment_v2.2*.csv > mutfam_enrichment_v2.2_g3d_link_UCEC_PVAL_10000_BATCH_25_TASKID_1-17.csv
#cat ../mutfam_v2.2_HDR.txt mutfam_enrichment_v2.2_g3d_link_UCEC_PVAL_10000_BATCH_25_TASKID_1-17.csv > mutfam_enrichment_v2.2_g3d_link_UCEC_PVAL_10000_BATCH_25_TASKID_1-17_HDR.csv
group_dir <- "~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP05";
runid <- "MFP033";
tcga <- "UCEC";
task_max <- 17;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep=""));
adjust_pval_bh(mfdat);

# THCA / I61
#tail -q -n +2 mutfam_enrichment_v2.2*.csv > mutfam_enrichment_v2.2_g3d_link_THCA_PVAL_10000_BATCH_25_TASKID_1-2.csv
#cat ../mutfam_v2.2_HDR.txt mutfam_enrichment_v2.2_g3d_link_THCA_PVAL_10000_BATCH_25_TASKID_1-2.csv > mutfam_enrichment_v2.2_g3d_link_THCA_PVAL_10000_BATCH_25_TASKID_1-2_HDR.csv
group_dir <- "~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP05";
runid <- "MFP034";
tcga <- "THCA";
task_max <- 2;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep=""));
adjust_pval_bh(mfdat);

# OV / I62
#tail -q -n +2 mutfam_enrichment_v2.2*.csv > mutfam_enrichment_v2.2_g3d_link_OV_PVAL_10000_BATCH_25_TASKID_1-2.csv
#cat ../mutfam_v2.2_HDR.txt mutfam_enrichment_v2.2_g3d_link_OV_PVAL_10000_BATCH_25_TASKID_1-2.csv > mutfam_enrichment_v2.2_g3d_link_OV_PVAL_10000_BATCH_25_TASKID_1-2_HDR.csv
group_dir <- "~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP05";
runid <- "MFP035";
tcga <- "OV";
task_max <- 2;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep=""));
adjust_pval_bh(mfdat);

# DLBC / I63
#tail -q -n +2 mutfam_enrichment_v2.2*.csv > mutfam_enrichment_v2.2_g3d_link_DLBC_PVAL_10000_BATCH_25_TASKID_1-1.csv
#cat ../mutfam_v2.2_HDR.txt mutfam_enrichment_v2.2_g3d_link_DLBC_PVAL_10000_BATCH_25_TASKID_1-1.csv > mutfam_enrichment_v2.2_g3d_link_DLBC_PVAL_10000_BATCH_25_TASKID_1-1_HDR.csv
group_dir <- "~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP05";
runid <- "MFP036";
tcga <- "DLBC";
task_max <- 1;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep=""));
adjust_pval_bh(mfdat);

# UCS / I64
#tail -q -n +2 mutfam_enrichment_v2.2*.csv > mutfam_enrichment_v2.2_g3d_link_UCS_PVAL_10000_BATCH_25_TASKID_1-1.csv
#cat ../mutfam_v2.2_HDR.txt mutfam_enrichment_v2.2_g3d_link_UCS_PVAL_10000_BATCH_25_TASKID_1-1.csv > mutfam_enrichment_v2.2_g3d_link_UCS_PVAL_10000_BATCH_25_TASKID_1-1_HDR.csv
group_dir <- "~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP05";
runid <- "MFP037";
tcga <- "UCS";
task_max <- 1;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_g3d_link_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep=""));
adjust_pval_bh(mfdat);


# ---------------------------------------------------------------------------------------------

# 05/04/17
# GMFP04
group_dir <- "~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP04";

# MFP028 CHOL
runid <- "MFP028";
tcga <- "CHOL";
task_max <-1;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep="")
                   );
adjust_pval_bh(mfdat);

# LIHC
runid <- "MFP029";
tcga <- "LIHC";
task_max <- 8;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep="")
);
adjust_pval_bh(mfdat);

# ESCA
runid <- "MFP030";
tcga <- "ESCA";
task_max <- 3;
mfdat <- file.path(group_dir, runid, 
                   paste("mutfam_enrichment_v2.2_",tcga,"_PVAL_10000_BATCH_25_TASKID_1-",task_max,"_HDR.csv", sep=""));
adjust_pval_bh(mfdat);

# *** ERROR FOUND IN ABOVE - SCRATCH GMFP04! *****  (IMPORTID missing from certain views refing CATH tables...)


# 23/03/16
humvar_poly_2.2 <- "/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP02/MFP010/mutfam_enrichment_v2.2_HUMVAR_POLY_PVAL_10000_BATCH_25_TASKID_1-3.csv";
adjust_pval_bh(humvar_poly_2.2);

humvar_2.2 <-"/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP02/MFP011/mutfam_enrichment_v2.2_HUMVAR_PVAL_10000_BATCH_25_TASKID_1-13.csv";
adjust_pval_bh(humvar_2.2);

GBM_2.2 <- "/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP02/MFP013/mutfam_enrichment_v2.2_GBM_PVAL_10000_BATCH_25_TASKID_1-2.csv"
adjust_pval_bh(GBM_2.2);

GLI_2.2 <- "/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP02/MFP015/mutfam_enrichment_v2.2_GLI_PVAL_10000_BATCH_25_TASKID_1-3.csv"
adjust_pval_bh(GLI_2.2);

# 04/04/16
BLCA_2.2 <- "/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP02/MFP016/mutfam_enrichment_v2.2_BLCA_PVAL_10000_BATCH_25_TASKID_1-3.csv"
adjust_pval_bh(BLCA_2.2);

# 05/04/16
# finally!
PANCAN_2.2 <- "/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP02/MFP012/mutfam_enrichment_v2.2_PANCANCER_PVAL_1000_BATCH_25_TASKID_1-249.csv"
adjust_pval_bh(PANCAN_2.2);

# GMFP03
#LUSC  Lung squamous cell carcinoma
#adjust_pval_bh('/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP018/mutfam_enrichment_v2.2_LUSC_PVAL_10000_BATCH_25_TASKID_1-6_HDR.csv');

# BRCA
#adjust_pval_bh('/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP019/mutfam_enrichment_v2.2_BRCA_PVAL_10000_BATCH_25_TASKID_1-6_HDR.csv')

# COAD
#adjust_pval_bh('/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP024/mutfam_enrichment_v2.2_COAD_PVAL_10000_BATCH_25_TASKID_1-17_HDR.csv')

# KIRC
#adjust_pval_bh('/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP020/mutfam_enrichment_v2.2_KIRC_PVAL_10000_BATCH_25_TASKID_1-3_HDR.csv')

# KIRP
#adjust_pval_bh('/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP021/mutfam_enrichment_v2.2_KIRP_PVAL_10000_BATCH_25_TASKID_1-1_HDR.csv')

#LAML
#adjust_pval_bh('/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP022/mutfam_enrichment_v2.2_LAML_PVAL_10000_BATCH_25_TASKID_1-1_HDR.csv')

#PRAD
#adjust_pval_bh('/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP023/mutfam_enrichment_v2.2_PRAD_PVAL_10000_BATCH_25_TASKID_1-1_HDR.csv')

# READ
#adjust_pval_bh('/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP025/mutfam_enrichment_v2.2_READ_PVAL_10000_BATCH_25_TASKID_1-3_HDR.csv')

# LUAD
#adjust_pval_bh('/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP017/mutfam_enrichment_v2.2_LUAD_PVAL_10000_BATCH_25_TASKID_1-17_HDR.csv')

# 23/05/16
# STAD
#adjust_pval_bh('/Users/ash/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP026/mutfam_enrichment_v2.2_STAD_PVAL_10000_BATCH_25_TASKID_1-12_HDR.csv')

# SKCM
#adjust_pval_bh('~/woofgit/fgfr-net-mutations/genome_wide/mutfams/run/GMFP03/MFP027/mutfam_enrichment_v2.2_SKCM_PVAL_10000_BATCH_25_TASKID_1-28_HDR.csv')




