# mutfam_enrichment_tsv.R
# 28/02/2019
# MutFam using supplied files containing e.g. FunFams and mutation counts
# Allows flexibility regarding database sources

# Standard read format
#            read.table( test_file, stringsAsFactors = F, header = T, sep="\t", fill = TRUE, quote = "\"" ) );


# ------------------------------------------------------------------------------
# mf.init
#  source utils & db connect
# ------------------------------------------------------------------------------
mf.init <- function() {

  library( data.table );
  library( plyr );
  
}

# ------------------------------------------------------------------------------
# mf.main
# ------------------------------------------------------------------------------
mf.main <- function(cancer_type, version){
  if (missing(cancer_type)){
    cancer_type="PANCANCER";
  }
  if (missing(version)){
    version="vXXX";
  }

  # Get data.tables needed
  # ----------------------
  # mf: muts in funfams
  dt_mf <- mf.calc_mf(cancer_type);

  # ng: muts in all containing genes
  dt_ng <- mf.calc_ng( cancer_type);

  # fm: amino fractions of funfams/genes
  dt_fm <- mf.calc_aa_ratio(dbConn);

  # Create main data.table
  dt_me <- mf.calc_dt_me(dt_mf,dt_ng,dt_fm);

  # Remove entries <?? mutations total
  # Maybe depends on cancer_type?
##  dt_me_small <- dt_me[dt_me$COS_MUTS_FUNFAMS >20,];
  dt_me_small <- dt_me[dt_me$FF_MUTS >10,];
  # sort
  dt_me_small <- dt_me_small[order(rank(-ef))];
  setkey(dt_me_small,SUPERFAMILY_ID, FUNFAM_NUMBER);

  # FunFam info
  dt_sfffi <- mf.get_sf_ff_info(dbConn);

  # The combo report
  dt_combo <- dt_sfffi[dt_me_small];
  dt_combo <- dt_combo[order(rank(-ef))];
  # No point if ef<1!
  dt_combo <- dt_combo[ef>1];

  # 21/03/2016 - also remove where fm>75 (as per Miller/Sander)
  dt_combo <- dt_combo[fm_ratio<=75];

  # Sfam info
  dt_sfi <- mf.get_sf_info(dbConn);

  # The combo report
  setkey(dt_combo,SUPERFAMILY_ID);
  dt_combo <- dt_combo[dt_sfi];
  setnames(dt_combo,"i.NAME","Superfamily_Name");
  # Remove nulls
  dt_combo <- dt_combo[!is.na(dt_combo$ef),];
  dt_combo <- dt_combo[order(rank(-ef))];

}

# ------------------------------------------------------------------------------
# mf.calc_dt_me : Build data.table combining FunFam muts, Gene muts and
# mf - the fraction of aminos contained in FunFams compared to genes.
# ------------------------------------------------------------------------------
mf.calc_dt_me <- function(dt_mf, dt_ng, dt_fm){

  # Join FunFam and gene mutations and rename columns
  dt_me <- dt_ng[dt_mf];
#  setnames(dt_me, "TOT_MUTS_COSMIC", "COS_MUTS_GENES");
#  setnames(dt_me, "i.TOT_MUTS_COSMIC", "COS_MUTS_FUNFAMS")

  # Join the fraction of funfam / gene aminos dt
  dt_me <- dt_me[dt_fm];

  # Calculate FunFam COSMIC mutation enrichment score
  # ef = mf / (ng x fm)
  dt_me <- dt_me[,.(SUPERFAMILY_ID,
                    FUNFAM_NUMBER,
                    FF_MUTS,
                    GENE_MUTS_IN_FFS,
                    SUM_FF,
                    FF_SUM_GENE_LEN,
                    fm_ratio,
                    ef = FF_MUTS / (GENE_MUTS_IN_FFS * fm_ratio)
                    )];

  # REMOVE NAs
  dt_me$ef <- sapply (dt_me$ef, na<-function(i){if (is.na(i)){0}else{i}});

  # Delete rows where FF_MUTS is NULL (none found -> not relevant)
  dt_me <- dt_me[!is.na(dt_me$FF_MUTS),];

  return(dt_me);
}

  mf.test <- function() {
    # Testing principle with single FFs
    # test few ffs
    # ----------------------
    ff <- 7;
    dt_mf[FUNFAM_NUMBER==7];
    #dt_fms <-  function (dbConn);

  #   # Create db connection
  #   result = tryCatch({
  #     dbConn <-mfu.connectFGFRlocal()
  #   }, error = function(e) {
  #     print("Couldn't connect to Oracle database")
  #    break;
  #   }, finally = {
  #    # cleanup-code
  #   })

  }

# ------------------------------------------------------------------------------
# mf.calc_ef( ff_num) :
# Enrichment score of mutations (for a given database such as COSMIC)
# occurring in FunFam f.
# This is the entry function to calculate ef given Oracle db connection (dbConn)
# Passes dbConn alongside ff_num to relevant functions.
# ------------------------------------------------------------------------------
mf.calc_ef <- function ( ff_num){
	me <- mf.calc_me(ff_num);
	if (me >0){
		return(mf.calc_mf( ff_num) / me);
	}else {
		return(0);
	}
}

# ------------------------------------------------------------------------------
# mf.calc_me( ff_num):
# Simple entry function to calculate me by calls to mf.calc_fm and mf.calc_ng.
# ------------------------------------------------------------------------------
mf.calc_me <- function ( ff_num){
	return(mf.calc_ng( ff_num) * mf.calc_fm( ff_num));
}


# ------------------------------------------------------------------------------
# mf.calc_mf( [ff_num]):
# Optional ff_num to specify single FunFam
# Sums observed number of mutations occurring in FunFam(s)
# across all genes referenced in the COSMIC dataset.
# Returns data.table
# ------------------------------------------------------------------------------
  #19/10/15 Updated SQL joins on
  # 1) MATERIALIZED view;
  # 2) Joins on UniProt ID mapped UP ID, not COSMIC gene name
  # sql=paste("SELECT FUNFAM_NUMBER, TOT_MUTS_COSMIC FROM FGFR.VW_MF_MUTS_COS_FUNFAMS", sep="");

  # 30/10/15
  # Combine all/spec FunFam functions with optional ff_num argument

  # 11/11/15
  # Fixed issue with UP mapping join (to MOD_AC)
  # FGFR.MVW_MF_MUTS_COS_FUNFAMS_UP now returns site and gene info, so need to
  # add group by depending on cancer_type, if passed.

  # 17/11/2015: Correct bug re: FunFamIDs - need to always use composite of Sfam-FunFam!
  # Have removed FunFam WHERE at present.

  # 23/11/2015: Add germline mutations from humsavar (via Oracle/SQL)
  # Allow cancer_type="humvar"

  # 16/02/2016: Updates based on TBL_CANCER_LOOKUP and HGNC_ID changes.

  # 21/03/2016: v2.2 updates re: FunFam multi-region issue and correcting SQL SELECT

#mf.calc_mf <- function ( sf_num, ff_num, cancer_type){
mf.calc_mf <- function ( cancer_type){
  require (data.table);

  # # 12/02/16 Also allow Polymorphism (MV changed too!)
  if ((cancer_type == "HUMVAR") | (cancer_type == "HUMVAR_POLY")){
    # germline mutations (special case)
    sql <- "SELECT SUPERFAMILY_ID, FUNFAM_NUMBER, SUM(TOT_MUTS_HUMVAR) FF_MUTS FROM FGFR.MVW_MF_GERMUTS_FUNFAMS_UP ";
    # FGFR.MVW_GERMLINE_HUMVAR GH
  } else{
    # Cancer/COSMIC
    # 12/02/2016 HGNC_ID chg made in MVW
    sql <- "SELECT SUPERFAMILY_ID, FUNFAM_NUMBER, SUM(TOT_MUTS_COSMIC) FF_MUTS FROM FGFR.MVW_MF_MUTS_COS_FUNFAMS_UP ";
  }
  # # Swanton - special case - need better way of doing this!
  # if (cancer_type == "SWAMUTS"){
  #   # 12/02/2016 HGNC_ID chg needed - NOT DONE
  #   ##sql <- "SELECT SUPERFAMILY_ID, FUNFAM_NUMBER, SUM(TOT_MUTS_HUMVAR) FF_MUTS FROM FGFR.VW_MF_SWAMUTS_FUNFAMS_UP ";
  # }

  # Need to build appropriate WHERE to count muts from the materialised view
  # Get WHERE clause
  sql_get <- paste("SELECT WHERE_CLAUSE FROM FGFR.VW_CANCER_LOOKUP_CURRENT WHERE ABBREV='", cancer_type,"'",sep="");
  sql_where <- dbGetQuery( sql_get);
  if (nrow(sql_where) ==1 ){
    if (!is.na(sql_where)){
      sql <- paste(sql, " WHERE ",sql_where)
    }
  }else{
    return(-1)
  }

  # Get GROUP BY clause
  sql_get <- paste("SELECT GROUP_BY FROM FGFR.VW_CANCER_LOOKUP_CURRENT WHERE ABBREV='", cancer_type,"'",sep="");
  sql_group_by <- dbGetQuery( sql_get);
  if (nrow(sql_group_by) ==1 ){
    sql <- paste(sql, " GROUP BY ",sql_group_by);
  }else{
    return(-1)
  }

  # Finally, execute to get MF results!
  dt_mf <- data.table(dbGetQuery( sql));
  if (is.na(dt_mf)) {
    return(-1)
  } else
  {
    setkey(dt_mf, SUPERFAMILY_ID, FUNFAM_NUMBER);
    return(dt_mf)
  };


  # sql_where <- "";
  # if (missing(cancer_type)){
  #   # All cancer or germline don't have WHERE
  #   cancer_type <- "PANCANCER";
  # } else
  # {
  #   if (str_to_upper(cancer_type)=="GBM"){
  #     sql_where <- " WHERE PRIMARY_SITE='central_nervous_system' AND PRIMARY_HISTOLOGY='glioma' "
  #   }
  #   if (str_to_upper(cancer_type)=="BLADDER"){
  #     sql_where <- " WHERE PRIMARY_SITE='urinary_tract' AND SITE_SUBTYPE='bladder' "
  #   }
  #   if (str_to_upper(cancer_type)=="HUMVAR_POLY"){
  #     sql_where <- " WHERE PRIMARY_HISTOLOGY='Polymorphism' "
  #   }
  # }

  # Add WHERE clause if specific Sfam/FunFam
#   if (!missing(ff_num)){
#     if (sql_where==""){
#       sql_where <- paste(" WHERE FUNFAM_NUMBER = ", ff_num, sep="");
#     }else{
#       sql_where <- paste(sql_where, " AND FUNFAM_NUMBER = ", ff_num, sep="");
#     }
#   }

  # GROUP BY FunFam
#  sql <- paste(sql, sql_where, " GROUP BY SUPERFAMILY_ID, FUNFAM_NUMBER", sep="");


}
# ------------------------------------------------------------------------------
# mf.calc_ng( [ff_num]):
# Optional ff_num to specify single FunFam
# Sums observed number of muts occurring in all genes containing FunFams
# Returns data.table
# ------------------------------------------------------------------------------
# Updates 12/11/15: Now uses new materialised view with fixes.
# See EN: "MutFams - p-values & testing"

# 17/11/2015: Correct bug re: FunFamIDs - need to always use composite of Sfam-FunFam!
# Have removed FunFam WHERE at present.

# 23/11/2015: Add germline mutations from humsavar (via Oracle/SQL)
# Allow cancer_type="humvar"

# 12/02/2016 HGNC_ID chg made in MVWs and TBL_CANCER_LOOKUP useage

# 21/03/2016: v2.2 updates re: FunFam multi-region issue and correcting SQL SELECT

#mf.calc_ng <- function ( ff_num, cancer_type){
mf.calc_ng <- function ( cancer_type){
  require (data.table);

  # # 12/02/16 Also allow Polymorphism (MV changed too!)
  if ((cancer_type == "HUMVAR") | (cancer_type == "HUMVAR_POLY")){
    # germline mutations (special case)
    sql <- "SELECT SUPERFAMILY_ID, FUNFAM_NUMBER, SUM(TOT_MUTS_HUMVAR) GENE_MUTS_IN_FFS FROM FGFR.MVW_MF_GERMUTS_GENE_BY_FUNFAM ";
    # FGFR.MVW_GERMLINE_HUMVAR GH
  } else{
    # Cancer/COSMIC
    # 12/02/2016 HGNC_ID chg made in MVW
    sql <- "SELECT SUPERFAMILY_ID, FUNFAM_NUMBER, SUM(TOT_MUTS_COSMIC) GENE_MUTS_IN_FFS FROM FGFR.MVW_MF_MUTS_COS_GENE_BY_FUNFAM ";
  }

  if (sql==""){return(-1)}


  # Need to build appropriate WHERE to count muts from the caner lookup view
  sql_get <- paste("SELECT WHERE_CLAUSE FROM FGFR.VW_CANCER_LOOKUP_CURRENT WHERE ABBREV='", cancer_type,"'",sep="");
  sql_where <- dbGetQuery( sql_get);
  if (nrow(sql_where) == 1 ){
    if (!is.na(sql_where)){
      sql <- paste(sql, " WHERE ",sql_where)
    }
  }else{
    return(-1)
  }

  # Get GROUP BY clause
  sql_get <- paste("SELECT GROUP_BY FROM FGFR.VW_CANCER_LOOKUP_CURRENT WHERE ABBREV='", cancer_type,"'",sep="");
  sql_group_by <- dbGetQuery( sql_get);
  if (nrow(sql_group_by) ==1 ){
    sql <- paste(sql, " GROUP BY ",sql_group_by);
  }else{
    return(-1)
  }

  # Finally, execute to get MF results!
  dt_ng <- data.table(dbGetQuery( sql));
  if (is.na(dt_ng)) {
    return(-1)
  } else
  {
    setkey(dt_ng, SUPERFAMILY_ID, FUNFAM_NUMBER);
    return(dt_ng)
  };

}

# ------------------------------------------------------------------------------
# mf.calc_aa_ratio( ff_num):
# Calculate ratio of number of amino acids in FunFam to number of amino acids
# in all genes that contain the FunFam.
# ------------------------------------------------------------------------------
# 17/11/2015: Correct bug re: FunFamIDs - need to always use composite of Sfam-FunFam!
# Have removed FunFam WHERE at present.
#mf.calc_aa_ratio <- function ( ff_num){

# 21/03/2016: v2.2 updates re: FunFam multi-region issue and correcting SQL SELECT
mf.calc_aa_ratio <- function (dbConn){
 require (data.table);

  # What is the total amino acid length (across all genes) of a given FunFam?
  # 21/03/16 sql_ff <- paste("SELECT SUPERFAMILY_ID, FUNFAM_NUMBER, SUM_FF FROM FGFR.VW_MF_SUM_GENE_FF_LEN", sep="");
  sql_ff <- paste("SELECT SUPERFAMILY_ID, FUNFAM_NUMBER, SUM_FF FROM FGFR.MVW_MF_SUM_GENE_FF_LEN", sep="");

  # What is the total amino acid length of all genes that contain a given FunFam?
  # 21/03/16 sql_gene <- paste("SELECT SUPERFAMILY_ID, FUNFAM_NUMBER, FF_SUM_GENE_LEN FROM FGFR.VW_MF_SUM_GENE_LEN", sep="");
  sql_gene <- paste("SELECT SUPERFAMILY_ID, FUNFAM_NUMBER, FF_SUM_GENE_LEN FROM FGFR.MVW_MF_SUM_GENE_LEN", sep="");

 fun_sum <- data.table(dbGetQuery(sql_ff));
 if (is.na(fun_sum)) {return(0)};

 gene_sum <- data.table(dbGetQuery( sql_gene));
 if (is.na(gene_sum)) {return(0)};

#if (missing(ff_num)){
	# Return ratio of FunFam lengths to total genes containing those FunFam's lengths.
  # Join the gene_len and ff_len data tables
  setkey(fun_sum, SUPERFAMILY_ID,FUNFAM_NUMBER);
  setkey(gene_sum, SUPERFAMILY_ID,FUNFAM_NUMBER);
  dt_ffgene <- gene_sum[fun_sum];
  # Calculate the ratio of FunFam to total length of all FunFam containing genes
  dt_ffgene <- dt_ffgene[,.(SUPERFAMILY_ID,FUNFAM_NUMBER, SUM_FF, FF_SUM_GENE_LEN, fm_ratio=SUM_FF/FF_SUM_GENE_LEN)];
  setkey(dt_ffgene,SUPERFAMILY_ID, FUNFAM_NUMBER);

  return(dt_ffgene);
 # } else{
 #   if (gene_sum > 0){
  #      return (fun_sum/gene_sum);
   #   } else{
        return(0);
    #  }
#  }
}


# ------------------------------------------------------------------------------
# mf.get_ff_aa_range(dbConn):
# What are the amino start / end points for each FF?
# Used for calculating FF mutation enrichment p-values
# Includes other fields for help/validation!
# ------------------------------------------------------------------------------
# 17/11/15 Add sfam field and key

mf.get_ff_aa_range <- function (dbConn){
 require (data.table);

# 21/03/16 ranges don't work this way...
return("Needs updating")

}

# ------------------------------------------------------------------------------
# mf.get_sf_ff_info(dbConn):
# FunFam names
# ------------------------------------------------------------------------------
mf.get_sf_ff_info <- function (dbConn){
 require (data.table);
 sql <- paste("SELECT SUPERFAMILY_ID, FUNFAM_NUMBER, NAME FROM STAGING.STA_CATH_FUNFAM WHERE IMPORTID=(SELECT IMPORTID FROM FGFR.IMPORT_CURR_VER WHERE TYPE='CATH')", sep="");
 sf_ff_info <- data.table(dbGetQuery(sql));
 if (is.na(sf_ff_info)) {
   return(0)
   } else {
     setkey(sf_ff_info,SUPERFAMILY_ID,FUNFAM_NUMBER);
     return(sf_ff_info);
   }
}
# ------------------------------------------------------------------------------
# mf.get_sf_ff_info(dbConn):
# Superfamily names
# ------------------------------------------------------------------------------
mf.get_sf_info <- function (dbConn){
 require (data.table);
 sql <- paste("SELECT SUPERFAMILY_ID, NAME, FF_COUNT FROM STAGING.STA_CATH_SUPERFAMILY_SUMMARY WHERE IMPORTID=(SELECT IMPORTID FROM FGFR.IMPORT_CURR_VER WHERE TYPE='CATH')", sep="");
 sf_sf_info <- data.table(dbGetQuery(sql));
 if (is.na(sf_sf_info)) {
   return(0)
   } else {
     setkey(sf_sf_info,SUPERFAMILY_ID);
     return(sf_sf_info);
   }
}
