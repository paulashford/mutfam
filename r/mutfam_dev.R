# mutfam_dev_2019.R
# Scratch / dev area

# [0] Init
  source('mutfam_enrichment_tsv.R');
  source('mutfam_enrichment.R')
  
  library( data.table );
  library( plyr );

# [1] FGFR Oracle conn
  #ssh -N -L1521:gene3d01:1521 gene3d01 -v
  library(RJDBC);
  drv <- JDBC("oracle.jdbc.OracleDriver",classPath = "../db/ojdbc6.jar");
  con <- dbConnect(drv, "jdbc:oracle:thin:@localhost:1521:FGFR","XXXXXX","XXXXXXXX");

# [2] Get MutFams (LUAD) using Oracle db (old data!)
    cancer_type <- "LUAD";
    dbConn <- con;

    mf_data_dir <- "~/woofgit/mutfam/data/test/mf4_0";
    mf_cancer_dir <- file.path( mf_data_dir, cancer_type );

  # (a) Muts in FunFam and write out tsv
    dt_mf_luad <- mf.calc_mf( con, cancer_type ); 
    write.table( dt_mf_luad,
                file = file.path( mf_cancer_dir, "mf.tsv" ),
                row.names = FALSE, 
                quote = FALSE, 
                sep = "\t", 
                qmethod = "escape" ); # fileEncoding = "UTF-8")

  # (b) Muts in overall gene
    dt_ng_luad <- mf.calc_ng( con, cancer_type ); 
    write.table( dt_ng_luad,
                file = file.path( mf_cancer_dir, "ng.tsv" ),
                row.names = FALSE, 
                quote = FALSE, 
                sep = "\t", 
                qmethod = "escape" ); # fileEncoding = "UTF-8")

  
# [3] Get 'auxilliary' datasets
  # (a) fm: amino fractions of funfams/genes
    dt_fm <- mf.calc_aa_ratio( con ); 
    write.table( dt_fm,
                file = file.path( mf_data_dir, "dt_fm.tsv" ),
                row.names = FALSE, 
                quote = FALSE, 
                sep = "\t", 
                qmethod = "escape" ); # fileEncoding = "UTF-8")

  # (b) FunFam info
    dt_sfffi <- mf.get_sf_ff_info( dbConn );
    write.table( dt_sfffi,
                file = file.path( mf_data_dir, "dt_sfffi.tsv" ),
                row.names = FALSE, 
                quote = FALSE, 
                sep = "\t", 
                qmethod = "escape" ); # fileEncoding = "UTF-8")

  # (c) Sfam info
    dt_sfi <- mf.get_sf_info( dbConn );
    write.table( dt_sfi,
                file = file.path( mf_data_dir, "dt_sfi.tsv" ),
                row.names = FALSE, 
                quote = FALSE, 
                sep = "\t", 
                qmethod = "escape" ); # fileEncoding = "UTF-8")

# [4] Process MutFam data
  dt_me <- mf.calc_dt_me( dt_mf_luad, dt_ng_luad, dt_fm );
  # Remove entries <?? mutations total
  dt_me_small <- dt_me[ dt_me$FF_MUTS > 10, ];
  # sort
  dt_me_small <- dt_me_small[ order( rank( -ef ) ) ];

  # The combo report (with FunFam info)
  setkey( dt_me_small, SUPERFAMILY_ID, FUNFAM_NUMBER );
  setkey( dt_sfffi, SUPERFAMILY_ID, FUNFAM_NUMBER );

  dt_combo <- dt_sfffi[ dt_me_small ];
  dt_combo <- dt_combo[ order( rank( -ef ) ) ];
  # No point if ef<1!
  dt_combo <- dt_combo[ ef>1 ];

  # 21/03/2016 - also remove where fm>75 (as per Miller/Sander)
  dt_combo <- dt_combo[ fm_ratio <= 75 ];

  # The combo report
  setkey( dt_combo, SUPERFAMILY_ID );
  setkey( dt_sfi, SUPERFAMILY_ID );
  dt_combo <- dt_combo[ dt_sfi ];
  setnames( dt_combo, "i.NAME", "Superfamily_Name" );
  # Remove nulls
  dt_combo <- dt_combo[ !is.na( dt_combo$ef ), ];
  dt_combo <- dt_combo[ order( rank( -ef ) ) ];

