#-----#
plink_reg <- function(bedfile=NULL, bimfile=NULL, famfile=NULL, yFile=NULL, covFile=NULL,  covar=NULL, fids=NULL, ids=NULL, rsids=NULL, ,fnSNPs=NULL, fnOut=NULL, wd=NULL, logistic=FALSE,nthreads=1){
  #bedfile: path to bedfile
  #bimfile: path to bimfile
  #famfile: path to famfile
  #yFile: path to phenotype file
  #fids: family IDs
  #id: ID
  #rsids: vector of SNPs used in GWAS
  #fnSNPs: if you already have a list stored on disk
  #covFile: matrix containing fids, ids, and covariates used in GWAS
  #covar: list which covariables used in the analysis - should be in one character string; fx "sex pc1 pc2 pc3", and not "sex", "pc1", "pc2" "pc3"
  #fnOut: path+name of GWAS results
  #wd: working directory where the job should put temporary files
  #logistic: whether the regression should use the link function for logistic regression
  
  fnY <- paste0(wd,"y")
  fnIDs <- paste0(wd, "id")
  fnCOV <- paste0(wd,"cov")

  write.table( yFile, row.names = F, col.names = T, quote = F, sep = "\t", file = fnY  )  
  write.table( cbind(fids,ids), row.names = F, col.names = T, quote = F, sep = "\t", file = fnIDs  )  
  write.table( covFile, row.names = F, col.names = T, quote = F, sep = "\t", file = fnCOV )
  
  if(!is.null(rsids)){
      fnSNPs <- paste0(wd, "rsids")
      write.table( rsids, row.names = F, col.names = F, quote = F, sep = "\t", file = fnSNPs )
  }
  
  if(logistic==TRUE){
    system( paste("/home/rohde/software/plink2/plink2 --bfile", bedfile,
                "--fam", famfile,
                "--bim", bimfile,
                "--keep", fnIDs,
                "--extract", fnSNPs,   
                "--logistic hide-covar",
                "--keep-allele-order",
                "--covar", fnCOV, 
                "--covar-name", covar,
                "--pheno", fnY,
                "--silent",
                "--threads", nthreads,
                "--out", fnOut ) ) 
    }

  if(logistic==FALSE){
    system( paste("/home/rohde/software/plink2/plink2 --bfile", bedfile,
                "--fam", famfile,
                "--bim", bimfile,
                "--keep", fnIDs,
                "--extract", fnSNPs,  
                "--linear hide-covar",
                "--keep-allele-order",
                "--covar", fnCOV, 
                "--covar-name", covar,
                "--pheno", fnY,
                "--silent",
				"--threads", nthreads,
                "--out", fnOut ) ) 
    }       
}



plink_regOLD <- function(bedfile=NULL, bimfile=NULL, famfile=NULL, yFile=NULL, covFile=NULL,  covar=NULL, fids=NULL, ids=NULL, MAF=NULL,rsids=NULL, fnSNPs=NULL, fnOut=NULL, wd=NULL, logistic=FALSE){
  #bedfile: path to bedfile
  #bimfile: path to bimfile
  #famfile: path to famfile
  #yFile: path to phenotype file
  #fids: family IDs
  #id: ID
  #rsids: vector of SNPs used in GWAS
  #fnSNPs: if you already have a list stored on disk
  #covFile: matrix containing fids, ids, and covariates used in GWAS
  #covar: list which covariables used in the analysis - should be in one character string; fx "sex pc1 pc2 pc3", and not "sex", "pc1", "pc2" "pc3"
  #fnOut: path+name of GWAS results
  #wd: working directory where the job should put temporary files
  #logistic: whether the regression should use the link function for logistic regression
  
  fnY <- paste0(wd,"y")
  fnIDs <- paste0(wd, "id")
  fnCOV <- paste0(wd,"cov")

  write.table( yFile, row.names = F, col.names = T, quote = F, sep = "\t", file = fnY  )  
  write.table( cbind(fids,ids), row.names = F, col.names = F, quote = F, sep = "\t", file = fnIDs  )  
  write.table( covFile, row.names = F, col.names = T, quote = F, sep = "\t", file = fnCOV )
  
  if(!is.null(rsids)){
      fnSNPs <- paste0(wd, "rsids")
      write.table( rsids, row.names = F, col.names = F, quote = F, sep = "\t", file = fnSNPs )
  }
  
  if(logistic==TRUE){
    system( paste("plink --bfile", bedfile,
                "--fam", famfile,
                "--bim", bimfile,
                "--keep", fnIDs,
                "--extract", fnSNPs,   
                "--logistic hide-covar",
                "--keep-allele-order",
                "--covar", fnCOV, 
                "--covar-name", covar,
                "--pheno", fnY,
                "--maf",MAF,
                "--silent",
                "--out", fnOut ) ) 
    }

  if(logistic==FALSE){
    system( paste("plink --bfile", bedfile,
                "--fam", famfile,
                "--bim", bimfile,
                "--keep", fnIDs,
                "--extract", fnSNPs,  
                "--linear hide-covar",
                "--keep-allele-order",
                "--covar", fnCOV, 
                "--covar-name", covar,
                "--pheno", fnY,
				"--maf",MAF,
                "--silent",
                "--out", fnOut ) ) 
    }       
}
