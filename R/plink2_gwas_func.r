#' Performs linear/logistic regression using PLINK 2.0 alpha
#'
#' @param bedfile Path to the PLINK bedfile
#' @param bimfile Path to the PLINK bimfile
#' @param famfile Path to the PLINK famfile
#' @param yFile path to phenotype file, which is a matrix containing fids, ids and the phenotypic values
#' @param covFile path to the covariate file, which is a matrix containing fids, ids, and the covariates used in the analysis
#' @param covar a character string of covaraites used in the analysis
#' @param fids family IDs
#' @param ids IDs
#' @param rsids vector of SNPs used in the association analysis
#' @param fnSNPs path to a list of SNPs if they already are stored on disk
#' @param fnOut path and file name of the association results
#' @param wd working directory of temporary files
#' @param logistic logical if TRUE the logic link function is used
#' @param firth logical if TRUE Firth's logistic regression will be performed.
#' @param ncores number of cores used in the analysis
#' @param dir path to PLINK 2.0 execuatble
#'
#' @export
plink_reg <- function(bedfile=NULL, bimfile=NULL, famfile=NULL, yFile=NULL, covFile=NULL,  covar=NULL, fids=NULL, ids=NULL,
                      rsids=NULL, fnSNPs=NULL, fnOut=NULL, wd=NULL, logistic=FALSE,firth=FALSE, ncores=1, dir=NULL){

  fnY <- paste0(wd,"y")
  fnIDs <- paste0(wd, "id")
  fnCOV <- paste0(wd,"cov")

  fwrite(yFile, file=fnY, sep="\t")
  fwrite(as.data.frame(cbind(fids,ids)), file=fnIDs, sep="\t")
  fwrite(covFile, file=fnCOV, sep="\t")

  if(!is.null(rsids)){
      fnSNPs <- paste0(wd, "rsids")
      fwrite( rsids, sep = "\t", file = fnSNPs )
  }

  if(logistic==TRUE){
    system( paste(dir, " --bfile", bedfile,
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
                "--threads", ncores,
                "--out", fnOut ) )
    }

  if(logistic==FALSE){
    system( paste(dir, " --bfile", bedfile,
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
				"--threads", ncores,
                "--out", fnOut ) )
    }
    if(firth==TRUE){
    system( paste(dir, " --bfile", bedfile,
                "--fam", famfile,
                "--bim", bimfile,
                "--keep", fnIDs,
                "--extract", fnSNPs,
                "--glm firth-fallback hide-covar",
                "--keep-allele-order",
                "--covar", fnCOV,
                "--covar-name", covar,
                "--pheno", fnY,
                "--silent",
                "--threads", ncores,
                "--out", fnOut ) )
    }
}



plink_regOLD <- function(bedfile=NULL, bimfile=NULL, famfile=NULL, yFile=NULL, covFile=NULL,  covar=NULL, fids=NULL, ids=NULL, rsids=NULL, fnSNPs=NULL, fnOut=NULL, wd=NULL, logistic=FALSE){
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
                "--silent",
                "--out", fnOut ) )
    }
}
