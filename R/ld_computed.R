# Function to computed LD
#' computed linkage disequilibrium statistics for a locus of interest
#'
#' This function allows the user to extract linkage disequilibrium statistics between a variant of interest and surrounding variants within a contiguous genomic region. This function uses the University of Michigan LocusZoom API (<https://portaldev.sph.umich.edu/>) to obtain LD information, and allows the user to specify genome-build and ancestry of interest.
#'
#' @param chrom Integer - chromosome of reference variant
#' @param pos Integer - position of reference variant
#' @param ref Character - reference allele (or effect allele) for reference variant
#' @param alt Character - alternate allele (or non-effect allele) for reference variant
#' @param start Integer - starting position of range of interest
#' @param stop Integer - ending position of range of interest
#' @param metric Character - one of "r", "rsquare", or "cov", referring to the correlation statistic of interest
#' @param bfile Character - contained binary fileset plink.bed + plink.bim + plink.fam to be referenced
#' @param plink Character - contained binary plink executable
#'
#' @return A tibble containing each variant within the supplied range surrounding the variant of interest, with the requested linkage disequilibrium information with respect to the variant of interest
#' @export
#'
#' @examples
#' \dontrun{
#' ld_extract_locuszoom(chrom = 16, pos = 53830055, ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5,  metric = "rsquare")
#' }
#'

run_plink <-function(plink,args, outf=NULL){
	 suppressWarnings(
           out<- system2(command = normalizePath(plink),args = args,stdout = TRUE,stderr = TRUE)
         )
         if(!is.null(attributes(out))){
            stop(
	      paste("error plink : \n",
              paste0(out, collapse = " \n")," \n", "commands plink : ",
              normalizePath(plink), " ", paste0(args, collapse = " "), "\n"
            ))
         }
	 if(!is.null(outf) & !file.exists(outf)){
		         stop(
              paste("plink not created output file :",outf," \n",
              paste0(out, collapse = " \n")," \n", "commands plink : ",
              normalizePath(plink), " ", paste0(args, collapse = " "), "\n"
            ))


	 }
	 return(out)
}

ld_computed <- function(chrom, pos, ref, alt, start, stop, bfile,beddata=NULL,metric='r2',  plink='plink') {
          #--ld-snp rs12345
          #--ld-window-kb 1000
          #--ld-window 99999
          #--ld-window-r2 0
	 tmpdir=tempdir()
	 filebed=paste(tmpdir,'/',chrom,'_',start,'_',stop,'.pos.bed',sep='')
	 bfiletmp=paste(tmpdir,'/',chrom,'_',start,'_',stop,'',sep='')
         if(is.null(beddata))writeLines(paste(chrom,start,stop,start,sep='\t'),con=filebed)
	 else write.table(beddata[,c(1,2,3,2)], sep='\t',quote=F, row.names=F,col.names=F,file=filebed)
	 args=paste("-bfile ",bfile, '--extract range ',filebed,' --make-bed --keep-allele-order -out ',bfiletmp, sep=' ')
	 out<-run_plink(plink,args, paste(bfiletmp,'.bim',sep=''))
	 # read bim
         bimf<- paste(bfiletmp,'.bim',sep='')
	 databim<-read.table(bimf)
	 databim$V2<-paste(databim$V1,':',databim$V4,'_',databim$V6,'/',databim$V5,sep='')
	 write.table(databim, file=bimf, row.names=F,col.names=F, quote=F,sep='\t')
         rsid<-databim[databim$V1==chrom &databim$V4==pos, 'V2']
	 if(length(rsid)!=1){
            stop(
              paste("rsid of :",chrom,':',pos, 'not found in plink file',bfile," \n",
            ))
	 }
	 args_ld<-paste("-bfile ",bfiletmp,' -out ', bfiletmp, ' --ld-snp ', rsid,' --ld-window-r2 ',0,'--ld-window-kb 100000','--ld-window 99999 --r2')
         out_ld<-run_plink(plink,args_ld, paste(bfiletmp,'.ld',sep=''))
	 ld<-read.table(paste(bfiletmp,'.ld',sep=''),header=T);
	 ld<-ld[,c('SNP_A','CHR_A','BP_A','SNP_B','CHR_B','BP_B','R2')];names(ld)<-c('variant1','chromosome1','position1','variant2','chromosome2','position2', 'correlation')
	 return(ld)
}

