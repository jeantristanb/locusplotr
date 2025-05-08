library(tidyr)
library(dplyr)
library(data.table)
fileout=paste(tempdir(),'lz_38',sep='/')
update_gc38=T
update_gc37=T
if(update_gc38){
download.file('https://www.ebi.ac.uk/gwas/api/search/downloads/alternative',fileout, method='wget')
gc_db<-fread(fileout,sep='\t')
gc_db<-gc_db %>% mutate(rsid = SNPS , chro=CHR_ID, bp=as.integer(CHR_POS), label=MAPPED_TRAIT) %>% select(rsid, chro,bp,label)%>%drop_na()
gc_hg38<-gc_db
save(gc_hg38, file='../data/gc_hg38.rda')
}

fileout=paste(tempdir(),'lz_37.gz',sep='/')
if(update_gc37){
download.file('http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz',fileout, method='wget')
gc_db<-fread(fileout,sep='\t')
names(gc_db)<-c("bin", "chrom", "chromStart", "chromEnd", "name", "pubMedID", "author", "pubDate", "journal", "title", "trait", "initSample","replSample","region", "genes", "riskAllele", "riskAlFreq", "pValue", "pValueDesc", "orOrBeta", "ci95","platform", "cnv")
gc_db<-gc_db %>% mutate(rsid = gsub('-[A-Za-z?]+','',riskAllele), chro=gsub('chr','',chrom), bp=as.integer(chromEnd), label=trait) %>% select(rsid, chro,bp,label)%>%drop_na()
gc_hg37<-unique(gc_db)
save(gc_hg37, file='../data/gc_hg37.rda')
}

