library(tidyr)
library(dplyr)
library(data.table)
fileout=paste(tempdir(),'lz',sep='/')
download.file('https://www.ebi.ac.uk/gwas/api/search/downloads/alternative',fileout, method='wget')
gc_db<-fread(fileout,sep='\t')
gc_db<-gc_db %>% mutate(rsid = SNPS , chro=as.integer(CHR_ID), bp=as.integer(CHR_POS), label=MAPPED_TRAIT) %>% select(rsid, chro,bp,label)%>%drop_na()
gc_hg38<-gc_db
save(gc_hg38, file='../data/gc_hg38.rda')


