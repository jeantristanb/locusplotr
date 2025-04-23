#!/usr/bin/env Rscript

#install.packages('remotes','optparse','data.table', 'dplyr','tidyr')
#remotes::install_github('jeantristanb/locusplotr',ref='test')
suppressMessages(library('locusplotr'))
suppressMessages(library("optparse"))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))

option_list = list(
  make_option(c("--sumstat"), type="character", default=NULL,
              help="sumstat file", metavar="character"),
  make_option(c( "--head_chr"), type="character", default=NULL,
              help="chromosome header of sumstat", metavar="character"),
  make_option(c( "--head_bp"), type="character", default=NULL,
              help="position header of sumstat", metavar="character"),
  make_option(c( "--head_pval"), type="character", default=NULL,
              help="pvalue header of sumstat", metavar="character"),
  make_option(c( "--head_af"), type="character", default=NULL,
              help="pvalue header of sumstat", metavar="character"),
  make_option(c( "--head_rs"), type="character", default=NULL,
              help="pvalue header of sumstat", metavar="character"),
  make_option(c( "--head_ref"), type="character", default=NULL,
              help="ref/effect header of sumstat", metavar="character"),
  make_option(c( "--head_alt"), type="character", default=NULL,
              help="alt/non effect header of sumstat", metavar="character"),
  make_option(c( "--head_beta"), type="character", default=NULL,
              help="alt/non effect header of sumstat", metavar="character"),
  make_option(c( "--head_se"), type="character", default=NULL,
              help="alt/non effect header of sumstat", metavar="character"),
  make_option(c( "--genome_build"), type="character", default="GRCh38",
              help="genome build : GRCh38 / GRCh37", metavar="character"),
  make_option(c( "--bfile"), type="character", default=NULL,
              help="head of bfile if need to compute ld", metavar="character"),
  make_option(c( "--wind_kb"), type="double", default=250,
              help="windows in kb", metavar="double"),
  make_option(c( "--bp"), type="integer", default=NULL,
              help="windows in kb or begin - end", metavar="integer"),
  make_option(c( "--begin"), type="double", default=NULL,
              help="windows in kb or begin - end", metavar="double"),
  make_option(c( "--end"), type="double", default=NULL,
              help="windows in kb or begin - end", metavar="double`"),
  make_option(c( "--chr"), type="character", default=NULL,
              help="chr", metavar="character"),
  make_option(c( "--maf"), type="double", default=0,
              help="chr", metavar="character"),
  make_option(c( "--ld_pop_ref"), type="character", default="ALL",
              help="chr", metavar="character"),
  make_option(c("--out"), type="character",default='out',
              help="output header", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

rs_col<-opt[['head_rs']]
bp_col<-opt[['head_bp']]
chr_col<-opt[['head_chr']]
ref_col<-opt[['head_ref']]
alt_col<-opt[['head_alt']]
effect_col=opt[['head_beta']]
stderr_col=opt[['head_se']]
p_col=opt[['head_p']]
genome_build=opt[['genome_build']]
around<-opt[['wind_kb']] * 1000

chr=opt[['chr']]
bp=opt[['bp']]
gwas<-fread(opt[['sumstat']])
bfile=opt[['bfile']]

compute_ld=T
out=opt[['out']]

if(!is.null(bp)){
gwassub<-gwas %>%  rename_with(~c('rsid', 'chrom', 'pos'), all_of(c(rs_col, chr_col, bp_col))) %>% filter(chrom== chr, pos>=bp-around, pos<= bp+around)
gwas_rs<-gwassub %>% filter(pos == bp) 
lead_snp = gwas_rs[,rsid]
}else{
begin <- opt[['begin']]
end <- opt[['end']]
gwassub<-gwas %>%  rename_with(~c('rsid', 'chrom', 'pos'), all_of(c(rs_col, chr_col, bp_col))) %>% filter(chrom== chr, pos>=begin, pos<= end)
lead_snp=NULL
}
population = NULL
if(is.null(bfile)){
compute_ld=F
population=opt[['ld_pop_ref']]
}
if(!is.null(opt[['head_af']]) & !is.null(opt[['maf']]) & opt[['maf']]>0){
	maf=opt[['maf']]
gwassub<-gwassub %>% rename_with(~c('af'),all_of(opt[['head_af']])) %>% filter(af>maf , af<1-maf)
}
write.table(gwassub, file='all.tsv', row.names=F, col.name=T,sep='\t',quote=F)

suppressWarnings(gg_locusplot(gwassub, lead_snp = lead_snp, rsid = 'rsid', chrom = 'chrom', pos = 'pos', ref = ref_col, alt = alt_col, effect = effect_col, std_err = stderr_col, p_value = p_col, trait = NULL, plot_pvalue_threshold = 1, plot_subsample_prop = 0.25, plot_distance = around, genome_build = genome_build, population = population, plot_genes = T, plot_recombination = T, plot_title = NULL, plot_subtitle = NULL, path = out,ld_extracted=NULL,plot_gc=T,gc_db=NULL,compute_ld=compute_ld,bfile=bfile))


