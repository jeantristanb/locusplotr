library(data.table)
library(ggplot2)
source('gg_gc.R')
#install.packages('gginnards')
library(ggtext)
library(dplyr)
library(tidyr)
source('gg_locusplot.R')
source('gg_geneplot.R')
source('recomb_extract_locuszoom.R')
source('ld_extract_locuszoom.R')
#source('ld_extract_locuszoom')
#install.packages('snpsettest')

ld<-fread('~/locuszoom2_test/plink.ld');ld<-ld[,c('SNP_A','CHR_A','BP_A','SNP_B','CHR_B','BP_B','R2')]
#   variant1     chromosome1 position1 variant2 chromosome2 position2 correlation
names(ld)<-c('variant1','chromosome1','position1','variant2','chromosome2','position2', 'correlation')
gwas<-fread("~/locuszoom2_test/summstats.sorted.tab.gz")
rsid<-'MarkerName'
pos<-'Position'
chrom<-'Chromosome'
ref<-'Allele1'
alt<-'Allele2'
effect='Effect'
std_err='StdErr'
p_value='P.value'
genome_build='GRCh38'
#      select(rsid = {{ rsid }}, chromosome = {{ chrom }}, position = {{ pos }}, ref = {{ ref }}, alt = {{ alt }}, log10_pval) %>%
#Chromosome Allele1 Allele2
#~c(rs_col, chr_col, bp_col,a1_col,a2_col)
chr<-13
rg<-range(c(ld[ld[[1]]==chr,2], ld[ld[[1]]==chr,5]))
gwassub<-gwas %>% filter(Chromosome == chr, Position>=rg[1]-500000, Position <= rg[2]+500000) 
gc_db<-fread('full', sep='\t');names(gc_db)[names(gc_db)=='DISEASE/TRAIT']<-'DISEASE'
gc_db<-gc_db %>% mutate(rsid = SNPS , chro=as.integer(CHR_ID), bp=as.integer(CHR_POS), label=DISEASE) %>% select(rsid, chro,bp,label)%>%drop_na()

#gwas_sub<-gwassub %>% rename_with(~c('rsid', 'chrom', 'pos','ref','alt'), c(rs_col, chr_col, bp_col,a1_col,a2_col))
##' @param df Dataframe containing columns with rsid, chromosome, position, reference/effect allele, alternate/non-effect allele, and p-value for all variants within the range of interest
#   variant1     chromosome1 position1 variant2 chromosome2 position2 correlation


gg_locusplot(gwassub, lead_snp = NULL, rsid = rsid, chrom = chrom, pos = pos, ref = ref, alt = alt, effect = effect, std_err = std_err, p_value = p_value, trait = NULL, plot_pvalue_threshold = 0.5, plot_subsample_prop = 0.25, plot_distance = 500000, genome_build = genome_build, population = "ALL", plot_genes = T, plot_recombination = T, plot_title = NULL, plot_subtitle = NULL, path = NULL,ld_extracted=ld)
