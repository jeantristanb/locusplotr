#' Plot gwas catalog located within a genomic region of interest
#'
#' Returns a ggplot containing the gwas catalog within a specified genomic region. The function uses database gave by user or gwas catalog result compile on the library
#'
#' @param chr Integer - chromosome
#' @param start Integer - starting position for region of interest
#' @param end Integer - ending position for region of interest
#' @param genome_build Character - genome build - one of "GRCh37" or "GRCh38"
#' @param gwas_cat_db data.frame - contained gwas catalog hit
#'
#' @return A ggplot object containing a plot of gwas catalog within the region of interest
#' @export
#'
#' @examples
#' \dontrun{
#' gg_gc(1, 170054349 - 1e6, 170054349 + 1e6, "GRCh38")
#' }
#'

gg_gc<- function(chr, start, end, genome_build = "GRCh38", gwas_cat_db=NULL) {
  checkmate::assert_numeric(chr)
  checkmate::assert_numeric(start)
  checkmate::assert_numeric(end)
  checkmate::assert_data_frame(gwas_cat_db,null.ok=T)
  checkmate::assert_choice(genome_build, choices = c("GRCh37", "GRCh38"))

  # Select the appropriate gene table based on the genome version
  if (genome_build == "GRCh38" & is.null(gwas_cat_db)) {
    gwas_cat_db<- gc_hg38
  } #else if (genome_build == "GRCh37") {
    #gene_table <- snpsettest::gene.curated.GRCh37
  #} else {
  #  stop("Invalid genome version. Use 'GRCh37' or 'GRCh38'.")
  #}
  if(is.null(gwas_cat_db)){
  stop('gwas_cat_db is null ')
  }
  chromosome <- chr
  filter_start <- start
  filter_end <- end
  # Filter genes within the specified region
  around<-(end - start)*0.05
  gc <- gwas_cat_db %>%
    filter(chr == chro,
           start +around <= bp,
           end -around >= bp) %>%
    select(rsid, bp,  label)


  # Check if any genes were found
  if (nrow(gc) == 0) {
    warning("No gwas catalog found in the specified region.")
    return(NULL)
  }

  # Assign y-levels to genes
  # Create the plot
  #  geom_text(aes(x = (start + end) / 2, label = gene), vjust = -0.5, size = 3) +
  gc<-unique(gc[,c('bp','label')])
  balise<-T
  gc2<-gc[order(gc$bp),c('bp','label')]
  nbdel=0
  while(any(balise)){
    tmp<-abs(gc2$bp[1:(nrow(gc2)-1)] - gc2$bp[2:(nrow(gc2))])/(end-start)
    balise<-tmp<0.02
    if(any(balise)){
	    pos_del<-which(balise)+1
	    nbdel=nbdel+1
            gc2<-gc2[-pos_del,]
    }
  }
  p <- ggplot(gc2, aes(x= bp)) + geom_segment(aes(x = bp, xend = bp, y=4.8,yend = 5), linewidth = 0.5, color = "darkblue") + geom_text(aes(x = bp, y = 4.7, label = label), angle = 45, size = 2.5, hjust='top') + ylim(1, 5) + labs(x = glue::glue("Position on Chromosome {chr} (Mb)"),y = paste(nbdel," hits omited")) +     scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6), limits = c(start, end)) + theme_bw(base_size = 16) +theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())
  return(p)
}
