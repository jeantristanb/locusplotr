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
  } else if (genome_build == "GRCh37") {
    gwas_cat_db<- gc_hg37
  } else {
    stop("Invalid genome version. Use 'GRCh37' or 'GRCh38'.")
  }
  if(is.null(gwas_cat_db)){
  stop('gwas_cat_db is null ')
  }
  chromosome <- chr
  filter_start <- start
  filter_end <- end
  # Filter genes within the specified region
  gc <- gwas_cat_db %>%
    filter(chr == chro,
           start  <= bp,
           end  >= bp) %>%
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
  balise_plot1=F
  gc2<-gc[order(gc$bp),c('bp','label')]
  if(balise_plot1){
    balise<-T
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
   p <- ggplot(gc2, aes(x= bp)) + geom_segment(aes(x = bp, xend = bp, y=4.8,yend = 6), linewidth = 0.5, color = "darkblue") + geom_text(aes(x = bp, y = 4.7, label = label), angle = 45, size = 2.5, hjust='top') + ylim(1, 5) + labs(x = glue::glue("Position on Chromosome {chr} (Mb)"),y = paste(nbdel," hits omited")) +     scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6), limits = c(start, end)) + theme_bw(base_size = 16) +theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())
  }else{
  gc2$y<- (0:(nrow(gc2)-1)%%10)/2
  ymax=max(gc2$y)+1
   #p <- ggplot(gc2, aes(x= bp)) + ylim(0,5)+geom_label_repel(aes(x=bp,y=y,label=label),force=1,ylim=c(0,5),size=3, max.iter = 1e5, max.time = 4,direction= "y",  hjust=0.5,vjust=1,fill = NA)  + labs(x = glue::glue("Position on Chromosome {chr} (Mb)"),y = "") + geom_segment(aes(x = bp, xend = bp, y=y+0.2,yend = 5), linewidth = 0.5, color = "darkblue") +     scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6), limits = c(start, end)) + theme_bw(base_size = 16) +theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())
   p <- ggplot(gc2, aes(x= bp,y=y,label=label))+geom_label_repel(size=3, max.iter = 1e5, max.time = 4,direction= "y",  hjust=0.5,vjust=1,fill = NA,segment.alpha = 0, nudge_y = 0.4)  #
  #+ labs(x = glue::glue("Position on Chromosome {chr} (Mb)"),y = "") + geom_segment(aes(y=y-0.2,xend = bp,yend = 5),linewidth = .4,arrow=arrow()) +     scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6), limits = c(start, end)) + theme_bw(base_size = 16) +theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())

  # p is your ggplot object that already contains geom_label_repel()
built    <- ggplot_build(p)
# index of the ggrepel layer (works for geom_label_repel and geom_text_repel)
lbl_idx  <- which(vapply(p$layers,
                         \(l) inherits(l$geom, c("GeomLabelRepel","GeomTextRepel")),
                         logical(1)))

label_df <- built$data[[lbl_idx]]   # <- data‑frame with the new coordinates
  write.table(label_df,file='tmp')
  print(head(label_df))
  maxy<-max(label_df$y)
  p<-ggplot(label_df, aes(x= x,y=y,label=label))  + ylim(c(0,maxy+0.5))+geom_segment(aes(x=x,y=y+0.2,xend = x_orig,yend = maxy+0.5),linewidth = .4,arrow=arrow(type = "closed",length = unit(4, "pt")),size = 0.2,color='blue',alpha=0.5) + geom_label(size=2.5,fill = NA, color='black', label.size = NA)   + labs(x = glue::glue("Position on Chromosome {chr} (Mb)"),y = "") +     scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6), limits = c(start, end)) + theme_bw(base_size = 16) +theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())


  }
  return(p)
}
