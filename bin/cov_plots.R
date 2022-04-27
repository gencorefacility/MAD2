#!/usr/bin/Rscript
library(ggplot2)
library(glue)

args<-commandArgs(TRUE)
prefix = args[1]

##################### COVERAGE PLOTS #########################################
Covplot = function(df, savedir, coverage_cut, percentcov, prefix){
	  
	  xsamp = length(levels(factor(df$name)))
  
  #x = ggplot(data = df, aes(x=ntpos, y=totalcount, colour=(totalcount>=coverage_cut))) +
  x = ggplot(data = df, aes(x=ntpos, y=totalcount)) +
	      
	      geom_hline(yintercept = coverage_cut, linetype = 2) +
	          
	          geom_line(aes(group=1)) +
		      
		      #scale_color_manual(values=c('red','black')) +
		      
		      theme_bw() +
		          
		          theme(legend.key = element_blank(),
				          strip.background = element_rect(colour="black", fill="white"),
					            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    
    #facet_grid(name~segment, scales ='free') +
    facet_grid(name~segment) +
        
        xlab("Nucleotide Position") +
	    
	    ylab("Raw Read Depth")
      
      ggsave(x,
	              filename = glue("{savedir}/{prefix}.CoveragePlot.{coverage_cut}.{percentcov}.pdf"),
		               width = 12, height = xsamp*1.2, limitsize=FALSE)
}


# generate coverage files using the filtered mixed
Logcovplot = function(df, savedir, cov, percentcov, prefix){
	  
	  xsamp = length(levels(factor(df$name)))
  
  #x = ggplot(data = df, aes(x=ntpos, y=log10(totalcount), colour=(totalcount>=log10(cov)))) +
  x = ggplot(data = df, aes(x=ntpos, y=log10(totalcount))) +
	      
	      #geom_hline(yintercept = log10(cov), linetype = 2) +
	      
	      geom_hline(yintercept = log10(cov), linetype = 2) +
	          
	          geom_line(aes(group=1)) +
		      
		      #scale_color_manual(values=c('red','black')) +
		      
		      theme_bw() +
		          
		          theme(legend.key = element_blank(),
				          strip.background = element_rect(colour="black", fill="white"),
					            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    
    #facet_grid(name~segment, scales ='free') +
    facet_grid(name~segment) +
        
        xlab("Nucleotide Position") +
	    
	    ylab("Log10 Read Depth")
      
      #print(x)
      ggsave(x,
	              filename = glue("{savedir}/{prefix}.CoveragePlot.log10.{cov}.{percentcov}.pdf"),
		               width = 12, height = xsamp*1.2, limitsize=FALSE)
}

################################################################################

oldvars <- read.table(file = 'cov_data.tsv', sep = '\t', header =TRUE)
savedir = "./"
Covplot(oldvars,savedir,100,50,prefix)
# This would show coverage plots above the thresholds of 50% of the segment having at least 100X coverage
Logcovplot(oldvars,savedir,100,50,prefix)
# Same thresholds but log(coverage) instead of raw coverage
