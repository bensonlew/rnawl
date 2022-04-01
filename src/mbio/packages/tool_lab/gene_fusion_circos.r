#!/usr/bin/env Rscript

# load package
library(getopt)
library(RCircos)


# set options
command <- matrix(c(
  "chr_infos", "c", 1, "character", "input file containing chr_infos data",
  "gene_infos", "g", 1, "character", "input file containing gene_infos data",
  "fusion_infos", "f", 1, "character", "input file containing fusion_infos data",
  "config", "e", 1, "character", "plot config info ",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

print("start")
# prepare data
chr_df = read.table(opts$chr_infos, sep='\t', header=TRUE)
gene_df = read.table(opts$gene_infos, sep='\t', header=TRUE)
fusion_df = read.table(opts$fusion_infos, sep='\t', header=TRUE)

config = read.table(opts$config, sep='\t', header=TRUE)
drop_chrs = config$drop_chrs

# if (drop_chrs == "None" ){
#    f_drop_chrs = NULL
#    }else {
#    f_drop_chrs = strsplit(drop_chrs,",")
#    }
if (drop_chrs[1] == "None" ){
   f_drop_chrs = NULL
   }else {
   f_drop_chrs = drop_chrs
   }


# plot
print("start plot outside circ")
chr.exclude <- f_drop_chrs
cyto.info <- chr_df
tracks.inside <- 10
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
RCircos.List.Plot.Parameters()
pdf(file="gene_fusion_circos.pdf", height=8, width=8, compress=TRUE)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

print("start plot gene line")
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(gene_df, track.num, side)
name.col <- 4
track.num <- 2
RCircos.Gene.Name.Plot(gene_df, name.col,track.num, side)
print("start plot fusion line")
track.num <- 5
RCircos.Link.Plot(fusion_df, track.num, TRUE)

dev.off()