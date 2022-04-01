#!/usr/bin/env Rscript

# load package
library(getopt)
library(ggplot2)


# set options
command <- matrix(c(
  "data", "i", 1, "character", "input file containing processed data",
  "x_axis_name", "x", 1, "character", "x_axis_name",
  "y_axis_name", "y", 1, "character", "y_axis_name",
  "method", "m", 1, "integer", "color schemes",
  "title_name", "t", 1, "character", "title",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

# prepare data
input = read.table(opts$data, sep='\t', header=TRUE)


# calculate
colo = c("blue", "grey","red")
if (opts$method == 1 ){
   colo = c("down" = "blue","no" = "grey","up" = "red")
   }else if (opts$method == 2){
   colo = c("down" = "green", "no" = "grey","up" = "red")
   }else if (opts$method == 3){
   colo = c("down" = "blue","no" =  "black","up" = "red")
   }else if(opts$method == 4){
   colo = c("down" = "green","no" =  "black","up" = "red")
   }

pdf(file = 'diff_ma.pdf',height = 5,width = 5)
ggplot(data = input,aes(x = log10fpkm,y=log2fc,colour=regulate))+ ggtitle(opts$title_name) +
 geom_point( size=2) + scale_color_manual(values=colo)+
 xlab(opts$x_axis_name) + ylab(opts$y_axis_name) + theme_bw() +
 theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
 legend.title = element_text(size = 15),legend.text = element_text(size = 12),
 legend.position = c(1, 0.7),legend.justification = c(1, 1),legend.background = element_rect(fill = NULL, colour = "black"))
dev.off()


