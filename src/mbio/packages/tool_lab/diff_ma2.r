#!/usr/bin/env Rscript

# load package
library(getopt)
library(ggplot2)


# set options
command <- matrix(c(
  "data", "i", 1, "character", "input file containing processed data",
   "config", "c", 1, "character", "plot config info ",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

# prepare data
input = read.table(opts$data, sep='\t', header=TRUE)
config = read.table(opts$config, sep='\t', header=TRUE)
method =  config$m
title =  config$t
x_axis_title = config$x
y_axis_title = config$y


# calculate
colo = c("blue", "grey","red")
if (method == 1 ){
   colo = c("down" = "blue","normal" = "grey","up" = "red")
   }else if (method == 2){
   colo = c("down" = "green", "normal" = "grey","up" = "red")
   }else if (method == 3){
   colo = c("down" = "blue","normal" =  "black","up" = "red")
   }else if(method == 4){
   colo = c("down" = "green","normal" =  "black","up" = "red")
   }

pdf(file = 'diff_ma.pdf',height = 5,width = 5)
ggplot(data = input,aes(x = log10fpkm,y=log2fc,colour=regulate))+ ggtitle(title) +
 geom_point( size=2) + scale_color_manual(values=colo)+
 xlab(x_axis_title) + ylab(y_axis_title) + theme_bw() +
 theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
 legend.title = element_text(size = 15),legend.text = element_text(size = 12),
 legend.background = element_rect(fill = NULL, colour = "black"))
dev.off()

png(file = 'diff_ma.png')
ggplot(data = input,aes(x = log10fpkm,y=log2fc,colour=regulate))+ ggtitle(title) +
 geom_point( size=2) + scale_color_manual(values=colo)+
 xlab(x_axis_title) + ylab(y_axis_title) + theme_bw() +
 theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
 legend.title = element_text(size = 15),legend.text = element_text(size = 12),
 legend.background = element_rect(fill = NULL, colour = "black"))
dev.off()

svg(file = 'diff_ma.svg',height = 5,width = 5)
ggplot(data = input,aes(x = log10fpkm,y=log2fc,colour=regulate))+ ggtitle(title) +
 geom_point( size=2) + scale_color_manual(values=colo)+
 xlab(x_axis_title) + ylab(y_axis_title) + theme_bw() +
 theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
 legend.title = element_text(size = 15),legend.text = element_text(size = 12),
 legend.background = element_rect(fill = NULL, colour = "black"))
dev.off()




# pdf(file = 'diff_ma.pdf',height = 5,width = 5)
# ggplot(data = input,aes(x = log10fpkm,y=log2fc,colour=regulate))+ ggtitle(title) +
#  geom_point( size=2) + scale_color_manual(values=colo)+
#  xlab(x_axis_title) + ylab(y_axis_title) + theme_bw() +
#  theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
#  legend.title = element_text(size = 15),legend.text = element_text(size = 12),
#  legend.position = c(1, 0.7),legend.justification = c(1, 1),legend.background = element_rect(fill = NULL, colour = "black"))
# dev.off()
#
# png(file = 'diff_ma.png')
# ggplot(data = input,aes(x = log10fpkm,y=log2fc,colour=regulate))+ ggtitle(title) +
#  geom_point( size=2) + scale_color_manual(values=colo)+
#  xlab(x_axis_title) + ylab(y_axis_title) + theme_bw() +
#  theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
#  legend.title = element_text(size = 15),legend.text = element_text(size = 12),
#  legend.position = c(1, 0.7),legend.justification = c(1, 1),legend.background = element_rect(fill = NULL, colour = "black"))
# dev.off()
#
# svg(file = 'diff_ma.svg',height = 5,width = 5)
# ggplot(data = input,aes(x = log10fpkm,y=log2fc,colour=regulate))+ ggtitle(title) +
#  geom_point( size=2) + scale_color_manual(values=colo)+
#  xlab(x_axis_title) + ylab(y_axis_title) + theme_bw() +
#  theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
#  legend.title = element_text(size = 15),legend.text = element_text(size = 12),
#  legend.position = c(1, 0.7),legend.justification = c(1, 1),legend.background = element_rect(fill = NULL, colour = "black"))
# dev.off()
