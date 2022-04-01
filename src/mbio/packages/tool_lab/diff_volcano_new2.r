#!/usr/bin/env Rscript

# ~/app/bioinfo/rna/miniconda2/bin/R
# install.packages("ggrepel")
# load package
library(getopt)
library(ggplot2)
library(ggrepel)


# set options
command <- matrix(c(
  "data", "i", 1, "character", "input file containing processed data",
  "config", "c", 1, "character", "plot config info ",
  "showname_num", "n", 2, "character", "show the names of the first N genes",
  "showname_str", "s", 2, "character", "show the names of the first genes names.Seperated by ','",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

print("start")
# prepare data
input = read.table(opts$data, sep='\t', header=TRUE,colClasses=c("character","numeric","numeric"))
config = read.table(opts$config, sep='\t', header=TRUE)
#
show_id = list()
if (!is.null(opts$showname_str)){
   show_id = c(show_id,as.list(strsplit(opts$showname_str, ',')[[1]]))
}
if (!is.null(opts$showname_num)){
   input = input[sort(input$padjust,index.return=TRUE)$ix,]
   show_id = c(show_id,as.list(input[1:as.numeric(opts$showname_num),1]))
}
print(show_id)
#
log2fc_cutoff = log2(config$l)
method =  config$m
padjust_cutoff = config$p
title =  config$t
x_axis_title = config$x
y_axis_title = config$y

# calculate
total_num = nrow(input)
up_data<- input[input$padjust < padjust_cutoff & input$log2fc > log2fc_cutoff, ]
up_num = nrow(up_data)
down_data<- input[input$padjust < padjust_cutoff & input$log2fc < -log2fc_cutoff, ]
down_num = nrow(down_data)
up_name = paste("up:",as.character(up_num))
down_name = paste("down:",as.character(down_num))
no_num = (total_num - up_num - down_num )
no_name = paste("no_diff:",as.character(total_num - up_num - down_num ))
colo = c("blue", "grey","red")
if (method == 1 & up_num!= 0  & down_num != 0 &  no_num != 0 ){
   colo = c("blue","grey","red")
   }else if (method == 2 & up_num!= 0  & down_num != 0 &  no_num != 0 ){
   colo = c("green", "grey","red")
   }else if (method == 3 & up_num!= 0  & down_num != 0 &  no_num != 0 ){
   colo = c("blue", "black","red")
   }else if(method == 4 & up_num!= 0  & down_num != 0 &  no_num != 0 ){
   colo = c("green", "black","red")
   }

if (method == 1 & up_num == 0  & down_num != 0 &  no_num != 0 ){
   colo = c("blue","grey")
   }else if (method == 2 & up_num == 0  & down_num != 0 &  no_num != 0 ){
   colo = c("green", "grey")
   }else if (method == 3 & up_num == 0  & down_num != 0 &  no_num != 0 ){
   colo = c("blue", "black")
   }else if(method == 4 & up_num == 0  & down_num != 0 &  no_num != 0 ){
   colo = c("green", "black")
   }

if (method == 1 & up_num != 0  & down_num == 0 &  no_num != 0 ){
   colo = c("grey","red")
   }else if (method == 2 & up_num != 0  & down_num == 0 &  no_num != 0 ){
   colo = c("grey","red")
   }else if (method == 3 & up_num != 0  & down_num == 0 &  no_num != 0 ){
   colo = c( "black","red")
   }else if(method == 4 &  up_num != 0  & down_num == 0 &  no_num != 0){
   colo = c("black","red")
   }

if (method == 1 & up_num != 0  & down_num != 0 &  no_num == 0 ){
   colo = c("blue","red")
   }else if (method == 2 & up_num != 0  & down_num != 0 &  no_num == 0 ){
   colo = c("green", "red")
   }else if (method == 3 & up_num != 0  & down_num != 0 &  no_num == 0 ){
   colo = c("blue","red")
   }else if(method == 4 & up_num != 0  & down_num != 0 &  no_num == 0){
   colo = c("green","red")
   }

if (method == 1 & up_num!= 0  & down_num == 0 &  no_num == 0 ){
   colo = c("red")
   }else if (method == 2 & up_num!= 0  & down_num == 0 &  no_num == 0  ){
   colo = c("red")
   }else if (method == 3 & up_num!= 0  & down_num == 0 &  no_num == 0  ){
   colo = c("red")
   }else if(method == 4 & up_num!= 0  & down_num == 0 &  no_num == 0  ){
   colo = c("red")
   }

if (method == 1 & up_num == 0  & down_num != 0 &  no_num == 0 ){
   colo = c("blue")
   }else if (method == 2 & up_num == 0  & down_num != 0 &  no_num == 0  ){
   colo = c("green")
   }else if (method == 3 & up_num == 0  & down_num != 0 &  no_num == 0  ){
   colo = c("blue")
   }else if(method == 4 & up_num == 0  & down_num != 0 &  no_num == 0  ){
   colo = c("green")
   }

if (method == 1 & up_num == 0  & down_num == 0 &  no_num != 0 ){
   colo = c("grey")
   }else if (method == 2  & up_num == 0  & down_num == 0 &  no_num != 0  ){
   colo = c("grey")
   }else if (method == 3  & up_num == 0  & down_num == 0 &  no_num != 0 ){
   colo = c("black")
   }else if(method == 4 & up_num == 0  & down_num == 0 &  no_num != 0 ){
   colo = c("black")
   }



the <- max(abs(range(input$log2fc)))
the <- ceiling(the)

input$Significant = as.factor(ifelse (input$padjust < padjust_cutoff & abs(input$log2fc)>= log2fc_cutoff,ifelse(input$log2fc > log2fc_cutoff,up_name,down_name),no_name))

# pdf(file = 'volcano.pdf',height = 5,width = 5)
# ggplot(data = input,aes(x = log2fc,y=-log10(padjust),colour=Significant,shape=Significant))+ ggtitle(title) +
#  geom_point( size=2) + scale_color_manual(values=colo)+
#  scale_shape_manual(values=c(17,16,15)) + xlab(x_axis_title) + ylab(y_axis_title) + theme_bw() +xlim(-the,the)+
#  theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
#  legend.title = element_text(size = 15),legend.text = element_text(size = 12),
#  legend.position = c(1, 0.7),legend.justification = c(1, 1),legend.background = element_rect(fill = NULL, colour = "black"))+
#  geom_vline(xintercept=c(-log2fc_cutoff,log2fc_cutoff),lty=2,col="black",lwd=0.4)+geom_hline(yintercept = -log10(padjust_cutoff),lty=2,col="black",lwd=0.4)
# dev.off()

pdf(file = 'volcano.pdf',height = 8,width = 8)
ggplot(data = input,aes(x = log2fc,y=-log10(padjust),colour=Significant,shape=Significant))+ ggtitle(title) +
 geom_point( size=2) + geom_text_repel(data=subset(input, seq_id %in% show_id), aes(label=seq_id)) + scale_color_manual(values=colo)+
 scale_shape_manual(values=c(17,16,15)) + xlab(x_axis_title) + ylab(y_axis_title) + theme_bw() +xlim(-the,the)+
 theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
 legend.title = element_text(size = 15),legend.text = element_text(size = 12),
 legend.background = element_rect(fill = NULL, colour = "black"))+
 geom_vline(xintercept=c(-log2fc_cutoff,log2fc_cutoff),lty=2,col="black",lwd=0.4)+geom_hline(yintercept = -log10(padjust_cutoff),lty=2,col="black",lwd=0.4)
dev.off()


png(file = 'volcano.png')
ggplot(data = input,aes(x = log2fc,y=-log10(padjust),colour=Significant,shape=Significant))+ ggtitle(title) +
 geom_point( size=2) + geom_text_repel(data=subset(input, seq_id %in% show_id), aes(label=seq_id)) + scale_color_manual(values=colo)+
 scale_shape_manual(values=c(17,16,15)) + xlab(x_axis_title) + ylab(y_axis_title) + theme_bw() +xlim(-the,the)+
 theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
 legend.title = element_text(size = 15),legend.text = element_text(size = 12),
 legend.background = element_rect(fill = NULL, colour = "black"))+
 geom_vline(xintercept=c(-log2fc_cutoff,log2fc_cutoff),lty=2,col="black",lwd=0.4)+geom_hline(yintercept = -log10(padjust_cutoff),lty=2,col="black",lwd=0.4)
dev.off()

svg(file = 'volcano.svg',height = 8,width = 8)
ggplot(data = input,aes(x = log2fc,y=-log10(padjust),colour=Significant,shape=Significant))+ ggtitle(title) +
 geom_point( size=2) + geom_text_repel(data=subset(input, seq_id %in% show_id), aes(label=seq_id)) + scale_color_manual(values=colo)+
 scale_shape_manual(values=c(17,16,15)) + xlab(x_axis_title) + ylab(y_axis_title) + theme_bw() +xlim(-the,the)+
 theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
 legend.title = element_text(size = 15),legend.text = element_text(size = 12),
 legend.background = element_rect(fill = NULL, colour = "black"))+
 geom_vline(xintercept=c(-log2fc_cutoff,log2fc_cutoff),lty=2,col="black",lwd=0.4)+geom_hline(yintercept = -log10(padjust_cutoff),lty=2,col="black",lwd=0.4)
dev.off()


