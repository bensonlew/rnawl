#!/usr/bin/env Rscript

# load package
library(getopt)
library(ggplot2)


# set options
command <- matrix(c(
  "data", "i", 1, "character", "input file containing processed data",
  "fc_cutoff", "l", 1, "double", "log2fc_cutoff for classify ",
  "padjust", "p", 1, "double", "padjust_cutoff for classify ",
  "method", "m", 1, "integer", "color schemes",
  "x_axis_name", "x", 1, "character", "x_axis_name",
  "y_axis_name", "y", 1, "character", "y_axis_name",
  "title_name", "t", 1, "character", "title",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

print("start")
print(opts$y_axis_name)
# prepare data
input = read.table(opts$data, sep='\t', header=TRUE)


# calculate
log2fc_cutoff = log2(opts$fc_cutoff)
padjust_cutoff = opts$padjust
total_num = nrow(input)
up_data<- input[input$padjust < padjust_cutoff & input$log2fc > log2fc_cutoff, ]
up_num = nrow(up_data)
down_data<- input[input$padjust < padjust_cutoff & input$log2fc > log2fc_cutoff, ]
down_num = nrow(down_data)
up_name = paste("up:",as.character(up_num))
down_name = paste("down:",as.character(down_num))
no_name = paste("no_diff:",as.character(total_num - up_num - down_num ))
colo = c("blue", "grey","red")
if (opts$method == 1 ){
   colo = c("blue","grey","red")
   }else if (opts$method == 2){
   colo = c("green", "grey","red")
   }else if (opts$method == 3){
   colo = c("blue", "black","red")
   }else if(opts$method == 4){
   colo = c("green", "black","red")
   }
the <- max(abs(range(input$log2fc)))
the <- ceiling(the)

input$Significant = as.factor(ifelse (input$padjust < padjust_cutoff & abs(input$log2fc)>= log2fc_cutoff,ifelse(input$log2fc >log2fc_cutoff,up_name,down_name),no_name))
pdf(file = 'volcano.pdf',height = 5,width = 5)
ggplot(data = input,aes(x = log2fc,y=-log10(padjust),colour=Significant,shape=Significant))+ ggtitle(opts$title_name) +
 geom_point( size=2) + scale_color_manual(values=colo)+
 scale_shape_manual(values=c(17,16,15)) + xlab(opts$x_axis_name) + ylab(opts$y_axis_name) + theme_bw() +xlim(-the,the)+
 theme(plot.title = element_text(hjust = 0.5) ,axis.title = element_text(size = 15),axis.text = element_text(size = 12),
 legend.title = element_text(size = 15),legend.text = element_text(size = 12),
 legend.position = c(1, 0.7),legend.justification = c(1, 1),legend.background = element_rect(fill = NULL, colour = "black"))+
 geom_vline(xintercept=c(-log2fc_cutoff,log2fc_cutoff),lty=2,col="black",lwd=0.4)+geom_hline(yintercept = -log10(padjust_cutoff),lty=2,col="black",lwd=0.4)
dev.off()


