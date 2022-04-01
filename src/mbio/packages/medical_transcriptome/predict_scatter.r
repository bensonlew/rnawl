# Title     : sub-script of table.kit.py
# Objective : use for scatter_plot for somatic_predict results
# Created by: fwy
# Created on: 2020/11/11

# load package
library(getopt)
library(ggplot2)
# set options
command <- matrix(c(
  "predict_matrix", "i", 1, "character", "predict_matrix after predeal",
  "output", "o", 1, "character", "output directory for results",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

# check options
if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status = -1)
}

# read data and predeal
raw_data <- as.matrix(read.table(opts$predict_matrix,sep ="\t",header = TRUE))
data_df  <- as.data.frame(raw_data)
data_df$Polyphen2_HDIV_score <-  as.numeric(levels(data_df$Polyphen2_HDIV_score)[data_df$Polyphen2_HDIV_score])
data_df$SIFT_score<-  as.numeric(levels(data_df$SIFT_score)[data_df$SIFT_score])
pdf(paste(opts$output,"pdf",sep = "."))
ggplot(data_df,aes(x=Polyphen2_HDIV_score,y= SIFT_score,colour = type))+geom_point()+
    scale_x_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1))+scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1))+
    geom_hline(aes(yintercept=0.05), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=0.453), colour="#990000", linetype="dashed")
dev.off()
png(paste(opts$output,"png",sep = "."))
ggplot(data_df,aes(x=Polyphen2_HDIV_score,y= SIFT_score,colour = type))+geom_point()+
    scale_x_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1))+scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1))+
    geom_hline(aes(yintercept=0.05), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=0.453), colour="#990000", linetype="dashed")
dev.off()
svg(paste(opts$output,"svg",sep = "."))
ggplot(data_df,aes(x=Polyphen2_HDIV_score,y= SIFT_score,colour = type))+geom_point()+
    scale_x_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1))+scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1))+
    geom_hline(aes(yintercept=0.05), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=0.453), colour="#990000", linetype="dashed")
dev.off()


