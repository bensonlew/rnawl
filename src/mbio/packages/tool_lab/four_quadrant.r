#!/usr/bin/env Rscript

library(ggpubr)
library(getopt)

command <- matrix(c(
  'dataframe', 'd', 1, 'character', 'user input file. A 5-column-dataframe for graph.',
  'x_lab', 'x', 2, 'character', 'x-axis label',
  'y_lab', 'y', 2, 'character', 'y-axis label',
  'help', 'h', 0, 'logical', 'show this help message and exit'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}
if (is.null(opts$dataframe)) {
  cat(getopt(command, usage = TRUE))
  q(status=-2)
}

## 设置默认值
if ( is.null(opts$x_lab) ) {
  opts$x_lab = 'Protein(Log2FC)'
}
if ( is.null(opts$y_lab) ) {
  opts$y_lab = "Gene(Log2FC)"
}


data<-read.table(opts$dataframe, header = TRUE)
data$Significance<-paste(data$sig2, data$sig1,sep='_')
lm_formula<-as.formula(paste(colnames(data)[2], colnames(data)[4], sep='~'))
result<-summary(lm(lm_formula, data = data))

p <- ggscatter(data,
               x = colnames(data)[4],
               y = colnames(data)[2],
               xlab = opts$x_lab,
               ylab = opts$y_lab,
               color = "Significance",
               size = 1.5,
               font.label = 7,
               repel = T,
               legend='right',
               title = paste(strsplit(opts$y_lab, '(', fixed = TRUE)[[1]][1], 'vs.', strsplit(opts$x_lab, '(', fixed = TRUE)[[1]][1], sep=' '),
               palette = c("#CC0000", "#2F5688", "#7594D0", '#BEBEBE'))+
  stat_cor(method = "pearson")+
  theme(plot.title = element_text( size = 12.5, hjust = 0.5))+
  geom_hline(yintercept=0,linetype=4)+
  geom_vline(xintercept=0,linetype=4) +
  geom_abline(slope=as.numeric(result$coefficients[2]),
              intercept=as.numeric(result$coefficients[1]))

pdf("four_quadrant.pdf", onefile=FALSE)
plot(p)
dev.off()
