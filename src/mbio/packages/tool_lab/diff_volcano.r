# load package
library(ggplot2)

args = commandArgs(T)
print(args)

if (length(args) > 2){
    p_type <- args[3]
}else{
    p_type <- 'pvalue'
}

for (i in Sys.glob(paste(args[1], "/*vs*",sep=''))){
    print(i)
    volcano_path = i
    volcano = read.table(volcano_path, sep='\t', header=TRUE)
    pdf_name = paste("/", basename(i), "_volcano", ".pdf", sep="")
    print(pdf_name)
    pdf_path = paste(args[1], pdf_name, sep="")
    print(pdf_path)
    if (args[2] == 'noiseq'){
    p = ggplot(volcano, aes(x = log2fc, y = D, colour=regulate)) + geom_point() +
    theme_bw() +
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
    } else{
    p = ggplot(volcano, aes(x = log2fc, y = log10pvalue, colour=regulate)) + geom_point() +
    theme_bw() + labs(y=paste0('-log10(', p_type, ')'), x='log2(fc)') +
    scale_colour_manual(name="regulate", values=c('forestgreen', 'grey41', 'red'))+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
    }
    ggsave(pdf_path, p)
}






