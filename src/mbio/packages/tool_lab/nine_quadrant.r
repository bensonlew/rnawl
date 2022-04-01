# load package
library(ggplot2)
library(RColorBrewer)
library(scales)

args = commandArgs(T)
print(args)
nine_quadrant = read.table(args[1], sep='\t', header=T)
y = max(abs(nine_quadrant$trans_log2fc))
x = max(abs(nine_quadrant$protein_log2fc))

pdf_path = args[4]
mycolors <- brewer.pal(9, args[5])
p = ggplot(data=nine_quadrant,aes(x=protein_log2fc, y=trans_log2fc, colour=group)) + geom_point() +
geom_vline(xintercept=c(-log2(as.numeric(args[2])),log2(as.numeric(args[2]))), linetype="dotted") +
geom_hline(yintercept=c(-log2(as.numeric(args[3])),log2(as.numeric(args[3]))), linetype="dotted") + ylim(-y,y) + xlim(-x,x) +
theme_bw() +
theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) +
scale_color_manual(values=mycolors) +
xlab(args[6]) +
ylab(args[7])
ggsave(pdf_path, p)



