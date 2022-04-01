# load package
library(ggplot2)

args = commandArgs(T)
print(args)
busco_path = args[1]
title = args[2]
x_lab = args[3]
output_file = args[4]
title_all = paste('BUSCO Assessment Results\n', title)
busco = read.table(busco_path, sep='\t', header=T)
busco$name = factor(busco$name, levels = busco$name)
p = ggplot(busco,aes(spe,persent,fill=name))+
geom_bar(stat="identity",position="stack",width = 0.4)+
theme_classic()+
coord_flip()+
xlab(x_lab) +
ylab('%Busco')+
theme(text=element_text(size=16,  family="serif"))+
theme(axis.text.y = element_blank())+
ggtitle(title_all)+
theme(legend.title = element_blank())
ggsave(output_file, p, height = 3.5)


