#!/usr/bin/R
# -*- coding: utf-8 -*-
# __author__ = konghualei, 20170418

fpkm = "/mnt/ilustre/users/sanger-dev/workspace/20170413/Single_rsem_stringtie_zebra_9/Express/output/oldrsem/genes.TMM.fpkm.matrix"
out_path = "/mnt/ilustre/users/sanger-dev/sg-users/konghualei"
library(gplots)
data = read.delim(fpkm, sep="\t",header=T,as.is=T,row.names=1)
outfile_path = out_path
data_venn = list()
for(s in seq(dim(data)[2])){
    data_venn[[s]] = data[,s]
}
venn_data=venn(data_venn, show.plot=F)
venn_info = attr(venn_data,"intersections")
venn_names = names(venn_info)
venn_table = data.frame(ncol = 3, nrow = length(venn_names))
for(i in venn_names){
    tmp = venn_info[venn_names]
    tmp_name = paste(data,sep=",")
    venn_table[i]=c(i,length(tmp),tmp_name)
}
print(venn_table)
write(venn_table,paste(outfile_path,"venn_table.xls",sep="/"),sep="\t",quote=F,rownames=NULL, colnames=NULL)
# venn_data = t(data)
# venn_graph = data.frame(ncol=2)
# for(j in seq(dim(data)[2])){
    # tmp1 = paste(data[,j],sep=",")
    # venn_graph = rbind(venn_graph, c(colnames(data)[j],tmp1))
# }
# colnames(venn_graph) = c("#group_name","gene_id")
# write(venn_graph, paste(outfile_path,"venn_graph.xls",sep="/"), sep="\t",quote=F,colnames=c("#group_name","gene_id"))

