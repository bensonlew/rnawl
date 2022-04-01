#load packages
library(getopt)
library(ggbio)
library(rtracklayer)

#set options
command <- matrix(c(
  'seq', 's', 1, 'character', 'consists of gene id/gene name/transcript id. If multiple items supplied, separating by comma',
  'annotation','a',1, 'character','user input annotation file of the species, format should be in BED/GTF/GFF', 
  'help', 'h', 0, 'logical', 'show this help message and exit'
), byrow = TRUE, ncol = 5)

opts <- getopt(command)

#check options
if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}

if (is.null(opts$seq) || is.null(opts$annotation)) {
  cat(getopt(command, usage = TRUE))
  q(status=-2)
}

#main
seq_list<-strsplit(opts$seq,",",fixed = TRUE)[[1]]
anno<-import(opts$annotation)
anno_filtered<-anno[tolower(anno$type) %in% c("cds", "exon", "utr", "gap","gene"),]

plot_list<-c()
file_list<-c()

for (i in 1:length(seq_list)){
    if (seq_list[i] %in% anno_filtered$gene_name){
        selected<-anno_filtered[anno_filtered$gene_name %in% c(seq_list[i])]
        title<-paste0("Gene Name: ",seq_list[i])
    }else if (seq_list[i] %in% anno_filtered$gene_id) {
        selected<-anno_filtered[anno_filtered$gene_id %in% c(seq_list[i])]
        title<-paste0("Gene ID: ",seq_list[i])
    }else if (seq_list[i] %in% anno_filtered$transcript_id){
        selected<-anno_filtered[anno_filtered$transcript_id %in% c(seq_list[i])]
        title<-paste0("Transcript ID: ",seq_list[i])
    }else{
        print(paste0(seq_list[i]," cannot be found in your annotation file."))
        next
    }

    if (!('cds' %in% tolower(selected$type)) & ('exon' %in% tolower(selected$type)) &
     ('utr' %in% tolower(selected$type))){
        colour<-c("black","white")
    }else{
        colour<-c("black", NA, NA)
    }

    names(selected)<-selected@elementMetadata@listData[["transcript_id"]]
    selected<-split(selected,names(selected))

    #plotting
    p<-ggbio()+geom_alignment(selected,color="black",gap.geom = "chevron",label.color = "black",
                              label.size = 5,aes(fill=type))+scale_fill_manual(values=colour)+
                              theme_bw()+theme(legend.position="none")+labs(title=title)

    plot_list=c(plot_list,p)
    file_list=c(file_list,paste0(seq_list[i],".pdf"))
}

# save as pdf
for (p in 1:length(plot_list)){
  pdf(file_list[p])
  print(plot_list[[p]])
  dev.off()
}






