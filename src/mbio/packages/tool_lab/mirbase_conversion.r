library(miRBaseConverter)
library(methods)
library(stringr)
library(getopt)

command <- matrix(c(
  'infile', 'f', 2, 'character', 'User-input file with original miRNA ID/Accessions in the first column.',
  'instr', 's', 2, 'character', 'A string of miRNA names/Accessions, separated by comma.',
  'target_version', 't', 1, 'character', 'One of the available version: v6,v7_1,v8,v8_1,v8_2,v9,v9_1,v9_2,v10,v10_1,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20,v21,v22.',
  'help', 'h', 0, 'logical', 'Show this help message and exit'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

#check options
if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}

if (is.null(opts$infile) & is.null(opts$instr)) {
  cat(getopt(command, usage = TRUE))
  q(status=-2)
}

# process inputs
t_version<-opts$target
if (!is.null(opts$infile)){
    origin_df<-read.table(opts$infile,header = F)
    origin<-origin_df[,1]
}else if (!is.null(opts$instr)){
    origin<-strsplit(opts$instr,',')
}


# conversion
if (grepl('-', origin[1])==TRUE){
    convert<-miRNAVersionConvert(origin,targetVersion = t_version,exact = TRUE,verbose = FALSE)
    o_version<-checkMiRNAVersion(origin, verbose = F)
    col_names<-c(paste0(o_version, '_OriginalName'), paste0(t_version,"_TargetName"), paste0(t_version,"_Accession"))
}else{
    convert<-miRNA_AccessionToName(origin, targetVersion = t_version)
    col_names<-c('OriginalName', paste0(t_version,"_TargetName"))
}

uniq_name<-str_split_fixed(convert$TargetName,"\\&",2)
convert[,2]<-uniq_name[,1]

## replace column1 in input file using converted miRNA Names
#if (!is.null(opts$infile)){
#    origin_df[,1]<-sapply(origin_df[,1], function(x) x=ifelse(convert[convert[[1]]==x,2]=='', as.character(x), convert[convert[[1]]==x,2]))
#    write.table(origin_df,paste0(t_version, "_ID_infile.xls"),sep = "\t",quote = F,row.names = F,col.names = F)
#}

# generate result file
colnames(convert)<-col_names
write.table(convert,'result.xls',sep = "\t",quote = F, row.names = F, na = 'NA')