# load package

library(getopt)

# set options
command <- matrix(c(
    "json", "i", 1, "character", "setting file in JSON format",
    "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (! is.null(opts$help)) {
    cat(getopt(command, usage = TRUE))
    q(status = 2)
}
if (is.null(opts$json)) {
    cat(getopt(command, usage = TRUE))
    q(status = - 1)
}

library(rjson)
library(org.Mm.eg.db)
library(dplyr)
library(biomaRt)
#library(clusterProfiler)
json.lst <- fromJSON(file = opts$json)
ensembl_id <-json.lst$ensembl_id
output.path <-json.lst$output
mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "asia",host='asia.ensembl.org')
#mart <- useMart("ensembl","mmusculus_gene_ensembl")

gene_id_name_description<-getBM(attributes=c("ensembl_gene_id","external_gene_name","description"),filters = "ensembl_gene_id",values = ensembl_id, mart = mart)
names(gene_id_name_description)=c('gene_id','gene_name', 'description')
write.table(new_dropNA,file=output.path,sep = '\t', row.names = F,col.names = TRUE)


#symbol<-mapIds(org.Mm.eg.db,keys = ensembl_id,column ="SYMBOL",keytype = "ENSEMBL", multiVals = 'first')
#symbol <- as.matrix(symbol)
#symbol= as.data.frame(symbol)
#colnames(symbol) = 'SYMBOL'
#ENSEMBL = row.names(symbol)
#new = data.frame(ENSEMBL,symbol)
#new['description']='NA'
#names(new)=c('gene_id','gene_name', 'description')
#new_dropNA = filter(new, gene_name != 'NA')



