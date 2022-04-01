# houshuang 20190926

library(qiimer)
library(biom)
library(Matrix)
library(Tax4Fun)
library(getopt)
	
command <- matrix(c("infile", "i", 1, "character",
					"outfile", "o", 1, "character",
					"database", "d", 1, "character",
					"kegg", 'k', 1, "character"), byrow = T, ncol = 4)

args <- getopt(command)

infile <- args$infile
outfile <- args$outfile
database <- args$database
kegg <- args$kegg

# Tax4Fun functional profiling
QIIMESingleData <- importQIIMEData(infile)
Tax4FunOutput <- Tax4Fun(QIIMESingleData, database, fctProfiling = TRUE, refProfile = "UProC", shortReadMode = TRUE, normCopyNo = TRUE)
tax4fun_gene <- as.data.frame(t(Tax4FunOutput$Tax4FunProfile))
temp <- tax4fun_gene$gene
temp$KO <- t(data.frame(strsplit(rownames(tax4fun_gene), ';')))[ ,1]
temp$Description <- t(data.frame(strsplit(rownames(tax4fun_gene), '; ')))[ ,2]
write.table(cbind(temp,tax4fun_gene), outfile, row.names = FALSE, sep = '\t', quote = FALSE)

# for check result
otu_table = QIIMESingleData$otuTable
colSums(otu_table)
write.table("ID\t", paste(outfile, '.data', sep=""), append = FALSE, quote = FALSE, sep="\t",eol = "", na = "NA", dec = ".", row.names = F, col.names = F)
write.table(otu_table, paste(outfile, '.data', sep=""), append = T, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

# Tax4Fun metabolic profiling
Tax4FunOutput <- Tax4Fun(QIIMESingleData, database, fctProfiling = FALSE, refProfile = 'UProC', shortReadMode = TRUE, normCopyNo = TRUE)
tax4fun_pathway <- as.data.frame(t(Tax4FunOutput$Tax4FunProfile))
pathway <- rownames(tax4fun_pathway)
write.table(cbind(pathway, tax4fun_pathway), paste(outfile, '.pathway', sep=""), row.names = FALSE, sep = '\t', quote = FALSE)

# KEGG annotation
kegg_anno <- read.delim(kegg, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
tax4fun_pathway$Pathway_level3 <- t(data.frame(strsplit(rownames(tax4fun_pathway), ';')))[ ,1]
tax4fun_pathway <- merge(kegg_anno, tax4fun_pathway, by = 'Pathway_level3')
write.table(tax4fun_pathway, paste(outfile, '.pathway.anno', sep=""), row.names = FALSE, sep = '\t', quote = FALSE)
