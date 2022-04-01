#!/usr/bin/env Rscript

# load package
library(getopt)
library(MEGENA)
library(ggplot2)


# set options
command <- matrix(c(
  "exp_table", "e", 1, "character", "input file express matrix",
  "tmp_dir", "t", 1, "character", "tmp work dir ",
  "config", "c", 1, "character", "plot config info ",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)


print("start set options")
config = read.table(opts$config, sep='\t', header=TRUE)

n.cores <- 2; # number of cores/threads to call for PCP
doPar <-TRUE; # do we want to parallelize?
method = config$method  # method for correlation. either pearson or spearman.
FDR.cutoff = config$fdr_cutoff  # method for correlation. either pearson or spearman.
module.pval = config$module_pval # module significance p-value. Recommended is 0.05.
min_size = config$min_size #min size for cluster
hub.pval = 0.05  # connectivity significance p-value based random tetrahedral networks
cor.perm = 10  # number of permutations for calculating FDRs for all correlation pairs.
hub.perm = 100  # number of permutations for calculating connectivity significance p-value.
tmp_dir = opts$tmp_dir
annot.table=NULL
id.col = 1
symbol.col= 2



print("start read_table")
# prepare data
datExpr <- read.table(opts$exp_table, sep='\t', header=TRUE)
# ijw <- calculate.correlation(datExpr,doPerm = cor.perm,method=method,output.corTable = FALSE,output.permFDR = FALSE)
print("end read_table")
print("start calculate correlation")
ijw <- calculate.correlation(datExpr,doPerm = cor.perm,method=method,output.corTable = TRUE,output.permFDR = TRUE,saveto = tmp_dir)
print("end calculate correlation")
print("start calculate PFN")
el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores,keep.track = FALSE)
g <- graph.data.frame(el,directed = FALSE)
print("end calculate PFN")
MEGENA.output <- do.MEGENA(g,
mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
min.size = min_size,max.size = vcount(g)/2,
doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
save.output = TRUE)

summary.output <- MEGENA.ModuleSummary(MEGENA.output,
mod.pvalue = module.pval,hub.pvalue = hub.pval,
min.size = min_size,max.size = vcount(g)/2,
annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
output.sig = TRUE)


if (!is.null(annot.table))
{
# update annotation to map to gene symbols
V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name,annot.table[[id.col]])],V(g)$name,sep = "|")
summary.output <- output[c("mapped.modules","module.table")]
names(summary.output)[1] <- "modules"
}

print(head(summary.output$modules,2))
print(summary.output$module.table)
dir.create(paste(tmp_dir, "/", "modue_network_plots", sep = ""))
result_dir <- paste(tmp_dir, "/", "modue_network_plots", sep = "")
for (module in summary.output$module.table$module.id){
    pdf(paste(result_dir,paste(module,"pdf", sep = "."), sep = "/"))
    a<-plot_module(output.summary = summary.output,PFN = g,subset.module = module,
                 layout = "kamada.kawai",label.hubs.only = TRUE,
                 gene.set = NULL,color.code =  "grey",
                 output.plot = FALSE,out.dir = "modulePlot",col.names = c("magenta","green","cyan"),label.scaleFactor = 20,
                 hubLabel.col = "black",hubLabel.sizeProp = 1,show.topn.hubs = Inf,show.legend = TRUE)
    print(a[1])
    dev.off()
}

warnings()

module.table <- summary.output$module.table
colnames(module.table)[1] <- "id" # first column of module table must be labelled as "id".
pdf(paste(result_dir,paste("module_hierarchy","pdf", sep = "."), sep = "/"))
plot_module_hierarchy(module.table = module.table,label.scaleFactor = 0.15,
                                    arrow.size = 0.03,node.label.color = "blue")
dev.off()



output.summary = summary.output
mdf= output.summary$module.table
pdf(paste(result_dir,paste("module_hill","pdf", sep = "."), sep = "/"))
draw_sunburst_wt_fill(module.df = mdf,feat.col = "module.size",log.transform = FALSE,
    fill.type = "continuous",
    fill.scale = scale_fill_gradient2(low = "white",high = "red",
    na.value = "white"),
    id.col = "module.id",parent.col = "module.parent")
dev.off()


