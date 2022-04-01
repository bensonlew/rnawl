otu_table <- read.table("${otu_file}",comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE)
otu_table_temp <- read.table("${otu_file}",comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE, colClasses = c("character"))
rownames(otu_table) <- row.names(otu_table_temp)
otu_table <- t(otu_table)
# fix(otu_table)
env_factor <- read.table("${env_file}", comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE)
env_factor_temp <- read.table("${env_file}", comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE, colClasses = c("character"))
rownames(env_factor) <- row.names(env_factor_temp)
# fix(env_factor)
centroids <- FALSE
biplot <- FALSE
for (i in env_factor) {
    if (is.factor(i)){
        centroids <- TRUE
    } else{
        biplot <- TRUE
    }
}
library(vegan)
capscale.result<- capscale(otu_table~.,as.data.frame(env_factor),dist="${distance_algorithm}",sqrt.dist=${sqrt_dist},add=${add})
cot<-summary(capscale.result)
pdf('${output_dir}/db_rda.pdf')
plot_values <- plot(capscale.result)
dev.off()
write.table(plot_values$sites, '${output_dir}/db_rda_sites.xls', sep = '\t', col.names = NA, quote = FALSE)
#write.table(plot_values$species, '${output_dir}/db_rda_species.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(cot$cont$importance, '${output_dir}/db_rda_cont.xls', sep = '\t', col.names = NA, quote = FALSE)##add by zhujuan for Proportion Explained 2017.08.21
if (centroids){
    write.table(plot_values$centroids, '${output_dir}/db_rda_centroids.xls', sep = '\t', col.names = NA, quote = FALSE)
}
if (biplot){
    write.table(plot_values$biplot * attr(plot_values$biplot,"arrow.mul"), '${output_dir}/db_rda_biplot.xls', sep = '\t', col.names = NA, quote = FALSE)
}
sink('${output_dir}/env_data.temp')
if (centroids) {
    print('centroids:TRUE')
}else{
    print('centroids:FALSE')
}
if (biplot) {
    print('biplot:TRUE')
}else{
    print('biplot:FALSE')
}
sink(NULL)
