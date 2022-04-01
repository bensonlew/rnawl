otu_table <- read.table("${species_table}",comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE)  # add 'species_table' by zhujuan for add db_rda_species.xls 20171010
otu_table_temp <- read.table("${species_table}",comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE, colClasses = c("character"))
rownames(otu_table) <- row.names(otu_table_temp)
otu_table <- t(otu_table)
dis_matrix <- read.table("${dis_matrix}",comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE)
dis_matrix_temp <- read.table("${dis_matrix}",comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE, colClasses = c("character"))
rownames(dis_matrix) <- row.names(dis_matrix_temp)
env_factor <- read.table("${env_file}", comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE)
env_factor_temp <- read.table("${env_file}", comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE, colClasses = c("character"))
rownames(env_factor) <- row.names(env_factor_temp)
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
capscale.result<- capscale(as.dist(dis_matrix)~.,as.data.frame(env_factor),comm=otu_table,sqrt.dist=${sqrt_dist},add=${add})
fit <- try(cot<-summary(capscale.result))
if ("try-error" %in% summary(fit)){
capscale.result<- dbrda(as.dist(dis_matrix)~.,as.data.frame(env_factor),comm=otu_table,sqrt.dist=${sqrt_dist},add=${add})
}
ef<-envfit(capscale.result,as.data.frame(env_factor),permu=999)
ef<-data.frame(ef$vectors$arrows,ef$vectors$r,ef$vectors$pvals)
colnames(ef)[length(colnames(ef))] = "p_values"
colnames(ef)[length(colnames(ef))-1] = "r2"
pdf('${output_dir}/db_rda.pdf')
plot_values <- plot(capscale.result)
dev.off()
write.table(plot_values$sites, '${output_dir}/db_rda_sites.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(ef,'${output_dir}/db_rda_envfit.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(summary(capscale.result)$cont$importance, '${output_dir}/db_rda_cont.xls', sep = '\t', col.names = NA, quote = FALSE)## add by zhujuan for Proportion Explained  2017.08.21
 write.table(plot_values$species, '${output_dir}/db_rda_species.xls', sep = '\t', col.names = NA, quote = FALSE)
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
