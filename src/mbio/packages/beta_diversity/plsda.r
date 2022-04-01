##############loadings file#########################
otu_table <- read.table("${otu_file}",comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE)
otu_table_temp <- read.table("${otu_file}",comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE, colClasses = c("character"))
rownames(otu_table) <- row.names(otu_table_temp)
otu_table <- t(otu_table)
env_factor <- read.table("${env_file}", comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE)
env_factor_temp <- read.table("${env_file}", comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE, colClasses = c("character"))
rownames(env_factor) <- row.names(env_factor_temp)
##############format table##########################
inter_samples <- sort(intersect(rownames(env_factor), rownames(otu_table)))
env_temp <- env_factor[inter_samples, ]
if (is.data.frame(env_temp)){
    env_factor <- env_temp
}else{
    env_temp <- data.frame(env_temp)
    rownames(env_temp) <- inter_samples
    colnames(env_temp) <- colnames(env_factor)
    env_factor <- env_temp
}

otu_temp <- otu_table[inter_samples, ]
if (is.data.frame(otu_temp)){
    otu_table <- otu_temp
}else{
    otu_temp <- data.frame(otu_temp)
    rownames(otu_temp) <- inter_samples
    colnames(otu_temp) <- colnames(otu_table)
    otu_table <- otu_temp
}
############calculate plsda##########################
library(mixOmics)
num <- length(unique(env_factor$${group_name}))
normal_run <- function(){
    print('start run plsda with maximal and near zero.');
    plsda_otu <<- plsda(otu_table,env_factor$${group_name},ncomp=num,near.zero.var=T)
}
error_run <- function(e=''){
    print('ERROR:');
    print(e);
    print('restart run plsda without maximal and near zero.');
    plsda_otu <<- plsda(otu_table,env_factor$${group_name})
}
tryCatch(normal_run(), error=function(e){error_run(e)})
plsda_sites<-plsda_otu$variates$X
plsda_rotat<-plsda_otu$loadings$X
plsda_impo<-plsda_otu$loadings$Y
plsda_pre<-plsda_otu$explained_variance$X
plsda_sites_file<-paste("${output_dir}", "/plsda_sites.xls", sep='')
plsda_rotat_file<-paste("${output_dir}", "/plsda_rotation.xls", sep='')
plsda_impo_file<-paste("${output_dir}", "/plsda_importance.xls", sep='')
plsda_pre_file<-paste("${output_dir}","/plsda_pre.xls", sep='')
write.table(plsda_sites, plsda_sites_file, sep="\t", col.names = NA, quote = FALSE)
write.table(plsda_rotat, plsda_rotat_file, sep="\t", col.names = NA, quote = FALSE)
write.table(plsda_impo, plsda_impo_file, sep="\t", col.names = NA, quote = FALSE)
write.table(plsda_pre, plsda_pre_file, sep="\t", col.names = NA, quote = FALSE)
