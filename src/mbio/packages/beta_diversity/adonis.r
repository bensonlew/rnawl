##############loadings file#########################
dist_matrix <- read.table("${dist_matrix}",comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE)
dist_matrix_temp <- read.table("${dist_matrix}",comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE, colClasses = c("character"))
rownames(dist_matrix) <- row.names(dist_matrix_temp)
group_factor <- read.table("${group_file}", comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE)
group_factor_temp <- read.table("${group_file}", comment.char = "",sep = "\t",row.names=1, header = TRUE, check.names=FALSE, colClasses = c("character"))
rownames(group_factor) <- row.names(group_factor_temp)
############calculate adonis##########################
library(vegan)
set.seed(20)
adonis_result <- adonis(as.dist(dist_matrix) ~ group_factor$"${flied}", permutations = ${permutations})
${flied} <- c('', "${flied}")
R2 <- c('R^2', adonis_result$aov.tab$R2[1])
pr <- c('Pr(>F)', adonis_result$aov.tab$'Pr(>F)'[1])
format_adonis <- data.frame(${flied}, R2, pr)
format_adonis_file <- paste("${output_dir}", "/adonis_format.txt", sep='')
write.table(format_adonis, format_adonis_file, sep='\t', , col.names=FALSE, quote=FALSE, row.names=FALSE)
result_file <- paste("${output_dir}", "/adonis_results.txt", sep='')
df <- c('Df', adonis_result$aov.tab$Df)
sum <- c('SumsOfSqs', adonis_result$aov.tab$SumsOfSqs)
mean <- c('MeanSqs', adonis_result$aov.tab$MeanSqs)
mean[is.na(mean)] <- '-'
model <- c('F_Model', adonis_result$aov.tab$F.Model)
model[is.na(model)] <- '-'
r <- c('R^2', adonis_result$aov.tab$R2)
r[is.na(r)] <- '-'
p <- c('Pr(>F)', adonis_result$aov.tab$'Pr(>F)')
p[is.na(p)] <- '-'
name <- c('#', "group_factor$${flied}", "Residuals", "Total")
adonis_result_file <- data.frame(name, df, sum, mean, model, r, p)
write.table(adonis_result_file, result_file, sep='\t', , col.names=FALSE, quote=FALSE, row.names=FALSE)
#sink(result_file)
#print(adonis_result)
#sink(NULL)
