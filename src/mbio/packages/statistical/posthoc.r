# Rscript post_hoc.r in_table group_file meth out_prefix
library("PMCMRplus")

args <- commandArgs(trailingOnly = TRUE)

in_table <- read.csv(args[1], row.names=1, header=1, sep='\t')
group_table <- read.csv(args[2], row.names=1, header=1, sep='\t', colClasses=c("character"))

gp_lab = NULL
for (g in colnames(in_table)) {
    gp_lab <- c(gp_lab, group_table[g, 1])
}

if (args[3] == "dunn") {
    my_t = kwAllPairsDunnTest
} else if (args[3] == "nemenyi") {
    my_t = kwAllPairsNemenyiTest
} else if (args[3] == "conover-iman") {
    my_t = kwAllPairsDunnTest
} else if (args[3] == "steel_dwass") {
    my_t = dscfAllPairsTest
}

gp_pair <- NULL
result_pvalue <- list()
result_statistic <- list()
for (i in rownames(in_table)){
    one_line <-  as.vector(t(in_table[i,]))
    ans <- my_t(one_line, as.factor(gp_lab))
    pvalue <- ans$p.value
    has_value <- !upper.tri(pvalue)
    colnames(has_value) <- colnames(pvalue)
    rownames(has_value) <- rownames(pvalue)
    statistic <- ans$statistic
    tmp_gpair <- NULL
    one_p <- NULL
    one_s <- NULL
    for (c in colnames(pvalue)) {
        for (r in rownames(pvalue)) {
            if (has_value[r, c]) {
                tmp_gpair <- c(tmp_gpair, paste(c, r, sep='-'))
                one_p <- c(one_p, pvalue[r, c])
                one_s <- c(one_s, statistic[r, c])
            }
        }
    }
    gp_pair <- tmp_gpair
    result_pvalue[[i]] = one_p
    result_statistic[[i]] = one_s
}
pvalue_head <- NULL
statistic_head <- NULL
for (i in 1:length(gp_pair)) {
    pvalue_head <- c(pvalue_head, paste(gp_pair[i], "pvalue", sep='_'))
    statistic_head <- c(statistic_head, paste(gp_pair[i], "statistic", sep='_'))
}
result_pvalue <- t(as.data.frame(result_pvalue))
colnames(result_pvalue) <- pvalue_head
result_statistic <- t(as.data.frame(result_statistic))
colnames(result_statistic) <- statistic_head
result <- cbind(rownames(result_statistic), result_statistic, result_pvalue)

write.table(result, args[4], sep='\t', row.names=F, col.names=T, quote=F)
