library(parallel)
system.time({
otu_data <- read.table("${inputfile}",sep = "\t",comment.char = '', colClasses="character")
samp <- t(otu_data[1,-1])
otu_data <- otu_data[-1,]
rownames(otu_data) <- otu_data[,1]
otu_data <- otu_data[,-1]
colnames(otu_data) <- samp
s1 <- "${sample1}"
s2 <- "${sample2}"
otu_data <- otu_data[,which(samp %in% c(s1,s2))]
# otu_data <- otu_data[apply(otu_data,1,function(x)any(x>0)),]
otu_data <- otu_data[apply(otu_data,1,function(x)any(x != "0.0")),]
otu_data <- otu_data[apply(otu_data,1,function(x)any(x != "0")),]
lendata <- nrow(otu_data)
test <- "${choose_test}"
sam1_sum <- sum(as.numeric(as.vector(otu_data$"${sample1}")))
sam2_sum <- sum(as.numeric(as.vector(otu_data$"${sample2}")))
})

two_sample_test <- function(i){
  c1 <- as.numeric(as.vector(otu_data$"${sample1}")[i])
  c2 <- sam1_sum - c1
  c3 <- as.numeric(as.vector(otu_data$"${sample2}")[i])
  c4 <- sam2_sum - c3
  data <- matrix(c(c1,c2,c3,c4),ncol = 2)
  if (test == "chi") {
    tt <- chisq.test(data)
    statlist <- tt$statistic
  }else{
    tt <- fisher.test(data,alternative = "${test_type}",conf.level = ${ci})
    statlist <- tt$estimate
  }
  statistic <- c(statlist[[1]])
  pvalue <- c(tt$p.value)
  pro1 <- c1 / sam1_sum
  pro2 <- c3 / sam2_sum
  name <- c(rownames(otu_data)[i])
  samp1 <- c(c1)
  samp2 <- c(c3)
  return(data.frame(name, samp1, samp2, statistic, pvalue))
}
print("cluster start")
system.time({
cl <- makeCluster(2)
clusterExport(cl, c("otu_data", "test", "sam1_sum", "sam2_sum"), envir=environment())
result_list <- parLapply(cl, 1:lendata, two_sample_test)
result <- do.call('rbind', result_list)
stopCluster(cl)
})
print("cluster end")
system.time({
pvalue <- p.adjust(as.numeric(result$pvalue),method = "none")
qvalue <- p.adjust(as.numeric(result$pvalue),method = "${mul_test}")
for(i in 1:length(pvalue)){
  pvalue[i] <- signif(pvalue[i],4)
  qvalue[i] <- signif(qvalue[i],4)
}
result$pvalue <- pvalue
result$qvalue <- qvalue
if(test == "chi"){
  colnames(result) <- c(" ",paste(s1,"",sep=''),paste(s2,"",sep=''),"statistic","pvalue","corrected_pvalue")
result_order <- result[order(-(as.numeric(result[,2])+as.numeric(result[,3]))),]
}else{
  colnames(result) <- c(" ",paste(s1,"",sep=''),paste(s2,"",sep=''),"odds_ratio","pvalue","corrected_pvalue")
result_order <- result[order(-(as.numeric(result[,2])+as.numeric(result[,3]))),]
}
if(lendata == 1){
  a <- data.frame(result_order)
  result_order <- t(a)
}
write.table(result_order,"${outputfile}",sep="\t",col.names=T,row.names=F,quote = F)
})