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
result <- matrix(nrow = nrow(otu_data),ncol = 3)
pvalue <- 1
for(i in 1:nrow(otu_data)){
  c1 <- as.numeric(as.vector(otu_data$"${sample1}")[i])
  c2 <- sum(as.numeric(as.vector(otu_data$"${sample1}"))) - c1
  c3 <- as.numeric(as.vector(otu_data$"${sample2}")[i])
  c4 <- sum(as.numeric(as.vector(otu_data$"${sample2}"))) - c3
  data <- matrix(c(c1,c2,c3,c4),ncol = 2)
  test <- "${choose_test}"
  if (test == "chi") {
    tt <- chisq.test(data)
  }else{
    tt <- fisher.test(data,alternative = "${test_type}",conf.level = ${ci})
  }
  pvalue <- c(pvalue,tt$p.value)
  # pro1 <- c1 / sum(as.numeric(as.vector(otu_data$"${sample1}")))
  # pro2 <- c3 / sum(as.numeric(as.vector(otu_data$"${sample2}")))
  
  result[i,] = c(rownames(otu_data)[i],c1,c3)
}
pvalue <- pvalue[-1]
pvalue <- p.adjust(as.numeric(pvalue),method = "none")
qvalue <- p.adjust(as.numeric(pvalue),method = "${mul_test}")
for(i in 1:length(pvalue)){
  pvalue[i] <- signif(pvalue[i],4)
  qvalue[i] <- signif(qvalue[i],4)
}
result <- cbind(result,pvalue)
result <- cbind(result,qvalue)
colnames(result) <- c(" ",paste(s1,"-propotion",sep=''),paste(s2,"-propotion",sep=''),"pvalue","corrected_pvalue")
result_order <- result[order(-(as.numeric(result[,2])+as.numeric(result[,3]))),]
if(lendata == 1){
  a <- data.frame(result_order)
  result_order <- t(a)
}
write.table(result_order,"${outputfile}",sep="\t",col.names=T,row.names=F,quote = F)
