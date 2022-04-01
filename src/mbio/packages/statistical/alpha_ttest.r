otu_data <- read.table("${inputfile}",sep = "\t",comment.char = '', colClasses="character")
samp <- t(otu_data[1,-1])
otu_data <- otu_data[-1,]
rownames(otu_data) <- otu_data[,1]
otu_data <- otu_data[,-1]
colnames(otu_data) <- samp
#lendata <- nrow(otu_data)
#read groupfile to make the dataframe for test
group <- read.table("${groupfile}",sep="\t", colClasses="character")
#group <- group[-1,]
gsamp <- group[,1]
g1 <- group[1,2]
g2 <- group[which(!group[,2] %in% g1),2]
g2 <- g2[1]
gsamp1=group[which(group[,2] %in% g1),1]
gsamp2=group[which(group[,2] %in% g2),1]
otu_data <- otu_data[,which(samp %in% gsamp)]
otu_data <- otu_data[apply(otu_data,1,function(x)length(unique(x))!=1),]
lendata <- nrow(otu_data)
samp <- samp[which(samp %in% gsamp)]
result <- matrix(nrow = nrow(otu_data),ncol = 6)
pvalue <- 1
for(i in 1:nrow(otu_data)){
  o1 <- as.numeric(as.vector(unlist(otu_data[i,which(samp %in% gsamp1)])))
  o2 <- as.numeric(as.vector(unlist(otu_data[i,which(samp %in% gsamp2)])))
  me1 <- mean(o1)
  me2 <- mean(o2)
  sd1 <- sd(o1)
  sd2 <- sd(o2)
  #tt <- t.test(o1,o2,var.equal = TRUE,alternative = "two.side",conf.level = 0.95)
  test <- "${choose_test}"
  if(test == "student"){
    tt <- t.test(o1,o2,var.equal = TRUE,alternative = "two.side",conf.level = 0.95)
  }else if(test == "welch"){
    tt <- t.test(o1,o2,var.equal = FALSE,alternative = "two.side",conf.level = 0.95)
  }else{
    tt <- wilcox.test(o1,o2,alternative = "two.side",exact = F,conf.level = 0.95)
  }
  statlist <- tt$statistic
  statistic <- statlist[[1]]
  pvalue <- c(pvalue,tt$p.value)
  result[i,] = c(rownames(otu_data)[i],me1,sd1,me2,sd2,statistic)
}
pvalue <- pvalue[-1]
pvalue <- p.adjust(as.numeric(pvalue),method = "none")
qv <- p.adjust(as.numeric(pvalue),method = "fdr")
for(i in 1:length(pvalue)){
  pvalue[i] <- signif(pvalue[i],4)
  qv[i] <- signif(qv[i],4)
}
result <- cbind(result,pvalue)
result <- cbind(result,qv)
colnames(result) <- c(" ",paste("mean(",g1,")",sep=''),paste("sd(",g1,")",sep=''),paste("mean(",g2,")",sep=''),paste("sd(",g2,")",sep=''),"statistic","pvalue","corrected_pvalue")
result_order <- result[order(-(as.numeric(result[,2])+as.numeric(result[,4]))),]
if(lendata == 1){
  a <- data.frame(result_order)
  result_order <- t(a)
}
write.table(result_order,"${outputfile}",sep="\t",col.names=T,row.names=F,quote = F)
