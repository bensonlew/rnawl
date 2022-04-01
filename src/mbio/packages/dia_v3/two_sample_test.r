otu_data <- read.table("${inputfile}",sep = "\t",comment.char = '', colClasses="character")
otu_data[otu_data=='-']='0'
samp <- t(otu_data[1,-1])
otu_data <- otu_data[-1,]
rownames(otu_data) <- otu_data[,1]
otu_data <- otu_data[,-1]
colnames(otu_data) <- samp
s1 <- "${sample1}"
s2 <- "${sample2}"
otu_data <- otu_data[,which(samp %in% c(s1,s2))]
#otu_data <- otu_data[apply(otu_data,1,function(x)any(x != "0.0")),]
#otu_data <- otu_data[apply(otu_data,1,function(x)any(x != "0")),]
lendata <- nrow(otu_data)
pvalue <- 1
exp_11 <- data.frame()
exp_01 <- data.frame()
exp_10 <- data.frame()
exp_00 <- data.frame()
for(i in 1:lendata){
    if(as.numeric(otu_data[i,s1]) != 0 && as.numeric(otu_data[i,s2]) != 0)  {
    exp_11=rbind(exp_11,otu_data[i,])
    }
    else if(as.numeric(otu_data[i,s1]) == 0 && as.numeric(otu_data[i,s2]) != 0)  {
    exp_01=rbind(exp_01,otu_data[i,])
    }
    else if(as.numeric(otu_data[i,s1]) != 0 && as.numeric(otu_data[i,s2]) == 0)  {
    exp_10=rbind(exp_10,otu_data[i,])
    }
    else if(as.numeric(otu_data[i,s1]) == 0 && as.numeric(otu_data[i,s2]) == 0)  {
    exp_00=rbind(exp_00,otu_data[i,])
    }
}

maxfc <- 0
result <- matrix(nrow = nrow(exp_11),ncol = 5)
if (nrow(exp_11)>0) {for(i in 1:nrow(exp_11)){
  c1 <- as.numeric(as.vector(exp_11$"${sample1}")[i])
  c2 <- sum(as.numeric(as.vector(exp_11$"${sample1}"))) - c1
  c3 <- as.numeric(as.vector(exp_11$"${sample2}")[i])
  c4 <- sum(as.numeric(as.vector(exp_11$"${sample2}"))) - c3
  data <- matrix(c(c1,c2,c3,c4),ncol = 2)
  test <- "${choose_test}"

  # added for log transformation
  log_trans<- "${log_trans}"
  if (log_trans == 'log10'){
    log_data = log10(data)
  }else if(log_trans == "log2"){
    log_data = log2(data)
  }else{
    log_data=data
  }

  if (test == "chi") {
    tt <- chisq.test(log_data)
  }else{
    tt <- fisher.test(log_data,alternative = "${test_type}",conf.level = ${ci})
  }
  pvalue <- c(pvalue,tt$p.value)
  # pro1 <- c1 / sum(as.numeric(as.vector(otu_data$"${sample1}")))
  # pro2 <- c3 / sum(as.numeric(as.vector(otu_data$"${sample2}")))

  fc <- signif(c3/c1,4)
  if (fc > maxfc) {maxfc <- fc}
  log2fc <- signif(log(fc,2),2)
  result[i,] = c(rownames(exp_11)[i],c1,c3,fc,log2fc)
}
pvalue <- pvalue[-1]
pvalue <- p.adjust(as.numeric(pvalue),method = "none")
qvalue <- p.adjust(as.numeric(pvalue),method = "${mul_test}")
for(i in 1:length(pvalue)){
  pvalue[i] <- signif(pvalue[i],4)
  qvalue[i] <- signif(qvalue[i],4)
}
}
result <- cbind(result,pvalue)
result <- cbind(result,qvalue)
result_order <- result[order(-(as.numeric(result[,2])+as.numeric(result[,3]))),]

if(maxfc <16) {maxfc <- 16}

if (nrow(exp_01)>0) {for(i in 1:nrow(exp_01)){
   c1 <- as.numeric(as.vector(exp_01$"${sample1}")[i])
   c2 <- as.numeric(as.vector(exp_01$"${sample2}")[i])
   fc <- maxfc*2
   log2fc <- log(fc,2)
   pvalue <- 0
   qv <- 0
   exp_01_1 <- c(rownames(exp_01)[i],c1,c2,fc,log2fc,pvalue,qv)
   exp_01_1 <- t(as.data.frame(exp_01_1))
   result_order <- rbind(result_order,exp_01_1)
   }
}

if (nrow(exp_10)>0) {for(i in 1:nrow(exp_10)){
   c1 <- as.numeric(as.vector(exp_10$"${sample1}")[i])
   c2 <- as.numeric(as.vector(exp_10$"${sample2}")[i])
   fc <- 0.00001
   log2fc <- signif(log(fc,2),4)
   pvalue <- 0
   qv <- 0
   exp_10_1 <- c(rownames(exp_10)[i],c1,c2,fc,log2fc,pvalue,qv)
   exp_10_1 <- t(as.data.frame(exp_10_1))
   result_order <- rbind(result_order,exp_10_1)
   }
}

if (nrow(exp_00)>0) {for(i in 1:nrow(exp_00)){
   c1 <- as.numeric(as.vector(exp_00$"${sample1}")[i])
   c2 <- as.numeric(as.vector(exp_00$"${sample2}")[i])
   fc <- 0
   log2fc <- 0
   pvalue <- 1
   qv <- 1
   exp_00_1 <- c(rownames(exp_00)[i],c1,c2,fc,log2fc,pvalue,qv)
   exp_00_1 <- t(as.data.frame(exp_00_1))
   result_order <- rbind(result_order,exp_00_1)
   }
}

colnames(result_order) <- c(" ",paste(s1,"-propotion",sep=''),paste(s2,"-propotion",sep=''),"fc","log2fc","pvalue","corrected_pvalue")

if(lendata == 1){
  a <- data.frame(result_order)
  result_order <- t(a)
}
write.table(result_order,"${outputfile}",sep="\t",col.names=T,row.names=F,quote = F)
