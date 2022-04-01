otu_data <- read.table("${inputfile}",sep = "\t",comment.char = '', colClasses="character")
otu_data[otu_data=='-']='0'
samp <- t(otu_data[1,-1])
otu_data <- otu_data[-1,]
rownames(otu_data) <- otu_data[,1]
otu_data <- otu_data[,-1]
colnames(otu_data) <- samp

#read groupfile to make the dataframe for test
group <- read.table("${groupfile}",sep="\t",colClasses="character")
# otu_data <- otu_data[,-1]
colnames(otu_data) <- samp
gsamp <- group[,1]
g1 <- group[1,2]
g2 <- group[which(!group[,2] %in% g1),2]
g2 <- g2[1]
gsamp1=group[which(group[,2] %in% g1),1]
gsamp2=group[which(group[,2] %in% g2),1]
tmp <- otu_data
otu_data <- otu_data[,which(samp %in% gsamp)]
otu_data1 <- tmp[which(samp %in% gsamp1)]
otu_data2 <- tmp[which(samp %in% gsamp2)]

cutoff_a = as.numeric("${cutoff_a}")
cutoff_b = as.numeric("${cutoff_b}")

lendata <- nrow(otu_data)
exp_11 <- data.frame()
exp_01 <- data.frame()
exp_10 <- data.frame()
exp_00 <- data.frame()
for (i in seq(lendata)) {
    # if (length(otu_data1[i,][which(as.numeric(otu_data1[i,])!=0)])>=length(otu_data1[i,][which(as.numeric(otu_data1[i,])==0)]) && length(otu_data2[i,][which(as.numeric(otu_data2[i,])!=0)])>=length(otu_data2[i,][which(as.numeric(otu_data2[i,])==0)])){
    # exp_11=rbind(exp_11,otu_data[i,])
    # }
    # else if(length(otu_data1[i,][which(as.numeric(otu_data1[i,])!=0)])<length(otu_data1[i,][which(as.numeric(otu_data1[i,])==0)]) && length(otu_data2[i,][which(as.numeric(otu_data2[i,])!=0)])>=length(otu_data2[i,][which(as.numeric(otu_data2[i,])==0)])) {exp_01=rbind(exp_01,otu_data[i,])}
    # else if (length(otu_data1[i,][which(as.numeric(otu_data1[i,])!=0)])>=length(otu_data1[i,][which(as.numeric(otu_data1[i,])==0)]) && length(otu_data2[i,][which(as.numeric(otu_data2[i,])!=0)])<length(otu_data2[i,][which(as.numeric(otu_data2[i,])==0)])) {exp_10=rbind(exp_10,otu_data[i,])}
    # else if (length(otu_data1[i,][which(as.numeric(otu_data1[i,])!=0)])<length(otu_data1[i,][which(as.numeric(otu_data1[i,])==0)]) && length(otu_data2[i,][which(as.numeric(otu_data2[i,])!=0)])<length(otu_data2[i,][which(as.numeric(otu_data2[i,])==0)])){exp_00=rbind(exp_00,otu_data[i,])}
    if (length(otu_data1[i,][which(as.numeric(otu_data1[i,])!=0)])>=cutoff_a && length(otu_data2[i,][which(as.numeric(otu_data2[i,])!=0)])>=cutoff_b){
    exp_11=rbind(exp_11,otu_data[i,])
    }
    else if(length(otu_data1[i,][which(as.numeric(otu_data1[i,])!=0)])<cutoff_a && length(otu_data2[i,][which(as.numeric(otu_data2[i,])!=0)])>=cutoff_b) {exp_01=rbind(exp_01,otu_data[i,])}
    else if (length(otu_data1[i,][which(as.numeric(otu_data1[i,])!=0)])>=cutoff_a && length(otu_data2[i,][which(as.numeric(otu_data2[i,])!=0)])<cutoff_b) {exp_10=rbind(exp_10,otu_data[i,])}
    else if (length(otu_data1[i,][which(as.numeric(otu_data1[i,])!=0)])<cutoff_a && length(otu_data2[i,][which(as.numeric(otu_data2[i,])!=0)])<cutoff_b){exp_00=rbind(exp_00,otu_data[i,])}
}

if(nrow(exp_11) > 1){
  da <- exp_11
  exp_11 <-apply(da,2,function(x) as.numeric(x))
  exp_11[is.na(exp_11)] <- 0
  rownames(exp_11)<-rownames(da)
  exp_11 <- exp_11[apply(exp_11,1,function(x)any(x>0)),]
}

if(nrow(exp_10) > 1){
  da <- exp_10
  exp_10 <-apply(da,2,function(x) as.numeric(x))
  exp_10[is.na(exp_10)] <- 0
  rownames(exp_10)<-rownames(da)
}

if(nrow(exp_01) > 1){
  da <- exp_01
  exp_01 <-apply(da,2,function(x) as.numeric(x))
  exp_01[is.na(exp_01)] <- 0
  rownames(exp_01)<-rownames(da)
}

if(nrow(exp_00) > 1){
  da <- exp_00
  exp_00 <-apply(da,2,function(x) as.numeric(x))
  exp_00[is.na(exp_00)] <- 0
  rownames(exp_00)<-rownames(da)
}
samp <- samp[which(samp %in% gsamp)]
#otu_data <- otu_data[apply(otu_data,1,function(x) (length(unique(x[which(samp %in% gsamp1)])) != 1 | length(unique(x[which(samp %in% gsamp2)])) != 1)),]
exp_11 <- exp_11[apply(exp_11,1,function(x) (length(unique(x[which(samp %in% gsamp1)])) != 1 | length(unique(x[which(samp %in% gsamp2)])) != 1)),]
#exp_10 <- exp_10[apply(exp_10,1,function(x) (length(unique(x[which(samp %in% gsamp1)])) != 1 | length(unique(x[which(samp %in% gsamp2)])) != 1)),]
#exp_01 <- exp_01[apply(exp_01,1,function(x) (length(unique(x[which(samp %in% gsamp1)])) != 1 | length(unique(x[which(samp %in% gsamp2)])) != 1)),]
#exp_00 <- exp_00[apply(exp_00,1,function(x) (length(unique(x[which(samp %in% gsamp1)])) != 1 | length(unique(x[which(samp %in% gsamp2)])) != 1)),]
result <- matrix(nrow = nrow(exp_11),ncol = 7)
box_result <- matrix(nrow = nrow(otu_data),ncol = 11)
pvalue <- 1
test <- "${choose_test}"
maxfc <- 0
for(i in 1:nrow(exp_11)){
  if(test != "signal"){
    o1 <- as.numeric(as.vector(unlist(exp_11[i,which(samp %in% gsamp1)])))
    o2 <- as.numeric(as.vector(unlist(exp_11[i,which(samp %in% gsamp2)])))
  }else{
    o1 = numeric()
    o2 = numeric()
    for(j in 1:length(gsamp1)){
      o1[j] <- exp_11[i, gsamp1[j]]
    }
    for(j in 1:length(gsamp2)){
      o2[j] <- exp_11[i, gsamp2[j]]
    }
  }
   if (length(o1) == 2) {
        if(o1[1]==o1[2]){
            o1[2] = o1[2]+(1e-10)
        }
    }
  if (length(o2) == 2) {
        if(o2[1]==o2[2]){
            o2[2] = o2[2]+(1e-10)
        }
    }
  if (length(o1) != 2){o1[o1==0]=NA}
  if (length(o2) != 2){o2[o2==0]=NA}
  # o1[o1==0]=NA
  # o2[o2==0]=NA
  if(test == "student"){
    tt <- try(t.test(o1,o2,var.equal = TRUE,alternative = "${test_type}",conf.level = ${ci}))
  }else if(test == "welch"){
    tt <- t.test(o1,o2,var.equal = FALSE,alternative = "${test_type}",conf.level = ${ci})
  }else if(test == "signal"){
    tt <- wilcox.test(o1,o2,alternative = "${test_type}",exact = F,conf.level = ${ci},paired = T)
  }else{
    tt <- wilcox.test(o1,o2,alternative = "${test_type}",exact = F,conf.level = ${ci})
  }
  if(nrow(otu_data) == 1){
     o1 = o1/o1
     o2 = o2/o2
  }
  sum1 <- as.numeric(summary(o1))[-4]
  sum2 <- as.numeric(summary(o2))[-4]
  for(l in 1:length(sum1)){
    sum1[l] <- signif(sum1[l],4)
    sum2[l] <- signif(sum2[l],4)
  }
  # me1 <- signif(mean(o1, na.rm=T)*100,4)
  # me2 <- signif(mean(o2, na.rm=T)*100,4)
  # sd1 <- signif(sd(o1)*100,4)
  # sd2 <- signif(sd(o2)*100,4)
  # me1 <- signif(mean(o1, na.rm=T),4)
  # me2 <- signif(mean(o2, na.rm=T),4)
  me1 <- round(mean(o1, na.rm=T),4)     # modified by zhangyitong on 20210924
  me2 <- round(mean(o2, na.rm=T),4)
  sd1 <- signif(sd(o1, na.rm=T),4)
  sd2 <- signif(sd(o2, na.rm=T),4)
  #pvalue <- c(pvalue,tt$p.value)
  if ('try-error' %in% class(tt)){
  pvalue <- c(pvalue,1)}else{pvalue <- c(pvalue,tt$p.value)}
  fc <- signif(me2/me1,4)
  if (fc > maxfc) {maxfc <- fc}
  log2fc <- signif(log(fc,2),4)
  result[i,] = c(rownames(exp_11)[i],me1,sd1,me2,sd2,fc,log2fc)
  # box_result[i,] = c(rownames(otu_data)[i],sum1[1],sum1[2],sum1[3],sum1[4],sum1[5],sum2[1],sum2[2],sum2[3],sum2[4],sum2[5])
}
pvalue <- pvalue[-1]
pvalue <- p.adjust(as.numeric(pvalue),method = "none")
qv <- p.adjust(as.numeric(pvalue),method = "${mul_test}")
for(i in 1:length(pvalue)){
  pvalue[i] <- signif(pvalue[i],4)
  qv[i] <- signif(qv[i],4)
}
result <- cbind(result,pvalue)
result <- cbind(result,qv)
result_order <- result[order(-(as.numeric(result[,2])+as.numeric(result[,4]))),]

if(maxfc <16) {maxfc <- 16}

if (nrow(exp_01)>0) {
    for(i in 1:nrow(exp_01)){
        o1 <- as.numeric(as.vector(unlist(exp_01[i,which(samp %in% gsamp1)])))
        o2 <- as.numeric(as.vector(unlist(exp_01[i,which(samp %in% gsamp2)])))
        me1 <- 0
        sd1 <- 0
        # me2 <- signif(mean(o2),4)
        me2 <- round(mean(o2),4)
        sd2 <- signif(sd(o2),4)
        fc <- maxfc*2
        log2fc <- log(fc,2)
        pvalue <- 0
        qv <- 0
        exp_01_1 <- c(rownames(exp_01)[i],me1,sd1,me2,sd2,fc,log2fc,pvalue,qv)
        exp_01_1 <- t(as.data.frame(exp_01_1))
        result_order <- rbind(result_order,exp_01_1)

}
}

if (nrow(exp_10)>0) {
    for(i in 1:nrow(exp_10)){
        o1 <- as.numeric(as.vector(unlist(exp_10[i,which(samp %in% gsamp1)])))
        o2 <- as.numeric(as.vector(unlist(exp_10[i,which(samp %in% gsamp2)])))
        # me1 <- signif(mean(o1),4)
        me1 <- round(mean(o1),4)
        sd1 <- signif(sd(o1),4)
        me2 <- 0
        sd2 <- 0
        fc <- 0.00001
        log2fc <- signif(log(fc,2),4)
        pvalue <- 0
        qv <- 0
        exp_10_1 <- c(rownames(exp_10)[i],me1,sd1,me2,sd2,fc,log2fc,pvalue,qv)
        exp_10_1 <- t(as.data.frame(exp_10_1))
        result_order <- rbind(result_order,exp_10_1)
}
}

if (nrow(exp_00)>0) {
    for(i in 1:nrow(exp_00)){
        o1 <- as.numeric(as.vector(unlist(exp_00[i,which(samp %in% gsamp1)])))
        o2 <- as.numeric(as.vector(unlist(exp_00[i,which(samp %in% gsamp2)])))
        me1 <- 0
        sd1 <- 0
        me2 <- 0
        sd2 <- 0
        fc <- 0
        log2fc <- 0
        pvalue <- 1
        qv <- 1
        exp_00_1 <- c(rownames(exp_00)[i],me1,sd1,me2,sd2,fc,log2fc,pvalue,qv)
        exp_00_1 <- t(as.data.frame(exp_00_1))
        result_order <- rbind(result_order,exp_00_1)
}
}

colnames(result_order) <- c(" ",paste(g1,"-mean",sep=''),paste(g1,"-sd",sep=''),paste(g2,"-mean",sep=''),paste(g2,"-sd",sep=''),"fc","log2fc","pvalue","corrected_pvalue")

if(lendata == 1){
  a <- data.frame(result_order)
  result_order <- t(a)
}
write.table(result_order,"${outputfile}",sep="\t",col.names=T,row.names=F,quote = F)
colnames(box_result) <- c(" ",paste("min(",g1,")",sep=''),paste("Q1(",g1,")",sep=''),paste("Median(",g1,")",sep=''),paste("Q3(",g1,")",sep=''),paste("max(",g1,")",sep=''),paste("min(",g2,")",sep=''),paste("Q1(",g2,")",sep=''),paste("Median(",g2,")",sep=''),paste("Q3(",g2,")",sep=''),paste("max(",g2,")",sep=''))