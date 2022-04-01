data = read.table("${exp}", header=T, row.names=1, sep="\t", com='', na.strings = 'NA')
col_ordering_a = c(${control})
col_ordering_b = c(${other})
test = "${method}"
alternative = "${alternative}"
ci = 0.95
p_correct = "${p_correct}"
cutoff_a = as.numeric("${cutoff_a}")
cutoff_b = as.numeric("${cutoff_b}")
average = "${average}"

up = ${up}
down = ${down}
p_cutoff = ${p_cutoff}

options(warn=-1)

data[data=='-']<-0
data[data=='_']<-0
data[is.na(data)] <- 0

otu_data <- data[,c(col_ordering_a, col_ordering_b)]
otu_data1 <- data[,col_ordering_a]
if (length(cutoff_a)==1){otu_data1 <- data[col_ordering_a]}
otu_data2 <- data[,col_ordering_b]
if (length(cutoff_b)==1){otu_data2 <- data[col_ordering_b]}

col_ordering_a <- seq(length(col_ordering_a))
col_ordering_b <- seq(length(col_ordering_a)+1, length(col_ordering_a) + length(col_ordering_b))

lendata <- nrow(otu_data)
exp_11 <- data.frame()
exp_01 <- data.frame()
exp_10 <- data.frame()
exp_00 <- data.frame()
for (i in seq(lendata)) {
    if (length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))!=0)])>=cutoff_a && length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))!=0)])>=cutoff_b){
    exp_11=rbind(exp_11,otu_data[i,])
    }
    else if (length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))!=0)])<cutoff_a && length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))!=0)])>=cutoff_b) {exp_01=rbind(exp_01,otu_data[i,])}
    else if (length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))!=0)])>=cutoff_a && length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))!=0)]) < cutoff_b) {exp_10=rbind(exp_10,otu_data[i,])}
    else if (length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))!=0)])<cutoff_a && length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))!=0)])<cutoff_b){exp_00=rbind(exp_00,otu_data[i,])}
}

#for (i in seq(lendata)) {
#    if (length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))!=0)])>=length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))==0)]) && length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))!=0)])>=length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))==0)])){
#    exp_11=rbind(exp_11,otu_data[i,])
#    }
#    else if(length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))!=0)])<length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))==0)]) && length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))!=0)])>=length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))==0)])) {exp_01=rbind(exp_01,otu_data[i,])}
#    else if (length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))!=0)])>=length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))==0)]) && length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))!=0)])<length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))==0)])) {exp_10=rbind(exp_10,otu_data[i,])}
#    else if (length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))!=0)])<length(otu_data1[i,][which(as.numeric(as.matrix(otu_data1[i,]))==0)]) && length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))!=0)])<length(otu_data2[i,][which(as.numeric(as.matrix(otu_data2[i,]))==0)])){exp_00=rbind(exp_00,otu_data[i,])}
#}

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

maxfc <- 0
outfile11 <- matrix(nrow=nrow(exp_11),ncol=8)
for(i in 1:nrow(exp_11)){
    acc <- as.character(rownames(exp_11)[i])
    if (average=='true'){exp_11[i,col_ordering_a][exp_11[i,col_ordering_a]==0]<-sum(exp_11[i,col_ordering_a])/(length(exp_11[i,col_ordering_a])-1)}
    if (average=='true'){exp_11[i,col_ordering_b][exp_11[i,col_ordering_b]==0]<-sum(exp_11[i,col_ordering_b])/(length(exp_11[i,col_ordering_b])-1)}

    exp_a <- as.numeric(format(exp_11[i,col_ordering_a]),scientific=TRUE)
    if (length(col_ordering_a) != 2){exp_a[exp_a==0]=NA}

    exp_b <- as.numeric(format(exp_11[i,col_ordering_b]),scientific=TRUE)
    if (length(col_ordering_b) != 2){exp_b[exp_b==0]=NA}

    avg_a <- mean(exp_a, na.rm=T)
    avg_b <- mean(exp_b, na.rm=T)
    fc <- avg_b/avg_a
    if (fc > maxfc) {maxfc <- fc}
    logfc <- log(fc,2)

    if (length(col_ordering_a) == 2) {
        if(exp_a[1]==exp_a[2]){
            exp_a[2] = exp_a[2]+(1e-10)
        }
    }

    if (length(col_ordering_b) == 2) {
        if(exp_b[1]==exp_b[2]){
            exp_b[2] = exp_b[2]+(1e-10)
        }
    }

    # t <- t.test(exp_a, exp_b, var.equal = T)
    if(test == "student"){
    tt <- t.test(exp_a,exp_b,var.equal = TRUE,alternative = alternative,conf.level = ci)
    }else if(test == "welch"){
    tt <- t.test(exp_a,exp_b,var.equal = FALSE,alternative = alternative,conf.level = ci)
    }else if(test == "wilcox"){
    tt <- wilcox.test(exp_a,exp_b,alternative = alternative,exact = F,conf.level = ci,paired = T)
    } else if (test == "chisq") {
    c1 <- as.numeric(as.vector(exp_a))
    c2 <- sum(as.numeric(as.vector(otu_data1[,1]))) - c1
    c3 <- as.numeric(as.vector(exp_b))
    c4 <- sum(as.numeric(as.vector(otu_data2[,1]))) - c3
    data1 <- matrix(c(c1,c2,c3,c4),ncol = 2)
    tt <- chisq.test(data1)
    }else if (test == "fisher") {
    c1 <- as.numeric(as.vector(exp_a))
    c2 <- sum(as.numeric(as.vector(otu_data1[,1]))) - c1
    c3 <- as.numeric(as.vector(exp_b))
    c4 <- sum(as.numeric(as.vector(otu_data2[,1]))) - c3
    data1 <- matrix(c(c1,c2,c3,c4),ncol = 2)
    tt <- fisher.test(data1)
    }

    #t <- t.test(exp_a, exp_b)
    p <- tt$p.value

    if (p < p_cutoff) {
        sig <- "yes"
    } else {
        sig <- "no"
    }

    if (fc > up) {
        reg <- "up"
    } else if (fc < down) {
        reg <- "down"
    } else {
        reg <- "no change"
    }
    outfile11[i,] <- c(acc,avg_a,avg_b,fc,logfc,p,sig,reg)
}

if (maxfc < 16) {maxfc <- 8}
outfile01 <- matrix(nrow=nrow(exp_01),ncol=8)
for(i in 1:nrow(exp_01)){
    acc <- as.character(rownames(exp_01)[i])
    if (average=='true'){exp_01[i,col_ordering_b][exp_01[i,col_ordering_b]==0]<-sum(exp_01[i,col_ordering_b])/(length(exp_01[i,col_ordering_b])-1)}

    exp_a <- as.numeric(format(exp_01[i,col_ordering_a]),scientific=TRUE)
    exp_a[exp_a==0]=NA

    exp_b <- as.numeric(format(exp_01[i,col_ordering_b]),scientific=TRUE)
    exp_b[exp_b==0]=NA

    avg_a <- 0
    avg_b <- mean(exp_b, na.rm=T)
    fc <- maxfc*2
    logfc <- log(fc,2)

    p <- 0.000000000000000000000001
    sig <- "yes"
    reg <- "up"

    outfile01[i,] <- c(acc,avg_a,avg_b,fc,logfc,p,sig,reg)
}

outfile10 <- matrix(nrow=nrow(exp_10),ncol=8)
for(i in 1:nrow(exp_10)){
    acc <- as.character(rownames(exp_10)[i])
    if (average=='true'){exp_10[i,col_ordering_a][exp_10[i,col_ordering_a]==0]<-sum(exp_10[i,col_ordering_a])/(length(exp_10[i,col_ordering_a])-1)}
    exp_a <- as.numeric(format(exp_10[i,col_ordering_a]),scientific=TRUE)
    exp_a[exp_a==0]=NA

    exp_b <- as.numeric(format(exp_10[i,col_ordering_b]),scientific=TRUE)
    exp_b[exp_b==0]=NA

    avg_a <- mean(exp_a, na.rm=T)
    avg_b <- 0
    fc <- 0.00001
    logfc <- signif(log(fc,2),4)

    p <- 0.000000000000000000000001
    sig <- "yes"
    reg <- "down"

    outfile10[i,] <- c(acc,avg_a,avg_b,fc,logfc,p,sig,reg)
}

outfile00 <- matrix(nrow=nrow(exp_00),ncol=8)
if (nrow(exp_00)>0){for(i in 1:nrow(exp_00)){
    acc <- as.character(rownames(exp_00)[i])


    avg_a <- 0
    avg_b <- 0
    fc <- 0
    logfc <- 0

    p <- 1
    sig <- "no"
    reg <- "no change"

    outfile00[i,] <- c(acc,avg_a,avg_b,fc,logfc,p,sig,reg)
}
}

header=c("Accession", "${control_name}", "${other_name}", "FC(${sympol})", "log2FC(${sympol})","Pvalue(${sympol})","significant","regulate")
colnames(outfile11)=header
corrected_pvalue<-p.adjust(outfile11[,6], method = p_correct)
outfile11<-cbind(outfile11,corrected_pvalue)
write.table(outfile11,file="${outxls}",col.names=T,row.names=F,quote=F,sep="\t")

colnames(outfile01)=header
corrected_pvalue<-rep('_',nrow(outfile01))
outfile01<-cbind(outfile01,corrected_pvalue)
write.table(outfile01,file="${other_name}_vs_${control_name}.01.xls",col.names=T,row.names=F,quote=F,sep="\t")

colnames(outfile10)=header
corrected_pvalue<-rep('_',nrow(outfile10))
outfile10<-cbind(outfile10,corrected_pvalue)
write.table(outfile10,file="${other_name}_vs_${control_name}.10.xls",col.names=T,row.names=F,quote=F,sep="\t")

colnames(outfile00)=header
corrected_pvalue<-rep('_',nrow(outfile00))
outfile00<-cbind(outfile00,corrected_pvalue)
write.table(outfile00,file="${other_name}_vs_${control_name}.00.xls",col.names=T,row.names=F,quote=F,sep="\t")
