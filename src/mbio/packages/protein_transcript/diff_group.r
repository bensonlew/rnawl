data = read.table("${exp}", header=T, row.names=1, sep="\t", com='', na.strings = 'NA')
col_ordering_a = c(${control})
col_ordering_b = c(${other})
test = "${method}"
alternative = "${alternative}"
ci = 0.95
p_correct = "${p_correct}"
up = ${up}
down = ${down}
p_cutoff = ${p_cutoff}


options(warn=-1)
outfile <- matrix(nrow=nrow(data),ncol=8)
for(i in 1:nrow(data)){
    acc <- as.character((rownames(data[i,])))

    exp_a <- as.numeric(format(data[i,col_ordering_a]),scientific=TRUE)
    exp_a[exp_a==0]=NA

    exp_b <- as.numeric(format(data[i,col_ordering_b]),scientific=TRUE)
    exp_b[exp_b==0]=NA

    avg_a <- mean(exp_a, na.rm=T)
    avg_b <- mean(exp_b, na.rm=T)
    fc <- avg_b/avg_a
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
    c2 <- sum(as.numeric(as.vector(data[,col_ordering_a]))) - c1
    c3 <- as.numeric(as.vector(exp_b))
    c4 <- sum(as.numeric(as.vector(data[,col_ordering_b]))) - c3
    data1 <- matrix(c(c1,c2,c3,c4),ncol = 2)
    tt <- chisq.test(data1)
    }else if (test == "fisher") {
    c1 <- as.numeric(as.vector(exp_a))
    c2 <- sum(as.numeric(as.vector(data[,col_ordering_a]))) - c1
    c3 <- as.numeric(as.vector(exp_b))
    c4 <- sum(as.numeric(as.vector(data[,col_ordering_b]))) - c3
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
    outfile[i,] <- c(acc,avg_a,avg_b,fc,logfc,p,sig,reg)

}
colnames(outfile)=c("Accession", "${control_name}", "${other_name}", "FC(${sympol})", "log2FC(${sympol})","Pvalue(${sympol})","significant","regulate")
corrected_pvalue<-p.adjust(outfile[,6], method = p_correct)
outfile<-cbind(outfile,corrected_pvalue)
write.table(outfile,file="${outxls}",col.names=T,row.names=F,quote=F,sep="\t")