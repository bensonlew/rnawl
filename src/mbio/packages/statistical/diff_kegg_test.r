protein_data <- read.table("${proteinfile}",sep = "\t",comment.char = '', colClasses="character")
enzyme_data <- read.table("${enzymefile}", sep="\t", comment.char = "", colClasses="character")

# two group test
two_group <- function(otu_data="", gsamp1="", gsamp2="", g1="", g2="", samp="") {
    otu_data <- otu_data[apply(otu_data,1,function(x)any(x>0)),]
	result <- matrix(nrow = nrow(otu_data), ncol=2)
	# pvalue <- 1
	test <- "${choose_test}"
	for (i in 1:nrow(otu_data)){
		if(test != "signal"){
			o1 <- as.numeric(as.vector(unlist(otu_data[i,which(samp %in% gsamp1)])))
			o2 <- as.numeric(as.vector(unlist(otu_data[i,which(samp %in% gsamp2)])))
		}else{
			o1 = numeric()
			o2 = numeric()
			for(j in 1:length(gsamp1)){
				o1[j] <- otu_data[i, gsamp1[j]]
			}
			for(j in 1:length(gsamp2)){
				o2[j] <- otu_data[i, gsamp2[j]]
			}
		}
		if(test == "student"){
			tt <- t.test(o1,o2, var.equal=TRUE, alternative = "${test_type}",conf.level = ${ci})
		}else if(test == "welch"){
			tt <- t.test(o1,o2,var.equal = FALSE,alternative = "${test_type}",conf.level = ${ci})
		}else if(test == "signal"){
			tt <- wilcox.test(o1,o2,alternative = "${test_type}",exact = F,conf.level = ${ci},paired = T)
		}else{
			tt <- wilcox.test(o1,o2,alternative = "${test_type}",exact = F,conf.level = ${ci})
		}
		# pvalue <- c(pvalue, tt$p.value)
		result[i,] <- c(rownames(otu_data)[i], tt$p.value)
	}
	# pvalue <- pvalue[-1]
	pvalue <- p.adjust(as.numeric(result[,2]),method = "none")
	qv <- p.adjust(as.numeric(pvalue),method = "${mul_test}")
	for(i in 1:length(pvalue)){
		pvalue[i] <- signif(pvalue[i],4)
		qv[i] <- signif(qv[i],4)
	}
	# result <- cbind(result,pvalue)
	# result <- cbind(result,qv)
	result[,2] <- pvalue
	result <- cbind(result, qv)
	colnames(result) <- c(" ", paste(g1, "-", g2, "_pvalue",sep=""), paste(g1, "-", g2, "_qvalue", sep=""))
	return(result)
}

# mean and sd
summary <- function(otu_data="", gsamp1="", gsamp2="", g1="", g2=""){
	result <- matrix(nrow = nrow(otu_data), ncol=5)
	for (i in 1:nrow(otu_data)){
		o1 = numeric()
		o2 = numeric()
		for(j in 1:length(gsamp1)){
			o1[j] <- otu_data[i, gsamp1[j]]
		}
		for(j in 1:length(gsamp2)){
			o2[j] <- otu_data[i, gsamp2[j]]
		}
		if(nrow(otu_data) == 1){
			o1 = o1/o1
			o2 = o2/o2
		}
		# sum1 <- as.numeric(summary(o1))[-4]
		# sum2 <- as.numeric(summary(o2))[-4]
		# for (l in 1:length(sum1)){
		# 	sum1[l] <- signif(sum1[l], 4)
		# 	sum2[l] <- signif(sum2[l], 4)
		# }
		me1 <- signif(mean(o1) * 100, 4)
		me2 <- signif(mean(o2) * 100, 4)
		sd1 <- signif(sd(o1) * 100, 4)
		sd2 <- signif(sd(o2) * 100, 4)
		result[i, ] = c(rownames(otu_data)[i], me1, sd1, me2, sd2)
	}
	colnames(result) <- c("", g1, paste(g1, "_sd", sep=""), g2, paste(g2, "_sd", sep=""))
	return(result)
}

# group test group by group
multiple_compare <- function(otu_data="", filter="") {
	samp <- t(otu_data[1,-1])
	otu_data <- otu_data[-1,]
	rownames(otu_data) <- otu_data[,1]
	otu_data <- otu_data[,-1]
	colnames(otu_data) <- samp
	data <- otu_data
	data <- data[, which(samp %in% gsamp)]
	samp <- samp[samp %in% gsamp]
	otu_data <- apply(data,2,function(x) as.numeric(x)/sum(as.numeric(x)))
	otu_data[is.na(otu_data)] <- 0
	rownames(otu_data) <- rownames(data)
	row_name <- rownames(data)
	otu_data <-  otu_data[which(row_name %in% filter),]
	# otu_data <- otu_data[apply(otu_data,1,function(x) any (x>0)),]
	summary_result <- FALSE
	test_result <- FALSE
	for(i in 1:(length(group_name)-1)){
		for(j in (i+1):length(group_name)){
			g1 <- group_name[i]
			g2 <- group_name[j]
			gsamp1 <- group[which(group[,2] %in% group_name[i]),1]
			gsamp2 <- group[which(group[,2] %in% group_name[j]),1]
			tmp_test <- two_group(otu_data=otu_data, gsamp1=gsamp1, gsamp2=gsamp2, g1=g1, g2=g2, samp=samp)
			tmp_summary <- summary(otu_data=otu_data, gsamp1=gsamp1, gsamp2=gsamp2, g1=g1, g2=g2)
			if(test_result != FALSE){
				test_result <- merge(x=test_result, y=tmp_test, all=TRUE)
			}else{
				test_result <- tmp_test
			}
			if(summary_result != FALSE){
				summary_result <- merge(x=summary_result, y=tmp_summary)
			}else{
				summary_result <- tmp_summary
			}
		}
	}
	if(test_result != FALSE){
		colnames(test_result)[1] <- c("Query")
	}
	if(summary_result != FALSE){
		colnames(summary_result)[1] <- c("Query")
	}
	result_data = merge(x=summary_result, y=test_result, by="Query", all=TRUE)
}

# read filterfile
# filterfile <- read.table("${filterfile}", sep = "\t", colClasses="character")
filterdata <- strsplit(readLines(file("${filterfile}", "r")), "\t")
profilter <- filterdata[[1]]
enzfilter <- filterdata[[2]]

#read groupfile to make the dataframe for test
group <- read.table("${groupfile}",sep="\t",colClasses="character")
gsamp <- group[,1]
group_name <- unique(group[,2])
enzyme_result = enzyme_have_result = FALSE
protein_result = protein_have_result = FALSE
enzyme_have_one = FALSE  # need add only one record condition
protein_have_one = FALSE
if(length(enzfilter) > 1){
    enzyme_have_result = TRUE
	enzyme_result <- multiple_compare(otu_data=enzyme_data, filter=enzfilter)
}else if(length(enzfilter) == 1){
	enzyme_have_result = TRUE
	enzyme_have_one = TRUE
}
if(length(profilter) > 1){
    protein_have_result = TRUE
	protein_result <- multiple_compare(otu_data=protein_data, filter=profilter)
}else if(length(profilter) ==1){
	protein_have_result = FALSE
	protein_have_one = TRUE
}
if (enzyme_have_result && protein_have_result){
	table_result <- rbind(enzyme_result, protein_result)
}else if(enzyme_have_result){
	table_result <- enzyme_result
}else if(protein_have_result){
	table_result <- protein_result
}
write.table(table_result, "${outputfile}", sep="\t", col.names=T, row.names=F, quote=F)
