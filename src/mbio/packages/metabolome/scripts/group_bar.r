otu <- read.table('${inputfile}', sep = '\t',comment.char='', colClasses="character")
# otu <- read.table('C:\\Users\\ping.qiu.MAJORBIO\\Desktop\\otu_file.xls', sep = '\t',comment.char='')
sam <- t(otu[1,-1])
otu <- otu[-1,]
rownames(otu) <- otu[,1]
otu <- otu[,-1]
colnames(otu) <- sam
# group <- read.table('C:\\Users\\ping.qiu.MAJORBIO\\Desktop\\group_file_input.group.xls',sep = '\t')
group <- read.table('${groupfile}',sep = '\t', colClasses="character")
gsamp <- group[,1]
otu <- otu[,which(sam %in% gsamp)]
length <- nrow(otu)
#data <- apply(otu,2,function(x)signif(as.numeric(x)/sum(as.numeric(x)),4))
data <- apply(otu,2,function(x)signif(as.numeric(x),4))  #20200522
data[is.na(data)] <- 0
if(length == 1){
    data <- t(as.data.frame(data))
}
rownames(data) <- rownames(otu)
col <- colnames(data)
for(i in 1:(length(col))){
    col[i] <- paste(group[which(gsamp %in% col[i]),2],col[i],sep='-')
}
colnames(data) <- col
# write.table(data,"C:\\Users\\ping.qiu.MAJORBIO\\Desktop\\test5",sep="\t",col.names=T,row.names=T,quote = F)
write.table(data,"${outputfile}",sep="\t",col.names=T,row.names=T,quote = F)
