arg<-commandArgs()
print(arg)
raw_data<-read.table(arg[6],header = T,row.names = 1,stringsAsFactors = F)
data_median<-sweep(raw_data,2,apply(raw_data,2,median,na.rm=T),FUN = "/")
FC<-data_median[,2]/data_median[,1]
log2FC <-log2(FC)
r0<-median(log2FC)
r1<-quantile(log2FC,0.8413)
r2<-quantile(log2FC,0.1587)
sd1<-r1-r0
sd2<-r0-r2
Z2<-(r0-log2FC)/sd2
Z1<-(log2FC-r0)/sd1
Z_score <- vector(mode="numeric",length= length(raw_data[,1]))
Z_P <- vector(mode = "numeric",length = length(raw_data[,1]))
for (i in 1:length(raw_data[,1])){
Z_score[i]<- max(Z1[i],Z2[i])
Z_P[i]<- 1-pnorm(Z_score[i])
}
diff_exp_tmp<-cbind(raw_data,FC,log2FC,Z_P)
df_re <- data.frame(c(rownames(diff_exp_tmp)),c(diff_exp_tmp[,3]))
colnames(df_re)<-c("Accession","FC_a")
df_up <- subset(df_re, (FC_a > as.numeric(arg[7])))
df_up <- cbind(df_up, rep("up", nrow(df_up)))
df_down <- subset(df_re, (FC_a < as.numeric(arg[8])))
df_down <- cbind(df_down, rep("down", nrow(df_down)))
df_no_change <- subset(df_re, (FC_a >= as.numeric(arg[8]) & FC_a <= as.numeric(arg[7])))
df_no_change <- cbind(df_no_change, rep("no_change", nrow(df_no_change)))
colnames(df_up)[3] <- "regulate"
colnames(df_down)[3] <- "regulate"
colnames(df_no_change)[3] <- "regulate"
df_sig <- data.frame(c(rownames(diff_exp_tmp)),c(diff_exp_tmp[,5]))
colnames(df_sig)<-c("Accession","Pvalue")
df_sig_y <- subset(df_sig, (Pvalue < arg[9]))
df_sig_y <- cbind(df_sig_y, rep("yes", nrow(df_sig_y)))
df_sig_n <- subset(df_sig, (Pvalue >= arg[9]))
df_sig_n <- cbind(df_sig_n, rep("no", nrow(df_sig_n)))
colnames(df_sig_y)[3] <- "significance"
colnames(df_sig_n)[3] <- "significance"

df_r <- rbind(df_up, df_down, df_no_change)
df_s <- rbind(df_sig_y,df_sig_n)
exp_re<-merge(df_r,df_s,by=c("Accession"),all=T)
new_diff_exp<- cbind(rownames(diff_exp_tmp),diff_exp_tmp)
rownames(new_diff_exp)<- NULL
colnames(new_diff_exp)[1]<- "Accession"
diff_exp_tmp2<- merge(new_diff_exp,exp_re,by=c("Accession"),all=T)
diff_exp<- diff_exp_tmp2[,c(1,2,3,4,5,6,10,8)]
colnames(diff_exp)<- c("Accession","control","sample","FC","log2FC","Pvalue","significance","regulate")
write.table(diff_exp,"sample_vs_control.diff.exp.xls",sep="\t",quote=FALSE,row.names = FALSE)