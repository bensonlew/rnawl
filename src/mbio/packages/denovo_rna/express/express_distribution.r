# __author__ = "qiuping"
# last_modify:20170110
# -*- coding: utf-8 -*-
fpkm = read.table("${input_matrix}", header=TRUE, check.names=0, row.names=1 )
fpkm2 = log(fpkm+1)/log(2)
fpkm10 = log(fpkm+1)/log(10)
fpkm = fpkm
filename = "${filename}"
#fpkm = log(fpkm)/log(2)
#maxfpkm = 20
#samples=colnames(fpkm)
fpkm_density<-function(fpkm_data, logValue){ #add functions by khl 20170403
    samples=colnames(fpkm_data)
    maxfpkm = as.integer(max(fpkm_data)) + 2
    dens=c()
    for(i in 1:length(samples)){
       a=fpkm_data[i]
       den=density(a[,1],n = 601,from = -10,to = maxfpkm)
       dens=cbind(dens,den$y)
    }
    dens_new=cbind(den$x,dens)
    if(logValue == "None"){
       samples = c("fpkm", samples)
    }else{
       samples = c(paste("log", as.character(logValue), "fpkm", sep=""), samples)
    }
    colnames(dens_new)=samples
    return(dens_new)
}
fpkmlog2=fpkm_density(fpkm2,"2")
fpkmlog10=fpkm_density(fpkm10,'10')
fpkm=fpkm_density(fpkm,'None') 
#choose: fpkm>=0.01
fpkmlog2_new = rbind(fpkmlog2[1,], fpkmlog2[-1,][fpkmlog2[-1,1]>=0.1,])
fpkmlog10_new = rbind(fpkmlog10[1,], fpkmlog10[-1,][fpkmlog10[-1,1]>=0.1,])
fpkm_new = rbind(fpkm[1,], fpkm[-1,][fpkm[-1,1]>=0.1,])
write.table(fpkmlog2_new, file=paste("${outputfile}","/log2",filename,"_distribution.xls", sep=""), sep='\t', quote=F, row.names=F, col.names=T)
write.table(fpkmlog10_new, file=paste("${outputfile}","/log10",filename, "_distribution.xls", sep=""), sep='\t', quote=F, row.names=F, col.names=T)
write.table(fpkm_new, file=paste("${outputfile}","/",filename,"_distribution.xls", sep=""), sep='\t', quote=F, row.names=F, col.names=T)
#dens_new<-cbind(den$x,dens)
#samples=c('log2fpkm',samples)
#colnames(dens_new)=samples
#write.table(dens_new, file="/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tofiles/", sep='\t', quote=F, row.names=F, col.names=T)
