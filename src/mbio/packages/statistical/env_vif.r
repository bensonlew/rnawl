library('vegan',quietly=TRUE)
viflim<-as.integer("${viflim}")
output<-"${output_dir}"
method<-"${method}"

# preprocess of the spe and env table
spe<-read.table(file="${abund_table}",header=T,check.names=FALSE,sep="\t")
rownames(spe) <-spe[,1]
rownames(spe) <-sapply(rownames(spe),function(x) gsub("_*{.+}"," ",x,perl = TRUE))
#write.table(rownames(spe),"rowname.xls",sep="\t") 
spe <-spe[,-1]
spe_head <-colnames(spe)

env<-read.table(file="${env_table}",header=T,check.names=FALSE,sep="\t")
rownames(env) <-env[,1]
rownames(env) <-sapply(rownames(env),function(x) gsub("_*{.+}"," ",x,perl = TRUE)) 
env <-env[,-1]
env <- env[spe_head,]
env_head <-colnames(env)
    
#spe<-spe[,rownames(env)]

#### decide rda or cca to use ####
   #spe<-log1p(spe)
   #decorana<-decorana(t(spe))
   dca_normal_run <- function(){dca_analysis <<- decorana(t(spe))}
   dca_error_run <- function(){
       scale_data <- apply(t(spe), 2, function(x) as.numeric(x)/sum(as.numeric(x)))
       dca_analysis <<- decorana(scale_data)
   }
   tryCatch(dca_normal_run(), error=function(e){print("error info");print(e);dca_error_run()})
   sink(paste(output,"/DCA.txt",sep=""),append=FALSE)
   dca_analysis
   sink(file=NULL)

### env select ###
#   bioenv <- summary(bioenv(t(spe),env))
#   sink("bioenv.txt",append=FALSE)
#   bioenv
#   sink(file=NULL)

dca1_length<-substr(readLines(paste(output,"/DCA.txt",sep=""))[11],17,21)

if (method=="cca" ||(method=="" && (dca1_length>=3.5))){
    spe.cca<-cca(t(spe)~.,data=env)
    env.vif<-vif.cca(spe.cca)
    #env.vif<-round(env.vif,2)
    #env.vif<-na.omit(env.vif)
    allvif<-as.data.frame(env.vif)
    write.table(t(allvif),paste(output,"/raw_cca_vif.txt",sep=""),sep='\t',quote=F,row.names=F)    
    env.vif[is.na(env.vif)] <- viflim + 100000
    allvif<-as.data.frame(env.vif)
    while (length(allvif[allvif>=viflim])>0 && length(allvif[,])>2){
        env<-env[ ,-which.max(env.vif)]
        spe.cca<-cca(t(spe)~.,data=env)
        env.vif<-vif.cca(spe.cca)
        #env.vif<-round(env.vif,2)
        allvif<-as.data.frame(env.vif)
    }
    write.table(t(allvif),paste(output,"/final_cca_vif.txt",sep=""),sep='\t',quote=F,row.names=F)    
    
}else if(method=="rda" ||(method=="" && (dca1_length<3.5))){
### rda/cca attributes ####
    #spe.rda <-rda(decostand(spe, "hellinger")~., data=env) ##
    spe.rda2<-rda(t(spe)~.,data=env)
#    sink(paste(output,"/rda.txt",sep=""),append=FALSE)
#    spe.rda2
#    sink(file=NULL)
    
### env vif attributes ####    
    env.vif<-vif.cca(spe.rda2)
    #env.vif<-round(env.vif,2)
    #env.vif<-na.omit(env.vif)
    allvif<-as.data.frame(env.vif)
    write.table(t(allvif),paste(output,"/raw_rda_vif.txt",sep=""),sep='	',quote=F,row.names=F)    

    #judge the raw_vif 
    env.vif[is.na(env.vif)] <- viflim + 100000
    allvif<-as.data.frame(env.vif)
    while (length(allvif[allvif>=viflim])>0 && length(allvif[,])>2){
        env<-env[ ,-which.max(env.vif)]
        spe.rda2<-rda(t(spe)~.,data=env)
        env.vif<-vif.cca(spe.rda2)
        #env.vif<-round(env.vif,2)
        allvif<-as.data.frame(env.vif)
    }
    write.table(t(allvif),paste(output,"/final_rda_vif.txt",sep=""),sep='\t',quote=F,row.names=F)    
}
