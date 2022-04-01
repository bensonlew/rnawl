library("ropls")
mymethod <- paste("${scr_dir}", "ropls_ellipse.r",sep="/")
source(mymethod)

df <- read.table("${inputfile}",sep="\t",header=T,check.names = F,row.names = 1)
gf <- read.table("${groupfile}",sep="\t",header=F,check.names = F,row.names = 1)
classFc <- as.factor(gf[,1])
gsamp <- rownames(gf)
samp <- colnames(df)
df <- df[,which(samp %in% gsamp)]

df <- t(df)
df <- df[rownames(gf),]
df <- df[,colSums(df)!=0]
ci <- strsplit("${ci}",";")[[1]]
permutation <- strsplit("${perm}",";")[[1]]
mul_type = strsplit("${mul_type}",";")[[1]]
scale = strsplit("${data_trans}",";")[[1]]

if(length(gsamp) <7){
crossvalI = length(gsamp)
}else{
crossvalI = 7 # package default
}
# for plsda ,oplsda result table
get_pls_result <- function(pls_result, mytype, out_prefix,confidence){
        pls <- pls_result
        out_prefix <- paste(out_prefix,mytype,sep="/")
        name.pls.x <- paste(out_prefix,"sites.xls",sep=".")
        name.pls.loading <- paste(out_prefix,"loading.xls",sep=".")
        name.pls.vip <- paste(out_prefix,"vip.xls",sep=".")
        name.pls.sum <- paste(out_prefix,"model.xls",sep=".")
        name.pls.permMN<- paste(out_prefix,"permMN.xls",sep=".")
        name.pls.intercept <- paste(out_prefix,"intercept.xls",sep=".")
        name.pls.ellipse <- paste(out_prefix,"ellipse.xls",sep=".")
        pls.perMN <- pls@suppLs$permMN[,c(2,3,7)]
        pls.all_model <- (pls@modelDF)[,1:6]
        if(mytype=="PLS-DA"){
            pls.x <- getScoreMN(pls)
        }else{
            pls.p <- getScoreMN(pls)
            pls.o <- pls@orthoScoreMN
            pls.x <- cbind(pls.p,pls.o)
        }

        if(mytype=="OPLS-DA"){
            pls.loading <- cbind(getLoadingMN(pls),pls@orthoLoadingMN)  #############20200316 zgq
            loadingMN = pls@loadingMN
            corlist = list()
            for(i in 1:nrow(pls@loadingMN)){
                 corlist[i] = cor(pls@scoreMN,pls@suppLs$xModelMN[,i])
            }
            splot = cbind(loadingMN, corlist)
            write.table(splot, paste(out_prefix,"splot.xls",sep="."), sep='\t',quote=F)
        }else{
            pls.loading <- getLoadingMN(pls)
        }

        pls.sum <- pls@modelDF
        pls.vip <- getVipVn(pls)
        pls.vip <- as.data.frame(pls.vip)
        colnames(pls.vip) <- "VIP"
        pls.z1 <- lm(pls.perMN[,1]~pls.perMN[,3])$coefficients[1]
        pls.z2 <- lm(pls.perMN[,2]~pls.perMN[,3])$coefficients[1]
        pls.z <- as.matrix(c(pls.z1,pls.z2))
        ellipse <- add_ellipse(pls,classFc,confidence,parCompVi=c(1,2))
        write.table(pls.x,name.pls.x,sep="\t",quote=F,col.names=NA)
        write.table(pls.loading,name.pls.loading,sep="\t",quote=F,col.names=NA)
        write.table(pls.vip,name.pls.vip,sep="\t",quote=F,col.names=NA)
        write.table(pls.sum,name.pls.sum,sep="\t",quote=F,col.names=NA)
        write.table(pls.z,name.pls.intercept,sep="\t",quote=F)
        write.table(ellipse,name.pls.ellipse,sep="\t",row.names=F,quote=F,col.names=F)
        write.table(pls.perMN,name.pls.permMN,sep="\t",row.names=F,quote=F)
}

# from mul_type get ci,perm,sacle value
get_function_var <- function(m_type,myvar,is_numeric=T){
    if(is_numeric){
    result <- as.numeric(myvar[which(mul_type == m_type )])
    }else{
     result <- myvar[which(mul_type == m_type )]
    }
    return(result)
}

# from scale abbreviation to scale method
get_scale <- function(abbreviation){
    if(abbreviation == "UV"){
        scale <- "standard"
    }else if(abbreviation =="Ctr"){
        scale <- "center"
    }else if(abbreviation == "Par"){
        scale <- "pareto"
    }else{
        scale <- "none"
    }
    return(scale)
}

if ("pca"  %in% mul_type){
        confidence <- get_function_var("pca",ci)
        trans <- get_function_var("pca",scale,is_numeric=F)
        pca <- opls(df,printL=F,plotL=F,predI=NA,scaleC=get_scale(trans),crossvalI=crossvalI, algoC="nipals")
        if(pca@summaryDF[["pre"]]=="1"){
                pca <- opls(df,printL=F,plotL=F,predI=2,scaleC=get_scale(trans),crossvalI=crossvalI, algoC="nipals")
        }
        name.pca.x <- paste("${output}","PCA.sites.xls",sep="/")
        name.pca.loading <- paste("${output}","PCA.loading.xls",sep="/")
        name.pca.sum <- paste("${output}","PCA.model.xls",sep="/")
        name.pca.ellipse <- paste("${output}","PCA.ellipse.xls",sep="/")
        pca.x <- getScoreMN(pca)
        pca.loading <- getLoadingMN(pca)
        pca.sum <- pca@modelDF
        ellipse <- add_ellipse(pca,classFc,confidence,parCompVi=c(1,2))
        write.table(pca.x,name.pca.x,sep="\t",quote=F,col.names=NA)
        write.table(pca.loading,name.pca.loading,sep="\t",quote=F,col.names=NA)
        write.table(pca.sum,name.pca.sum,sep="\t",quote=F,col.names=NA)
        write.table(ellipse,name.pca.ellipse,sep="\t",row.names=F,quote=F,col.names=F)
}
if("plsda" %in% mul_type){
        confidence <- get_function_var("plsda",ci)
        trans <- get_function_var("plsda",scale,is_numeric=F)
        perm <- get_function_var("plsda",permutation)
        plsda <- opls(df,classFc,printL=F,plotL=F,predI=2,scaleC=get_scale(trans),permI=perm,crossvalI=crossvalI)
        if(plsda@modelDF$Signif.[2] == "NS" | plsda@modelDF$Signif.[2] == "N4"| plsda@modelDF$Signif.[1] == "NS"|plsda@modelDF$Signif.[1] == "N4"){
            plsda <- opls(df,classFc,printL=F,plotL=F,predI=2,scaleC=get_scale(trans),permI=perm,crossvalI=crossvalI)
        }else{
        plsda <- opls(df,classFc,printL=F,plotL=F,predI=NA,scaleC=get_scale(trans),permI=perm,crossvalI=crossvalI)
        }
        get_pls_result(plsda,"PLS-DA","${output}",confidence)
}
if("oplsda" %in% mul_type){
        confidence <- get_function_var("oplsda",ci)
        trans <- get_function_var("oplsda",scale,is_numeric=F)
        perm <- get_function_var("oplsda",permutation)
        oplsda <- opls(df, classFc, predI=1, orthoI=1,printL=F,plotL=F,scaleC=get_scale(trans),permI=perm,crossvalI=crossvalI)
        if(oplsda@modelDF[[1,"Signif."]] !="NS" & oplsda@modelDF[[2,"Signif."]] !="NS" & oplsda@modelDF[[1,"Signif."]] !="N4" & oplsda@modelDF[[2,"Signif."]] !="N4" ){
                oplsda <- opls(df,classFc,predI=1,orthoI=NA,printL=F,plotL=F,scaleC=get_scale(trans),permI=perm,crossvalI=crossvalI)
        }
        get_pls_result(oplsda,"OPLS-DA","${output}",confidence)
}


