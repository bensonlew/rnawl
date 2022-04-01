options(warn=-1)
library(methods)
data = read.table("${exp}", header=T, row.names=1, sep="\t", quote='\"', na.strings = 'NA')
feature_num = nrow(data)
feature_names = vector()
for(i in 1:feature_num){
feat_name = paste0("Feature_",i)
feature_names <- append(feature_names,feat_name)}

data_des = data.frame(Feature=feature_names,Accession=rownames(data))
rownames(data) = NULL
rownames(data) = feature_names

group = read.table("${group}",sep="\t",quote='\"',na.strings="NA")
names(group)<- c("Sample","Group")
na_ratio = as.numeric("${na_ratio}")
specific = "${specific}"
specific_threshold = as.numeric("${specific_threshold}")
fill_type = '${fill_type}'
interact = '${interactive}'

method = "${method}"

minConc<-min(data[data>0], na.rm=T)/2;

Merge_func <- function(x,y) {
    df <- merge(x, y, by = "Feature", all.x = T, all.y = T)
    return(df)
}

Merge_rname <- function(x,y){
    x[["Feature"]] <- rownames(x)
    rownames(x)<-NULL
    y[["Feature"]] <- rownames(y)
    rownames(y)<- NULL
    df <- merge(x,y,by="Feature",all=FALSE)
    return(df)
}

nafunctions<-function(x,method="zero"){
    df<-df1<-as.data.frame(x)
    method<-tolower(method)
    if(method=="zero"){
      df[is.na(df)]<-0
    }
    else if(method=="min"){
      df[is.na(df)]<-min(df1,na.rm = TRUE)
    }
    else if(method=="colmedian"){
      library(e1071)
      df<-impute(df1,what ="median")
    }
    else if(method=="rowmedian"){
      library(e1071)
      dfx<-impute(t(df1),what ="median")
      df<-t(dfx)
    }
    else if(method=="mean"){
      library(e1071)
      dfx<-impute(t(df1),what ="mean")
      df<-t(dfx)
    }
    else if(method=="knnmethod"){
      library(impute)
      data_zero1<-impute.knn(as.matrix(df1),k = 10, rowmax = 0.9, colmax = 0.9)
      df<-data_zero1$data
    }
    else if(method=="seqknn"){
      library(SeqKnn)
      df <- SeqKNN(df1,k = 10)
    }
    else if(method=="bpca"){
      library(pcaMethods)
      data_zero1<-pcaMethods::pca(as.matrix(df1), nPcs = ncol(df1)-1, method = "bpca", maxSteps =100)
      df<-completeObs(data_zero1)
    }
    else if(method=="svdmethod"){
      library(pcaMethods)
      data_zero1<-pcaMethods::pca(as.matrix(df1), nPcs = ncol(df1)-1, method = "svdImpute")
      df<-completeObs(data_zero1)
    }
    else if(method=="lls"){
      library(pcaMethods)
      data_zero1<-llsImpute(t(df1), k = 10)
      df<-t(completeObs(data_zero1))
    }
    else if(method=="mle"){
      library(norm)
      xxm<-as.matrix(df1)
      ss <- norm::prelim.norm(xxm)
      thx <- norm::em.norm(ss)
      norm::rngseed(123)
      df <- norm::imp.norm(ss, thx, xxm)
    }
    else if(method=="qrilc"){
      library(imputeLCMD)
      xxm<-t(df1)
      data_zero1 <- imputeLCMD::impute.QRILC(xxm, tune.sigma = 1)[[1]]
      df<-t(data_zero1)
    }
    else if(method=="mindet"){
      library(imputeLCMD)
      xxm<-as.matrix(df1)
      df <- imputeLCMD::impute.MinDet(xxm, q = 0.01)
    }
    else if(method=="minprob"){
      library(imputeLCMD)
      xxm<-as.matrix(df1)
      df <- imputeLCMD::impute.MinProb(xxm, q = 0.01, tune.sigma = 1)
    }
    else if(method=="irm"){
      library(VIM)
      df <- irmi(df1, trace = TRUE,imp_var=FALSE)
      rownames(df)<-rownames(df1)
    }
    else if(method=="impseq"){
      library(rrcovNA)
      df <- impSeq(df1)
    }
    else if(method=="impseqrob"){
      library(rrcovNA)
      data_zero1 <- impSeqRob(df1, alpha=0.9)
      df<-data_zero1$x
    }
    else if(method=="mice-norm"){
      library(mice)
      minum<-5
      datareadmi<-mice(df1,m=minum,seed = 1234, method ="norm")
      newdatareadmi<-0
      for (i in 1:minum) {
        newdatareadmi<-complete(datareadmi,action = i)+newdatareadmi
      }
      df<-newdatareadmi/minum
      rownames(df)<-rownames(df1)
    }
    else if(method=="mice-cart"){
      library(mice)
      minum<-5
      datareadmi<-mice(df1,m=minum,seed = 1234, method ="cart")
      newdatareadmi<-0
      for (i in 1:minum) {
        newdatareadmi<-complete(datareadmi,action = i)+newdatareadmi
      }
      df<-newdatareadmi/minum
      rownames(df)<-rownames(df1)
    }
    else if(method=="trknn"){
      source('Trunc_KNN/Imput_funcs.r')
      sim_trKNN_wrapper <- function(data) {
        result <- data %>% as.matrix %>% t %>% imputeKNN(., k=10, distance='truncation', perc=0) %>% t
        return(result)
      }
      df1x <- sim_trKNN_wrapper(t(df1))
      df<-as.data.frame(t(df1x))
    }
    else if(method=="rf"){
      library(missForest)
      data_zero1 <- missForest(t(df1), maxiter =10,ntree = 500,mtry=floor(nrow(df1)^(1/3)),verbose = TRUE)
      df<-t(data_zero1$ximp)
    }
    else if(method=="pi"){
      width <- input$piwidth
      downshift <- input$pidownshift
      for(i in 1:ncol(df1)){
        temp <- df1[[i]]
        if(sum(is.na(temp))>0){
          temp.sd <- width * sd(temp[!is.na(temp)], na.rm = TRUE)
          temp.mean <- mean(temp[!is.na(temp)], na.rm = TRUE) - downshift * sd(temp[!is.na(temp)], na.rm = TRUE)
          n.missing <- sum(is.na(temp))
          temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
          df[[i]]<-temp
        }
      }
      df
    }
    else if(method=="grr"){
      library(DreamAI)
      df<-impute.RegImpute(data=as.matrix(df1), fillmethod = "row_mean", maxiter_RegImpute = 10,conv_nrmse = 1e-03)
    }
    else if(method=="gms"){
      library(GMSimpute)
      df<-GMS.Lasso(df1,nfolds=3,log.scale=FALSE,TS.Lasso=TRUE)
    }
    else if (method =='none' | method == ''){
        df<-as.data.frame(df)
    }
    else{
      stop("Unspported methods so far~~")
    }
    df<-as.data.frame(df)
    df
  }

extract_suijidataout<-function(data){
  dfx<-data ##为根据CV过滤后的数据；
  naratiox<-dim(dfx[!complete.cases(dfx),])[1]/nrow(dfx)  ##有缺失值的行占总行数的比值；
  nacolratio<-apply(dfx[!complete.cases(dfx),],2,function(x){sum(is.na(x))})/nrow(dfx[!complete.cases(dfx),])   ##每一列的缺失值比例
  nacolratio[is.na(nacolratio)]<-0  # added by zhangyitong on 20211221 allow complete exp matrix
  datanonaoutx<-dfx[complete.cases(dfx),]  ##没有na的数据矩阵；
  nanum<-round(naratiox*nrow(datanonaoutx))  #有缺失值的行数； 他是按照行和列分别来进行的随机替换和分析的；
  set.seed(123)
  samplenaindex<-sample(1:nrow(datanonaoutx),nanum)
  datanonaoutxx1<-datanonaoutx
  for(i in 1:ncol(datanonaoutx)){
    set.seed(i)
    samplenaindexi<-sample(samplenaindex,round(nacolratio[i]*nanum))
    datanonaoutxx1[samplenaindexi,i]<-NA
  }
  if(nrow(datanonaoutxx1)!=sum(complete.cases(datanonaoutxx1))){
    datanonaoutxx1_1<-datanonaoutx[samplenaindex,]
    datanonaoutxx1_2<-datanonaoutx[-samplenaindex,]
    for(ii in 1:nrow(datanonaoutxx1_1)){
      set.seed(ii)
      eachrownaratio<-round(naratiox*ncol(datanonaoutx))
      samplenaindexi<-sample(1:ncol(datanonaoutx),ifelse(eachrownaratio==0,round(0.3*ncol(datanonaoutx)),eachrownaratio))
      datanonaoutxx1_1[ii,samplenaindexi]<-NA
    }
    datanonaoutxx1<-rbind(datanonaoutxx1_1,datanonaoutxx1_2)
    datanonaoutxx1<-datanonaoutxx1[match(rownames(datanonaoutx), rownames(datanonaoutxx1)),]
  }
  rowcolindex1<-which(is.na(datanonaoutxx1),arr.ind = TRUE)   ##缺失值的位置；
  phosdata_na_sim1<-datanonaoutxx1
  suijidataresult <- list(suijinadatadf=phosdata_na_sim1,rowcolindex=rowcolindex1)
  return(suijidataresult)
}

calculate_nrmsedfout<-function(data,method="bpca"){  # data 是含有缺失值的；计算缺失值，填充后与真实值的差距；这个data是原始去掉Na过多和CV 阈值的data;
    dfx<-data	#dim(dfx) :3640   20
	fillmethod <- method
	suijidataout <- extract_suijidataout(data)
    datanonaoutx<-dfx[complete.cases(dfx),]   #完整数据集； dim(datanonaoutx): 3330   20  -不包含缺失值；
    suijinadatadfx<-suijidataout$suijinadatadf  #随机抽样的数据集； dim(suijinadatadfx): 3330   20  --包含缺失值；
    rowcolindexx<-suijidataout$rowcolindex   #缺失值的位置索引
    suijidatajieguo<- nafunctions(suijinadatadfx,method=fillmethod)   #dim(suijidatajieguo): 3330   20 --缺失值填充后的数据矩阵
    impdata<-as.numeric(suijidatajieguo[rowcolindexx])  #填充的缺失值
    truedata<-as.numeric(datanonaoutx[rowcolindexx])   #相同位置的真实值；
    nrmsejisuan<-sqrt(mean((impdata - truedata)^{2})/var(truedata))
	NRMSE=round(nrmsejisuan,digits=5)
	return(NRMSE)
    }


calculate_cvdata<-function(data,group){
    dfx<- data
    samplesdf <- group
    grnames<-unique(samplesdf$Group)
    cvdf<-NULL
    cvdfnames<- vector()
    for(i in 1:length(grnames)){
        sample<-samplesdf[samplesdf$Group==grnames[i],1]
        datai<-dfx[,names(dfx)%in%sample]
        #datai<-dfx[,samplesdf$Group==grnames[i]]
        cvi<-apply(datai,1,function(x){
            #if(all(as.numeric(x)==0)|is.na(all(as.numeric(x)==0))){
            x<-as.numeric(x)
            x[x==0]=NA
            x<-as.numeric(na.omit(x))
            if(length(x)<2){
                NA
            }else{
                raster::cv(as.numeric(x),na.rm = TRUE,aszero=TRUE)
            }
        })
        cvdf<-cbind(cvdf,cvi)
        cvdfnames<- append(cvdfnames,paste0('cv_',grnames[i]))
    }
	cvg <- apply(data,1,calculate_groupcv)
	cvdf <- cbind(cvdf,cvg)
	cvdfnames<-append(cvdfnames,'cv_group')
    round(cvdf,digits=5)
    cvdf <- as.data.frame(cvdf)
    names(cvdf) <- cvdfnames
    return(cvdf)
}

summarise_cv<-function(data){   # summarising cv distribution within each group
  cv<-data
  result<-data.frame()
  frequency_name<-c()
  percent_name<-c()
  cumulative_name<-c()
  for (i in c(2:ncol(cv)-1)){
    all_protein<-cv[,i]
    na_count<-length(all_protein[is.na(all_protein)==TRUE])
    protein<-na.omit(all_protein)
    mean<-mean(protein)
    #range<-range(protein)
    quant<-seq(10,100,10)
    #quant<-data.frame(quantile(range,c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)))
    #quant<-as.data.frame(t(quant))
    freq<-c()
    percent<-c()
    cumulative<-c()
    for (j in c(1:length(quant))){
      if (j == 1){
        freq1<-length(protein[protein<quant[j]])
      }else if (j == length(quant)){
        freq1<-length(protein[protein>=quant[j-1]&protein<=quant[j]])
      }else{
        freq1<-length(protein[protein>=quant[j-1]&protein<quant[j]])
      }
      freq<-c(freq,freq1)
      percent<-c(percent,freq1/length(all_protein))
      cumulative<-c(cumulative,sum(percent))
      if (i == 2){
        frequency_name<-c(frequency_name,paste0(quant[j],"%_frequency"))
        percent_name<-c(percent_name,paste0(quant[j],"%_percent"))
        cumulative_name<-c(cumulative_name,paste0(quant[j],"%_cumulative"))
      }
    }
    exceed<-length(protein[protein>100])
    freq<-c(freq, exceed, na_count)
    percent<-c(percent,exceed/length(all_protein), na_count/length(all_protein))
    cumulative<-c(cumulative,sum(percent[1:length(percent)-1]), sum(percent))
    result_sample<-c(freq,percent,cumulative,mean)
    result<-rbind(result,result_sample)
  }
  rownames(result)<-colnames(cv)[2:ncol(cv)-1]
  frequency_name<-c(frequency_name,">100%_frequency","na_frequency")
  percent_name<-c(percent_name,">100%_percent", "na_percent")
  cumulative_name<-c(cumulative_name,">100%_cumulative", "na_cumulative")
  names<-c(frequency_name,percent_name,cumulative_name,"mean_cv")
  colnames(result)<-names
  return(result)
}

calculate_groupcv <- function(x){
	exp = as.numeric(x) #运用apply时，每一行是一个向量；
	#tmp_group = group[as.factor(names(x)) == as.vector(group[["Sample"]]),2]
	tmp_group = sapply(names(x), function(x) group[x==group$Sample,2])
	#print(tmp_group)
	#print(names(x))
	var_data= data.frame(exp=exp,group = tmp_group)
    var_data = var_data[complete.cases(var_data),]
    if (length(unique(var_data[['group']]))<2){
        cv_grp = NA
    }else{
        fit = aov(exp~group,data = var_data)
        result = summary(fit)[[1]]
        mean_exp = mean(exp,na.rm =T)
        mean_sq = result[1,][,"Mean Sq"]
        cv_grp = sqrt(mean_sq)/mean_exp*100
    }
	return(cv_grp)
}

count_nanumber<-function(data,group){
    dfx<- data
    samplesdf <- group
    grnames<-unique(samplesdf$Group)
    nadf<-NULL
    nadfnames<- vector()
    for(i in 1:length(grnames)){
        sample<-samplesdf[samplesdf$Group==grnames[i],1]
        datai<-dfx[,names(dfx)%in%sample]
        #datai<-dfx[,samplesdf$Group==grnames[i]]
        datai[is.na(datai)] = 0
        nai = rowSums(datai==0)
        nadf<-cbind(nadf,nai)
        nadfnames<- append(nadfnames,paste0(as.character(grnames[i])))
    }
    nadf <- as.data.frame(nadf)
    names(nadf) <- nadfnames
    dfx[is.na(dfx)] =0
    nadf[["NA_num"]] <- rowSums(dfx==0)
    return(nadf)
}

my_impute<- function(data_df,group_df,na_ratio_num,method,specific='yes',specific_threshold=0.5,fill_type='group'){
    df_list = list()
    raw_data <- data_df
    group <- group_df
    na_ratio <- na_ratio_num
    group_specific <- specific
    group_specific_threshold <- specific_threshold
    method = tolower(method)
    #首先计算缺失值的个数与CV值的大小；然后根据缺失值的个数进行判断；
    na_result <- count_nanumber(raw_data,group)
    na_result[["Feature"]] <- rownames(na_result)
    rownames(na_result)<-NULL
    #raw_data[is.na(raw_data)] = 0
    #raw_data[["data_na"]] = rowSums(raw_data==0)
    data <- raw_data[na_result["NA_num"] < na_ratio * (ncol(raw_data)),]    #去掉初始缺失值的完整数据矩阵；data 有 0 无缺失；
    print(ifelse(dim(data)[1]==0, 'OOPS', 'ALLGOOD'))
    # cv_result<- calculate_cvdata(data,group)
    # cv_result[["Feature"]] <- rownames(cv_result)
    # rownames(cv_result)<-NULL
    # cv_summary<-summarise_cv(cv_result)
    #cv_na_statis_result <- Merge_rname(cv_result,na_result)
    if (group_specific =="yes"){	#需要判断是否为特异的，因此需要按组进行填充；对数据进行处理；
        nrmses = vector()
        nrmses_basic = vector()
        for(i in unique(group[,2])){
            samples <- group[group[,2] == i,1]
            exp_grp <- data[,names(data) %in% samples]
            exp_grp[is.na(exp_grp)] = 0
            exp_grp[["data_na"]] = rowSums(exp_grp==0)  #计算了这个组缺失值个数；大于存在标准为existed; 小于等于存在标准为：non_existed
            na_num <- (ncol(exp_grp)-1) * (1-specific_threshold)
            exp_grp_keep <- exp_grp[exp_grp["data_na"]<na_num,-ncol(exp_grp)]  #为这个组存在的表达量
            exp_grp_del <- exp_grp[exp_grp["data_na"]>=na_num,-ncol(exp_grp)]  #这个组无效表达量；
            exp_grp_keep[exp_grp_keep==0] <- NA
            if (sum(is.na(exp_grp_keep)) == 0){
                method = 'zero'
            }else{
                method = method
            }
            exp_grp_keep_fill <- nafunctions(exp_grp_keep,method= method)
            exp_grp_fill_nrmse <- calculate_nrmsedfout(exp_grp_keep,method=method)
            nrmses <- append(nrmses,exp_grp_fill_nrmse)
            exp_grp_keep_fill_0 <- nafunctions(exp_grp_keep,method= "zero")
            exp_grp_fill_nrmse_0 <- calculate_nrmsedfout(exp_grp_keep,method="zero")  ##创建一个基础比较对象；
            nrmses_basic <- append(nrmses_basic,exp_grp_fill_nrmse_0)
            exp_grp_new <- rbind(exp_grp_del,exp_grp_keep_fill)
            exp_grp_new[["Feature"]] <- rownames(exp_grp_new)
            rownames(exp_grp_new ) <- NULL
            exp_grp_new <- exp_grp_new[,c(ncol(exp_grp_new ),seq(1,ncol(exp_grp_new)-1,1))]
            exp_grp_new[is.na(exp_grp_new)] = 0
            grp_exp_name = paste0("exp_grp_",i)
            #print(grp_exp_name)
            df_list[[grp_exp_name]]<- exp_grp_new}
        result_df <- Reduce(Merge_func,df_list)
        result_nrmse <- mean(nrmses)
        result_nrmse_0 <- mean(nrmses_basic)
        imputeresult <- list(imputedatadf=result_df,nrmse=result_nrmse,nrmse_base= result_nrmse_0, na_df= na_result)
    return(imputeresult)
    }else if(group_specific =="no" && fill_type=='group') {
        nrmses = vector()
        nrmses_basic = vector()
        for(i in unique(group[,2])){
            samples <- group[group[,2] == i,1]
            exp_grp <- data[,names(data) %in% samples]
            exp_grp_keep <- exp_grp
            exp_grp_keep[exp_grp_keep==0] <- NA
            if (sum(is.na(exp_grp_keep)) == 0){
                method = 'zero'
            }else{
                method = method
            }
            exp_grp_keep_fill <- nafunctions(exp_grp_keep,method= method)
            exp_grp_fill_nrmse <- calculate_nrmsedfout(exp_grp_keep,method=method)
            nrmses <- append(nrmses,exp_grp_fill_nrmse)
            exp_grp_keep_fill_0 <- nafunctions(exp_grp_keep,method= "zero")
            exp_grp_fill_nrmse_0 <- calculate_nrmsedfout(exp_grp_keep,method="zero")  ##创建一个基础比较对象；
            nrmses_basic <- append(nrmses_basic,exp_grp_fill_nrmse_0)

            exp_grp_keep_fill[["Feature"]] <- rownames(exp_grp_keep_fill)
            rownames(exp_grp_keep_fill) <- NULL
            exp_grp_keep_fill <- exp_grp_keep_fill[,c(ncol(exp_grp_keep_fill),seq(1,ncol(exp_grp_keep_fill)-1,1))]
            exp_grp_keep_fill[is.na(exp_grp_keep_fill)] = 0
            grp_exp_name = paste0("exp_grp_",i)
            #print(grp_exp_name)
            df_list[[grp_exp_name]]<- exp_grp_keep_fill}
        result_df <- Reduce(Merge_func,df_list)
        result_nrmse <- mean(nrmses)
        result_nrmse_0 <- mean(nrmses_basic)
        imputeresult <- list(imputedatadf=result_df,nrmse=result_nrmse,nrmse_base= result_nrmse_0, na_df= na_result)
    return(imputeresult)
    }else{     #不关心组特异的蛋白，直接将整个数据矩阵填充为完整的矩阵；
        samples <- names(data)
        exp_all <- data[,names(data) %in% samples]
        exp_all_keep <- exp_all
        exp_all_keep[exp_all_keep==0] <- NA
        if (sum(is.na(exp_all_keep)) == 0){
            method = 'zero'
        }else{
            method = method
        }
        exp_all_keep_fill <- nafunctions(exp_all_keep,method = method)
        exp_all_fill_nrmse <- calculate_nrmsedfout(exp_all_keep,method=method)
        exp_all_keep_fill_0 <- nafunctions(exp_all_keep,method= "zero")
        exp_all_fill_nrmse_0 <- calculate_nrmsedfout(exp_all_keep,method="zero")  ##创建一个基础比较对象；

        exp_all_keep_fill[["Feature"]] <- rownames(exp_all_keep_fill)
        rownames(exp_all_keep_fill) <- NULL
        exp_all_keep_fill <- exp_all_keep_fill[,c(ncol(exp_all_keep_fill),seq(1,ncol(exp_all_keep_fill)-1,1))]
        exp_all_keep_fill[is.na(exp_all_keep_fill)] = 0
        result_df <- exp_all_keep_fill
        imputeresult <- list(imputedatadf=result_df,nrmse=exp_all_fill_nrmse,nrmse_base=exp_all_fill_nrmse_0,
                             na_df= na_result)
    return(imputeresult)
    }
}

run_cv_stats<-function(data, group, type=''){
    cv_result<-calculate_cvdata(data, group)
    cv_result[["Feature"]] <- rownames(cv_result)
    rownames(cv_result)<-NULL
    cv_summary<-summarise_cv(cv_result)
    result_cv_add<-merge(data_des,cv_result,by = "Feature", all.y = T)
    result_cv_add <- result_cv_add[,-1]
    write.table(result_cv_add,file= paste0("${outprefix}_cv_stats", type, ".xls"),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(cv_summary,file=paste0("${outprefix}_cv_summary", type, ".xls"),col.names=T,row.names=T,quote=F,sep="\t")
}

if (interact != 'yes'){
    run_cv_stats(data, group, type='_raw')
}
result_data <- my_impute(data,group,na_ratio,method,specific=specific,specific_threshold=specific_threshold, fill_type=fill_type)
result_data_add <- merge(data_des,result_data$imputedatadf,by = "Feature", all.y = T)
result_data_add <- result_data_add[,-1]
processed_data <- result_data$imputedatadf[,-1]
rownames(processed_data) <- result_data$imputedatadf[,1]
run_cv_stats(processed_data, group)
# result_cv_add<-merge(data_des,result_data$cv_df,by = "Feature", all.y = T)
# result_cv_add <- result_cv_add[,-1]
write.table(result_data_add,file="${outprefix}_fillna.txt",col.names=T,row.names=F,quote=F,sep="\t")
# write.table(result_cv_add,file="${outprefix}_cv_stats.xls",col.names=T,row.names=F,quote=F,sep="\t")
# write.table(result_data$cv_summary,file="${outprefix}_cv_summary.xls",col.names=T,row.names=T,quote=F,sep="\t")
nrmse_result <- data.frame(fill= result_data$nrmse,zero = result_data$nrmse_base)
names(nrmse_result) <- c(method,"zero")
write.table(nrmse_result,file="${outprefix}_nrmse.txt",col.names=T,row.names=F,quote=F,sep="\t")

if (interact != 'yes'){
    result_na_add<-merge(data_des,result_data$na_df,by = "Feature", all.y = T)
    result_na_add <- result_na_add[,-1]
    write.table(result_na_add,file="${outprefix}_na_stats.xls",col.names=T,row.names=F,quote=F,sep="\t")
}


