library(broom)
library(survival)
library(methods)
library(getopt)

command <- matrix(c(
  'clinical', 'c', 1, 'character', 'User-input file of clinical data, columns of survival time and survival status must be included.',
  'exp', 'e', 1, 'character', 'User-input file of expression matrix.',
  'geneset', 'g', 2, 'character', 'A string of gene symbols, separated by comma. The default is an empty string.',
  'method', 'm', 1, 'character', 'One of the available selection methods: lasso, stepwise, backward, forward',
  'help', 'h', 0, 'logical', 'Show this help message and exit'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

#check options
if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}

#data processing, merge clinical and exp files
processing<-function(clinical,exp,geneset=''){
  clinical_df<-read.table(clinical,header=T)
  exp_df<-read.table(exp,row.names = 1,header=T)

  if('sample'%in%colnames(clinical_df)){
    rownames(clinical_df)<-clinical_df$sample
  }else{
    rownames(clinical_df)<-clinical_df[[1]]
  }
  os<-subset(clinical_df,select =c(time,status))
  os<-na.omit(os)
  os<-os[os$time!=0,]
  
  sample<-intersect(colnames(exp_df),rownames(os))
  exp_filtered<-exp_df[,colnames(exp_df)%in%sample]
  clinical_filtered<-os[rownames(os)%in%sample,]
  
  clinical_filtered$sample<-rownames(clinical_filtered)
  exp_filtered_t<-data.frame(t(exp_filtered))
  if (geneset!=''){
    geneset<-strsplit(geneset,',')[[1]]
    exp_filtered_t<-exp_filtered_t[,colnames(exp_filtered_t) %in% geneset]
    if(!is.null(!(geneset %in% colnames(exp_filtered_t)))){
        print('The following genes were not found in the expression profile:')
        print(geneset[!(geneset %in% colnames(exp_filtered_t))])
    }
    geneset<-geneset[geneset %in% colnames(exp_filtered_t)]
  }else{
    geneset<-colnames(exp_filtered_t)
  }
  exp_filtered_t$sample<-row.names(exp_filtered_t)
  merged_df<-merge(exp_filtered_t,clinical_filtered,by="sample")
  
  merged_df$s_status<-''
  status<-sort(unique(merged_df$status))
  status_1<-as.integer(c(0,1))
  status_2<-as.integer(c(1,2))
  if (identical(status,status_1)==FALSE & identical(status,status_2)==FALSE){
    deceased<-c('DEAD','DECEASED')
    alive<-c('LIVING','ALIVE')
    merged_df[toupper(merged_df$status) %in% alive,]$s_status<-1
    merged_df[toupper(merged_df$status) %in% deceased,]$s_status<-2
  }else if (identical(status,status_1)==TRUE){
    merged_df[merged_df$status==0,]$s_status<-1
    merged_df[merged_df$status==1,]$s_status<-2
  }else if (identical(status,status_2)==TRUE){
    merged_df$s_status<-merged_df$status
  }
  merged_df$s_status<-as.numeric(merged_df$s_status)

  processed<-list(geneset=geneset,merged_df=merged_df,exp_filtered=exp_filtered)
  return(processed)
}


# add legend for LASSO plot1
lbs_fun <- function(fit) {
  L <- length(fit$lambda)
  y <- fit$beta[, L]
  labs <- names(y)
  if (length(labs)<25){
    legend('bottomright', legend=labs[-1], col=1:6, lty=1,cex=.5)
  }
}


# LASSO
lasso<-function(geneset,merged_df){
  library(glmnet)
  library(caret)
  
  set.seed(1)
  formula<-as.formula(paste('~ ',paste(geneset, collapse= "+")))
  x<-model.matrix(formula,merged_df)
  y <- Surv(merged_df$time, merged_df$s_status)
  
  fit <- glmnet(x, y, family="cox")
  fit_cv <- cv.glmnet(x, y, family="cox")

  pdf("Lasso_result1.pdf")
  plot(fit,xvar = 'lambda',xlab="log(lambda)")
  lbs_fun(fit)
  abline(v=log(fit_cv$lambda.min), lty=3)
  dev.off()

  pdf("Lasso_result2.pdf")
  plot(fit_cv, xlab="log(lambda)")
  dev.off()
  
  coefficients<-coef(fit, s = fit_cv$lambda.min)
  result<-tidy(coefficients)[,c(1,3)]
  colnames(result)<-c('Gene','Coefficient')
  return(result)
}

#Stepwise, Backward, Forward
step_selection<-function(geneset,merged_df,method){
  set.seed(1)
  formula<-as.formula(paste('Surv(time,s_status)~ ',paste(geneset, collapse= "+")))
  model_cox<-coxph(formula,data=merged_df)
  if (method=='stepwise'){
    tstep<-step(model_cox, direction = "both")
  }else if (method =='forward'){
    tstep<-step(model_cox, direction = "forward")
  }else if (method == 'backward'){
    tstep<-step(model_cox, direction = "backward")
  }

  aic<-drop1(tstep)
  rownames(aic)[1]<-'Overall'
  aic$Gene<-rownames(aic)
  result<-aic[,c(3,2)]
  return(result)
}


# Run
if(is.null(opts$geneset)){
    processed = processing(opts$clinical, opts$exp)
}else{
    processed = processing(opts$clinical, opts$exp, geneset=opts$geneset)
}

if(opts$method == 'lasso'){
  result = lasso(processed$geneset, processed$merged_df)
}else{
  result = step_selection(processed$geneset, processed$merged_df, opts$method)
}
write.table(result, paste0(opts$method,'_filtered_genelist.xls'), col.names = T, sep = '\t', quote = F, row.names = F)

# Generate an exp file for filtered genes
selected_exp<-processed$exp_filtered[rownames(processed$exp_filtered)%in%result$Gene,]
write.table(selected_exp, "genelist_exp.xls", col.names = T, sep = '\t', quote = F, row.names = T)