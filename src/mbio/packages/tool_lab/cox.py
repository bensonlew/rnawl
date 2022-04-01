import os
from biocluster.config import Config
from multiprocessing import Pool
import subprocess
import pandas as pd


# def run_script(codes):
#     subprocess.call("Rscript {}".format(codes), shell=True)

class cox(object):
    def __init__(self, surv_table, exp_table, factor_list=None, gene_list=None, output=None):
        self.surv_table = surv_table
        self.exp_table = exp_table
        self.factor_list = factor_list.strip().split(';')
        self.output = output
        self.gene_list = gene_list
        if self.gene_list:
            gene = pd.read_table(self.gene_list, header=0, sep='\t')
            self.genes = list(gene['seq_id'])
        software_dir = Config().SOFTWARE_DIR
        self.rscript = software_dir + "/program/R-3.3.1/bin/Rscript"
        # self.rscript = '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/miniconda/bin/Rscript'
    # def run(self, script_list):
    #     pool = Pool(self.pool_size)
    #     pool.map(run_script, script_list)
    #     pool.close()
    #     pool.join()

    def prepare(self):

        meta = pd.read_table(self.surv_table, header=0, index_col=0, sep='\t')
        if self.gene_list:
            count = pd.read_table(self.exp_table, header=0, index_col=0, sep='\t')
            gene = pd.read_table(self.gene_list, header=0, sep='\t')
            median_total = count.loc[gene['seq_id']].median(axis=1)
            for i in gene['seq_id']:
                meta[i] = ['high' if j >= median_total[i] else 'low' for j in count.loc[i]]
            meta.to_csv(os.path.join(self.output, 'meta_new.txt'), sep='\t', header=True, index=True)

    def cox(self, output=None):
        if output is None:
            output = os.getcwd()
        if self.gene_list:
            meta_new = os.path.join(self.output, 'meta_new.txt')
        else:
            meta_new = self.surv_table
        script_name = os.path.join(output, 'cox.r')
        if self.factor_list and self.gene_list:
            f = open(script_name, 'w')
            f.write('library(survival)\n')
            f.write('library(survminer)\n')
            f.write('library(dplyr)\n')
            f.write('library(caret)\n')
            f.write('library(Hmisc)\n')

            f.write('surv_cox <- read.table("{}", sep="\t", header=TRUE)\n'.format(meta_new))
            f.write('expr <- read.table("{}", sep="\t", header=TRUE, row.names=1)\n'.format(self.exp_table))
            f.write('gene_list <- read.table("{}", sep="\t", header=TRUE)\n'.format(self.gene_list))
            factor_all = ["'{}'".format(factor) for factor in self.factor_list+self.genes]
            f.write('covariates <- c({})\n'.format(','.join(factor_all)))
            f.write('univ_formulas <- sapply(covariates, function(x) as.formula(paste("Surv(time, status) ~", x)))\n')
            f.write('univ_models <- lapply(univ_formulas, function(x) {coxph(x, data = surv_cox)})\n')
            f.write('univ_results <- lapply(univ_models, function(x){\n')
            f.write('x <- summary(x)\n')
            f.write('p.value <- signif(x$wald["pvalue"], digits=2)\n')
            f.write('wald.test <- signif(x$wald["test"], digits=2)\n')
            f.write('beta <- signif(x$coef[1], digits=2)\n')
            f.write('HR <- signif(x$coef[2], digits=2)\n')
            f.write('HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)\n')
            f.write('HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)\n')
            f.write('HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")\n')
            f.write('res<-c(beta, HR, wald.test, p.value)\n')
            f.write('names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", "p.value")\n')
            f.write('return(res)\n')
            f.write('})\n')
            f.write('res <- t(as.data.frame(univ_results, check.names = FALSE))\n')
            f.write('res <- as.data.frame(res)\n')
            f.write('write.table(t(c("feature", "beta", "HR (95% CI for HR)", "wald_test", "pvalue")), file="{}", sep="\t", quote=F, col.names=F,row.names=F)\n'.format(self.output+'/cox.txt'))
            f.write('write.table(file="{}", res, sep="\t", quote=F, row.names=T, col.names=F, append=T)\n'.format(self.output+'/cox.txt'))

            f.write('meta<-read.table("{}", sep="\t", head=T, check.names=F)\n'.format(self.surv_table))
            f.write('expr<-read.table("{}", sep="\t", head=T, row.names=1, check.names=F)\n'.format(self.exp_table))
            f.write('gene=read.table("{}", head=T)\n'.format(self.gene_list))
            f.write('meta=meta[meta$seq_id %in% as.factor(colnames(expr)),]\n')
            f.write('expr=expr[,colnames(expr) %in% as.character(meta$seq_id)]\n')
            f.write('coxphs<-function(g){\n')
            f.write('meta$group=ifelse( as.numeric(expr[as.character(g),]) > median(as.numeric(expr[as.character(g),])) ,"high" ,"low" )\n') # modified by zhangyitong on 20211011
            f.write('formula1=as.formula(paste("Surv(time, status) ~ {}","group",sep="+"))\n'.format('+'.join(self.factor_list)))
            f.write('m=coxph(formula1, data =  meta)\n')
            f.write('beta <- coef(m)\n')
            f.write('se <- sqrt(diag(vcov(m)))\n')
            f.write('HR <- exp(beta)\n')
            f.write('HRse <- HR * se\n')
            f.write('tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1), HR = HR, HRse = HRse, '
                    'HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1), HRCILL = exp(beta - qnorm(.975, 0, 1) * se), '
                    'HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)\n')
            f.write('return(tmp["grouplow",])\n')
            f.write('}\n')
            f.write('cox_results <-lapply(gene[,1],coxphs)\n')
            f.write('cox<-t(as.data.frame(cox_results))\n')
            f.write('rownames(cox)=gene[,1]\n')
            f.write('cox=as.data.frame(cox)\n')
            f.write('coxe=rownames(cox[cox$p<1,])\n')
            f.write('coxm=as.data.frame(log2(t(expr[coxe,])+1))\n')
            f.write('coxm$seq_id=rownames(coxm)\n')
            f.write('dat=left_join(meta,coxm,by="seq_id")\n')
            f.write('formula=as.formula(paste("Surv(time, status) ~ {}",paste(coxe,sep="+",collapse="+"),sep="+"))\n'.format('+'.join(self.factor_list)))
            f.write('model <- coxph(formula, data = dat )\n')
            f.write('p <- ggforest(model,data=dat,main = "Hazard ratio", cpositions = c(0.01, 0.22, 0.4), fontsize = 1, refLabel = "reference", noDigits = 3)\n')
            # f.write('ggsave(file="{}",print(p),device="png")\n'.format(self.output + '/cox.png'))
            f.write('ggsave(file="{}",p,device="pdf")\n'.format(self.output + '/cox.pdf'))
        if self.factor_list and not self.gene_list:
            f = open(script_name, 'w')
            f.write('library(survival)\n')
            f.write('library(survminer)\n')
            f.write('library(dplyr)\n')
            f.write('library(caret)\n')
            f.write('library(Hmisc)\n')

            f.write('surv_cox <- read.table("{}", sep="\t", header=TRUE)\n'.format(meta_new))
            factor_all = ["'{}'".format(factor) for factor in self.factor_list]
            f.write('covariates <- c({})\n'.format(','.join(factor_all)))
            f.write('univ_formulas <- sapply(covariates, function(x) as.formula(paste("Surv(time, status) ~", x)))\n')
            f.write('univ_models <- lapply(univ_formulas, function(x) {coxph(x, data = surv_cox)})\n')
            f.write('univ_results <- lapply(univ_models, function(x){\n')
            f.write('x <- summary(x)\n')
            f.write('p.value <- signif(x$wald["pvalue"], digits=2)\n')
            f.write('wald.test <- signif(x$wald["test"], digits=2)\n')
            f.write('beta <- signif(x$coef[1], digits=2)\n')
            f.write('HR <- signif(x$coef[2], digits=2)\n')
            f.write('HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)\n')
            f.write('HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)\n')
            f.write('HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")\n')
            f.write('res<-c(beta, HR, wald.test, p.value)\n')
            f.write('names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", "p.value")\n')
            f.write('return(res)\n')
            f.write('})\n')
            f.write('res <- t(as.data.frame(univ_results, check.names = FALSE))\n')
            f.write('res <- as.data.frame(res)\n')
            f.write(
                'write.table(t(c("feature", "beta", "HR (95% CI for HR)", "wald_test", "pvalue")), file="{}", sep="\t", quote=F, col.names=F,row.names=F)\n'.format(
                    self.output + '/cox.txt'))
            f.write('write.table(file="{}", res, sep="\t", quote=F, row.names=T, col.names=F, append=T)\n'.format(
                self.output + '/cox.txt'))

            f.write('meta<-read.table("{}", sep="\t", head=T, check.names=F)\n'.format(self.surv_table))
            f.write(
                'formula=as.formula("Surv(time, status) ~ {}")\n'.format(
                    '+'.join(self.factor_list)))
            f.write('model <- coxph(formula, data = meta )\n')
            f.write(
                'p <- ggforest(model,data=meta,main = "Hazard ratio", cpositions = c(0.01, 0.22, 0.4), fontsize = 1, refLabel = "reference", noDigits = 3)\n')
            # f.write('ggsave(file="{}",print(p),device="png")\n'.format(self.output + '/cox.png'))
            f.write('ggsave(file="{}",p,device="pdf")\n'.format(self.output + '/cox.pdf'))
        if not self.factor_list and self.gene_list:
            pass
        f.close()
        os.system('{} cox.r'.format(self.rscript))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to identify pirna from piRbase database.')
    parser.add_argument('-m', type=str, help='meta file')
    parser.add_argument('-e', type=str, help='exp file')
    parser.add_argument('-f', type=str, help='factor list')
    parser.add_argument('-g', type=str, help='gene list', default=None)
    parser.add_argument('-o', type=str, help='output path')
    args = parser.parse_args()
    cox = cox(args.m, args.e, args.f, args.g, args.o)
    cox.prepare()
    cox.cox()
