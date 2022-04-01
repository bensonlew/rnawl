# load package
library(survival)
library(survminer)
library(getopt)
command = matrix(c(
    'help', 'h',  0, 'logical',
    'km_table', 'k', 1, 'character',
    'risk_table', 'r', 1, 'logical',
    'conf', 'c', 1, 'logical',
    'output', 'o', 1, 'character'
  ),byrow=TRUE,ncol =4);
opt = getopt(command);
pdf_path = opt$output
kmsurvival <- read.table(opt$km_table, sep='\t', header=TRUE)
km <- survfit(Surv(time,status)~group,data=kmsurvival)
pdf.options(reset=TRUE, onefile=FALSE)
pdf(pdf_path)
ggsurvplot(km,
           pval = TRUE,
           conf.int = opt$conf,
           #pval.method = TRUE,
           risk.table = opt$risk_table, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF")
)
dev.off()
#ggsave(pdf_path, plot=print(p), device = 'pdf')








