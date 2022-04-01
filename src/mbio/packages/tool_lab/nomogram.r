# load package
library(survival)
library(rms)

args = commandArgs(T)
surv_path = args[1]
surv_path
print(args)
pdf_path = args[4]

surv_cox <- read.table(surv_path, sep='\t', header=TRUE)
dd = datadist(surv_cox)
options(datadist='dd')
res.cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data =  lung)



km <- survfit(Surv(time,status)~group,data=kmsurvival)
p <- ggsurvplot(km,
           pval = TRUE,
           conf.int = args[2],
           #pval.method = TRUE,
           risk.table = args[3], # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF")
)
ggsave(pdf_path, p)








