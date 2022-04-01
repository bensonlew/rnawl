library(RTCGAToolbox)
library(survival)

BRCA.data <- getFirehoseData("BRCA", clinical = TRUE, Mutation = TRUE)

clinical.data <- biocExtract(BRCA.data, "clinical")

survival.data <- data.frame(Samples = rownames(clinical.data),
                            Time = as.numeric(clinical.data[, 4]),
                            Censor = as.numeric(clinical.data[, 3]))

surv.obj <- Surv(survival.data$Time, survival.data$Censor)
survfit.list <- survfit(surv.obj ~ 1)
plot(survfit.list)
