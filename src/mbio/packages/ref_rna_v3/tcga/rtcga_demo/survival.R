library(dplyr)
library(RTCGA.mutations)
library(RTCGA.clinical)

mutationsTCGA(BRCA.mutations, OV.mutations) %>%
	filter(Hugo_Symbol == "TP53") %>%
	filter(substr(bcr_patient_barcode, 14, 15) == "01") %>%
	mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> BRCA_OV.mutations

survivalTCGA(BRCA.clinical, OV.clinical, extract.cols = "admin.disease_code") %>%
	rename(disease = admin.disease_code) -> BRCA_OV.clinical

BRCA_OV.clinical %>%
	left_join(BRCA_OV.mutations, by = "bcr_patient_barcode") %>%
	mutate(TP53 = ifelse(!is.na(Variant_Classification), "Mut", "WILDorNOINFO")) -> BRCA_OV.clinical_mutations

BRCA_OV.clinical_mutations %>%
	select(times, patient.vital_status, disease, TP53) -> BRCA_OV.2plot

kmTCGA(BRCA_OV.2plot, explanatory.names = c("TP53", "disease"), break.time.by = 400, xlim = c(0,4000))
kmTCGA(BRCA_OV.2plot, explanatory.names = "disease", break.time.by = 400, xlim = c(0, 4000))
