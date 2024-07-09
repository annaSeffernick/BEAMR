###################################
# Create Example Data for
# BEAMR R Package
# AES
# 2023-03-16
######################################

ddir <- "C:/Users/aseffern/Box/BookChapterApril2022"
#source(paste0(ddir, "/Xueyuan-BEAMR/prep-beam-data.R"))
library(survival)

lesion<-readRDS(paste0(ddir, "/Abdel-GRIN/TALL-binary-lsn-mtx-atleast5patients.RDS"))
#lesion <- cbind(lesion.ID=rownames(lesion), lesion)
lesion.id <- order(rowSums(lesion), decreasing = T)
lesion.sm <- lesion[lesion.id[1:20],]
# Extract ensembl IDs to get matching RNA if it exists
les.rw <- rownames(lesion.sm)
les.rw.list <- strsplit(les.rw, "_")
les.rw.df <- do.call(rbind.data.frame, les.rw.list)
colnames(les.rw.df) <- c("ID", "Type")

rdat<-"C:/Users/aseffern/Box/BookChapterApril2022/DataSet/v1/TALL-dataset.Rdata"
load(rdat)
geneann<-read.csv(paste0(ddir, "/DataSet/v1/ensembl.annotation.csv"))
RNA<-RNA[!is.na(RNA$ensembl.ID), ]
rownames(RNA)<-RNA$ensembl.ID
match.id <- which(rownames(RNA) %in% les.rw.df$ID)
RNA.sm <- RNA[c(1:14, match.id),]
geneann.sm <- geneann[which(geneann$ensembl.ID %in% rownames(RNA.sm)),]
#RNA<-RNA[rownames(RNA)%in%substring(rownames(lesion), 1, 15), ]
#lesion<-lesion[substring(rownames(lesion), 1, 15)%in%rownames(RNA), ]

omicdat<-list(Lesion=as.matrix(lesion.sm),
              RNA=as.matrix(RNA.sm[, -1]))
omicann<-list(Lesion=cbind.data.frame(prob.id=rownames(lesion.sm),
                                      gene.id=substring(rownames(lesion.sm), 1, 15)),
              RNA=geneann.sm)


clin$EFScensor[clin$First_Event%in%c("None", "Censored")]<-0
clin$EFScensor[clin$First_Event%in%c("Death","Progression", "Relapse" , "Second Malignant Neoplasm")]<-1

clin$OScensor[clin$Vital_Status%in%"Alive"]<-0
clin$OScensor[clin$Vital_Status%in%"Dead"]<-1

rownames(clin)<-clin$ID
clin$"RNA.clm" <- clin$ID
clin$"Lesion.clm" <- clin$ID
clin$RNA.id<-clin$Lesion.id<-clin$ID
clin$EFS<-Surv(clin$Event_Days, clin$EFScensor)
clin$OS<-Surv(clin$OS_Days, clin$OScensor)

clinf<-clin[, c(1, 8, 15:ncol(clin))]

omicann<-list(Lesion=cbind.data.frame(id=rownames(lesion.sm),
                                      gene=substring(rownames(lesion.sm), 1, 15)),
              RNA=cbind.data.frame(id=rownames(RNA.sm),
                                   gene=rownames(RNA.sm)))
setdat<-cbind.data.frame(set.id=c(rownames(RNA.sm), substring(rownames(lesion.sm), 1, 15)),
                         mtx.id=c(rep("RNA", nrow(RNA.sm)), rep("Lesion", nrow(lesion.sm))),
                         row.id=c(rownames(RNA.sm), rownames(lesion.sm)))


usethis::use_data(omicdat, overwrite = TRUE)
usethis::use_data(clinf, overwrite = TRUE)
usethis::use_data(omicann, overwrite = TRUE)
usethis::use_data(setdat, overwrite = TRUE)

# Keep intermediate steps to save time
beam_dat <- prep_beam_data(main.data=clinf, mtx.data=omicdat,
                           mtx.anns=omicann, set.data=setdat,
                           set.anns=NULL, n.boot=10, seed=123)
beam_specs <- prep_beam_specs(beam.data=beam_dat, endpts=c("MRD29", "EFS", "OS"),
                              firth=TRUE)
beam_stats <- compute_beam_stats(beam.data=beam_dat, beam.specs=beam_specs)

usethis::use_data(beam_dat, overwrite = TRUE)
usethis::use_data(beam_specs, overwrite = TRUE)
usethis::use_data(beam_stats, overwrite = TRUE)

# make very small dataset to compute_beam_stats with in less than 5s
beam_dat_sm <- prep_beam_data(main.data=clinf, mtx.data=omicdat,
                           mtx.anns=omicann, set.data=setdat,
                           set.anns=NULL, n.boot=2, seed=123)
beam_specs_sm <- prep_beam_specs(beam.data=beam_dat, endpts=c("MRD29", "EFS"),
                              firth=TRUE)
beam_stats_sm <- compute_beam_stats(beam.data=beam_dat_sm, beam.specs=beam_specs_sm)

usethis::use_data(beam_dat_sm, overwrite = TRUE)
usethis::use_data(beam_specs_sm, overwrite = TRUE)
usethis::use_data(beam_stats_sm, overwrite = TRUE)
