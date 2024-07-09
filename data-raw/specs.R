## Prepare `specs` dataset goes here
beam.name=paste0(rep(c("Lesion","RNA"),times=3),
                 rep(c(".MRD29",".EFS",".OS"),each=2))
beam.mtx=rep(c("Lesion","RNA"),times=3)
beam.mdl=rep(c("lm(MRD29~mtx.row,data=main.data,model=T)",
               "coxphf2(EFS~mtx.row,data=main.data)",
               "coxphf2(OS~mtx.row,data=main.data)"),each=2)

specs=cbind.data.frame(name=beam.name,
                       mtx=beam.mtx,
                       mdl=beam.mdl)

usethis::use_data(specs, overwrite = TRUE)
