#' Print summary information about beam.stats object
#'
#' @param x An object of class "beam.stats"
#' @param ... Other arguments passed to or from other methods
#'
#' @return Messages about the beam.data object
#' @export
#'
#' @examples
#' \dontrun{
#' data(clinf)
#' data(omicdat)
#' data(omicann)
#' data(setdat)
#' test.beam.data <- prep_beam_data(main.data=clinf, mtx.data=omicdat,
#'                                  mtx.anns=omicann, set.data=setdat,
#'                                  set.anns=NULL, n.boot=10, seed=123)
#' specs <- prep_beam_specs(beam.data=test.beam.data, endpts=c("MRD29", "EFS", "OS"))
#' test.beam.stats <- compute_beam_stats(beam.data=test.beam.data, beam.specs=specs)
#' print(test.beam.stats)
#' }
#' data(beam_stats)
#' print(beam_stats)
print.beam.stats <- function(x,...)
{
  beam.stats <- unclass(x)
  beam.stats <- as.list(beam.stats)
  cat(paste0("Contains ",
             length(beam.stats$beam.stats)," association estimate matrices with ",
             nrow(beam.stats$beam.data$boot.index)-1," bootstraps: \n"))

  splt <- strsplit(names(beam.stats$beam.stats), split="\\.")
  splt.df <- do.call(rbind.data.frame,splt)
  colnames(splt.df) <- c("Omics", "Endpoints")

  for(i in 1:length(splt)){
    cat(paste0("  Association Estimate Matrix of ", splt[[i]][1], " with ", splt[[i]][2],
               " has dimensions ", nrow(beam.stats$beam.stats[[i]]), " x ",
               ncol(beam.stats$beam.stats[[i]]), ". \n"))
  }

  for(i in 1:length(splt)){
    cat(paste0("\nExample Association Estimate Matrix for ", splt[[i]][1], " with ",
               splt[[i]][2],": \n"))
    print(beam.stats$beam.stats[[i]][1:5,1:5])
  }

  cat(paste0("\nExample Endpoint Data: \n"))
  un.eps <- unique(splt.df$Endpoints)
  print(beam.stats$beam.data$main.data[1:5,which(colnames(beam.stats$beam.data$main.data) %in% un.eps)])

  cat(paste0("\nBEAM Model Specifications: \n"))
  print(beam.stats$beam.specs)

  cat(paste0("\nBEAM data used to create Association Estimate Matrices: \n"))
  print(beam.stats$beam.data)


}
