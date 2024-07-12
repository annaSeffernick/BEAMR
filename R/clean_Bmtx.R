#' Clean up bootstrap coefficient matrix
#'
#' @param B Matrix of bootstrap coefficients
#'
#' @return Matrix of cleaned bootstrap coefficients
#' @export
#'
#' @examples
#' \donttest{
#' data(clinf)
#' data(omicdat)
#' data(omicann)
#' data(setdat)
#' test.beam.data <- prep_beam_data(main.data=clinf, mtx.data=omicdat,
#'                                  mtx.anns=omicann, set.data=setdat,
#'                                  set.anns=NULL, n.boot=10, seed=123)
#' specs <- prep_beam_specs(beam.data=test.beam.data,
#'                          endpts=c("MRD29", "EFS", "OS"))
#' test.beam.stats <- compute_beam_stats(beam.data=test.beam.data,
#'                                       beam.specs=specs)
#' B.mtx <- test.beam.stats$beam.stats[[1]]
#' B.cln <- clean_Bmtx(B.mtx)
#' }
#'
#'data(beam_stats)
#'B.mtx <- beam_stats$beam.stats[[1]]
#'B.cln <- clean_Bmtx(B.mtx)
clean_Bmtx=function(B)
{
  m=nrow(B)
  res=apply(B,1,replace_value)
  res=t(res)
}

############################
# replace values

replace_value=function(x)
{
  inf=!is.finite(x)
  na=is.na(x)
  mn=mean(x[na&!inf])
  x[na]=mn
  x[inf]=sign(x[inf])*max(abs(x[!inf]))
  return(x)
}
