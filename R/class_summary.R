#' Summarising the results of the five main NOSM functions
#'
#' @description S3 method for class 'NOSM'. summary.NOSM creates summary
#'   statistics for objects of class NOSM. The exact summary statistics computed
#'   depends on the 'Type' attribute (e.g. 'bip') of the NOSM object (see
#'   below). The summary method generates more useful information for the user
#'   than the standard NOSM functions. Another S3 method (print.summary.NOSM;
#'   not documented) is used to print the output.
#' @param object An object of class 'NOSM'.
#' @param y (default of 3) The adjustment value for the computation of the z
#'   value (see Strona & Veech, 2015).
#' @param ...	further arguments passed to or from other methods.
#' @return Returns object of class 'summary.NOSM' with a Type attribute (e.g.
#'   'bip') which is inherited. For NOSM objects of Type 'Pot_dir', 'bip' or
#'   'Dir', the summary.NOSM method returns the mean of the overlap values for
#'   the "in nodes" (NOS_In), the mean of the overlap values for the "out nodes"
#'   (NOS_Out), the mean of Nos In and Nos Out (NOS), the standard
#'   deviation of the overlap values for the "in nodes" (MOD_In), the SD of the
#'   overlap values for the "out nodes" (MOD_Out), and the SD of the combined
#'   set of overlap values (MOD; network modularity).
#'
#'   For NOSM objects of Type 'Dir' and 'Undir', the summary.NOSM method returns
#'   just the NOS and MOD values (network modularity).
#'
#'   For all types of NOSM object, the z value and associated p value are also
#'   provided (see Strona & Veech, 2015).
#' @seealso \code{\link{NOSM_bip}}, \code{\link{NOSM_POT_dir}},
#'   \code{\link{NOSM_POT_undir}}, \code{\link{NOSM_dir}},
#'   \code{\link{NOSM_undir}}
#' @references Strona, G. & Veech, J. A. (2015). A new measure of ecological
#'   network structure based on node overlap and segregation. Methods in Ecology
#'   and Evolution, 6(8), 907-915.
#' @examples
#' data(boreal)
#' z <- boreal[sample(rownames(boreal), 200, FALSE),] #subset for speed
#' x <- NOSM_bip(z, perc = 1, sl = 1)
#' summary(x, y = 3)
#' @importFrom dplyr %>%
#' @export


summary.mmSAR2 <- function(object, ..., y = 3){

  if (attributes(object)$Type == "lin_pow"){
    object2 <- object$Model
    c = object2$coefficients[1, 1]
    z = object2$coefficients[2, 1]
    z.sig = object2$coefficients[2, 4]
    r2 = object2$r.squared
    md_res <- c(c, z, z.sig, r2) %>% round(2)
    names(md_res) <- c("c", "z", "z.sig", "r2")
    fit_df <- data.frame(Area = object$Area, Fitted = object$Fitted) %>% round(2)
    res <- list(Summary = md_res, df = fit_df)
  }
  class(res) <- "summary.mmSAR2"
  attr(res, "Type") <- attributes(object)$Type
  attr(res, "Dataset") <- attributes(object)$Dataset
  return(res)
}

