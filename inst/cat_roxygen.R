

cat_roxygen <- function(model, funName, fileName){
  
  cat1 <- function(...){cat(..., file = fileName, append = T)}
  #cat1 <- function(...){cat(...)}
  
  cat1(paste0("#' ","Fit the ", model$name," model", "\n"))
  cat1("\n")
  cat1(paste0("#' @description ","Fit the ", model$name," model to SAR data", "\n"))
  cat1(paste0("#' @usage ", funName,"(data, custstart = NULL, normtest = 'lillie')", "\n"))
  cat1(paste0("#' @param ", "data ", "A dataset in the form of a dataframe with two columns: 
                    the first with island/site areas, and the second with the species richness 
                    of each island/site.", "\n"))
  cat1(paste0("#' @return "))

  cat1(paste0("#' @examples", "\n", "#' data(galap)", "\n", "#' fit <- ", funName,
              "(galap)", "\n", "#' summary(fit)","\n", "#' plot(fit)", "\n"))

  cat1(paste0("#' @export"))
  
}#eo cat_roxygen





