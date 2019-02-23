## Try to get a function that build that creates sars.

sars_builder <- function(data, name, start = NULL, grid_start = NULL, normaTest = "lillie",
  homoTest = "cor.fitted", verb = TRUE) {
  data <- check_data(data)
  model <- switch(name,
    asymp = model_asymp(data),
    betap = model_betap(data),
    chapman = model_chapman(data),
    epm1 = model_epm1(data)
  )
  model <- compmod(model)
  fit <- get_fit(model = model, data = data, start = start, grid_start = grid_start, algo = 'Nelder-Mead',
    normaTest =  normaTest, homoTest = homoTest, verb = verb)
  #
  return_fit(fit)
}



check_data <- function(data) {
  #
  if (!(is.matrix(data) | is.data.frame(data)))
    stop('data must be a matrix or dataframe')
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop('NAs present in data')
  data <- data[order(data[,1L]), ]
  #
  colnames(data) <- c('A', 'S')
  data
}

return_fit <- function(fit) {
  if (is.na(fit$value)) {
    fit <- list(value = NA)
  } else {
    # return(list(value = NA))
    obs <- obs_shape(fit)
    fit$observed_shape <- obs$fitShape
    fit$asymptote <- obs$asymp
    class(fit) <- 'sars'
    attr(fit, 'type') <- 'fit'
  }
  fit
}
