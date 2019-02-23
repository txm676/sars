## Try to get a function that build that creates sars.

sars_builder <- function(data, name, normaTest = "lillie",
  homoTest = "cor.fitted", start = NULL, grid_start = NULL, verb = TRUE) {

  if (name == "linear") return(sar_linear(data, normaTest, homoTest))

  data <- check_data(data)
  model <- switch(name,
    asymp = model_asymp(data),
    betap = model_betap(data),
    chapman = model_chapman(data),
    epm1 = model_epm1(data),
    epm2 = model_epm2(data),
    gompertz = model_gompertz(data),
    heleg = model_heleg(data),
    koba = model_koba(data),
    loga = model_loga(data),
    mmf = model_mmf(data),
    monod = model_monod(data),
    negexpo = model_negexpo(data),
    p1 = model_p1(data),
    p2 = model_p2(data),
    power = model_power(data),
    powerR = model_powerR(data),
    ratio = model_ratio(data),
    weibull3 = model_weibull3(data),
    weibull4 = model_weibull4(data),
    message_sars()
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

message_sars <- function() {
  stop(
    paste0("Valid values for argument `name` are: ",
    paste(sars_models(), collapse = ", ")
    )
  )
}
