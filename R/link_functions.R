#### INTERNAL FUNCTION(S)


#  LINK FUNCTIONS FOR CONSTRAINT OPTIMIZATION (internal)

#We use transLink for Rplus and unif parameters to provide starting parameter
#values that then constrain the parameter values optim generates to be within a
#more sensible range (i.e. mimicking boundary conditions). Otherwise the scale
#of the parameter is too large and optim generates unsuitable parameter
#estimates to test that can be negative etc. The transformed starting estimates
#are back-transformed within the model rss function anyway so it is the proper
#parameter estimates that are used to generate the rss.

transLink <- function(x,boundType) {

  # fonction logit
  logit <- function(y) {
    log(y/(1-y))
  }#end of logit

  if (length(x) == 1) {
    res <- switch(boundType,
                 R = x,
                 Rplus = log(x),
                 unif = logit(x)
    )#end of switch
  } else {

    res <- vector(length = length(x))

    for (i in seq_along(x)) {

      res[i] <- switch(boundType[i],
                      R = x[i],
                      Rplus = log(x[i]),
                      unif = logit(x[i])
      )#end of switch
    }#end of for
    res
  }#end of if/else
}#end of transLink

##

backLink <- function(x,boundType) {

  # fonction reciproque de la fonction logit
  invlogit <- function(x) {
    1/(1+exp(-x))
  }#end of invlogit

  if (length(x) ==1) {
    res <- switch(boundType,
                 R = x,
                 Rplus = exp(x),
                 unif = invlogit(x)
    )#end of switch
  } else {

    res <- vector(length = length(x))

    for (i in seq_along(x)) {

      res[i] <- switch(boundType[i],
                      R = x[i],
                      Rplus = exp(x[i]),
                      unif = invlogit(x[i])
      )#end of switch
    }#end of for
    res
  }#end of if/else
}#end of backLink
