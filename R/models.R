#### INTERNAL FUNCTION(S)

# Asymptotic Regression
model_asymp <- function(data) {
  list(
    name = "Asymptotic regression",
    formula = expression(S == d - c*z^A),
    exp = expression(d - c*z^A),
    shape = "convex",
    asymp = function(pars)pars["d"],
    parLim  = c("Rplus","R","R"),
    #initial values function
    init = function(data) {#Ratkowsky 1983 p178
      #d determination (asymptote)
      d = max(data$S) + .25*max(data$S)
      #Intermediate variable calculation
      Z = log(d-data$S)
      #we have also Z=log(c)+Xlog(z) -> linear regression
      dat = data.frame("a" = data$A, "Z" = Z)
      zf = stats::lm(Z~a,dat)$coefficients
      c(d, exp(zf[1]), exp(zf[2]))
    }
  )
}

# Beta-P function (cumulative)
model_betap <- function(data) {
  list(
    name = c("Beta-P cumulative"),
    formula = expression(S == d*(1-(1+(A/c)^z)^-f)),
    exp = expression(d*(1-(1+(A/c)^z)^-f)),
    shape = "sigmoid",
    asymp = function(pars) pars["d"],
    parLim  =  c("Rplus", "R", "R", "R"),
    #initial values function
    init = function(data) c(max(data$S)+1, .5, .5, .5)
  )
}

#Chapman–Richards 3 S = a [1 − exp(−bA)]c Flather (1996)
model_chapman <- function(data) {
  list(
    name = c("Chapman Richards"),
    formula = expression(S == d * (1 - exp(-z*A)^c )),
    exp = expression(d * (1 - exp(-z*A)^c )),
    shape = "sigmoid",
    asymp = function(pars)pars["d"],
    #limits for parameters
    parLim  =  c("Rplus","R","R"),
    #initials values function
    init = function(data) {
      d = max(data$S)
      Z = (-log((-data$S/(max(data$S)+1))+1))/data$A
      z = mean(Z)
      c(d,z,1)
    }
  )
}

# EXTENDED POWER MODEL 1 (TJORVE 2009)
model_epm1 <- function(data) {
  list(
    name = c("Extended Power model 1"),
    formula = expression(S==c*A^(z*A^-d)),
    exp = expression(c*A^(z*A^-d)),
    shape = "sigmoid",
    asymp = function(pars) FALSE,
    parLim = c("Rplus","R","R"),
    custStart = function(data) c(5, .25, .15),
    #initials values function
    init = function(data) {
      if (any(data$S==0)) {
        log.data <- data.frame(A=log(data$a), S=log(data$S+.5))
      } else log.data <- log(data)
      res <- stats::lm(S~A,log.data)$coefficients
      res <- c(exp(res[1]),res[2],.15)
      res
    }
  )
}

# EXTENDED POWER MODEL 2 (TJORVE 2009)
model_epm2 <- function(data) {
  list(
    name = c("Extended Power model 2"),
    formula = expression(S==c*A^(z-(d/A))),
    exp = expression(c*A^(z-(d/A))),
    shape = "sigmoid",
    asymp = function(pars) FALSE,
    custStart = function(data) c(5, .25, 1),
    #limits for parameters
    parLim = c("Rplus","unif","R"),
    #initials values function
    init = function(data) c(0,0,0)
  )
}

# Gompertz model
model_gompertz <- function(data) {
  list(
    name = c("Gompertz"),
    formula = expression(S==d*e^(-e^(-z*(A-c)))),
    exp = expression(d*exp(-exp(-z*(A-c)))),
    shape= "sigmoid",
    asymp = function(pars)pars["d"],
    parLim = c("Rplus","R","R"),
    custStart = function(data) {
      if (any(data$S==0)) data=data[data$S!=0,]
      #d determination (asymptote)
      d <- max(data$S)+1
      #t.0 determination (min obs. age)
      t.0 <- min(data$A)-1
      #Intermediate variable calculation
      Z <- log(-log(data$S/d))
      #we have also Z=-kT + kt0 -> linear regression
      dat <- data.frame("a"=data$A,"Z"=Z)
      reg <- stats::lm(Z~a,dat)$coefficients
      #transformations of coeficients
      k.first <- -reg[2L]
      k.second <- reg[1L]/t.0
      k.final <- mean(c(k.first,k.second))
      #estimates return
      c(d,t.0,k.final)
    },
    #initials values function
    init = function(data) {
      if(any(data$S==0)) data <- data[data$S!=0,]
      # d determination (asymptote)
      d <- max(data$S)+1
      # t.0 determination (min obs. age)
      t.0 <- min(data$A)-1
      # Intermediate variable calculation
      Z = log(-log(data$S/d))
      # we have also Z=-kT + kt0 -> linear regression
      dat = data.frame("a" = data$A, "Z" = Z)
      reg = stats::lm(Z~a,dat)$coefficients
      #transformations of coeficients
      k.first <- -reg[2L]
      k.second <- reg[1L]/t.0
      k.final <- mean(c(k.first,k.second))
      #estimates return
      c(d,t.0,k.final)
    }
  )
}

# LOGISTIC FUNCTION (HE & LEGENDRE 1996)
model_heleg <- function(data) {
  list(
    name = c("Heleg(Logistic)"),
    formula = expression(S == c/(f + A^(-z))),
    exp = expression(c/(f + A^(-z))),
    shape = "sigmoid",
    asymp = function(pars)pars["c"]/pars["f"],
    parLim = c("Rplus","Rplus","Rplus"),
    custStart = function(data)c(max(data$S),10,.01),
    #initial values function
    init = function(data){
      if(any(data$S==0)) data <- data[data$S!=0,]
      # c calculation (asymptote)
      c <- max(data$S)+1
      # Intermediate variable calculation
      # long=length(data[[2]])
      Z <- log((c/data$S) -1)
      # We have the z and f init values by linear regression of Z on data[[1]]
      dat <- data.frame("a"=data$A,"Z"=Z)
      zf <- stats::lm(Z~a,dat)$coefficients
      c(max(data$S)*exp(-zf[1]), exp(-zf[1]), -zf[2])
    }
  )
}

# Kobayashi logarithmic (KOBAYASHI 1975), convex upward, no asymptote
model_koba <- function(data) {
  list(
    name = "Kobayashi",
    formula = expression(S == c*log(1+ A/z)),
    exp = expression(c*log(1 + A/z)),
    shape = "convex",
    asymp = function(pars) FALSE,
    #limits for parameters
    parLim = c("R","Rplus"),
    custStart = function(data) c(5,0.1),
    #initials values function
    init = function(data) c(max(data$S), max(data$S))
  )
}



# Logarithmic MODEL (GLEASON 1922)
model_loga <- function(data) {
  list(
    name = "Logarithmic",
    formula = expression(S==c+z*log(A)),
    exp = expression(c+z*log(A)),
    shape = "convex",
    asymp = function(pars) FALSE,
    #limits for parameters
    parLim = c("R","R"),
    custStart = function(data) c(5,.25),
    #initials values function
    init = function (data) {
      semilog.data = data.frame(log(data$A),data$S)
      names(semilog.data) = c("A","S")
      par=stats::lm(S~A, semilog.data)$coefficients
      names(par) = c("c", "z")
      par
    }
  )
}



# "Morgan Mercier Family" curve (Williams et al. 2009 formula)
model_mmf <- function(data) {
  list(
    name = "MMF",
    formula = expression(S==d/(1+c*A^(-z))),
    exp = expression(d/(1+c*A^(-z))),
    shape = "sigmoid",
    asymp = function(pars) pars["d"],
    #limits for parameters
    parLim = c("Rplus","Rplus","Rplus"),
    custStart = function(data) c(max(data$S), 5, .25),
    init = function(data) {
      if (any(data$S==0)) data <- data[data$S!=0,]
      d = max(data$S)*4
      newVar = log((d/data$S) - 1)
      reg = stats::lm(newVar~log(data$A))
      c = exp(reg$coefficients[1L])
      z = -reg$coefficients[2L]
      c(d,c,z)
    }
  )
}


# MONOD CURVE (MONOD 1950, Willimas et al. 2009 formula)
model_monod <- function(data) {
  list(
    name = "Monod",
    formula = expression(S==d/(1+c*A^(-1))),
    exp = expression(d/(1+c*A^(-1))),
    shape = "convex",
    asymp = function(pars)pars["d"],
    #limits for parameters
    parLim = c("Rplus", "Rplus"),
    custStart = function(data) c(stats::quantile(data$A,c(0.25)), max(data$S)),
    #initials values function
    init = function(data) {
      if (any(data$S==0)) data=data[data$S!=0,]
      d <- as.double(max(data$A)+max(data$S)/4)
      c <- data[[1]]*(d/data$S - 1)
      c(d,stats::quantile(c, c(0.25)))
    }
  )
}


# NEGATIVE EXPONENTIAL (Holdridge et al. 1971)
model_negexpo <- function(data) {
  list(
    name = c("Negative exponential"),
    formula = expression(S == d*(1 - exp(-z*A) )),
    exp = expression(d*(1-exp(-z*A))),
    shape = "convex",
    asymp = function(pars) pars["d"],
    #limits for parameters
    parLim = c("Rplus","unif"),
    custStart = function(data) c(max(data$S), .01),
    #initials values function
    init = function(data) {
      d <- max(data$S)
      Z <- (-log( (-data$S/(max(data$S)+1))+1))/data$A
      z <- mean(Z)
      c(d,z)
    }
  )
}


# PERSISTENCE FUNCTION 1 (TJORVE 2009)
model_p1 <- function(data) {
  list(
    name = c("Persistence function 1"),
    formula = expression(S == c*A^z * exp(-d*A)),
    exp = expression(c*A^z * exp(-d*A)),
    mod = s~c*a^z * exp(-d*a),
    shape = "convex",
    asymp = function(pars) FALSE,
    custStart = function(data)c(5,.25,.15),
    #limits for parameters
    parLim = c("Rplus","Rplus","Rplus"),
    #initials values function
    init = function(data) {
      if (any(data$S == 0)) {
        log.data <- data.frame(A = log(data$A), S = log(data$S + .5))
      } else log.data <- log(data)
      res <- stats::lm(S~A,log.data)$coefficients
      res <- c(res[1],res[2],1)
      names(res) <- c("c","z","d")
      res
    }
  )
}


# PERSISTENCE FUNCTION 2 (TJORVE 2009)
model_p2 <- function(data) {
  list(
    name = c("Persistence function 2"),
    formula = expression(S == c*A^z * exp(-d/A)),
    exp = expression(c*A^z * exp(-d/A)),
    mod = s~c*a^z * exp(-d/a),
    shape = "sigmoid",
    asymp = function(pars) FALSE,
    custStart = function(data)c(5,.25,.15),
    #limits for parameters
    parLim = c("Rplus","Rplus","R"),
    #initials values function
    init = function(data){
      if (any(data$S==0)){
        log.data = data.frame(A=log(data$A),S=log(data$S+.5))
      } else log.data = log(data)
      res = stats::lm(S~A,log.data)$coefficients
      c(res[1L], res[2L], .1)
    }
  )
}

# POWER MODEL (ARRHENIUS 1921)
model_power <- function(data) {
  list(
    name = c("Power"),
    formula = expression(S == c * A ^ z),
    exp = expression(c * A ^ z),
    shape="convex",
    asymp= function(pars) FALSE,
    parLim = c("R", "R"),
    custStart = function(data) c(5,0.25),
    init = function(data) {
      if (any(data$S == 0)){
        log.data = data.frame(A = log(data$A), S = log(data$S + .5))
      } else log.data = log(data)
    res = stats::lm(S ~ A, log.data)$coefficients
    res = c(exp(res[1]), res[2])
    names(res) = c("c", "z")
    res
    }
  )
}

# POWER MODEL BIS (ROSENSWEIG 1995)
model_powerR <- function(data) {
  list(
    name = "PowerR",
    formula = expression(S == f + c*A^z),
    exp = expression(f + c*A^z),
    shape = "convex",
    asymp = function(pars) FALSE,
    custStart = function(data) c(5,.25,0),
    #limits for parameters
    parLim = c("R", "R", "R"),
    #initials values function
    init = function(data) {
      if (any(data$S == 0)) {
        log.data <- data.frame(A = log(data$A), S = log(data$S + .5))
      } else log.data <- log(data)
      res <- stats::lm(S~A,log.data)$coefficients
      res <- c(exp(res[1L]), res[2L])
      res <- c(stats::lm(S~A,data)$coefficients[1], res)
      names(res) <- c("f","c","z")
      return(res)
    }
  )
}

#RATIONAL FUNCTION ratkowski (1990)
model_ratio <- function(data) {
  list(
    name = c("Rational function"),
    formula = expression(S == (c + z*A)/(1+d*A)),
    exp = expression((c + z*A)/(1+d*A)),
    shape = "convex",
    asymp = function(pars)pars["z"]/pars["d"],
    parLim = c("R","Rplus","unif"),
    custStart = function(data) c(1,1,0.000001),
    #initial values function
    init = function(data) c(0,0,.5)
  )
}

#CUMULATIVE WEIBULL DISTRIBUTION with 3 parameters
model_weibull3 <- function(data) {
  list(
    name = "Cumulative Weibull 3 par.",
    formula = expression(S == d(1 - exp(-c*A^z)) ),
    exp = expression(d*(1 - exp(-c*A^z)) ),
    shape = "sigmoid",
    asymp = function(pars) pars["d"],
    parLim = c("Rplus", "Rplus", "Rplus"),
    custStart = function(data) c(10, .01, max(data$S)),
    init = function(data) {
      data <- data[data$S!=0,]
      #c calculation (asymptote)
      d <- max(data$S)+max(data$S)/4
      #z calculation
      Z <- log(-log((d-data$S)/d)) # = log z + flogX
      Z[][is.infinite(Z)] <- NA
      c <- exp(min(Z))
      dat <- data.frame("A"=log(data$A), "S"=Z)
      c <- exp(stats::lm(S~A,dat)$coefficients[[1L]])
      #f calculation
      z <- stats::lm(S~A,dat)$coefficients[[2L]]
      c(d,c,z)
    }
  )
}

# CUMULATIVE WEIBULL DISTRIBUTION with 4 parameters
model_weibull4 <- function(data) {
  list(
    name = "Cumulative Weibull 4 par.",
    formula = expression(S == d * (1 - exp(-c*A^z))^f ),
    exp = expression(d * (1 - exp(-c*A^z))^f ),
    shape = "sigmoid",
    asymp = function(pars) pars["d"],
    parLim = c("Rplus", "Rplus", "Rplus", "Rplus"),
    custStart = function(data)c(10,.01,max(data$S),.01),
    #initial values function
    init = function(data) {
      data <- data[data$S != 0,]
      #c calculation (asymptote)
      d <- max(data$S) + .25*max(data$S)
      #z calculation
      Z <- log(-log((d-data$S)/d)) # = log z + flogX
      Z[][is.infinite(Z)] <- NA
      c <- exp(min(Z))
      dat <- data.frame("A" = log(data[[1]]),"S"=Z)
      c <- exp(stats::lm(S~A,dat)$coefficients[[1L]])
      #f calculation
      z <- stats::lm(S~A,dat)$coefficients[[2L]]
      c(d, c, z, 1)
    }
  )
}
