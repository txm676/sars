
#### INTERNAL FUNCTION(S)


###complete models
### add derivatives, residual sum of square function ...

compmod <- function(model){

  #get the model expression
  mod.exp <- model$exp

  #intern variables (not parameters)
  modChars <- c(" ","+","-","*","/","^","(",")","l","o","g","s","i","n",
                "e","x","p","A","1")

  #extracting characaters from the equation
  chars <- unlist(strsplit(as.character(mod.exp),split=""))

  #parameters extraction
  model$parNames <- unique(chars[!is.element(chars,modChars)])
  names(model$parNames) <- seq_along(model$parNames)

  #1st & 2nd derivatives
  model$d1.exp <- stats::D(mod.exp,"A")
  model$d2.exp <- stats::D(model$d1.exp,"A")

  ##########################################################################
  ######### creating functions
  model$mod.fun <- create_fun_mod(model$exp, model$parNames)
  # first derivative function
  model$d1.fun <- create_fun_mod(model$d1.exp, model$parNames)
  # second derivative function
  model$d2.fun <- create_fun_mod(model$d2.exp, model$parNames)

  # rss function
  #you'd use opt=FALSE if you were actually providing this function with your own
  #estimates that you thus wouldn't want transformed before calculating rss. Is a
  #legacy option and not used at all now.
  model$rss.fun <- function(par, data, parLim = model$parLim, opt = TRUE){
    #cat(".....",opt,"\n")
    #cat("bef trans : ",par,"\n")
    if (opt) par <- backLink(par,parLim)
    #cat("aft trans : ",par,"\n")
    S <- data$S
    A <- data$A
    res <- sum((S - model$mod.fun(A,par))^2)
    #cat("RSS : ",res,"\n")
    res
  }#eo rss.fun

  #cat("# model constructed ...\n")
  #------------> here end the "model constructor"

  invisible(model)
}

create_fun_mod <- function(expr, nam) {
  function(A, par) {
    for (i in seq_along(par)) assign(nam[i], par[i])
    eval(expr)
  }
}




