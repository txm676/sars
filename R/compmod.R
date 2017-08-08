###complete models
### add derivatives, residual sum of square function ...

compmod <- function(model){

  #get the model expression
  mod.exp <- model$exp

  #intern variables (not parameters)
  modChars <- c(" ","+","-","*","/","^","(",")","l","o","g","s","i","n","e","x","p","A","1")

  #extracting characaters from the equation
  chars <- unlist(strsplit(as.character(mod.exp),split=""))

  #parameters extraction
  model$parNames <- unique(chars[!is.element(chars,modChars)])
  names(model$parNames) <- 1:length(model$parNames)

  #1rst & 2nd derivatives
  model$d1.exp <- D(mod.exp,"A")
  model$d2.exp <- D(model$d1.exp,"A")

  ##########################################################################
  ######### creating functions
  stringStartFun <- " <- function(A,par){ eval("
  stringStartFunRss <- " <- function(data,par,opt){ eval("
  stringStartList <- "list(A=A,"
  stringPars <- paste(model$parNames,"=par[",names(model$parNames),"]",sep="")
  stringEnd <- ")) }"

  #model function
  eval(parse(text=paste("model$mod.fun",stringStartFun,"model$exp,",stringStartList,paste(stringPars,collapse=","),stringEnd ,sep="")))

  #first derivative function
  eval(parse(text=paste("model$d1.fun",stringStartFun,"model$d1.exp,",stringStartList,paste(stringPars,collapse=","),stringEnd ,sep="")))

  #second derivative function
  eval(parse(text=paste("model$d2.fun",stringStartFun,"model$d2.exp,",stringStartList,paste(stringPars,collapse=","),stringEnd ,sep="")))

  #rss function
  model$rss.fun <- function(par,data,parLim=model$parLim,opt=TRUE){
    #cat(".....",opt,"\n")
    #cat("bef trans : ",par,"\n")
    if(opt) { par <- backLink(par,parLim) }
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
