######### Multiple starting values optimization function

grid_start_fit <- function(model, data, n = 100, verb = FALSE) {
  
  if(length(model$parNames)<4){
    ns <- 100
  }else{
    ns <- 10
  }
  
  start.list <- lapply(model$parLim,function(x){
    res = switch(x,
                 R = sample(seq(-500,500),ns),
                 Rplus = seq(.1,500,length.out = ns),
                 unif = runif(ns)
                 )
    return(res)
  })
  
  names(start.list) <- model$parNames
  
  grid.start <- expand.grid(start.list)
  
  grid.start <- grid.start[sample.int(dim(grid.start)[1],n),]
  
  if (verb) cat("- running grid optim: ")
  
  fit.list <- apply(grid.start,1,function(x){
    if (verb) cat(".")
    tryCatch(rssoptim(model,data,custStart=x),error=function(e) list(value=NA))
    })
  
  fit.list <- as.list(fit.list)
  
  #fit.anal <- tryCatch(rssoptim(model,data),error=function(e) list(value=NA))
  
  #fit.list[[length(fit.list)+1]] <- fit.anal
  
  values <- unlist(lapply(fit.list,function(x){x$value}))
  
  min <- which.min(values)
  
  if(length(min) != 0) {
    fit <- fit.list[[min]]
  }else{
    fit <- NA
  }
  
  return(fit)
  
}#eo grid_start_fit
