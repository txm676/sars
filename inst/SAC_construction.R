

##PA matrix, with sites as cols and sp as rows

dat <- as.data.frame(matrix(round(runif(25, 0, 1),0), ncol = 5, nrow = 5))


##function work out cumulative richness of a set of sites

cum_sp <- function (dat){
  
  dat <- t(dat)
  if (nrow(dat) <= 1  || ncol(dat) <= 1) stop("Only 1 site or species present in dataframe")

  sr <- vector(mode = "numeric", length(nrow(dat)))
  sr[1] <- sum(dat[1,]) 
  add <- c()
  
  for (i in 2:nrow(dat)){
    for (j in 1:ncol(dat)){
      if ((dat[i, j] == 1) && (sum(dat[, j][1:(i - 1)]) == 0)){
        add[j] <- 1
      } else {
        add[j] <- 0
      }
    } #eo j
    sr[i] <- sum(add)
    add <-c()
  } #eo i

sr <- cumsum(sr)
return (sr)
  
}


##########function to construct the SAC



SAC_construct <- function(dat, method = c("PR")){
  if (!method %in% c("PR", "RP", "SL", "LS")) stop("Incorrect method argument provided")
  
  dplyr::arrange()
  
}



data2 <- data2[(order(data2$sp.r,decreasing=FALSE)),]

###for and if loop to accumulate area if vector is ordered in terms of sp.r

acar_pr <- c()

for (i in 1:length(data2$area)){
  
  if (i == 1)
    acar_pr[i] = data2$area[i]
  
  else acar_pr[i] <- acar_pr[i-1]+data2$area[i]
  
}


e=data[,2:ncol(data_fc)]

e2 <- e[,(order(e[spn,],decreasing=FALSE))]

e2 <- cbind(data_fc[,1],e2)

prc <- cum_sp(e2)