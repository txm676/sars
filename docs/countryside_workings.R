
#add starting parameters for sar_habitat

#sort out mes and printing of extrap function

#Add in logarithmic version and maybe logistic?

#some sort of plot function that looks at extinction
#rates if you change amount of habitat!?

library(sars)
data <- read.table("D:\\documents\\Work\\On-going projects\\Habitat SAR Models\\Countryside files\\Vania\\MatrixMHUb.txt", 1)
data <- data[,c(2:4, 6,8,7,9)]


f <- sar_countryside(data, ubiSp = TRUE,
                     habNam = c("AG", "SH", "F"),
                     grid_start = "partial")

countryside_extrap(f, 1:3)




##user startPar

SP <- matrix(rep(c(3.061161e+08, 2.104536e-01,
               1.074748e+00, 1.223700e-01),4),
             ncol = 4) %>% t()

f <- sar_countryside(data, ubiSp = TRUE,
                     habNam = c("AG", "SH", "F"),
                     startPar  = SP)






#define CSAR 2 - used to estimate total species
CountrysideSAR2<-function(a1,a2,a3,par)
{
  species<-(a1*par[1]+par[2]*a2+par[3]*a3)^par[4]
  names(species)<-NULL
  species
}

# determine total species
TotalCSAR <- function (a1,a2,a3, f)
{
  sp.AG<-CountrysideSAR2(a1,a2,a3,coef(f[[1]]))
  sp.UB<-CountrysideSAR2(a1,a2,a3,coef(f[[4]]))	
  sp.SH<-CountrysideSAR2(a1,a2,a3,coef(f[[2]]))
  sp.QF<-CountrysideSAR2(a1,a2,a3,coef(f[[3]]))
  sp.AG+sp.SH+sp.QF+sp.UB
}

TotalCSAR(1,2,3,f)



