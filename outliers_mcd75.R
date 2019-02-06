### multivariate outliers detection based one paper, 
### link: https://www.sciencedirect.com/science/article/pii/S0022103117302123
###############################################################################
library(MASS)
  datapath='D:/fMRI_data/Rscripts/data'
  setwd(datapath)
  listcsv<-dir(pattern='*.csv')
  
  for(file in 1:length(listcsv)){
    alldata<-read.csv(listcsv[file],header=TRUE)
    #alldata=read.csv('sub1_nodata.csv',header=TRUE)
    data<-na.exclude(alldata)
    output75<-cov.mcd(data,quantile.used=nrow(data)*.75)
    md<-mahalanobis(data,colMeans(data),cov(data))
    mhmcd75<-mahalanobis(data,output75$center,output75$cov)
    alpha=.01
    cutoff<-(qchisq(p=1-alpha,df=2))
    names_outliers_MCD75<-which(mhmcd75>cutoff)
    exclued<-names_outliers_MCD75
    data2<-data[-exclued,]
    write.csv(data2,paste0('outliers_stripping_',listcsv[file]),row.names = FALSE)
  }
