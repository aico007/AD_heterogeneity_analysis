install.packages("glmnet")
library(glmnet)

data = read.table("skin_PBMC_clinical_data_combined_.txt",sep="\t",header=T,row.names=1) 
colnames(data)[1:3]
# "EASI_total", "EASI_erythema", "EASI_papule"
data_EASI = data[ ,!colnames(data) %in% c("EASI_erythema","EASI_papule")]
colnames(data_EASI)[1]="trait"

regression = function(x){
  set.seed(x)
  id.random = sample(nrow(data_EASI),round(nrow(data_EASI)*0.8)) 
  d.train = data_EASI[id.random,]
  d.test = data_EASI[setdiff(1:nrow(data_EASI),id.random),]
  
  m = cv.glmnet(as.matrix(d.train[,-1]), d.train[,1], alpha=0.5, family="gaussian", maxit=100000, nfolds=10) 
  best.lambda = m$lambda.min
  en.train = glmnet(as.matrix(d.train[,-1]), d.train[,1], lambda=best.lambda, alpha=0.5, family="gaussian")
  en.coef=as.matrix(drop(en.train$beta))
  en.coef=en.coef[en.coef[,1]!=0,,drop=F]#[order(abs(en.coef[,1]),decreasing=T),]
  d.train.select = d.train[,colnames(d.train) %in% c("trait",rownames(en.coef))]
  en.selected.lm = lm(d.train.select$trait~., data=d.train.select) 
  selected.coef = summary(en.selected.lm)$coefficients
  d.test$predict.en = predict(en.train, newx=as.matrix(d.test[,-1]), s=best.lambda)
  SS.test.total = sum((d.test$trait - mean(d.test$trait))^2)
  SS.test.residual.en = sum((d.test$trait - d.test$predict.en)^2) 
  test.rsq.en = 1 - SS.test.residual.en/SS.test.total
  
  lm.train = lm(d.train$trait~., data=d.train)
  d.test$predict.lm = predict(lm.train, as.data.frame(d.test[,2:ncol(d.train)]))
  SS.test.residual.lm   <- sum((d.test$trait - d.test$predict.lm)^2)
  test.rsq.lm <- 1 - SS.test.residual.lm/SS.test.total  #lm
  return(list(selected.coef, 
              data.frame(lambda=best.lambda,
                         train.rsq.en=summary(en.selected.lm)$adj.r.squared,
                         train.rsq.lm=summary(lm.train)$adj.r.squared,
                         test.rsq.en=test.rsq.en,
                         test.rsq.lm=test.rsq.lm)))
}

res1 = NULL
for (i in 1:100){
  res1 = c(res1,regression(i))
}

res2 = NULL
for (i in seq(2,200,2)){
  res2 = rbind(res2,res1[[i]])
}
#rsq
colMeans(res2) 


#longitudinal dataset
d.cross = read.table("cross-sectional_PBMC_clinical_data_combined.txt",sep="\t",header=T,row.names=1) 
d.longi = read.table("longitudinal_PBMC_clinical_data_combined.txt",sep="\t",header=T,row.names=1) 

regression_longi = function(x){
  set.seed(x) 
  m = cv.glmnet(as.matrix(d.cross[,-1]), d.cross[,1], alpha=0.5, family="gaussian", maxit=100000, nfolds=10) 
  best.lambda = m$lambda.min
  en.train = glmnet(as.matrix(d.cross[,-1]), d.cross[,1], lambda=best.lambda, alpha=0.5, family="gaussian") 
  
  d.longi$predict.en = predict(en.train, newx=as.matrix(d.longi[,colnames(d.longi)!="trait"]), s=best.lambda, alpha=0.5, family="gaussian")
  SS.test.total = sum((d.longi$trait - mean(d.longi$trait))^2)
  SS.test.residual.en = sum((d.longi$trait - d.longi$predict.en)^2) 
  test.rsq.en = 1 - SS.test.residual.en/SS.test.total
  
  lm.train = lm(d.cross$trait~., data=d.cross)
  d.longi$predict.lm = predict(lm.train, as.data.frame(d.longi[,2:ncol(d.cross)]))
  SS.test.residual.lm   <- sum((d.longi$trait - d.longi$predict.lm)^2)
  test.rsq.lm <- 1 - SS.test.residual.lm/SS.test.total  
  
  return(data.frame(lambda=best.lambda,
                         train.rsq.en=summary(en.selected.lm)$adj.r.squared,
                         train.rsq.lm=summary(lm.train)$adj.r.squared,
                         test.rsq.en=test.rsq.en,
                         test.rsq.lm=test.rsq.lm))
}

res3 = NULL
for (i in 1:100){
  res3 = rbind(res3,regression_longi(i))
}
#rsq
colMeans(res3)
