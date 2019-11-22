#### This is inference function withou seting torlance, BP,Backpropagation,sae:Stacked AutoEncoder,DBN,Deep belief network
rm(nntrain,target,param,sumstat,nss,Nor_sumstat,Nor_param,Nor_target, ltransf,numparam,cond1,cond2,transf,gwt,paramnames,statnames,scaled.sumstat,sum1,nacc,ds,wt1,aux,ss,unadj.values,statvar,param.mad,pred,fv0,fit0,pred.med,fitted.values,fitted.values0,residuals,pred1,fv1,fit1,pred1,pred.sd,adj.values,fv.sd,unadj.values)

encodernet1=encodernet(as.matrix(test_paramSSq0[1,-(1:4)]), as.matrix(trainingSSq0[,2:4]), as.matrix(trainingSSq0[,-(1:4)]), iteration=10, encoder="sae",hidden = c(40,40), activationfun = "sigm", learningrate = 0.8,
  momentum = 0.5, learningrate_scale=1, output = "sigm", sae_output = "sigm",
  numepochs = 100, batchsize = 400, hidden_dropout = 0.1, visible_dropout = 0.1)


encodernet2=encodernet(as.matrix(test_paramSSq0[1,-(1:4)]), as.matrix(trainingSSq0[,2:4]), as.matrix(trainingSSq0[,-(1:4)]), iteration=10, encoder="db",hidden = c(40,40), activationfun = "sigm", learningrate = 0.8,
                       momentum = 0.5, learningrate_scale=1, output = "sigm", 
                       numepochs = 100, batchsize = 400, hidden_dropout = 0.1, visible_dropout = 0.1)


encodernet3=encodernet(as.matrix(test_paramSSq0[1,-(1:4)]), as.matrix(trainingSSq0[,2:4]), as.matrix(trainingSSq0[,-(1:4)]), iteration=10, encoder="None",hidden = c(40,40), activationfun = "sigm", learningrate = 0.8,
                       momentum = 0.5, learningrate_scale=1, output = "sigm", 
                       numepochs = 100, batchsize = 400, hidden_dropout = 0.1, visible_dropout = 0.1)

nn=nn.train(as.matrix(trainingSSq0[,2:4]), as.matrix(trainingSSq0[,-(1:4)]),hidden = c(40,40), activationfun = "sigm", learningrate = 0.8,
         momentum = 0.5, learningrate_scale=1, output = "sigm", 
         numepochs = 100, batchsize = 400, hidden_dropout = 0.1, visible_dropout = 0.1)



target=as.matrix(test_paramSSq0[1,-(1:4)])
param=as.matrix(trainingSSq0[,2:4])
sumstat=as.matrix(trainingSSq0[,-(1:4)])


##### Note: Adding nn.tarin directly and do not normolize the parameter if they range from 0-1

### we can set the batchsize = sum(wt1) for if there is tol,

encodernet=function(target, param, sumstat,iteration,encoder,hidden , activationfun , learningrate,
         momentum,  learningrate_scale,output , sae_output,numepochs, batchsize, hidden_dropout, visible_dropout,cd){

Nor_sumstat=matrix(data=0,ncol=ncol(sumstat),nrow=nrow(sumstat))
Nor_target=matrix(data=0,ncol=ncol(target),nrow=nrow(target))
Nor_param=matrix(data=0,ncol=ncol(param),nrow=nrow(param))
for (i in 1:ncol(sumstat)) {
  Nor_sumstat[,i] <- normalise(sumstat[,i],sumstat[,i])
  Nor_target[,i] <- normalise(target[,i],target[,i])
}
param.mad=c()
for (j in 1:ncol(param)){
     param.mad[j] <- mad(param[,j])
     Nor_param[,j]=normalise(param[,j],param[,j]) 
}
colnames(Nor_sumstat)=colnames(sumstat)
colnames(Nor_target)=colnames(target)
colnames(Nor_param)=colnames(param)



#min=min(param) 
#max=max(param) 
#scale=max-min 
#Nor_param<- (param -min)/scale
#Nor_param=as.data.frame(Nor_param)


### batchsize = length(Nor_sumstat[,1])
residuals<- array(dim=c(length(Nor_sumstat[,1]), length(Nor_param[1,]), iteration))
fv0 <- array(dim=c(length(Nor_sumstat[,1]), length(Nor_param[1,]), iteration)) #length(Nor_param[1,]) is number of parameters ncol(param)
pred<- matrix(nrow=length(Nor_param[1,]), ncol=iteration)
if (encoder=="sae") {
  for(i in 1:iteration){
   
  #  saetrain=sae.dnn.train(as.matrix(Nor_sumstat), as.matrix(Nor_param),hidden = c(40,40), activationfun = "sigm", learningrate = 0.8,
   # momentum = 0.5,  output = "sigm", sae_output = "sigm", numepochs = 100, batchsize = 400, hidden_dropout = 0, visible_dropout = 0)
   saetrain=sae.dnn.train(as.matrix(Nor_sumstat), as.matrix(Nor_param), hidden , activationfun , learningrate,
                  momentum, learningrate_scale, output, sae_output,numepochs, batchsize, hidden_dropout, visible_dropout)
cat(i)
 residuals[,,i]=saetrain$e
 fv0[,,i]=saetrain$post[[4]]
pred[,i] <- nn.predict(saetrain, data.frame(rbind(Nor_target)))
}
pred.med <- apply(pred, 1, median)
pred.med <- matrix(pred.med, nrow=length(param[,1]), ncol=ncol(param), byrow=T)
fitted.values0 <- apply(fv0, c(1,2), median)
residuals=apply(residuals, c(1,2), median)
residuals1 <- Nor_param - fitted.values0
adj.values <- pred.med + residuals
adj.values1 <- pred.med + residuals1

}
if (encoder=="db") {
  
  for(i in 1:iteration){
    
  dbtrain=dbn.dnn.train(as.matrix(Nor_sumstat), as.matrix(Nor_param),  hidden , activationfun , learningrate,
                        momentum, learningrate_scale, output,numepochs, batchsize, hidden_dropout, visible_dropout)
  
 
 #fv0[,,i]=as.matrix(rowMeans(do.call(cbind, lapply(saetrain$net.result, data.frame))))
 cat(i)
  residuals[,,i]=dbtrain$e
 fv0[,,i]=dbtrain$post[[4]]
 pred[,i] <- nn.predict(dbtrain, data.frame(rbind(Nor_target)))
}
pred.med <- apply(pred, 1, median)
pred.med <- matrix(pred.med, nrow=nrow(param), ncol=ncol(param), byrow=T)
fitted.values0 <- apply(fv0, c(1,2), median)
residuals=apply(residuals, c(1,2), median)
residuals1 <- Nor_param - fitted.values0
adj.values <- pred.med + residuals
adj.values1 <- pred.med + residuals1
}

if (encoder=="None") {
  
  for(i in 1:iteration){
    
   # nntrain=nn.train(Nor_sumstat, Nor_param,  hidden , activationfun , learningrate,
    #                    momentum, learningrate_scale, output,numepochs, batchsize, hidden_dropout, visible_dropout)
    
   nntrain=nn.train(Nor_sumstat, Nor_param, hidden = c(40,40), activationfun = "sigm", learningrate = 0.8,
             momentum = 0.5, learningrate_scale=1, output = "sigm", 
            numepochs = 100, batchsize = 400, hidden_dropout = 0.1, visible_dropout = 0.1)
    #fv0[,,i]=as.matrix(rowMeans(do.call(cbind, lapply(saetrain$net.result, data.frame))))
    cat(i)
    residuals[,,i]=nntrain$e
    fv0[,,i]=nntrain$post[[4]]
    pred[,i] <- nn.predict(nntrain, data.frame(rbind(Nor_target)))
  }
  pred.med <- apply(pred, 1, median)
  pred.med <- matrix(pred.med, nrow=nrow(param), ncol=ncol(param), byrow=T)
  fitted.values0 <- apply(fv0, c(1,2), median)
  residuals=apply(residuals, c(1,2), median)
  residuals1 <- Nor_param - fitted.values0
  adj.values <- pred.med + residuals
  adj.values1 <- pred.med + residuals1
}
for(i in 1:ncol(param)){
  ##       residuals[,i] <- residuals[,i]*param.mad[i] not much sense...
  adj.values[,i] <- adj.values[,i]*param.mad[i] 
  adj.values1[,i] <- adj.values1[,i]*param.mad[i] 
  pred.med[,i]=pred.med[,i]*param.mad[i]
}
 return(list(adj.values=adj.values,adj.values1=adj.values1,pred.med=pred.med))

}

