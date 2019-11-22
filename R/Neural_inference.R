
rm(target,param,sumstat,nss,ltransf,numparam,cond1,cond2,transf,gwt,paramnames,statnames,scaled.sumstat,sum1,nacc,ds,wt1,aux,ss,unadj.values,statvar,param.mad,pred,fv0,fit0,pred.med,fitted.values,fitted.values0,residuals,pred1,fv1,fit1,pred1,pred.sd,adj.values,fv.sd,unadj.values)
target=a
param=b
sumstat=c

trainSSq0=cbind(b,c)

tol=0.5
iteration=20

neurg=neuralgression(a, b,c,tol=0.1,iteration=10, hcorr = TRUE,
                     transf = "log", kernel = "epanechnikov", 
                    hidden = c(40,40), activationfun = "sigm", learningrate = 0.8,learningrate_scale=1,
                     momentum = 0.5,  output = "sigm", sae_output = "sigm",
                     numepochs = 100, hidden_dropout = 0.1, visible_dropout = 0.1)


deepnetgression <- function(target, param, sumstat, tol, iteration, hcorr = TRUE,
                           transf = "none", logit.bounds = c(0,0), subset = NULL,hidden, activationfun, learningrate,learningrate_scale,
                           momentum, output, sae_output, numepochs,  hidden_dropout, visible_dropout,kernel ) {
  
  call <- match.call()
  
  ## general checks that the function is used correctly
  ## ###################################################
  
  if(missing(target)) stop("'target' is missing")
  if(missing(param)) stop("'param' is missing")
  if(missing(sumstat)) stop("'sumstat' is missing")
  if(!is.matrix(param) && !is.data.frame(param) && !is.vector(param)) stop("'param' has to be a matrix, data.frame or vector.")
  if(!is.matrix(sumstat) && !is.data.frame(sumstat) && !is.vector(sumstat)) stop("'sumstat' has to be a matrix, data.frame or vector.")
  if(missing(tol)) stop("'tol' is missing")
  
  
  if(!any(kernel == c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine"))){
    kernel <- "epanechnikov"
    warning("Kernel is incorrectly defined. Setting to default (Epanechnikov)")
  }
  
  if(is.data.frame(param)) param <- as.matrix(param)
  if(is.data.frame(sumstat)) sumstat <- as.matrix(sumstat)
  if(is.list(target)) target <- unlist(target)
  if(is.vector(sumstat)) sumstat <- matrix(sumstat, ncol=1)
  if(length(target)!=dim(sumstat)[2]) stop("Number of summary statistics in 'target' has to be the same as in 'sumstat'.")
  
  ## stop if zero var in sumstat
  ## #########################
  nss <- length(sumstat[1,])
  cond1 <- !any(as.logical(apply(sumstat, 2, function(x) length(unique(x))-1)))
  if(cond1) stop("Zero variance in the summary statistics.")
  
  
  ## transformations
  ## ################
  ltransf <- length(transf)
  if(is.vector(param)){
    numparam <- 1
    param <- matrix(param, ncol=1)
  }
  else numparam <- dim(param)[2]
  for (i in 1:ltransf){
    if(sum(transf[i] == c("none","log","logit")) == 0){
      stop("Transformations must be none, log, or logit.")
    }
    if(transf[i]=="logit"){
      if(logit.bounds[i,1] >= logit.bounds[i,2]){
        stop("Logit bounds are incorrect.")       
      }
    }
  }
  
  if (numparam != ltransf) {
    if (length(transf) == 1) {
      transf <- rep(transf[1], numparam)
      warning("All parameters are \"", transf[1], "\" transformed.", 
              sep = "", call. = F)
    }
    else stop("Number of parameters is not the same as number of transformations.", 
              sep = "", call. = F)
  }
  ## if no logit change logit.bounds to NULL
  ## if(any(transf) == "logit") logit.bounds <- NULL
  
  ## parameter and/or sumstat values that are to be excluded 
  ## ####################################################### 
  gwt <- rep(TRUE,length(sumstat[,1])) 
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE 
  if(missing(subset)) subset <- rep(TRUE,length(sumstat[,1])) 
  gwt <- as.logical(gwt*subset) 
  
  ## extract names of parameters and statistics if given 
  ## ################################################### 
  if(!length(colnames(param))){ 
    warning("No parameter names are given, using P1, P2, ...") 
    paramnames <- paste("P", 1:numparam, sep="") 
  } 
  else paramnames <- colnames(param) 
  
  if(!length(colnames(sumstat))){ 
    warning("No summary statistics names are given, using S1, S2, ...") 
    statnames <- paste("S", 1:nss, sep="") 
  } 
  else statnames <- colnames(sumstat) 
  
  normalise <- function(x,y){
    if(mad(y) == 0 || is.na(mad(y)))
      return (x)
    else 
      return (x/mad(y))
  }
  
  ## scale everything
  ## #################
  
  scaled.sumstat <- sumstat
  for(j in 1:nss){
    scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
  }
  
  for(j in 1:nss){
    target[j] <- normalise(target[j],sumstat[,j][gwt])
  }
  ## calculate euclidean distance
  ## ############################
  sum1 <- 0
  for(j in 1:nss){
    sum1 <- sum1 + (scaled.sumstat[,j]-target[j])^2
  }
  dist <- sqrt(sum1)
  
  ## includes the effect of gwt in the tolerance
  dist[!gwt] <- floor(max(dist[gwt])+10)
  
  ## wt1 defines the region we're interested in
  
  nacc=ceiling(length(dist)*tol)
  ds=sort(dist)[nacc]
  wt1 <- (dist <= ds)
  aux<-cumsum(wt1)
  wt1 <- wt1 & (aux<=nacc)
  if(kernel == "gaussian") 
  {
    wt1 <- rep(TRUE, length(dist))
  }
  ## transform parameters
  ## ######################
  for (i in 1:numparam){
    if(transf[i] == "log"){
      if(min(param[,i]) <= 0){
        cat("log transform: values out of bounds - correcting...")
        x.tmp <- ifelse(param[,i] <= 0,max(param[,i]),param[,i])
        x.tmp.min <- min(x.tmp)
        param[,i] <- ifelse(param[,i] <= 0, x.tmp.min,param[,i])
      }
      param[,i] <- log(param[,i])
    }
    else if(transf[i] == "logit"){
      if(min(param[,i]) <= logit.bounds[i,1]){
        x.tmp <- ifelse(param[,i] <= logit.bounds[i,1],max(param[,i]),param[,i])
        x.tmp.min <- min(x.tmp)
        param[,i] <- ifelse(param[,i] <= logit.bounds[i,1], x.tmp.min,param[,i])
      }
      if(max(param[,i]) >= logit.bounds[i,2]){
        x.tmp <- ifelse(param[,i] >= logit.bounds[i,2],min(param[,i]),param[,i])
        x.tmp.max <- max(x.tmp)
        param[,i] <- ifelse(param[,i] >= logit.bounds[i,2], x.tmp.max,param[,i])
      }
      param[,i] <- (param[,i]-logit.bounds[i,1])/(logit.bounds[i,2]-logit.bounds[i,1])
      param[,i] <- log(param[,i]/(1-param[,i]))
    }
  } # end of parameter transformations
  
  
  ## select summary statistics in region
  ## ###################################
  ss <- sumstat[wt1,]
  unadj.values <- param[wt1,]
  statvar <- as.logical(apply(cbind(sumstat[wt1, ]), 2, function(x) length(unique(x)) -  1))
  cond2 <- !any(statvar)
  if(cond2) cat("Warning messages:\nStatistic(s)", statnames[!statvar], "has/have zero variance in the selected region.\nConsider using larger tolerance or the rejection method or discard this/these statistics, which might solve the collinearity problem in 'lsfit'.\n", sep=", ")
  
  ## weights
  if(kernel == "epanechnikov") weights <- 1 - (dist[wt1]/ds)^2
  if(kernel == "rectangular") weights <- dist[wt1]/ds
  if(kernel == "gaussian") weights <- 1/sqrt(2*pi)*exp(-0.5*(dist/(ds/2))^2)
  if(kernel == "triangular") weights <- 1 - abs(dist[wt1]/ds)
  if(kernel == "biweight") weights <- (1 - (dist[wt1]/ds)^2)^2
  if(kernel == "cosine") weights <- cos(pi/2*dist[wt1]/ds)
  
  ## neural network regression
  ## ##########################
  
  linout <- TRUE
  
  ## normalise parameters
  param.mad <- c()
  for(i in 1:numparam){
    param.mad[i] <- mad(param[,i][gwt]) # save for renormalisation
    param[,i] <- normalise(param[,i],param[,i][gwt])
  }
  
  #lambda <- sample(lambda, iteration, replace=T)
  # fv0 <- array(dim=c(sum(wt1), numparam, iteration)) ## this is the original, may be wrong
  # fv <- array(dim=c(sum(wt1), length(param), iteration)) ## fv[,col,row]
  fv0 <- array(dim=c(length(param[,1]), numparam, iteration)) #
  pred <- matrix(nrow=numparam, ncol=iteration)
  
  for (i in 1:iteration) {
    
    fit0=sae.dnn.train(scaled.sumstat, param, hidden, activationfun, learningrate,
                       momentum, learningrate_scale, output, sae_output, numepochs, batchsize = length(param[,1]), hidden_dropout , visible_dropout)
    
    # fit0=sae.dnn.train(scaled.sumstat[wt1,], param[wt1,], hidden, activationfun, learningrate,
    #                   momentum, learningrate_scale = 1, output, sae_output, numepochs, batchsize = length(param[wt1,1]), hidden_dropout , visible_dropout)
    cat(i)
    fv0[,,i]=fit0$e
    #  fv01[,,i]=fit0$L
    pred[,i] <- nn.predict(fit0, data.frame(rbind(target)))
  }
  
  cat("\n")
  pred.med <- apply(pred, 1, median)
  pred.med <- matrix(pred.med, nrow=length(param[,1]), ncol=numparam, byrow=T)
  #fitted.values0 <- apply(fv0[,,], c(1,2), median)
  fitted.values0 <- apply(fv0, c(1,2), median)
 # residuals <- param[wt1,] - fitted.values0 # median of fitted values par nnets for each accepted point and parameter
  residuals <- param - fitted.values0
  if(hcorr == TRUE){
    pred1 <- matrix(nrow=numparam, ncol=iteration)
    # fv1 <- array(dim=c(sum(wt1), numparam, iteration))
    fv1 <- array(dim=c(length(param[,1]), numparam, iteration))
    # fv2 <- array(dim=c(sum(wt1), length(param), iteration))
    for (i in 1:iteration){
      #fit1=sae.dnn.train(scaled.sumstat, log(residuals^2), hidden=c(40,40), activationfun="sigm", learningrate=0.8, momentum=3, learningrate_scale = 1, output="sigm", sae_output="sigm", numepochs=100, batchsize = length(param[,1]), hidden_dropout=0.2 , visible_dropout=0.2)
    fit1=sae.dnn.train(scaled.sumstat, log(residuals^2), hidden, activationfun, learningrate ,
                        momentum , learningrate_scale, output, sae_output, numepochs, batchsize =length(param[,1]) , hidden_dropout, visible_dropout)
      # fit1=sae.dnn.train(scaled.sumstat[wt1,], log(residuals^2), hidden, activationfun, learningrate ,
      #                   momentum ,  output, sae_output, numepochs, batchsize = sum(wt1), hidden_dropout, visible_dropout) 
     # fit0=sae.dnn.train(scaled.sumstat[wt1,], param[wt1,], hidden=c(40,40), activationfun="sigm", learningrate=0.8, momentum=3, learningrate_scale = 1, output="sigm", sae_output="sigm", numepochs=100, batchsize = length(param[wt1,1]), hidden_dropout=0.2 , visible_dropout=0.2)
      cat(i)
      fv1[,,i]=fit1$e
      #  fv02[,,i]=fit0$L
      pred1[,i] <- nn.predict(fit1, data.frame(rbind(target)))
    }
    cat("\n")
    pred.sd <- sqrt(exp(apply(pred1, 1, median)))
    pred.sd <- matrix(pred.sd, nrow=length(param[,1]), ncol=numparam, byrow=T)
    fv.sd <- sqrt(exp(apply(fv1, c(1,2), median)))
    # fv.sd <- sqrt(exp(apply(fv1, c(1), median)))
    # adj.values <- as.vector(pred.med) + (as.vector(pred.sd)*as.vector(residuals))/fv.sd # correction heteroscedasticity
    # residuals<-(as.vector(pred.sd)*as.vector(residuals))/fv.sd
   ## adj.values <- pred.med + (as.vector(pred.sd) *as.vector(residuals))/as.vector(fv.sd)
  #  residuals <- (as.vector(pred.sd) * as.vector(residuals))/as.vector(fv.sd)
    adj.values <- pred.med + (pred.sd *residuals)/fv.sd
    residuals <- (pred.sd * residuals)/fv.sd
    }
  else{
    adj.values <- pred.med + residuals
  }
  colnames(adj.values) <- colnames(unadj.values)
  
  ## renormalise
  for(i in 1:numparam){
    ##       residuals[,i] <- residuals[,i]*param.mad[i] not much sense...
    adj.values[,i] <- adj.values[,i]*param.mad[i] ## this is only apply to one parameters
  }
  
  # end of neuralnet
  
  
  
  
  ## back transform parameter values
  ## ################################
  if(numparam == 1){
    unadj.values <- matrix(unadj.values, ncol=1)
  }
  
  for (i in 1:numparam){
    if(transf[i] == "log"){
      unadj.values[,i] <- exp(unadj.values[,i])
      adj.values[,i] <- exp(adj.values[,i])
    }
    else if(transf[i] == "logit"){
      unadj.values[,i] <- exp(unadj.values[,i])/(1+exp(unadj.values[,i]))
      unadj.values[,i] <- unadj.values[,i]*(logit.bounds[i,2]-logit.bounds[i,1])+logit.bounds[i,1]
      adj.values[,i] <- exp(adj.values[,i])/(1+exp(adj.values[,i]))
      adj.values[,i] <- adj.values[,i]*(logit.bounds[i,2]-logit.bounds[i,1])+logit.bounds[i,1]
    }
  }
  
  
  
  
  neuralregression.return <- function(transf, logit.bounds, method, call, numparam, nss, paramnames, statnames,
                                      unadj.values, adj.values, ss, weights, residuals, dist, wt1, gwt,  hcorr,aic,bic){
    
    
    
    out <- list(adj.values=adj.values, unadj.values=unadj.values,
                ss=ss, weights=weights, residuals=residuals, dist=dist,
                call=call, na.action=gwt, region=wt1, transf=transf, logit.bounds = logit.bounds,
                hcorr = hcorr,  numparam=numparam, numstat=nss,
                names=list(parameter.names=paramnames, statistics.names=statnames))
    
    
    class(out) <- "neuralregression"
    out
  }
  neuralregression.return(transf, logit.bounds, method, call, numparam, nss, paramnames, statnames,
                          unadj.values, adj.values, ss, weights, residuals, dist, wt1, gwt,  hcorr,aic,bic)
  
  
}








