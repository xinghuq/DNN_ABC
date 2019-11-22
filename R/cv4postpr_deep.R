
target=as.matrix(test_paramSSq0[2,-(1:4)])
index=as.vector(trainingSSq0[,1])
sumstat=as.matrix(trainingSSq0[,-(1:4)])

cv4postpr_deep0=cv4postpr_deep(index, sumstat, postpr.out = postpr_deep1, nval=60, tols=0.1,
method="mnlogistic", subset = NULL, kernel = "epanechnikov",numnet = 10, sizenet = 5, lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit = 50,hidden=c(30,30), activationfun="sigm", learningrate=0.5,
momentum=0.5, batchsize=100,learningrate_scale=1, output="sigm", sae_output="sigm", numepochs=100,  hidden_dropout=0.1 , visible_dropout=0.1,threshold = 0.01,
stepmax = 1e+05, rep = 1, startweights = NULL,
learningrate.limit = NULL, learningrate.factor = list(minus = 0.5,plus = 1.2),  lifesign = "none",
lifesign.step = 1000, algorithm = "rprop+", err.fct = "sse",
act.fct = "logistic", linear.output = TRUE, exclude = NULL,
constant.weights = NULL, likelihood = FALSE)

summary(cv4postpr_deep0)

plot(cv4postpr_deep0)

cv4postpr_deep <- function(index, sumstat, postpr.out = NULL, nval, tols,method, subset = NULL, kernel = "epanechnikov",numnet = 10, sizenet = 5, lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit = 500,hidden, activationfun, learningrate,
                      momentum, batchsize,learningrate_scale, output, sae_output, numepochs, hidden_dropout , visible_dropout,threshold,
                      stepmax, rep , startweights,
                      learningrate.limit , learningrate.factor,lifesign,
                      lifesign.step, algorithm, err.fct,
                      act.fct, linear.output, exclude,
                      constant.weights, likelihood){
  
  linout <- TRUE
  ## checks:
  ## ######
  if(missing(nval)) stop("'nval' must be supplied.", call.=F)
  if(is.null(postpr.out) && missing(method)) stop("Method must be supplied when 'postpr.out' is NULL.", call.=F)
  if(length(index)!=na.omit(length(index))) stop("'index' contains missing values. Models must be specified for each simulation.", call.=F)
  
  if(!prod(table(index)>nval)) stop("'nval' has to be smaller or equal to number of simulations for any of the models. Choose a smaller 'nval'.", call.=F)
  
  ## set random seeds
  ## ################
  if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  
  ## define defaults:
  ## #################
  
  if(!is.null(postpr.out)){
    subset <- postpr.out$na.action
    method <- postpr.out$method
    kernel <- "epanechnikov"
  }
  
  ## checks, numbers of stats, params, sims
  if(is.vector(sumstat)){
    numstat <- 1
    sumstat <- matrix(sumstat, ncol=1)
  }
  else numstat <- dim(sumstat)[2]
  numsim <- length(index)
  
  ## names
  if(!is.null(postpr.out)){ # indexnames & statnames from postpr.out
    if(numstat != postpr.out$numstat || numsim != length(postpr.out$na.action)){
      stop("The number of summary statistics, or simulations provided in 'sumstat' are not the same as in 'postpr.out'.", call.=F)
    }
    else if(!prod(unique(index) %in% postpr.out$names$models)){
      stop("Models in 'index' are not the same as in 'postpr.out', or different names are used.", call.=F)
    }
    else if(!prod(colnames(sumstat) %in% postpr.out$names$statistics.names)){
      stop("Summary statistics in 'sumstat' are not the same as in 'postpr.out', or different names are used.", call.=F)
    }
    else{
      mymodels <- postpr.out$names$models
      statnames <- postpr.out$names$statistics.names
    }
  }
  else{ # statnames o/w
    mymodels <- levels(factor(index))
    if(length(colnames(sumstat))){
      statnames <- colnames(sumstat)
    }
    else{
      warning("No statistics names are given, using S1, S2, ...", call.=F, immediate=T)
      statnames <- paste("S", 1:numstat, sep="")
    }
  }
  
  ## ## order data by models
  ## index <- factor(index)
  ## sumstat <- sumstat[order(index), ]
  ## index <- index[order(index)]
  
  ## indices for the CV sample and check that the sample is not actually an NA
  gwt <- rep(TRUE,length(sumstat[,1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if(is.null(subset)) subset <- rep(TRUE,length(sumstat[,1]))
  gwt <- as.logical(gwt*subset)
  ## CV samples
  cvsamp <- unlist(tapply(c(1:length(index))[gwt], index[gwt], sample, nval))
  
  ## if tols is a vector have to loop through all values
  tols <- sort(tols)
  allprobs <- list()
  mycall <- list()
  for(mytol in tols){
    res <- matrix(ncol = length(unique(index)), nrow = length(cvsamp))
    for(i in 1:length(cvsamp)){
      ## things to over-write from original call: tolerances, target, index, sumstat
      mysamp <- cvsamp[i]
      mytrue <- index[mysamp]
      mytarget <- sumstat[mysamp,]
      myindex <- index[-mysamp]
      mysumstat <- sumstat[-mysamp,]
      mysubset <- subset[-mysamp]
      subres <- postpr_deep(target = mytarget, index = myindex, sumstat = mysumstat, tol=mytol,subset = mysubset, method = method, kernel = kernel,hidden=hidden, activationfun=activationfun, 
                            learningrate=learningrate, momentum=momentum, learningrate_scale=learningrate_scale, output=output, sae_output=sae_output, 
                            numepochs=numepochs, batchsize=length(myindex), hidden_dropout=hidden_dropout, visible_dropout=visible_dropout,threshold=threshold,
                            stepmax=stepmax, rep=rep , startweights=startweights,
                            learningrate.limit=learningrate.limit , learningrate.factor=learningrate.factor,lifesign=lifesign,
                            lifesign.step=lifesign.step, algorithm=algorithm, err.fct=err.fct,
                            act.fct=act.fct, linear.output=linear.output, exclude=exclude,
                            constant.weights=constant.weights, likelihood=likelihood)
      
      if(subres$method=="rejection") res[i,] <- summary.postpr(subres, print = F)$Prob
      if(subres$method=="mnlogistic") res[i, ] <- summary.postpr(subres, print = F)$mnlogistic$Prob
      if(subres$method=="neuralnet") res[i, ] <- summary.postpr(subres, print = F)$neuralnet$Prob
      if(subres$method=="deepnet") res[i, ] <- summary.postpr(subres, print = F)$deepnet$Prob
      if(subres$method=="deepneuron") res[i, ] <- summary.postpr(subres, print = F)$deepneuron$Prob
    }
    colnames(res) <- mymodels
    rownames(res) <- index[cvsamp]
    allprobs[[paste("tol", mytol, sep="")]] <- res
    allnames <- lapply(allprobs, apply, 1, function(xx) mymodels[which(xx==max(xx))])
    
    mycall[[paste("tol", mytol, sep="")]] <-
      call("postpr_deep", target = quote(target), index = quote(index), sumstat = quote(sumstat), tol= mytol,
           subset = quote(subset), method = subres$method, kernel = subres$kernel,hidden=hidden, activationfun=activationfun, 
           learningrate=learningrate, momentum=momentum, batchsize=batchsize,learningrate_scale=learningrate_scale, output=output, sae_output=sae_output, 
           numepochs=numepochs, hidden_dropout=hidden_dropout, visible_dropout=visible_dropout,threshold=threshold,
           stepmax=stepmax, rep=rep , startweights=startweights,
           learningrate.limit=learningrate.limit , learningrate.factor=learningrate.factor,lifesign=lifesign,
           lifesign.step=lifesign.step, algorithm=algorithm, err.fct=err.fct,
           act.fct=act.fct, linear.output=linear.output, exclude=exclude,
           constant.weights=constant.weights, likelihood=likelihood)
  }
  
  cv4postpr.out <-list(calls = mycall, cvsamples = cvsamp, tols = tols, true = index[cvsamp],estim = allnames,
         model.probs = allprobs,
         method = method, names = list(models = mymodels, statistics.names = statnames), seed = seed)
  
  class(cv4postpr.out) <- "cv4postpr_deep"
  invisible(cv4postpr.out)
  
}
