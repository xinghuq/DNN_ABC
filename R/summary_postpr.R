

summary.postpr <- function(object, rejection = TRUE, print = TRUE, digits = max(3, getOption("digits")-3)){
  
  if (!inherits(object, "postpr_deep")) 
    stop("Use only with objects of class \"postpr_deep\".", call.=F)
  
  postpr.out <- object
  cl <- postpr.out$call
  npost <- length(postpr.out$values)
  pred <- postpr.out$pred
  allvals <- postpr.out$values
  postmod <- levels(postpr.out$values)
  nmod <- length(postmod)
  method <- postpr.out$method
  corr <- postpr.out$corr
  nmodels <- postpr.out$nmodels
  
  if(print){
    cat("Call: \n")
    dput(cl, control=NULL)  
    cat(paste("Data:\n postpr.out$values (",npost," posterior samples)\n", sep=""))
    cat(paste("Models a priori:\n "))
    cat(postpr.out$names$models, sep=", ")
    cat(paste("\nModels a posteriori:\n "))
    cat(postmod, sep=", ")
    if(corr & length(unique(nmodels))>1){
      cat("\n")
      warning("Posterior model probabilities are corrected for unequal number of simulations per models.", immediate.=T, call.=F)
    }
    cat("\n\n")
  }
  
  if(rejection || method == "rejection"){
    
    if(print) cat("Proportion of accepted simulations (rejection):\n")
    allpr <- table(allvals)/length(allvals)
    if(corr){
      ratio <- (allpr*npost) / nmodels
      allpr <- ratio/sum(ratio)
    }
    prnames <- dimnames(allpr)$allvals
    allpr <- c(allpr); names(allpr) <- prnames
    if(print) print(round(allpr, digits=digits))
    
    if(nmod>1){
      pr.rej <- table(allvals)/length(allvals)
      bf.rej <- t(matrix(pr.rej, nmod, nmod, byrow=T)/matrix(pr.rej, nmod, nmod, byrow=F))
      colnames(bf.rej) <- postmod
      rownames(bf.rej) <- postmod
      bf.rej <- as.table(bf.rej)
      if(print){
        cat("\nBayes factors:\n")
        print(round(bf.rej, digits=digits))
        cat("\n\n")
      }
    }
    else bf.rej <- NA
    
  }
  
  if(method == "mnlogistic" | method == "neuralnet"| method=="deepnet" | method=="deepneuron"){
    
    if(print){
      cat(paste("Posterior model probabilities (", method, "):\n", sep=""))
      print(round(pred, digits=digits))
    }
    if(nmod>1){
      bf.reg <- t(matrix(pred[pred!=0], nmod, nmod, byrow=T)/matrix(pred[pred!=0], nmod, nmod, byrow=F))
      colnames(bf.reg) <- postmod
      rownames(bf.reg) <- postmod
      bf.reg <- as.table(bf.reg)
      if(print){
        cat("\nBayes factors:\n")
        print(round(bf.reg, digits=digits))
        cat("\n")
      }
    }
    else bf.reg <- NA
    
    if(rejection){
      if(method == "mnlogistic")
        out <- list(rejection=list(Prob=allpr, BayesF=bf.rej), mnlogistic=list(Prob=pred, BayesF=bf.reg))
      if(method == "neuralnet")
        out <- list(rejection=list(Prob=allpr, BayesF=bf.rej), neuralnet=list(Prob=pred, BayesF=bf.reg))
      if(method == "deepnet")
        out <- list(rejection=list(Prob=allpr, BayesF=bf.rej), deepnet=list(Prob=pred, BayesF=bf.reg))
      if(method == "deepneuron")
        out <- list(rejection=list(Prob=allpr, BayesF=bf.rej), deepneuron=list(Prob=pred, BayesF=bf.reg))

      }
    
    else{
      out <- list(Prob=pred, BayesF=bf.reg)
    }
  }
  else out <- list(Prob=allpr, BayesF=bf.rej)
  invisible(out)
}
