is.cv4postpr <- function(x){
  if (inherits(x, "cv4postpr_deep")) TRUE
  else FALSE
}

summary.cv4postpr_deep <- function(object, probs=TRUE, print = TRUE, digits = max(3, getOption("digits")-3)){
  
  if (!inherits(object, "cv4postpr_deep")) 
    stop("Use only with objects of class \"cv4postpr_deep\".", call.=F)
  
  cv4postpr.out <- object
  tols <- cv4postpr.out$tols
  numtols <- length(tols)
  true <- cv4postpr.out$true
  estim <- cv4postpr.out$estim
  method <- cv4postpr.out$method
  model.probs <- cv4postpr.out$model.probs
  nmodels <- length(cv4postpr.out$names$models)
  nval <- length(true)/nmodels
  
  if(print) cat("Confusion matrix based on ", nval, " samples for each model.\n\n", sep="")
  cm <- lapply(estim, function(x) table(true, x))
  cm <- lapply(cm, function(x) {attributes(dimnames(x))$names <- NULL; x})
  if(print) print(cm); cat("\n")
  
  if(probs){
    if(print) cat(paste("Mean model posterior probabilities (", method ,")\n\n", sep="")) 
    myprs <- lapply(model.probs, apply, 2, tapply, true, mean)
    if(print) print(lapply(myprs, round, digits=digits))
    out <- list(conf.matrix=cm, probs=myprs)
  }
  else out <- cm
  
  invisible(out)
}

plot.cv4postpr_deep <- function(x, probs=FALSE, file = NULL, postscript = FALSE, onefile = TRUE, ask = !is.null(deviceIsInteractive()), caption = NULL, ...){
  
  if (!inherits(x, "cv4postpr_deep")) 
    stop("Use only with objects of class \"cv4postpr_deep\".", call.=F)
  
  cv4postpr.out <- x
  tols <- cv4postpr.out$tols
  numtols <- length(tols)
  true <- cv4postpr.out$true
  estim <- cv4postpr.out$estim
  method <- cv4postpr.out$method
  model.probs <- cv4postpr.out$model.probs
  nmodels <- length(cv4postpr.out$names$models)
  nval <- length(true)/nmodels
  
  if(is.null(caption)) caption <- "Confusion matrix"
  
  ## Devices
  save.devAskNewPage <- devAskNewPage()
  if(!is.null(file)){
    file <- substitute(file)
    if(!postscript) pdf(file = paste(file, "pdf", sep="."), onefile=onefile)
    if(postscript) postscript(file = paste(file, "ps", sep="."), onefile=onefile)
  }
  else{
    if (ask && 1 < numtols) {
      devAskNewPage(TRUE)
    }
  }
  
  par(cex = 1, cex.main = 1.2, cex.lab = 1.1)
  for(i in 1:numtols){
    if(probs){
      mym <- lapply(model.probs, apply, 2, tapply, true, mean)
      barplot(t(mym[[paste("tol", tols[i], sep="")]]), ylab="Mean model probability", ...)
    }
    else{
      mym <- lapply(estim, table, true)
      barplot(mym[[paste("tol", tols[i], sep="")]], ylab="Frequency", ...)
    }
    title(caption, sub=paste("Tolerance rate = ", tols[i], sep=""))
  }
  
  if(!is.null(file)){
    dev.off()
  }
  else devAskNewPage(save.devAskNewPage)
  invisible()
  
}
