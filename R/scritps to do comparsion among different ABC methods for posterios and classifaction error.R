fv0[,,i]=as.matrix(rowMeans(do.call(cbind, lapply(saetrain$net.result, data.frame))))
library(abc)

# use the remaining 99% of data to training data
trainSS_all=as.matrix(SSall_dataset[,-(1:4)])
trainSS_q0= as.matrix(SSq0_dataset[,-(1:4)])
trainSS_q1= as.matrix(SSq1_dataset[,-(1:4)])
trainSS_q2= as.matrix(SSq2_dataset[,-(1:4)])

train_paramSSq0=as.matrix(SSq0_dataset[,2:4])
train_paramSSq1=as.matrix(SSq1_dataset[,2:4])
train_paramSSq2=as.matrix(SSq2_dataset[,2:4])
train_paramSSall=as.matrix(SSall_dataset[,2:4])
##five test data

target=as.matrix(validationSSq0[2,-(1:4)])
target=as.matrix(validationSSq1[1,-(1:4)])
target=as.matrix(validationSSq2[1,-(1:4)])
target=as.matrix(validationSSall[1,-(1:4)])

index=as.vector(SSall_dataset[,1])

a=as.matrix(validationSSall[1,-(1:4)])
b=as.matrix(SSall_dataset[,2:4])
c=as.matrix(SSall_dataset[,-(1:4)])

write.csv(validationSSq1,file="observed_SSq1.csv")
write.csv(validationSSall,file="observed_SSall.csv")

##### now doing inference using different methods and different SS

abcSSq1_rejection=abc1(as.matrix(validationSSall[1,-(1:4)]), as.matrix(SSall_dataset[,2:4]), as.matrix(SSall_dataset[,-(1:4)]), tol=0.1, method="rejection", hcorr = TRUE, transf = c("log", "none", "log" ))
abcSSq1_rejection$unadj.values
abcSSq1_rejection$unadj.values[,1]-validationSSall[1,2]
abcSSq1_loclinear=abc1(as.matrix(validationSSall[1,-(1:4)]), as.matrix(SSall_dataset[,2:4]), as.matrix(SSall_dataset[,-(1:4)]), tol=0.1, method="loclinear", hcorr = TRUE, transf = c("log", "none", "log" ),
                       subset = NULL, kernel = "epanechnikov", numnet =
                        10, sizenet = 5, lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit =
                        50)
abcSSq1_loclinear=abc1(a,b,c,tol=0.1,method="loclinear",transf = c("log", "none", "log" ), hcorr = TRUE,kernel = "epanechnikov")
abcSSq1_loclinear$adj.values
abcSSq1_neuralnet=abc1(a,b,c, tol=0.1, method="neuralnet", hcorr = TRUE, transf = c("none"),
                       subset = NULL, kernel = "epanechnikov", numnet =10, sizenet = 10, lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit =500)

abcSSq1_neuralnet=abc1(a,b,c,tol=0.1,method="neuralnet",transf = c("none"),numnet =10, sizenet = 10)
abcSSq1_neuralnet$adj.values
abcSSq1_ridge=abc1(as.matrix(validationSSq1[2,-(1:4)]), SSq1_dataset[,2:4], as.matrix(SSq1_dataset[,-(1:4)]), tol=0.1, method="ridge", hcorr = TRUE, transf = c("none", "none", "log" ),
                      subset = NULL, kernel = "epanechnikov", numnet =
                        10, sizenet = 10, lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit =
                        500)

abcSSq1_ridge=abc1(a,b,c,tol=0.1,method="ridge",transf = c("log", "none", "log" ))
abcSSq1_ridge$adj.values
abcinferenceSSq1_hiersp=cbind(abcSSq1_rejection$unadj.values,abcSSq1_loclinear$adj.values,abcSSq1_neuralnet$adj.values,abcSSq1_ridge$adj.values)
write.csv(abcinferenceSSq1_hiersp,file="abcinferenceSSall_Panmixia_only_rejection__ridge_work.csv")

abcinferenceSSq1_Island=cbind(abcSSq1_rejection$unadj.values,abcSSq1_loclinear$adj.values,abcSSq1_neuralnet$adj.values)
write.csv(abcinferenceSSq1_Island,file="abcinferenceSSq1_Island_new.csv")
write.csv(abcSSq1_ridge$adj.values,file="abcinferenceSSall_stepping_ridge.csv")


### this is deep learning
deepregSSall_none=deepnetgression(as.matrix(validationSSall[1,-(1:4)]), as.matrix(SSall_dataset[,2:4]), as.matrix(SSall_dataset[,-(1:4)]),tol=0.1,iteration=10, hcorr =TRUE,
                             transf = c("log", "none", "log" ), kernel = "epanechnikov", encoder="none",
                             hidden = c(40,40), activationfun = "sigm", learningrate = 0.8,learningrate_scale=1,
                             momentum = 0.5,  output = "sigm",numepochs = 100, hidden_dropout = 0.1, visible_dropout = 0.1)
colMeans(deepregSSall_none$adj.values)
deepregSSall_sae=deepnetgression(as.matrix(validationSSall[1,-(1:4)]), as.matrix(SSall_dataset[,2:4]), as.matrix(SSall_dataset[,-(1:4)]),tol=0.1,iteration=10, hcorr =TRUE,
                                 transf = c("log", "none", "log" ), kernel = "epanechnikov", encoder="sae",
                                 hidden = c(40,40), activationfun = "sigm", learningrate = 0.8,learningrate_scale=1,
                                 momentum = 0.5,  output = "sigm",sae_output="sigm",numepochs = 100, hidden_dropout = 0.1, visible_dropout = 0.1)
colMeans(deepregSSall_sae$adj.values)
deepregSSall_db=deepnetgression(as.matrix(validationSSall[1,-(1:4)]), as.matrix(SSall_dataset[,2:4]), as.matrix(SSall_dataset[,-(1:4)]),tol=0.1,iteration=10, hcorr =TRUE,
                                 transf = c("log", "none", "log" ), kernel = "epanechnikov", encoder="db",
                                 hidden = c(40,40), activationfun = "sigm", learningrate = 0.8,learningrate_scale=1,
                                 momentum = 0.5,  output = "sigm",numepochs = 100, hidden_dropout = 0.1, visible_dropout = 0.1)

colMeans(deepregSSall_db$adj.values)

deepinferenceSSall=cbind(deepregSSall_none$adj.values,deepregSSall_sae$adj.values,deepregSSall_db$adj.values)

deepinferenceSSall_value1=cbind(deepregSSall_none$adj.values1,deepregSSall_sae$adj.values1,deepregSSall_db$adj.values1)

deepinferenceSSall_value2=cbind(deepregSSall_none$adj.values2,deepregSSall_sae$adj.values2,deepregSSall_db$adj.values2)

deepinferenceSSall_value3=cbind(deepregSSall_none$adj.values3,deepregSSall_sae$adj.values3,deepregSSall_db$adj.values3)

write.csv(deepinferenceSSall,file="deepinferenceSS_all_Panmixia.csv")
write.csv(deepinferenceSSall_value1,file="deepinferenceSS_all_value1_Hier_steppingstone.csv")
write.csv(deepinferenceSSall_value2,file="deepinferenceSS_all_value2_Hier_steppingstone.csv")
write.csv(deepinferenceSSall_value3,file="deepinferenceSS_all_value3_Hier_steppingstone.csv")


abline(v=0.0656827, col="blue")

plot(density(deepregSSall_none$adj.values[,1],weights=deepregSSall_none$weights/sum(deepregSSall_none$weights)), col = "red", main="m within")
lines(density(abcSSq1_rejection$unadj.values[,1]))
lines(density(abcSSq1_ridge$adj.values[,1]))
lines(density(deepregSSall_sae$adj.values[,1]))
lines(density(deepregSSall_db$adj.values[,1]))
lines(density(deepregSSall_sae$adj.values1[,1]))
lines(density(deepregSSall_sae$adj.values3[,1]))
plot(deepregSSall_db$unadj.values[,1],deepregSSall_db$adj.values[,1])

lines(density(abcSSq1_rejection$unadj.values[,1]),cor="black")
lines(density(abcSSq1_rejection$unadj.values[,1]),cor="black")

plot(density(abcSSq1_rejection$unadj.values[,1]))


library(reshape2)
library(ggplot2)

Scenarios1:   Panmixia

adjectvaluem_within_SSq1=read.excel()
adjectvaluem_within_SSq1=melt(adjectvaluem_within_SSq1)
p=ggplot(adjectvaluem_within_SSq1,aes(x=value, color=variable)) + geom_density(aes(x = value, colour = variable)) + labs(x = "m_within",y="Density") +labs(colour = "Methods") +labs(title = "Panmixia")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=0.0656827, linetype="solid", color = "black")


adjectvalues_between_SSq1=read.excel()
adjectvalues_between_SSq1=melt(adjectvalues_between_SSq1)
p=ggplot(adjectvalues_between_SSq1,aes(x=value, color=variable)) + geom_density() + labs(x = "m_between",y="Density") +labs(colour = "Methods") +labs(title = "Panmixia")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=0, linetype="solid", color = "blue")+theme(legend.background = element_rect(fill="transparent",size=0.5,  colour ="transparent"))+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                                                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

p=ggplot(adjectvaluem_within_SSq1,aes(x=value, color=variable)) + geom_density() + labs(x = "m_within",y="Density") +labs(colour = "Methods") +labs(title = "Panmixia")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=0.0656827, linetype="solid", color = "black")
p=ggplot(adjectvaluem_within_SSq1,aes(x=value, color=variable)) + geom_density() + labs(x = "m_within",y="Density") +labs(colour = "Methods") +labs(title = "Panmixia")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=0.0656827, linetype="solid", color = "blue")+theme(legend.background = element_rect(fill="transparent",size=0.5,  colour ="transparent"))+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                                                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))


adjectvaluem_NPOP_SSq1=read.excel()
adjectvaluem_NPOP_SSq1=melt(adjectvaluem_NPOP_SSq1)
p=ggplot(adjectvaluem_NPOP_SSq1,aes(x=value, color=variable)) + geom_density() + labs(x = "NPOPt",y="Density") +labs(colour = "Methods") +labs(title = "Panmixia")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=19964, linetype="solid", color = "blue")+theme(legend.background = element_rect(fill="transparent",size=0.5,  colour ="transparent"))+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                                                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))





Scenarios2:   Island
Islandadjectvalues_within_SSq1=read.excel()
Islandadjectvalues_within_SSq1=melt(Islandadjectvalues_within_SSq1)
p=ggplot(Islandadjectvalues_within_SSq1,aes(x=value, color=variable)) + geom_density() + labs(x = "m_within",y="Density") +labs(colour = "Methods") +labs(title = "Island")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=0.0609778, linetype="solid", color = "blue")+theme(legend.background = element_rect(fill="transparent",size=0.5,  colour ="transparent"))+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                                                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

Islandadjectvalues_between_SSq1=read.excel()
Islandadjectvalues_between_SSq1=melt(Islandadjectvalues_between_SSq1)
p=ggplot(Islandadjectvalues_between_SSq1,aes(x=value, color=variable)) + geom_density() + labs(x = "m_between",y="Density") +labs(colour = "Methods") +labs(title = "Island")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=0, linetype="solid", color = "blue")+theme(legend.background = element_rect(fill="transparent",size=0.5,  colour ="transparent"))+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                                                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

Islandadjectvalues_NPOPt_SSq1=read.excel()
Islandadjectvalues_NPOPt_SSq1=melt(Islandadjectvalues_NPOPt_SSq1)
p=ggplot(Islandadjectvalues_NPOPt_SSq1,aes(x=value, color=variable)) + geom_density() + labs(x = "NPOPt",y="Density") +labs(colour = "Methods") +labs(title = "Island")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=17360, linetype="solid", color = "blue")+theme(legend.background = element_rect(fill="transparent",size=0.5,  colour ="transparent"))+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                                                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))


Scenarios3:   Hier_Hier_Island
Hier_Islandadjectvalues_within_SSq1=read.excel()
Hier_Islandadjectvalues_within_SSq1=melt(Hier_Islandadjectvalues_within_SSq1)
p=ggplot(Hier_Islandadjectvalues_within_SSq1,aes(x=value, color=variable)) + geom_density() + labs(x = "m_within",y="Density") +labs(colour = "Methods") +labs(title = "Hier_Island")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=0.0612868, linetype="solid", color = "blue")+theme(legend.background = element_rect(fill="transparent",size=0.5,  colour ="transparent"))+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                                                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

Hier_Islandadjectvalues_between_SSq1=read.excel()
Hier_Islandadjectvalues_between_SSq1=melt(Hier_Islandadjectvalues_between_SSq1)
p=ggplot(Hier_Islandadjectvalues_between_SSq1,aes(x=value, color=variable)) + geom_density() + labs(x = "m_between",y="Density") +labs(colour = "Methods") +labs(title = "Hier_Island")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=0.0049995, linetype="solid", color = "blue")+theme(legend.background = element_rect(fill="transparent",size=0.5,  colour ="transparent"))+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                                                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

Hier_Islandadjectvalues_NPOPt_SSq1=read.excel()
Hier_Islandadjectvalues_NPOPt_SSq1=melt(Hier_Islandadjectvalues_NPOPt_SSq1)
p=ggplot(Hier_Islandadjectvalues_NPOPt_SSq1,aes(x=value, color=variable)) + geom_density() + labs(x = "NPOPt",y="Density") +labs(colour = "Methods") +labs(title = "Hier_Hier_Island")+ theme(legend.position=c(0.9, 0.8))+geom_vline(xintercept=16666, linetype="solid", color = "blue")+theme(legend.background = element_rect(fill="transparent",size=0.5,  colour ="transparent"))+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                                                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))




### model choice and 
prod(table(as.vector(SSq0_dataset[,1])))
table(index)
postpr_deepSSq1=postpr_deep(target, index, trainSS_all, tol=0.1, subset=NULL, method="neuralnet", corr=TRUE, kernel="epanechnikov",
                         numnet = 10, sizenet = 5, lambda = c(0.0001,0.001,0.01), trace = TRUE, maxit = 500, hidden=c(30,30), activationfun="sigm", learningrate=0.5,
                         momentum=0.5, learningrate_scale=1, output="sigm", sae_output="sigm", numepochs=100, batchsize=100, hidden_dropout=0.1 , visible_dropout=0.1,threshold = 0.01,
                         stepmax = 1e+05, rep = 1, startweights = NULL,
                         learningrate.limit = NULL, learningrate.factor = list(minus = 0.5,plus = 1.2),  lifesign = "none",
                         lifesign.step = 1000, algorithm = "rprop+", err.fct = "sse",
                         act.fct = "logistic", linear.output = TRUE, exclude = NULL,
                         constant.weights = NULL, likelihood = FALSE)


summary.postpr(postpr_deepSSq1)

cv4postpr_deep1=cv4postpr_deep(index, trainSS_all, postpr.out = NULL, nval=10, tols=0.1,
                               method="rejection", subset = NULL, kernel = "epanechnikov",numnet = 5, sizenet = 5, lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit = 50, hidden=c(30,30), activationfun="sigm", learningrate=0.5,
                               momentum=0.5, batchsize=50,learningrate_scale=1, output="sigm", sae_output="sigm", numepochs=100,  hidden_dropout=0.1 , visible_dropout=0.1,threshold = 0.01,
                               stepmax = 1e+05, rep = 1, startweights = NULL,
                               learningrate.limit = NULL, learningrate.factor = list(minus = 0.5,plus = 1.2),  lifesign = "none",
                               lifesign.step = 1000, algorithm = "rprop+", err.fct = "sse",
                               act.fct = "logistic", linear.output = TRUE, exclude = NULL,
                               constant.weights = NULL, likelihood = FALSE)

summary(cv4postpr_deep1)

plot(cv4postpr_deep0)


cv.modsel <- cv4postpr_deeep(models, stat.3pops.sim, nval=5, tol=.01, method="mnlogistic")
modsel.ha <- postpr_deep(stat.voight["hausa",], models, stat.3pops.sim, tol=.05, method="mnlogistic")
modsel.it <- postpr_deep(stat.voight["italian",], models, stat.3pops.sim, tol=.05, method="mnlogistic")
modsel.ch <- postpr_deep(stat.voight["chinese",], models, stat.3pops.sim, tol=.05, method="mnlogistic")
summary(modsel.ha)