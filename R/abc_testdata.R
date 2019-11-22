library(abc.data)
data(human)    stat.voight stat.3pops.sim models  
View(par.italy.sim)
prior_param=subset(stat.3pops.sim,  subset=models=="bott")

View(prior_param)
library(randomForest)
data(musigma2)
data(ppc)
