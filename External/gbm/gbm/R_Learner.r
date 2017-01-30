library(randomForest) ;
data.all <- read.table("R_LrnData",header=T) ;
data <- subset(data.all , select = -Label) ;
labels <- data.all$Label ;
labels.factor <- as.factor(labels) ;
mdl <- randomForest(data,labels.factor,sampsize=c(836,418),ntree=1000)
save(mdl,file="R_Model") ;
sink("R_Trees") ;
print(mdl$forest$ntree) ;
for (i in 1:mdl$forest$ntree) {
	tree <- getTree(mdl,i) ;
	print(attr(tree,"dim")[[1]]) ;
	print(tree) ;
}
sink() ;
