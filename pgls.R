library(caper)
library(ggplot2)
library(ggrepel)
library(treeio)

ppp <- function(cdat,tag){
  
  mod <- pgls(V3 ~ V2,cdat,lambda = "ML")
  
  sum <- summary(mod)
  
  write.table(sum$coefficients,file=paste(args[3],tag,"pgls.txt",sep=".") ,sep="\t",quote = F)
  
  preds <- data.frame( fit = as.matrix(mod$fitted),r = mod$x[,2])
  
  pdf(paste(args[3],tag,"pgls.pdf",sep="."),width = 6 ,height = 4)
  da <- cdat$data
  gp <- ggplot(da,aes(x=V2,y=V3,label = rownames(da) )) + geom_point() +
    geom_text_repel() + 
    geom_line(aes(x=r,y=fit), data = preds,colour = "red") +
    annotate("text",x=(range(da$V2)[2] - range(da$V2)[1])/8,y=0.9  ,label= paste("p-value = ",round(sum$coefficients[2,4],digits = 3),sep="")) +
    annotate("text",x=(range(da$V2)[2] - range(da$V2)[1])/8,y=0.8  ,label= paste("lambda = ",round(sum$param[2],digit = 3),sep="")) +
    theme_bw()
  print(gp)
  par(mfrow=c(2,2))
  plot(mod)
  dev.off()
}
## trait   tree   prefix
args <- commandArgs(trailingOnly = T)
#args <- c("akr","pla")

da <- read.table(args[1])
da$V3 <- as.numeric(da$V3)
da <- da[!da$V1 == "Talpa_occidentalis",]
da <- na.omit(da)
da <- da[da$V2 != 0,]


t <- read.codeml_mlc(args[2])
phy <- as.phylo(t)


cdat <- comparative.data(data <- da,phy = phy,names.col = "V1")
#rownames(da) <- da$V1


mod <- pgls(V3  ~ V2,cdat,lambda = "ML")
res <- residuals(mod, phylo = TRUE)
res<- res/sqrt(var(res))[1]
#rownames(res)[(abs(res)>3)]
cdat <- cdat[abs(res) < 3,]


ppp(cdat,"raw")

if(F){
  library(Routliers)
  ol <- outliers_mad(da$V2)
  cdats <- subset(cdat, subset=!cdat$data$V2 %in% ol$outliers)
  ppp(cdats,"rmout")
  
  
  cdate <- subset(cdat, subset= cdat$data$V3 != 0 )
  ppp(cdate,"ext")
  
}