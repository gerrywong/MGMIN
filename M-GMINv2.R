###INPUT:
###df: a data.frame consisting of 
###(1st column) beta values for a given sample and 
###(column name: type ) corresponding probe types (1=type1, 2=type2).
###sampleID: TCGA sample ID.

require(mixtools)

MGMIN <- function(df, sampleID){
flag <- any(is.na(df[,1]))
if(flag){
idx <- !is.na(df[,1])
df2 <- df[idx,]
}
else{df2<-df}

## transit beta-values to M-values

obeta <- as.numeric(df2[,1])

index1<-which(df2$type==1)
index2<-which(df2$type==2)
obeta1<-obeta[index1]
obeta2<-obeta[index2]

if(min(obeta1)==0){
	obeta1[obeta1==0] <- min(setdiff(obeta1,0))
}
if(min(obeta2)==0){
	obeta2[obeta2==0] <- min(setdiff(obeta2,0))
}
if(max(obeta1)==1){
	obeta1[obeta1==1] <- max(setdiff(obeta1,1))
}
if(max(obeta2)==1){
	obeta2[obeta2==1] <- max(setdiff(obeta2,1))
}

mvaluev1<-log2(obeta1/(1-obeta1))
mvaluev2<-log2(obeta2/(1-obeta2))
em1<-normalmixEM(mvaluev1,mu=c(-4,-1,4),sigma=c(1,1,1))
em2<-normalmixEM(mvaluev2,mu=c(-4,0,4),sigma=c(1,1,1))
#plot(em1,whichplots=2)
#plot(em2,whichplots=2)

#added 0729, mvaluev1 fit four mixed normal distributions
flag1<-(em1$mu[3]-em1$sigma[3])>(em2$mu[3]-em2$sigma[3])
if (flag1){
	em1<-em1
}
else{
	em1<-normalmixEM(mvaluev1,mu=c(-4,0,3,5),sigma=c(1,1,1,1))
}


type1class<-apply(em1$posterior, 1, which.max)
type1.v <- rep(2,length(mvaluev1));
type1.v[which(mvaluev1 < max(mvaluev1[type1class==1]))] <- 1;
type1.v[which(mvaluev1 > min(mvaluev1[type1class==3]))] <- 3;
type1class <- type1.v
#type1thres <-c(mean(c(max(mvaluev1[type1class==1]),min(mvaluev1[type1class==2]))),mean(c(max(mvaluev1[type1class==2]),min(mvaluev1[type1class==3]))))
type2class<-apply(em2$posterior, 1, which.max)
type2.v <- rep(2,length(mvaluev2));
type2.v[which(mvaluev2 < max(mvaluev2[type2class==1]))] <- 1;
type2.v[which(mvaluev2 > min(mvaluev2[type2class==3]))] <- 3;
type2class <- type2.v
#type2thres <-c(mean(c(max(mvaluev2[type2class==1]),min(mvaluev2[type2class==2]))),mean(c(max(mvaluev2[type2class==2]),min(mvaluev2[type2class==3]))))


type1mu <- em1$mu
type2mu <- em2$mu

nmvaluev2 <- mvaluev2

type2U.idx <- which(type2class==1)
selUR.idx <- type2U.idx[which(mvaluev2[type2U.idx] > type2mu[1])]
selUL.idx <- type2U.idx[which(mvaluev2[type2U.idx] < type2mu[1])]
# p.v <- pnorm(mvaluev2[selUR.idx], em2$mu[1], em2$sigma[1], lower.tail=FALSE)
# q.v <- qnorm(p.v, em1$mu[1], em1$sigma[1], lower.tail=FALSE)
# nmvaluev2[selUR.idx] <- q.v
p.v <- pnorm(mvaluev2[selUL.idx], em2$mu[1], em2$sigma[1], lower.tail=TRUE)
q.v <- qnorm(p.v, em1$mu[1], em1$sigma[1], lower.tail=TRUE)
nmvaluev2[selUL.idx] <- q.v

type2M.idx <- which(type2class==3)
selMR.idx <- type2M.idx[which(mvaluev2[type2M.idx] > type2mu[3])]
selML.idx <- type2M.idx[which(mvaluev2[type2M.idx] < type2mu[3])]
p.v <- pnorm(mvaluev2[selMR.idx], em2$mu[3], em2$sigma[3], lower.tail=FALSE)
if(flag1){
	q.v <- qnorm(p.v, em1$mu[3], em1$sigma[3], lower.tail=FALSE)

}
else{
	q.v <- qnorm(p.v, em1$mu[4], em1$sigma[4], lower.tail=FALSE)

}
nmvaluev2[selMR.idx] <- q.v
#p.v <- pnorm(mvaluev2[selML.idx], em2$mu[3], em2$sigma[3], lower.tail = TRUE)
#q.v <- qnorm(p.v, em1$mu[3], em1$sigma[3], lower.tail = TRUE)
#nmvaluev2[selML.idx] <- q.v

type2H.idx <- c(selUR.idx, which(type2class==2), selML.idx)

minH <- min(mvaluev2[type2H.idx])
#maxH <- max(mvaluev2[type2H.idx])
maxH <- min(mvaluev2[selMR.idx])
delH<- maxH-minH
delUH <- min(mvaluev2[type2H.idx]) - max(mvaluev2[selUL.idx])
#delHM <- min(mvaluev2[selMR.idx]) - max(mvaluev2[type2H.idx])
delHM <- min(mvaluev2[selMR.idx]) - max(mvaluev2[selML.idx])

nmaxH<- min(nmvaluev2[selMR.idx]) - delHM
nminH <- max(nmvaluev2[selUL.idx]) + delUH
ndelH <- nmaxH - nminH

hf <- ndelH/delH

nmvaluev2[type2H.idx] <- nminH + hf*(mvaluev2[type2H.idx] - minH)

####plot typeI
if(flag1){
	n=3
}
else{
	n=4
}
prob1 = colSums(em1$posterior)/length(mvaluev1)
tmpL.v <- as.vector(rmultinom(1:n,length(mvaluev1),prob=prob1));
tmpB.v <- vector();
for(l in 1:n){
tmpB.v <- c(tmpB.v,rnorm(tmpL.v[l],em1$mu[l],em1$sigma[l]));
}

pdf(paste0("Type1fit-",sampleID,"-",
	paste(round(em1$mu,3),collapse=" "),".pdf"),
	pointsize=10.46094, family="Times", width=4, height=3);
dv1.o <- density(mvaluev1)
d.o <- density(tmpB.v);
ymax <- max(d.o$y, dv1.o$y)
plot(density(mvaluev1),main="",ylim=c(0,ymax));
points(d.o$x,d.o$y,col="green",type="l")
legend("top",legend=c("typeI-obs","typeI-fit"),bty="n",fill=c("black","green"));
plot(em1,whichplots=2)
dev.off();

##plot type2
prob2 = colSums(em2$posterior)/length(mvaluev2)
tmpL.v <- as.vector(rmultinom(1:3,length(mvaluev2),prob=prob2));
tmpB2.v <- vector();
for(l in 1:3){
tmpB2.v <- c(tmpB2.v,rnorm(tmpL.v[l],em2$mu[l],em2$sigma[l]));
}
pdf(paste0("Type2fit-",sampleID,"-",
	paste(round(em2$mu,3),collapse=" "),".pdf"),
	pointsize=10.46094, family="Times", width=4, height=3);
dv2.o <-density(mvaluev2)
d.o <- density(tmpB2.v);
ymax <- max(d.o$y, dv2.o$y)
plot(density(mvaluev2),main="",ylim=c(0,ymax));
points(d.o$x,d.o$y,col="green",type="l")
legend("top",legend=c("typeII-obs","typeII-fit"),bty="n",fill=c("black","green"));
plot(em2,whichplots=2)
dev.off();

mvalue<-df2[,1]
mvalue[index1] <- mvaluev1
mvalue[index2] <- nmvaluev2
nbeta <- 2^mvalue/(1+2^mvalue)
df2[,1]<-nbeta

####2 plot
d1.o <- density(nbeta[index1])
d2.o <- density(obeta[index2])
d2n.o<- density(nbeta[index2])
ymax <- max(d1.o$y, d2.o$y, d2n.o$y)
pdf(paste0("fitting_curves_after_normalization_of_",sampleID, ".pdf"), 
	pointsize=10.46094, family="Times", width=4.724, height=3.7792)
plot(density(obeta[index2]), type="l", ylim=c(0,ymax), main="")
points(d1.o$x, d1.o$y, col="blue", type="l")
points(d2n.o$x, d2n.o$y, col="red", type="l")
legend(x=0.3,y=ymax,legend=c("typeI","typeII","typeII-MGMIN"),fill=c("blue","black","red"), bty="n")
dev.off()

if(flag){
	df[idx,1] <- nbeta
}
else{df<-df2}

print(paste0("mu of type1 M-values:", paste(round(type1mu,3),collapse=' ')))
print(paste0("sigma of type1 M-values:", paste(round(em1$sigma,3),collapse=' ')))
print(em1$lambda)
print(paste0("mu of type2 M-values:", paste(round(type2mu,3),collapse=' ')))
print(paste("sigma of type2 M-values:",paste(round(em2$sigma,3),collapse=" ")))
print(em2$lambda)
print(hf)
print(paste0("Finish normalization of ", sampleID,"."))
return(df)
}
