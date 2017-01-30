library(ggplot2)
library(data.table)
library(splines)

#sparse spline design
splinedegree=3 #Cubic spline
Krow=100
Ncs=1000
cutsites = sort(runif(n=Ncs, min=1, max=10))
dx = 1.01*(max(cutsites)-min(cutsites))/(Krow-splinedegree)
t = min(cutsites) - dx*0.01 + dx * seq(-splinedegree,Krow-splinedegree+3)
Bsp = spline.des(cutsites, knots = t, outer.ok = T, sparse=T)$design
X = rbind(cbind(Bsp/2,Bsp/2),bdiag(Bsp,Bsp),bdiag(Bsp,Bsp))
Sm2=Diagonal(5*Ncs,x=c(rnorm(Ncs)^2*10,rnorm(2*Ncs)^2*5,rnorm(2*Ncs)^2))
ggplot(as.data.table(melt(as.matrix(X))))+geom_raster(aes(Var1,Var2,fill=value))
#
liota=cos(2*cutsites)+cutsites/10
lrho=(2+sin(cutsites))*cutsites^2/20
liota=liota-mean(liota)
lrho=lrho-mean(lrho)
y=matrix(c(10+(liota+lrho)/2,5+liota,5+lrho,1+liota,1+lrho)+rnorm(5*Ncs,sd=0.1),ncol=1)
ggplot(data.table(x=rep(cutsites,5),cat=rep(1:5,each=Ncs),y=y[,1]))+geom_point(aes(x,y))+facet_wrap(~cat)
data.table(x=rep(cutsites,5),cat=rep(1:5,each=Ncs),y=y[,1])[,mean(y),by=cat]
#
diags = list(rep(1,Krow), rep(-2,Krow))
D = bandSparse(Krow-2, Krow, k=0:2, diagonals=c(diags, diags[1]))
D = bdiag(D,D)
ggplot(as.data.table(melt(as.matrix(D))))+geom_raster(aes(Var1,Var2,fill=value))
#
A=crossprod(crossprod(Sm2,X),X)+Krow*0.1*crossprod(D)
cholA=Cholesky(A)
ggplot(as.data.table(melt(as.matrix(A))))+geom_raster(aes(Var1,Var2,fill=value))#A
ggplot(as.data.table(melt(as.matrix(A)==0)))+geom_raster(aes(Var1,Var2,fill=value))#A==0
ggplot(as.data.table(melt(as.matrix(crossprod(expand(cholA)[[1]])==0))))+geom_raster(aes(Var1,Var2,fill=value))#permuted A
ggplot(as.data.table(melt(as.matrix(expand(cholA)[[2]]==0))))+geom_raster(aes(Var1,Var2,fill=value))#factor
ggplot(as.data.table(melt(as.matrix(expand(cholA)[[1]]==0))))+geom_raster(aes(Var1,Var2,fill=value))#permutation

#
W=cbind(c(rep(0,3*Ncs),rep(1,Ncs),rep(0,Ncs)), c(rep(0,4*Ncs),rep(1,Ncs)))
ggplot(as.data.table(melt(as.matrix(W==0))))+geom_raster(aes(Var1,Var2,fill=value))
#
U=cbind(c(rep(1,Ncs),rep(0,4*Ncs)),c(rep(0,Ncs),rep(1,2*Ncs),rep(0,2*Ncs)),c(rep(0,3*Ncs),rep(1,2*Ncs)))
ggplot(as.data.table(melt(as.matrix(U==0))))+geom_raster(aes(Var1,Var2,fill=value))
#
Gamma=t(X)%*%W
Gamma=Diagonal(2*Krow) - Gamma%*%solve(t(Gamma)%*%solve(cholA,Gamma), t(solve(cholA,Gamma)))
all(abs(Gamma%*%Gamma-Gamma)<1e-5)
#
beta_y=solve(cholA, Gamma%*%t(X)%*%Sm2%*%y)
beta_U=solve(cholA, Gamma%*%t(X)%*%Sm2%*%U)
#
e=solve(t(U)%*%Sm2%*%(U-X%*%beta_U),t(U)%*%Sm2%*%(y-X%*%beta_y))
#beta = beta_y - beta_U%*%e
beta=solve(cholA, Gamma%*%t(X)%*%Sm2%*%(y-U%*%e))
t(W)%*%X%*%beta
#
mu=t(X)%*%W
mu=solve(t(mu)%*%solve(cholA,mu), t(mu)%*%solve(cholA, t(X)%*%Sm2%*%(y-U%*%e)))
#
lambda=(Krow-2)/(Krow^2*crossprod(D%*%beta))

###
result = data.table(x=rep(cutsites,5), cat=rep(1:5,each=Ncs), y=y[,1],
                    mu=as.matrix(X%*%beta)[,1],
                    mid=as.matrix(U%*%e)[,1],
                    muy=as.matrix(X%*%beta_y)[,1],
                    muU=as.matrix(X%*%beta_U%*%e)[,1])
ggplot(result)+
  geom_point(aes(x,y))+facet_wrap(~cat)+geom_line(aes(x,mu+mid),colour="green")
ggplot(result)+
  geom_point(aes(x,y))+facet_wrap(~cat)+geom_line(aes(x,mu,colour="mu"))+
  geom_line(aes(x,muU,colour="mu U"))+
  geom_line(aes(x,muy,colour="mu y"))+
  geom_line(aes(x,mid,colour="mid"))


