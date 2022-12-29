 
 
 
############################################################################
############################################################################
##############Sequential boostrap####################
############################################################################
############################################################################
library(nleqslv)  


FindRoot = function(x,gamma1,gamma2,alpha) {
return(pnorm(x-gamma2)-exp(2*gamma2*x)*pnorm(-x-gamma2)-1+alpha)
}
MAX<-100
max<-100
Qn1<-rep(0,max)
Times1<-rep(0,max)
n1<-50
n2<-50
n<-100
sigma1<-1
sigma2<-1
mu1<-0.5
mu2<--0.5
mu<--0.5
c<-0.5
Z_alpha<-nleqslv(1, FindRoot, control=list(btol=.000001),jacobian=TRUE,method="Newton",gamma1=gamma1,gamma2=2*mu-2*c,alpha= 0.05)$x
hatsigma1<-sigma1
hatsigma2<-sigma2

for(k in 1:max){
    Tn1<-rep(0,MAX)
    Tn2<-rep(0,MAX)
    Tn3<-rep(0,MAX)
    Tn4<-rep(0,MAX)
     X<-rnorm(n,mu,sigma1)
     hatsigma<-sd(X)
     x<-c-X
    y<-X-c
    hatsigma1<-sd(x)
    hatsigma2<-sd(y)
    for(j in 1:MAX){
    if(((sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*sqrt((hatsigma1^2/n+hatsigma2^2/n))))<=0){
   xx<-sample(x,n,replace=TRUE)
   yy<-sample(y,n,replace=TRUE) 
    #hatsigma1<-sd(xx)
     #  hatsigma2<-sd(yy)
    Tn1[j]<-(mean(xx)-mean(yy)-(2*c-2*mean(X)))
    Tn2[j]<-mean(xx)-mean(yy)
    }else{
   xx<-sample(x,n,replace=TRUE)
   yy<-sample(y,n,replace=TRUE) 
    # hatsigma1<-sd(xx)
     #  hatsigma2<-sd(yy)
    Tn3[j]<-(mean(yy)-mean(xx)-(2*mean(X)-2*c))
    Tn4[j]<-mean(yy)-mean(xx)}      
   }
   Qn1[k]<-sum(Tn2+Tn4)/(MAX)+sum(Tn1+Tn3)/(sqrt(MAX)*sqrt((hatsigma1^2/n+hatsigma2^2/n)))
   Times1[k]<-(abs(Qn1[k])>Z_alpha)+0  
 }

densityf<-function(y,alpha,c){
exp(-(y^2-2*alpha*(abs(y-c)-abs(c))+alpha^2)/2)/sqrt(2*pi)-
alpha*exp(2*alpha*abs(y-c))*(1-pnorm(abs(c)+abs(y-c)+alpha,0,1))
}   
YY<-seq(-5,5,length=100)
Alpha<-2*mu-2*c
res<-array(0,dim=c(length(YY),length(Alpha)))
for(j in 1:length(YY)){
for(k in 1:length(Alpha))
res[j,k]<-densityf(YY[j],Alpha[k],0)
}
plot(YY,res,type="l",ylab="density",xlab="",col="blue")
lines(density(Qn1),ylab="density",col="red",xlab="")
mean(Times1)

############################################################################
##################Sequential data/streaming data########################
############################################################################

library(nleqslv)  
FindRoot = function(x,gamma1,gamma2,alpha) {
return(pnorm(x-gamma2)-exp(2*gamma2*x)*pnorm(-x-gamma2)-1+alpha)
}
MAX<-10
max<-500
Qn1<-rep(0,max)
Times1<-rep(0,max)
Mean<-matrix(data=0,max, MAX)
n<-100
mu<-1
c<-1.1
sigma1<-0.68
Z_alpha<-nleqslv(1, FindRoot, control=list(btol=.000001),jacobian=TRUE,method="Newton",gamma1=gamma1,gamma2=-c,alpha= 0.05)$x
hatsigma1<-sigma1
rm<-matrix(data=0,max,MAX)
for(k in 1:max){
    Tn1<-rep(0,MAX)
    Tn2<-rep(0,MAX)
    Tn3<-rep(0,MAX)
    Tn4<-rep(0,MAX)
    print(k)
   for(j in 1:MAX){
       
    X<-rnorm(n,mu,sigma1) 
    mean<-mean(X)
    rm[k,j]=mean
    if(j<2){hatsigma1=0.3
    h<-rnorm(1,0,1) }
    else{hatsigma1<-sd(X)
    h=-1}

    if(((sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*sqrt((hatsigma1^2/n))))<=0 & h<0){
    Tn1[j]<-mean(X)-c
    Tn2[j]<-mean(X)
    }else{
    Tn3[j]<-(c-mean(X))
    Tn4[j]<--mean(X)}      
   }
   Qn1[k]<-sum(Tn2+Tn4)/(MAX)+sum(Tn1+Tn3)/(sqrt(MAX)*sqrt((hatsigma1^2/n)))
   Times1[k]<-(abs(Qn1[k])>Z_alpha)+0  
 }
densityf<-function(y,alpha,c){
exp(-(y^2-2*alpha*(abs(y-c)-abs(c))+alpha^2)/2)/sqrt(2*pi)-
alpha*exp(2*alpha*abs(y-c))*(1-pnorm(abs(c)+abs(y-c)+alpha,0,1))
}   
YY<-seq(-5,5,length=1000)
Alpha<--c
res<-array(0,dim=c(length(YY),length(Alpha)))
for(j in 1:length(YY)){
for(k in 1:length(Alpha))
res[j,k]<-densityf(YY[j],Alpha[k],0)
}
plot(YY,res,type="l",ylab="density",xlab="",col="blue",panel.first=grid(10,10,col="gray70"), xlim=c(-20,20))#
lines(density(Qn1),ylab="density",col="red",xlab="")
mean(Times1)
sd(apply(rm,1,sd))
############################################################################
##################ÏÔÖøÐÔ¼ìÑé########################
############################################################################

library(nleqslv)  
FindRoot = function(x,gamma1,gamma2,alpha) {
    return(pnorm(x-gamma2)-exp(2*gamma2*x)*pnorm(-x-gamma2)-1+alpha)
}
MAX<-500
max<-1000
Qn1<-rep(0,max)
Times1<-rep(0,max)
Mean<-matrix(data=0,max, MAX)
#n1<-50
#n2<-50
n<-100
#sigma1<-1
#sigma2<-1
#mu1<-0.5
#mu2<--0.5
mu<--0.02
c<-0
Z_alpha<-nleqslv(1, FindRoot, control=list(btol=.000001),jacobian=TRUE,method="Newton",gamma1=gamma1,gamma2=0,alpha= 0.05)$x
hatsigma1<-sigma1
hatsigma2<-sigma2

for(k in 1:max){
    Tn1<-rep(0,MAX)
    Tn2<-rep(0,MAX)
    Tn3<-rep(0,MAX)
    Tn4<-rep(0,MAX)
    #hatsigma1<-sigma1
    
    print(k)
    
    for(j in 1:MAX){
        #if(j < 500){
        #X<-rnorm(n,mu,sigma1)}
        #   else{X<-rnorm(n,mu-1,sigma1)}
        
        X<-rnorm(n,mu,sigma1)   
        hatsigma1<-sd(X)
        #Mean[k,j]<-mean(X)
        #row_mean = mean(Mean[k,1:j])
        
        
        
        if(((sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*sqrt((hatsigma1^2/n))))<=0){
            
            Tn1[j]<-mean(X)-c
            Tn2[j]<-mean(X)
        }else{
            
            Tn3[j]<-(c-mean(X))
            Tn4[j]<--mean(X)}      
    }
    Qn1[k]<-sum(Tn2+Tn4)/(MAX)+sum(Tn1+Tn3)/(sqrt(MAX)*sqrt((hatsigma1^2/n)))
    Times1[k]<-(abs(Qn1[k])>Z_alpha)+0  
}

densityf<-function(y,alpha,c){
    exp(-(y^2-2*alpha*(abs(y-c)-abs(c))+alpha^2)/2)/sqrt(2*pi)-
        alpha*exp(2*alpha*abs(y-c))*(1-pnorm(abs(c)+abs(y-c)+alpha,0,1))
}   
YY<-seq(-5,5,length=1000)
Alpha<-0
res0<-array(0,dim=c(length(YY),length(Alpha)))
for(j in 1:length(YY)){
    for(k in 1:length(Alpha))
        res0[j,k]<-densityf(YY[j],0,0)
}

res1<-array(0,dim=c(length(YY),length(Alpha)))
for(j in 1:length(YY)){
    for(k in 1:length(Alpha))
        res1[j,k]<-densityf(YY[j],1,0)
}

res2<-array(0,dim=c(length(YY),length(Alpha)))
for(j in 1:length(YY)){
    for(k in 1:length(Alpha))
        res2[j,k]<-densityf(YY[j],2,0)
}

res11<-array(0,dim=c(length(YY),length(Alpha)))
for(j in 1:length(YY)){
    for(k in 1:length(Alpha))
        res11[j,k]<-densityf(YY[j],-1,0)
}

res22<-array(0,dim=c(length(YY),length(Alpha)))
for(j in 1:length(YY)){
    for(k in 1:length(Alpha))
        res22[j,k]<-densityf(YY[j],-2,0)
}

#ress<-array(0,dim=c(length(YY),length(Alpha)))
#for(j in 1:length(YY)){
#for(k in 1:length(Alpha))
#        ress[j,k]<-densityf(YY[j],-Alpha[k],0)
#}
TeX("$\\kappa$")

plot(YY,res0,type="l",ylab="density",xlab="",col="skyblue1",panel.first=grid(10,10,col="gray70"),ylim=c(0,2))#xlim=c(-15,15)
#lines(YY,ress,type="l",ylab="density",xlab="",col="red")
legend("topleft",inset=0.05,  col=c("red","pink","skyblue1","gray","black"),c(TeX("$\\kappa$=2"),TeX("$\\kappa$=1"),TeX("$\\kappa$=0"),TeX("$\\kappa$=-1"),TeX("$\\kappa$=-2")),lty=c(1,1),bty="n")
lines(YY,res1,ylab="density",col="pink",xlab="")
lines(YY,res2,ylab="density",col="red",xlab="")
lines(YY,res11,ylab="density",col="gray",xlab="")
lines(YY,res22,ylab="density",col="black",xlab="")
mean(Times1)
