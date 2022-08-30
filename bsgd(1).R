####################################################################
#基于bootstrap的带有扰动的随机梯度下降�?    ##################### 
#################################################################
result<-function(){
  W<-rnorm(150,1,1)
  Wb<-matrix(data = 0,150,30)
  for (b in 1:30) {
    Wb[,b]<-sample(W,150,replace = TRUE)
  }
  #生成bootstrap数据
  btsgd<-function(x,y,error,maxiter,step=0.001,B,b,Wb)#B是bootsrap总次�?,b是第几次bootstrap
  { 
    m<-nrow(x)
    
    n<-ncol(x) 
    theta<-matrix(rep(0,n),n,1)  #ktheta初始值都设置�?0
    iter<-0   #迭代次数
    k<-0  #第k个样�?
    newerror<-1 
    
    while(iter<maxiter|newerror>error){
      
      iter<-iter+1
      k<-sample(1:m, 1, replace = TRUE)
      xk<-x[k,,drop=FALSE]
      yk<-y[k,,drop=FALSE]
      hk<-xk%*%theta
      des<-t((hk-yk)%*%xk)
      new_theta<-theta-Wb[k,b]*step*des
      #newerror<-t(new_theta-theta)%*%(new_theta-theta)
      theta<-new_theta #/(iter)+(theta*(iter-1))/(iter)   
      
      
    }
    costfunction<-(t(x%*%theta-y)%*%(x%*%theta-y))/m
    result<-c(theta,iter,costfunction)
    #names(result)<-c('系数','系数','系数','系数','系数','迭代次数','误差')
    result 
  }
  #########################
  #生成数据  100个样�?#####
  #########################
  
  x<-array(runif(500,0,30),c(100,5))#协变�?
  thet<-c(1,2,3,4,10) #要估计的参数
  y<-x%*%thet+rnorm(100)#模型
  #btsgd(x,y,error=1,maxiter=1000,step=0.001,B=500,b=3,Wb)
  
  
  
  
  ############################
  ##   500次bootsrap      ###
  ############################
  
  r1<-matrix(data=0,5,30)
  #r2<-matrix(data=0,1,500)
  #r3<-matrix(data=0,1,30)
  for (b in 1:30) {
    r<-btsgd(x,y,error=1,maxiter=150,step=0.001,B=30,b,Wb)
    r1[,b]<-r[1:5]
    #r2[,b]<-do.call(cbind,r1[,b])
    #r3[b]<-(t(r1[,b]-thet))%*%(r1[,b]-thet)
  }
  
  #plot(density(r3),type="l",ylab="密度函数",xlab="误差",col="red",xlim=c(-5,5))
  #mean(r3)
  #r1是我们需要的数据,但是这里r1是向�?
  #r1
  #dim(r1)
  return(r1)
}


############################################################################
############################################################################
##################prove rhe spike distribution Theorem 1####################
##################prove rhe spike distribution Theorem 1####################
############################################################################
############################################################################
library(nleqslv)  
FindRoot = function(x,gamma1,gamma2,alpha) {
  return(pnorm(x-gamma2)-exp(2*gamma2*x)*pnorm(-x-gamma2)-1+alpha)
}
MAX<-150
max<-200
Qn1<-rep(0,max)
Times1<-rep(0,max)
n1<-50
n2<-50
sigma1<-1
sigma2<-1
mu1<-0.5
mu2<--0.5
Z_alpha<-nleqslv(1, FindRoot, control=list(btol=.000001),jacobian=TRUE,method="Newton",gamma1=gamma1,gamma2=mu2-mu1,alpha= 0.05)$x
hatsigma1<-sigma1
hatsigma2<-sigma2



for(k in 1:max){
  Tn1<-rep(0,MAX)
  Tn2<-rep(0,MAX)
  Tn3<-rep(0,MAX)
  Tn4<-rep(0,MAX)
  print(k)
  hatsigma1<-1
  for(j in 1:MAX){
    
    if(((sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*hatsigma1))<=0){
     # xx<-rnorm(n1,mu1,sigma1)
      #yy<-rnorm(n2,mu2,sigma2)
      #hatsigma2<-sd(yy)
      
      
      #r1=rchisq(50,3)
      #r1=rnorm(50,10,1)
      xx=r1[,5]
      yy=11
      Tn1[j]<-yy-mean(xx)-(mu1-mu2)
      Tn2[j]<-yy-mean(xx)
      hatsigma1<-(sd(xx)/sqrt(30))
    }else{
      #xx<-rnorm(n1,mu1,sigma1)
      #yy<-rnorm(n2,mu2,sigma2)
      #hatsigma2<-sd(yy)
      # hatsigma1<-0.2
      
      r2=result()
      #r2=rchisq(50,3)
      #r2=rnorm(50,10,1)
      xxx=r2[,5]
      Tn3[j]<-mean(xxx)-yy-(mu2-mu1)
      Tn4[j]<-mean(xxx)-yy 
      hatsigma1<-(sd(xxx)/sqrt(30))
    }}
  Qn1[k]<-(sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*hatsigma1)
  Times1[k]<-(abs(Qn1[k])>Z_alpha) 
}

densityf<-function(y,alpha,c){
  exp(-(y^2-2*alpha*(abs(y-c)-abs(c))+alpha^2)/2)/sqrt(2*pi)-
    alpha*exp(2*alpha*abs(y-c))*(1-pnorm(abs(c)+abs(y-c)+alpha,0,1))
}   
YY<-seq(-5,5,length=100)
Alpha<-mu2-mu1
res<-array(0,dim=c(length(YY),length(Alpha)))
c<-0
for(j in 1:length(YY)){
  for(k in 1:length(Alpha))
    res[j,k]<-densityf(YY[j],Alpha[k],0)
}
plot(YY,res,type="l",ylab="density",xlab="",col="blue")
lines(density(Qn1),ylab="density",col="red",xlab="")
mean(Times1)




############################################################################
############################################################################
############################################################################
############################################################################
##################power vary along with different H1########################
############################################################################
############################################################################
############################################################################
############################################################################

library(nleqslv)  
FindRoot = function(x,gamma1,gamma2,alpha) {
  return(pnorm(x-gamma2)-exp(2*gamma2*x)*pnorm(-x-gamma2)-1+alpha)
}

MU<-seq(0.499,0.488,by=-0.001)
Power<-array(0,dim=c(length(MU),2))
##under Ho######
mu11<-0.5
mu22<--0.5
n1<-20
n2<-20
sigma1<-1
sigma2<-1
Z_alpha<-nleqslv(1, FindRoot, control=list(btol=.000001),jacobian=TRUE,method="Newton",gamma1=gamma1,gamma2=mu22-mu11,alpha= 0.05)$x
for(l in 1:length(MU)){
  MAX<-5000
  max<-800
  Qn1<-rep(0,max)
  Times1<-rep(0,max)
  Times2<-rep(0,max)
  #Times3<-rep(0,max)
  ####under H1########
  mu1<- MU[l]
  mu2<--MU[l]
  ####under H1########
  
  hatsigma1<-sigma1
  hatsigma2<-sigma2
  
  for(k in 1:max){
    Tn1<-rep(0,MAX)
    Tn2<-rep(0,MAX)
    Tn3<-rep(0,MAX)
    Tn4<-rep(0,MAX)
    XX<-array(0,dim=c(MAX,n1))
    YY<-array(0,dim=c(MAX,n2))
    for(j in 1:MAX){
      if(((sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*sqrt((hatsigma1^2/n1+hatsigma2^2/n2))))<=0){
        xx<-rnorm(n1,mu1,sigma1)
        yy<-rnorm(n2,mu2,sigma2)
        XX[j,]<-xx
        YY[j,]<-yy
        hatsigma1<-sd(xx)
        hatsigma2<-sd(yy)
        Tn1[j]<-(mean(xx)-mean(yy)-(mu11-mu22))
        Tn2[j]<-mean(xx)-mean(yy)
      }else{
        XX[j,]<-xx
        YY[j,]<-yy
        xx<-rnorm(n1,mu1,sigma1)
        yy<-rnorm(n2,mu2,sigma2)
        hatsigma1<-sd(xx)
        hatsigma2<-sd(yy)
        Tn3[j]<-(mean(yy)-mean(xx)-(mu22-mu11))
        Tn4[j]<-mean(yy)-mean(xx)}      
    }
    Qn1[k]<-sum(Tn2+Tn4)/(MAX)+sum(Tn1+Tn3)/(sqrt(MAX)*sqrt((hatsigma1^2/n1+hatsigma2^2/n2)))
    Times1[k]<-(abs(Qn1[k])>Z_alpha)+0
    #Times2[k]<-(mean(XX)-mean(YY)-(mu11-mu22))/(sqrt(var(c(XX))/(n1*MAX)+var(c(YY))/(n2*MAX)))<qnorm(0.05)
    Times2[k]<-(mean(XX)-mean(YY)-(mu11-mu22))/(sqrt(var(c(XX))/(n1*MAX)+var(c(YY))/(n2*MAX)))<qt(0.05,n1*MAX+n2*MAX-2)
  }
  Power[l,1]<-mean(Times1)
  Power[l,2]<-mean(Times2)
  print(l)
}
Power
indx<-1:l
plot(MU[indx],Power[indx,1],type="l",ylab="power",xlab="",col="blue",ylim=c(0,1))
lines(MU[indx],Power[indx,2],col="red",xlab="") 
lines(MU[indx],Power[indx,3],col="green",xlab="") 





