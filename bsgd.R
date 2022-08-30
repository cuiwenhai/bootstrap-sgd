 


  
#����ÿ������
R<-function(maxiter,B){
x<-array(runif(500,-10,10),c(100,5))#�Ա���
theta<-c(0.2,0.4,0.6,0.8,1) #Ҫ���ƵĲ���
y<-x%*%theta+rnorm(100)#ģ��

#��������ֵ
btsgd<-function(x,y,maxiter,B)#B��bootsrap�ܴ���,b�ǵڼ���bootstrap
  { step=0.5/((maxiter)^0.5)
    m<-nrow(x)
    n<-ncol(x) 
    theta<-matrix(rep(0,n*B),n,B)#ktheta��ʼֵ������Ϊ0
    new_theta<-theta
    iter<-0   #��������
    k<-0  #��k������
    #����bootstrap����
    W<-rnorm(maxiter,1,1)
    Wb<-matrix(data = 0,maxiter,B)
    for (b in 1:B) {
      Wb[,b]<-sample(W,maxiter,replace = TRUE)
    }
    #���²���
    while(iter<maxiter){
      iter<-iter+1
      k<-sample(1:m, 1, replace = TRUE)
      xk<-x[k,,drop=FALSE]
      yk<-y[k,,drop=FALSE]
      for (b in 1:B) {
        hk<-xk%*%theta[,b]
        des<-t((hk-yk)%*%xk)
        new_theta[,b]<-theta[,b]-Wb[iter,b]*step*des
        #newerror<-t(new_theta-theta)%*%(new_theta-theta)
        theta[,b]<-(new_theta[,b] /iter)+theta[,b]*(1-1/iter)
      }
       }
    #costfunction<-(t(x%*%theta-y)%*%(x%*%theta-y))/m
    result<-theta[n,]
    sresult<-sort(result)
    r<-sresult[-c(1,2,3,4,9,10,11,12)]#ɾ���쳣��
    #names(result)<-c('ϵ��','ϵ��','ϵ��','ϵ��','ϵ��','��������','���')
    return(r)
}
return(btsgd(x,y,maxiter,B))
}
R(1000,12)
sd(R(1000,12))


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
MAX<-30
max<-300
Qn1<-rep(0,max)
Times1<-rep(0,max)
#n1<-50
#n2<-50
#sigma1<-1
#sigma2<-1
mu<-1
d0<-1.01
rm<-matrix(data=0,max,MAX)

Z_alpha<-nleqslv(1, FindRoot, control=list(btol=.000001),jacobian=TRUE,method="Newton",gamma1=gamma1,gamma2=d0,alpha= 0.05)$x
#hatsigma1<-sigma1
#hatsigma2<-sigma2
for(k in 1:max){
Tn1<-rep(0,MAX)
Tn2<-rep(0,MAX)
Tn3<-rep(0,MAX)
Tn4<-rep(0,MAX)
print(k)
for(j in 1:MAX){
  B=12
  xx=R(1000,B)#��������
  mean<-mean(xx)
  rm[k,j]=mean
  if(j<2){hatsigma1=0.048
  h<-rnorm(1,0,1) }
  else{hatsigma1<-sd(rm[k,1:j])
  h<--1}
  #hatsigma1=0.006
 if(((sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*hatsigma1))<=0 & h<0){
  Tn1[j]<-mean-d0
  Tn2[j]<-mean
  }else{
Tn3[j]<-d0-mean
    Tn4[j]<--mean
  }}
Qn1[k]<-(sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*hatsigma1)
Times1[k]<-(abs(Qn1[k])>Z_alpha) 

    }
densityf<-function(y,alpha,c){
exp(-(y^2-2*alpha*(abs(y-c)-abs(c))+alpha^2)/2)/sqrt(2*pi)-
alpha*exp(2*alpha*abs(y-c))*(1-pnorm(abs(c)+abs(y-c)+alpha,0,1))
}   
YY<-seq(-5,5,length=100)
Alpha<--d0
res<-array(0,dim=c(length(YY),length(Alpha)))
c<-0
for(j in 1:length(YY)){
for(k in 1:length(Alpha))
res[j,k]<-densityf(YY[j],Alpha[k],0)
}
plot(YY,res,type="l",ylab="density",xlab="",col="blue",panel.first=grid(12,12,col="gray70"),xlim = c(-5,5))
lines(density(Qn1),ylab="density",col="red",xlab="")
mean(Times1)














#######                                                    #######
##################    ��α�����bsgd�µ���ʽ     ####################
########################            ################################
############################     ####################################

library(nleqslv)  
FindRoot = function(x,gamma1,gamma2,alpha) {
  return(pnorm(x-gamma2)-exp(2*gamma2*x)*pnorm(-x-gamma2)-1+alpha)
}
MAX<-100
max<-200
Qn1<-rep(0,max)
Times1<-rep(0,max)
#n1<-50
#n2<-50
#sigma1<-1
#sigma2<-1
mu<-0.2
d0<-0.2
B<-10

rm<-matrix(data=0,max,MAX)
Z_alpha<-nleqslv(1, FindRoot, control=list(btol=.000001),jacobian=TRUE,method="Newton",gamma1=gamma1,gamma2=d0,alpha= 0.05)$x
####����bootstrap����
W<-rnorm(maxiter,1,1)
Wb<-matrix(data = 0,maxiter,B)
for (b in 1:B) {
  Wb[,b]<-sample(W,maxiter,replace = TRUE)
}
####����ͳ����
for(k in 1:max){
  Tn1<-rep(0,MAX)
  Tn2<-rep(0,MAX)
  Tn3<-rep(0,MAX)
  Tn4<-rep(0,MAX)
  print(k)
  for(j in 1:MAX){
  ###����������
      x<-array(runif(100,10,100),c(100,1))#�Ա���
      theta<-c(0.2) #Ҫ���ƵĲ���
      y<-x%*%theta+rnorm(100)#ģ��
  ###���ɹ���ֵ 
      step=0.001
      m<-nrow(x)
      n<-ncol(x) 
      theta<-matrix(rep(0,n*B),n,B)#ktheta��ʼֵ������Ϊ0
      new_theta<-theta
      iter<-0   #��������
      k<-0  #��k������
      #���²���
      while(iter<maxiter){
          iter<-iter+1
          k<-sample(1:m, 1, replace = TRUE)
          xk<-x[k,,drop=FALSE]
          yk<-y[k,,drop=FALSE]
          for (b in 1:B) {
            hk<-xk%*%theta[,b]
            des<-t((hk-yk)%*%xk)
            new_theta[,b]<-theta[,b]-Wb[iter,b]*step*des
            theta[,b]<-(new_theta[,b] /iter)+theta[,b]*(1-1/iter)
          }
        }
        xx<-c(theta)
#######��������ֵ   
    
    mean<-mean(xx)
    rm[k,j]=mean
    if(j<2){hatsigma1=0.006}
    else{hatsigma1<-sd(rm[k,1:j])}
    #hatsigma1=0.006
    if(((sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*hatsigma1))<=0){
      Tn1[j]<-mean-d0
      Tn2[j]<-mean
    }else{
      Tn3[j]<-d0-mean
      Tn4[j]<--mean
    }}
  Qn1[k]<-(sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*hatsigma1)
  Times1[k]<-(abs(Qn1[k])>Z_alpha) 
}
densityf<-function(y,alpha,c){
  exp(-(y^2-2*alpha*(abs(y-c)-abs(c))+alpha^2)/2)/sqrt(2*pi)-
    alpha*exp(2*alpha*abs(y-c))*(1-pnorm(abs(c)+abs(y-c)+alpha,0,1))
}   
YY<-seq(-5,5,length=100)
Alpha<--d0
res<-array(0,dim=c(length(YY),length(Alpha)))
c<-0
for(j in 1:length(YY)){
  for(k in 1:length(Alpha))
    res[j,k]<-densityf(YY[j],Alpha[k],0)
}
plot(YY,res,type="l",ylab="density",xlab="",col="blue",panel.first=grid(10,10,col="gray70"))
lines(density(Qn1),ylab="density",col="red",xlab="")
mean(Times1)

