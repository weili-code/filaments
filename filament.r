# This script implement the Bayesian estimation for filaments on the earthquake data
# Author: Wei Li
# 
# ref: Posterior Contraction and Credible Sets for Filaments of Regression Function
# Li, W. and Ghosal, S., Electronic Journal of Statistics, 14, 1707-1743, 2020.
#    

 
library(mapproj)
library(mapdata)
library(scales)
#library(ggplot2)
set.seed(1355)

#----------------------load data -------------------------------------------#

earth_dat <- read.csv(file="earthquake_dat2.csv", header=TRUE, sep=",")
row.names(earth_dat) <- NULL
earth_dat
nrow(earth_dat)

# subset data for california only
earth_dat<-subset(earth_dat, !(   latitude>39.174 & latitude<42.163   & longitude> -119  & longitude < -113.774 )  )
earth_dat<-subset(earth_dat, !(   latitude>38.4 & latitude<39.174   & longitude> -118.103  & longitude < -113.774 )  )
earth_dat<-subset(earth_dat, !(   latitude>37.276& latitude<38.4   & longitude> -116.78  & longitude < -113.774 )  )
earth_dat<-subset(earth_dat, !(   latitude>36 & latitude<37.276   & longitude> -116.038 & longitude < -113.774 )  )

coord <- mapproject(earth_dat[,2], earth_dat[,1] , proj="gilbert", orientation=c(90, 0, 225))  #convert points to projected lat/long

#map(database= "state", regions="California",col="grey80", fill=TRUE, projection="gilbert", orientation= c(90,0,225))
#points(coord$x,coord$y , pch=21, cex=1, bg=4, col="red", lwd=.4)  #plot converted points
#points(coord$x,coord$y , pch=21, cex=rescale(earth_dat[,3], to = c(.5, 3)) , col="red", lwd=1)  #plot converted points


X<-cbind( coord$x,coord$y  )
y<-earth_dat[,3]      
dat<-cbind(X,y)
nobs<-length(y)
N<-nobs


library("splines")
library("splines2")
#library(rgl)
library(MASS)
library(pracma)
library("gtools")
library(parallel)

#source("compute_hausdorff.R")


###########################################
###### define a function ##################
###########################################

Bdata<-function(nobs,x1bs,x2bs,J1,J2){
  out<-matrix(,0,J1*J2)
  for (i in 1:nobs ){
    out<-rbind(out, kronecker(x1bs[i,], x2bs[i,]))
  } 
  return(out)
}


#########################################
# find the optimal choice of number of J's
#########################################

min(X[,1])
max(X[,1])
min(X[,2])
max(X[,2])




###############################################################
# parameters (optimal choice)
optJ<-c(32,32)
q1=4; k1=optJ[1]-q1   # k1 is number of interior knots, q1=order
J1=q1+k1;

# parameters (optimal choice)
q2=4; k2=optJ[2]-q2   # k2 is number of interior knots, q2=order
J2=q2+k2;

knots1<-seq(min(X[,1])-.05, max(X[,1])+.05,length=k1+2)
bdknots1<-knots1[c(1,length(knots1))]
inknots1<-knots1[2:(length(knots1)-1)]
knots2<-seq(min(X[,2])-.05, max(X[,2])+.05,length=k2+2)
bdknots2<-knots2[c(1,length(knots2))]
inknots2<-knots2[2:(length(knots2)-1)]

x1<-dat[,1]
x2<-dat[,2]
x1bs<-bs(x1,knots=inknots1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)
x2bs<-bs(x2,knots=inknots2, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)

B<-Bdata(nobs,x1bs,x2bs,J1,J2)  # the capital B matrix which is n by J1*J2



##############################################################################
## constructing function, derivatives, Hessian 
##############################################################################

############# regression function ####################
f_theta<-function(x1,x2,theta) {
# this function evaluates the B-spline approximation given by theta, at point (x1,x2)
# input: x1, x2 should be scalars,theta coefficients ordered in dictionary order
	res1<-bs(x1, knots=inknots1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,] 
	res2<-bs(x2, knots=inknots2, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
	#res1<-bbase2(x1, min(x_test),max(x_test),k1+1,q1-1)
	#res2<-bbase2(x2, min(x_test),max(x_test),k2+1,q2-1) 
	out1<- kronecker(res1, res2)%*%as.vector(theta)
	output<-as.numeric(out1)
	return(output)
}


f10<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
	temp1<-dbs(x1, knots=inknots1, derivs=1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,] 
    temp2<-bs(x2, knots=inknots2, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))
 
}	
	

f20<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
 	temp1<-dbs(x1, knots=inknots1, derivs=2, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,]  
    temp2<-bs(x2, knots=inknots2, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))
 
}


f30<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
 	temp1<-dbs(x1, knots=inknots1, derivs=3, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,]  
    temp2<-bs(x2, knots=inknots2, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))
 
}


f21<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
 	temp1<-dbs(x1, knots=inknots1, derivs=2, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,]  
    temp2<-dbs(x2, knots=inknots2, derivs=1, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))
 
}




f01<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane,
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
    
    temp1<-bs(x1, knots=inknots1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,] 
    temp2<-dbs(x2, knots=inknots2, derivs=1, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))


}


f02<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  

    temp1<-bs(x1, knots=inknots1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,] 
    temp2<-dbs(x2, knots=inknots2, derivs=2, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))
 
}


f12<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  

    temp1<-dbs(x1, knots=inknots1, derivs=1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,] 
    temp2<-dbs(x2, knots=inknots2, derivs=2, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))
 
}

f03<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  

    temp1<-bs(x1, knots=inknots1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,] 
    temp2<-dbs(x2, knots=inknots2, derivs=3, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))
 
}



f11<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  

    temp1<-dbs(x1, knots=inknots1,derivs=1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,] 
    temp2<-dbs(x2, knots=inknots2,derivs=1, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))

}



# filament points equation for regression function (with B spline approximation)
filamentEqn<-function(x,index=NULL,theta) {
# input: x is a vector, index,theta
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
  
  if(is.null(index)) {
    x1<-x[1]
    x2<-x[2]
    temp1<-2*f10(x1,x2,theta=theta)*(f20(x1,x2,theta=theta)-f02(x1,x2,theta=theta)+f11(x1,x2,theta=theta)-sqrt((f02(x1,x2,theta=theta)-f20(x1,x2,theta=theta))^2+4*(f11(x1,x2,theta=theta))^2 ) )
    temp2<-f01(x1,x2,theta=theta)*(f02(x1,x2,theta=theta)-f20(x1,x2,theta=theta)+4*f11(x1,x2,theta=theta)-sqrt((f02(x1,x2,theta=theta)-f20(x1,x2,theta=theta))^2+4*(f11(x1,x2,theta=theta))^2 ) ) 
    out<-temp1+temp2
    return(out)
    } else {
    x1<-x[1]
    x2<-x[2]
    temp1<-2*f10(x1,x2,index,theta=theta)*(f20(x1,x2,index,theta=theta)-f02(x1,x2,index,theta=theta)+f11(x1,x2,index,theta=theta)-sqrt((f02(x1,x2,index,theta=theta)-f20(x1,x2,index,theta=theta))^2+4*(f11(x1,x2,index,theta=theta))^2 ) )
    temp2<-f01(x1,x2,index,theta=theta)*(f02(x1,x2,index,theta=theta)-f20(x1,x2,index,theta=theta)+4*f11(x1,x2,index,theta=theta)-sqrt((f02(x1,x2,index,theta=theta)-f20(x1,x2,index,theta=theta))^2+4*(f11(x1,x2,index,theta=theta))^2 ) ) 
    out<-temp1+temp2
    return(out)
  }  
}


########  gradient ##############
gradBs<-function(x,index=NULL,theta) {
# to-pass values: 
# input: x is a vector, index,theta
# index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
# output: 2-dim vector
  if(is.null(index)) {
  return(c(f10(x[1],x[2],theta=theta), f01(x[1],x[2],theta=theta)))
  } else {
  return(c(f10(x[1],x[2],index,theta=theta), f01(x[1],x[2],index,theta=theta)))
  }
}

#######  second eigenvector ##############
Veign<- function(x,index=NULL,theta){
# to-pass values: 
# input: x is a vector, index,theta
# index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
# output: 2-dim vector
  if(is.null(index)) {
  x1<-x[1]
  x2<-x[2]
  a1<- 2*(f20(x1,x2,theta=theta)-f02(x1,x2,theta=theta)+f11(x1,x2,theta=theta)-sqrt((f02(x1,x2,theta=theta)-f20(x1,x2,theta=theta))^2+4*(f11(x1,x2,theta=theta))^2 ) )
  a2<- (f02(x1,x2,theta=theta)-f20(x1,x2,theta=theta)+4*f11(x1,x2,theta=theta)-sqrt((f02(x1,x2,theta=theta)-f20(x1,x2,theta=theta))^2+4*(f11(x1,x2,theta=theta))^2 ) ) 
  return( c(a1,a2))
  } else {
  x1<-x[1]
  x2<-x[2]
  a1<- 2*(f20(x1,x2,index,theta=theta)-f02(x1,x2,index,theta=theta)+f11(x1,x2,index,theta=theta)-sqrt((f02(x1,x2,index,theta=theta)-f20(x1,x2,index,theta=theta))^2+4*(f11(x1,x2,index,theta=theta))^2 ) )
  a2<- (f02(x1,x2,index,theta=theta)-f20(x1,x2,index,theta=theta)+4*f11(x1,x2,index,theta=theta)-sqrt((f02(x1,x2,index,theta=theta)-f20(x1,x2,index,theta=theta))^2+4*(f11(x1,x2,index,theta=theta))^2 ) ) 
  return( c(a1,a2))
  }
}


#######  second eigenvalue ##############
secondeign<- function(x,index=NULL,theta){
# to-pass values: 
# input: x is a vector, index,theta
# index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
# output: 2-dim vector
  
  if(is.null(index)) {
  return(  (f20(x[1],x[2],theta=theta)+f02(x[1],x[2],theta=theta)- sqrt( (f20(x[1],x[2],theta=theta)-f02(x[1],x[2],theta=theta))^2+ 4*(f11(x[1],x[2],theta=theta))^2)  )/2)
  } else {
  return(  (f20(x[1],x[2],index,theta=theta)+f02(x[1],x[2],index,theta=theta)- sqrt( (f20(x[1],x[2],index,theta=theta)-f02(x[1],x[2],index,theta=theta))^2+ 4*(f11(x[1],x[2],index,theta=theta))^2)  )/2)
  }
}

#######  first eigenvalue ##############
firsteign<- function(x,index=NULL,theta){
# to-pass values: 
# input: x is a vector, index
# index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
# output: 2-dim vector
  if(is.null(index)) {
    return(  (f20(x[1],x[2],theta=theta)+f02(x[1],x[2],theta=theta)+ sqrt( (f20(x[1],x[2],theta=theta)-f02(x[1],x[2],theta=theta))^2+ 4*(f11(x[1],x[2],theta=theta))^2)  )/2)
  } else {
    return(  (f20(x[1],x[2],index,theta=theta)+f02(x[1],x[2],index,theta=theta)+ sqrt( (f20(x[1],x[2],index,theta=theta)-f02(x[1],x[2],index,theta=theta))^2+ 4*(f11(x[1],x[2],index,theta=theta))^2)  )/2)
  }
}





################### priors ###########################

theta_0<-matrix(rep(0,J1*J2))
Lambda_0<-diag(1,J1*J2)
# empirical bayes estiamte for sigma2
sigma2_hat<-nobs^{-1}*t(y-B%*%theta_0)%*%solve(B%*%Lambda_0%*%t(B)+diag(1,nobs))%*%(y-B%*%theta_0)

################## posterior sampling ################

# posterior samples for theta --------------------------------------------
a<-6
b<-2
nsim<-100 # number of posterior samples
#sigma2.post<-1/rgamma(nsim, (a+nobs)/2, rate=(b+nobs*sigma2_hat)/2 ) # 1*nsim vector of posterior sample for sigma2
sigma2_hat # as alternative for sigma2.post

BB_lambda<- (solve(Lambda_0)+t(B)%*%B)
BB_lambda_inv<-solve(BB_lambda)
theta_post_EXPmean<- BB_lambda_inv%*%(t(B)%*%y+solve(Lambda_0)%*%theta_0)

#E<-matrix(rnorm(nsim*J1*J2,0,sqrt(sigma2.post)), nsim, J1*J2,byrow=FALSE ) 
# E is nsim by J1J2 matrix. each row corresponds to nsim iid simulated normals with sd corresponding to one element of sigma2.post

#using empirical sigma2_hat. E is nsim by J1J2 matrix. each row corresponds to nsim iid simulated normals with sd=sigma2_hat
E<-matrix(rnorm(nsim*J1*J2,0,sqrt(sigma2_hat)), nsim, J1*J2,byrow=FALSE ) 

chol_BB_lambda_inv<-chol(BB_lambda_inv)
theta.post<-t(t(chol_BB_lambda_inv)%*%t(E)+ matrix(theta_post_EXPmean,ncol=nsim, nrow=J1*J2, byrow=FALSE))  # nsim*J1J2 matrix of posterior sample for theta
theta_postmean<-as.matrix(colMeans(theta.post)) #monte carlo mean of the regression function
#thetaMat_postmean<-matrix(theta_postmean,nrow=J1,byrow = TRUE)  # matrix of theta, rows=first index, columns=second index
#as.vector(t(thetaMat_postmean))-theta_postmean=0

#########################################################################
#-------------compute the 1-gamma quantile for GP(0,Sigma)-----------------------------
#########################################################################


######################################################
#### plot the posterior mean regression function #####
######################################################
x1s<-seq(bdknots1[1],bdknots1[2],length.out=100)
x2s<-seq(bdknots2[1],bdknots2[2],length.out=100)

f_post<-function(x1,x2){
	f_theta(x1,x2,as.vector(theta_postmean))
}
fz_post<-outer(x1s,x2s,Vectorize(f_post))
#write.table(fdenz_post, file="fdenz_post.txt", row.names=FALSE, col.names=FALSE)
#fdenz_post<- matrix(scan("fdenz_post.txt"), nrow=1000, byrow=TRUE)

#contour(x1s,x2s,fz_post,xlab="x1",ylab="x2")
#persp3d(x1s,x2s,fz_post,col="skyblue",aspect = c(1, 1, 1))

######################################################
######### search for filaments #######################
######################################################


niter<-500; step_len<-.000005; xdat<-X; tol<-3; stop<-1e-6

########### a wrapper for sc_mean_shift for all observations ########
sc_mean_shift<-function(k){
          theta<-theta_postmean
          out<-matrix(,0,2)
          for (i in 1:(N)){
             if(abs(f_theta(xdat[i,1],xdat[i,2],theta_postmean ))>tol){
                x_c<-as.matrix(xdat[i,])
                for (i in 1:niter){
                  V<-Veign(as.vector(x_c),theta=theta)
                  V<-V/sqrt(sum(V^2))  # normalized V
                  if ( abs(t(as.matrix(V))%*%as.matrix(gradBs(as.vector(x_c),theta=theta))) <stop ) {
                          break }
                  x_c<-x_c+step_len*as.matrix(V)%*%t(as.matrix(V))%*%as.matrix(gradBs(as.vector(x_c),theta=theta))
                }
                  out<-rbind(out,t(x_c))   
            }
          }
    return(out)
}

detectCores()
filament_all<-mclapply(1:1,mc.cores=detectCores(), sc_mean_shift)
class(filament_all)
length(filament_all)


################## plot #################################


pdf("filament.pdf",width=6.8,height=6.8)
	map(database= "state", regions="California",col="grey80", fill=TRUE, projection="gilbert", orientation= c(90,0,225))
    points(coord$x,coord$y , pch=21, cex=rescale(earth_dat[,3], to = c(.5, 3)) , col="red", lwd=1)  #plot converted points
    #contour(x1s,x2s,fz_post, xlab=expression(X[1]), ylab=expression(X[2]))    
    for (k in 1:1){
			for (i in 1: nrow(filament_all[[k]])){
                  if (k==1){ # the filament using posterior mean of parameters
                    points(filament_all[[k]][i,1],filament_all[[k]][i,2] ,xlab="x1",ylab="x2",,col="blue",pch=20,cex=1)
                  } else {
                    #points(filament_all[[k]][i,1],filament_all[[k]][i,2] ,xlab="x1",ylab="x2",col=k+2,lwd=.1,lty=2,pch=20,cex=.2)
				    #points(x_c[1,1],x_c[2,1] ,xlab="x1",ylab="x2",col=k+2,lwd=.1,lty=2,pch=20,cex=.1)
                  }
                 
         }
	}
dev.off() 


##########################################################
save.image(file="filament.RData") # nsim=100, J=32, rho=1.2, gamma=.1
