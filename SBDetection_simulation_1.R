
rm(list=ls(all=TRUE))

library("glmnet")
library("matrixcalc")
library("fields")
library("vars")
library("MTS")
library("mvtnorm")
library("xtable")
library("lattice")
library(ggplot2)

######## FUNCTIONS ##########################

source("functions_SBDetection.R")

#############################################

#############################################

#############################################
######## DATA GENERATION ####################
#############################################

# for real data application, this part will be replaced by loading the data!
T <- 300; # number of time points
k <- 20; # number of time series

# TRUE BREAK POINTS WITH T+1 AS THE LAST ELEMENT
brk <- c(floor(T/3),floor(2*T/3),T+1)

m <- length(brk)
p.t <- 1; ## the true AR order
########## PHI GENERATION ###################
phi.full <- matrix(0,k,k*p.t*m)
aa <- 0.75; bb <- 0.75
for (j in 1:(k-1)){
  phi.full[j,((1-1)*p.t*k+j+1)] <- -0.6;
  phi.full[j,((2-1)*p.t*k+j+1)] <- 0.75;
  phi.full[j,((3-1)*p.t*k+j+1)] <- -0.8;
}

print(plot.matrix(abs(phi.full),p=m,name="TRUTH"))

for(i in 2:m){
  phi.1 <- phi.full[, ((i-2)*k*p.t+1):((i-1)*k*p.t)];
  phi.2 <- phi.full[, ((i-1)*k*p.t+1):((i-0)*k*p.t)];
  print(sqrt( sum((phi.1-phi.2)^2)   ))
}


#############################################
######## DATA ANALYSIS   ####################
#############################################


###### METHODS ##############################
method.all <- c("VAR","LASSO","HVARC","HVAROO","HVARELEM","SSLASSO","SSHVARC","SSHVAROO","SSHVARELEM",
                "DHVAR","DHVARC","SSDHVAR","SSDHVARC");
method.all <- c("LASSO");
final.check <- matrix(0,length(method.all)+1,6)
final.check[2:(length(method.all)+1),1] <- method.all; phi.hat.all <- vector("list",length(method.all));
final.check[1,2] <- c("MSPE"); final.check[1,3] <- c("SD of MSPE"); final.check[1,4] <- c("MEDIAN of MSPE"); 
final.check[1,5] <- c("MEAN of RPE"); final.check[1,6] <- c("SD of RPE");
count.non <- 0; count.zero <- 0; p <- p.t;  
for (i in 1:k){
  for (j in 1:(k*p)){
    if ( phi[i,j] != 0  ){count.non <- count.non + 1;}
    if ( phi[i,j] == 0  ){count.zero <- count.zero + 1;}
  }
}
est.error <-  matrix(0,length(method.all)+2,9)
est.error[3:(length(method.all)+2),1] <- method.all; est.error[1,2]<-c("l2 est error"); est.error[1,3]<-c("SD of l2"); 
est.error[1,4] <- c("TZ"); est.error[2,4] <- c(count.zero); est.error[2,1] <- c("TRUTH");
est.error[1,5] <- c("TNZ"); est.error[2,5] <- c(count.non);
est.error[1,6] <- c("FZ"); est.error[1,7] <- c("FNZ"); 
est.error[1,8] <- c("lag.L.mean"); est.error[1,9] <- c("lag.L.sd"); 
#############################################
#############################################
N <- 1; # number of replicates
iter.final <- rep(0,length(method.all));
pts <- vector("list",N); data.full <- vector("list",N); pts.final <- vector("list",N);
phi.final <- vector("list",N); pts.final.sbs <- vector("list",N);


########### GENERAL PARAMETERS ##############
tol <- 4*10^(-2); # tolerance 
step.size <- 2*10^(-4); # step size 
max.iteration <- 100; # max number of iteration for the LASSo solution
p <- p.t; ## the selected AR order
sig <- matrix(0,k,k);
rho <- 0.5;
for(i in 1:k){for(j in 1:k){sig[i,j] <- rho^(abs(i-j))}}


for ( j.1 in 1:N){
  
  #############################################
  set.seed(123456*j.1)

  e.sigma <- as.matrix(0.01*diag(k));
  try=var.sim.break(T,arlags=seq(1,p.t,1),malags=NULL,phi=phi.full,sigma=e.sigma,brk = brk)
  data <- try$series
  data <- as.matrix(data)
  data.full[[j.1]] <- data;
  #############################################
######################################################################
######## FIRST STEP : INITIAL BRK POINT SELECTION ####################
######################################################################
method<- c("LASSO"); # other types of penalization has not been coded yet, so the only option as of now is just LASSO!
lambda.1 <- 0.50; #The first tuning parameter 
lambda.1.cv <- seq(0.10,0.50,0.05); #The first tuning parameter
cv.length <- 10;
cv.index <- seq(p+cv.length,T-1,cv.length); # cv index
temp.first <- first.step.cv.new(method,data, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4),cv.index)
mspe.plot(method,temp.first$cv,lambda.1.cv,temp.first$cv.final, "1")
fisrt.brk.points <- temp.first$brk.points;
# plotting the data with the candidate break points from initial step.
MTSplot(data)
abline(v=fisrt.brk.points)
print(fisrt.brk.points);
pts[[j.1]] <- fisrt.brk.points;

######################################################################
######## SECOND STEP : SCREENING                  ####################
######################################################################
n <- T - p;
omega <- (1/1.75)*(((log(n))^1)*log(k)); # the penalty term in the information criterion -- the higher omega, the smaller number of break points selected.
lambda.2 <- (1/1)*(log(n)*log(k))/n # the second tuning parameter. This default number seems to be working for many simulation and real data examples!
temp <- second.step(data, lambda = lambda.2, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), fisrt.brk.points, omega)
# final break points selected
final.brk.points <- temp$pts;
# plotting the data with selected break points
MTSplot(data)
abline(v=final.brk.points)
pts.final[[j.1]] <- final.brk.points;

######################################################################
######## THIRD STEP : AR PARAMETER ESTIMATION     ####################
######################################################################
lambda.final <- seq(0.0025,0.25, 0.0025); # the range of tuning parameter for AR parameter estimation 
r.n <- floor(omega)+1; # the length of time interval removal

par.fit <-  ar.est(method, data, weight = NULL, lambda.final, p, final.brk.points, r.n, max.iteration = 1000, tol = tol, step.size = 10^(-3))

# AR parameter estimates
phi.est <- par.fit$phi.hat;
# checkinh the plot of MSPE for the LASSO
mspe.plot(method,par.fit$pred.error,lambda.final,par.fit$tune.final, "1")
# plotting the AR parameter estimates
print(plot.matrix(abs(par.fit$phi.hat),length(final.brk.points)+1,name = ""))

phi.final[[j.1]] <- phi.est;

}



############################################
############# DETECTION CHECK     ##########
############################################
############################################
m <- length(brk);
pts.final.full <- vector("list",N);
for(i in 1:N){
  if ( length(pts.final[[i]]) > (m-1)   ){print("OVER-ESTIAMTING"); break;}
  if ( length(pts.final[[i]]) == (m-1)   ){pts.final.full[[i]] <- pts.final[[i]];}
  if ( length(pts.final[[i]]) == 0 ) {  pts.final.full[[i]] <- rep(0,m-1);  }
  if ( length(pts.final[[i]]) > 0 && length(pts.final[[i]]) < (m-1) ){
    ll <- length(pts.final[[i]]); pts.final.full[[i]] <- rep(0,m-1);
    for(j in 1:ll){
      if (  pts.final[[i]][j] < (brk[1] + (1/2)*(brk[2] - brk[1]) )   ) {pts.final.full[[i]][1] <- pts.final[[i]][j];}
      for(kk in 2:(m-1)){
        if (  pts.final[[i]][j] >=  (brk[(kk-1)] + (1/2)*(brk[kk] - brk[(kk-1)]) ) && pts.final[[i]][j] <  (brk[(kk)] + (1/2)*(brk[kk+1] - brk[(kk)]) )   ){
          pts.final.full[[i]][kk] <- pts.final[[i]][j];
        }
      }
    }
  }
}

detection <- matrix(0,m,5)

detection[1,1] <- c("break points"); detection[1,2] <- c("truth"); detection[1,3] <- c("mean");
detection[1,4] <- c("std"); detection[1,5] <- c("selection rate");
for(i in 1:(m-1)){
  detection[(i+1),1] <- c(i); detection[(i+1),2] <- c(brk[i]/T); loc <- rep(0,N);
  for(j in 1:N){
    temp <- pts.final.full[[j]]; l <- length(temp); loc[j] <- temp[i];
  }
  loc <- loc[which(loc!=0)]; T.new <- length(loc); detection[(i+1),3] <- mean(loc/T); detection[(i+1),4] <- sd(loc/T); detection[(i+1),5] <- T.new/N;
}

for(i in 2:(m)){
  for(j in 2:5){
    detection[i,j] <- round(as.numeric(detection[i,j]),digits = 4)
  }
}

table.1 <- xtable(detection, hline.after = c(1,2,7))
align(table.1) <- "r|r"
print(table.1,include.rownames = FALSE, include.colnames = FALSE)


############################################
############# PHI SELECTION CHECK ##########
############################################

count.non <- 0; count.zero <- 0; p <- p.t; m <- length(brk) 
for (i in 1:k){
  for (j in 1:(k*m*p)){
    if ( phi.full[i,j] != 0  ){count.non <- count.non + 1;}
    if ( phi.full[i,j] == 0  ){count.zero <- count.zero + 1;}
  }
}
est.error <-  matrix(0,length(method.all)+2,7)
est.error[3:(length(method.all)+2),1] <- method.all; est.error[1,2]<-c("REE"); est.error[1,3]<-c("SD(REE)"); 
est.error[1,4] <- c("TZ"); est.error[2,4] <- c(count.zero); est.error[2,1] <- c("TRUTH");
est.error[1,5] <- c("TNZ"); est.error[2,5] <- c(count.non);
est.error[1,6] <- c("FZ"); est.error[1,7] <- c("FNZ"); 


error <- estimation.check.new(phi.full,phi.final,k,p,m,pts.final.full)
est.error[1+2,2] <- error$l2.error.mean; est.error[1+2,3] <- error$l2.error.sd; 
est.error[1+2,4] <- error$true.zero.median; 
est.error[1+2,5] <- error$true.non.zero.median; 
est.error[1+2,6] <- error$false.zero.median; est.error[1+2,7] <- error$false.non.zero.median; 

est.error[2:3,2:7] <- round(as.numeric(est.error[2:3,2:7]),digits = 4)


##################################################################
##################################################################
##################################################################

