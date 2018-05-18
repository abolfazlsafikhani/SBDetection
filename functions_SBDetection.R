

######## FUNCTIONS #################

var.sim <- function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL,theta = NULL, skip = 200, sigma) {
  if (!is.matrix(sigma)) 
    sigma = as.matrix(sigma)
  k = nrow(sigma)
  nT = nobs + skip
  at = rmvnorm(nT, rep(0, k), sigma)
  nar = length(arlags)
  p = 0
  if (nar > 0) {
    arlags = sort(arlags)
    p = arlags[nar]
  }
  q = 0
  nma = length(malags)
  if (nma > 0) {
    malags = sort(malags)
    q = malags[nma]
  }
  ist = max(p, q) + 1
  zt = matrix(0, nT, k)
  if (length(cnst) == 0) 
    cnst = rep(0, k)
  for (it in ist:nT) {
    tmp = matrix(at[it, ], 1, k)
    if (nma > 0) {
      for (j in 1:nma) {
        jdx = (j - 1) * k
        thej = theta[, (jdx + 1):(jdx + k)]
        atm = matrix(at[it - malags[j], ], 1, k)
        tmp = tmp - atm %*% t(thej)
      }
    }
    if (nar > 0) {
      for (i in 1:nar) {
        idx = (i - 1) * k
        phj = phi[, (idx + 1):(idx + k)]
        ztm = matrix(zt[it - arlags[i], ], 1, k)
        tmp = tmp + ztm %*% t(phj)
      }
    }
    zt[it, ] = cnst + tmp
  }
  zt = zt[(1 + skip):nT, ]
  at = at[(1 + skip):nT, ]
  VARMAsim <- list(series = zt, noises = at)
}

var.sim.break <- function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL,theta = NULL, skip = 200, sigma, brk = nobs+1) {
  if (!is.matrix(sigma)) 
    sigma = as.matrix(sigma)
  k = nrow(sigma)
  m <- length(brk)
  nT = nobs + skip
  at = rmvnorm(nT, rep(0, k), sigma)
  nar = length(arlags)
  p = 0
  if (nar > 0) {
    arlags = sort(arlags)
    p = arlags[nar]
  }
  q = 0
  nma = length(malags)
  if (nma > 0) {
    malags = sort(malags)
    q = malags[nma]
  }
  ist = max(p, q) + 1
  zt = matrix(0, nT, k)
  if (length(cnst) == 0) 
    cnst = rep(0, k)
  if (m == 1){
    for (it in ist:nT) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
  }
  
  if (m > 1){
    
    for (it in ist:(skip+brk[1]-1)) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
    for ( mm in 1:(m-1)){
      for (it in (skip+brk[mm]):(skip+brk[mm+1]-1) ) {
        tmp = matrix(at[it, ], 1, k)
        if (nma > 0) {
          for (j in 1:nma) {
            jdx = (j - 1) * k
            thej = theta[, (jdx + 1):(jdx + k)]
            atm = matrix(at[it - malags[j], ], 1, k)
            tmp = tmp - atm %*% t(thej)
          }
        }
        if (nar > 0) {
          for (i in 1:nar) {
            idx = (i - 1) * k
            phj = phi[, ((mm)*p*k+idx + 1):((mm)*p*k+idx + k)]
            ztm = matrix(zt[it - arlags[i], ], 1, k)
            tmp = tmp + ztm %*% t(phj)
          }
        }
        zt[it, ] = cnst + tmp
      }
    }
  }
  
  zt = zt[(1 + skip):nT, ]
  at = at[(1 + skip):nT, ]
  VARMAsim <- list(series = zt, noises = at)
}

pred <- function(Y,phi,p,T,k,h){
  concat.Y <- matrix(0,k,T+h); concat.Y[,1:T] <- Y[,1:T];
  for ( j in 1:h){
    temp <- matrix(0,k,1);
    for (i in 1:p){temp <- temp +  phi[,((i-1)*k+1):(i*k)]%*%concat.Y[,T+j-i];}
    concat.Y[,T+j] <- temp; 
  }
  return(as.matrix(concat.Y[,T+h]))
}

soft <- function(L,weight,lambda){
  for (i in 1:length(L[1,])){
    lambda <- lambda*(1+weight[i])
    if ( L[i] > lambda){L[i] <- L[i] - lambda}
    if ( L[i] < -lambda){L[i] <- L[i] + lambda}
    if ( abs(L[i]) <= lambda){L[i] <- 0}
  }
  return(L)
}

soft.full <- function(L,lambda,k,p,n){
  
  # for(kk in 1:n){
  #   temp <- L[,((kk-1)*k*p+1):(kk*k*p)];
  #   nrm <- sum(abs(temp))
  #   if ( nrm <= lambda){ L[,((kk-1)*k*p+1):(kk*k*p)] <- matrix(0,k,k*p)  }
  #   if ( nrm > lambda) { L[,((kk-1)*k*p+1):(kk*k*p)] <- L[,((kk-1)*k*p+1):(kk*k*p)] - matrix((lambda/(p*k^2)),k,k*p); }
  # }
  # nrm <- sum(abs(L))
  # if ( nrm <= lambda){ L <- matrix(0,k,k*p)  }
  # if ( nrm > lambda) { L <- L - matrix((lambda/(p*k^2)),k,k*p);  }
  
  
  for (i in 1:length(L[1,])){
    for(j in 1:length(L[,1])){
      if ( L[j,i] > lambda){L[j,i] <- L[j,i] - lambda}
      if ( L[j,i] < -lambda){L[j,i] <- L[j,i] + lambda}
      if ( abs(L[j,i]) <= lambda){L[j,i] <- 0}
    }
  }
  return(L)
}

soft.group <- function(L,weight,group.lag,lambda){
  L <- L/(1+weight);
  for (i in 1:length(group.lag[,1])){
    temp <- group.lag[i,]; temp <- temp[which(temp!=0)];
    L[temp] <- (max(0,1-lambda/(sqrt(sum((L[temp])^2)))))*L[temp]
  }
  return(L*(1+weight))
}

myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  # ColorRamp <- rgb( seq(0,1,length=256),  # Red
  #                   seq(0,1,length=256),  # Green
  #                   seq(1,0,length=256))  # Blue
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- grey(seq(1, 0, length = 256))
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}

shvar.fit <- function(method,data, weight = NULL, lambda, p,T.1,T.2,max.iteration = 1000, tol = 10^(-4)){
  method.full <- c("VAR","LASSO","HVARC","HVAROO","HVARELEM","SLASSO","SHVARC","SHVAROO","SHVARELEM","SSHVARC",
                   "SSHVAROO","SSHVARELEM","SSLASSO","DHVAR","DHVARC","SSDHVAR","SSDHVARC");
  if ( !(method %in% method.full) ){print("ERROR"); break; }
  
  k <- length(data[1,]); T <- length(data[,1]);
  iter <- matrix(0,k,length(lambda));
  phi.hat <- matrix(0,k,k*p); phi.hat.fista <- matrix(0,max.iteration,k*p);
  pred.error <- rep(0,length(lambda)); phi.hat.temp <- matrix(0,k,k*p*length(lambda));
  Y <- as.matrix(t(data)); 
  Y <- as.matrix(Y[,(p+1):T.1]); 
  # Y <- Y%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  Z <- matrix(0,k*p,T.1-p); 
  # Z <- Z%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  for ( i in 1:(T.1-p)){
    for ( j in 1:p){
      Z[((j-1)*k+1):(j*k),i] <- t(data[i+p-j,])
    }
  }
  step.size <- 1/( max(svd(Z)$d)  )^2; 
  # Z.tilde <- vector("list",k);
  # for ( i in 1:k){
  #   Z.new <- matrix(0,k*p,T.1-p);
  #   for (j in 1:(T.1-p)){Z.new[,j] <- Z[,j]/(rep(1,k*p)+rep(weight[i,],p))}
  #   Z.tilde[[i]] <- Z.new;
  # }
  # step.size.new <- rep(0,k)
  # for (i in 1:k){step.size.new[i] <- 1/( max(svd(Z.tilde[[i]])$d)  )^2; }
  # 
  
  if ( method == 'VAR'){
    ll <- 1;
    for ( ii in 1:k){
      l <- 2;
      while( l < max.iteration){
        l <- l+1;
        phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
        phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
        if ( max(abs(phi.new - phi.temp)) < tol) {break;}
        if (max(abs(phi.new - phi.temp)) > tol ) {
          phi.hat.fista[l,] <- phi.new; }
      }
      iter[ii,ll] <- l;
      phi.hat[ii,] <- phi.new;
    }
    ll.final <- 1;
  }
  
  # if ( method == 'VAR'){
  #   ll <- 1;
  #   for ( ii in 1:k){
  #     phi.temp <- matrix(0,1,k*p); l <- 0;
  #     while( l < max.iteration){
  #       l <- l+1; 
  #       phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
  #       if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
  #       if (max(abs(phi.new - phi.temp)) > tol ) {
  #         phi.temp <- phi.new; }
  #     }
  #     iter[ii,ll] <- l;
  #     phi.hat[ii,] <- phi.new;
  #   }
  #   ll.final <- 1;
  # }
  
  # else if (method == 'LASSO'){
  #   for (ll in 1:length(lambda)){
  #     for ( ii in 1:k){
  #       phi.temp <- matrix(0,1,k*p); l <- 0;
  #       while( l < max.iteration){
  #         l <- l+1; 
  #         # print(k)
  #         phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
  #         phi.new <- soft(phi.new,rep(0,k*p),lambda[ll]);
  #         if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
  #         if (max(abs(phi.new - phi.temp)) > tol ) {
  #           # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
  #           phi.temp <- phi.new; }
  #       }
  #       # print("l="); print(l)
  #       iter[ii,ll] <- l;
  #       phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
  #     }
  #     #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
  #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
  #     pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
  #   }
  #   ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
  #   phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
  # }
  
  else if (method == 'LASSO'){
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        l <- 2;
        while( l < max.iteration){
          l <- l+1; 
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft(phi.new,rep(0,k*p),lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.1-p)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,p+jjj,k,1)  )
      res <- t(data[(p+1):(T.1),])-forecast;
      temp <- AIC.BIC(res,phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)]);
      pred.error[ll] <- temp$BIC;
      # forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      # pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
  }
  
  else if (method == 'SLASSO'){
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        l <- 2;
        while( l < max.iteration){
          l <- l+1; 
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          # phi.new <- soft(phi.new,rep(weight[ii,],p)-rep(1,k*p),lambda[ll]);
          phi.new <- soft(phi.new,rep(weight[ii,],p),lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
  }
  
  else if (method == 'SSLASSO'){
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        l <- 2; Z.new <- Z.tilde[[ii]];
        while( l < max.iteration){
          l <- l+1; 
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size.new[ii]*t(Z.new%*%(t(Y[ii,]-phi.temp%*%Z.new)));
          # phi.new <- soft(phi.new,rep(weight[ii,],p)-rep(1,k*p),lambda[ll]);
          phi.new <- soft(phi.new,rep(0,k*p),lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new/(rep(1,k*p)+rep(weight[ii,],p));
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
  }
  
  else if (method == 'HVARC'){
    group.lag <- matrix(0,p,p*k);
    for ( l in 1:p){
        group.lag[l,(((p-l)*k+1):(p*k))] <- c((((p-l)*k+1):(p*k)));
    }
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        l <- 2;
        while( l < max.iteration){
          l <- l+1; 
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft.group(phi.new,rep(0,k*p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
  }
  
  else if (method == 'SHVARC'){
    group.lag <- matrix(0,p,p*k);
    for ( l in 1:p){
      group.lag[l,(((p-l)*k+1):(p*k))] <- c((((p-l)*k+1):(p*k)));
    }
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        l <- 2; 
        while( l < max.iteration){
          l <- l+1; 
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft.group(phi.new,rep(weight[ii,],p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
  }
  
  else if (method == 'SSHVARC'){
    group.lag <- matrix(0,p,p*k);
    for ( l in 1:p){
      group.lag[l,(((p-l)*k+1):(p*k))] <- c((((p-l)*k+1):(p*k)));
    }
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        l <- 2; Z.new <- Z.tilde[[ii]];
        while( l < max.iteration){
          l <- l+1;
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size.new[ii]*t(Z.new%*%(t(Y[ii,]-phi.temp%*%Z.new)));
          phi.new <- soft.group(phi.new,rep(0,k*p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new/(rep(1,k*p)+rep(weight[ii,],p));
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
  }
  
  else if (method == 'HVAROO'){
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        ## OO ###
        group.lag <- matrix(0,2*p,p*k);
        for ( l.new in 1:(2*p)){
          if ( (l.new %% 2) == 1 ){
            l.temp <- floor(l.new/2)+1;
            lag.temp <- c(((p-l.temp)*k+1):((p-l.temp+1)*k)); lag.temp <- lag.temp[-c(ii)]; group.lag[l.new,1:(k-1)] <- lag.temp;
            l.temp <- floor(l.new/2);
            if (l.temp != 0) {group.lag[l.new,(((p-l.temp)*k+1):(p*k))] <- c((((p-l.temp)*k+1):(p*k)));}
          }
          if ( (l.new %% 2) == 0 ) {l.temp <- l.new/2; group.lag[l.new,(((p-l.temp)*k+1):(p*k))] <- c((((p-l.temp)*k+1):(p*k)));}
        }
        l <- 2;
        while( l < max.iteration){
          l <- l+1; 
          # print(k)
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft.group(phi.new,rep(0,k*p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
    
  }
  
  else if (method == 'SHVAROO'){
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        ## OO ###
        group.lag <- matrix(0,2*p,p*k);
        for ( l.new in 1:(2*p)){
          if ( (l.new %% 2) == 1 ){
            l.temp <- floor(l.new/2)+1;
            lag.temp <- c(((p-l.temp)*k+1):((p-l.temp+1)*k)); lag.temp <- lag.temp[-c(ii)]; group.lag[l.new,1:(k-1)] <- lag.temp;
            l.temp <- floor(l.new/2);
            if (l.temp != 0) {group.lag[l.new,(((p-l.temp)*k+1):(p*k))] <- c((((p-l.temp)*k+1):(p*k)));}
          }
          if ( (l.new %% 2) == 0 ) {l.temp <- l.new/2; group.lag[l.new,(((p-l.temp)*k+1):(p*k))] <- c((((p-l.temp)*k+1):(p*k)));}
        }
        l <- 2;
        while( l < max.iteration){
          l <- l+1; 
          # print(k)
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft.group(phi.new,rep(weight[ii,],p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
    
  }
  
  else if (method == 'SSHVAROO'){
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        ## OO ###
        group.lag <- matrix(0,2*p,p*k);
        for ( l.new in 1:(2*p)){
          if ( (l.new %% 2) == 1 ){
            l.temp <- floor(l.new/2)+1;
            lag.temp <- c(((p-l.temp)*k+1):((p-l.temp+1)*k)); lag.temp <- lag.temp[-c(ii)]; group.lag[l.new,1:(k-1)] <- lag.temp;
            l.temp <- floor(l.new/2);
            if (l.temp != 0) {group.lag[l.new,(((p-l.temp)*k+1):(p*k))] <- c((((p-l.temp)*k+1):(p*k)));}
          }
          if ( (l.new %% 2) == 0 ) {l.temp <- l.new/2; group.lag[l.new,(((p-l.temp)*k+1):(p*k))] <- c((((p-l.temp)*k+1):(p*k)));}
        }
        l <- 2; Z.new <- Z.tilde[[ii]];
        while( l < max.iteration){
          l <- l+1; 
          # print(k)
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size.new[ii]*t(Z.new%*%(t(Y[ii,]-phi.temp%*%Z.new)));
          phi.new <- soft.group(phi.new,rep(0,k*p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new/(rep(1,k*p)+rep(weight[ii,],p));
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
    
  }
  
  else if (method == 'HVARELEM'){
    ## ELEM ###
    group.lag <- matrix(0,k*p,p);
    for(j in 1:k){
      for ( l in 1:p){
        temp <- j + seq(0,((p-1)*k),k);
        group.lag[((j-1)*p+l),(1:l)] <- temp[(p-l+1):(p)]
      }
    }
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        l <- 2;
        while( l < max.iteration){
          l <- l+1; 
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft.group(phi.new,rep(0,k*p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
    
  }
  
  else if (method == 'SHVARELEM'){
    ## ELEM ###
    group.lag <- matrix(0,k*p,p);
    for(j in 1:k){
      for ( l in 1:p){
        temp <- j + seq(0,((p-1)*k),k);
        group.lag[((j-1)*p+l),(1:l)] <- temp[(p-l+1):(p)]
      }
    }
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        l <- 2;
        while( l < max.iteration){
          l <- l+1; 
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft.group(phi.new,rep(weight[ii,],p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
    
  }
  
  else if (method == 'SSHVARELEM'){
    ## ELEM ###
    group.lag <- matrix(0,k*p,p);
    for(j in 1:k){
      for ( l in 1:p){
        temp <- j + seq(0,((p-1)*k),k);
        group.lag[((j-1)*p+l),(1:l)] <- temp[(p-l+1):(p)]
      }
    }
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        l <- 2; Z.new <- Z.tilde[[ii]];
        while( l < max.iteration){
          l <- l+1; 
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size.new[ii]*t(Z.new%*%(t(Y[ii,]-phi.temp%*%Z.new)));
          phi.new <- soft.group(phi.new,rep(0,k*p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new/(rep(1,k*p)+rep(weight[ii,],p));
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
    
  }
  
  else if (method == 'DHVAR'){
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        ## OO ###
        group.lag <- matrix(0,p*k,k);
        for ( l.new in 1:(p)){
          for(j1 in 1:k){
            lag.temp <- rank(weight[ii,], ties.method = "first"); 
            group.lag[((l.new-1)*k+j1),(1):(j1)] <- rep((l.new-1)*k,j1) + lag.temp[(k-j1+1):(k)]
          }
        }

        l <- 2;
        while( l < max.iteration){
          l <- l+1; 
          # print(k)
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft.group(phi.new,rep(0,k*p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
    
  }
  
  else if (method == 'DHVARC'){
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        ## OO ###
        lag.temp <- rank(weight[ii,], ties.method = "first"); 
        group.lag <- matrix(0,k,p*k);
          for(j1 in 1:k){
            for ( l.new in 1:(p)){
            group.lag[(j1),((l.new-1)*k+1):((l.new-1)*k+j1)] <- rep((l.new-1)*k,j1) + lag.temp[(k-j1+1):(k)]
          }
        }
        
        l <- 2;
        while( l < max.iteration){
          l <- l+1; 
          # print(k)
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft.group(phi.new,rep(0,k*p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
    
  }
  
  else if (method == 'SSDHVAR'){
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        ## OO ###
        group.lag <- matrix(0,p*k,k);
        for ( l.new in 1:(p)){
          for(j1 in 1:k){
            lag.temp <- rank(weight[ii,], ties.method = "first"); 
            group.lag[((l.new-1)*k+j1),(1):(j1)] <- rep((l.new-1)*k,j1) + lag.temp[(k-j1+1):(k)]
          }
        }
        
        l <- 2; Z.new <- Z.tilde[[ii]];
        while( l < max.iteration){
          l <- l+1; 
          # print(k)
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size.new[ii]*t(Z.new%*%(t(Y[ii,]-phi.temp%*%Z.new)));
          phi.new <- soft.group(phi.new,rep(0,k*p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new/(rep(1,k*p)+rep(weight[ii,],p));
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
    
  }
  
  else if (method == 'SSDHVARC'){
    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        ## OO ###
        lag.temp <- rank(weight[ii,], ties.method = "first"); 
        group.lag <- matrix(0,k,p*k);
        for(j1 in 1:k){
          for ( l.new in 1:(p)){
            group.lag[(j1),((l.new-1)*k+1):((l.new-1)*k+j1)] <- rep((l.new-1)*k,j1) + lag.temp[(k-j1+1):(k)]
          }
        }
        
        l <- 2; Z.new <- Z.tilde[[ii]];
        while( l < max.iteration){
          l <- l+1; 
          # print(k)
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size.new[ii]*t(Z.new%*%(t(Y[ii,]-phi.temp%*%Z.new)));
          phi.new <- soft.group(phi.new,rep(0,k*p),group.lag,lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp))); 
            phi.hat.fista[l,] <- phi.new; }
        }
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new/(rep(1,k*p)+rep(weight[ii,],p));
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
    
  }
  
  
  return(list(phi.hat = phi.hat, iter = iter, pred.error = pred.error, tune.final = lambda[ll.final]))
}

var.lasso.brk <- function(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = 10^(-4)){

  k <- length(data[1,]); T <- length(data[,1]); T.1 <- T;
  iter <- matrix(0,k,length(lambda));
  phi.hat <- matrix(0,k,k*p); phi.hat.fista <- matrix(0,max.iteration,k*p);
  pred.error <- rep(0,length(lambda)); phi.hat.temp <- matrix(0,k,k*p*length(lambda));
  Y <- as.matrix(t(data)); 
  Y <- as.matrix(Y[,(p+1):T.1]); 
  # Y <- Y%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  Z <- matrix(0,k*p,T.1-p); 
  # Z <- Z%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  for ( i in 1:(T.1-p)){
    for ( j in 1:p){
      Z[((j-1)*k+1):(j*k),i] <- t(data[i+p-j,])
    }
  }
  step.size <- 1/( max(svd(Z)$d)  )^2; 
  # Z.tilde <- vector("list",k);
  # for ( i in 1:k){
  #   Z.new <- matrix(0,k*p,T.1-p);
  #   for (j in 1:(T.1-p)){Z.new[,j] <- Z[,j]/(rep(1,k*p)+rep(weight[i,],p))}
  #   Z.tilde[[i]] <- Z.new;
  # }
  # step.size.new <- rep(0,k)
  # for (i in 1:k){step.size.new[i] <- 1/( max(svd(Z.tilde[[i]])$d)  )^2; }
  
  

  
  

    for (ll in 1:length(lambda)){
      for ( ii in 1:k){
        l <- 2;
        while( l < max.iteration){
          l <- l+1; 
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft(phi.new,rep(0,k*p),lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
          if (max(abs(phi.new - phi.temp)) > tol ) {
            phi.hat.fista[l,] <- phi.new; 
            # print(l);  
            }
        }
        
        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- matrix(0,k,T.1-p);
      forecast <- sapply(c((p):(T.1-1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,jjj,k,1)  )
      pred.error[ll] <- sum((t(data[(p+1):(T.1),])-forecast)^2) + lambda[ll]*sum(abs(phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)]));
    }
    ll.final <- 1
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
  
  


  
  
  return(list(phi.hat = phi.hat, iter = iter, pred.error = pred.error, tune.final = lambda[ll.final]))
}

var.break.fit <- function(method,data, weight = NULL, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3),initial.phi = NULL){
  method.full <- c("LASSO");
  if ( !(method %in% method.full) ){print("ERROR"); break; }
  
  k <- length(data[1,]); T <- length(data[,1]);
  iter <- matrix(0,k,length(lambda));
  n <- T - p;
  Y <- matrix(0,k*p,n);
  for( i in p:(T-1)){
    Y[,(i-p+1)] <- sapply(c(1:p), function(jjj) data[i-jjj+1,]  )
  }
  
  C <- vector("list",n);
  # for(jjj in 1:n){C[[jjj]] <- as.matrix(Y[,jjj],k*p,1)%*%t(as.matrix(data[jjj+p,],1,k)); }
  C <- lapply(c(1:n), function(jjj) as.matrix(Y[,jjj],k*p,1)%*%t(as.matrix(data[jjj+p,],1,k))    );
  C.sum <- matrix(0,k*p*n,k);
  C.sum[1:(k*p),] <- C[[1]];
  for(i in 2:n){C.sum[((i-1)*k*p+1):(i*k*p),] <- C.sum[((i-2)*k*p+1):((i-1)*k*p),] + C[[i]] }
  C.sum.new <- matrix(0,k*p*n,k);
  C.sum.new[1:(k*p),] <- C.sum[((n-1)*k*p+1):(n*k*p),];
  for(i in 2:n){C.sum.new[((i-1)*k*p+1):(i*k*p),] <- C.sum[((n-1)*k*p+1):(n*k*p),] - C.sum[((i-2)*k*p+1):((i-1)*k*p),] }
  
  D <- vector("list",n);
  D <- lapply(c(1:n), function(jjj) as.matrix(Y[,jjj],k*p,1)%*%t(as.matrix(Y[,jjj],k*p,1))    );
  D.sum <- matrix(0,k*p*n,k*p);
  D.sum[1:(k*p),] <- D[[1]];
  for(i in 2:n){D.sum[((i-1)*k*p+1):(i*k*p),] <- D.sum[((i-2)*k*p+1):((i-1)*k*p),] + D[[i]] }
  D.sum.new <- matrix(0,k*p*n,k*p);
  D.sum.new[1:(k*p),] <- D.sum[((n-1)*k*p+1):(n*k*p),];
  # D.sum.new[(k*p+1):(n*k*p),] <- sapply(c(2:n), function(jjj) D.sum[((n-1)*k*p+1):(n*k*p),] - D.sum[((jjj-2)*k*p+1):((jjj-1)*k*p),]    )
  for(i in 2:n){D.sum.new[((i-1)*k*p+1):(i*k*p),] <- D.sum[((n-1)*k*p+1):(n*k*p),] - D.sum[((i-2)*k*p+1):((i-1)*k*p),] }
  D.sum.new.inv <- matrix(0,k*p*n,k*p);
  for(jjj in 1:n){D.sum.new.inv[((jjj-1)*k*p+1):(jjj*k*p),] <- solve(D.sum.new[((jjj-1)*k*p+1):(jjj*k*p),] +  (tol)*diag(k*p) );   }
  # D.sum.new.inv <- sapply(c(1:n), function(jjj)  solve(D.sum.new[((jjj-1)*k*p+1):(jjj*k*p),])   )
  
  
  phi.hat <- matrix(0,k,k*p*n);
  if (!is.null(initial.phi)){phi.hat <- initial.phi;}
  active <- rep(0,n);
  active <- sapply(c(1:n), function(jjj) if ( sum((phi.hat[,((jjj-1)*k*p+1):(jjj*k*p)])^2) != 0  ) {jjj} else {0}   )
  active <- active[which(active!=0)]
  phi.new <- matrix(0,k,k*p*n);
  # phi.temp <- phi.hat;

  
  
  # X <- matrix(0,n,n*k*p);
  # for ( i in 1:n){
  #   for(j in 1:i){
  #     X[i,((j-1)*p*k+1):(j*p*k)] <-  t(as.matrix(Y[,i],k*p,1));
  #   }
  # }
  # step.size <- 1/( max(svd(X)$d)  )^2;
  # tol <- 0.5*step.size;
  # step.size <- (0.25)*tol;
  # step.size <- 0.01;
  
  # print(step.size)
  
  
  
  

  
   if (method == 'LASSO'){
    for (ll in 1:length(lambda)){
        l <- 2;

        while( l < max.iteration){
          l <- l+1; 
          phi.compare <- phi.hat;
          # cl <- makePSOCKcluster(2);
          # registerDoParallel(cl);
          # foreach (ii=1:n) %dopar%  {
          for (ii in 1:n){
            
              E <- vector("list",n);
              E <- lapply(c(1:n), function(jjj) D.sum.new[((max(jjj,ii)-1)*k*p+1):(max(jjj,ii)*k*p),]%*%t(phi.hat[,((jjj-1)*k*p+1):(jjj*k*p)])   )
              E <- Reduce("+",E);
              E <- E - D.sum.new[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%t(phi.hat[,((ii-1)*k*p+1):(ii*k*p)]);
            
            S <- C.sum.new[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]  - E;
              
              
            # B <- t(S)%*%D.sum.new[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%S
            # B <- sqrt(sum(B^2));
            # if( lambda[ll] >= B) { phi.temp <- matrix(0,k,k*p);   }
            # if( lambda[ll] < B) { phi.temp <- (1 - lambda[ll]/B)*D.sum.new.inv[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%S   }
              
              S <- soft.full(S,lambda[ll],k,p,n)
              phi.temp <- D.sum.new.inv[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%( S )
              phi.temp <- t(phi.temp);
              
              phi.hat[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;
              phi.new[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;
            
            
            # phi.temp <- D.sum.new.inv[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%( S )
            # phi.temp <- t(phi.temp);
            # phi.temp <- soft.full(phi.temp,lambda[ll],k,p,n)
            # phi.hat[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;
            # phi.new[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;
          }

          
          if ( max(abs(phi.new - phi.compare)) < tol) {
            print(max(abs(phi.new - phi.compare)));
            break;} 
          if (max(abs(phi.new - phi.compare)) > tol ) {
              phi.hat <- phi.new;    
              print(max(abs(phi.new - phi.compare)))
              }
        }
        # print("l="); print(l);

    }

  }
  
  # stopCluster(cl)
  
  # for(i in 1:(n*p*k)){
  #   for( j in 1:k){
  #     if ( abs (phi.hat[j,i] <= tol) ){phi.hat[j,i] <- 0;}
  #   }
  # }
  return(list(phi.hat = phi.hat, iter = iter))
}

break.var <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts ){
  k <- length(data[1,]); T <- length(data[,1]); m <- length(pts); L.n <- rep(0,m+1);
  if ( m == 0) { pts.temp <- c(1,T+1);}
  if ( m > 0){  pts.temp <- rep(0,m+2); pts.temp[1] <- 1; pts.temp[m+2] <- T+1; pts.temp[(2):(m+1)] <- pts;}
  for(mm in 1:( m+1 )){
    # print("mmmmm"); print(mm)
    data.temp <- data[(pts.temp[mm]):(pts.temp[(mm+1)]-1),];
    try <- var.lasso.brk(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = 10^(-4))
    L.n[mm] <- try$pred.error;
  }
  return(list(L.n = sum(L.n)))
}

backward.sel <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts){
  k <- length(data[1,]); T <- length(data[,1]);
  m <- length(pts); L.n <- rep(0,m); L.n.curr <- rep(0,1);
  
  ###### SHVAR FUN FULL PTS ##################
  try <- break.var(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts )
  L.n.curr <- try$L.n;
  
  for( mm in 1:m){
    pts.temp <- pts[-mm];
    
    ### SHVAR FIT FUNCTION ###################
    try <- break.var(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts.temp );
    L.n[mm] <- try$L.n;
    
    
  }
  return(list(L.n = L.n, L.n.curr = L.n.curr  ))
  
} 

second.step <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts, omega){
  m <- length(pts); if( m == 0){break;}
  mm <- 0;
  while(mm < m){
    mm <- mm + 1;
    try <- backward.sel(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts);
    L.n <- try$L.n; L.n.curr <- try$L.n.curr;
    if ( min(L.n) + (m-1)*omega >= L.n.curr + m*omega   ) {
      ic <- L.n.curr + (m - mm + 1)*omega;
      break;}
    if ( min(L.n) + (m-1)*omega < L.n.curr + m*omega   ){
      pts <- pts[-which(L.n == min(L.n))];
      print(pts);
      }
    
  }
  
  
  return(list(pts = pts, ic = ic ))
}

plot.new <- function (data, caltime = NULL){
  if (!is.matrix(data)) 
    data = as.matrix(data)
  if (is.ts(data)) {
    plot(data)
  }
  else {
    nT = dim(data)[1]
    tdx = c(1:nT)
    if (length(caltime) > 1) 
      tdx = caltime
    k = dim(data)[2]
    # if (k < 4) {
    #   par(mfcol = c(k, 1))
    #   for (j in 1:k) {
    #     plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
    #          type = "l")
    #   }
    # }

    if (k >= 1) {
      par(mfcol = c(1, 1))
      yl = range(data) * 1.05;
      # plot(tdx, data[, 1], xlab = "time", ylab = " ", type = "l", 
      #      ylim = yl)
      plot(tdx, data[, 1], xlab = "seconds", ylab = " ", type = "l", 
           ylim = yl, xaxt="n")
      axis(1, c(1,500,1000,1500,2000), c(0,50,100,150,200))
      axis(2)
      for (j in 2:k) {
        lines(tdx, data[, j], lty = j, col = j)
      }
    }
  }
  par(mfcol = c(1, 1))
}

MTSplot.new <- function (data, caltime = NULL){
  if (!is.matrix(data)) 
    data = as.matrix(data)
  if (is.ts(data)) {
    plot(data)
  }
  else {
    nT = dim(data)[1]
    tdx = c(1:nT)
    if (length(caltime) > 1) 
      tdx = caltime
    k = dim(data)[2]
    if (k < 0) {
      par(mfcol = c(k, 1))
      for (j in 1:k) {
        plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
             type = "l")
      }
    }
    if (k == 0) {
      par(mfcol = c(2, 2))
      for (j in 1:k) {
        plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
             type = "l")
      }
    }
    if ((k > 0) && (k < 1)) {
      par(mfcol = c(3, 2), mai = c(0.3, 0.3, 0.3, 0.3))
      k1 = 6
      jcnt = 0
      for (j in 1:k) {
        plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
             type = "l", cex.axis = 0.8)
        jcnt = jcnt + 1
        if ((jcnt == k1) && (k > 6)) {
          jcnt = 0
          cat("Hit return for more plots: ", "\n")
          readline()
        }
      }
    }
    if (k > 0) {
      par(mfcol = c(1, 1))
      yl = range(data) * 1.05
      plot(tdx, data[, 1], xlab = "time", ylab = " ", type = "l", 
           ylim = yl)
      for (j in 2:k) {
        lines(tdx, data[, j], lty = j, col = j)
      }
    }
  }
  par(mfcol = c(1, 1))
}

ar.est <- function(method, data, weight = NULL, lambda, p, break.pts, r.n, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3)){
  method.full <- c("LASSO");
  if ( !(method %in% method.full) ){print("ERROR"); break; }
  
  k <- length(data[1,]); T <- length(data[,1]); T.1 <- T; m.hat <- length(break.pts) + 1;
  ind.remain <- rep(0,(2+2*length(break.pts))); ind.remain[1] <- p; ind.remain[(2+2*length(break.pts))] <- T.1;
  for(i in 1:length(break.pts)){ind.remain[(2*i)] <- break.pts[i] - r.n - 1; ind.remain[(2*i+1)] <- break.pts[i] + r.n + 1;  }
  iter <- matrix(0,k,length(lambda));
  phi.hat <- matrix(0,k,k*m.hat*p); phi.hat.fista <- matrix(0,max.iteration,k*m.hat*p);
  pred.error <- rep(0,length(lambda)); phi.hat.temp <- matrix(0,k,k*m.hat*p*length(lambda)); std.res <- rep(0,length(lambda));
  
  Y <- as.matrix(t(data)); Y <- Y[,-seq(1,p,1)];
  Z <- matrix(0,k*p,T.1-p); 
  # Z <- Z%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  for ( i in 1:(T.1-p)){
    for ( j in 1:p){
      Z[((j-1)*k+1):(j*k),i] <- t(data[i+p-j,])
    }
  }
  
  for (i in 1:length(break.pts)){
    Y <- Y[,-seq(break.pts[i]-r.n,break.pts[i]+r.n,1)];
    # Z <- Z[,-seq(break.pts[i]-r.n,break.pts[i]+r.n,1)];
  }
  
  n <- length(Y[1,]);
  Z.new <- matrix(0,T.1-p,k*m.hat*p);
  Z.new[(1:(break.pts[1]-r.n-1)),1:(k*p)] <- t(Z[1:(k*p),(1:(break.pts[1]-r.n-1))]);
  if( m.hat > 2 ){
    for(i in 1:(m.hat-2)){
      # ind <- break.pts[i]-r.n-1;
      Z.new[((break.pts[i]+r.n+1):(break.pts[i+1]-r.n-1)),((i)*k*p+1):((i+1)*k*p)] <- t(Z[1:(k*p),((break.pts[i]+r.n+1):(break.pts[i+1]-r.n-1))]);
    }
  }

  Z.new[((break.pts[m.hat-1]+r.n+1):(T.1-p)),((m.hat-1)*k*p+1):((m.hat)*k*p)] <- t(Z[1:(k*p),((break.pts[m.hat-1]+r.n+1):(T.1-p))])
  del.ind <- c();
  for(i in 1:(T.1-p)){
    if (sum(Z.new[i,]^2) == 0){del.ind <- c(del.ind,i)}
  }
  Z.new <- Z.new[-del.ind,]
  Z <- t(Z.new);

  step.size <- 1/( max(svd(Z)$d)  )^2;
  # print(step.size)
  # step.size <- 10^(-3);
  
  
  for (ll in 1:length(lambda)){
    for ( ii in 1:k){
      l <- 2;
      while( l < max.iteration){
        l <- l+1; 
        phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
        phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
        phi.new <- soft(phi.new,rep(0,k*m.hat*p),lambda[ll]);
        if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
        if (max(abs(phi.new - phi.temp)) > tol ) {
          phi.hat.fista[l,] <- phi.new; 
          # print(l);  
        }
      }
      
      # print("l="); print(l)
      iter[ii,ll] <- l;
      phi.hat.temp[ii,((ll-1)*k*m.hat*p+1):(ll*k*m.hat*p)] <- phi.new;
    }
    
    forecast <- matrix(0,k,T.1);
    for(i in 1:m.hat){
      len <- ind.remain[(2*i)] - ind.remain[(2*i-1)];
      delay <- floor(2*len/3);
      delay <- 0;
      l.b <-  ind.remain[(2*i-1)] + delay +  1; u.b <- ind.remain[(2*i)];
      forecast[,(l.b):(u.b)] <- sapply(c((l.b-1):(u.b-1)), 
                                       function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*m.hat*p+(i-1)*k*p+1):((ll-1)*k*m.hat*p+(i)*k*p)],p,jjj,k,1)  )
    }
    
    if( ll == 1){
      del.ind <- c();
      for(i in 1:(T.1)){
        if (sum(forecast[,i]^2) == 0){del.ind <- c(del.ind,i)}
      }
    }
    
    # foo <- t(forecast);
    # plot(c(1:300),foo[,1],type="l",col="red")
    # lines(c(1:300),data[,1],col="green")
    
    residual <- t(data) - forecast;
    BIC.temp <- 0;
    for(i in 1:m.hat){
      l.b <-  ind.remain[(2*i-1)] + 0 +  1; u.b <- ind.remain[(2*i)];
      temp <- AIC.BIC(residual[,(l.b):(u.b)],phi.hat.temp[,((ll-1)*k*m.hat*p+(i-1)*k*p+1):((ll-1)*k*m.hat*p+(i)*k*p)]);
      BIC.temp <- BIC.temp + temp$BIC;
    }
    
    
    residual <- residual[,-del.ind]; nn <- length(residual[1,]); 
    # check <- AIC.BIC(residual,phi.hat.temp[,((ll-1)*k*m.hat*p+1):(ll*k*m.hat*p)]);
    # pred.error[ll] <- sum(residual^2)/(k*nn); std.res[ll] <- sd(residual[,]^2);
    # pred.error[ll] <- check$BIC;
    pred.error[ll] <- BIC.temp;
    # pred.error[ll] <- check$AIC;
  # print(ll)
  }
  ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
  # sd.final <- 1*median(std.res)/(sqrt(k*nn));
  # # sd.final <- 1*max(std.res)/(sqrt(k*nn));
  # # sd.final <- 3*sd(pred.error);
  # # sd.final <- quantile(pred.error[1:ll.final],0.975) + pred.error[ll.final];
  # # sd.final <- 1*median(std.res)/(sqrt(nn));
  # p.error <- abs(pred.error[ll.final:length(lambda)] - (pred.error[ll.final] + sd.final));
  # ll.final <- which(p.error==min(p.error)); ll.final <- min(ll.final);
  phi.hat <- phi.hat.temp[,((ll.final-1)*k*m.hat*p+1):(ll.final*k*m.hat*p)];
  
  # for(i in 1:k){
  #   for(j in 1:(k*m.hat*p)){
  #     if(abs(phi.hat[i,j]) < lambda[ll.final]  ){phi.hat[i,j] <- 0;}
  #   }
  # }
  
  return(list(phi.hat = phi.hat, iter = iter, pred.error = pred.error, tune.final = lambda[ll.final]))
  
  
  
}

mspe.plot <- function(method,pred.error,lambda,tune.final,jj){
  plot( lambda, pred.error, type = 'o', col = "blue", main=c(jj,method))
  abline(v=tune.final)
}

estimation.check <- function(phi,phi.hat){
  k <- length(phi[,1]); p <- (length(phi[1,]))/k; N <- (length(phi.hat[1,]))/(k*p); L <- matrix(0,k,k);
  l2.error <- rep(0,N); true.zero <- rep(0,N); true.non.zero <- rep(0,N); L.hat <- matrix(0,k,k*N); L.hat.final <- rep(0,N);
  false.zero <- rep(0,N); false.non.zero <- rep(0,N); 
  
  count.non <- 0; count.zero <- 0;
  for (i in 1:k){
    for (j in 1:(k*p)){
      if ( phi[i,j] != 0  ){count.non <- count.non + 1;}
      if ( phi[i,j] == 0  ){count.zero <- count.zero + 1;}
    }
  }
  for (i in 1:k){
    for (j in 1:k){
      temp <- 0; 
      for (l in 1:p){if ( phi[i,((l-1)*k+j)] !=0   ) {temp <- temp+1;}}
      L[i,j] <- temp;
    }
  }
  for (jj in 1:N){
    phi.temp <- phi.hat[,((jj-1)*k*p+1):(jj*k*p)];
    count.false.zero <- 0; count.false.non.zero <- 0; count.true.non.zero <- 0; count.true.zero <- 0;
    for (i in 1:k){
      for (j in 1:(k*p)){
        if ( phi[i,j] != 0 && phi.hat[i,((jj-1)*k*p+j)] == 0   ){count.false.zero <- count.false.zero + 1;}
        if ( phi[i,j] == 0 && phi.hat[i,((jj-1)*k*p+j)] != 0   ){count.false.non.zero <- count.false.non.zero + 1;}
        if ( phi[i,j] == 0 && phi.hat[i,((jj-1)*k*p+j)] == 0   ){count.true.zero <- count.true.zero + 1;}
        if ( phi[i,j] != 0 && phi.hat[i,((jj-1)*k*p+j)] != 0   ){count.true.non.zero <- count.true.non.zero + 1;}
        if ( phi[i,j] == 0 ){phi.temp[i,j] <- 0;}
      }
    }
    l2.error[jj] <- sum((phi.temp-phi)^2);
    true.zero[jj] <- count.true.zero; true.non.zero[jj] <- count.true.non.zero;
    false.zero[jj] <- count.false.zero; false.non.zero[jj] <- count.false.non.zero;
  }
  for (jj in 1:N){
    for (i in 1:k){
      for (j in 1:k){
        temp <- 0; 
        for (l in 1:p){if ( phi.hat[i,(((jj-1)*k*p)+((l-1)*k)+j)] !=0   ) {temp <- temp+1;}}
        L.hat[i,(((jj-1)*k)+j)] <- temp;
      }
    }
    
  }
  L.hat.final <- sapply(c(1:N), function(jj) (sum((L.hat[,((jj-1)*k+1):(jj*k)]-L)^2))/(sum(L))  )
  
  return( list(l2.error.mean = mean(l2.error), l2.error.sd = sd(l2.error), true.zero.median = median(true.zero),
               true.non.zero.median = median(true.non.zero), false.zero.median = median(false.zero),
               false.non.zero.median = median(false.non.zero), lag.error.mean = mean(L.hat.final),
               lag.error.sd = sd(L.hat.final)) )
}

estimation.check.new <- function(phi,phi.final,k,p,m.hat,pts.final.full){
  N <- length(phi.final); 
  l2.error <- rep(0,N); true.zero <- rep(0,N); true.non.zero <- rep(0,N); 
  false.zero <- rep(0,N); false.non.zero <- rep(0,N); 
  
  count.non <- 0; count.zero <- 0;
  for (i in 1:k){
    for (j in 1:(k*m.hat*p)){
      if ( phi[i,j] != 0  ){count.non <- count.non + 1;}
      if ( phi[i,j] == 0  ){count.zero <- count.zero + 1;}
    }
  }

  for (jj in 1:N){
    phi.temp <- phi.final[[jj]]; pt.temp <- pts.final.full[[jj]];
    if( length(phi.temp[1,]) < k*m.hat*p  ){
      phi.temp.new <- matrix(0,k,k*m.hat*p); phi.temp.new[,(1):(k*p)] <- phi.temp[,(1):(k*p)]; ind <- 0;
      for(i in 2:m.hat){
        if( pt.temp[i-1] == 0  ){ind <- ind+1; phi.temp.new[,((i-1)*k*p+1):(i*k*p)] <- phi.temp.new[,((i-2)*k*p+1):((i-1)*k*p)];    }
        if( pt.temp[i-1] != 0  ){phi.temp.new[,((i-1)*k*p+1):(i*k*p)] <- phi.temp[,((i-1-ind)*k*p+1):((i-ind)*k*p)];    }
      }
      phi.temp <- phi.temp.new;
    }
    
    count.false.zero <- 0; count.false.non.zero <- 0; count.true.non.zero <- 0; count.true.zero <- 0;
    for (i in 1:k){
      for (j in 1:(k*m.hat*p)){
        if ( phi[i,j] != 0 && phi.temp[i,j] == 0   ){count.false.zero <- count.false.zero + 1;}
        if ( phi[i,j] == 0 && phi.temp[i,j] != 0   ){count.false.non.zero <- count.false.non.zero + 1;}
        if ( phi[i,j] == 0 && phi.temp[i,j] == 0   ){count.true.zero <- count.true.zero + 1;}
        if ( phi[i,j] != 0 && phi.temp[i,j] != 0   ){count.true.non.zero <- count.true.non.zero + 1;}
        if ( phi[i,j] == 0 ){phi.temp[i,j] <- 0;}
      }
    }
    l2.error[jj] <- sqrt(sum((phi.temp-phi)^2)/sum(phi^2));
    true.zero[jj] <- count.true.zero; true.non.zero[jj] <- count.true.non.zero;
    false.zero[jj] <- count.false.zero; false.non.zero[jj] <- count.false.non.zero;
  }

  
  return( list(l2.error.mean = mean(l2.error), l2.error.sd = sd(l2.error), true.zero.median = median(true.zero),
               true.non.zero.median = median(true.non.zero), false.zero.median = median(false.zero),
               false.non.zero.median = median(false.non.zero)) )
}

AIC.BIC <- function(residual,phi){
  k <- length(phi[,1]); k.lam <- length(phi[1,]); T.new <- length(residual[1,]); count <- 0;
  for (i in 1:k){for (j in 1:k.lam){if(phi[i,j] != 0){count <- count + 1;}}}
  # for ( i in 1:T.new){residual[,i] <- residual[,i] - as.matrix(rowMeans(residual))}
  # sigma.hat <- (1/(T.new-1))*(residual%*%t(residual)) + 10^(-10)*diag(k);
  sigma.hat <- 0*diag(k);
  for(i in 1:T.new){sigma.hat <- sigma.hat +  residual[,i]%*%t(residual[,i]);  }
  sigma.hat <- (1/(T.new))*sigma.hat;
  ee.temp <- min(eigen(sigma.hat)$values); 
  if(ee.temp <= 0){sigma.hat <- sigma.hat + (2.0)*(abs(ee.temp)+10^(-4))*diag(k);}
  # sigma.hat <- (1/(T.new-1))*(t(residual)%*%(residual)) + 0*10^(-10)*diag(T.new);
  log.det <- log(det(sigma.hat)+ 0*10^(-10)); 
  # print(det(sigma.hat));
  return(list(AIC = log.det + 2*count/T.new, BIC = log.det + log(T.new)*count/T.new))
}

AIC.BIC.CV <- function(residual,phi){
  k <- length(phi[,1]); k.lam <- length(phi[1,]); T.new <- length(residual[1,]); count <- 0;
  for (i in 1:k){for (j in 1:k.lam){if(phi[i,j] != 0){count <- count + 1;}}}
  for ( i in 1:T.new){residual[,i] <- residual[,i] - as.matrix(rowMeans(residual))}
  sigma.hat <- (1/(T.new-1))*(residual%*%t(residual)) 
  log.det <- log(det(sigma.hat)); if (abs(det(sigma.hat)) < 10^(-8) ){log.det <- -8^(1)}
  return(list(AIC = log.det + 2*count/T.new, BIC = log.det + log(T.new)*count/T.new))
}


plot.matrix <- function (phi,p,name = NULL) {
  B <- phi
  if (nrow(B) == 1) {
    B <- matrix(B[, 1:ncol(B)], nrow = 1)
  }
  else {
    B <- B[, 1:ncol(B)]
  }
  k <- nrow(B)
  s1 <- 0
  m <- 0
  s <- 0
  s <- s + s1
  text <- c()
  for (i in 1:p) {
    text1 <- as.expression(bquote(bold(Phi)^(.(i))))
    text <- append(text, text1)
  }
  if (m > 0) {
    for (i in (p + 1):(p + s + 1)) {
      text1 <- as.expression(bquote(bold(beta)^(.(i - p - 
                                                    s1))))
      text <- append(text, text1)
    }
  }
  f <- function(m) t(m)[, nrow(m):1]
  rgb.palette <- colorRampPalette(c("white", "blue"), space = "Lab")
  at <- seq(k/2 + 0.5, p * (k) + 0.5, by = k)
  if (m > 0) {
    at2 <- seq(p * k + m/2 + 0.5, p * k + s * m + 0.5, by = m)
  }
  else {
    at2 = c()
  }
  at <- c(at, at2)
  se2 = seq(1.75, by = k, length = k)
  L2 <- levelplot(as.matrix(f(B)), col.regions = rgb.palette, 
                  colorkey = NULL, xlab = NULL, ylab = NULL, main = list(label = name, 
                                                                         cex = 1), panel = function(...) {
                                                                           panel.levelplot(...)
                                                                           panel.abline(a = NULL, b = 1, h = seq(1.5, m * s + 
                                                                                                                   p * k + 0.5, by = 1), v = seq(1.5, by = 1, length = p * 
                                                                                                                                                   k + m * s), lwd = 0.5)
                                                                           bl1 <- seq(k + 0.5, p * k + 0.5, by = k)
                                                                           b23 <- seq(p * k + 0.5, p * k + 0.5 + s * m, by = m)
                                                                           b1 <- c(bl1, b23)
                                                                           panel.abline(a = NULL, b = 1, v = p * k + 0.5, lwd = 3)
                                                                           panel.abline(a = NULL, b = 1, v = b1, lwd = 2)
                                                                         }, scales = list(x = list(alternating = 1, labels = text, 
                                                                                                   cex = 1, at = at, tck = c(0, 0)), y = list(alternating = 0, 
                                                                                                                                              tck = c(0, 0))))
  return(L2)
}



first.step <- function(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4)){
  
  test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  phi.hat.full <- test$phi.hat;
  ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  temp.lam <- ll;
  # temp.lam <- quantile(phi.hat.full,ll)
  ll <- c(0);
  brk.points.list <- vector("list",length(ll));
  
  for(j in 1:length(ll)){
    
    phi.hat <- phi.hat.full;
    # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
    
    n <- T - p;
    m.hat <- 0; brk.points <- rep(0,n);
    
    for (i in 2:n)
    {
      if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
        m.hat <- m.hat + 1; brk.points[m.hat] <- i;
      }
    }
    
    loc <- rep(0,m.hat);
    brk.points <- brk.points[1:m.hat];
    brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
    m.hat <- length(brk.points);
    if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
    loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
    brk.points.list[[j]] <- brk.points;
  }
  
  
  return(list(brk.points = brk.points.list, phi.hat = phi.hat.full))
}


first.step.cv <- function(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4)){
  
  
  kk <- length(lambda);
  cv <- rep(0,kk);
  phi.final <- vector("list",kk);
  T <- length(data.temp[,1]); k <- length(data.temp[1,]);
  brk.points.final <- vector("list",kk);
  
  for (i in 1:kk) {
    test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;
    
    # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
    # temp.lam <- ll;
    # temp.lam <- quantile(phi.hat.full,ll)
    ll <- c(0);
    brk.points.list <- vector("list",length(ll));
    
    for(j in 1:length(ll)){
      
      phi.hat <- phi.hat.full;
      # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
      
      n <- T - p;
      m.hat <- 0; brk.points <- rep(0,n);
      
      for (iii in 2:n)
      {
        if ( sum((phi.hat[,((iii-1)*k*p+1):(iii*k*p)] )^2 ) != 0   ){
          m.hat <- m.hat + 1; brk.points[m.hat] <- iii;
        }
      }
      
      loc <- rep(0,m.hat);
      brk.points <- brk.points[1:m.hat];
      brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
      m.hat <- length(brk.points);
      if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
      loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      brk.points.list[[j]] <- brk.points;
    }
    
    brk.points.final[[i]] <- brk.points;
    m.hat <- length(brk.points);
    
    
    
    phi.full.all <- vector("list",T-p);
    phi.temp.cv <- matrix(0,k,k*p);
    forecast <- matrix(0,k,T-p);
    for(j in (p+1):T){
      phi.temp.cv <- phi.temp.cv + phi.hat.full[,((j-p-1)*k*p+1):((j-p)*k*p)];
      phi.full.all[[(j-p)]] <- phi.temp.cv;
      forecast[,(j-p)] <- pred(t(data),phi.temp.cv,p,j-1,k,1)
    }
    
    phi.hat.temp <- matrix(0,k,(k*p*(m.hat+1)));
    if( m.hat > 1){
      for(jj in 1:m.hat){
        phi.hat.temp[,((jj-1)*k*p+1):((jj)*k*p)] <- phi.full.all[[(brk.points[jj]-p-1)]];
      }
    }
    phi.hat.temp[,((m.hat)*k*p+1):((m.hat+1)*k*p)] <- phi.full.all[[(T-p)]];
    
    residual <- t(data[(p+1):T,]) - forecast;
    BIC.temp <- 0;
    if (m.hat == 0){
      temp <- AIC.BIC.CV(residual[,(1):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
      BIC.temp <- BIC.temp + temp$BIC;
    }
    if(m.hat >=1){
      temp <- AIC.BIC.CV(residual[,(1):(brk.points[1]-1)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
      BIC.temp <- BIC.temp + temp$BIC;
      if ( m.hat >= 2){
        for(ii in 2:m.hat){
          l.b <-  brk.points[(ii-1)]; u.b <- brk.points[ii]-1;
          temp <- AIC.BIC.CV(residual[,(l.b):(u.b)],phi.hat.temp[,((1-1)*k*m.hat*p+(ii-1)*k*p+1):((1-1)*k*m.hat*p+(ii)*k*p)]);
          BIC.temp <- BIC.temp + temp$BIC;
        }
      }
      temp <- AIC.BIC.CV(residual[,(brk.points[m.hat]):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(m.hat)*k*p+1):((1-1)*k*m.hat*p+(m.hat+1)*k*p)]);
      BIC.temp <- BIC.temp + temp$BIC;
    }
    
    
    cv[i] <- sum( (forecast - t(data[(p+1):T,])  )^2 );
    # cv[i] <- BIC.temp;
    # temp <- AIC.BIC.CV(residual,phi.hat.temp);
    # cv[i] <- temp$BIC;
    print("====================================")
  }
  
  lll <- min(which(cv==min(cv)));
  phi.hat.full <- phi.final[[lll]];
  
  
  # test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  # phi.hat.full <- test$phi.hat;
  # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  # temp.lam <- ll;
  # # temp.lam <- quantile(phi.hat.full,ll)
  # ll <- c(0);
  # brk.points.list <- vector("list",length(ll));
  # 
  # for(j in 1:length(ll)){
  #   
  #   phi.hat <- phi.hat.full;
  #   # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
  #   
  #   n <- T - p;
  #   m.hat <- 0; brk.points <- rep(0,n);
  #   
  #   for (i in 2:n)
  #   {
  #     if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
  #       m.hat <- m.hat + 1; brk.points[m.hat] <- i;
  #     }
  #   }
  #   
  #   loc <- rep(0,m.hat);
  #   brk.points <- brk.points[1:m.hat];
  #   brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
  #   m.hat <- length(brk.points);
  #   if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
  #   loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
  #   brk.points.list[[j]] <- brk.points;
  # }
  
  
  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv.final = lambda[lll]))
}



first.step.cv.new <- function(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4),cv.index){
  
  cv.l <- length(cv.index); data.org <- data.temp; T.org <- length(data.temp[,1]); k.org <- length(data.temp[1,]);
  data.temp <- data.temp[-cv.index,];
  kk <- length(lambda);
  cv <- rep(0,kk); cv.var <- rep(0,kk);
  phi.final <- vector("list",kk);
  T <- length(data.temp[,1]); k <- length(data.temp[1,]);
  brk.points.final <- vector("list",kk);
  
  for (i in 1:kk) {
    if ( i == 1){
      test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
    }
    if ( i > 1 ){
      initial.phi <- phi.final[[(i-1)]]
      test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = initial.phi)
    }
    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;
    
    # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
    # temp.lam <- ll;
    # temp.lam <- quantile(phi.hat.full,ll)
    ll <- c(0);
    brk.points.list <- vector("list",length(ll));
    
    for(j in 1:length(ll)){
      
      phi.hat <- phi.hat.full;
      # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
      
      n <- T - p;
      m.hat <- 0; brk.points <- rep(0,n);
      
      for (iii in 2:n)
      {
        if ( sum((phi.hat[,((iii-1)*k*p+1):(iii*k*p)] )^2 ) != 0   ){
          m.hat <- m.hat + 1; brk.points[m.hat] <- iii;
        }
      }
      
      loc <- rep(0,m.hat);
      brk.points <- brk.points[1:m.hat];
      brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
      m.hat <- length(brk.points);
      if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
      loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      brk.points.list[[j]] <- brk.points;
    }
    
    brk.points.final[[i]] <- brk.points;
    m.hat <- length(brk.points);
    
    
    
    phi.full.all <- vector("list",T-p);
    phi.temp.cv <- matrix(0,k,k*p);
    forecast <- matrix(0,k,T-p); forecast.new <- matrix(0,k,cv.l);
    for(j in (p+1):T){
      phi.temp.cv <- phi.temp.cv + phi.hat.full[,((j-p-1)*k*p+1):((j-p)*k*p)];
      phi.full.all[[(j-p)]] <- phi.temp.cv;
      forecast[,(j-p)] <- pred(t(data.temp),phi.temp.cv,p,j-1,k,1)
    }
    
    for(j in (1):cv.l){
      forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j]-1-p-j+1)]],p,cv.index[j]-1,k,1)
    }
    
    forecast.all <- matrix(0,k,T.org-p); forecast.all[,cv.index] <- forecast.new; forecast.all[,-cv.index] <- forecast;
    residual <- (forecast.all - t(data.org[((p+1):T.org),]))^2;
    var.matrix <-  matrix(0,k,m.hat+1);
    if ( m.hat == 0){var.matrix <- sapply(c(1:k), function(jjj)  var(residual[jjj,])  );}
    if ( m.hat >= 1){
      var.matrix[,1] <- t(as.vector(sapply(c(1:k), function(jjj)  var(residual[jjj,(1):(brk.points[1]-1)])  )));
      var.matrix[,(m.hat+1)] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(m.hat)]):(T.org-p)])  );
      }
    if ( m.hat >=2 ){
      for(mm in 2:m.hat){
        var.matrix[,mm] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(mm-1)]):(brk.points[mm]-1)])  );
      }
    }
    
    
    if ( m.hat == 0) {cv.var[i] <- sum(var.matrix);}
    if ( m.hat >= 1){
      for (i.1 in 1:cv.l){
        ind <- 0;
        for (i.2 in 1:m.hat){
          if (cv.index[i.1] > brk.points[i.2]){ind <- ind + 1;}
        }
        cv.var[i] <- cv.var[i] + sum(var.matrix[,(ind+1)]);
      }
    }
    
    cv.var[i] <- (1/(k*cv.l))*sqrt(cv.var[i]);
    
    # phi.hat.temp <- matrix(0,k,(k*p*(m.hat+1)));
    # if( m.hat > 1){
    #   for(jj in 1:m.hat){
    #     phi.hat.temp[,((jj-1)*k*p+1):((jj)*k*p)] <- phi.full.all[[(brk.points[jj]-p-1)]];
    #   }
    # }
    # phi.hat.temp[,((m.hat)*k*p+1):((m.hat+1)*k*p)] <- phi.full.all[[(T-p)]];
    # 
    # residual <- t(data.temp[(p+1):T,]) - forecast;
    # BIC.temp <- 0;
    # if (m.hat == 0){
    #   temp <- AIC.BIC.CV(residual[,(1):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    # }
    # if(m.hat >=1){
    #   temp <- AIC.BIC.CV(residual[,(1):(brk.points[1]-1)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    #   if ( m.hat >= 2){
    #     for(ii in 2:m.hat){
    #       l.b <-  brk.points[(ii-1)]; u.b <- brk.points[ii]-1;
    #       temp <- AIC.BIC.CV(residual[,(l.b):(u.b)],phi.hat.temp[,((1-1)*k*m.hat*p+(ii-1)*k*p+1):((1-1)*k*m.hat*p+(ii)*k*p)]);
    #       BIC.temp <- BIC.temp + temp$BIC;
    #     }
    #   }
    #   temp <- AIC.BIC.CV(residual[,(brk.points[m.hat]):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(m.hat)*k*p+1):((1-1)*k*m.hat*p+(m.hat+1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    # }
    
    
    cv[i] <- (1/(k*cv.l))*sum( (forecast.new - t(data.org[cv.index,])  )^2 );
    # cv[i] <- BIC.temp;
    # temp <- AIC.BIC.CV(residual,phi.hat.temp);
    # cv[i] <- temp$BIC;
    print("====================================")
  }
  
  lll <- min(which(cv==min(cv)));
  ind.new <- 0;
  if (lll < kk){
    for(i.3 in (lll+1):(kk)){
      if ( cv[i.3] < (cv[lll] + cv.var[lll] )   ){ind.new <- ind.new + 1;}
    }
  }
  lll <- lll + ind.new;
  phi.hat.full <- phi.final[[lll]];
  print("CV"); print(cv);
  print("CV.VAR"); print(cv.var);
  
  # test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  # phi.hat.full <- test$phi.hat;
  # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  # temp.lam <- ll;
  # # temp.lam <- quantile(phi.hat.full,ll)
  # ll <- c(0);
  # brk.points.list <- vector("list",length(ll));
  # 
  # for(j in 1:length(ll)){
  #   
  #   phi.hat <- phi.hat.full;
  #   # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
  #   
  #   n <- T - p;
  #   m.hat <- 0; brk.points <- rep(0,n);
  #   
  #   for (i in 2:n)
  #   {
  #     if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
  #       m.hat <- m.hat + 1; brk.points[m.hat] <- i;
  #     }
  #   }
  #   
  #   loc <- rep(0,m.hat);
  #   brk.points <- brk.points[1:m.hat];
  #   brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
  #   m.hat <- length(brk.points);
  #   if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
  #   loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
  #   brk.points.list[[j]] <- brk.points;
  # }
  
  
  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv.final = lambda[lll]))
}



generate.phi.new <- function(k,p.t,l.bound,max.iter,seed){
  
  for ( try in 1:max.iter){
    phi <- matrix(0,k,k*p.t)
    set.seed(try+seed)
    base <- matrix(2*(runif((k^2)*p.t)-1/2),k,k*p.t);
    for(i in 1:k){
      for(j in 1:(k*p.t)){
        d <- (abs(abs(i-j)-1)+1);
        if ( abs(base[i,j]/d) > l.bound   ){phi[i,j] <- base[i,j]/d; }
        if ( abs(base[i,j]/d) <= l.bound   ){phi[i,j] <- 0; }
      }
    }
    
    companion.phi <- matrix(0,k*p.t,k*p.t);
    companion.phi[1:k,] <- phi;
    if( p.t > 1){
      for(i in 1:(p.t-1)){companion.phi[(i*k+1):((i+1)*k),((i-1)*k+1):((i)*k)] <- diag(k)       }
    }
    aaa <- eigen(companion.phi)$values; aaa <- Mod(aaa);
    # print("TRY=="); print(try)
    if(max(aaa) < 1){break;}
  }
  
  return(list(phi = phi, try = try))
}






