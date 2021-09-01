####################################
# Needed packages:
# library(statmod)
# library(mvtnorm)
# library(Matrix)
# library(splines)

################################
#
# Y    =  Response vector encoded 0 and 1
# X    =  Design matrix corresponding to Y
# T    =  Vector of times at which the observations in Y were taken
# U    =  Vector which identifies user IDs for the obsevrations in Y, X, and T
# Tmax =  Max time to be included in analysis
# knot.seq = Interior knot set to be used to model subject specific trajectories
# Iter =  Numnber of MCMC iterates

Mixed.SVM<-function(Y, X, T, U, Tmax, knot.seq=c(.5), Iter){

  ############################################################
  # Structure the data
  user<-unique(U)    # Identify unique subjects
  J<-length(user)

  n.l<-lapply(1:J, function(j){ id<-which(U==user[j]); 
                              time<-T[id]/Tmax ;
                              ct<-time<=1;
                              time<-time[ct==1];
                              length(time)
                            })

  Z.l<-lapply(1:J, function(j){ id<-which(U==user[j]); 
                              time<-T[id]/Tmax ;
                              ct<-time<=1;
                              time<-time[ct==1];  
                              bs(time, knots=knot.seq, Boundary.knots = c(0,1))   
                            })

  X.l<-lapply(1:J, function(j){ id<-which(U==user[j]); 
                              time<-T[id]/Tmax ;
                              ct<-time<=1;
                              time<-time[ct==1];  
                              Xt<-X[id,];
                              Xt<-Xt[ct==1,]                                                       
                              cbind(Xt)
                            })
     
  Y.l<-lapply(1:J, function(j){ id<-which(U==user[j]); 
                              Yt<-Y[id]; 
                              time<-T[id]/Tmax ;
                              ct<-time<=1;
                              idy<-which(Yt==0); 
                              Yt[idy]<--1;
                              Yt<-Yt[ct==1] 
                              return(Yt)
                            })

  ##################################################
  # Structure design matrices specifically for SVM 

  Xt.l<-lapply(1:J,function(j){X.l[[j]]*as.vector(Y.l[[j]])})
  Zt.l<-lapply(1:J,function(j){Z.l[[j]]*as.vector(Y.l[[j]])})

  ##################################################################################################################
  ##################################################################################################################
  # Buidling the fixed effects design matrix 

  Xt<-NULL
  X<-NULL
  Y<-NULL
  for(j in 1:J){
    Xt<-rbind(Xt,Xt.l[[j]])
    X<-rbind(X,X.l[[j]])
    Y<-c(Y,Y.l[[j]])
  }

  ###########################################
  # Compute once 
  P<-dim(X.l[[1]])[2]    # Number of fixed effects covariates
  Q<-dim(Z.l[[1]])[2]    # Number of random effects

  #######################################
  # Prior parameters
  Di<-penalty.create(knots=knot.seq,Boundary.knots=c(0,1))   #Smoothing penalty   

  ###############################################
  # Initializations
  sig2<-apply(X,2,var)   
  sig2[1]<-1
  beta<-SVM.fit(Y,X,h=0.1,nu=10,alpha=2,sig2)
  gamma.l<-lapply(1:J,function(i,Q){rep(0,Q)},Q=Q)
  T.i<-diag(rep(1,P))
  lamj<-diag(rep(1,P))
  alpha<-1
  eta<-1


  nu<-10
  tau2.i<-0.0000001

  ######################################
  # MCMC setup 
  beta.save<-matrix(-99,nrow=Iter,ncol=P)
  gamma.save<-array(-99,c(Iter,Q,J))

  monitor<-rep(-99, Iter)

  pb <- txtProgressBar(0,Iter,initial=0,style=3)
  
  for(g in 1:Iter){
    setTxtProgressBar(pb,g)
    ############################################################################################
    # Sample and structure the lambda
    lam.i.l<-lapply(X=1:J, lam.i.samp, Xt.l=Xt.l,Zt.l=Zt.l,n.l=n.l,beta=beta, gamma.l=gamma.l)
    lam.i<-unlist(lam.i.l,use.names = FALSE)


    ############################################################################################
    # Sample beta
    ZtG.l<-lapply(1:J,function(j){Zt.l[[j]]%*%gamma.l[[j]]})
    ZtG<-unlist(ZtG.l,use.names = FALSE)


    Bi<-T.i+t(Xt*lam.i)%*%Xt
    B<-solve(Bi)

    b<-B%*%t(Xt)%*%(1+lam.i-lam.i*ZtG)
    beta<-as.vector(rmvnorm(1, mean = b, sigma = B))

    ############################################################################################
    # Sample omgj

    mu.omgi<- sqrt(lamj^2/(beta^2))
    omg.i<-rinvgauss(n=P, mean=mu.omgi, shape=lamj^2)
    T.i<-diag(omg.i)

    ############################################################################################
    # Sample lamj

    lamj<-rgamma(P,  (alpha+1),(abs(beta)+eta)      )


    ################################################
    ######## Sample alpha (hyper parameter) ########

    a.seq<-seq(0.001,0.999,.01)

    pmat<-log(matrix((1+abs(beta)/eta),ncol=P,nrow=length(a.seq),byrow=TRUE))
    pmat<-pmat*(-1/a.seq)
    pvec1<-apply(pmat,1,sum)
    pvec2<-  P*log((1-a.seq)/a.seq)
    
    const<- -max(pvec1 +pvec2)
    wk<-exp(pvec1 + pvec2 +const)
    swk<-sum(wk)
    wk<-wk/swk
    a<-sample(a.seq[wk!=0],1,prob=(wk[wk!=0]))
    alpha<-1/a-1
    
    ##############################################
    ######## Sample eta (hyper parameter) ########

    e.seq<-seq(0.001,0.999,.1)

    dmat<-matrix(abs(beta),ncol=P,nrow=length(e.seq),byrow=TRUE)
    dmat<- -(alpha+1)*log(1+dmat*(e.seq/(1-e.seq)))
    dvec1<-apply(dmat,1,sum)
    dvec2<-  P*log((e.seq)/(1-e.seq))

    const<- -max(dvec1 +dvec2)
    wk<-exp(dvec1 + dvec2 +const)
    swk<-sum(wk)
    wk<-wk/swk
    e<-sample(e.seq[wk!=0],1,prob=(wk[wk!=0]))
    eta<-1/e-1
   

    ############################################################################################
    # Sample gamma
    gamma.l<-lapply(X=1:J, gamma.samp, Xt.l=Xt.l, Zt.l=Zt.l, beta=beta, lam.i.l=lam.i.l, Di=Di, tau2.i=tau2.i)

    ############################################################################################
    # Sample tau2

    gamma2.l<-lapply(1:J,function(j,Di){t(gamma.l[[j]])%*%Di%*%gamma.l[[j]]}, Di=Di)
    app<-1 + (J*Q)/2 
    bpp<-1 + sum(unlist(gamma2.l))/2
    tau2.i<-rgamma(1,shape=app,rate=bpp)

    beta.save[g,]<-beta

    for(j in 1:J){
      gamma.save[g,,j]<-gamma.l[[j]]
    }

  
    monitor[g] <- sum(abs(Y-unlist(lapply(1:J,function(j){ifelse( (X.l[[j]]%*%beta+Z.l[[j]]%*%gamma.l[[j]])>0,1,-1)})) ))

    close(pb)
  }

  return(list("Beta"=beta.save,"Gamma"=gamma.save,"Misclassified"=monitor))
}



#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#
# Support functions 
#

#############################################################
# Continuous approx to the hinge loss
Lh<-function(g,h){
  res<- 0+ (1-g)*((1-g)>h) + ((1+h-g)^2/(4*h))*(abs(1-g)<= h) 
  return(res)
}

#############################################################
# SVM objective function

SVM.obj<-function(beta,Y,X,h,nu,alpha,sig){
  gi<-Y*(X%*%beta)
  res<-sum(Lh(gi,h)) + nu^(-alpha)*sum(abs(beta/sig)^alpha)
  return(res)
}


##############################################################
# Completes model fitting for standard SVM  (Note: Look for continuous/differentiable approx. of |\cdot| for future implementation)
SVM.fit<-function(Y,X,h,nu,alpha,sig2){
  res<-optim(rep(0,dim(X)[2]),SVM.obj,Y=Y,X=X,h=h,nu=nu,alpha=alpha,sig=sqrt(sig2))$par
  return(res)
}


##############################################################
# Functions to do sampling via lapply

# Sample lambda inverse 
lam.i.samp<-function(j,Xt.l,Zt.l,n.l,beta,gamma.l){
  mu.lami<-1/abs(1-Xt.l[[j]]%*%beta-Zt.l[[j]]%*%gamma.l[[j]])
	lam.i<-rinvgauss(n=n.l[[j]], mean=mu.lami, shape=1)
  return(lam.i)
}

# Sample gamma
gamma.samp<-function(j, Xt.l, Zt.l, beta, lam.i.l, Di, tau2.i){
	Si.g<-tau2.i*Di+t(Zt.l[[j]]*lam.i.l[[j]])%*%Zt.l[[j]]
	S.g<-solve(Si.g)
	mu.g<-S.g %*%t(Zt.l[[j]])%*%(1+lam.i.l[[j]]-lam.i.l[[j]]*(Xt.l[[j]]%*%beta))
	gamma.s<-as.vector(rmvnorm(1, mean = mu.g, sigma = S.g))
  return(gamma.s)
}


##############################################################
# Constructs smoothing penalty matrix
penalty.create<-function(knots,Boundary.knots){
  h<-0.001
  t<-seq(0,1,by=h)
  m<-length(t)
  Bf<-bs(t,knots=knots,Boundary.knots=Boundary.knots)
  dBf<-(Bf[-1,]-Bf[-m,])/h
  m<-m-1
  d2Bf<-(dBf[-1,]-dBf[-m,])/h
  m<-m-1
  L<-dim(Bf)[2]
  R<-matrix(-99,nrow=L,ncol=L)
  for(i in 1:L){
    for(j in 1:L){
      temp<-d2Bf[,i]*d2Bf[,j]
      R[i,j]<-sum((temp[-1]+temp[-m])/2*h)  # Need to change .001 based on spacing
    }
  }
  return(R)
}





#################################################################
# Regression coefficient summary
# 
# MCMC  = output from Mixed.SVM
# burn  = burn in to use
# names = covariate names 

Reg.summary<-function(MCMC,burn,names,knot.seq){

  lower<-burn
  upper<-dim(MCMC$Beta)[1]

  s1<-apply(MCMC$Beta[lower:upper,],2,mean)
  s2<-apply(MCMC$Beta[lower:upper,],2,sd)
  s3<-apply(MCMC$Beta[lower:upper,],2,quantile,prob=0.025)
  s4<-apply(MCMC$Beta[lower:upper,],2,quantile,prob=0.975)

  s<-round(cbind(s1,s3,s4),4)
  s<-cbind(s,ifelse( (s3*s4>0),"significant","not"))
  s<-cbind(names,s)
  colnames(s)<-c("Variable","Point Est.", "Lower:CI95","Upper:CI95", "")
  return(s)
  #print(s[-1:(length(knot.seq)+3)])
}

#################################################################
# Plotting population trajectory

Pop.trajectory<-function(MCMC,burn,knot.seq,Tmax){
  lower<-burn
  upper<-dim(MCMC$Beta)[1]

  tp<-seq(0,1,length.out=500)
  Xp<-bs(tp, knots=knot.seq)
  Q<-dim(Xp)[2]

  f.est<-matrix(-99,nrow=length(lower:upper),ncol=length(tp))
  tick<-1
  for(i in lower:upper){
    f.est[tick,]<-Xp%*%MCMC$Beta[i,1:Q]
    tick<-tick+1
  }

  plot(tp*Tmax, apply(f.est,2,mean), xlab="Days",ylab="f(Days)",ylim=c(-.5,.25),type="l")
  lines(tp*Tmax, apply(f.est,2,quantile,prob=0.25),col="red")
  lines(tp*Tmax, apply(f.est,2,quantile,prob=0.75),col="red")
  lines(tp*Tmax, apply(f.est,2,quantile,prob=0.025),col="green")
  lines(tp*Tmax, apply(f.est,2,quantile,prob=0.975),col="green")

  abline(0,0,col="blue")

}


#################################################################
# Functions to plot subject specific trajectory

Ind.trajectory<-function(subj, X, T, U, MCMC, burn, knot.seq, Tmax){
    user<-unique(U)    # Identify unique subjects
    lower<-burn
    upper<-dim(MCMC$Beta)[1]

    sel<-which(user==subj)

    id<-which(U==subj)
    time<-T[id]/Tmax ;
    ct<-time<=1;
    time<-time[ct==1];  
    Xt<-X[id,];
    Xt<-Xt[ct==1,]
    B<-bs(time, knots=knot.seq, Boundary.knots = c(0,1))                           
    F<-cbind(Xt)

    f.est<-matrix(-99,nrow=length(lower:upper),ncol=length(time))
    tick<-1
    for(i in lower:upper){
        f.est[tick,]<-F%*%MCMC$Beta[i,] + B%*%MCMC$Gamma[i,,sel]
        tick<-tick+1
    }
    data <- data.frame(x=time*Tmax, y= apply(f.est,2,mean))
    fig <- ggplot(data, aes(x=x,y=y)) + geom_point(size=1) + theme(text=element_text(size=11, family="Courier")) + labs(title=paste0(subj)) + xlab("Days") + ylab("f(Days)")

    return(fig)
}



# color version
Ind.trajectory.col<-function(subj, X, T, U, MCMC, burn, knot.seq, Tmax, col.sel){
    user<-unique(U)    # Identify unique subjects
    lower<-burn
    upper<-dim(MCMC$Beta)[1]

    sel<-which(user==subj)

    id<-which(U==subj)
    time<-T[id]/Tmax ;
    ct<-time<=1;
    time<-time[ct==1];  
    Xt<-X[id,];
    Xt<-Xt[ct==1,]
    B<-bs(time, knots=knot.seq, Boundary.knots = c(0,1))                           
    F<-cbind(Xt)

    f.est<-matrix(-99,nrow=length(lower:upper),ncol=length(time))
    tick<-1
    for(i in lower:upper){
        f.est[tick,]<-F%*%MCMC$Beta[i,] + B%*%MCMC$Gamma[i,,sel]
        tick<-tick+1
    }
    data <- data.frame(x=time*Tmax, y= apply(f.est,2,mean))
    #plot(time*Tmax, apply(f.est,2,mean), xlab="Days",ylab="f(Days)",col=(Xt[,col.sel]+1))
    #abline(0,0,col="blue")
    fig <- plot_ly(data, type="scatter",mode="markers",x=~x,y=~y, color=(Xt[,col.sel]+1))
    fig <- fig %>% layout(
        title = paste0(subj),
        xaxis = list(title = "Days"),
        yaxis = list(title = "f(Days)"),
        showlegend = FALSE
    )

    fig <- fig %>% hide_colorbar()
    
    return(fig)
}


# Pretty version, specific to time-varying covariates in ABQ DrinQ
Ind.trajectory.pretty<-function(subj, X, T, U, MCMC, burn, knot.seq, Tmax){
    user<-unique(U)    # Identify unique subjects
    lower<-burn
    upper<-dim(MCMC$Beta)[1]

    sel<-which(user==subj)

    id<-which(U==subj)
    time<-T[id]/Tmax ;
    ct<-time<=1;
    time<-time[ct==1];  
    Xt<-X[id,];
    Xt<-Xt[ct==1,]
    B<-bs(time, knots=knot.seq, Boundary.knots = c(0,1))                           
    F<-cbind(Xt)

    f.est<-matrix(-99,nrow=length(lower:upper),ncol=length(time))
    tick<-1
    for(i in lower:upper){
        f.est[tick,]<-F%*%MCMC$Beta[i,] + B%*%MCMC$Gamma[i,,sel]
        tick<-tick+1
    }
    data <- data.frame(x=time*Tmax, y= apply(f.est,2,mean))

    #heatcolor <- heat.colors(8)
    heatcolor  <- rainbow(8)   

    col <- array(dim=nrow(Xt))
    shape <- array(dim=nrow(Xt))
    for (i in 1:nrow(Xt)) {
        col[i] <- heatcolor[1]
        #shape[i] <- 18
        if(Xt[i,"THCUSEDAY"]==1 & Xt[i,"NICUSEDAY"]==0 & Xt[i,"OTHERDRUGUSE"]==0) col[i] <- heatcolor[2]
        if(Xt[i,"THCUSEDAY"]==0 & Xt[i,"NICUSEDAY"]==1 & Xt[i,"OTHERDRUGUSE"]==0) col[i] <- heatcolor[3]
        if(Xt[i,"THCUSEDAY"]==0 & Xt[i,"NICUSEDAY"]==0 & Xt[i,"OTHERDRUGUSE"]==1) col[i] <- heatcolor[4]
        if(Xt[i,"THCUSEDAY"]==1 & Xt[i,"NICUSEDAY"]==1 & Xt[i,"OTHERDRUGUSE"]==0) col[i] <- heatcolor[5]
        if(Xt[i,"THCUSEDAY"]==1 & Xt[i,"NICUSEDAY"]==0 & Xt[i,"OTHERDRUGUSE"]==1) col[i] <- heatcolor[6]
        if(Xt[i,"THCUSEDAY"]==0 & Xt[i,"NICUSEDAY"]==1 & Xt[i,"OTHERDRUGUSE"]==1) col[i] <- heatcolor[7]
        if(Xt[i,"THCUSEDAY"]==1 & Xt[i,"NICUSEDAY"]==1 & Xt[i,"OTHERDRUGUSE"]==1) col[i] <- heatcolor[8]
    }
                  

    fig <- ggplot(data, aes(x=x,y=y,col=col)) + geom_point(size=1) + 
    theme(text=element_text(size=11, family="Courier"),legend.position="bottom") + labs(title=paste0(subj)) + xlab("Days") + 
    ylab("f(Days)") + #xlim(0,700) + ylim(-3,3) +
    scale_colour_discrete(name  ="Substance Use",
      breaks=c('#FF0000FF','#FFBF00FF','#80FF00FF','#00FF40FF','#00FFFFFF','#0040FFFF','#8000FFFF','#FF00BFFF'),
        labels=c("None", "THC", "NICOTINE", "OTHER DRUGS", "THC + NICOTINE", "THC + OTHER DRUGS", "NICOTINE + OTHER DRUGS", "All"))                    
    return(fig)
}



############################################
# Standard drinking summaries
log.summaries <- function(iid, timeframe,visitday,includevisitday) {
    beginday <- visitday-timeframe

    results <- data.frame()

    for( i in 1:length(iid)) {

        cur <- subset(drinks,ursi==iid[i]) # pull for current subject
        cur <- cur[order(cur$Day),] # sort by date
        if(includevisitday==FALSE) {
            cur <- subset(cur,Day >= beginday & Day < visitday)
        } else {
            cur <- subset(cur,Day > beginday & Day <= visitday)
        }
        pdd <- NA
        phdd <- NA
        ddd <- NA
        dpd <- NA
        totaldrinks <- NA
        totaldays <- NA
        totaldrinkdays <- NA
        totalheavydays <- NA
        phwd <- NA
        mxd <- NA
        
        if (nrow(cur)>0) {
          pdd <- sum(cur$AD)/nrow(cur) # percent drinking days
          phdd <- sum(cur$HD)/nrow(cur) # percent heavy drinking days

          if (sum(cur$AD)==0) {ddd <- 0} else {ddd <- sum(cur$SD)/sum(cur$AD)} # ave std. drinks per drinking day
          dpd <- sum(cur$SD)/nrow(cur) # ave. std. drinks per day
          totaldrinks <- sum(cur$SD) # total std. drinks consumed
          totaldays <- nrow(cur) # total days
          totaldrinkdays <- sum(cur$AD) # total number of drinking days
          totalheavydays <- sum(cur$HD) # total number of heavy drinking days
          if (sum(cur$AD)==0) {phwd <- 0} else {phwd <- sum(cur$HD)/sum(cur$AD)} # percent heavy drinking days while drinking
          mxd <- max(cur$SD,na.rm=TRUE) # maximum standard drinks
        }
        df <- data.frame(usri=iid[i],pdd, phdd, ddd, dpd, totaldrinks, totaldays, totaldrinkdays, totalheavydays, phwd, mxd)
        
        results <- rbind(results,df)
    }
    return(results)
}

#######################################################################
# Computes model performance metrics for the Bayesian mixed effects SVM
Pred.summary.MSVM<-function(MCMC.res,Y, X, T, U, Tmax, knot.seq=c(.5), Iter, burn.in){

  beta<-MCMC.res$Beta
  gamma<-MCMC.res$Gamma


  ############################################################
  # Structure the data
  user<-unique(U)    # Identify unique subjects
  J<-length(user)

  n.l<-lapply(1:J, function(j){ id<-which(U==user[j]); 
                              time<-T[id]/Tmax ;
                              ct<-time<=1;
                              time<-time[ct==1];
                              length(time)
                            })

  Z.l<-lapply(1:J, function(j){ id<-which(U==user[j]); 
                              time<-T[id]/Tmax ;
                              ct<-time<=1;
                              time<-time[ct==1];  
                              bs(time, knots=knot.seq, Boundary.knots = c(0,1))   
                            })

  X.l<-lapply(1:J, function(j){ id<-which(U==user[j]); 
                              time<-T[id]/Tmax ;
                              ct<-time<=1;
                              time<-time[ct==1];  
                              Xt<-X[id,];
                              Xt<-Xt[ct==1,]                                                       
                              cbind(Xt)
                            })
     
  Y.l<-lapply(1:J, function(j){ id<-which(U==user[j]); 
                              Yt<-Y[id]; 
                              time<-T[id]/Tmax ;
                              ct<-time<=1;
                              idy<-which(Yt==0); 
                              Yt[idy]<--1;
                              Yt<-Yt[ct==1] 
                              return(Yt)
                            })

  ##################################################
  # Structure design matrices specifically for SVM 

  Xt.l<-lapply(1:J,function(j){X.l[[j]]*as.vector(Y.l[[j]])})
  Zt.l<-lapply(1:J,function(j){Z.l[[j]]*as.vector(Y.l[[j]])})

  ##################################################################################################################
  ##################################################################################################################
  # Buidling the fixed effects design matrix 

  Xt<-NULL
  X<-NULL
  Y<-NULL

  for(j in 1:J){
    Xt<-rbind(Xt,Xt.l[[j]])
    X<-rbind(X,X.l[[j]])
    Y<-c(Y,Y.l[[j]])
  }


  TP<-rep(-99,Iter)
  FP<-rep(-99,Iter)
  TN<-rep(-99,Iter)
  FN<-rep(-99,Iter)

  for(g in 1:Iter){
  
    ZG.l<-lapply(1:J,function(j){Z.l[[j]]%*%gamma[g,,j]})
    ZG<-unlist(ZG.l,use.names = FALSE)
    XB<-X%*%beta[g,] 
    Yp<-XB+ZG >=0
    Yp[Yp==0]<- -1 

    TP[g]<- mean(  Yp[Y==1]==Y[Y==1]   )
    FP[g]<- mean(  Yp[Y==-1]!=Y[Y==-1]   )
    TN[g]<- mean(  Yp[Y==-1]==Y[Y==-1]   )
    FN[g]<- mean(  Yp[Y==1]!=Y[Y==1]   )
  
  }

  TPm <- mean(TP[burn.in:Iter])
  FPm <- mean(FP[burn.in:Iter])
  TNm <- mean(TN[burn.in:Iter])
  FNm <- mean(FN[burn.in:Iter])


  return(list("TP"=TPm,"FP"=FPm,"TN"=TNm,"FN"=FNm))
}



#######################################################################
# Computes model performance metrics for the Bayesian fixed effects SVM
# i.e., standard SVM
Pred.summary.SVM<-function(MCMC.resY, X, T, U, Sim=10000, burn.in){

  beta<-MCMC.res$Beta

  ############################################################
  # Structure the data
  user<-unique(U)    # Identify unique subjects
  J<-length(user)

  X.l<-lapply(1:J, function(j){ id<-which(U==user[j]); 
                              time<-T[id]/Tmax ;
                              ct<-time<=1;
                              time<-time[ct==1];  
                              Xt<-X[id,];
                              Xt<-Xt[ct==1,]                                                       
                              cbind(Xt)
                            })
     
  Y.l<-lapply(1:J, function(j){ id<-which(U==user[j]); 
                              Yt<-Y[id]; 
                              time<-T[id]/Tmax ;
                              ct<-time<=1;
                              idy<-which(Yt==0); 
                              Yt[idy]<--1;
                              Yt<-Yt[ct==1] 
                              return(Yt)
                            })


  ##################################################################################################################
  ##################################################################################################################
  # Buidling the fixed effects design matrix 

 
  X<-NULL
  Y<-NULL
  for(j in 1:J){
    X<-rbind(X,X.l[[j]])
    Y<-c(Y,Y.l[[j]])
  }



  TP<-rep(-99,Iter)
  FP<-rep(-99,Iter)
  TN<-rep(-99,Iter)
  FN<-rep(-99,Iter)

  for(g in 1:Iter){
  
    XB<-X%*%beta[g,] 
    Yp<-XB >=0
    Yp[Yp==0]<- -1 

    TP[g]<- mean(  Yp[Y==1]==Y[Y==1]   )
    FP[g]<- mean(  Yp[Y==-1]!=Y[Y==-1]   )
    TN[g]<- mean(  Yp[Y==-1]==Y[Y==-1]   )
    FN[g]<- mean(  Yp[Y==1]!=Y[Y==1]   )
  
  }

  mean(TP[burn.in:Iter])
  mean(FP[burn.in:Iter])
  mean(TN[burn.in:Iter])
  mean(FN[burn.in:Iter])


  return(list("TP"=TP,"FP"=FP,"TN"=TN,"FN"=FN))
}

########################################################################
# Functions to plot subject specific trajectory with credible 
# intervals around functon values
Ind.trajectory.CI<-function(subj, T, U, MCMC, burn, knot.seq, Tmax){
    user<-unique(U)    # Identify unique subjects
    lower<-burn
    upper<-dim(MCMC$Beta)[1]

    sel<-which(user==subj)

    id<-which(U==subj)
    time<-T[id]/Tmax ;
    ct<-time<=1;
    time<-time[ct==1];  
 
    B<-bs(time, knots=knot.seq, Boundary.knots = c(0,1))                           
    
    f.est<-matrix(-99,nrow=length(lower:upper),ncol=length(time))
    tick<-1
    for(i in lower:upper){
        f.est[tick,]<- B%*%MCMC$Gamma[i,,sel]
        tick<-tick+1
    }
    fhat<-apply(f.est,2,mean)
    fl<-apply(f.est,2, quantile, prob=0.025)
    fu<-apply(f.est,2, quantile, prob=0.975)
    fmax<-max(fu)
    fmin<-min(fl)
   
    data <- data.frame(x=time*Tmax, fhat,fu,fl)
    fig <- ggplot(data, aes(x=x,y=fhat)) + geom_line() + geom_hline(yintercept=0,color="red") + xlab("Days") + ylab("f(Days)") + geom_line(aes(x=x,y=fu),linetype="dashed",color="blue") + geom_line(aes(x=x,fl), linetype="dashed", color="blue") + labs(title=paste0(subj)) 
    
    return(fig)
}



