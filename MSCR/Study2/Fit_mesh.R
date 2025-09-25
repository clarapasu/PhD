

# Let's assume that the dataset comes in continuous time like we simulated, i.e. for each individual
# we have a list of capture times and traps where the animal was seen. 
# The discretize functions creates a discetization of the time from 0 to T in r segments.
# Segments in which the individual is unseen gets the value 0, while segments where the individual is
# observed get cut in two at the event time. The segment until the event time gets the value of the trap where 
# it is seen, and so on. 
# This function does this for 1 individual. 

#In this file the Likelihood function takes "mesh" as a vector of meshes (although it should work for one mesh only too)
#which makes it work in the blue/red mesh setting 



seen<-function(T,trap,s,theta){
  h<-sum(apply(trap,1,halfnormal,theta=theta,s=s))
  U<-1-exp(-T*h)
  return(U)
}


discretize<-function(df,T,r){
  t<-round(seq(0, T, T/r),digits=3)
  y<-rep(0,r+1)
  for (j in 1:length(df$Time)){
    time<-round(df$Time[j],digits=3)
    i<-findInterval(time, t)
    if (round(time,digits=2)==0.3){  # This is because of the 0.3 issue in computers, still weird.... 
      y[i+1]<-df$y[j]
    }
    else if (round(time,digits=2)==0.6){  # Also weird
      y[i+1]<-df$y[j]
    }
    else if (round(time,digits=2)==0.7){  # I don't understand this
      y[i+1]<-df$y[j]
    }
    else if(round(time,digits=3)==round(t[i],digits=3)){
      y[i]<-df$y[j]
    }
    else{
      t<-append(t,time,after=i)
      y<-append(y,df$y[j],after=i)
    }
  }
  data<-data.frame(t,y)
  return(data)
}


Likelihood_ind<-function(theta,trap,data,s,memory){
  m<-dim(data)[1]
  L<-1
  t<-data$t
  y<-data$y
  seen<-0
  capt_time<-t[1]
  start_hazards<-sum(apply(trap,1,halfnormal,theta=theta,s=unlist(s)))
  first_capt<-which(data$y!=0)[1]
  capt_time<-t[first_capt]
  z<-trap[y[first_capt],]
  L<-exp(-capt_time*start_hazards)*halfnormal(z,theta,s)
  j<-first_capt
  for (i in (first_capt+1):m){
    L<-L*Survival(trap,theta,t,i,j,z,s,memory)
    if (y[i]!=0){
      k<-trap[y[i],]
      L<-L*hazard(k,theta,(t[i]-capt_time),z,s,memory)
      z<-k
      capt_time<-t[i]
      j<-i
    }
  }
  #seen_prob<-seen(t[m],trap,s,theta)
  return(L)
} 



Likelihood_test<-function(theta,trap,df,ac,T,r,memory){
  n<-length(unique(df$id))
  L<-0 
  for (i in 1:n){ 
    data<-discretize(df[df$id==i,],T,r)
    seen_prob<-seen(T,trap,ac[i,],theta)
    L<-L+log(Likelihood_ind(theta,trap,data,ac[i,],memory))-log(seen_prob)
    # This is not working anymore, it did when I divided by seen_prob inside Likelihood_ind
    # so I am a bit confused
  }
  return(-L)
}



Likelihood_integrate<-function(theta,trap,data,mesh,memory){
  Likeli<-0
  A<-0
  for (mask in mesh){
    if (dim(mask)[1]!=0){
      a = attr(mask,"a") # grid cell area
      D<-dim(mask)[1]
      A = A+a*D
      L<-0
      mask<-matrix(c(mask[,1],mask[,2]),ncol=2,nrow=D)
      #mask<-as.matrix(mask)
      for (i in 1:D){
        L<-L+Likelihood_ind(theta,trap,data,mask[i,],memory)
      }

      Likeli<-Likeli+(L*a)
    }
    }
    
  return(log(Likeli/A))
}




#Likelihood_integrate_test<-function(theta,trap,data,mesh,memory){
#  a = attr(mesh,"a") # grid cell area
#  D<-dim(mesh)[1]
#  A = a*D
#  mesh<-matrix(c(mesh[,1],mesh[,2]),ncol=2,nrow=D)
#  L<-0
#  for (i in 1:D){
#    L<-L+Likelihood_ind(theta,trap,data,mesh[i,],memory)
#  }
#  return(c(L*a,A))
#}




#Seen_int<-function(T,trap,mesh,theta){
#  a = attr(mesh,"a") # grid cell area
#  D<-dim(mesh)[1]
#  A = a*D
#  mesh<-matrix(c(mesh[,1],mesh[,2]),ncol=2,nrow=D)
#  S<-0
#  for (i in 1:D){
#    S<-S+seen(T,trap,as.double(mesh[i,]),theta)
#  }
#  return(log(S*a/A))
#}


Seen_int<-function(T,trap,mesh,theta){
  Integral<-0
  A<-0
  if (is.vector(mesh)==FALSE){
    mesh<-list(mesh)
  }
  for (mask in mesh){
    a = attr(mask,"a") # grid cell area
    D<-dim(mask)[1]
    A = A +a*D
    S<-0
    mask<-matrix(c(mask[,1],mask[,2]),ncol=2,nrow=D)
    for (i in 1:D){
      S<-S+seen(T,trap,as.double(mask[i,]),theta)
    }
    Integral<-Integral+(S*a)
  }
  return(log(Integral/A))
}
#   

#mesh can be a mesh in the normal sense OR it can be a list of different grids covering the area


Likelihood<-function(theta,trap,df,mesh,mask,T,r,memory){
  n<-length(unique(df$id))
  L<-0
  for (i in 1:n){ 
    data<-df[df$id==i,]
    L<-L+Likelihood_integrate(theta,trap,data,mesh[[i]],memory)
  }
  L<-L-n*Seen_int(T,trap,mask,theta)
  return(-L)
}



# Likelihood<-function(theta,trap,df,mesh,T,r,memory){
#   n<-length(unique(df$id))
#   L<-0
#   if (is.vector(mesh)==FALSE){
#     mesh<-list(mesh)
#   }
#   for (i in 1:n){ 
#     data<-discretize(df[df$id==i,],T,r)
#     L<-L+Likelihood_integrate(theta,trap,data,mesh,memory)
#   }
#   L<-L-n*Seen_int(T,trap,mesh,theta)
#   #L<-L-Seen_int(T,trap,mesh,theta)
#   return(-L)
# }



Full_likelihood<-function(par,n,trap,df,mesh,T,r,memory){
  N<-exp(par[1])+n
  theta<-par[2:4]
  p<-exp(Seen_int(T,trap,mesh,theta))
  L<-log(factorial(N))-log(factorial(N-n))+n*log(p)+(N-n)*log(1-p) # use gamma function instead of factorial 
  Lik<-0
  for (i in 1:n){ 
    data<-discretize(df[df$id==i,],T,r)
    Lik<-Lik+Likelihood_integrate(theta,trap,data,mesh,memory)
  }
  L<-L+Lik-n*log(p)
  return(-L)
}




Likelihood_onefunc<-function(theta,trap,df,mesh,T,memory){
  
  n<-length(unique(df$id))
  
  A = (max(mesh[,1]) - min(mesh[,1])) * (max(mesh[,2]) - min(mesh[,2]))
  a = A / nrow(mesh) # grid cell area
  D <- nrow(mesh)
  
  S<-0
  start_hazards <- c()
  for (i in 1:D){
    start_hazards[i] <- 0
    for(j in 1:nrow(trap)){
      start_hazards[i] <- start_hazards[i] + halfnormal(trap[j,],theta=theta,s=mesh[i,])
    }
    U <- 1 - exp(-T* start_hazards[i])
    S <- S + U
  }
  LS <- log(S*a/A)
  
  LL <- 0 
  # sums log-likelihood over individuals
  for (h in 1:n){
    data <- df[df$id == h, 1:2]
    m <- dim(data)[1]
    
    L <- 0 
    # cat("L is ", LL, "\n")
    # sums up likelihood over possible AC locations (integrating out latent ACs)
    for (j in 1:D){
      t<-data$t
      y<-data$y
      s <- mesh[j, ]
      first_capt<-which(data$y!=0)[1]
      capt_time<-t[first_capt]
      z<-trap[y[first_capt],]
      Lj<-exp(-capt_time*start_hazards[j])*halfnormal(z,theta,s)
      seen_ind <- first_capt
      
      for (i in (first_capt+1):m){
        total_hazard <- 0
        for(jj in 1:nrow(trap)){
          #total_hazard <- total_hazard + hazard(trap[jj,], theta=theta, (t[i] - t[i-1])/2, z, s, memory)
          total_hazard <- total_hazard + hazard(trap[jj,], theta=theta, t[i-1]-t[seen_ind]+((t[i]-t[i-1])/2), z, s, memory)
        }
        survival <- exp(-(t[i]-t[i-1]) * total_hazard)
        Lj <- Lj * survival
        if (y[i]!=0){
          k<-trap[y[i],]
          Lj <- Lj * hazard(k,theta,(t[i]-capt_time),z,s,memory)
          z<-k
          capt_time<-t[i]
          seen_ind <- i
        }
        
        
      }
      L <- L + Lj
    }
    LL <- LL + log(L * a/A)
  }
  
  LL <- LL - n * LS
  return(-LL)
}


confint_param<-function(par,hessian){
  Hessian<-hessian
  cov<-solve(Hessian)
  se<-sqrt(diag(cov)) 
  upper<-par+1.96*se
  lower<-par-1.96*se
  interval<-data.frame(value=par, lower=lower,upper=upper)
  return(interval)
}

confint_pop<-function(theta_est,Hessian,T,trap,mesh,n){
  cov<-solve(Hessian)
  Seen_est<-exp(Seen_int(T,trap,mesh,theta_est)) # the estimated probability of being observed
  Pop_size<-function(theta,T,trap,mesh,n,Seen_est){
    N_est<-n/Seen_est # Horvitz-Thompson estimator of the population size
    return(N_est)
  }
  N_est<-Pop_size(theta_est,T,trap,mesh,n,Seen_est=Seen_est)
  d<- grad(Pop_size,theta_est, T=T,trap=trap,mesh=mesh,n=n,Seen_est=Seen_est)
  Var<- d%*% cov %*%d + n*(1-Seen_est)/Seen_est
  upper<-N_est+1.96*sqrt(Var)
  lower<-N_est-1.96*sqrt(Var)
  interval<-data.frame(N_est,lower,upper)
  return(interval)
}
