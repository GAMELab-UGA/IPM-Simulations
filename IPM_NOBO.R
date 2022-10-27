######################################################################
## IPM of simulated radiotelemetry and CMR data from NOBO
#######################################################################


library(rjags)



##########################################################################
########################################################################
## Part 1. Parameters
#######################################################################
#########################################################################


## State Parameters
##############################################################################

# Using two age classes (juvenile, adult) and females only

# Assuming 100% site fidelity (no immigration or emigration)

# Assuming two periods within the annual cycle (winter/non-breeding) 
# Capture/telemetry occurs in fall (October - early November) and 
#   spring (late February - March)(Wann et al. 2020). Corresponds to post-breeding, pre-hunting
#   (fall) and post-hunting, pre-breeding (spring).
# Winter period is 5 months while breeding period is 7.
# Juveniles become adults in first spring.

N.init.a <- 200 # Initial adult female abundance
N.init.j <- 470 # Initial juvenile female abundance

FemFecundity <- 1.85 # Per bird female fecundity from Yeiser et al. 2021 (not-trapped)
Malebrooding <- 0.28 # Rate of male brooding (Yeiser et al. 2021, Sandercock et al. 2008)
relFec <- FemFecundity + FemFecundity*Malebrooding # Relative per/female fecundity to use in model



phi.w.a <- 0.5 # Winter survival of adults (Yeiser et al.2021)
phi.b.a <- 0.5 # Breeding survival of adults (Yeiser et al. 2021), longer breeding season than winter season means monthly survival is higher in the summer (no hunting)
phi.b.j <- 0.5 # Assuming similar breeding survival between juveniles and adults
phi.chick <- 0.41 # Three month chick survival (Sandercock et al. 2008)
phi.w.j <- round((phi.w.a^(1/5))^2 * phi.chick,2) # Juvenile survival over the first winter is product of chick survival (3 months) and 2 months of adult winter survival




## Detection Parameters
##################################################################################

ptelem <- 1 # Probability of relocating a radiotagged bird should be near 1 for the entirety of sampling
tagfail <- 0.1 # Rate of tag failure
taglife <- 2 # Tags are likely to only last for <1 year (2 periods), so assuming can only be detected in deployment period and next period
pcap <- 0.34 # Band-only recapture probability (Terhune et al. 2007)
############################################################################




########################################################################




###########################################################################
#########################################################################
### Part 2. Simulating data and analyzing
#########################################################################
############################################################################

## Simulation parameters
###############################################################################
Nyears <- 2
Nperiods <- Nyears * 2 # Spring and fall sampling periods in each year
periods <- rep(1:0,times=Nyears) # 1 represents fall, 0 represents spring
Ntags <- 40 # Number of transmitters per year
tagsplit <- 0.5 # Proportion of total tags deployed in first season (fall)
Nsim <- 50 # Number of simulations


# Empty matrices to fill with simulation data
pcap.est <- phi.w.j.est <- phi.est <- tagfail.est <-relFec.est <- matrix(NA,nrow=Nsim,ncol=3)
N.act <- N.a.est <- N.j.est <- N.t.est <- vector(length=Nsim, mode="list")
###################################################################












#######################################################################
# Simulations
#####################################################################


# Population projection
##############################################################################


for (x in 1:Nsim){
  
  # Starting in fall
  
  # Separate projection matrices for winter and breeding
  proj.w <- matrix(c(0,0,phi.w.j,phi.w.a), ncol=2, byrow=T) # No breeding in winter
  proj.b <- matrix(c(relFec,relFec,phi.b.j,phi.b.a), ncol=2, byrow=T)
  proj <- list(proj.w, proj.b)
  
  N <- surviv <- chicks <- matrix(NA, ncol=2, nrow=Nperiods)
  
  
  N[1,] <- c(N.init.j, N.init.a)
  
  for (t in 2:Nperiods){
    
    # Survival
    surviv[t,1] <- rbinom(1,N[t-1,1],proj[[periods[t]+1]][2,1])
    surviv[t,2] <- rbinom(1,N[t-1,2],proj[[periods[t]+1]][2,2])
    
    # Reproduction of surviving birds, reproduction occurs at beginning of season
    chicks[t,1] <- rpois(1,N[t-1,1]*proj[[periods[t]+1]][1,1])
    chicks[t,2] <- rpois(1,N[t-1,2]*proj[[periods[t]+1]][1,2])
    
    # Abundance calculation differs between spring and fall
    # Juveniles/Adults survivors are abundance in spring
    # Juveniles are newly-hatch chicks and adults are surviving adults and juveniles in fall
    N[t,] <- (1-periods[t]) * c(surviv[t,1],surviv[t,2]) + 
      (periods[t]) * c(sum(chicks[t,]),sum(surviv[t,])) 
  }
  
  N_data <- as.data.frame(N)
  colnames(N_data) <- c("N_juv","N_adult")
  N_data$N_tot <- N_data$N_juv + N_data$N_adult
  
  N_data <- cbind(data.frame(Year=rep(1:Nyears,each=2), Period=ifelse(periods==1,"Fall","Spring")),N_data)
  
  
  N.act[[x]] <- N_data
  
  
  ## Simulating capture/recapture data
  ####################################################################################
  
  # Age and fate matrix, NA = not banded, 1 = marked/alive/juvenile, 2 = marked/alive/adult, 3 = dead 
  proj.matrix.fs <- matrix(c(phi.w.j,0,1-phi.w.j,0,phi.w.a,1-phi.w.a,0,0,1), nrow=3, ncol=3, byrow=T)
  proj.matrix.sf <- matrix(c(0,phi.b.j,1-phi.b.j,0,phi.b.a,1-phi.b.a,0,0,1), nrow=3, ncol=3, byrow=T)
  proj.matrix <- list(proj.matrix.fs, proj.matrix.sf)
  
  # Detection is either captured juvenile (1), captured adult (2), or not captured (3), will never be detected if dead
  # Only need one detection matrix, since based on current state
  cmr.det.matrix <- matrix(c(pcap,0,1-pcap,0,pcap,1-pcap,0,0,1), nrow=3, ncol=3, byrow=T)
  
  
  # Initial number banded
  N.band.init.j <- rbinom(1,N[1,1],pcap)
  N.band.init.a <- rbinom(1,N[1,2],pcap)
  
  # Initializing age-state matrix for ecological state and detection state
  cmr.state <- cmr.det <- matrix(NA, nrow=sum(N.band.init.j,N.band.init.a), ncol=Nperiods)
  cmr.state[,1] <- cmr.det[,1] <- c(rep(1, times=N.band.init.j), rep(2, times=N.band.init.a))
  
  
  
  for (t in 2:Nperiods){
    
    for (i in 1:nrow(cmr.state)){
      
      # Ecological state of marked birds
      cmr.state[i,t] <- which(rmultinom(1,1,proj.matrix[[periods[t]+1]][cmr.state[i,t-1],])==1)
      
      # Detection state of marked birds, conditional on ecological state and cap probability
      cmr.det[i,t] <- which(rmultinom(1,1,cmr.det.matrix[cmr.state[i,t],])==1)
    } # i
    
    # Need to determine number of recaptures
    recaps.j <- length(cmr.det[,t][cmr.det[,t]==1])
    recaps.a <- length(cmr.det[,t][cmr.det[,t]==2])
    
    # How many new bands in each age class, if <0, setting to 0
    new.mark.j <- max(rbinom(1,N[t,1],pcap) - recaps.j, 0)
    new.mark.a <- max(rbinom(1,N[t,2],pcap) - recaps.a, 0)
    
    # Adding newly banded birds
    cmr.state.new <- cmr.det.new <- matrix(NA, nrow=sum(new.mark.j,new.mark.a), ncol=Nperiods)
    cmr.state.new[,t] <- cmr.det.new[,t] <- c(rep(1, times=new.mark.j), rep(2, times=new.mark.a))
    
    # Combining
    cmr.state <- rbind(cmr.state, cmr.state.new)
    cmr.det <- rbind(cmr.det, cmr.det.new)
  } # t
  
  
  
  # Total number captured
  tot.cap <- nrow(cmr.det)
  
  # First period and age of capture for each individual
  age.first.cap <- period.first.cap <- period.last.cap <- rep(NA,times=tot.cap)
  for (i in 1:tot.cap){
    period.first.cap[i] <- min(which(!is.na(cmr.det[i,])))
    period.last.cap[i] <- max(which(cmr.det[i,]<3))
    age.first.cap[i] <- cmr.state[i,period.first.cap[i]]
  } #i
  
  # Number of juvenile,adult captures in each period
  j.caps <- a.caps <- rep(NA,times=Nperiods)
  for(t in 1:Nperiods){
    j.caps[t] <- length(which(cmr.det[,t]==1))
    a.caps[t] <- length(which(cmr.det[,t]==2))
  }
  
  
  # Age determines parameters such as first-winter survival, but can
  #  be calculated deterministically from age and period of first capture
  # Will be computationally faster to model captures as Bernoulli rather
  #  than age-structured categorical
  cmr.age.state <- cmr.state
  cmr.age.state[cmr.age.state==3] <- 0
  cmr.state <- cmr.age.state
  cmr.state[cmr.state>0] <- 1
  cmr.det[cmr.det==3] <- 0
  cmr.det[cmr.det>0] <- 1
  
  
  
  
  
  ## Simulating telemetry data
  #############################################################################################
  
  # Simulating number of total tags deployed in the fall, then number deployed in the spring
  #   is the remainder
  
  tags.deploy <- rep(0,times=Nperiods)
  
  tags.deploy[1] <- rbinom(1,Ntags,tagsplit)
  
  for (i in 2:Nperiods){
    tags.deploy[i] <- periods[i] * rbinom(1,Ntags,tagsplit) +
      (1-periods[i]) * (Ntags - tags.deploy[i-1])
  }
  tags.cum <- cumsum(tags.deploy)
  
  # Age state matrix is the same as for CMR, but detection matrix is different
  # Can be detected alive as juvenile (1), detected alive as adult (2), or not detected
  # Detection is 1 as long as transmitter has not failed
  tel.det.matrix <- matrix(c(1-tagfail,0,tagfail,0,1-tagfail,tagfail,0,0,1), nrow=3, ncol=3, byrow=T)
  
  
  # Empty matrices to fill out 
  tel.state <- tel.det <- matrix(NA, nrow=Ntags*Nyears, ncol=Nperiods)
  
  # Initializing in each period based on likelihood of capturing juvenile vs adult birds
  tags.juv <- rbinom(1,tags.deploy[1],N[1,1]/sum(N[1,]))
  tel.state[1:tags.cum[1],1] <- tel.det[1:tags.cum[1],1] <- c(rep(1,times=tags.juv), rep(2,times=tags.deploy[1]-tags.juv))
  for (k in 2:Nperiods){
    tags.juv <- rbinom(1,tags.deploy[k],N[k,1]/sum(N[k,]))
    tel.state[(tags.cum[k-1]+1):tags.cum[k],k] <- tel.det[(tags.cum[k-1]+1):tags.cum[k],k] <- c(rep(1,times=tags.juv), rep(2,times=tags.deploy[k]-tags.juv))
  }
  
  
  
  for (t in 2:Nperiods){
    
    for (i in 1:nrow(tel.state)){
      
      if (!is.na(tel.state[i,t-1])){
        
        # Ecological state of marked birds
        tel.state[i,t] <- which(rmultinom(1,1,proj.matrix[[periods[t]+1]][tel.state[i,t-1],])==1)
        
        # Detection state of marked birds, conditional on ecological state and probability of radio failure
        tel.det[i,t] <- which(rmultinom(1,1,tel.det.matrix[tel.state[i,t],])==1)
      } # if (!is.na(tel.state[i,t]))
    } # i
  } # t
  
  
  # Removing records for transmitters deployed in last period, not telling anything
  tel.state <- tel.state[!is.na(tel.state[,Nperiods-1]),]
  tel.det <- tel.det[!is.na(tel.det[,Nperiods-1]),]
  
  
  
  # Age of capture and deployment period for each tracked individual
  age.first.tel <- period.first.tel <- rep(NA,times=nrow(tel.state))
  for (i in 1:nrow(tel.state)){
    period.first.tel[i] <- min(which(!is.na(tel.state[i,])))
    age.first.tel[i] <- tel.state[i,period.first.tel[i]]
  } #i
  
  
  
  
  # Transmitters only last a certain time, after which birds are not detectable
  # Changing to NA when captured birds are no longer detectable
  # With a taglife of 2, will be able to track birds in deployment period and next one, but not after
  # Making a matrix of detectability
  tagactive <- matrix(NA, ncol=ncol(tel.state), nrow=nrow(tel.state))
  for (i in 1:length(period.first.tel)){
    tagactive[i,(period.first.tel[i])] <- 1
    tagactive[i,(period.first.tel[i]+1)] <- 1
  } # i
  
  # Multiplying with detection matrix to mask by possible detections
  tel.det <- tel.det * tagactive
  
  # Again, changing to binary state and detection data
  tel.det[tel.det==3] <- 0
  tel.det[tel.det>0] <- 1
  tel.state.age <- tel.state
  tel.state.age[tel.state.age==3] <- 0 
  tel.state <- tel.state.age
  tel.state[tel.state>0] <- 1
  
  # Number of detections per period
  perdet <- colSums(tel.det, na.rm=T)
  
  
  
  
  
  
  ## Simulating reproduction data
  #############################################
  
  # Only have reproduction data from radiotagged
  #  birds and only for breeding season
  # Cutting out the last year, since ending in the
  #  spring before breeding
  
  # Number of available breeders
  avail.per <- ((1-periods)*perdet)
  avail.per <- avail.per[which(avail.per>0)]
  avail <- avail.per[-length(avail.per)]

  
  # Empty matrix to store data
  reproduction.data <- rep(NA,times=0)
  
  for (p in 1:(Nyears-1)){
    reproduction.data <- c(reproduction.data,rpois(avail[p],relFec)) 
  }
  # relFec not indexed by year, so just making a vector
  
###########################################################################
  
  
  
  
  
  
  ## IPM
  #########################################################################
  
  
  
  # Model
  #######################################################################
  
  sink("NOBO_IPM_sim.jag")
  cat("
model{

  # Priors
  phi.w.j ~ dunif(0,1)
  phi ~ dunif(0,1)
  pcap ~ dunif(0,1)
  tagfail ~ dunif(0,1)
  #tagfail ~ dunif(0.08,0.12) # Informative prior, likely known with some certainty
  relFec ~ dunif(0,10)
  
  
  # Survival matrix, columns are juvenile/adult, rows are overwinter/overbreeding
  phi.mat[1,1:2] <- c(phi.w.j, phi)
  phi.mat[2,1:2] <- c(phi, phi)
  

  # CMR Likelihood
  ####################################################################
  for (i in 1:tot.cap){
  
    # Define latent state and age at first capture, age and state known for captured birds
    z.cmr[i,period.first.cap[i]] <- cmr.det[i,period.first.cap[i]] # define state at first capture
    age.cmr[i,period.first.cap[i]] <- age.first.cap[i]
    
    # State at subsequent captures
    for (t in (period.first.cap[i]+1):Nperiods){
      
      # Age. will transition from 1 to 2 if age in t-1 is 1 and period is 1 (fall)
      age.cmr[i,t] <- ifelse(age.cmr[i,t-1]==1,
                         age.cmr[i,t-1] + periods[t],
                         2)
      
      # State Process
      z.cmr[i,t] ~ dbern(z.cmr[i,t-1] * phi.mat[(periods[t]+1),age.cmr[i,t-1]])
      
      # Observation Process
      cmr.det[i,t] ~ dbern(z.cmr[i,t] * pcap)
    }
  }
  
  
  # Telemetry Likelihood
  #######################################################################
  for (i in 1:tag.count){
  
    # Define latent state and age at first capture for telemetry
    z.tel[i,period.first.tel[i]] <- tel.det[i,period.first.tel[i]] # define state at first capture
    age.tel[i,period.first.tel[i]] <- age.first.tel[i]
    
    # State in subsequent years
    for (t in (period.first.tel[i]+1):(period.first.tel[i]+1)){
      
      # Age. will transition from 1 to 2 if age in t-1 is 1 and period is 1 (fall)
      age.tel[i,t] <- ifelse(age.tel[i,t-1]==1,
                         age.tel[i,t-1] + periods[t],
                         2)
      
      # State Process
      z.tel[i,t] ~ dbern(z.tel[i,t-1] * phi.mat[(periods[t]+1),age.tel[i,t-1]])
      
      # Observation Process
      tel.det[i,t] ~ dbern(z.tel[i,t] * (1-tagfail))
    }
  }
  
  
  # Reproduction Likelihood
  ###############################################

    for (i in 1:Nnests){
        repro[i] ~ dpois(relFec)
    }
      
  
  # Abundance Estimation
  ####################################################################
  
  # Have to set initial values for N
  N.j[1] <- N.init.j
  N.a[1] <- N.init.a
  N.t[1] <- N.j[1] + N.a[1]
  
  
  for (t in 2:Nperiods){
    
    # Surviving birds
    S.j[t] ~ dbinom(phi.mat[periods[t]+1,1],N.j[t-1])
    S.a[t] ~ dbinom(phi, N.a[t-1]) # Phi is similar throughout year for adults
    
    # Recruiting juveniles, will be zero if spring (period=0), otherwise based on rel.fecund and total N in t-1
    R[t] ~ dpois(periods[t] * relFec * (N.j[t-1] + N.a[t-1]))
    
    # Abundance of juveniles is R[t] in fall and S.j[t] in spring
    N.j[t] <- R[t] + (S.j[t] * (1 - periods[t])) # first part of line will be 0 in spring, second part will be 0 in fall
    # Abundance of adults is surving adults in spring, sum of surviving adults and juveniles in fall
    N.a[t] <- S.a[t] + (S.j[t] * periods[t])
    
    N.t[t] <- N.j[t] + N.a[t]
  }
  
  # Estimating recaptures based on age-structured abundance
  for (k in 1:Nperiods){
    j.caps[k] ~ dbinom(pcap, N.j[k])
    a.caps[k] ~ dbinom(pcap, N.a[k])
  }
  
  
  
    }", fill = TRUE)
  
  sink()
  
  
  
  # JAGS data
  #####################################################################
  jd <- list(periods=periods, # Vector of period ID, 1=fall, 0=winter
             Nperiods=Nperiods, # Number of periods
             cmr.det=cmr.det, # CMR detection data
             tot.cap=tot.cap, # Total individuals captured
             period.first.cap=period.first.cap, # Period of first capture for each individual
             age.first.cap=age.first.cap, # Age at first capture for each individual
             tel.det=tel.det, # Telemetry detection data
             tag.count=nrow(tel.det), # Number of total transmitters deployed
             period.first.tel=period.first.tel, # Period of first telemetry for each bird
             age.first.tel=age.first.tel, # Age at first telemetry for each bird
             j.caps=j.caps, # Number of juvenile captures each period
             a.caps=a.caps, # Number of adult captures each period
             N.init.a=N.init.a, # Initial number of adults
             N.init.j=N.init.j, # Initial number of juveniles
             repro=reproduction.data, # Reproductive data
             Nnests=length(reproduction.data)) # Number telemetried per period
  
  
  
  # Initial values
  ########################################################################
  
  # Need to set initial values of state for CMR for periods after period.first.cap
  # Setting to alive between first and last capture
  z.cmr.init <- cmr.det
  for (i in 1:nrow(z.cmr.init)){
    for (t in 1:ncol(z.cmr.init)){
      if (t>period.first.cap[i] & t<=period.last.cap[i]){
        z.cmr.init[i,t] <- 1
      }
    }
    z.cmr.init[i,period.first.cap[i]] <- NA
  }
  
  # Need to set initial values of state for periods after period.first.tel
  # Setting to alive between first and last detection
  # Simple when taglife is 2
  z.tel.init <- tel.det
  for (i in 1:nrow(z.tel.init)){
    z.tel.init[i,period.first.tel[i]] <- NA
  }
  
  # Initial values for S.a, S.j, and R 
  # Need to ensure that they are higher than the captures in each period,
  #   but not higher than the abundance in the previous period
  # Initializing with a population projection
  N.init <- matrix(0,ncol=2,nrow=Nperiods)
  R.init <-S.a.init <- S.j.init <- rep(NA,times=Nperiods)
  N.init[1,1:2] <- c(N.init.j,N.init.a) 
  
  for (i in 2:Nperiods){
    if (periods[i]==0){
      R.init[i] <- 0
      S.j.init[i] <- round(runif(1,j.caps[i],N.init[i-1,1]))
      S.a.init[i] <- round(runif(1,a.caps[i],N.init[i-1,2]))
      N.init[i,] <- c(S.j.init[i],S.a.init[i])
    } else {
      R.init[i] <- round(runif(1,j.caps[i],sum(N[i-1,])*2))
      N.init[i,1] <- R.init[i]
      while(N.init[i,2] < a.caps[i]){
        S.j.init[i] <- round(runif(1,0,N.init[i-1,1]))
        S.a.init[i] <- round(runif(1,0,N.init[i-1,2]))
        N.init[i,2] <- S.j.init[i] + S.a.init[i]
      }
    }
  }
  
  
  ji <- function(){
    list(z.cmr=z.cmr.init,
         z.tel=z.tel.init,
         R=R.init,
         S.j=S.j.init,
         S.a=S.a.init,
         phi.w.j=runif(1,0,1),
         phi=runif(1,0,1),
         pcap=runif(1,0,1),
         tagfail=runif(1,0,1),
         relFec=runif(1,0,10))
  }
  
  
  # Params
  ######################################################################
  jp <- c("pcap","phi.w.j","phi","N.t","N.a","N.j","tagfail","relFec")
  
  
  
  # Model
  ########################################################################
  jm <- jags.model(file="NOBO_IPM_sim.jag",
                   data=jd, inits=ji,
                   n.adapt=1000, n.chains=1)
  jc.sz <- coda.samples(jm, jp, n.iter=10000)
  
  
  
  # Post-processing
  mc <- as.matrix(jc.sz)
  
  
  phi.w.j.est[x,1:3] <- quantile(mc[,colnames(mc)=="phi.w.j"],c(0.025,0.5,0.975))
  phi.est[x,1:3] <- quantile(mc[,colnames(mc)=="phi"],c(0.025,0.5,0.975))
  pcap.est[x,1:3] <- quantile(mc[,colnames(mc)=="pcap"],c(0.025,0.5,0.975))
  tagfail.est[x,1:3] <- quantile(mc[,colnames(mc)=="tagfail"],c(0.025,0.5,0.975))
  relFec.est[x,1:3] <- quantile(mc[,colnames(mc)=="relFec"],c(0.025,0.5,0.975))
  N.a.est[[x]] <- t(apply(mc[,grep("N.a",colnames(mc))],2,quantile,probs=c(0.025,0.5,0.975)))
  N.j.est[[x]] <- t(apply(mc[,grep("N.j",colnames(mc))],2,quantile,probs=c(0.025,0.5,0.975)))
  N.t.est[[x]] <- t(apply(mc[,grep("N.t",colnames(mc))],2,quantile,probs=c(0.025,0.5,0.975)))
  
}


data.list <- list(Nyears=Nyears,
                  Ntags=Ntags,
                  tagsplit=tagsplit,
                  N.init.a=N.init.a,
                  N.init.j=N.init.j,
                  N.act=N.act,
                  phi.w.j=phi.w.j,
                  phi=phi.b.a,
                  pcap=pcap,
                  tagfail=tagfail,
                  relFec=relFec,
                  phi.w.j.est=phi.w.j.est,
                  phi.est=phi.est,
                  pcap.est=pcap.est,
                  tagfail.est=tagfail.est,
                  relFec.est=relFec.est,
                  N.a.est=N.a.est,
                  N.j.est=N.j.est,
                  N.t.est=N.t.est)

filename <- paste("NOBO_IPM_simresults_",Nyears,"Years_",Ntags,"Tags.gzip",sep="")
save(data.list,file=filename)
