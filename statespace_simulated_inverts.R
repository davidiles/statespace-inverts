# Required libraries
library(ggplot2)
library(cowplot)
library(jagsUI)
library(reshape2)

set.seed(1)

rm(list=ls())

#-------------------------------------------------------------
# GENERATE DATA
#-------------------------------------------------------------

time = 90
sites = 20

logN0_GRANDMEAN = 7
logN0_SPATIALSD = 0.7

mean_emergence_GRANDMEAN = 65
mean_emergence_SPATIALSD = 10

logsd_emergence_GRANDMEAN = log(6)
logsd_emergence_SPATIALSD = 0.35

survival_grandmean = 0.95
survival_sd_logit = 0.6
p = 0.1 # Probability of capturing

summary_of_sites = data.frame()

for (s in 1:sites){

    #Generate a value of N0, mean emergence, sd emergence, and survival for this site

    #------------------------------------------------
    #NOTE 1: Site-abundances are drawn from a lognormal distribution
    #------------------------------------------------

    N0 = rlnorm(1,meanlog = logN0_GRANDMEAN, sdlog = logN0_SPATIALSD)
    mean_emergence = rnorm(1,mean_emergence_GRANDMEAN, mean_emergence_SPATIALSD)
    #sd_emergence = sd_emergence_GRANDMEAN
    sd_emergence = rlnorm(1,meanlog = logsd_emergence_GRANDMEAN, sdlog = logsd_emergence_SPATIALSD)

    survival = plogis(qlogis(survival_grandmean) + rnorm(1,0,survival_sd_logit))

    summary_of_sites = rbind(summary_of_sites,
                             data.frame( site = s,
                                         N0,
                                         mean_emergence,
                                         sd_emergence,
                                         survival))
}

# True parameter values for each site
summary_of_sites

#Construct a matrix to store abundance dynamics for all sites (columns are sites, rows are days of season)
N = matrix(NA,nrow=time, ncol = sites)
Nobs = matrix(NA,nrow=time, ncol = sites) #Observed dynamics (removal sampling)

#Days each site was actually visited
visit_days = round(seq(2,time,length.out = 8))

#Generate abundances
for (s in 1:sites){

    sdat = summary_of_sites[s,]

    #Emergence dynamics
    t = 1:time
    t_stand = (t-sdat$mean_emergence)/sdat$sd_emergence
    pdf = pnorm(t_stand)

    recruits = rep(0,time)
    for (t in 2:time){

        #------------------------------------------------
        #NOTE 2: Emergence includes demographic stochasticity
        #------------------------------------------------

        recruits[t] = rbinom(1,round(sdat$N0), pdf[t]-pdf[t-1])

    }

    #Survival
    surv = sdat$survival

    Nt = rep(0,time)
    for (t in 2:time){

        #------------------------------------------------
        #NOTE 3: Survival includes demographic stochasticity
        #------------------------------------------------

        Nt[t] = rbinom(1,Nt[t-1],surv) + recruits[t]

        if (t %in% visit_days){
            observed = rbinom(1,Nt[t],p) #Removal sampling
            Nobs[t,s] = observed
            Nt[t] = Nt[t] - observed
        }
    }
    N[,s] = Nt

}


N_df = as.data.frame(N)
colnames(N_df)[1:sites] = paste("Site",1:sites)
N_df$Day = seq(1,nrow(N))

Nobs_df = as.data.frame(Nobs)
colnames(Nobs_df)[1:sites] = paste("Site",1:sites)
Nobs_df$Day = seq(1,nrow(Nobs))

Ndata <- melt(N_df, id=c("Day"), value.name = "true_count")
colnames(Ndata) = c("Day","Site","true_count")

Nobsdata <- melt(Nobs_df, id=c("Day"), value.name = "obs_count")
colnames(Nobsdata) = c("Day","Site","obs_count")

#Plot simulated data (and observations)
plot = ggplot(data = Ndata) +

    geom_line(data = Ndata, aes(x = Day, y = true_count*p), col = "dodgerblue", size = 1)+
    geom_point(data = Nobsdata, aes(x = Day, y = obs_count))+
    #geom_line(data = na.omit(mdata), aes(x = Day, y = obs_count), col = "black")+

    facet_grid(Site~., scales = "free")+
    theme_bw()+
    ggtitle("Simulated Inverts")+
    ylab("Inverts")

print(plot)

#write.csv(mdata, file = "simulated_Invert_data.csv", row.names=FALSE)

#-------------------------------------------------------------
# BAYESIAN ANALYSIS
#-------------------------------------------------------------

sink("statespace_simulated_invert.jags")
cat("
    model {

    #--------------------------------------
    # Priors and constraints
    #--------------------------------------

    # Population abundance (hyperparameters)
    logNtotal_spatialmean ~ dunif(0,10)        # Total population size
    logNtotal_sigma ~ dunif(0,2)               # Variance in popsize across sites (random effect for popsize)
    logNtotal_tau <- pow(logNtotal_sigma,-2)

    # Population emergence (hyperparameters)
    emergence_mean_grandmean ~ dunif(0,T) # Mean emergence grand mean across sites
    emergence_mean_sigma ~ dunif(0,T)     # Mean emergence sd across sites (random effect for mean emergence)
    emergence_mean_tau <- pow(emergence_mean_sigma,-2)

    log_emergence_sd_grandmean ~ dunif(0,2)  #SD emergence grand mean across sites
    log_emergence_sd_sigma ~ dunif(0,2)      #SD emergence sd across sites (random effect for variance in emergence)
    log_emergence_sd_tau <- pow(log_emergence_sd_sigma,-2)

    # Survival
    surv_grandmean ~ dunif(0.01,0.99) # survival mean across sites
    surv_sigma_logit ~ dunif(0,2) # survival sd across sites (random effect for mean survival)
    surv_tau_logit <- pow(surv_sigma_logit,-2)

    #--------------------------------------
    # Likelihood
    #--------------------------------------
    for (i in 1:S){

    # Population Size
    logNtotal[i] ~ dnorm(logNtotal_spatialmean,logNtotal_tau)
    Ntotal[i] <- exp(logNtotal[i])

    # Emergence
    #mean_emerge[i] <- emergence_mean_grandmean
    mean_emerge[i] ~ dnorm(emergence_mean_grandmean,emergence_mean_tau)
    log_sd_emerge[i] ~ dnorm(log_emergence_sd_grandmean,log_emergence_sd_tau)
    sd_emerge[i] <- exp(log_sd_emerge[i])

    # Survival
    surv_logit[i] ~ dnorm(log(surv_grandmean/(1-surv_grandmean)),surv_tau_logit)
    surv[i] <- 1/(1+exp(-surv_logit[i]))

    # Fix initial population sizes to 0
    N_true[1,i] <- 0
    t_stand[1,i] <- 0
    pdf_emerged[1,i] <- 0
    cdf_emerged[1,i] <- 0
    recruits[1,i] <- 0

    # Observation Process
    for (t in 1:T){

    y[t,i] ~ dpois(N_true[t,i]+0.1)  #If binomial detection

    }

    # Population Dynamics
    for (t in 2:T){

    # Emergence Dynamics
    t_stand[t,i] <- (t - mean_emerge[i])/sd_emerge[i]
    cdf_emerged[t,i] <- phi(t_stand[t,i])
    pdf_emerged[t,i] <- cdf_emerged[t,i] - cdf_emerged[t-1,i] + 0.001 # Addition to improve stability

    #binomial recruitment
    #recruits[t,i] ~ dbin(pdf_emerged[t,i], round(Ntotal[i])) #if recruitment is binomial
    recruits[t,i] <- round(pdf_emerged[t,i] * Ntotal[i]) #non-binomial recruitment

    # binomial survival
    #survivors[t,i] ~ dbin(surv[i],N_true[t-1,i]) # binomial survival
    survivors[t,i] <- round(surv[i]*N_true[t-1,i]) #non-binomial (smoothed) survival

    N_true[t,i] <-  survivors[t,i] + recruits[t,i]
    }
    }


    }
    ", fill = TRUE)
sink()

y = dcast(data = Nobsdata,formula = Day~Site, fun.aggregate = sum, value.var = "obs_count")
y = y[,-1]

# Compile data
bugs.data <- list(y = as.matrix(y),
                  T = nrow(y),
                  S = ncol(y))

# Initial values
inits <- list( list(logNtotal_spatialmean = 7,
                    logNtotal_sigma = 0.4,
                    emergence_mean_grandmean = 60,
                    emergence_mean_sigma = 10,
                    log_emergence_sd_grandmean = 2,
                    log_emergence_sd_sigma = 0.3,
                    surv_grandmean = 0.8,
                    surv_sigma_logit = 0.7),

               list(logNtotal_spatialmean = 6,
                    logNtotal_sigma = 0.6,
                    emergence_mean_grandmean = 40,
                    emergence_mean_sigma = 10,
                    log_emergence_sd_grandmean = 1,
                    log_emergence_sd_sigma = 0.5,
                    surv_grandmean = 0.9,
                    surv_sigma_logit = 0.5)

)

# Parameters to be monitored
parameters <- c("logNtotal_spatialmean","logNtotal_sigma",
                "emergence_mean_grandmean","emergence_mean_sigma",
                "log_emergence_sd_grandmean","log_emergence_sd_sigma",
                "surv_grandmean","surv_sigma_logit",

                "mean_emerge","Ntotal","surv",
                "N_true")

# MCMC settings
ni <- 50000  # Number of iterations
nt <- 10      # Thinning rate
nb <- 45000  # Burn-in period
nc <- 2       # Number of chains

out <- jags(data  = bugs.data,
            inits,
            parameters,
            "statespace_simulated_invert.jags",
            nc,
            nt,
            ni,
            nb)

#Save output (this file is huge - so this is commented out right now)
#save(out, file = "Invert_simulated_basic.RData")
#rm(out)

#-------------------------------------------------------------
# SUMMARIZE AND PLOT RESULTS
#-------------------------------------------------------------

# Calculate posterior summary statistics for estimated N at each site
N_median = N_0.025 = N_0.975 = as.data.frame(y*NA)

for (s in 1:ncol(y)){
    N_median[,s] = apply(out$sims.list$N_true[,,s],2,median)
    N_0.025[,s] = apply(out$sims.list$N_true[,,s],2,function(x) quantile(x,0.025))
    N_0.975[,s] = apply(out$sims.list$N_true[,,s],2,function(x) quantile(x,0.975))
    print(s)
}

N_median$Day = N_0.025$Day = N_0.975$Day = 1:nrow(N_median)

N_median_melt <- melt(N_median, id=c("Day"), value.name = "est_count")
N_0.025_melt <- melt(N_0.025, id=c("Day"), value.name = "est_count")
N_0.975_melt <- melt(N_0.975, id=c("Day"), value.name = "est_count")

colnames(N_median_melt) = colnames(N_0.025_melt) = colnames(N_0.975_melt) = c("Day","Site","est_count")

# Plot model estimates (black) compared to truth (blue)
plot = ggplot(data = Ndata) +

    #Estimated abundance
    geom_ribbon(aes(x = N_median_melt$Day,
                    ymin = N_0.025_melt$est_count,
                    ymax = N_0.975_melt$est_count), alpha = 0.3)+
    # True values
    geom_line(aes(x = Day, y = true_count*p, col = "dodgerblue"), size = 1)+

    geom_line(aes(x = N_median_melt$Day,
                  y = N_median_melt$est_count,
                  col = "black"),size = 0.5)+

    # Observed Data
    geom_point(data = Nobsdata, aes(x = Day, y = obs_count,
                                    shape = "19"))+


    facet_grid(Site~., scales = "free")+
    theme_bw()+
    xlab("Day")+
    ylab("Count")+
    ggtitle(paste("Simulated Inverts; ",ni," iterations", sep=""))+

    # legend
    scale_shape_manual(name = '',values =c(19), labels = c('Observed Count'))+
    scale_color_manual(name = '',values =c("black","dodgerblue"), labels = c('Estimated Abundance','True Abundance'))


print(plot)

#Save plot in working directory
ggsave("simulated_model_fit.pdf", width = 6, height = 15, units = "in")
