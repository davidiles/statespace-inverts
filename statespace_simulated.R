# Required libraries
library(ggplot2)
library(cowplot)
library(jagsUI)
library(reshape2)

set.seed(354)

rm(list=ls())

#-------------------------------------------------------------
# GENERATE DATA
#-------------------------------------------------------------

time = 120
sites = 20

logN0_GRANDMEAN = 6
logN0_SPATIALSD = 0.8

mean_emergence_GRANDMEAN = 70
mean_emergence_SPATIALSD = 10

sd_emergence_GRANDMEAN = 4

survival = 0.95

summary_of_sites = data.frame()
for (s in 1:sites){

    #Generate a value of N0, mean emergence, sd emergence, and survival for this site

    #------------------------------------------------
    #NOTE 1: Site-abundances are drawn from a lognormal distribution
    #------------------------------------------------

    N0 = rlnorm(1,meanlog = logN0_GRANDMEAN, sdlog = logN0_SPATIALSD)
    mean_emergence = rnorm(1,mean_emergence_GRANDMEAN, mean_emergence_SPATIALSD)
    sd_emergence = sd_emergence_GRANDMEAN

    summary_of_sites = rbind(summary_of_sites,
                             data.frame( site = s,
                                         N0,
                                         mean_emergence,
                                         sd_emergence,
                                         survival))
}

# True parameter values for each site
summary_of_sites
#write.csv(summary_of_sites, file = "summary_of_sites_simulated.csv",row.names=FALSE)

#Construct a matrix to store abundance dynamics for all sites (columns are sites, rows are days of season)
N = matrix(NA,nrow=time, ncol = sites)

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

    #survival
    surv = sdat$survival

    Nt = rep(0,time)
    for (t in 2:time){

        #------------------------------------------------
        #NOTE 3: Survival includes demographic stochasticity
        #------------------------------------------------

        Nt[t] = rbinom(1,Nt[t-1],surv) + recruits[t]
    }
    N[,s] = Nt

}

sim_df = as.data.frame(N)
colnames(sim_df)[1:sites] = paste("Site",1:sites)
sim_df$Day = seq(1,nrow(N))

mdata <- melt(sim_df, id=c("Day"), value.name = "true_count")
colnames(mdata) = c("Day","Site","true_count")

###########################
# Observation process (currently poisson)
# Poisson error
mdata$obs_count = rpois(nrow(mdata),mdata$true_count)
#Days each site was actually visited the field
visit_days = round(seq(1,time,length.out = 8))
missed_days = which(seq(1:time) %in% visit_days == FALSE)
mdata$obs_count[which(mdata$Day %in% missed_days == TRUE) ] = NA
###########################

#Plot simulated data (and observations)
plot = ggplot(data = mdata) +

    geom_line(aes(x = Day, y = true_count), col = "dodgerblue", size = 1)+
    geom_point(aes(x = Day, y = obs_count))+
    #geom_line(data = na.omit(mdata), aes(x = Day, y = obs_count), col = "black")+

    facet_grid(Site~., scales = "free")+
    theme_bw()+
    ggtitle("Simulated Inverts")

print(plot)

write.csv(mdata, file = "simulated_invert_data.csv", row.names=FALSE)

#-------------------------------------------------------------
# BAYESIAN ANALYSIS
#-------------------------------------------------------------

sink("invert_statespace_simulated.bug")
cat("
    model {

    #--------------------------------------
    # Priors and constraints
    #--------------------------------------

    # Population abundances
    logN0_grandmean ~ dnorm(0,0.001)
    logN0_sd ~ dunif(0,5)
    logN0_tau <- pow(logN0_sd,-2)

    # Emergence
    emergence_mean_grandmean ~ dunif(1,T) # grand mean emergence date across sites
    emergence_mean_sigma ~ dunif(0,T) # variation in mean emergence date across sites
    emergence_mean_tau <- pow(emergence_mean_sigma,-2)

    emergence_sd ~ dunif(0,T) # SD within a site (assumed to be the same for all sites)

    # Survival
    surv_grandmean ~ dunif(0.01,0.99) # survival is assumed to be same across sites and through time

    # Detection prob
    #p ~ dunif(0.01,0.99) # only required if observation model is binomial

    #--------------------------------------
    # Likelihood
    #--------------------------------------

    for (i in 1:S){

        # Fix initial population sizes to 0
        N_true[1,i] <- 0
        t_stand[1,i] <- 0
        pdf_emerged[1,i] <- 0
        cdf_emerged[1,i] <- 0
        recruits[1,i] <- 0

        N_total[i] ~ dlnorm(logN0_grandmean,logN0_tau)

        # Emergence
        mean_emerge[i] ~ dnorm(emergence_mean_grandmean, emergence_mean_tau)

        sd_emerge[i]   <- emergence_sd #all sites have same sd

        surv[i] <- surv_grandmean

        # Observation Process
        for (t in 1:T){
            #y[t,i] ~ dnorm(N_true[t,i], tau.obs[i])
            y[t,i] ~ dpois(N_true[t,i])
            #y[t,i] ~ dbin(p,N_true[t,i]) #If binomial detection

        }

        # Population Dynamics
        for (t in 2:T){

            # Emergence Dynamics
            t_stand[t,i] <- (t - mean_emerge[i])/sd_emerge[i]
            cdf_emerged[t,i] <- phi(t_stand[t,i])
            pdf_emerged[t,i] <- cdf_emerged[t,i] - cdf_emerged[t-1,i] # normal PDF for daily survival

            #binomial recruitment
            recruits[t,i] ~ dbin(pdf_emerged[t,i], round(N_total[i])) #if recruitment is binomial
            #recruits[t,i] <- round(pdf_emerged[t,i] * N_total[i]) #non-binomial recruitment

            # binomial survival
            survivors[t,i] ~ dbin(surv[i],N_true[t-1,i]) # binomial survival
            #survivors[t,i] <- round(surv[i]*N_true[t-1,i]) #non-binomial (smoothed) survival

            N_true[t,i] <-  survivors[t,i] + recruits[t,i]
        }

    }


    }
    ", fill = TRUE)
sink()

y = dcast(data = mdata,formula = Day~Site, fun.aggregate = sum, value.var = "obs_count")
y = y[,-1]

# Compile data
bugs.data <- list(y = as.matrix(y), T = nrow(y), S = ncol(y))

# Initial values
inits <- list( list(logN0_grandmean = 6,
                    logN0_sd = 0.6,
                    emergence_mean_grandmean = 50,
                    emergence_mean_sigma = 15,
                    emergence_sd = 10,

                    surv_grandmean = 0.8),

               list(logN0_grandmean = 7,
                    logN0_sd = 0.8,
                    emergence_mean_grandmean = 70,
                    emergence_mean_sigma = 12,
                    emergence_sd = 15,

                    surv_grandmean = 0.9)

)

# Parameters to be monitored
parameters <- c(

    # Hyper-parameters
    "logN0_grandmean","logN0_sd",
    "emergence_mean_grandmean","emergence_mean_sigma", "emergence_sd",
    "surv_grandmean",

    # Site-specific values
    "N_total","mean_emerge","sd_emerge","N_true")

# MCMC settings
ni <- 10000  # Number of iterations
nt <- 10      # Thinning rate
nb <- 8000  # Burn-in period
nc <- 2       # Number of chains

out <- jags(data = bugs.data,inits,parameters,"invert_statespace_simulated.bug",
            nc,nt,ni,nb)

#Save output (this file is huge)
save(out, file = "invert_simulated_basic.RData")
#rm(out)

#-------------------------------------------------------------
# SUMMARIZE AND PLOT RESULTS
#-------------------------------------------------------------
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

# Plot model estimates, compared to truth (blue)
plot = ggplot(data = mdata) +

    #Estimated abundance
    geom_ribbon(aes(x = N_median_melt$Day,
                    ymin = N_0.025_melt$est_count,
                    ymax = N_0.975_melt$est_count), alpha = 0.3)+
    # True values
    geom_line(aes(x = Day, y = true_count, col = "dodgerblue"), size = 1)+

    geom_line(aes(x = N_median_melt$Day,
                  y = N_median_melt$est_count,
                  col = "black"),size = 0.5)+

    # Observed Data
    geom_point(aes(x = Day, y = obs_count,
                   shape = "19"))+


    facet_grid(Site~., scales = "free")+
    theme_bw()+
    xlab("Day")+
    ylab("Count")+
    ggtitle("Simulated Inverts")+

    #legend
    scale_shape_manual(name = '',values =c(19), labels = c('Observed Count'))+
    scale_color_manual(name = '',values =c("black","dodgerblue"), labels = c('Estimated Abundance','True Abundance'))


print(plot)

#Save plot in working directory
ggsave("simulated_model_fit.pdf", width = 5, height = 20, units = "in")

