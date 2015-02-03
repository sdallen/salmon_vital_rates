################################################################################
################################################################################
# This program uses Model MLB-4 to estimate conditional natural mortality rates, 
# exploitation rates, maturation rates and 3 compound probabilities (optional) that maximizes
# the likelihood of the given catch and maturation at age data
# for SIMULATED datasets.
################################################################################
# For analysis of multiple release groups per brood year and/or the inclusion of 
# incomplete brood years, contact author (shanae.allen@noaa.gov). 
################################################################################
# This model formulation (MLB-4) is time/age variant for 
# natural mortality rates (age class and calendar years are additive on complementary log-log scale), 
# exploitation rates (age class and calendar years areadditive on complementary log-log scale), 
# maturation rates (ages and cohorts are additive on complementary log-log scale).
# Natural mortality age classes are ages 1 and ages 2-4 and age classes for
# exploitation rates are age 2, age 3, and age 4-5.
# This model is similar to that in Mohr, Hankin, and Speed (unpublished, 1994)
###############################################################################
# Convergence Criteria - standard error file must have been created but
# this script does NOT screen out estimated standard errors that are equal to -1 or NA
###############################################################################
# this script allows for missing brood years
###############################################################################
# Data set file must have first column specifiying i.data (in order to use this script for both simulated and single datasets)
# 

# If SIMULATED/DETERMINISTIC
# data file column descriptions
#  i.data - Dataset number
#  BY - Brood year
#  AGE - age
#  R - release size
#  i - impact generating rate (is not dependent on population size)
#  m - maturation generating rate
#  d - natural mortality generating rate 
#  i.rate - realized impact rate
#  m.rate - realized maturation rate
#  d.rate - realized natural mortality rate
#  N - Beginning of age abundance 
#  S.i - number of impacts sampled 
#  S.m - number of fish that matured sampled
#  I - number of impacts
#  M - number of maturation deaths
#  D - number of natural mortality deaths
#  I.est - estimated number of impacts (expanded for sampling)
#  M.est - estimated number of maturation deaths (expanded for sampling)



# If REAL
# data file column descriptions
#  i.data - Dataset number
#  BY - Brood year
#  AGE - age
#  R - release size
#  i.cr.est - CR estimated realized impact rate (is dependent on population size)
#  m.cr.est - CR estimated realized maturation rate
#  d.cr.est - CR estimated realized natural mortality rate 
#  N.cr.est - CR estimated beginning-of-age abundance 

#  I.est - estimated number of impacts (expanded for sampling)
#  M.est - estimated number of maturation deaths (expanded for sampling)

# ==============================================================================


rm(list = ls(all = TRUE))

# ==============================================================================
# Libraries/ Settings

library(R2admb)


# ==============================================================================
# Directories/Files.


ADMB_HOME <- "C:/ADMB"


setup_admb(ADMB_HOME)




# =============================================================================
# functions

#  source( paste( Sharedfn.dir,                   "Mod.SAT.LL.r", sep="/" ) )

CLL <- function(x){
  
  log(-log(1 - ( x ) ))
}

CLLinv <- function(x){
  
  1 - exp( -exp(x))
}

# ==============================================================================
# User defined file names and quantities

# Number of data sets to fit 
  # max number = 100 if using stochastic simulated data sets  
  #            = 1 if using deterministic data set

n.data.sets <- 1 #100

# Number of times to refit the model each time using different starting points
n.iter <- 100

# Specify model
model <- "mlb4"

#  admb file name
admb.file <- "mlb4_no_s_c4" #changed convergence criteria back to default

#  age 5 maturation rate 
m5constant <- 1  



  
  # Specify the set of generating rates 
  
 gen.rates <- "CR_tv_GenRates_IGH_Y_no_missBY" # This set of generating rate does NOT satisfy model assumptions
 # gen.rates <- "mu_tv_GenRates_IGH_Y_no_missBY" # This set of generating rates satisfies all model assumptions

  # Specify the dataset file (no extension)

  #data <- "sim.data.no.missing.cohorts" #stochastic data sets
  data <- "det.data.no.missing.cohorts" #deterministic (noise-free) data sets

  # Specify filename of dataset file 
  
  data.fn <- paste( gen.rates, paste(data, ".txt", sep=""), sep="/")
  
  # Specify output file
  mlb4.fn <-paste(model, ".output.", gen.rates, ".", data, ".txt", sep="")
  




# ==============================================================================
# Filenames

# ADMB files

admb.genfile <-   paste(admb.file, "_gen", sep="")
admb.gen.stdfile <- paste(admb.file, "_gen.std", sep="")
admb.gen.evafile <- paste(admb.file, "_gen.eva", sep="")  
admb.gen.parfile <- paste(admb.file, "_gen.par", sep="")  


# =============================================================================
# Read data file


dat_all <- read.table( data.fn, header=TRUE, sep=",", as.is=TRUE, colClasses="numeric")
colnames(dat_all)


#obtain ages
ages <- sort( unique( dat_all$AGE))
n.ages <- length( ages)


# =============================================================================
# Fixed quantities


ages_c <- ages[-1]
n.ages_c <- length(ages_c) # this cannot be changed

ages_d <- ages_c - 1
n.ages_d <- length(ages_d) # this cannot be changed


# Specify (unique) brood years and calculate calendar years

broodyrs <-  sort( unique(dat_all$BY))
n.broods <- length(broodyrs)

calyrs <- sort( unique(c(  broodyrs+2,  broodyrs +3,  broodyrs +4,  broodyrs +5)))
n.calyrs <- length(calyrs)   

calyrs_c <- sort( unique(c(  broodyrs+2,  broodyrs +3,  broodyrs +4,  broodyrs +5)))
n.calyrs_c <- length(calyrs_c)   

mincalyr_c <- min( calyrs_c)

calyrs_d <- sort( unique(c(  broodyrs +1,  broodyrs +2,  broodyrs +3,  broodyrs +4)))
n.calyrs_d <- length(calyrs_d)


# calculate number of release groups per brood year per data set

n.release <-  as.numeric( table( dat_all$BY[dat_all$AGE==2 & dat_all$i.data==1],  dat_all$AGE[dat_all$AGE==2 & dat_all$i.data==1]))
release.size <- dat_all$R[dat_all$AGE==2 & dat_all$i.data==1]


# put out a warning if n.release is greater than 1 

if ( sum( n.release > 1 ) > 0 ){
  warning( "Number of release groups exceeds 1 per brood year, Run alternate program")
}


# =============================================================================
# Model requirements

# ages or age classes to estimate age effect for exploitation rate
age.class_c <- c("2", "3", "4-5")

# number of ages or age classes to estimate age effect for exploitation rate
n.age.class_c <- length(age.class_c)

# ages or age classes to estimate age effect for natural mortality
age.class_d <- c("1", "2-4")

# number of ages or age classes to estimate age effect for natural mortality
n.age.class_d <- length(age.class_d)

# ages or age classes to estimate age effect for maturation rate
age.class_m <- c("2", "3")

# ages to estimate age effect for maturation rate 

n.age.class_m <- length(age.class_m)

# number of unconditional probabilities to estimate
n.s <- 0 

# cal years to estimate year effect for exploitation rate (depends on unconditional probabilities)
est.calyrs_c <- calyrs_c 
# first year to estimate exploitation rates
min.est.calyr_c <- min(est.calyrs_c)
# number of cal years to estimate year effect for exploitation rate (depends on unconditional probabilities)
n.est.calyrs_c <- length(est.calyrs_c)



# cohorts to estimate cohort effect for maturation rate (depends on unconditional probabilities)
est.bys <- broodyrs
# number of cohorts to estimate cohort effect for maturation rate (depends on unconditional probabilities)
n.est.bys <- length(est.bys)


# cal years to estimate year effect for age 1 natural mortality rate (depends on unconditional probabilities)
est.calyrs_d1 <- broodyrs + 1
min.est.calyr_d1 <- min(est.calyrs_d1)
# number of cal years to estimate year effect for natural mortalit rate (depends on unconditional probabilities)
n.est.calyrs_d1 <- length(est.calyrs_d1)

# cal years to estimate age 2 natural mortality rate (depends on unconditional probabilities)
est.calyrs_d2 <- broodyrs + 2
min.est.calyr_d2 <- min(est.calyrs_d2)
# number of cal years to estimate year effect for age 2 natural mortality rate  (depends on unconditional probabilities)
n.est.calyrs_d2 <- length(est.calyrs_d2)

# cal years to estimate age 3 natural mortality rate (depends on unconditional probabilities)
est.calyrs_d3 <- broodyrs + 3
min.est.calyr_d3 <- min(est.calyrs_d3)
# number of cal years to estimate year effect for age 3 natural mortality rate  (depends on unconditional probabilities)
n.est.calyrs_d3 <- length(est.calyrs_d3)

# cal years to estimate age 4 natural mortality rate (depends on unconditional probabilities)
est.calyrs_d4 <- broodyrs + 4
min.est.calyr_d4 <- min(est.calyrs_d4)
# number of cal years to estimate year effect for age 4 natural mortality rate  (depends on unconditional probabilities)
n.est.calyrs_d4 <- length(est.calyrs_d4)

# cal years to estimate year effects for age 2-4 natural mortality rate (depends on unconditional probabilities)
est.calyrs_d24 <- sort( unique( c( est.calyrs_d2, est.calyrs_d3, est.calyrs_d4)) )
min.est.calyr_d24 <- min(est.calyrs_d24)
# number of cal years to estimate year effect for age 2 natural mortality rate  (depends on unconditional probabilities)
n.est.calyrs_d24 <- length(est.calyrs_d24)

# cal years to estimate year effects for age 1-4 natural mortality rate (depends on unconditional probabilities)
est.calyrs_d <- sort( unique( c(est.calyrs_d1, est.calyrs_d2, est.calyrs_d3, est.calyrs_d4)) )
min.est.calyr_d <- min(est.calyrs_d)
# number of cal years to estimate year effect for ages 1-4 natural mortality rate  (depends on unconditional probabilities)
n.est.calyrs_d <- length(est.calyrs_d)


n.params = n.est.calyrs_d + n.age.class_d + n.est.calyrs_c + n.age.class_c + n.est.bys + n.age.class_m + 1 + n.s 


# the assigned baseline brood year - the first brood year OR the first brood year with ages 2 and 3 estimated
baseline.BY <- broodyrs[1]    

# the assigned baseline cal year - the first cal year with effects of all age classes estimated
baseline.CY_d <- est.calyrs_d[1]   


################################################################################
################################################################################
################################################################################
# START
################################################################################
################################################################################
################################################################################

k = 1

j = k


mlb4 <- NULL

for( k in j:n.data.sets){
  
  # =================================================================================
  # Vector of catches per age and release group followed by number matured per age and release group
  # used by ADMB
  # One vector must have Catch/Impacts (C2, C3, C4, ...) for each brood year followed by the same for escapement
  
  yvec <- c(dat_all$I.est[ dat_all$i.data==k & dat_all$AGE > 1], dat_all$M.est[ dat_all$i.data==k & dat_all$AGE > 1] ) *1.0 # To ensure yvec is treated as type double
  
  nobs <- length(yvec)
  
  # Vector of inital tagged releases
  nvec <-  rep(dat_all$R[ dat_all$i.data==k & dat_all$AGE > 1], 2 )*1.0 # To ensure nvec is treated as type double 
  
  # if dealing with incomplete cohorts, but this code isn't really set up for this?  
  nvec <- nvec[!is.na(yvec)]	
  
  yvec <- yvec[!is.na(yvec)]	
  nobs <- length(yvec)
  
  nobs_c <- sum( !is.na(dat_all$I.est[ dat_all$i.data==k ]))
  # =============================================================================
  # build data.frame that will store model results for all iterations and for each dataset
  
  
  mlb4_iter <- dat_all[ dat_all$i.data==k, ]
  
  # add columns to store model results
  
  mlb4_iter$i.iter <- NA
  
  mlb4_iter$d.est <- NA
  mlb4_iter$i.est <- NA
  mlb4_iter$m.est <- NA
  
  mlb4_iter$d.sd <- NA
  mlb4_iter$i.sd <- NA
  mlb4_iter$m.sd <- NA
  
  
  mlb4_iter$LL <- NA
  mlb4_iter$maxgrad <- NA
  
  
  
  # add columns to base starting points on
  
    mlb4_iter$m.start <- mlb4_iter$m.rate
    mlb4_iter$d.start <- mlb4_iter$d.rate
    mlb4_iter$i.start <- mlb4_iter$i.rate
  
  
  
  ################################################################################
  ################################################################################
  ################################################################################
  # GENERATE INITS - Select randomly each iteration
  ################################################################################
  ################################################################################
  ################################################################################
  
  
  iter <- 1
  
  iter.current <- iter
  
  for( iter in iter.current:n.iter) { 
    
    
    mlb4_iter$i.iter <- iter
    
    # =============================================================================
    # Draw initial maturation rates - starting with brood year 1 + age 3 = calendar year 4
    
    m.effect.names <- c( paste( "Age 2/", baseline.BY), paste( "Age 3/", baseline.BY), 
                         paste("BYr ", est.bys[est.bys!=baseline.BY]) )
    
    
    logbm_jk <- rep(0, length(m.effect.names ) ) 
    
    names(logbm_jk) <- m.effect.names 
    
    # intercept <- age 2, baseline brood year	= 1
    mu <- CLL( mean(mlb4_iter$m.start[ mlb4_iter$AGE==2], na.rm=TRUE) ) 
    logbm_jk[paste( "Age 2/", baseline.BY)] <- rnorm( 1, mean = mu, sd = 0.40* abs( mu) ) 
    
    # logbm_jk[paste( "Age 2/", baseline.BY)] must not be too negative (m2 > 0.00033) 
    logbm_jk[paste( "Age 2/", baseline.BY)] <- max( logbm_jk[paste( "Age 2/", baseline.BY)], -8)
    
    # Effect of age 3, baseline brood year
    mu <- CLL( mean( mlb4_iter$m.start[  mlb4_iter$AGE==3], na.rm=TRUE) ) - logbm_jk[1]
    logbm_jk[paste( "Age 3/", baseline.BY)] <- rnorm( 1, mean = mu, sd = 0.40* abs( mu) ) 
    
    # Assume no cohort effects for inits
    logbm_jk[paste("BYr ", est.bys[est.bys!=baseline.BY])] <- 0 
    
    
    # logbm_jk[paste( "Age 3/", baseline.BY)] must be greater than zero (i.e., maturation rate of age 3 is greater than maturation rate of age 2
    logbm_jk[paste( "Age 3/", baseline.BY)] <- max( 0.05, logbm_jk[paste( "Age 3/", baseline.BY)] )
    
    # logbm_jk[paste( "Age 2/", baseline.BY)] + logbm_jk[paste( "Age 3/", baseline.BY)] must be less than 2 (i.e., maturation rate of age 3 is not too close to 1)
    logbm_jk[paste( "Age 3/", baseline.BY)] <- min( 2, logbm_jk[paste( "Age 2/", baseline.BY)] + logbm_jk[paste( "Age 3/", baseline.BY)] ) - logbm_jk[paste( "Age 2/", baseline.BY)]
    
    
    # Inits for age 4 constant maturation rate
    mu <- CLL( mean(mlb4_iter$m.start[  mlb4_iter$AGE==4], na.rm=TRUE) )
    cll.m4 <- rnorm( 1, mean = mu, sd = 0.40* abs( mu) ) 
    
    # cll.m4 must be greater than logbm_jk[1] + logbm_jk[2] (i.e., maturation rate of age 4 is greater than or equal to age 3)
    cll.m4 <- max( logbm_jk[1] + logbm_jk[2], cll.m4 )  
    
    # cll.m4 must be less than 2 (i.e., maturation rate of age 4 is not too close to 1)
    cll.m4 <- min( 2, cll.m4 ) 
    
    
    
    # =============================================================================
    # Draw initial Exploitation Rates - indexed by age (j) - starting with 1nd brood year + age 2 = 3rd calendar year 
    
    c.effect.names <-  c( "Age 2", "Age 3", 
                          paste("Age 45/", est.calyrs_c) )
    
    
    logbc_jk <-  rep(0, length( c.effect.names) ) 
    names(logbc_jk) <- c.effect.names
    
    
    # age 45, all estimated exploitation rate years, assume logbc_jk for all age 45 are equal
    mu <- CLL( mean( mlb4_iter$i.start[mlb4_iter$AGE==4], na.rm=TRUE) )
    logbc_jk[ paste("Age 45/", est.calyrs_c)  ] <- rnorm( 1, mean=mu, sd = 0.40* abs( mu  ))
    
    # age 45, must not be too positive or negative (between -8 and 2)***cannot specify a higher lower limit because fisheries can be closed
    
    logbc_jk[ paste("Age 45/", est.calyrs_c)  ] <- pmin( 	pmax( logbc_jk[ paste("Age 45/", est.calyrs_c)] , -8), 2)
    
    # Effect of age 2
    mu <- CLL( mean(mlb4_iter$i.start[mlb4_iter$AGE==2], na.rm=TRUE) ) - mean(logbc_jk[paste("Age 45/", est.calyrs_c)])
    logbc_jk["Age 2"] <- rnorm( 1, mean=mu, sd = 0.40* abs( mu  ))
    
    # Effect of age 3
    mu <- CLL( mean(mlb4_iter$i.start[mlb4_iter$AGE==3], na.rm=TRUE) ) - mean(logbc_jk[paste("Age 45/", est.calyrs_c)])
    logbc_jk["Age 3"] <- rnorm( 1, mean=mu, sd = 0.40* abs( mu  ))
    
    # logbc_jk["Age 3"] must be less than or equal to 0 but logbc_jk["Age 3"] + logbc_jk["Age45"] >= -8 ( c_y,3 >=.0003 and c_y,3 <= c_y,45)
    # Thus 0 <= log(q3) <= -log(fCY) - 8 
    
    logbc_jk["Age 3"] <- max( min( 0, logbc_jk["Age 3"]), -logbc_jk[ paste("Age 45/", est.calyrs_c)]  - 8)
    
    # logbc_jk["Age 2"] must be less than logbc_jk["Age 3"] (i.e., exploitation rate of age 2 is less than exploitation rate of age 3) 
    logbc_jk["Age 2"] <- min( logbc_jk["Age 3"] - 0.05, logbc_jk["Age 2"] )
    
    # logbc_jk["Age 2"] + logbc_jk["Age 45"] must be greater than -9 ( c_y,2 >=.0001 )
    
    logbc_jk["Age 2"] <- max( logbc_jk["Age 2"], -logbc_jk[ paste("Age 45/", est.calyrs_c)] - 9 )
    
    
    # =============================================================================
    # Draw initial Natural mortality rates - each age class is treated independently 
    
    d.effect.names <- c( paste( "Age 1/", baseline.CY_d), paste( "Age 24/", baseline.CY_d), 
                         paste("CYr ", est.calyrs_d[est.calyrs_d!=baseline.CY_d]) )
    
    
    logbd_jk <- rep(0, length(d.effect.names ) ) 
    
    names(logbd_jk) <- d.effect.names 
    
    # intercept <- age 1, baseline cal year	= 1
    mu <- CLL( mean(mlb4_iter$d.start[ mlb4_iter$AGE==1], na.rm=TRUE) ) 
    logbd_jk[paste( "Age 1/", baseline.CY_d)] <- rnorm( 1, mean = mu, sd = 0.40* abs( mu) ) 
    
    # logbd_jk[paste( "Age 1/", baseline.CY_d)] must not be too positive
    logbd_jk[paste( "Age 1/", baseline.CY_d)] <- min( logbd_jk[paste( "Age 1/", baseline.CY_d)], 2)
    
    # Effect of age 24, baseline cal year
    mu <- CLL( mean( mlb4_iter$d.start[  mlb4_iter$AGE %in% c(2,3)], na.rm=TRUE) ) - logbd_jk[1]
    logbd_jk[paste( "Age 24/", baseline.CY_d)] <- rnorm( 1, mean = mu, sd = 0.40* abs( mu) ) 
    
    # Assume no cohort effects for inits
    logbd_jk[paste("CYr ", est.calyrs_d[est.calyrs_d!=baseline.CY_d])] <- 0 
    
    
    # logbd_jk[paste( "Age 24/", baseline.CY_d] must be less than zero (i.e., natural mortality rate of age 24 is less than natural mortality rate of age 1
    logbd_jk[paste( "Age 24/", baseline.CY_d)] <- min( -0.05, logbd_jk[paste( "Age 24/", baseline.CY_d)] )
    
    # logbd_jk[paste( "Age 1/", baseline.CY_d)] + logbd_jk[paste( "Age 24/", baseline.CY_d)] must be greater than  -3.676247 (i.e., natural mortality rate of age 24 is greater than .025)
    logbd_jk[paste( "Age 24/", baseline.CY_d)] <- max( -3.7, logbd_jk[ paste( "Age 1/", baseline.CY_d)] + logbd_jk[paste( "Age 24/", baseline.CY_d)] ) - logbd_jk[paste( "Age 1/", baseline.CY_d)]
    
    # =============================================================================
    # Set starting values for parameters
    
    
    logbc_jkstart <- logbc_jk
    logbd_jkstart <- logbd_jk
    logbm_jkstart <- logbm_jk
    
    cllm4start <- cll.m4
    
    
    ##################################################################################
    ##################################################################################
    # MAXIMIZE LOG - LIKELIHOOD
    ##################################################################################
    ##################################################################################
    
    # 1. List data - change variable names with '.' to '_' or remove
    
    z <- list( yvec = yvec, nvec = nvec, nobs = nobs,
               nestcalyrs_c = n.est.calyrs_c, ages_c = ages_c, nages_c = n.ages_c,  
               nestcalyrs_d = n.est.calyrs_d, ages_d = ages_d, nages_d = n.ages_d,  
               estbys = est.bys, nestbys = n.est.bys, 
               m5constant = m5constant,
               lengthbd_jk = length(logbd_jk), 
               lengthbm_jk = length(logbm_jk), 
               
               minestcalyr_d = min.est.calyr_d,
               lengthbc_jk = length(logbc_jk), minestcalyr_c = min.est.calyr_c, nparams = n.params
    )
    
    
    
    # 3. List inits
    
    
    inits <- list(
      logbm_jk = logbm_jkstart,
      logbc_jk = logbc_jkstart,
      logbd_jk = logbd_jkstart,
      
      cllm4 = cllm4start
    )
    
    params <- c( names(logbm_jk),
                 names(logbc_jk),	names(logbd_jk), "m4"  )
    
    
    
    
    ##################################################################################
    ##################################################################################
    # Perform UNBOUNDED Optimization
    ##################################################################################
    ##################################################################################
    
    clean_admb(admb.file, which="all")
    
    
    
    # This is to prevent breaking out of the loop if the hessian is not positive definite
    # bounds - cannot specify a higher lower limit on exploitation rates because fisheries can be closed
    
    if(   
      class( try(
        admb.out <- do_admb(admb.file,
                            data=z, params=inits, 
                            bounds=list(logbd_jk=c(-9,2), logbc_jk=c(-9,8),logbm_jk=c(-9,8)),  #, m_jk = c(.001, .999), clld1=c(0,2), clld24=c(-3,1), logbm_jk = c(-15, 15)), #d1=c(.63, .999), d2=c(.05, .93), d3=c(.05, .93), d4=c(.05, .93)
                            run.opts=run.control(checkparam="write",
                                                 checkdata="write", clean_files=FALSE)) ,  #, extra.args="-rs"
        silent= FALSE) ) == "admb"  &
        file.exists(  admb.gen.stdfile ) 
    )
      
      
      
    {
      
      
      # --------------------------------------------------------------------------------
      # Do it MANUALLY if do_admb crashes - can be used to debug or convergence 
      
      # compile_admb(admb.genfile, safe = TRUE, admb_errors="warn", verbose = TRUE)
      
      # run_admb(admb.genfile, verbose = TRUE)  #, extra.args="-rs"
      
      # read_pars(admb.file)
      
      # admb.out <- read_admb(admb.genfile,verbose=TRUE)
      
      # clean_admb(admb.file,which="all")
      
      # --------------------------------------------------------------------------------
      # Extract/view Results
      
      # summary(admb.out)
      
      # variance - covariance matrix: look out for high collinearity (1, -1)
      # vcov <- vcov(admb.out) 
      
      #eigvals <- read.table( admb.gen.evafile)[1,-1]
      
      # --------------------------------------------------------------------------------
      # Read .par file manually
      
      raw.parfile <-  scan( admb.gen.parfile , what="raw", quiet=TRUE)
      
      maxgrad <- as.numeric( raw.parfile[16])
      
      # --------------------------------------------------------------------------------
      # READ STD FILE to get natural and fishing mortality on the proportion scale with standard errors 
      
      vars.table <- read.table(admb.gen.stdfile, header=FALSE, skip = 1, as.is=TRUE, 
                               col.names = c("index", "name", "value", "std.dev"))
      
      c.est.vec <- vars.table[ vars.table[, "name"]=="c_jk", "value" ] 
      c.sd.vec <- vars.table[ vars.table[, "name"]=="c_jk", "std.dev" ] 
      
      c.est <- t(matrix(c.est.vec, nrow=n.ages_c, ncol=n.est.calyrs_c, dimnames = list(ages_c,est.calyrs_c)))
      c.sd <- t(matrix(c.sd.vec, nrow=n.ages_c, ncol=n.est.calyrs_c, dimnames = list(ages_c,est.calyrs_c)))
      
      #natural mortality
      
      d.est.vec <- vars.table[ vars.table[, "name"]=="d_jk", "value" ] 
      d.sd.vec <- vars.table[ vars.table[, "name"]=="d_jk", "std.dev" ] 
      
      d.est <- t(matrix(d.est.vec, nrow=n.ages_d, ncol=n.est.calyrs_d, dimnames = list(ages_d,est.calyrs_d)))
      d.sd <- t(matrix(d.sd.vec, nrow=n.ages_d, ncol=n.est.calyrs_d, dimnames = list(ages_d,est.calyrs_d)))
      
      m.est.vec <- vars.table[ vars.table[, "name"]=="m_jk", "value" ] 
      m.sd.vec <- vars.table[ vars.table[, "name"]=="m_jk", "std.dev" ] 
      
      m.est <- t(matrix(m.est.vec, nrow=n.ages_c-1, ncol=n.est.bys, dimnames = list(c("2", "3", "4"), est.bys)))
      m.sd <- t(matrix(m.sd.vec, nrow=n.ages_c-1, ncol=n.est.bys, dimnames = list(c("2", "3", "4"),   est.bys)))
      
      
      
      objfuncval <- -logLik(admb.out)
      
      
      # --------------------------------------------------------------------------------
      # Assign estimates to dataframe
      # Make sure to assign to correct brood year!! CHECK THIS
      
      
      mlb4_iter$i.est[mlb4_iter$AGE==1] <- 0
      mlb4_iter$i.est[mlb4_iter$AGE==2] <- c.est[as.character(broodyrs+2), "2"]
      mlb4_iter$i.est[mlb4_iter$AGE==3] <- c.est[as.character(broodyrs+3), "3"]
      mlb4_iter$i.est[mlb4_iter$AGE==4] <- c.est[as.character(broodyrs+4), "4"]
      mlb4_iter$i.est[mlb4_iter$AGE==5] <- c.est[as.character(broodyrs+5), "5"]
      
      
      mlb4_iter$i.sd[mlb4_iter$AGE==2] <- c.sd[as.character(broodyrs+2), "2"]
      mlb4_iter$i.sd[mlb4_iter$AGE==3] <- c.sd[as.character(broodyrs+3), "3"]
      mlb4_iter$i.sd[mlb4_iter$AGE==4] <- c.sd[as.character(broodyrs+4), "4"]
      mlb4_iter$i.sd[mlb4_iter$AGE==5] <- c.sd[as.character(broodyrs+5), "5"]
      
      
      mlb4_iter$m.est[mlb4_iter$AGE==1] <- 0
      mlb4_iter$m.est[mlb4_iter$AGE==2] <- m.est[, "2"]
      mlb4_iter$m.est[mlb4_iter$AGE==3] <- m.est[, "3"]
      mlb4_iter$m.est[mlb4_iter$AGE==4] <- m.est[, "4"]
      mlb4_iter$m.est[mlb4_iter$AGE==5] <- 1
      
      
      mlb4_iter$m.sd[mlb4_iter$AGE==2] <- m.sd[, "2"]
      mlb4_iter$m.sd[mlb4_iter$AGE==3] <- m.sd[, "3"]
      mlb4_iter$m.sd[mlb4_iter$AGE==4] <- m.sd[, "4"]
      
      
      
      mlb4_iter$d.est[mlb4_iter$AGE==1] <- d.est[as.character(broodyrs+1), "1"]
      mlb4_iter$d.est[mlb4_iter$AGE==2] <- d.est[as.character(broodyrs+2), "2"]
      mlb4_iter$d.est[mlb4_iter$AGE==3] <- d.est[as.character(broodyrs+3), "3"]
      mlb4_iter$d.est[mlb4_iter$AGE==4] <- d.est[as.character(broodyrs+4), "4"]
      
      
      mlb4_iter$d.sd[mlb4_iter$AGE==1] <- d.sd[as.character(broodyrs+1), "1"]
      mlb4_iter$d.sd[mlb4_iter$AGE==2] <- d.sd[as.character(broodyrs+2), "2"]
      mlb4_iter$d.sd[mlb4_iter$AGE==3] <- d.sd[as.character(broodyrs+3), "3"]
      mlb4_iter$d.sd[mlb4_iter$AGE==4] <- d.sd[as.character(broodyrs+4), "4"]
      
      
      
      mlb4_iter$LL <- objfuncval
      mlb4_iter$maxgrad <-  maxgrad
      
      
    }# end of if try ==TRUE
    
    
    clean_admb(admb.file, which="all")
    
    
    mlb4 <- rbind( mlb4, mlb4_iter)	 
    
    # ===================================================================================================
    # END of FOR ITERATIONS
  } 
  
  
  
  # ===================================================================================================
  # END of FOR DATASETS
} 



# ===================================================================================================
# Write to ASCII file 

write.table( mlb4, file = mlb4.fn, row.names=FALSE, col.names=TRUE, sep=", ",
             append=FALSE, quote=FALSE)

# ==================================================================================================


