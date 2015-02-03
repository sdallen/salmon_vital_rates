# This function simulates data assuming given generating rates
# Allows for sampling noise in escapement data
# Allows for missing brood years
# For 1 release group per year
# broodyrs and release.size are vectors of the same size
# Vital rates include s - survival, mu - exploitation, sigma - maturation.
# Survival (s) and exploitation (mu) are indexed by calendar year and age
# while maturation (sigma) is indexed by brood year and age. 
# ocean.sample.fraction/river.sample.fraction is a single number 
# specifying the percent of total catch/escapement that's sampled
# Estimated catch/escapement is calculated by expanding the sampled catch/escapement by
# the assumed sampling fraction (can be different than simulated)
# Assumes no river harvest and no river natural mortality


# column descriptions
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
#  I - estimated number of impacts (expanded for sampling)
#  M - estimated number of maturation deaths (expanded for sampling)

# ==============================================================================
  sim.data <- function( gen.rate.dat,  ocean.sample.fraction, 
			 river.sample.fraction, broodyrs,calyrs, ages_c,
			  numbroods){


      sim.dat <- gen.rate.dat

 	sim.dat$d.rate <- NA
 	sim.dat$i.rate <- NA
 	sim.dat$m.rate <- NA
      
 	sim.dat$N <- NA
 	sim.dat$S.i <- NA
 	sim.dat$S.m <- NA
 	sim.dat$I <- NA
 	sim.dat$M <- NA
 	sim.dat$D <- NA

 	sim.dat$I.est <- NA  # this can be estimated catch or estimated impacts in ocean fisheries 
 	sim.dat$M.est <- NA


	# Age 1
 	sim.dat$N[sim.dat$AGE==1]	 <- sim.dat$R[sim.dat$AGE==1]
 	sim.dat$I[sim.dat$AGE==1]	 <- 0
      sim.dat$S.i[sim.dat$AGE==1]	 <- 0
 	sim.dat$I.est[sim.dat$AGE==1]	 <- 0
      sim.dat$i.rate[sim.dat$AGE==1] <- 0
 	sim.dat$M[sim.dat$AGE==1]	 <- 0
      sim.dat$S.m[sim.dat$AGE==1]	 <- 0
 	sim.dat$M.est[sim.dat$AGE==1]	 <- 0
      sim.dat$m.rate[sim.dat$AGE==1] <- 0
 	sim.dat$D[sim.dat$AGE==1]	 <-  rbinom( numbroods, size = sim.dat$N[sim.dat$AGE==1] - sim.dat$I[sim.dat$AGE==1] - sim.dat$M[sim.dat$AGE==1], 	
										prob = gen.rate.dat$d[gen.rate.dat$AGE==1])
      
      sim.dat$d.rate[sim.dat$AGE==1] <- sim.dat$D[sim.dat$AGE==1]/(sim.dat$N[sim.dat$AGE==1] - sim.dat$I[sim.dat$AGE==1]- sim.dat$M[sim.dat$AGE==1])


      a = ages_c[3]

  	for( a in ages_c){

	 	sim.dat$N[sim.dat$AGE==a] <- sim.dat$N[sim.dat$AGE==(a-1)] - sim.dat$I[sim.dat$AGE==(a-1)]- sim.dat$M[sim.dat$AGE==(a-1)] - sim.dat$D[sim.dat$AGE==(a-1)]

	 	sim.dat$I[sim.dat$AGE==a] <- rbinom( numbroods, size = sim.dat$N[sim.dat$AGE==a], 	
										prob = gen.rate.dat$i[gen.rate.dat$AGE==a])

	      sim.dat$i.rate[sim.dat$AGE==a] <- sim.dat$I[sim.dat$AGE==a]/sim.dat$N[sim.dat$AGE==a]



	 	sim.dat$S.i[sim.dat$AGE==a] <- rbinom( numbroods, size = sim.dat$I[sim.dat$AGE==a], 	
										prob = ocean.sample.fraction)

	 	sim.dat$I.est[sim.dat$AGE==a] <- sim.dat$S.i[sim.dat$AGE==a]/ocean.sample.fraction



	 	sim.dat$M[sim.dat$AGE==a] <- rbinom( numbroods, size = sim.dat$N[sim.dat$AGE==a] - sim.dat$I[sim.dat$AGE==a], 	
										prob = gen.rate.dat$m[gen.rate.dat$AGE==a])

	      sim.dat$m.rate[sim.dat$AGE==a] <- sim.dat$M[sim.dat$AGE==a]/(sim.dat$N[sim.dat$AGE==a] - sim.dat$I[sim.dat$AGE==a])


	 	sim.dat$S.m[sim.dat$AGE==a] <- rbinom( numbroods, size = sim.dat$M[sim.dat$AGE==a], 	
										prob = river.sample.fraction)


	 	sim.dat$M.est[sim.dat$AGE==a] <- sim.dat$S.m[sim.dat$AGE==a]/river.sample.fraction


	 	sim.dat$D[sim.dat$AGE==a] <- rbinom( numbroods, size = sim.dat$N[sim.dat$AGE==a] - sim.dat$I[sim.dat$AGE==a] - sim.dat$M[sim.dat$AGE==a], 	
										prob = gen.rate.dat$d[gen.rate.dat$AGE==a])

	      sim.dat$d.rate[sim.dat$AGE==a] <- sim.dat$D[sim.dat$AGE==a]/(sim.dat$N[sim.dat$AGE==a] - sim.dat$I[sim.dat$AGE==a]- sim.dat$M[sim.dat$AGE==a])

      }

   return(sim.dat)

  } # end of function


