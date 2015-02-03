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
# ==============================================================================
  det.data <- function( gen.rate.dat,  ocean.sample.fraction, 
			 river.sample.fraction, broodyrs,calyrs, ages_c,
			  numbroods){


      det.dat <- gen.rate.dat

 	det.dat$d.rate <- NA
 	det.dat$i.rate <- NA
 	det.dat$m.rate <- NA
      

 	det.dat$N <- NA
 	det.dat$S.i <- NA
 	det.dat$S.m <- NA
 	det.dat$I <- NA # this can be estimated catch or estimated impacts
 	det.dat$M <- NA
 	det.dat$I.est <- NA  # this can be estimated catch or estimated impacts
 	det.dat$M.est <- NA
 	det.dat$D <- NA


	# Age 1
 	det.dat$N[det.dat$AGE==1]	 <- det.dat$R[det.dat$AGE==1]
 	det.dat$I[det.dat$AGE==1]	 <- 0
      det.dat$S.i[det.dat$AGE==1]	 <- 0
      det.dat$i.rate[det.dat$AGE==1] <- 0

 	det.dat$M[det.dat$AGE==1]	 <- 0
      det.dat$S.m[det.dat$AGE==1]	 <- 0
 	det.dat$I.est[det.dat$AGE==1]	 <- 0
 	det.dat$M.est[det.dat$AGE==1]	 <- 0
      det.dat$m.rate[det.dat$AGE==1] <- 0

 	det.dat$D[det.dat$AGE==1]	 <-  (det.dat$N[det.dat$AGE==1] - det.dat$I[det.dat$AGE==1] - det.dat$M[det.dat$AGE==1])* 	
										(gen.rate.dat$d[gen.rate.dat$AGE==1])
      
      det.dat$d.rate[det.dat$AGE==1] <- det.dat$D[det.dat$AGE==1]/(det.dat$N[det.dat$AGE==1] - det.dat$I[det.dat$AGE==1]- det.dat$M[det.dat$AGE==1])



  	for( a in ages_c){

	 	det.dat$N[det.dat$AGE==a] <- det.dat$N[det.dat$AGE==(a-1)] - det.dat$I[det.dat$AGE==(a-1)]- det.dat$M[det.dat$AGE==(a-1)] - det.dat$D[det.dat$AGE==(a-1)]

	 	det.dat$I[det.dat$AGE==a] <- det.dat$N[det.dat$AGE==a]*gen.rate.dat$i[gen.rate.dat$AGE==a]

	      det.dat$i.rate[det.dat$AGE==a] <- det.dat$I[det.dat$AGE==a]/det.dat$N[det.dat$AGE==a]


	 	det.dat$S.i[det.dat$AGE==a] <- det.dat$I[det.dat$AGE==a]*ocean.sample.fraction

	 	det.dat$I.est[det.dat$AGE==a] <- det.dat$S.i[det.dat$AGE==a]/ocean.sample.fraction   


	 	det.dat$M[det.dat$AGE==a] <- (det.dat$N[det.dat$AGE==a] - det.dat$I[det.dat$AGE==a])*gen.rate.dat$m[gen.rate.dat$AGE==a]

	 	det.dat$S.m[det.dat$AGE==a] <- det.dat$M[det.dat$AGE==a]*river.sample.fraction

	      det.dat$m.rate[det.dat$AGE==a] <- det.dat$M[det.dat$AGE==a]/(det.dat$N[det.dat$AGE==a] - det.dat$I[det.dat$AGE==a])


	 	det.dat$M.est[det.dat$AGE==a] <- det.dat$S.m[det.dat$AGE==a]/river.sample.fraction


	 	det.dat$D[det.dat$AGE==a] <-  (det.dat$N[det.dat$AGE==a] - det.dat$I[det.dat$AGE==a] - det.dat$M[det.dat$AGE==a])*gen.rate.dat$d[gen.rate.dat$AGE==a]

	      det.dat$d.rate[det.dat$AGE==a] <- det.dat$D[det.dat$AGE==a]/(det.dat$N[det.dat$AGE==a] - det.dat$I[det.dat$AGE==a]- det.dat$M[det.dat$AGE==a])

      }

   return(det.dat)

  } # end of function


