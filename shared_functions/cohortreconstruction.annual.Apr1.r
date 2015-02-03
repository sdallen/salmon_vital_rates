# This ANNUAL cohort reconstruction assumes catch -> maturation -> natural mortality

cohortreconstruction.annual.Apr1 <- function( v, dat){

  
  dat$N.cr.est[dat$AGE==5] <- dat$I.est[dat$AGE==5] + dat$M.est[dat$AGE==5] + 0

  dat$N.cr.est[dat$AGE==4] <- dat$I.est[dat$AGE==4] + dat$M.est[dat$AGE==4] + dat$N.cr.est[dat$AGE==5]/(1- v["4"])

  dat$N.cr.est[dat$AGE==3] <-  dat$I.est[dat$AGE==3] + dat$M.est[dat$AGE==3] + dat$N.cr.est[dat$AGE==4]/(1- v["3"])

  dat$N.cr.est[dat$AGE==2] <-  dat$I.est[dat$AGE==2] + dat$M.est[dat$AGE==2] + dat$N.cr.est[dat$AGE==3]/(1- v["2"])

  dat$N.cr.est[dat$AGE==1] <- dat$R[dat$AGE==1]

  # conditional exploitation rates and maturation rates
  
  dat$i.cr.est[dat$AGE==1] <- 0

  dat$i.cr.est[dat$AGE==2] <- dat$I.est[dat$AGE==2]/dat$N.cr.est[dat$AGE==2]

  dat$i.cr.est[dat$AGE==3] <- dat$I.est[dat$AGE==3]/dat$N.cr.est[dat$AGE==3]

  dat$i.cr.est[dat$AGE==4] <- dat$I.est[dat$AGE==4]/dat$N.cr.est[dat$AGE==4]

  dat$i.cr.est[dat$AGE==5] <- dat$I.est[dat$AGE==5]/dat$N.cr.est[dat$AGE==5]


  # conditional maturation rates

  dat$m.cr.est[dat$AGE==1] <- 0

  dat$m.cr.est[dat$AGE==2] <- dat$M.est[dat$AGE==2]/(dat$M.est[dat$AGE==2] + dat$N.cr.est[dat$AGE==3]/(1 - v["2"]))

  dat$m.cr.est[dat$AGE==3] <- dat$M.est[dat$AGE==3]/(dat$M.est[dat$AGE==3] + dat$N.cr.est[dat$AGE==4]/(1 - v["3"]))

  dat$m.cr.est[dat$AGE==4] <- dat$M.est[dat$AGE==4]/(dat$M.est[dat$AGE==4] + dat$N.cr.est[dat$AGE==5]/(1 - v["4"]))

  dat$m.cr.est[dat$AGE==5] <- dat$M.est[dat$AGE==5]/(dat$M.est[dat$AGE==5] + 0)


  # conditional natural mortality rates

  dat$d.cr.est[dat$AGE==1] <- 1 - dat$N.cr.est[dat$AGE==2]/dat$R[dat$AGE==2]
  dat$d.cr.est[dat$AGE==2] <- v["2"]
  dat$d.cr.est[dat$AGE==3] <- v["3"]
  dat$d.cr.est[dat$AGE==4] <- v["4"]  
  

  return( dat )
}


