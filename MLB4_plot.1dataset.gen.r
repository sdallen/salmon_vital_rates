
 rm(list = ls(all = TRUE))

# ==============================================================================
# Directories/Files.


# =============================================================================
# functions



# ==============================================================================
# User defined file names and quantities


# Specify model
  model <- "mlb4"

#  admb file name
  admb.file <- "mlb4_no_s_c4" 

#  age 5 maturation rate 
  m5constant <- 1  

 
 # Specify the set of generating rates 
 
  gen.rates <- "CR_tv_GenRates_IGH_Y_no_missBY" # This set of generating rate does NOT satisfy model assumptions
# gen.rates <- "mu_tv_GenRates_IGH_Y_no_missBY" # This set of generating rates satisfies all model assumptions
 
# Specify the dataset file (no extension)
  
 #data <- "sim.data.no.missing.cohorts" #stochastic data sets
  data <- "det.data.no.missing.cohorts" #deterministic (noise-free) data sets

# Specify location of dataset file 

  data.fn <- paste(  gen.rates, paste(data, ".txt", sep=""), sep="/")

# Specify inpu file
  mlb4.fn <- paste(model, ".output.", gen.rates, ".", data, ".txt", sep="")

 




# =============================================================================
# Read data file


  mlb4 <- read.table( mlb4.fn, header=TRUE, sep=",", as.is=TRUE, colClasses="numeric")
  colnames(mlb4)


# =============================================================================
# PLOT all possible solutions for one particular dataset and GENERATING Rates 

 # find datasets that have estimates for at least one iteration

  
  sum.perdata <- aggregate( m.est ~ i.data, data=mlb4, subset = AGE==2, sum) 

  data.est <- unique( mlb4$i.data)[ !(is.na(sum.perdata$m.est))]

  k = sample(data.est, 1) # 4 #94 # pick dataset number randomly (will this code work when there is only 1 dataset? like when fitting a deterministic dataset?)
  k

 # specify which iterations to plot


 #plot only iterations that have a maximum gradient component that's less than tol (per dataset)
    # and out of those select iterations with the minimum -LL value 

 tol <- .0001  #.001 #- TSP default mccullough  #

 iter.crit <- unique( mlb4$i.iter[ mlb4$maxgrad <= tol & !is.na(mlb4$maxgrad) & mlb4$i.data==k & 
						mlb4$LL == min( mlb4$LL[ mlb4$maxgrad <=  tol & mlb4$i.data==k  ]) ] )
 length( iter.crit) 


 #plot all iterations

 iter.crit <- unique(mlb4$i.iter[ mlb4$i.data==k])
 



 # Pick colors for each iteration

 col.itercrit = rainbow(length(iter.crit))

								
 par( mfrow=c( 3, 1 ),
  		mar=c(2,4.5,1.7,.5), #inner margins
	     	mgp=c(3,.7,0), #how close axis labels are to graph
	      oma=c(3,0,3,1),  #outer-most margins #c(bottom, left, top, right)c(4,4,4,1)
	    	err=-1, bg="white") 
		#"din" )  

 tf <- mlb4$i.data==k & mlb4$i.iter== iter.crit[1] 
 plot( mlb4$BY[tf & mlb4$AGE==2] + 2, mlb4$i.est[tf & mlb4$AGE==2], pch=16, type="b", ylim=c(0,1),
		xlim=range( c( mlb4$BY + min(mlb4$AGE), mlb4$BY + max(mlb4$AGE)) ),  xlab="", ylab="Exploitation rate", col=col.itercrit[1])
 lines( mlb4$BY[tf & mlb4$AGE==3] + 3, mlb4$i.est[tf & mlb4$AGE==3], pch=16, type="b", col=col.itercrit[1])
 lines( mlb4$BY[tf & mlb4$AGE==4] + 4, mlb4$i.est[tf & mlb4$AGE==4], pch=16, type="b", col=col.itercrit[1])	
# lines( mlb4$BY[tf & mlb4$AGE==5] + 5, mlb4$i.est[tf & mlb4$AGE==5], pch=16, type="b", col=col.itercrit[1])	

# mtext(Model.title, 
#				side = 3, outer=TRUE, line=0)	
  title( main=paste( 'MLB4', 'dataset #', k, sep= " "))

  legend('topright', c('Age specific generating rate'), lty=2, col=c('black'), cex=.8)
				
 for ( i in iter.crit[-1]){
 	tf <- mlb4$i.data==k & mlb4$i.iter==i      
  	lines( mlb4$BY[tf & mlb4$AGE==2] + 2, mlb4$i.est[tf & mlb4$AGE==2], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])		
 	lines( mlb4$BY[tf & mlb4$AGE==3] + 3, mlb4$i.est[tf & mlb4$AGE==3], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])
 	lines( mlb4$BY[tf & mlb4$AGE==4] + 4, mlb4$i.est[tf & mlb4$AGE==4], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])
 #	lines( mlb4$BY[tf & mlb4$AGE==5] + 5, mlb4$i.est[tf & mlb4$AGE==5], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])

 }

 	tf <- mlb4$i.data==k & mlb4$i.iter==1   
  	lines( mlb4$BY[tf & mlb4$AGE==2] + 2, mlb4$i[tf & mlb4$AGE==2], pch=16, lty=2)		
 	lines( mlb4$BY[tf & mlb4$AGE==3] + 3, mlb4$i[tf & mlb4$AGE==3], pch=16, lty=2)
 	lines( mlb4$BY[tf & mlb4$AGE==4] + 4, mlb4$i[tf & mlb4$AGE==4], pch=16, lty=2)
 #	lines( mlb4$BY[tf & mlb4$AGE==5] + 5, mlb4$i[tf & mlb4$AGE==5], pch=16, lty=2)




# Plot Maturation rates

 tf <- mlb4$i.data==k & mlb4$i.iter== iter.crit[1] 
 plot( mlb4$BY[tf & mlb4$AGE==2] + 2, mlb4$m.est[tf & mlb4$AGE==2], pch=16, type="b", ylim=c(0,1),
		xlim=range( c( mlb4$BY + min(mlb4$AGE), mlb4$BY + max(mlb4$AGE)) ),  xlab="", ylab="Maturation rate", col=col.itercrit[1])
 lines( mlb4$BY[tf & mlb4$AGE==3] + 3, mlb4$m.est[tf & mlb4$AGE==3], pch=16, type="b", col=col.itercrit[1])
 lines( mlb4$BY[tf & mlb4$AGE==4] + 4, mlb4$m.est[tf & mlb4$AGE==4], pch=16, type="b", col=col.itercrit[1])	
# lines( mlb4$BY[tf & mlb4$AGE==5] + 5, mlb4$m.est[tf & mlb4$AGE==5], pch=16, type="b", col=col.itercrit[1])	

				
 for ( i in iter.crit[-1]){
 	tf <- mlb4$i.data==k & mlb4$i.iter==i      
  	lines( mlb4$BY[tf & mlb4$AGE==2] + 2, mlb4$m.est[tf & mlb4$AGE==2], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])		
 	lines( mlb4$BY[tf & mlb4$AGE==3] + 3, mlb4$m.est[tf & mlb4$AGE==3], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])
 	lines( mlb4$BY[tf & mlb4$AGE==4] + 4, mlb4$m.est[tf & mlb4$AGE==4], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])
# 	lines( mlb4$BY[tf & mlb4$AGE==5] + 5, mlb4$m.est[tf & mlb4$AGE==5], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])

 }

 	tf <- mlb4$i.data==k & mlb4$i.iter==1  
  	lines( mlb4$BY[tf & mlb4$AGE==2] + 2, mlb4$m[tf & mlb4$AGE==2], pch=16, lty=2)		
 	lines( mlb4$BY[tf & mlb4$AGE==3] + 3, mlb4$m[tf & mlb4$AGE==3], pch=16, lty=2)
 	lines( mlb4$BY[tf & mlb4$AGE==4] + 4, mlb4$m[tf & mlb4$AGE==4], pch=16, lty=2)
# 	lines( mlb4$BY[tf & mlb4$AGE==5] + 5, mlb4$m[tf & mlb4$AGE==5], pch=16, lty=2)



 tf <- mlb4$i.data==k & mlb4$i.iter== iter.crit[1] 
 plot( mlb4$BY[tf & mlb4$AGE==1] + 1, mlb4$d.est[tf & mlb4$AGE==1], pch=16, type="b", ylim=c(0,1),
		xlim=range( c( mlb4$BY + min(mlb4$AGE), mlb4$BY + max(mlb4$AGE)) ),  xlab="", ylab="Natural mortality rate", col=col.itercrit[1])

 lines( mlb4$BY[tf & mlb4$AGE==2] + 2, mlb4$d.est[tf & mlb4$AGE==2], pch=16, type="b", col=col.itercrit[1])
 lines( mlb4$BY[tf & mlb4$AGE==3] + 3, mlb4$d.est[tf & mlb4$AGE==3], pch=16, type="b", col=col.itercrit[1])	
# lines( mlb4$BY[tf & mlb4$AGE==4] + 4, mlb4$d.est[tf & mlb4$AGE==4], pch=16, type="b", col=col.itercrit[1])	 

				
 for ( i in iter.crit[-1]){
 	tf <- mlb4$i.data==k & mlb4$i.iter==i      
  	lines( mlb4$BY[tf & mlb4$AGE==1] + 1, mlb4$d.est[tf & mlb4$AGE==1], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])		
 	lines( mlb4$BY[tf & mlb4$AGE==2] + 2, mlb4$d.est[tf & mlb4$AGE==2], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])
 	lines( mlb4$BY[tf & mlb4$AGE==3] + 3, mlb4$d.est[tf & mlb4$AGE==3], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])
# 	lines( mlb4$BY[tf & mlb4$AGE==4] + 4, mlb4$d.est[tf & mlb4$AGE==4], pch=16, type="b", col=col.itercrit[match(i,iter.crit)])
 }

 	tf <- mlb4$i.data==k & mlb4$i.iter== 1  
  	lines( mlb4$BY[tf & mlb4$AGE==1] + 1, mlb4$d[tf & mlb4$AGE==1], pch=16, lty=2)		
 	lines( mlb4$BY[tf & mlb4$AGE==2] + 2, mlb4$d[tf & mlb4$AGE==2], pch=16, lty=2)
 	lines( mlb4$BY[tf & mlb4$AGE==3] + 3, mlb4$d[tf & mlb4$AGE==3], pch=16, lty=2)
 #	lines( mlb4$BY[tf & mlb4$AGE==4] + 4, mlb4$d[tf & mlb4$AGE==4], pch=16, lty=2)




