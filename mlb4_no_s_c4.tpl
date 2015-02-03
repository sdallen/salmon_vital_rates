  //This model assumes fully vulnerable fishing mortality is constant over years 
  //This model assumes age specific natural mortality rates share the same calendar year effects
  //This model assumes age specific maturation rates share the same cohort effect
	
  //mnow - unconditional maturation probability
  //cnow - unconditional exploitation probability
	//pivec - a vector of unconditional probabilities (mnow, cnow)
	//expected - a vector of expected catch and escapement numbers (pivec*N) 
  //d1 - a vector of estimated age 1 natural mortality 
  //d2 - a vector of estimated age 2-4 natural mortality 		
	//c_jk - a matrix of estimated age k and calendar year j exploitation rate
  //m2 - a vector of estimated age 2 maturation rates
  //m3 - a vector of estimated age 3 maturation rates 		
  //m4 - a vector of estimated age 4 maturation rates 	
	//s - a vector of estimated confounding probabilities
	
  //npostd_k - estimated abundance post period of natural mortality during age k
  //npostc_k - estimated abundance post period of exploitation during age k
  //npostm_k - estimated abundance post period of maturation during age k

  //allows for missing brood years
	
PARAMETER_SECTION

  sdreport_matrix d_jk(1,nestcalyrs_d,1,nages_d)

  sdreport_matrix c_jk(1,nestcalyrs_c,1,nages_c)
  sdreport_matrix m_jk(1,nestbys,1,nages_c-1) //don't include age 5 because standard errors don't exist, watch spaces

  vector pivec(1,nobs)    // unconditional catch/maturation probabilities
  vector expected(1,nobs)
  number cnow
  number mnow


  vector npostd_1(1,nestbys)
  vector npostd_2(1,nestbys)
  vector npostd_3(1,nestbys) 
  vector npostd_4(1,nestbys) 

  vector npostc_2(1,nestbys)
  vector npostc_3(1,nestbys)
  vector npostc_4(1,nestbys) 
  vector npostc_5(1,nestbys) 

  vector npostm_2(1,nestbys)
  vector npostm_3(1,nestbys)
  vector npostm_4(1,nestbys) 
  vector npostm_5(1,nestbys) 

PROCEDURE_SECTION

    // penalty variable if probabilities are outside of some range
	dvariable fpen=0.0;         
	
	
	// use complementary log-log transformation to obtain maturation rates,exploitation and natural mortality rates


	//maturation rates
	// 1st brood year/index 1 

        m_jk(1,1) = 1 - mfexp(-mfexp( logbm_jk(1) ) );  //age 2
        m_jk(1,2) = 1 - mfexp(-mfexp( logbm_jk(1) + logbm_jk(2)  ) );  //age3
        m_jk(1,3) = 1 - mfexp(-mfexp( cllm4 ) );   //age4


	//for 2nd and onward brood years:  

	for (int j=2;j<=nestbys;j++){
		int m = j-1;
		// age 2
		m_jk(j,1)= 1 - mfexp(-mfexp( logbm_jk(1) + logbm_jk(m+2) ) ); 
		// age 3
		m_jk(j,2)= 1 - mfexp(-mfexp( logbm_jk(1) + logbm_jk(2) + logbm_jk(m+2) ) ); 
		// age 4
		m_jk(j,3)= 1 - mfexp(-mfexp( cllm4 ) );   //age4
		
	}	






		
      //exploitation rate
      // year effects are fully vulnerable fishing mortality =  CLL(age 4/5)
      // additive age 2 and age 3 effects (catchability)

	//age 4/5
	for (int j=1;j<=nestcalyrs_c;j++){
		
		c_jk(j,3) = 1 - mfexp(-mfexp(  logbc_jk(j+2) ) );
		c_jk(j,4) = c_jk(j,3);

	}	


	//age 2

	for (int j=1;j<=nestcalyrs_c;j++){

		c_jk(j,1) = 1 - mfexp(-mfexp( logbc_jk(1) +  logbc_jk(j+2) ) );
		
		
	}	

	//age 3

	for (int j=1;j<=nestcalyrs_c;j++){

		c_jk(j,2) = 1 - mfexp(-mfexp( logbc_jk(2) +  logbc_jk(j+2) ) );
		
		
	}	


		
		
    //natural mortality rate


    //baseline calendar year - index 1

	d_jk(1,1) = 1 - mfexp(-mfexp( logbd_jk(1)  ) );   //age 1
	d_jk(1,2) = 1 - mfexp(-mfexp( logbd_jk(1) + logbd_jk(2)  ) );   //age 2
	d_jk(1,3) = 1 - mfexp(-mfexp( logbd_jk(1) + logbd_jk(2)  ) );   //age 3
	d_jk(1,4) = 1 - mfexp(-mfexp( logbd_jk(1) + logbd_jk(2)  ) );   //age 4

	for (int l=2;l<=nestcalyrs_d;l++){
		int m = l-1;

		d_jk(l,1) = 1 - mfexp(-mfexp( logbd_jk(1) + logbd_jk(m+2) ) );   //age 1
		d_jk(l,2) = 1 - mfexp(-mfexp( logbd_jk(1) + logbd_jk(2) + logbd_jk(m+2) ) );   //age 2
		d_jk(l,3) = 1 - mfexp(-mfexp( logbd_jk(1) + logbd_jk(2) + logbd_jk(m+2) ) );   //age 3
		d_jk(l,4) = 1 - mfexp(-mfexp( logbd_jk(1) + logbd_jk(2) + logbd_jk(m+2) ) );   //age 4
		
	}



	
	
	// Start 
	
	for (int bindex=1;bindex<=nestbys;bindex++){ 
		for(int aindex=1;aindex<=nages_c;aindex++){

		//determine calendar year
		int a = ages_c(aindex);
		int calyr_at = ages_c(aindex) + estbys(bindex);
		int calyrc_ind = calyr_at - minestcalyr_c +1;
		int calyrd_ind = calyr_at-1 - minestcalyr_d +1;
		
		//determine c, record in row.at of pi
		//determine m, record in row.at+num.rows of pi
			if (a==2){
			
			        //survive to 2           get fished at 2
			cnow= (1-d_jk(calyrd_ind,1)) * c_jk(calyrc_ind,aindex);
			        //survive to 2           don't get fished at 2      mature at 2
			mnow= (1-d_jk(calyrd_ind,1)) * (1-c_jk(calyrc_ind,aindex)) *  m_jk(bindex,1);

			npostd_1(bindex) = nvec((bindex-1)*nages_c + aindex)*(1-d_jk(calyrd_ind,1)) ;
			npostc_2(bindex) = npostd_1(bindex) * (1-c_jk(calyrc_ind,aindex));
			npostm_2(bindex) = npostc_2(bindex) *  (1 - m_jk(bindex,1));
 
			}
			if (a==3){
			
		                  //survive to 2     don't get fished at 2                don't mature at 2   survive to 3          get fished at 3
			cnow= (1-d_jk(calyrd_ind-1,1))  * (1-c_jk(calyrc_ind-1,aindex-1)) * (1-m_jk(bindex,1)) * (1-d_jk(calyrd_ind,2)) * c_jk(calyrc_ind,aindex);
			        //survive to 2         don't get fished at 2   don't mature at 2     survive to 3     don't get fished at 3      mature at 3
			mnow= (1-d_jk(calyrd_ind-1,1))  * (1-c_jk(calyrc_ind-1,aindex-1)) * (1-m_jk(bindex,1)) * (1-d_jk(calyrd_ind,2))  * (1-c_jk(calyrc_ind,aindex)) * m_jk(bindex,2);
	
			npostd_2(bindex) = npostm_2(bindex) * (1-d_jk(calyrd_ind,2));
			npostc_3(bindex) = npostd_2(bindex) * (1-c_jk(calyrc_ind,aindex));
			npostm_3(bindex) = npostc_3(bindex) * (1-m_jk(bindex,2));


			}
			if (a==4){
			
			        //survive to 2                 don't get fished at 2        don't mature at 2     survive to 3            don't get fished at 3        don't mature at 3      survive to 4            get fished at 4
			cnow= (1-d_jk(calyrd_ind-2,1))  * (1-c_jk(calyrc_ind-2,aindex-2)) * (1-m_jk(bindex,1)) * (1-d_jk(calyrd_ind-1,2))  *  (1-c_jk(calyrc_ind-1,aindex-1)) * (1-m_jk(bindex,2)) * (1-d_jk(calyrd_ind,3)) * c_jk(calyrc_ind,aindex);
			        //survive to 2           don't get fished at 2        don't mature at 2     survive to 3                     don't get fished at 3        don't mature at 3      survive to 4          don't get fished at 4      mature at 4
			mnow= (1-d_jk(calyrd_ind-2,1))  * (1-c_jk(calyrc_ind-2,aindex-2)) * (1-m_jk(bindex,1)) * (1-d_jk(calyrd_ind-1,2))  *  (1-c_jk(calyrc_ind-1,aindex-1)) * (1-m_jk(bindex,2)) * (1-d_jk(calyrd_ind,3))* (1-c_jk(calyrc_ind,aindex)) * m_jk(bindex,3);
	
			npostd_3(bindex) = npostm_3(bindex) * (1-d_jk(calyrd_ind,3)) ;
			npostc_4(bindex) = npostd_3(bindex) * (1-c_jk(calyrc_ind,aindex));
			npostm_4(bindex) = npostc_4(bindex) * (1 - m_jk(bindex,3));
							
			}
			if (a==5){
			
			        //survive to 2           don't get fished at 2        don't mature at 2             survive to 3           don't get fished at 3        don't mature at 3      survive to 4           don't get fished at 4          don't mature at 4      survive to 5                  get fished at 5
			cnow= (1-d_jk(calyrd_ind-3,1))  * (1-c_jk(calyrc_ind-3,aindex-3)) * (1-m_jk(bindex,1)) * (1-d_jk(calyrd_ind-2,2))  *  (1-c_jk(calyrc_ind-2,aindex-2)) * (1-m_jk(bindex,2)) * (1-d_jk(calyrd_ind-1,3))* (1-c_jk(calyrc_ind-1,aindex-1)) * (1-m_jk(bindex,3)) * (1-d_jk(calyrd_ind,4))* c_jk(calyrc_ind,aindex);
			        //survive to 2           don't get fished at 2        don't mature at 2     survive to 3                     don't get fished at 3        don't mature at 3      survive to 4                  don't get fished at 4          don't mature at 4      survive to 5                  don't get fished at 5      mature at 5
			mnow= (1-d_jk(calyrd_ind-3,1)) * (1-c_jk(calyrc_ind-3,aindex-3)) * (1-m_jk(bindex,1))  * (1-d_jk(calyrd_ind-2,2))  *  (1-c_jk(calyrc_ind-2,aindex-2)) * (1-m_jk(bindex,2)) * (1-d_jk(calyrd_ind-1,3))* (1-c_jk(calyrc_ind-1,aindex-1)) * (1-m_jk(bindex,3)) * (1-d_jk(calyrd_ind,4)) * (1-c_jk(calyrc_ind,aindex)) * m5constant;
	
			npostd_4(bindex) = npostm_4(bindex) * (1-d_jk(calyrd_ind,4)) ;
			npostc_5(bindex) = npostd_4(bindex) * (1-c_jk(calyrc_ind,aindex));
			npostm_5(bindex) = npostc_5(bindex) * (1 - m5constant);
		
			}

		
		//penalties to ensure pivec is greater than 0 
		fpen = 0;
		cnow = posfun(cnow,1.0e-10,fpen);
		f += 1.0*fpen;
		
		fpen = 0;
		mnow = posfun(mnow,1.0e-10,fpen);
		f += 1.0*fpen;

					
		pivec((bindex-1)*nages_c + aindex)=cnow;
		pivec(nages_c*nestbys + (bindex-1)*nages_c + aindex)=mnow;

		
	} // END OF AGE
    } // END OF BROODYR


	expected = elem_prod(nvec,pivec);  //this must be greater than zero because it's later logged
	// poisson negative log-likelihood
	f -= sum( elem_prod(yvec,log(expected) )  - expected - gammln(yvec+1.0) ); //because gamma function is defined as(n-1)!
  
	
	//list added penalties here (e.g.,  constrain logit transformation of extreme probablities to avoid derivatives near 0)
	
		
		
RUNTIME_SECTION
    maximum_function_evaluations 20000
    convergence_criteria 0.0001  // 0.000001  //default value = 0.0001, this is the same as setting the maximum gradient component criterion
	
	
TOP_OF_MAIN_SECTION

	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);	//need this when nestbys > 8
