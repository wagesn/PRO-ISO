library(Iso)

onetitecdp<-function (PIc, PIp, targetc, targetp, gam, n, obswinc, obswinp, rate, 
	accrual, minwin, cl) 
{

    	ndose<-length(PIc)
    	ycvec <- ypvec <- nj <- mcj <- mpj <- ecj <- epj <- dose.select <- ptox.hatc <- ptox.hatp <- rep(0,ndose)
    	len <- length(ycvec)
    	userind <- seq(1, len, 1)
    	stopc <- stopp <-0

	###Prior for pi_c
	xc<-2*targetc
	muc<-targetc
	uc<-0.95
	f<-function(b){
		pbeta(xc,muc*b/(1-muc),b)-uc
	}
	b0c<-uniroot(f,c(0.0001,100))$root
	a0c<-muc*b0c/(1-muc)

	###Prior for pi_p
	xp<-2*targetp
	mup<-targetp
	up<-0.95
	f<-function(b){
		pbeta(xp,mup*b/(1-mup),b)-up
	}
	b0p<-uniroot(f,c(0.0001,100))$root
	a0p<-mup*b0p/(1-mup)

	post.tox<-function(p,a,b,y,m,n,w,delta){
        term=1
	  prior=(((p^(a-1))*(1-p)^(b-1))/(beta(a,b)))
	  for(i in 1:n){
		term=term*(1-w[i]*p)**(1-delta[i])
	  }
	  return(term*(p^y)*((1-p)^(m))*prior);
      }
      
      #the posterior mean of ptox
      posttoxf <- function(p,a,b,y,m,n,w,delta) {p*post.tox(p,a,b,y,m,n,w,delta); }

    	if (accrual == "fixed") {
        next.arrival <- obswinc/rate
    	} else if (accrual == "poisson") {
        next.arrival <- rexp(1, rate/obswinc)
    	}
        
	  uc <- up <- yc <- yp <- level <- arrival <- weightsc <- weightsp <- NULL
        cur <- 1
        
        for (i in 1:(n - 1)) {
            arrival <- c(arrival, next.arrival)
            level <- c(level, cur)
         
            ycnew <- rbinom(1, 1, PIc[cur])
            if (ycnew){ 
                	ucnew <- runif(1, 0, obswinc)
            	} else {
				ucnew <- Inf
			}
            	yc <- c(yc, ycnew)
            	uc <- c(uc, ucnew)
            	uctox <- uc + arrival

           ypnew <- rbinom(1, 1, PIp[cur])
            if (ypnew) {
                	upnew <- runif(1, 0, obswinp)
            	} else {
				upnew <- Inf
			}
            	yp <- c(yp, ypnew)
            	up <- c(up, upnew)
            	uptox <- up + arrival
        	        
            if (accrual == "fixed") {
                next.arrival <- next.arrival + obswinc/rate
            } else if (accrual == "poisson") {
                next.arrival <- next.arrival + rexp(1, rate/obswinc)
            }

            Bc <- rep(0, length(yc))
            Bc[uctox <= next.arrival] <- 1
            censorc <- pmin(next.arrival, uctox) - arrival
            followupc <- pmin(censorc, obswinc)

            Bp <- rep(0, length(yp))
            Bp[uptox <= next.arrival] <- 1
            censorp <- pmin(next.arrival, uptox) - arrival
            followupp <- pmin(censorp, obswinp)

            weightsc <- followupc/obswinc
            weightsp <- followupp/obswinp

   		weightsc[Bc == 1] <- 1
    		y1c <- Bc[1:length(Bc)]
    		w1c <- weightsc[1:length(Bc)]

		weightsp[Bp == 1] <- 1
    		y1p <- Bp[1:length(Bp)]
    		w1p <- weightsp[1:length(Bp)]
		
   	 for(j in 1:length(ycvec)){
		ycvec[j]=sum(y1c[level==j])	
		nj[j]=sum(level==j)
		mcj[j] <- sum(w1c[level==j])#sum(w1c[level==j]==1 & y1c[level==j]==0)
		ecj[j] <- sum(w1c[level==j]>=(minwin/obswinc) & y1c[level==j]==0)
	}

	for(j in 1:length(ypvec)){
		ypvec[j]=sum(y1p[level==j])	
		mpj[j] <- sum(w1p[level==j])#==1 & y1p[level==j]==0)
		epj[j] <- sum(w1p[level==j]>=(minwin/obswinp) & y1p[level==j]==0)
	}


	tried=unique(level)
	#for(j in 1:length(tried)){
	#	marginal.toxc = integrate(post.tox,lower=0,upper=1,a0c,b0c,ycvec[j],mcj[j],nj[j],w1c[level==j],as.numeric(w1c[level==j]==1))$value;
	#	marginal.toxp = integrate(post.tox,lower=0,upper=1,a0p,b0p,ypvec[j],mpj[j],nj[j],w1p[level==j],as.numeric(w1p[level==j]==1))$value;
#
#		ptox.hatc[j] = integrate(posttoxf,lower=0,upper=1,a0c,b0c,ycvec[j],mcj[j],nj[j],w1c[level==j],as.numeric(w1c[level==j]==1))$value/marginal.toxc; 
#		ptox.hatp[j] = integrate(posttoxf,lower=0,upper=1,a0p,b0p,ypvec[j],mpj[j],nj[j],w1p[level==j],as.numeric(w1p[level==j]==1))$value/marginal.toxp; 
#
#	}
	ptox.hatc=qbeta(0.5,ycvec[tried]+a0c,mcj[tried]-ycvec[tried]+b0c)#(ycvec[tried]+a0c)/(mcj[tried]+a0c+b0c)
	ptox.hatp=qbeta(0.5,ypvec[tried]+a0p,mpj[tried]-ypvec[tried]+b0p)#(ypvec[tried]+a0p)/(mpj[tried]+a0p+b0p)

	#marginal.toxc = integrate(post.tox,lower=0,upper=1,a0c,b0c,ycvec[1],mcj[1],nj[1],w1c[level==1],as.numeric(w1c[level==1]==1))$value;
	#marginal.toxp = integrate(post.tox,lower=0,upper=1,a0p,b0p,ypvec[1],mpj[1],nj[1],w1p[level==1],as.numeric(w1p[level==1]==1))$value;
	#unsafe1c <-(integrate(post.tox,lower=targetc,upper=1,a0c,b0c,ycvec[1],mcj[1],nj[1],w1c[level==1],as.numeric(w1c[level==1]==1))$value/marginal.toxc)>cl
	#unsafe1p <-(integrate(post.tox,lower=targetp,upper=1,a0p,b0p,ypvec[1],mpj[1],nj[1],w1p[level==1],as.numeric(w1p[level==1]==1))$value/marginal.toxp)>cl

		if(1 - pbeta(targetc, ycvec[1] + a0c, mcj[1] - ycvec[1] + b0c) > cl){
			stopc<-3 
			break
		}
		if(1 - pbeta(targetp, ypvec[1] + a0p, mpj[1] - ypvec[1] + b0p) > cl){
			stopp<-3 
			break
		}
	#ptox.hatc<-ptox.hatc[tried]
	pipostc=pava(ptox.hatc,w=nj[tried])

	#ptox.hatp<-ptox.hatp[tried]
	pipostp=pava(ptox.hatp,w=nj[tried])
		
		lossvecc=ifelse(pipostc > targetc, (1-gam)*(pipostc-targetc), gam*(targetc-pipostc))  
		minlossc <- min(lossvecc)
		haveminc <- lossvecc==minlossc	
		Tc=lossvecc==min(lossvecc)
			possc=which(Tc)
			if(sum(Tc)==1){
				sugglevc=possc
			} else {
				if(all(pipostc[possc]>targetc)){
					sugglevc=min(possc)
				} else {
					sugglevc=max(possc)
				}
			}
			if(pipostc[sugglevc]<targetc & length(tried)<ndose){
				curc<-ifelse(nj[sugglevc+1]==0&ecj[sugglevc]>0,sugglevc+1,sugglevc)
#&ecj[sugglevc]>0,sugglevc+1,sugglevc)					
			} else {
				curc<-sugglevc
			}         

		lossvecp=ifelse(pipostp > targetp, (1-gam)*(pipostp-targetp), gam*(targetp-pipostp))  
		minlossp <- min(lossvecp)
		haveminp <- lossvecp==minlossp	
		Tp=lossvecp==min(lossvecp)
			possp=which(Tp)
			if(sum(Tp)==1){
				sugglevp=possp
			} else {
				if(all(pipostp[possp]>targetp)){
					sugglevp=min(possp)
				} else {
					sugglevp=max(possp)
				}
			}
			if(pipostp[sugglevp]<targetp & length(tried)<ndose){
				curp<-ifelse(nj[sugglevp+1]==0&epj[sugglevp]>0,sugglevp+1,sugglevp)
#,sugglevp+1,sugglevp)
#&epj[sugglevp]>0,sugglevp+1,sugglevp)					
			} else {
				curp<-sugglevp
			}     
			cur <- min(curc,curp) 
        } #end of for loop


        arrival <- c(arrival, next.arrival)
        level <- c(level, cur)

  		ycnew <- rbinom(1, 1, PIc[cur])
            if (ycnew) 
                	ucnew <- runif(1, 0, obswinc)
            	else ucnew <- Inf
            	yc <- c(yc, ycnew)
            	uc <- c(uc, ucnew)
            	uctox <- uc + arrival

           ypnew <- rbinom(1, 1, PIp[cur])
            if (ypnew) 
                	upnew <- runif(1, 0, obswinp)
            	else upnew <- Inf
            	yp <- c(yp, ypnew)
            	up <- c(up, upnew)
            	uptox <- up + arrival

        	
          	       	weightsc=weightsp=rep(1,length(uctox))	
    		y1c <- yc
    		y1p <- yp

    		w1c <- weightsc
		w1p <- weightsp
    		 
	for(j in 1:length(ycvec)){
		ycvec[j]=sum(y1c[level==j])	
		nj[j]=sum(level==j)
		mcj[j] <- sum(w1c[level==j])#sum(w1c[level==j]==1 & y1c[level==j]==0)
		#ecj[j] <- sum(w1c[level==j]>(minwin/obswinc) & y1c[level==j]==0)
	}

	for(j in 1:length(ypvec)){
		ypvec[j]=sum(y1p[level==j])	
		mpj[j] <- sum(w1p[level==j])#==1 & y1p[level==j]==0)
		#epj[j] <- sum(w1p[level==j]>(minwin/obswinp) & y1p[level==j]==0)
	}


	tried=unique(level)
	#for(j in 1:length(tried)){
	#	marginal.toxc = integrate(post.tox,lower=0,upper=1,a0c,b0c,ycvec[j],mcj[j],nj[j],w1c[level==j],as.numeric(w1c[level==j]==1))$value;
	#	marginal.toxp = integrate(post.tox,lower=0,upper=1,a0p,b0p,ypvec[j],mpj[j],nj[j],w1p[level==j],as.numeric(w1p[level==j]==1))$value;
#
#		ptox.hatc[j] = integrate(posttoxf,lower=0,upper=1,a0c,b0c,ycvec[j],mcj[j],nj[j],w1c[level==j],as.numeric(w1c[level==j]==1))$value/marginal.toxc; 
#		ptox.hatp[j] = integrate(posttoxf,lower=0,upper=1,a0p,b0p,ypvec[j],mpj[j],nj[j],w1p[level==j],as.numeric(w1p[level==j]==1))$value/marginal.toxp; 
#
#	}
	ptox.hatc=qbeta(0.5,ycvec[tried]+a0c,mcj[tried]-ycvec[tried]+b0c)#(ycvec[tried]+a0c)/(mcj[tried]+a0c+b0c)
	ptox.hatp=qbeta(0.5,ypvec[tried]+a0p,mpj[tried]-ypvec[tried]+b0p)#(ypvec[tried]+a0p)/(mpj[tried]+a0p+b0p)

	#marginal.toxc = integrate(post.tox,lower=0,upper=1,a0c,b0c,ycvec[1],mcj[1],nj[1],w1c[level==1],as.numeric(w1c[level==1]==1))$value;
	#marginal.toxp = integrate(post.tox,lower=0,upper=1,a0p,b0p,ypvec[1],mpj[1],nj[1],w1p[level==1],as.numeric(w1p[level==1]==1))$value;
	#unsafe1c <-(integrate(post.tox,lower=targetc,upper=1,a0c,b0c,ycvec[1],mcj[1],nj[1],w1c[level==1],as.numeric(w1c[level==1]==1))$value/marginal.toxc)>cl
	#unsafe1p <-(integrate(post.tox,lower=targetp,upper=1,a0p,b0p,ypvec[1],mpj[1],nj[1],w1p[level==1],as.numeric(w1p[level==1]==1))$value/marginal.toxp)>cl

		if(1 - pbeta(targetc, ycvec[1] + a0c, mcj[1] - ycvec[1] + b0c) > cl){
			stopc<-3 
			
		}
		if(1 - pbeta(targetp, ypvec[1] + a0p, mpj[1] - ypvec[1] + b0p) > cl){
			stopp<-3 
			
		}
	#ptox.hatc<-ptox.hatc[tried]
	pipostc=pava(ptox.hatc,w=nj[tried])

	#ptox.hatp<-ptox.hatp[tried]
	pipostp=pava(ptox.hatp,w=nj[tried])

		
		lossvecc=ifelse(pipostc > targetc, (1-gam)*(pipostc-targetc), gam*(targetc-pipostc))  
		minlossc <- min(lossvecc)
		haveminc <- lossvecc==minlossc	
		Tc=lossvecc==min(lossvecc)
			possc=which(Tc)
			if(sum(Tc)==1){
				sugglevc=possc
			} else {
				if(all(pipostc[possc]>targetc)){
					sugglevc=min(possc)
				} else {
					sugglevc=max(possc)
				}
			}
			
		lossvecp=ifelse(pipostp > targetp, (1-gam)*(pipostp-targetp), gam*(targetp-pipostp))  
		minlossp <- min(lossvecp)
		haveminp <- lossvecp==minlossp	
		Tp=lossvecp==min(lossvecp)
			possp=which(Tp)
			if(sum(Tp)==1){
				sugglevp=possp
			} else {
				if(all(pipostp[possp]>targetp)){
					sugglevp=min(possp)
				} else {
					sugglevp=max(possp)
				}
			}

   			mtd=ifelse(stopc==3|stopp==3,0,min(sugglevc,sugglevp))
			dose.select[mtd]=dose.select[mtd]+1;    
     return(list(dose.select=dose.select,tox.datac=ycvec,tox.datap=ypvec,pt.allocation=nj,duration=max(arrival)+max(obswinc,obswinp)))
}

###Load the function 'titecdp.sim' 
titecdp.sim<-function(ntrial, PIc, PIp, targetc, targetp, gam, n, obswinc, obswinp, rate, 
	accrual, minwin, cl){
	ndose=length(PIc)
	
	d.select<-yoc<-yop<-no<-matrix(nrow=ntrial,ncol=ndose)
	td<-mtdp<-pcs<-rep(0,ntrial)
	
	for(i in 1:ntrial){
		result<-onetitecdp(PIc, PIp, targetc, targetp, gam, n, obswinc, obswinp, rate, accrual, minwin, cl) 
		d.select[i,]=result$dose.select
		yoc[i,]=result$tox.datac
		yop[i,]=result$tox.datap
		no[i,]=result$pt.allocation
		#pcs[i]=result$dose.select[which.min(abs(PI-target))]
		#mtdp[i]=result$pt.allocation[which.min(abs(PI-target))]
		td[i]=result$duration
	}
 cat("Simulation results for the TITE-PRO-CDP design (Wages et al., 2022+)\n");
 cat("targeting a NCI-CTCAE DLT rate of", targetc,"\n");
 cat("and targeting a PRO-CTCAE DLT rate of", targetp,"\n\n");
	cat("True NCI-CTCAE DLT probability:\n");
	cat(round(PIc,3),  sep="\t",  "\n")
	cat("True PRO-CTACE DLT probability:\n");
	cat(round(PIp,3),  sep="\t",  "\n")
	cat("MTD Selection percentage:\n")
	cat(formatC(colMeans(d.select)*100, digits=1, format="f"), sep="\t",  "\n")
	cat("Average number of C-DLTs:\n")
	cat(formatC(colMeans(yoc), digits=1, format="f"), sep="\t",   "\n")
	cat("Average number of P-DLTs:\n")
	cat(formatC(colMeans(yop), digits=1, format="f"), sep="\t",   "\n")
	cat("Average number of patients:\n")
	cat(formatC(colMeans(no), digits=1, format="f"), sep="\t",   "\n")
	cat("Average trial duration:\n")
	cat(round(mean(td),2), sep="\t",   "\n")
	cat("Percent stopped for safety:\n");
      cat(formatC(100-sum(colMeans(d.select)*100), digits=1, format="f"), sep="\t",   "\n");
		}
##########'titecdp.sim' end here
c1<-c(0.05,0.05,0.25,0.40,0.55)
p1<-c(0.17,0.18,0.35,0.50,0.65)

c2<-c(0.05,0.25,0.40,0.55,0.70)
p2<-c(0.10,0.15,0.35,0.50,0.65)

c3<-c(0.01,0.02,0.05,0.10,0.25)
p3<-c(0.04,0.09,0.17,0.20,0.35)

c4<-c(0.02,0.05,0.10,0.25,0.40)
p4<-c(0.09,0.17,0.20,0.35,0.50)

c5<-c(0.05,0.10,0.16,0.25,0.40)
p5<-c(0.05,0.20,0.35,0.50,0.65)

c6<-c(0.05,0.18,0.20,0.25,0.40)
p6<-c(0.17,0.35,0.50,0.65,0.80)

c7<-c(0.01,0.05,0.10,0.16,0.25)
p7<-c(0.04,0.05,0.20,0.35,0.50)

##new randomly generated scenarios
c8<-c(0.11,0.12,0.35,0.39,0.71)
p8<-c(0.08,0.25,0.57,0.60,0.87)

c9<-c(0.03,0.11,0.15,0.54,0.82)
p9<-c(0.003,0.01,0.09,0.29,0.32)

c10<-c(0.03,0.06,0.17,0.19,0.24)
p10<-c(0.07,0.08,0.09,0.10,0.13)

c11<-c(0.03,0.05,0.05,0.10,0.13)
p11<-c(0.08,0.15,0.19,0.25,0.36)

c12<-c(0.10,0.13,0.14,0.34,0.77)
p12<-c(0.14,0.14,0.16,0.30,0.34)

c13<-c(0.08,0.11,0.17,0.17,0.44)
p13<-c(0.01,0.01,0.05,0.24,0.35)

c14<-c(0.05,0.12,0.13,0.63,0.77)
p14<-c(0.004,0.01,0.06,0.07,0.27)

c15<-c(0.15,0.23,0.50,0.54,0.70)
p15<-c(0.07,0.19,0.79,0.85,0.87)

c16<-c(0.12,0.12,0.34,0.51,0.71)
p16<-c(0.02,0.40,0.63,0.68,0.99)


targetc<-0.25      ##target c toxicity rate 
targetp<-0.35      ##target p toxicity rate 
n <- 18
obswinc <- obswinp <- 3
accrual <- "fixed"
rate <- 6
minwin <- 1
ntrial <- 10000
gam <- 0.5
prob <- 0.5
cl <- 0.99999


PIc<-c5
PIp<-p5
set.seed(6508580)
titecdp.sim(ntrial, PIc, PIp, targetc, targetp, gam, n, obswinc, obswinp, rate, 
	accrual, minwin, cl)

