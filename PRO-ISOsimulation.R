
library(Iso)

###Load the function 'incyte' 
procdp<-function(pc,pp,targetc,targetp,ac,bc,ap,bp,gam,cohortsize,ncohort,n.stop,start,cs){

### run a trial 	
    ndose = length(pc);   #number of combos
    yc=yp=n=numeric(ndose);  #number of toxicity/responses at each dose level
    n=numeric(ndose);  #number of treated patients at each dose level
    dose.select=codose.select=rep(0,ndose); # a vector of indicators for dose selection
    stop=stopc=stopp=0; #indicate if trial stops early
    currc = currp = start
    curr = min(currc,currp)
    unsafe=0
    i=1
while(i <= ncohort){
		
		if(any(n>n.stop)){
			stop<-0
			break
		}

		yc[curr] = yc[curr] + rbinom(1,cohortsize,pc[curr]);
		yp[curr] = yp[curr] + rbinom(1,cohortsize,pp[curr]);
		n[curr] = n[curr] + cohortsize;
		tried=which(n>0)
		

		if(1 - pbeta(targetc, yc[1] + ac, n[1] - yc[1] + bc) > cs){
			stopc=1
			break
		}

		if(1 - pbeta(targetp, yp[1] + ap, n[1] - yp[1] + bp) > cs){
			stopp=1
			break
		}


		if(length(tried)<ndose & all(yc==0)){
			currc<-sugglevc<-curr+1
		} else {
			uc=(yc[tried]+ac-(1/3))/(n[tried]+ac+bc-(2/3))
#qbeta(0.5,yc[tried]+ac,n[tried]-yc[tried]+bc)#(yc[tried]+ac)/(n[tried]+ac+bc)
			pipostc=pava(uc,w=n[tried])
			lossvecc=ifelse(pipostc > targetc, (1-gam)*(pipostc-targetc),
				gam*(targetc-pipostc))  
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
				currc=ifelse(n[sugglevc+1]==0,sugglevc+1,sugglevc)					
			} else {
				currc=sugglevc
			}
		}
		if(length(tried)<ndose & all(yp==0)){
			currp<-sugglevp<-curr+1
		} else {
			up=(yp[tried]+ap-(1/3))/(n[tried]+ap+bp-(2/3))
#qbeta(0.5,yp[tried]+ap,n[tried]-yp[tried]+bp)#(yp[tried]+ap)/(n[tried]+ap+bp)
			pipostp=pava(up,w=n[tried])
			lossvecp=ifelse(pipostp > targetp, (1-gam)*(pipostp-targetp),
				gam*(targetp-pipostp))  
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
				currp=ifelse(n[sugglevp+1]==0,sugglevp+1,sugglevp)					
			} else {
				currp=sugglevp
			}
		} 
		
		curr <- min(currc,currp)#sugglev <- min(currc,currp)
		#curr <- ifelse(sugglev==1|sum(n)<=(ncohort/2),sugglev,sample(c(sugglev,sugglev-1),1,prob=c(0.5,0.5)))
	
		i<-i+1
	}		
	if(stopc==0 & stopp==0){
			mtdc=sugglevc
			mtdp=sugglevp
			mtd=min(mtdc,mtdp)
			comtd=ifelse(mtd==1,0,mtd-1)
			dose.select[mtd]=dose.select[mtd]+1;
			codose.select[comtd]=codose.select[comtd]+1;		
		}
	return(list(dose.select=dose.select,codose.select=codose.select,tox.datac=yc,tox.datap=yp,
		pt.allocation=n,stopc=stopc,stopp=stopp))
}
##########'newcdp' end here


###Load the function 'bpocrm.sim' 
procdp.sim<-function(pc,pp,targetc,targetp,ac,
	bc,ap,bp,gam,cohortsize,ncohort,n.stop,start,cs,ntrial){


	ndose=length(pc)
	
	dose.select<-codose.select<-yc<-yp<-n<-matrix(nrow=ntrial,ncol=ndose)
	nstopc=nstopp=0
	
	for(i in 1:ntrial){
		result<-procdp(pc,pp,targetc,targetp,ac,bc,ap,bp,gam,cohortsize,ncohort,n.stop,start,cs)
		dose.select[i,]=result$dose.select
		codose.select[i,]=result$codose.select
		yc[i,]=result$tox.datac
		yp[i,]=result$tox.datap
		n[i,]=result$pt.allocation
		nstopc=nstopc+result$stopc
		nstopp=nstopp+result$stopp
	}
	cat("Simulation results for PRO-CDP design (Wages and Lin, 2023+)\n");
      cat("targeting a NCI-CTCAE DLT rate of", targetc,"\n");
      cat("and targeting a PRO-CTCAE DLT rate of", targetp,"\n\n");

	cat("True NCI-CTCAE DLT probability:\n");
	cat(round(pc,3), sep="\t",  "\n");
	cat("True PRO-CTACE DLT probability:\n");
      cat(round(pp,3), sep="\t",  "\n");
	cat("MTD selection percentage:\n");
	cat(formatC(colMeans(dose.select)*100, digits=1, format="f"), sep="\t",  "\n");
	cat("Average number of NCI-CTCAE DLTs:\n");
      cat(formatC(colMeans(yc), digits=1, format="f"), sep="\t",   "\n");
	cat("Average number of PRO-CTCAE DLTs:\n");
	cat(formatC(colMeans(yp), digits=1, format="f"), sep="\t",   "\n");
	cat("Average number of patients treated:\n");
	cat(formatC(colMeans(n), digits=1, format="f"), sep="\t",   "\n");
	cat("Percentage of trials stopped for NCI-CTCAE safety:\n");
	cat(nstopc/ntrial*100, "\n");
	cat("Percentage of trials stopped for PRO-CTCAE safety:\n");
	cat(nstopp/ntrial*100, "\n");
}

##########'procdp.sim' end here



start=1  		##starting dose
targetc=0.25      ##target c toxicity rate 
targetp=0.35      ##target p toxicity rate 

##True DLT probability scenarios
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
##new randomly generated scenarios
c1<-c(0.1445, 0.1921, 0.2783, 0.4747, 0.5066)
p1<-c(0.0614, 0.1007, 0.1716, 0.2224, 0.2992)

c2<-c(0.0031, 0.0885, 0.1529, 0.2779, 0.7576)
p2<-c(0.0297, 0.1033, 0.2252, 0.2364, 0.3766)

c3<-c(0.1103, 0.1135, 0.1266, 0.1782, 0.1866)
p3<-c(0.0548, 0.0622, 0.2390, 0.7226, 0.7394)

c4<-c(0.0393, 0.0902, 0.1029, 0.1456, 0.6634)
p4<-c(0.0058, 0.1493, 0.7544, 0.8593, 0.9743)

c5<-c(0.0426, 0.1438, 0.2421, 0.3182, 0.7565)
p5<-c(0.1611, 0.5437, 0.6156, 0.6454, 0.8405)

c6<-c(0.1074, 0.1091, 0.1343, 0.2545, 0.8185)
p6<-c(0.0089, 0.0858, 0.0994, 0.3783, 0.3961)


cohortsize=1    	##cohort size for each inclusion
ncohort=18      	##number of cohorts
n.stop=2455       	##Number of patients needed on one combination to stop the trial
ntrial=10000      ##number of simulated trials 

###calculate prior C
x<-2*targetc
mu<-targetc
u<-0.95

f<-function(b){
	pbeta(x,mu*b/(1-mu),b)-u
}
bc<-uniroot(f,c(0.0001,100))$root
ac<-mu*bc/(1-mu)
round(c(ac,bc),2)

###calculate prior P
x<-1.5*targetp
mu<-targetp
u<-0.9

f<-function(b){
	pbeta(x,mu*b/(1-mu),b)-u
}
bp<-uniroot(f,c(0.0001,100))$root
ap<-mu*bp/(1-mu)
round(c(ap,bp),2)
gam<-0.5
cs<-0.999999

pc<-c6
pp<-p6
procdp.sim(pc,pp,targetc,targetp,ac,bc,ap,bp,gam,cohortsize,ncohort,n.stop,start,cs,ntrial)


#
# Phase I study evaluating adjuvant hypofractionated
# whole pelvis radiation therapy (WPRT) in
# endometrial cancer (NCT04458402)
#
####################################################
start=1  		##starting dose
targetc=0.2      ##target c toxicity rate 
targetp=0.55      ##target p toxicity rate 

##True DLT probability scenarios

c1<-c(0.0821, 0.4096)
p1<-c(0.2707, 0.3176)

c2<-c(0.0891, 0.4725)
p2<-c(0.1259, 0.1876)

c3<-c(0.0021, 0.0966)
p3<-c(0.1744, 0.6809)

c4<-c(0.0534, 0.4548)
p4<-c(0.3586, 0.4450)

c5<-c(0.0053, 0.1731)
p5<-c(0.1233, 0.7996)

c6<-c(0.0752, 0.1500)
p6<-c(0.4259, 0.6118)

c7<-c(0.1368, 0.4834)
p7<-c(0.1190, 0.6251)

cohortsize=1    	##cohort size for each inclusion
ncohort=15      	##number of cohorts
n.stop=2455       	##Number of patients needed on one combination to stop the trial
ntrial=10000      ##number of simulated trials 

###calculate prior C
x<-2*targetc
mu<-targetc
u<-0.95

f<-function(b){
	pbeta(x,mu*b/(1-mu),b)-u
}
bc<-uniroot(f,c(0.0001,100))$root
ac<-mu*bc/(1-mu)
round(c(ac,bc),2)

###calculate prior P
x<-1.5*targetp
mu<-targetp
u<-0.9

f<-function(b){
	pbeta(x,mu*b/(1-mu),b)-u
}
bp<-uniroot(f,c(0.0001,100))$root
ap<-mu*bp/(1-mu)
round(c(ap,bp),2)
gam<-0.5
cs<-0.999999

pc<-c7
pp<-p7
procdp.sim(pc,pp,targetc,targetp,ac,bc,ap,bp,gam,cohortsize,ncohort,n.stop,start,cs,ntrial)
