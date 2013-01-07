# Mike Tiernay
# 3/22/2011
# This code shows how to use WinBugs to estimate bivariate probit and ordered probit
# models.  The results from the models can be compared to running the same models
# on the dataset in stata.  Also, the code shows how to recover parameter estimates
# of a bivariate probit model that are similar to estimates obtained in stata when the 
# cutpoints of the model are manually inserted.
# Finally, code that 'should' implement a bivariate probit model fails to produce
# reliable estimates of the parameters and their standard errors.

# Note that to implement this code you need to have the full version of Winbugs
# For insructions see: http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/WinBUGS14_immortality_key.txt

# input: tiernay_bugs.dta

# Set the working directory
setwd("C:/Users/mike/Documents/bugs")


# Load the relevant packages
library(arm)
library(foreign)
library(R2WinBUGS)




# Load the data
bugs <- read.dta("tiernay_bugs.dta")
attach(bugs)

# Values of the number of observations in the dataset, to be used later
N <- nrow(bugs)
N2<- 2*N







###############################################################################

# Estimate an ordered probit model

###############################################################################


# Create a data matrix
# These are the only variables we need from the dataset to run the ordered probit
data.b <- list("N", "y1o", "x1")

                         
# Load the initial values
initials <- function(){
  list(beta1=rnorm(1,0,.5), beta2=rnorm(1,0,.5),cut1=runif(1,0,1))
}




# This is the likelihood function for the ordered probit    
oprobit_model <- function(){
	for(i in 1:N){	 #Loop over all the observations		

		mu[i] <- beta1*x1[i]+beta2 # mu here is our X*Beta
		   z[i] ~ dnorm(mu[i], 1)  # Get the probabily values associated with X*Beta
			# Define y1o to be categorical with probability value p[i,]
			#p[i,1:3] is a vector with 0's and one 1
      		#corresponding to whether draw of z[i,j] falls in partition
      		#of standard normal				
			y1o[i] ~ dcat(p[i,])
	
		# The step() function returns 1 if the argument is greater than zero, 0 otherwise
		p[i,1] <- step(-z[i])
		p[i,2] <- step(z[i])*step(gam1-z[i])
		p[i,3] <- step(z[i]-gam1)	
	}
	
# prior for  unknown gamma (cutpoint)
# there is a constant in the model so the first cut point is zero
# parameterized this way to ensure it is greater than zero
	gam1 <- exp(cut1)
	cut1 ~ dnorm(1,1)


#prior for the parameters to be estimated
	beta1 ~ dnorm(0,.01)
	beta2 ~ dnorm(0,.01)
}

# End of the likelihood function



# For some reason we need to do the following step so Bugs can read the model
## write model file:
filename <- file.path(tempdir(), "oprobit_model.bug")
write.model(oprobit_model, filename)

## and let's take a look:
file.show(filename)



# This is the call to Bugs.  Bugs will appear, and one must close the Bugs window
# in order for R to receive results back from Bugs.

# The first line tells Bugs what data, initial values, and model to use                                                                    
results <- bugs(data.b, initials, "oprobit_model.bug",

# List the parameters we want to save
 parameters = c("beta1","beta2","cut1"), 
# Run 2 chains, with a specific burn in, thinning, and number of simulations to save 
 n.chains=2, n.burnin=5000, n.iter=20000,n.thin=1,
debug=TRUE, bugs.seed=12345,

# Tell R where bugs is located 
 bugs.directory = "C:/WinBUGS14")

#End call to Bugs





# Print the results
# Note that in rhat < 1.1, chains have mixed and we have convergence (Gelman and Hill, p. 352)
print(results, digits.summary = 3)

attach.bugs(results)
                                                                       
                                    























###############################################################################

# Estimate the bivariate probit model

###############################################################################

# Create a data matrix
data.bi <- list("N", "y1b", "x1", "y2b", "x2")

                         
#Load the initial values
initials <- function(){
  list(beta11=rnorm(1,0,.5), beta12=rnorm(1,0,.5), beta21=rnorm(1,0,.5),
	 beta22=rnorm(1,0,.5),rho=runif(1,0,1)    )
}








biprobit_model <- function(){

for(i in 1:N){

   #linear link (x1, x2 observed; beta[1:2] are parameters of interest)
   mu.beta[i,1] <- x1[i]*beta11+beta12
   mu.beta[i,2] <- x2[i]*beta21+beta22
	

   #drawing from bivariate standard normal
   z[i,1:2] ~ dmnorm(mu.beta[i,1:2], varmu[1:2,1:2])
        #looping through two covariates
	for(j in 1:2){
		
			
				#Note that dcat requires are dummy's to be coded 1,2 instead of 0,1
				y1b[i] ~ dcat(pis[i,1,])	
				y2b[i] ~ dcat(pis[i,2,])	
				
			
			
                #pis[i,j,1:2] is a vector with a 0 and a 1
                #corresponding to whether draw of z[i,j] falls in partition
                #of standard normal
		pis[i,j,1] <- step(-z[i,j])
		pis[i,j,2] <- step(z[i,j])


	}
   }

#transforming var-cov matrix to precision matrix, because we are bayesian
varmu2[1,1] <- 1
varmu2[2,2] <- 1
varmu2[1,2] <- rho
varmu2[2,1] <- rho
varmu[1:2,1:2] <- inverse(varmu2[,])

#prior for rho
rho ~ dunif(-1,1)



#priors for betas
beta11~ dnorm(0,0.01)
beta21 ~ dnorm(0,0.01)
beta12~ dnorm(0,0.01)
beta22 ~ dnorm(0,0.01)

}


## write model file:
filename <- file.path(tempdir(), "biprobit_model.bug")
write.model(biprobit_model, filename)
## and let's take a look:
file.show(filename)



#send to bugs                                                                       
results <- bugs(data.bi, initials, "biprobit_model.bug",
 parameters = c("beta11","beta12","beta21","beta22","rho"),
 n.chains=2, n.burnin=1000, n.iter=4000,n.thin=1,
debug=TRUE, bugs.seed=12345,
 bugs.directory = "C:/WinBUGS14")


#print the results
print(results, digits.summary = 3)










































###############################################################################

# Estimate the bivariate ordered probit model with cutpoints
# Note we are trying to recover similar parameter estimates.  
# The standard errors seem to be underestimated, probably due to
# the fact that we are giving cutpoints intsead of estimating them.

###############################################################################



# Create a data matrix
data.bio <- list("N", "y1o", "x1", "y2o", "x2")

                         
#Load the initial values
initials <- function(){
  list(beta11=rnorm(1,0,.5), beta21=rnorm(1,0,.5),
	 rho=runif(1,0,1)   )
}








bioprobit_model <- function(){

for(i in 1:N){

   #linear link (x1, x2 observed; beta[1:2] are parameters of interest)
   mu.beta[i,1] <- x1[i]*beta11+beta12
   mu.beta[i,2] <- x2[i]*beta21+beta22
	

   #drawing from bivariate standard normal
   z[i,1:2] ~ dmnorm(mu.beta[i,1:2], varmu[1:2,1:2])
        #looping through two covariates
	for(j in 1:2){
		
			
				
				y1o[i] ~ dcat(pis[i,1,])	
				y2o[i] ~ dcat(pis[i,2,])	
				
			
			
                #pis[i,j,1:5] is a vector with 0's and one 1
                #corresponding to whether draw of z[i,j] falls in partition
                #of standard normal
		pis[i,j,1] <- step(-z[i,j])
		pis[i,j,2] <- step(z[i,j])*step(gam1[j]-z[i,j])
		pis[i,j,3] <- step(z[i,j]-gam1[j])

	}
   }

#transforming var-cov matrix to precision matrix
varmu2[1,1] <- 1
varmu2[2,2] <- 1
varmu2[1,2] <- rho
varmu2[2,1] <- rho
varmu[1:2,1:2] <- inverse(varmu2[,])

#prior for rho
rho ~ dunif(-1,1)



#priors for betas
beta11~ dnorm(0,0.01)
beta21 ~ dnorm(0,0.01)

#These cutpoints will be given and not estimated in this model
#beta12~ dnorm(0,0.01)
#beta22 ~ dnorm(0,0.01)

#Give the cutpoints, taken from stata estimation
beta12 <- -.197
beta22 <- -.604


g1[1] <- log(2.12 + beta12) 
gam1[1] <- exp(g1[1]) 

g1[2] <- log(2.64 + beta22) 
gam1[2] <- exp(g1[2]) 



}


## write model file:
filename <- file.path(tempdir(), "bioprobit_model.bug")
write.model(bioprobit_model, filename)
## and let's take a look:
file.show(filename)



#send to bugs                                                                       
results <- bugs(data.bio, initials, "bioprobit_model.bug",
 parameters = c("beta11","beta21","rho"),
 n.chains=2, n.burnin=1000, n.iter=4000,n.thin=1,
debug=TRUE, bugs.seed=12345,
 bugs.directory = "C:/WinBUGS14")


#print the results
print(results, digits.summary = 3)









































###############################################################################

# Estimate the bivariate ordered probit model
# This time estimating the cutpoints ourselves

###############################################################################



# Create a data matrix
data.bio <- list("N", "y1o", "x1", "y2o", "x2")

                         
#Load the initial values
initials <- function(){
  list(beta11=rnorm(1,0,1), beta21=rnorm(1,0,1),beta12=rnorm(1,0,1), 
		beta22=rnorm(1,0,1),rho=runif(1,0,1),g1=runif(2,0,1)   )
}
  



bioprobit_model <- function(){

for(i in 1:N){

   #linear link (x1, x2 observed; beta[1:2] are parameters of interest)
   mu.beta[i,1] <- x1[i]*beta11+beta12
   mu.beta[i,2] <- x2[i]*beta21+beta22
	

   #drawing from bivariate standard normal
   z[i,1:2] ~ dmnorm(mu.beta[i,1:2], varmu[1:2,1:2])
        #looping through two covariates
	for(j in 1:2){
		
			
				
				y1o[i] ~ dcat(pis[i,1,])	
				y2o[i] ~ dcat(pis[i,2,])	
				
			
			
                #pis[i,j,1:3] is a vector with 0's and one 1
                #corresponding to whether draw of z[i,j] falls in partition
                #of standard normal
		pis[i,j,1] <- step(-z[i,j])
		pis[i,j,2] <- step(z[i,j])*step(gam1[j]-z[i,j])
		pis[i,j,3] <- step(z[i,j]-gam1[j])

	}
   }

#transforming var-cov matrix to precision matrix
varmu2[1,1] <- 1
varmu2[2,2] <- 1
varmu2[1,2] <- rho
varmu2[2,1] <- rho
varmu[1:2,1:2] <- inverse(varmu2[,])

#prior for rho
rho ~ dunif(-1,1)



#priors for betas
beta11~ dnorm(0,0.01)
beta21 ~ dnorm(0,0.01)

beta12~ dnorm(0,0.01)
beta22 ~ dnorm(0,0.01)


#prior for  unknown gammas (cutpoints)
for(j in 1:2){
   gam1[j] <- exp(g1[j])
   g1[j] ~ dnorm(1,1)
}




}


## write model file:
filename <- file.path(tempdir(), "bioprobit_model.bug")
write.model(bioprobit_model, filename)
## and let's take a look:
file.show(filename)



#send to bugs                                                                       
results <- bugs(data.bio, initials, "bioprobit_model.bug",
 parameters = c("beta11","beta21","beta12","beta22","rho","g1"),
 n.chains=2, n.burnin=1000, n.iter=4000,n.thin=1,
debug=TRUE, bugs.seed=12345,
 bugs.directory = "C:/WinBUGS14")


#print the results
print(results, digits.summary = 3)






