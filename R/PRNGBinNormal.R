

#' Pseudo-Random Number Generation of Correlated Binary (Bernoulli-Normal) or Count (Binomial-Normal) Outcomes
#'
#' This function produces pseudo-random numbers for responses for correlated (longitudinal) data.  The responses are generated as either binary (Bernoulli) data or count (Binomial), using a normally distributed random effect to induce autocorrelation.  The predictors must be provided as a vector.  The function returns a response vector, labeled "Outcomes".  
#' @param M The number of groups in the resulting data.  
#' @param mvec A vector indicating the length of each group in the data.  
#' @param mu The mean of the normal random effect.  
#' @param s The standard deviation of the normal random effect.  
#' @param D A vector of denominators for each response; a vector of 1's indicates Bernoulli data.  
#' @param beta0 The true intercept in the systematic component of the model.  
#' @param beta1 The true slope in the systematic component of the model.  
#' @param x A vector of predictor values.  
#' @param seed The seed for data generation.  
#' @keywords Data Generation Simulation Correlated
#' @export 
#' @examples
#' PRNGBinNormal()



PRNGBinNormal=function(M,mvec,mu,s,D,beta0,beta1,x,seed){

set.seed(seed)
seeds = rpois(M,125)

y = rep(0,length(D))

for(j in 1:M)
{
	set.seed(seeds[j])
	vj = rnorm(1,mu,s)
	for(k in 1:mvec[j])
	{
		set.seed(k*seeds[j])
		index = sum(mvec[0:(j-1)])+k
		pk = (exp(beta0+beta1*x[index]+vj))/(1+exp(beta0+beta1*x[index]+vj))
		y[index] = rbinom(1,D[index],pk)
	}
}

list(Outcomes=y)

} # END PRNGBinNormal

