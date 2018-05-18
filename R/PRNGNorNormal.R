

#' Pseudo-Random Number Generation of Correlated Continuous (Normal-Normal) Outcomes
#'
#' This function produces pseudo-random numbers for responses for correlated (longitudinal) data.  The responses are generated as continuous (Normal), using a normal random effect to induce autocorrelation.  The predictors must be provided as a vector.  The function returns a response vector, labeled "Outcomes".  
#' @param M The number of groups in the resulting data.  
#' @param mvec A vector indicating the length of each group in the data.  
#' @param mu_r The mean of the normal random effect.  
#' @param s_r The standard deviation of the normal random effect.  
#' @param s The standard deviation of the random error associated with the response.  
#' @param beta0 The true intercept in the systematic component of the model.  
#' @param beta1 The true slope in the systematic component of the model.  
#' @param x A vector of predictor values.  
#' @param seed The seed for data generation.  
#' @keywords Data Generation Simulation Correlated
#' @export 
#' @examples
#' PRNGNorNormal()



PRNGNorNormal=function(M,mvec,mu_r,s_r,s,beta0,beta1,x,seed){

set.seed(seed)
seeds = rpois(M,125)

y = rep(0,length(x))

for(j in 1:M)
{
	set.seed(seeds[j])
	vj = rnorm(1,mu_r,s_r)
	for(k in 1:mvec[j])
	{
		set.seed(k*seeds[j])
		index = sum(mvec[0:(j-1)])+k
		y[index] = rnorm(1,(beta0+beta1*x[index]+vj),s)
	}
}

list(Outcomes=y)

} # END PRNGNorNormal

