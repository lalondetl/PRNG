

#' Pseudo-Random Number Generation of Overdispersed Correlated Binary (Bernoulli-Beta) or Count (Binomial-Beta) Outcomes
#'
#' This function produces pseudo-random numbers for responses for correlated (longitudinal) data.  The responses are generated as either binary (Bernoulli) data or count (Binomial), using a transformed Beta random effect to induce autocorrelation.  Once the Beta-distributed random effects are generated, the values are transformed according to a logit transformation to ensure they interact with the response similarly to otherpredictors.  These random transformed effects are then passed through the systematic component to a Binomial-Beta pseudo-random number generator from the VGAM package.  Therefore there are two instances of random effects: in the transformed Beta random effects, and through the Binomial-Beta pseudo-random number generation.  The predictors must be provided as a vector.  The function returns a response vector, labeled "Outcomes".  
#' @param M The number of groups in the resulting data.  
#' @param mvec A vector indicating the length of each group in the data.  
#' @param D A vector of denominators for each response; a vector of 1's indicates Bernoulli data.  
#' @param beta0 The true intercept in the systematic component of the model.  
#' @param beta1 The true slope in the systematic component of the model.  
#' @param x A vector of predictor values.  
#' @param seed The seed for data generation.  
#' @keywords Data Generation Simulation Correlated
#' @export 
#' @examples
#' PRNGODBinBeta()



PRNGODBinBeta=function(M,mvec,D,beta0,beta1,x,seed){

set.seed(seed)
seeds = rpois(M,125)

y = rep(0,length(x))

for(i in 1:M)
{
	set.seed(seeds[i])
	ui = rbeta(1,2,3)
	vi = log(ui/(1-ui))
	seeds2 = rpois(mvec[i],75)

	for(j in 1:mvec[i])
	{
		set.seed(seeds2[j])
		index = sum(mvec[0:(i-1)])+j
		alphaij = 4*exp(beta0+beta1*x[index]+vi)
		y[index] = rbetabin.ab(1,D[index],alphaij,4)
	}
}

list(Outcomes=y)

} # END PRNGODBinBeta



