

#' Pseudo-Random Number Generation of Correlated Binary (Bernoulli) Outcomes
#'
#' This function produces pseudo-random numbers for both predictors and responses for correlated (longitudinal) data.  The responses are generated as binary (Bernoulli) data.  The predictors can be either continuous or binary, and can include both time-independent and time-dependent covariates.  Time-dependent covariates can be controlled as either Type I, Type II, Type III, or Type IV time-dependent covariates.  The function returns a response vector, a design matric for the time-dependent covariates (with intercept column), and a design matrix for the time-independent covariates.  
#' @param seed The seed used to control pseudo-random number generation.  
#' @param S The number of subjects or groups to be generated.  These are assumed to be collections of auto-correlated data.  
#' @param Tvec The vector indicating the number of times (for each subject) or the number of units (for each group) to be generated.  
#' @param rhoyy The auto-regressive weight, a correlation between responses from the same subject or group.  
#' @param rhoxy The weight representing the correlation between "previous" time-dependent covariate values on response values.  This is associated with Type II time-dependent covariates.  
#' @param rhoyx The weight representing the correlation between "previous" responses and time-dependent covariate values.  This is associated with Type IV time-dependent covariates.  
#' @param TDCTypes The vector indicating the type (I, II, III, or IV) of each time-dependent covariate.  
#' @param dataTypes The vector indicating the type of each predictor.  Either 'b' for binary or 'c' for continuous.  
#' @param beta The vector of "true" parameter values according to which data will be randomly generated.  This includes coefficients for all time-independent covariates, all time-dependent covariates, and an intercept.  
#' @param pred A vector of parameters describing the distributions of predictors.  Each binary predictor should have an associated probability of success, while each continuous predictors should have an associated standard deviation.  
#' @keywords Data Generation Simulation Correlated
#' @export 
#' @examples
#' PRNGCorrBer()



PRNGCorrBer = function(seed,S,Tvec,rhoyy,rhoxy,rhoyx,TDCTypes,dataTypes,beta,pred){

set.seed(seed)
seeds = rnorm(S,0,50)

# NUMBER OF TDC'S, EXCLUDE INTERCEPT #
q = length(beta) - 1

beta_x = 3
X_TIC = matrix(0,sum(Tvec),1)

for(i in 1:S)
{
	T_i = Tvec[i]
	set.seed(seeds[i])
	seeds_i = rnorm(T_i,0,25)

	# TIME-INDEPENDENT COVARIATE #
	TIC_i = rnorm(1,5,3)
	X_TIC[(sum(Tvec[1:i-1])+1):sum(Tvec[1:i]),] = TIC_i

	#initialize y, X and mu
	mu<-rep(0,T_i)
	y<-rep(0,T_i)
	X = matrix(0,T_i,q+1)
	px = matrix(0,T_i,q)
	X[,1] = rep(1,nrow(X)) # INTERCEPT COLUMN #

	#loop running through T_i "time points"
	for(t in 1:T_i)
	{
		set.seed(seeds_i[t])
		
		#first time point uses no previous info
		if(t == 1)
		{
			#randomly generate x-values for time 1
			for(j in 1:q)
			{
				if(dataTypes[j] == 'c')
				{
					px[t,j] = beta_x
					X[t,(j+1)] = rnorm(1,px[t,j],pred[j])
				}
				else # DATA TYPE == 'b' #
				{
					px[t,j] = pred[j]
					X[t,(j+1)] = rbinom(1,1,px[t,j])
				}
			}

			#store probs of success (mu) and outcomes (y)
			mu[t]<-exp(X[t,]%*%beta+0.2*X_TIC[i,1])/(1+exp(X[t,]%*%beta+0.2*X_TIC[i,1]))
			y[t]<-rbinom(1,1,mu[t])
		}

		#2+ time points uses previous info
		else # t > 1 #
		{
			#randomly generate x-values for time > 1
			for(j in 1:q)
			{
				if(dataTypes[j] == 'c')
				{
					if(TDCTypes[j] == 1 | TDCTypes[j] == 2)
					{
						px[t,j] = beta_x
						X[t,(j+1)] = rnorm(1,px[t,j],pred[j])
					}
					else # TYPE III OR TYPE IV #
					{
						px[t,j] = (beta_x+rhoyx*log(((mu[t-1]^(y[t-1]))*(1-mu[t-1])^(1-y[t-1]))/(1-(mu[t-1]^(y[t-1]))*(1-mu[t-1])^(1-y[t-1]))))
						X[t,(j+1)] = rnorm(1,px[t,j],pred[j])
					}
				}
				else # DATA TYPE == 'b' #
				{
					if(TDCTypes[j] == 1 | TDCTypes[j] == 2)
					{
						px[t,j] = pred[j]
						X[t,(j+1)] = rbinom(1,1,px[t,j])
					}
					else # TYPE III OR TYPE IV #
					{
						eta_j = log(pred[j]/(1-pred[j]))+rhoyx*log(((mu[t-1]^(y[t-1]))*(1-mu[t-1])^(1-y[t-1]))/(1-(mu[t-1]^(y[t-1]))*(1-mu[t-1])^(1-y[t-1])))
						px[t,j] = exp(eta_j)/(1+exp(eta_j))
						X[t,(j+1)] = rbinom(1,1,px[t,j])
					}
				}
			}

			#store probs of success (mu) and outcomes (y)
			types23 = ifelse((TDCTypes>1 & TDCTypes<4),1,0)
			Xnoint = as.matrix(X[,-1]) # ELIMINATE INTERCEPT TO MATCH DIMENSIONS OF px, dataTypes #
			values = ifelse((dataTypes=='c'),log((pnorm(Xnoint[t-1,],px[t-1,],pred))/(1-(pnorm(Xnoint[t-1,],px[t-1,],pred)))),log((px[t-1,]^(Xnoint[t-1,])*(1-px[t-1,])^(1-Xnoint[t-1,]))/(1-(px[t-1,]^(Xnoint[t-1,]))*(1-px[t-1,])^(1-Xnoint[t-1,]))))
			x_prev = t(types23)%*%values

			eta_t = X[t,]%*%beta+rhoxy*(x_prev)+rhoyy*log((mu[t-1]^(y[t-1])*(1-mu[t-1])^(1-y[t-1]))/(1-(mu[t-1]^(y[t-1]))*(1-mu[t-1])^(1-y[t-1])))
			mu[t]<-exp(eta_t+0.2*X_TIC[i,1])/(1+exp(eta_t+0.2*X_TIC[i,1]))

			y[t]<-rbinom(1,1,mu[t])
		}

	} # END T_i-LOOP #

	if(i == 1){
		Y<-y
		Xmat<-X
	}
	else{
		Y<-c(Y,y)
		Xmat<-rbind(Xmat,X)
	}

} # END S-LOOP #


list(yvec=Y,Xmat=Xmat,Zmat=X_TIC)

} # END PRNGCorrBer FUNCTION #





