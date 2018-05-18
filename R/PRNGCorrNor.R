

#' Pseudo-Random Number Generation of Correlated Normal Outcomes
#'
#' This function produces pseudo-random numbers for both predictors and responses for correlated (longitudinal) data.  The responses are generated as normal data.  The predictors can be either continuous or binary, and can include both time-independent and time-dependent covariates.  Time-dependent covariates can be controlled as either Type I, Type II, Type III, or Type IV time-dependent covariates.  The function returns a response vector, a design matrix for the time-dependent covariates (with intercept column), and a design matrix for the time-independent covariates.  
#' @param seed The seed used to control pseudo-random number generation.  
#' @param S The number of subjects or groups to be generated.  These are assumed to be collections of auto-correlated data.  
#' @param Tvec The vector indicating the number of times (for each subject) or the number of units (for each group) to be generated.  
#' @param rhoyy The auto-regressive weight, a correlation between responses from the same subject or group.  
#' @param rhoxy The weight representing the correlation between "previous" time-dependent covariate values on response values.  This is associated with Type II time-dependent covariates.  
#' @param rhoyx The weight representing the correlation between "previous" responses and time-dependent covariate values.  This is associated with Type IV time-dependent covariates.  
#' @param TDCTypes The vector indicating the type (0, 1, 2, 3, or 4) of each time-dependent covariate.  A type of "0" indicates a time-independent covariate.  This vector should be ordered by time-independent covariates first.  
#' @param dataTypes The vector indicating the type of each predictor.  Either 'b' for binary or 'c' for continuous.  This includes all time-independent covariates, all time-dependent covariates, in that order, but not the intercept.  
#' @param beta The vector of "true" parameter values according to which data will be randomly generated.  This includes coefficients for an intercept, all time-independent covariates, and all time-dependent covariates, in that order.  
#' @param pred A vector of parameters describing the distributions of predictors.  Each binary predictor should have an associated probability of success, while each continuous predictor should have an associated standard deviation.  This includes all time-dependent covariates, all time-independent covariates, but not the intercept.  Continuous predictors will be assumed to be "centered" (mean of 0) but not "standardized" (variance can be other than 1).  
#' @param sigma The error standard deviation associated with the normal outcomes generated.  
#' @keywords Data Generation Simulation Correlated
#' @export 
#' @examples
#' PRNGCorrNor()




PRNGCorrNor = function(seed,S,Tvec,rhoyy,rhoxy,rhoyx,TDCTypes,dataTypes,beta,pred,sigma){

set.seed(seed)
seeds = rnorm(S,0,50)

# NUMBER OF TDC'S, EXCLUDE INTERCEPT #
q = length(beta) - 1


# GENERATE DATA FOR EACH SUBJECT #
for(i in 1:S)
{
	T_i = Tvec[i]
	set.seed(seeds[i])
	seeds_i = rnorm(T_i,0,25)

	# INITIALIZE SYSTEMATIC COMPONENTS #
	mu = rep(0,T_i)
	y = rep(0,T_i)
	X = matrix(0,T_i,q+1)
	px = matrix(0,T_i,q)
	X[,1] = rep(1,nrow(X)) # INTERCEPT COLUMN #

	# GENERATE DATA FOR SUBJECT i FOR ALL T_i "TIME POINTS" #
	for(t in 1:T_i)
	{
		set.seed(seeds_i[t])
		
		# FIRST TIME POINT USES NO PREVIOUS INFO #
		if(t == 1)
		{
			# GENERATE X-VALUES FOR TIME 1 #
			for(j in 1:q)
			{
				if(dataTypes[j] == 'c')
				{
					X[t,(j+1)] = rnorm(1,0,pred[j])
				}
				else # DATA TYPE == 'b' #
				{
					X[t,(j+1)] = rbinom(1,1,pred[j])
				}
			}

			# MEANS (mu) AND OUTCOMES (y) #
			mu[t]<-X[t,]%*%beta
			y[t]<-rnorm(1,mu[t],sigma)
			
		} # END t=1 X-VALUES AND Y-VALUES #

		# TIME POINTS BEYOND 1 MAY USE PREVIOUS INFO #
		else # t > 1 #
		{
			# GENERATE X-VALUES FOR TIME > 1 #
			for(j in 1:q)
			{
				if(dataTypes[j] == 'c')
				{
					# TIME-INDEPENDENT COVARIATE: SAME AS BEFORE #
					if(TDCTypes[j] == 0)
					{
						X[t,(j+1)] = X[(t-1),(j+1)]						
					}

					# NO PRIOR EFFECT FROM TYPES I OR II #
					if(TDCTypes[j] == 1 | TDCTypes[j] == 2)
					{
						X[t,(j+1)] = rnorm(1,pred[j])
					}

					# PRIOR EFFECT FROM RESPONSES FOR TYPES III OR IV #
					else # TYPE III OR TYPE IV #
					{
						px[t,j] = (0+rhoyx*log((pnorm(y[t-1],mu[t-1],sigma))/(1-(pnorm(y[t-1],mu[t-1],sigma)))))
						X[t,(j+1)] = rnorm(1,px[t,j],pred[j])
					}
				}

				# BINARY PREDICTORS #
				else # DATA TYPE == 'b' #
				{
					# TIME-INDEPENDENT COVARIATE: SAME AS BEFORE #
					if(TDCTypes[j] == 0)
					{
						X[t,(j+1)] = X[(t-1),(j+1)]						
					}

					# NO PRIOR EFFECT FROM TYPES I OR II #
					if(TDCTypes[j] == 1 | TDCTypes[j] == 2)
					{
						X[t,(j+1)] = rbinom(1,1,pred[j])
					}

					# PRIOR EFFECT FROM RESPONSES FOR TYPES III OR IV #
					else # TYPE III OR TYPE IV #
					{
						eta_j = log(pred[j]/(1-pred[j]))+rhoyx*log((pnorm(y[t-1],mu[t-1],sigma))/(1-(pnorm(y[t-1],mu[t-1],sigma))))
						px[t,j] = exp(eta_j)/(1+exp(eta_j))
						X[t,(j+1)] = rbinom(1,1,px[t,j])
					}
				}
			} # END t>1 X-VALUES #

			## MEANS (mu) AND OUTCOMES (y) ##
			
			# FIRST IDENTIFY TYPES II AND III TO INCLUDE EFFECTS FROM PRIOR PREDICTORS #			
			types23 = ifelse((TDCTypes>1 & TDCTypes<4),1,0)
			Xnoint = as.matrix(X[,-1]) # ELIMINATE INTERCEPT TO MATCH DIMENSIONS OF px, dataTypes #
			values = ifelse((dataTypes=='c'),log((pnorm(Xnoint[t-1,],px[t-1,],pred))/(1-(pnorm(Xnoint[t-1,],px[t-1,],pred)))),log((px[t-1,]^(Xnoint[t-1,])*(1-px[t-1,])^(1-Xnoint[t-1,]))/(1-(px[t-1,]^(Xnoint[t-1,]))*(1-px[t-1,])^(1-Xnoint[t-1,]))))

			# x_prev WILL INCLUDE PRIOR EFFECTS FOR TYPES II AND III, ZEROS FOR ALL OTHERS #
			x_prev = t(types23)%*%values

			eta_t = X[t,]%*%beta+rhoxy*(x_prev)+rhoyy*log((pnorm(y[t-1],mu[t-1],sigma))/(1-(pnorm(y[t-1],mu[t-1],sigma))))
			mu[t] = eta_t

			y[t] = rnorm(1,mu[t],sigma)
			
		} # END t>1 X-VALUES AND Y-VALUES #

	} # END T_i-LOOP #

	if(i == 1){
		Y = y
		Xmat = X
	}
	else{
		Y = c(Y,y)
		Xmat = rbind(Xmat,X)
	}

} # END S-LOOP #


# SEPARATE TIC'S FROM TDC'S, ACCOUNT FOR INTERCEPT COLUMN TO BE INCLUDED IN TDC DESIGN #
X = Xmat[,c(1,(which(TDCTypes != 0) + 1))]
Z = Xmat[,(which(TDCTypes == 0) + 1)]


list(yvec=Y,Xmat=X,Zmat=Z)

} # END PRNGCorrNor FUNCTION #





