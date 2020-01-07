library(MASS)
library(ROCR)
library(randomForest)
library(binom)

load("C:\\U\\GHIT\\NatMed_Paper_OrigBrazData\\Figure5_algorithm\\Algorithm\\AB_LOGIT.RData")

load("C:\\U\\GHIT\\NatMed_Paper_OrigBrazData\\Figure5_algorithm\\Algorithm\\AB_LDA.RData")

load("C:\\U\\GHIT\\NatMed_Paper_OrigBrazData\\Figure5_algorithm\\Algorithm\\AB_QDA.RData")

load("C:\\U\\GHIT\\NatMed_Paper_OrigBrazData\\Figure5_algorithm\\Algorithm\\AB_Dtree.RData")

load("C:\\U\\GHIT\\NatMed_Paper_OrigBrazData\\Figure5_algorithm\\Algorithm\\AB_Rforest.RData")

load("C:\\U\\GHIT\\NatMed_Paper_OrigBrazData\\Figure5_algorithm\\Algorithm\\AB_DYN1.RData")

load("C:\\U\\GHIT\\NatMed_Paper_OrigBrazData\\Figure5_algorithm\\Algorithm\\AB_DYN2.RData")


##################################################################
##################################################################
##                                                              ##
##  ##### #   ## #####  ####  ##### #     # #####  ##    #####  ## 
##  ##    ##  ## ##    ##     ##    ##   ## ##  ## ##    ##     ##
##  ####  ### ## ####   ####  ####  ####### #####  ##    ####   ##
##  ##    ## ### ##        ## ##    ## # ## ##  ## ##    ##     ##
##  ##### ##  ## #####  ####  ##### ##   ## #####  ##### #####  ##
##                                                              ##
################################################################## 
##################################################################


STEP_CHULL = function( XX_VEC, YY_VEC )
{
	x_step = 0	
	y_step = 0	

	x_chull = x_step
	y_chull = y_step

	while( (x_step < 1) && (y_step < 1) )
	{
		index = intersect( which(XX_VEC >= x_step), which(YY_VEC > y_step) )
		index = index[-intersect( which(XX_VEC == x_step), which(YY_VEC == y_step) )]	

		next_step = index[which.min(XX_VEC[index])]
	
		x_step = XX_VEC[next_step]
		y_step = YY_VEC[next_step]

		x_chull = c(x_chull, x_step)
		y_chull = c(y_chull, y_step)
	}

	x_chull = c(x_chull, 1)
	y_chull = c(y_chull, 1)

	cbind( x_chull, y_chull)
	
}



######################################
######################################
##                                  ##
##  PART 1                          ##
##  Thailand predicting Thailand    ## 
##                                  ##                                  
######################################
######################################


AB_ENSEMBLE_thai_thai = list()

for(i in 1:7)
{
	AB_ENSEMBLE_thai_thai[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_thai_thai[[i]]))
	{
		AB_ENSEMBLE_thai_thai[[i]][[count]] = AB_LOGIT_thai_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_thai_thai[[i]]))
	{
		AB_ENSEMBLE_thai_thai[[i]][[count]] = AB_LDA_thai_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_thai_thai[[i]]))
	{
		AB_ENSEMBLE_thai_thai[[i]][[count]] = AB_QDA_thai_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_thai_thai[[i]]))
	{
		AB_ENSEMBLE_thai_thai[[i]][[count]] = AB_Dtree_thai_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_thai_thai[[i]]))
	{
		AB_ENSEMBLE_thai_thai[[i]][[count]] = AB_Rforest_thai_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_thai_thai[[i]]))
	{
		AB_ENSEMBLE_thai_thai[[i]][[count]] = AB_DYN1_thai_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_thai_thai[[i]]))
	{
		AB_ENSEMBLE_thai_thai[[i]][[count]] = AB_DYN2_thai_thai[[i]][[j]]

		count = count + 1
	}	
}


AB_ENSEMBLE_thai_thai_chull = list()

AB_ENSEMBLE_thai_thai_chull[[1]] = AB_thai_1

XX_VEC = AB_thai_1[,1]
YY_VEC = AB_thai_1[,2]

for(i in 1:length(AB_ENSEMBLE_thai_thai))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_thai_thai[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_ENSEMBLE_thai_thai[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_thai_thai[[i]][[j]][,2] )
	}

	if( length(which(XX_VEC < 0)) )
	{
		XX_VEC[which(XX_VEC < 0)] = 0
	}
	if( length(which(XX_VEC > 1)) )
	{
		XX_VEC[which(XX_VEC > 1)] = 1
	}

	if( length(which(YY_VEC < 0)) )
	{
		YY_VEC[which(YY_VEC < 0)] = 0
	}
	if( length(which(YY_VEC > 1)) )
	{
		YY_VEC[which(YY_VEC > 1)] = 1
	}



	XY_hull_plot = STEP_CHULL( XX_VEC, YY_VEC )


	AB_ENSEMBLE_thai_thai_chull[[1+i]] = XY_hull_plot 
}


######################################
######################################
##                                  ##
##  PART 6                          ##
##  Brazil predicting Brazil        ## 
##                                  ##                                  
######################################
######################################

AB_ENSEMBLE_braz_braz= list()

for(i in 1:7)
{
	AB_ENSEMBLE_braz_braz[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_braz_braz[[i]]))
	{
		AB_ENSEMBLE_braz_braz[[i]][[count]] = AB_LOGIT_braz_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_braz_braz[[i]]))
	{
		AB_ENSEMBLE_braz_braz[[i]][[count]] = AB_LDA_braz_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_braz_braz[[i]]))
	{
		AB_ENSEMBLE_braz_braz[[i]][[count]] = AB_QDA_braz_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_braz_braz[[i]]))
	{
		AB_ENSEMBLE_braz_braz[[i]][[count]] = AB_Dtree_braz_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_braz_braz[[i]]))
	{
		AB_ENSEMBLE_braz_braz[[i]][[count]] = AB_Rforest_braz_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_braz_braz[[i]]))
	{
		AB_ENSEMBLE_braz_braz[[i]][[count]] = AB_DYN1_braz_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_braz_braz[[i]]))
	{
		AB_ENSEMBLE_braz_braz[[i]][[count]] = AB_DYN2_braz_braz[[i]][[j]]

		count = count + 1
	}	
}


AB_ENSEMBLE_braz_braz_chull = list()

AB_ENSEMBLE_braz_braz_chull[[1]] = AB_braz_1

XX_VEC = AB_braz_1[,1]
YY_VEC = AB_braz_1[,2]

for(i in 1:length(AB_ENSEMBLE_braz_braz))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_braz_braz[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_ENSEMBLE_braz_braz[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_braz_braz[[i]][[j]][,2] )
	}

	if( length(which(XX_VEC < 0)) )
	{
		XX_VEC[which(XX_VEC < 0)] = 0
	}
	if( length(which(XX_VEC > 1)) )
	{
		XX_VEC[which(XX_VEC > 1)] = 1
	}

	if( length(which(YY_VEC < 0)) )
	{
		YY_VEC[which(YY_VEC < 0)] = 0
	}
	if( length(which(YY_VEC > 1)) )
	{
		YY_VEC[which(YY_VEC > 1)] = 1
	}



	XY_hull_plot = STEP_CHULL( XX_VEC, YY_VEC )


	AB_ENSEMBLE_braz_braz_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 11                         ##
##  Solomons predicting Solomons    ## 
##                                  ##                                  
######################################
######################################

AB_ENSEMBLE_sol_sol = list()

for(i in 1:7)
{
	AB_ENSEMBLE_sol_sol[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_sol_sol[[i]]))
	{
		AB_ENSEMBLE_sol_sol[[i]][[count]] = AB_LOGIT_sol_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_sol_sol[[i]]))
	{
		AB_ENSEMBLE_sol_sol[[i]][[count]] = AB_LDA_sol_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_sol_sol[[i]]))
	{
		AB_ENSEMBLE_sol_sol[[i]][[count]] = AB_QDA_sol_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_sol_sol[[i]]))
	{
		AB_ENSEMBLE_sol_sol[[i]][[count]] = AB_Dtree_sol_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_sol_sol[[i]]))
	{
		AB_ENSEMBLE_sol_sol[[i]][[count]] = AB_Rforest_sol_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_sol_sol[[i]]))
	{
		AB_ENSEMBLE_sol_sol[[i]][[count]] = AB_DYN1_sol_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_sol_sol[[i]]))
	{
		AB_ENSEMBLE_sol_sol[[i]][[count]] = AB_DYN2_sol_sol[[i]][[j]]

		count = count + 1
	}	
}


AB_ENSEMBLE_sol_sol_chull = list()

AB_ENSEMBLE_sol_sol_chull[[1]] = AB_sol_1

XX_VEC = AB_sol_1[,1]
YY_VEC = AB_sol_1[,2]

for(i in 1:length(AB_ENSEMBLE_sol_sol))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_sol_sol[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_ENSEMBLE_sol_sol[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_sol_sol[[i]][[j]][,2] )
	}

	if( length(which(XX_VEC < 0)) )
	{
		XX_VEC[which(XX_VEC < 0)] = 0
	}
	if( length(which(XX_VEC > 1)) )
	{
		XX_VEC[which(XX_VEC > 1)] = 1
	}

	if( length(which(YY_VEC < 0)) )
	{
		YY_VEC[which(YY_VEC < 0)] = 0
	}
	if( length(which(YY_VEC > 1)) )
	{
		YY_VEC[which(YY_VEC > 1)] = 1
	}


	XY_hull_plot = STEP_CHULL( XX_VEC, YY_VEC )


	AB_ENSEMBLE_sol_sol_chull[[1+i]] = XY_hull_plot 
}




######################################
######################################
##                                  ##
##  PART 16                         ##
##  All predicting All              ## 
##                                  ##                                  
######################################
######################################


AB_ENSEMBLE_all_all= list()

for(i in 1:7)
{
	AB_ENSEMBLE_all_all[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_all_all[[i]]))
	{
		AB_ENSEMBLE_all_all[[i]][[count]] = AB_LOGIT_all_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_all_all[[i]]))
	{
		AB_ENSEMBLE_all_all[[i]][[count]] = AB_LDA_all_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_all_all[[i]]))
	{
		AB_ENSEMBLE_all_all[[i]][[count]] = AB_QDA_all_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_all_all[[i]]))
	{
		AB_ENSEMBLE_all_all[[i]][[count]] = AB_Dtree_all_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_all_all[[i]]))
	{
		AB_ENSEMBLE_all_all[[i]][[count]] = AB_Rforest_all_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_all_all[[i]]))
	{
		AB_ENSEMBLE_all_all[[i]][[count]] = AB_DYN1_all_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_all_all[[i]]))
	{
		AB_ENSEMBLE_all_all[[i]][[count]] = AB_DYN2_all_all[[i]][[j]]

		count = count + 1
	}	
}


AB_ENSEMBLE_all_all_chull = list()


AB_ENSEMBLE_all_all_chull[[1]] = AB_all_1

XX_VEC = AB_all_1[,1]
YY_VEC = AB_all_1[,2]

for(i in 1:length(AB_ENSEMBLE_all_all))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_all_all[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_ENSEMBLE_all_all[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_all_all[[i]][[j]][,2] )
	}

	if( length(which(XX_VEC < 0)) )
	{
		XX_VEC[which(XX_VEC < 0)] = 0
	}
	if( length(which(XX_VEC > 1)) )
	{
		XX_VEC[which(XX_VEC > 1)] = 1
	}

	if( length(which(YY_VEC < 0)) )
	{
		YY_VEC[which(YY_VEC < 0)] = 0
	}
	if( length(which(YY_VEC > 1)) )
	{
		YY_VEC[which(YY_VEC > 1)] = 1
	}

	XY_hull_plot = STEP_CHULL( XX_VEC, YY_VEC )


	AB_ENSEMBLE_all_all_chull[[1+i]] = XY_hull_plot 
}







############################
############################
##                        ##
##  ###### #####  #####   ## 
##    ##   ##  ## ##  ##  ##
##    ##   #####  #####   ##
##    ##   ##     ##      ##
##    ##   ##     ##      ##
##                        ##
############################
############################


spec = 1 - AB_ENSEMBLE_thai_thai_chull[[8]][,1]
sens = AB_ENSEMBLE_thai_thai_chull[[8]][,2]


Thai_SS = matrix(NA, nrow=3, ncol=2)
colnames(Thai_SS) = c("sensitivity", "specificity")

Thai_SS[1,] = 0.5*( sens[which.min(abs( sens - spec ))] + spec[which.min(abs( sens - spec ))] )

Thai_SS[2,1] = 0.5
Thai_SS[2,2] = spec[which.min(abs(sens - 0.5))] 

Thai_SS[3,1] = sens[which.min(abs(spec - 0.5))]
Thai_SS[3,2] = 0.5






spec = 1 - AB_ENSEMBLE_braz_braz_chull[[8]][,1]
sens = AB_ENSEMBLE_braz_braz_chull[[8]][,2]


Braz_SS = matrix(NA, nrow=3, ncol=2)
colnames(Braz_SS) = c("sensitivity", "specificity")

Braz_SS[1,] = 0.5*( sens[which.min(abs( sens - spec ))] + spec[which.min(abs( sens - spec ))] )


Braz_SS[2,1] = 0.5
Braz_SS[2,2] = spec[which.min(abs(sens - 0.5))] 


##Braz_SS[3,1] = sens[which.min(abs(spec - 0.5))]
##Braz_SS[3,2] = 0.5


Braz_SS[3,1] = 1 
Braz_SS[3,2] = max(spec[which(sens == 1)])


spec = 1 - AB_ENSEMBLE_sol_sol_chull[[8]][,1]
sens = AB_ENSEMBLE_sol_sol_chull[[8]][,2]


Sol_SS = matrix(NA, nrow=3, ncol=2)
colnames(Sol_SS) = c("sensitivity", "specificity")

Sol_SS[1,] = 0.5*( sens[which.min(abs( sens - spec ))] + spec[which.min(abs( sens - spec ))] )


Sol_SS[2,1] = 0.5
Sol_SS[2,2] = spec[which.min(abs(sens - 0.5))] 

Sol_SS[3,1] = sens[which.min(abs(spec - 0.5))]
Sol_SS[3,2] = 0.5









spec = 1 - AB_ENSEMBLE_all_all_chull[[8]][,1]
sens = AB_ENSEMBLE_all_all_chull[[8]][,2]



All_SS = matrix(NA, nrow=3, ncol=2)
colnames(All_SS) = c("sensitivity", "specificity")

All_SS[1,] = 0.5*( sens[which.min(abs( sens - spec ))] + spec[which.min(abs( sens - spec ))] )


All_SS[2,1] = 0.5
All_SS[2,2] = spec[which.min(abs(sens - 0.5))] 

All_SS[3,1] = sens[which.min(abs(spec - 0.5))]
All_SS[3,2] = 0.5


###################################
###################################
##                               ##
##  ####    ####  ######  ####   ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ## ######   ##   ######  ## 
##  ## ##  ##  ##   ##   ##  ##  ##
##  ####   ##  ##   ##   ##  ##  ##
##                               ##
###################################
###################################

par(mfrow=c(1,3))


################
################
##            ##
##  Thailand  ##
##            ##
################
################


#########################################
## Create a function for data conversion

thai_date_convert = function( xxx )
{
	num_date = NA

	if( nchar(as.vector(xxx))==9 )
	{
		if( substr( xxx, 2, 2 ) == "/" )
		{
			day   = as.numeric(substr(xxx, 1, 1))
			month = as.numeric(substr(xxx, 3, 4))
			year  = as.numeric(substr(xxx, 6, 9))
		}

		if( year==2013 ){ num_date <- 0 }
		if( year==2014 ){ num_date <- 1*365 }

		if( month==1  ){ num_date <- num_date + 0 }
		if( month==2  ){ num_date <- num_date + 31 }
		if( month==3  ){ num_date <- num_date + 59 }
		if( month==4  ){ num_date <- num_date + 90 }
		if( month==5  ){ num_date <- num_date + 120 }
		if( month==6  ){ num_date <- num_date + 151 }
		if( month==7  ){ num_date <- num_date + 181 }
		if( month==8  ){ num_date <- num_date + 212 }
		if( month==9  ){ num_date <- num_date + 243 }
		if( month==10  ){ num_date <- num_date + 273 }
		if( month==11 ){ num_date <- num_date + 304 }
		if( month==12 ){ num_date <- num_date + 334 }
	
		num_date <- num_date + day
	}

	if( nchar(as.vector(xxx))==10 )
	{
		if( substr( xxx, 3, 3 ) == "/" )
		{
			day   = as.numeric(substr(xxx, 1, 2))
			month = as.numeric(substr(xxx, 4, 5))
			year  = as.numeric(substr(xxx, 7, 10))
		}

		if( year==2013 ){ num_date <- 0 }
		if( year==2014 ){ num_date <- 1*365 }

		if( month==1  ){ num_date <- num_date + 0 }
		if( month==2  ){ num_date <- num_date + 31 }
		if( month==3  ){ num_date <- num_date + 59 }
		if( month==4  ){ num_date <- num_date + 90 }
		if( month==5  ){ num_date <- num_date + 120 }
		if( month==6  ){ num_date <- num_date + 151 }
		if( month==7  ){ num_date <- num_date + 181 }
		if( month==8  ){ num_date <- num_date + 212 }
		if( month==9  ){ num_date <- num_date + 243 }
		if( month==10  ){ num_date <- num_date + 273 }
		if( month==11 ){ num_date <- num_date + 304 }
		if( month==12 ){ num_date <- num_date + 334 }
	
		num_date <- num_date + day
	}
	
	num_date
} 



#########################################
## Read in raw data from Thailand

Thai_data_read <- read.csv( "C:\\U\\GHIT\\2_All_Data\\2_Data\\data\\raw\\Thailand\\thai_epi_data.csv" )

Thai_PvPCR <- as.data.frame( Thai_data_read[,c(53,6)] )

Thai_PvPCR <- cbind( Thai_PvPCR, rep(NA, nrow(Thai_PvPCR)) )
colnames(Thai_PvPCR)[3] <- "studyday"

for(i in 1:nrow(Thai_PvPCR))
{
	Thai_PvPCR[i,3] <- thai_date_convert( Thai_PvPCR[i,2] )
}

Thai_PvPCR[,3] <- Thai_PvPCR[,3] - min(Thai_PvPCR[,3])


Thai_bin_edges <- seq(from=-14, by=28, length=15) 
N_thai_bins <- length(Thai_bin_edges) - 1

Thai_bin_edges[N_thai_bins+1] <- max(Thai_PvPCR[,3])



hist( Thai_PvPCR[,3], breaks=200, col="grey",
main="Thai samples", xlab="time (days)" )

for(j in 1:length(Thai_bin_edges))
{
	points(x=rep(Thai_bin_edges[j],2), y=c(0,2000), type='l', col="red", lwd=2)	
}

Thai_PvPCR_bins <- matrix(NA, nrow=3, ncol=N_thai_bins)
for(j in 1:N_thai_bins)
{
	index_j <- which( (Thai_PvPCR[,3] > Thai_bin_edges[j]) & (Thai_PvPCR[,3] <= Thai_bin_edges[j+1]) ) 

	Thai_PvPCR_bins[,j] <- as.numeric(as.vector(
            	                	 binom.confint( sum(Thai_PvPCR[index_j,1]), length(Thai_PvPCR[index_j,1]), method="wilson")[1,4:6]
            	                ))
}

Thai_PvPCR_times <- 0.5*( Thai_bin_edges[-1] + Thai_bin_edges[-length(Thai_bin_edges)] ) 

Thai_PvPCR_times <- 2019.95 - max(Thai_PvPCR_times/365) + Thai_PvPCR_times/365



################
################
##            ##
##  Brazil    ##
##            ##
################
################


#########################################
## Create a function for data conversion

braz_date_convert = function( xxx )
{
	num_date = NA

	NN <- nchar(as.vector(xxx))

	slash <- c()

	for(k in 1:NN)
	{
		if( substr(xxx, k, k) == "/" )
		{	
			slash <- c( slash, k)
		}
	}

	if( (NN==8) && (slash[1]==2) && (slash[2]==4) )
	{
		day   = as.numeric(substr(xxx, 1, 1))
		month = as.numeric(substr(xxx, 3, 3))
		year  = as.numeric(substr(xxx, 5, 8))

		if( year==2013 ){ num_date <- 0 }
		if( year==2014 ){ num_date <- 1*365 }
		if( year==2015 ){ num_date <- 2*365 }

		if( month==1  ){ num_date <- num_date + 0 }
		if( month==2  ){ num_date <- num_date + 31 }
		if( month==3  ){ num_date <- num_date + 59 }
		if( month==4  ){ num_date <- num_date + 90 }
		if( month==5  ){ num_date <- num_date + 120 }
		if( month==6  ){ num_date <- num_date + 151 }
		if( month==7  ){ num_date <- num_date + 181 }
		if( month==8  ){ num_date <- num_date + 212 }
		if( month==9  ){ num_date <- num_date + 243 }
		if( month==10  ){ num_date <- num_date + 273 }
		if( month==11 ){ num_date <- num_date + 304 }
		if( month==12 ){ num_date <- num_date + 334 }
	
		num_date <- num_date + day
	}

	if( (NN==9) && (slash[1]==2) && (slash[2]==5) )
	{
		day   = as.numeric(substr(xxx, 1, 1))
		month = as.numeric(substr(xxx, 3, 4))
		year  = as.numeric(substr(xxx, 6, 9))

		if( year==2013 ){ num_date <- 0 }
		if( year==2014 ){ num_date <- 1*365 }
		if( year==2015 ){ num_date <- 2*365 }

		if( month==1  ){ num_date <- num_date + 0 }
		if( month==2  ){ num_date <- num_date + 31 }
		if( month==3  ){ num_date <- num_date + 59 }
		if( month==4  ){ num_date <- num_date + 90 }
		if( month==5  ){ num_date <- num_date + 120 }
		if( month==6  ){ num_date <- num_date + 151 }
		if( month==7  ){ num_date <- num_date + 181 }
		if( month==8  ){ num_date <- num_date + 212 }
		if( month==9  ){ num_date <- num_date + 243 }
		if( month==10  ){ num_date <- num_date + 273 }
		if( month==11 ){ num_date <- num_date + 304 }
		if( month==12 ){ num_date <- num_date + 334 }
	
		num_date <- num_date + day
	}

	if( (NN==9) && (slash[1]==3) && (slash[2]==5) )
	{
		day   = as.numeric(substr(xxx, 1, 2))
		month = as.numeric(substr(xxx, 4, 4))
		year  = as.numeric(substr(xxx, 6, 9))

		if( year==2013 ){ num_date <- 0 }
		if( year==2014 ){ num_date <- 1*365 }
		if( year==2015 ){ num_date <- 2*365 }

		if( month==1  ){ num_date <- num_date + 0 }
		if( month==2  ){ num_date <- num_date + 31 }
		if( month==3  ){ num_date <- num_date + 59 }
		if( month==4  ){ num_date <- num_date + 90 }
		if( month==5  ){ num_date <- num_date + 120 }
		if( month==6  ){ num_date <- num_date + 151 }
		if( month==7  ){ num_date <- num_date + 181 }
		if( month==8  ){ num_date <- num_date + 212 }
		if( month==9  ){ num_date <- num_date + 243 }
		if( month==10  ){ num_date <- num_date + 273 }
		if( month==11 ){ num_date <- num_date + 304 }
		if( month==12 ){ num_date <- num_date + 334 }
	
		num_date <- num_date + day
	}

	if( (NN==10) && (slash[1]==3) && (slash[2]==6) )
	{
		day   = as.numeric(substr(xxx, 1, 2))
		month = as.numeric(substr(xxx, 4, 5))
		year  = as.numeric(substr(xxx, 7, 10))

		if( year==2013 ){ num_date <- 0 }
		if( year==2014 ){ num_date <- 1*365 }
		if( year==2015 ){ num_date <- 2*365 }

		if( month==1  ){ num_date <- num_date + 0 }
		if( month==2  ){ num_date <- num_date + 31 }
		if( month==3  ){ num_date <- num_date + 59 }
		if( month==4  ){ num_date <- num_date + 90 }
		if( month==5  ){ num_date <- num_date + 120 }
		if( month==6  ){ num_date <- num_date + 151 }
		if( month==7  ){ num_date <- num_date + 181 }
		if( month==8  ){ num_date <- num_date + 212 }
		if( month==9  ){ num_date <- num_date + 243 }
		if( month==10  ){ num_date <- num_date + 273 }
		if( month==11 ){ num_date <- num_date + 304 }
		if( month==12 ){ num_date <- num_date + 334 }
	
		num_date <- num_date + day
	}

	if( (NN==8) && (slash[1]==3) && (slash[2]==6) )
	{
		day   = as.numeric(substr(xxx, 1, 2))
		month = as.numeric(substr(xxx, 4, 5))
		year  = as.numeric(substr(xxx, 7, 8)) + 2000

		if( year==2013 ){ num_date <- 0 }
		if( year==2014 ){ num_date <- 1*365 }
		if( year==2015 ){ num_date <- 2*365 }

		if( month==1  ){ num_date <- num_date + 0 }
		if( month==2  ){ num_date <- num_date + 31 }
		if( month==3  ){ num_date <- num_date + 59 }
		if( month==4  ){ num_date <- num_date + 90 }
		if( month==5  ){ num_date <- num_date + 120 }
		if( month==6  ){ num_date <- num_date + 151 }
		if( month==7  ){ num_date <- num_date + 181 }
		if( month==8  ){ num_date <- num_date + 212 }
		if( month==9  ){ num_date <- num_date + 243 }
		if( month==10  ){ num_date <- num_date + 273 }
		if( month==11 ){ num_date <- num_date + 304 }
		if( month==12 ){ num_date <- num_date + 334 }
	
		num_date <- num_date + day
	}

	
	num_date
} 



#########################################
## Read in raw data from Brazil

Braz_data_read <- read.csv( "C:\\U\\GHIT\\2_All_Data\\2_Data\\data\\raw\\Brazil\\Manaus_finalbleed_SK.csv" )

Braz_PvPCR <- as.data.frame( Braz_data_read[,c(24,5)] )

Braz_PvPCR <- cbind( Braz_PvPCR, rep(NA, nrow(Braz_PvPCR)) )
colnames(Braz_PvPCR)[3] <- "studyday"

for(i in 1:nrow(Braz_PvPCR))
{
	Braz_PvPCR[i,3] <- braz_date_convert( Braz_PvPCR[i,2] )
}

Braz_PvPCR[,3] <- Braz_PvPCR[,3] - min(Braz_PvPCR[,3])

Braz_PvPCR <- Braz_PvPCR[-which(Braz_PvPCR[,3] > 400),]


Braz_bin_edges <- seq(from=-14, by=28, length=15) 
N_braz_bins <- length(Braz_bin_edges) - 1

Braz_bin_edges[N_braz_bins+1] <- max(Braz_PvPCR[,3])



hist( Braz_PvPCR[,3], breaks=200, col="grey",
main="Brazilian samples", xlab="time (days)" )


for(j in 1:length(Braz_bin_edges))
{
	points(x=rep(Braz_bin_edges[j],2), y=c(0,2000), type='l', col="red", lwd=2)	
}

Braz_PvPCR_bins <- matrix(NA, nrow=3, ncol=N_thai_bins)
for(j in 1:N_braz_bins)
{
	index_j <- which( (Braz_PvPCR[,3] > Braz_bin_edges[j]) & (Braz_PvPCR[,3] <= Braz_bin_edges[j+1]) ) 

	Braz_PvPCR_bins[,j] <- as.numeric(as.vector(
            	                	 binom.confint( sum(Braz_PvPCR[index_j,1]), length(Braz_PvPCR[index_j,1]), method="wilson")[1,4:6]
            	                ))
}

Braz_PvPCR_times <- 0.5*( Braz_bin_edges[-1] + Braz_bin_edges[-length(Braz_bin_edges)] ) 

Braz_PvPCR_times <- 2019.95 - max(Braz_PvPCR_times/365) + Braz_PvPCR_times/365



################
################
##            ##
##  Solomons  ##
##            ##
################
################


#########################################
## Create a function for data conversion

sol_date_convert = function( xxx )
{
	num_date = NA

	if( nchar(as.vector(xxx))==10 )
	{
		if( substr( xxx, 3, 3 ) == "/" )
		{
			day   = as.numeric(substr(xxx, 1, 2))
			month = as.numeric(substr(xxx, 4, 5))
			year  = as.numeric(substr(xxx, 7, 10))
		}

		if( year==2013 ){ num_date <- 0 }
		if( year==2014 ){ num_date <- 1*365 }

		if( month==1  ){ num_date <- num_date + 0 }
		if( month==2  ){ num_date <- num_date + 31 }
		if( month==3  ){ num_date <- num_date + 59 }
		if( month==4  ){ num_date <- num_date + 90 }
		if( month==5  ){ num_date <- num_date + 120 }
		if( month==6  ){ num_date <- num_date + 151 }
		if( month==7  ){ num_date <- num_date + 181 }
		if( month==8  ){ num_date <- num_date + 212 }
		if( month==9  ){ num_date <- num_date + 243 }
		if( month==10  ){ num_date <- num_date + 273 }
		if( month==11 ){ num_date <- num_date + 304 }
		if( month==12 ){ num_date <- num_date + 334 }
	
		num_date <- num_date + day
	}
	
	num_date
} 



#########################################
## Read in raw data from Solomons

Sol_data_read <- read.csv( "C:\\U\\GHIT\\2_All_Data\\2_Data\\data\\raw\\Solomon_Islands\\SI_ACD_data.csv" )

Sol_PvPCR <- as.data.frame( Sol_data_read[,c(37,3)] )

Sol_PvPCR <- cbind( Sol_PvPCR, rep(NA, nrow(Sol_PvPCR)) )
colnames(Sol_PvPCR)[3] <- "studyday"

for(i in 1:nrow(Sol_PvPCR))
{
	Sol_PvPCR[i,3] <- sol_date_convert( Sol_PvPCR[i,2] )
}

Sol_PvPCR[,3] <- Sol_PvPCR[,3] - min(Sol_PvPCR[,3])


Sol_bin_edges <- seq(from=-14, by=28, length=14) 
N_sol_bins <- length(Sol_bin_edges) - 1

Sol_bin_edges[N_sol_bins+1] <- max(Sol_PvPCR[,3])



hist( Sol_PvPCR[,3], breaks=200, col="grey",
main="Solomon Islands samples", "time (days)" )

for(j in 1:length(Sol_bin_edges))
{
	points(x=rep(Sol_bin_edges[j],2), y=c(0,2000), type='l', col="red", lwd=2)	
}

Sol_PvPCR_bins <- matrix(NA, nrow=3, ncol=N_sol_bins)
for(j in 1:N_sol_bins)
{
	index_j <- which( (Sol_PvPCR[,3] > Sol_bin_edges[j]) & (Sol_PvPCR[,3] <= Sol_bin_edges[j+1]) & (is.na(Sol_PvPCR[,1])==FALSE) ) 

	Sol_PvPCR_bins[,j] <- as.numeric(as.vector(
            	                	 binom.confint( sum(Sol_PvPCR[index_j,1]), length(Sol_PvPCR[index_j,1]), method="wilson")[1,4:6]
            	                ))
}

Sol_PvPCR_times <- 0.5*( Sol_bin_edges[-1] + Sol_bin_edges[-length(Sol_bin_edges)] ) 

Sol_PvPCR_times <- 2019.95 - max(Sol_PvPCR_times/365) + Sol_PvPCR_times/365




Sol_age_check <- read.csv( "C:\\U\\GHIT\\2_All_Data\\2_Data\\data\\raw\\Solomon_Islands\\20170124_SI_epi_ab_combined.csv" )

quantile( Sol_age_check[,8] )





##################################
##################################
##                              ## 
##  #####  ##     ####  ######  ##
##  ##  ## ##    ##  ##   ##    ##
##  #####  ##    ##  ##   ##    ##
##  ##     ##    ##  ##   ##    ##
##  ##     #####  ####    ##    ##
##                              ##
##################################
##################################


ant_cols = c("red", "orange", "gold", "yellowgreen", "forestgreen", "dodgerblue", "royalblue", "blue4")



MDA_cols = c("grey39", "sienna1", "forestgreen", "deepskyblue", "purple3" )



tiff(file="Figure5_ROC_curves_and_TPP_big_labels.tif", width=32, height=22, units="cm", res=500)


lay.mat <- rbind( c( 1, 2, 3, 4 ),
                  c( 5, 5, 5, 5 ), 
                  c( 6, 8,10,12 ),
                  c( 7, 9,11,13 ) )
layout(lay.mat, heights=c(10,1,8,8))
layout.show(13)


par(mar=c(5,5,2,1))
par(mgp=c(3.3, 1,0))


line_seq <- c(0.2, 0.4, 0.6, 0.8)


line.size = 2
lab.size  = 2
axis.size = 1
main.size = 1.1

pc_axis_size <- 1.5


######################################
######################################
##                                  ##
##  PANEL 1                         ##
##  Thailand predicting Thailand    ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="(A) training = Thailand; testing = Thailand",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=0.9*pc_axis_size ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=pc_axis_size ) 


for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

AB_ENSEMBLE_thai_thai_chull[[2]] = AB_ENSEMBLE_thai_thai_chull[[2]][-62,]

for(i in 1:8)
{
	points( x=AB_ENSEMBLE_thai_thai_chull[[i]][,1], 
      	  y=AB_ENSEMBLE_thai_thai_chull[[i]][,2], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}


points(x=1-Thai_SS[1,2], y=Thai_SS[1,1], pch=17, cex=2, col=MDA_cols[3])
points(x=1-Thai_SS[2,2], y=Thai_SS[2,1], pch=17, cex=2, col=MDA_cols[4])
points(x=1-Thai_SS[3,2], y=Thai_SS[3,1], pch=17, cex=2, col=MDA_cols[5])



######################################
######################################
##                                  ##
##  PANEL 2                         ##
##  Brazil predicting Brazil        ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="(B) training = Brazil; testing = Brazil",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=0.9*pc_axis_size ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=pc_axis_size ) 



for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 1:8)
{
	points( x=AB_ENSEMBLE_braz_braz_chull[[i]][,1], 
      	  y=AB_ENSEMBLE_braz_braz_chull[[i]][,2], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}




points(x=1-Braz_SS[1,2], y=Braz_SS[1,1], pch=17, cex=2, col=MDA_cols[3])
points(x=1-Braz_SS[2,2], y=Braz_SS[2,1], pch=17, cex=2, col=MDA_cols[4])
points(x=1-Braz_SS[3,2], y=0.99*Braz_SS[3,1], pch=17, cex=2, col=MDA_cols[5])






######################################
######################################
##                                  ##
##  PANEL 3                         ##
##  Solomons predicting Solomons    ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="(C) training = Solomons; testing = Solomons",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=0.9*pc_axis_size ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=pc_axis_size ) 



for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

for(i in 1:8)
{
	points( x=AB_ENSEMBLE_sol_sol_chull[[i]][,1], 
      	  y=AB_ENSEMBLE_sol_sol_chull[[i]][,2], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}





points(x=1-Sol_SS[1,2], y=Sol_SS[1,1], pch=17, cex=2, col=MDA_cols[3])
points(x=1-Sol_SS[2,2], y=Sol_SS[2,1], pch=17, cex=2, col=MDA_cols[4])
points(x=1-Sol_SS[3,2], y=Sol_SS[3,1], pch=17, cex=2, col=MDA_cols[5])



######################################
######################################
##                                  ##
##  PANEL 4                         ##
##  All predicting All              ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="(D) training = all regions; testing = all regions",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=0.9*pc_axis_size ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=pc_axis_size ) 


for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 1:8)
{
	points( x=AB_ENSEMBLE_all_all_chull[[i]][,1], 
      	  y=AB_ENSEMBLE_all_all_chull[[i]][,2], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}




points(x=1-All_SS[1,2], y=All_SS[1,1], pch=17, cex=2, col=MDA_cols[3])
points(x=1-All_SS[2,2], y=All_SS[2,1], pch=17, cex=2, col=MDA_cols[4])
points(x=1-All_SS[3,2], y=All_SS[3,1], pch=17, cex=2, col=MDA_cols[5])



######################################
######################################
##                                  ##
##  LEGEND                          ##
##                                  ##                                  
######################################
######################################

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("1 antigen", 
                  "2 antigens", 
                  "3 antigens",
                  "4 antigens",
                  "5 antigens",
                  "6 antigens",
                  "7 antigens", 
                  "8 antigens"), 
       fill = ant_cols,
       border = ant_cols,
       ncol=8, cex=1.75, bty="n" )



par(mar=c(3,4.5,2,1))
par(mgp=c(3, 0.55,0))



######################################
######################################
##                                  ## 
##  PANEL 5                         ##
##  PQ targeting: Thailand          ##
##                                  ##
######################################   
######################################

width = 1
edge = 0.1



plot(x=10000, y=10000, 
xlim=c(0,5*(width+edge)+edge), ylim=c(0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="", ylab="sensitivity",
main="(E) Thailand: hypnozoite carriers targeted",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)



PPP = 1

polygon( x=edge*1+c((1-1)*width,1*width,1*width,(1-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[1], border=MDA_cols[1] )


PPP = length(which(inf_cat_thai=="thai_current"))/length(which(bin_cat_thai=="new"))

polygon( x=edge*2+c((2-1)*width,2*width,2*width,(2-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[2], border=MDA_cols[2] )


PPP = Thai_SS[1,1]

polygon( x=edge*3+c((3-1)*width,3*width,3*width,(3-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[3], border=MDA_cols[3] )


PPP = Thai_SS[2,1]

polygon( x=edge*4+c((4-1)*width,4*width,4*width,(4-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[4], border=MDA_cols[4] )


PPP = Thai_SS[3,1]

polygon( x=edge*5+c((5-1)*width,5*width,5*width,(5-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[5], border=MDA_cols[5] )


for(i in 1:4)
{
	points(x=c(0,10), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
}

##points(x=c(0,100), y=c(0.8, 0.8), 
##type='l', lty="dashed", lwd=2)

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=0.9*pc_axis_size ) 

axis(1, 
at = c(0.6, 1.7, 2.8, 3.9, 5.0),
labels=c("MDA", "MSAT: PCR",
                  "SeroTAT_80_80", "SeroTAT_50_98", "SeroTAT_98_50"), cex.axis=0.5 )





######################################
######################################
##                                  ## 
##  PANEL 6                         ##
##  PQ targeting: Thailand          ##
##                                  ##
######################################   
######################################



plot(x=10000, y=10000, 
xlim=c(0,5*(width+edge)+edge), ylim=c(0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="", ylab="1 - specificity",
main="(I) Thailand: over-treatment",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


PPP = 1 - length(which(inf_cat_thai %in% c("thai_current", "thai_recent")))/length(inf_cat_thai)

polygon( x=edge*1+c((1-1)*width,1*width,1*width,(1-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[1], border=MDA_cols[1] )


PPP = 0

polygon( x=edge*2+c((2-1)*width,2*width,2*width,(2-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[2], border=MDA_cols[2] )


PPP = 1-Thai_SS[1,2]

polygon( x=edge*3+c((3-1)*width,3*width,3*width,(3-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[3], border=MDA_cols[3] )


PPP = 1-Thai_SS[2,2]

polygon( x=edge*4+c((4-1)*width,4*width,4*width,(4-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[4], border=MDA_cols[4] )


PPP = 1-Thai_SS[3,2]

polygon( x=edge*5+c((5-1)*width,5*width,5*width,(5-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[5], border=MDA_cols[5] )

for(i in 1:4)
{
	points(x=c(0,10), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
}


##points(x=c(0,100), y=c(0.8, 0.8), 
##type='l', lty="dashed", lwd=2)

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=0.9*pc_axis_size ) 

axis(1, 
at = c(0.6, 1.7, 2.8, 3.9, 5.0),
labels=c("MDA", "MSAT: PCR",
                  "SeroTAT_80_80", "SeroTAT_50_98", "SeroTAT_98_50"), cex.axis=0.5 )






######################################
######################################
##                                  ## 
##  PANEL 7                         ##
##  PQ targeting: Brazil            ##
##                                  ##
######################################   
######################################



plot(x=10000, y=10000, 
xlim=c(0,5*(width+edge)+edge), ylim=c(0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="", ylab="sensitivity",
main="(F) Brazil: hypnozoite carriers targeted",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)



PPP = 1

polygon( x=edge*1+c((1-1)*width,1*width,1*width,(1-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[1], border=MDA_cols[1] )


PPP = length(which(inf_cat_braz=="braz_current"))/length(which(bin_cat_braz=="new"))

polygon( x=edge*2+c((2-1)*width,2*width,2*width,(2-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[2], border=MDA_cols[2] )


PPP = Braz_SS[1,1]

polygon( x=edge*3+c((3-1)*width,3*width,3*width,(3-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[3], border=MDA_cols[3] )


PPP = Braz_SS[2,1]

polygon( x=edge*4+c((4-1)*width,4*width,4*width,(4-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[4], border=MDA_cols[4] )


PPP = Braz_SS[3,1]

polygon( x=edge*5+c((5-1)*width,5*width,5*width,(5-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[5], border=MDA_cols[5] )


for(i in 1:4)
{
	points(x=c(0,10), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
}

##points(x=c(0,100), y=c(0.8, 0.8), 
##type='l', lty="dashed", lwd=2)

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=0.9*pc_axis_size ) 

axis(1, 
at = c(0.6, 1.7, 2.8, 3.9, 5.0),
labels=c("MDA", "MSAT: PCR",
                  "SeroTAT_80_80", "SeroTAT_50_98", "SeroTAT_98_50"), cex.axis=0.5 )




######################################
######################################
##                                  ## 
##  PANEL 8                         ##
##  PQ targeting: Brazil            ##
##                                  ##
######################################   
######################################
 

plot(x=10000, y=10000, 
xlim=c(0,5*(width+edge)+edge), ylim=c(0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="", ylab="1 - specificity",
main="(J) Brazil: over-treatment",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)



PPP = 1 - length(which(inf_cat_braz %in% c("braz_current", "braz_recent")))/length(inf_cat_braz)


polygon( x=edge*1+c((1-1)*width,1*width,1*width,(1-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[1], border=MDA_cols[1] )


PPP = 0

polygon( x=edge*2+c((2-1)*width,2*width,2*width,(2-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[2], border=MDA_cols[2] )


PPP = 1-Braz_SS[1,2]

polygon( x=edge*3+c((3-1)*width,3*width,3*width,(3-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[3], border=MDA_cols[3] )


PPP = 1-Braz_SS[2,2]

polygon( x=edge*4+c((4-1)*width,4*width,4*width,(4-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[4], border=MDA_cols[4] )


PPP = 1-Braz_SS[3,2]

polygon( x=edge*5+c((5-1)*width,5*width,5*width,(5-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[5], border=MDA_cols[5] )


for(i in 1:4)
{
	points(x=c(0,10), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
}

##points(x=c(0,100), y=c(0.8, 0.8), 
##type='l', lty="dashed", lwd=2)

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=0.9*pc_axis_size ) 

axis(1, 
at = c(0.6, 1.7, 2.8, 3.9, 5.0),
labels=c("MDA", "MSAT: PCR",
                  "SeroTAT_80_80", "SeroTAT_50_98", "SeroTAT_98_50"), cex.axis=0.5 )




######################################
######################################
##                                  ## 
##  PANEL 9                         ##
##  PQ targeting: Solomons          ##
##                                  ##
######################################   
######################################

plot(x=10000, y=10000, 
xlim=c(0,5*(width+edge)+edge), ylim=c(0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="", ylab="sensitivity",
main="(G) Solomons: hypnozoite carriers targeted",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


PPP = 1

polygon( x=edge*1+c((1-1)*width,1*width,1*width,(1-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[1], border=MDA_cols[1] )


PPP = length(which(inf_cat_sol=="sol_current"))/length(which(bin_cat_sol=="new"))

polygon( x=edge*2+c((2-1)*width,2*width,2*width,(2-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[2], border=MDA_cols[2] )


PPP = Sol_SS[1,1]

polygon( x=edge*3+c((3-1)*width,3*width,3*width,(3-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[3], border=MDA_cols[3] )


PPP = Sol_SS[2,1]

polygon( x=edge*4+c((4-1)*width,4*width,4*width,(4-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[4], border=MDA_cols[4] )


PPP = Sol_SS[3,1]

polygon( x=edge*5+c((5-1)*width,5*width,5*width,(5-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[5], border=MDA_cols[5] )


for(i in 1:4)
{
	points(x=c(0,10), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
}

##points(x=c(0,100), y=c(0.8, 0.8), 
##type='l', lty="dashed", lwd=2)

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=0.9*pc_axis_size ) 

axis(1, 
at = c(0.6, 1.7, 2.8, 3.9, 5.0),
labels=c("MDA", "MSAT: PCR",
                  "SeroTAT_80_80", "SeroTAT_50_98", "SeroTAT_98_50"), cex.axis=0.5 )






######################################
######################################
##                                  ## 
##  PANEL 10                        ##
##  PQ targeting: Solomons          ##
##                                  ##
######################################   
######################################


plot(x=10000, y=10000, 
xlim=c(0,5*(width+edge)+edge), ylim=c(0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="", ylab="1 - specificity",
main="(K) Solomons: over-treatment",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


PPP = 1 - length(which(inf_cat_sol %in% c("sol_current", "sol_recent")))/length(inf_cat_sol)


polygon( x=edge*1+c((1-1)*width,1*width,1*width,(1-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[1], border=MDA_cols[1] )


PPP = 0

polygon( x=edge*2+c((2-1)*width,2*width,2*width,(2-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[2], border=MDA_cols[2] )


PPP =  1-Sol_SS[1,2]

polygon( x=edge*3+c((3-1)*width,3*width,3*width,(3-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[3], border=MDA_cols[3] )


PPP = 1-Sol_SS[2,2]

polygon( x=edge*4+c((4-1)*width,4*width,4*width,(4-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[4], border=MDA_cols[4] )


PPP = 1-Sol_SS[3,2]

polygon( x=edge*5+c((5-1)*width,5*width,5*width,(5-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[5], border=MDA_cols[5] )


for(i in 1:4)
{
	points(x=c(0,10), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
}

##points(x=c(0,100), y=c(0.8, 0.8), 
##type='l', lty="dashed", lwd=2)

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=0.9*pc_axis_size ) 

axis(1, 
at = c(0.6, 1.7, 2.8, 3.9, 5.0),
labels=c("MDA", "MSAT: PCR",
                  "SeroTAT_80_80", "SeroTAT_50_98", "SeroTAT_98_50"), cex.axis=0.5 )







######################################
######################################
##                                  ## 
##  PANEL 11                        ##
##  PQ targeting: All               ##
##                                  ##
######################################   
######################################


plot(x=10000, y=10000, 
xlim=c(0,5*(width+edge)+edge), ylim=c(0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="", ylab="sensitivity",
main="(H) All regions: hypnozoite carriers targeted",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)



PPP = 1

polygon( x=edge*1+c((1-1)*width,1*width,1*width,(1-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[1], border=MDA_cols[1] )


PPP = length(which(inf_cat_sol=="sol_current"))/length(which(bin_cat_sol=="new"))

polygon( x=edge*2+c((2-1)*width,2*width,2*width,(2-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[2], border=MDA_cols[2] )


PPP = All_SS[1,1]

polygon( x=edge*3+c((3-1)*width,3*width,3*width,(3-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[3], border=MDA_cols[3] )


PPP = All_SS[2,1]

polygon( x=edge*4+c((4-1)*width,4*width,4*width,(4-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[4], border=MDA_cols[4] )


PPP = All_SS[3,1]

polygon( x=edge*5+c((5-1)*width,5*width,5*width,(5-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[5], border=MDA_cols[5] )


for(i in 1:4)
{
	points(x=c(0,10), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
}

##points(x=c(0,100), y=c(0.8, 0.8), 
##type='l', lty="dashed", lwd=2)

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=0.9*pc_axis_size ) 

axis(1, 
at = c(0.6, 1.7, 2.8, 3.9, 5.0),
labels=c("MDA", "MSAT: PCR",
                  "SeroTAT_80_80", "SeroTAT_50_98", "SeroTAT_98_50"), cex.axis=0.5 )






######################################
######################################
##                                  ## 
##  PANEL 12                        ##
##  PQ targeting: All               ##
##                                  ##
######################################   
######################################


plot(x=10000, y=10000, 
xlim=c(0,5*(width+edge)+edge), ylim=c(0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="", ylab="1 - specificity",
main="(L) All regions: over-treatment",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)



PPP = 1 - length(which(inf_cat_all %in% c("thai_current", "thai_recent", "braz_current", "braz_recent", "sol_current", "sol_recent")))/length(inf_cat_all)

polygon( x=edge*1+c((1-1)*width,1*width,1*width,(1-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[1], border=MDA_cols[1] )


PPP = 0

polygon( x=edge*2+c((2-1)*width,2*width,2*width,(2-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[2], border=MDA_cols[2] )


PPP = 1-All_SS[1,2]

polygon( x=edge*3+c((3-1)*width,3*width,3*width,(3-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[3], border=MDA_cols[3] )


PPP = 1-All_SS[2,2]

polygon( x=edge*4+c((4-1)*width,4*width,4*width,(4-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[4], border=MDA_cols[4] )


PPP = 1-All_SS[3,2]

polygon( x=edge*5+c((5-1)*width,5*width,5*width,(5-1)*width), 
         y=c(0,0,PPP,PPP), 
         col=MDA_cols[5], border=MDA_cols[5] )


for(i in 1:4)
{
	points(x=c(0,10), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
}

##points(x=c(0,100), y=c(0.8, 0.8), 
##type='l', lty="dashed", lwd=2)

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), las=2, cex.axis=0.9*pc_axis_size ) 

axis(1, 
at = c(0.6, 1.7, 2.8, 3.9, 5.0),
labels=c("MDA", "MSAT: PCR",
                  "SeroTAT_80_80", "SeroTAT_50_98", "SeroTAT_98_50"), cex.axis=0.5 )





dev.off()










