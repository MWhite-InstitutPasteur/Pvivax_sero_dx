
load("C:\\U\\GHIT\\NatMed_Paper\\Figure5_algorithm\\Algorithm\\AB_LOGIT.RData")

load("C:\\U\\GHIT\\NatMed_Paper\\Figure5_algorithm\\Algorithm\\AB_LDA.RData")

load("C:\\U\\GHIT\\NatMed_Paper\\Figure5_algorithm\\Algorithm\\AB_QDA.RData")

load("C:\\U\\GHIT\\NatMed_Paper\\Figure5_algorithm\\Algorithm\\AB_Dtree.RData")

load("C:\\U\\GHIT\\NatMed_Paper\\Figure5_algorithm\\Algorithm\\AB_Rforest.RData")

load("C:\\U\\GHIT\\NatMed_Paper\\Figure5_algorithm\\Algorithm\\AB_DYN1.RData")

load("C:\\U\\GHIT\\NatMed_Paper\\Figure5_algorithm\\Algorithm\\AB_DYN2.RData")


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
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_thai_thai[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_thai_thai[[i]][[j]][,1] )
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
##  PART 2                          ##
##  Thailand predicting Brazil      ## 
##                                  ##                                  
######################################
######################################


AB_ENSEMBLE_thai_braz = list()

for(i in 1:7)
{
	AB_ENSEMBLE_thai_braz[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_thai_braz[[i]]))
	{
		AB_ENSEMBLE_thai_braz[[i]][[count]] = AB_LOGIT_thai_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_thai_braz[[i]]))
	{
		AB_ENSEMBLE_thai_braz[[i]][[count]] = AB_LDA_thai_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_thai_braz[[i]]))
	{
		AB_ENSEMBLE_thai_braz[[i]][[count]] = AB_QDA_thai_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_thai_braz[[i]]))
	{
		AB_ENSEMBLE_thai_braz[[i]][[count]] = AB_Dtree_thai_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_thai_braz[[i]]))
	{
		AB_ENSEMBLE_thai_braz[[i]][[count]] = AB_Rforest_thai_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_thai_braz[[i]]))
	{
		AB_ENSEMBLE_thai_braz[[i]][[count]] = AB_DYN1_thai_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_thai_braz[[i]]))
	{
		AB_ENSEMBLE_thai_braz[[i]][[count]] = AB_DYN2_thai_braz[[i]][[j]][[1]]

		count = count + 1
	}	
}


AB_ENSEMBLE_thai_braz_chull = list()



AB_ENSEMBLE_thai_braz_chull[[1]] = AB_braz_1

XX_VEC = AB_braz_1[,1]
YY_VEC = AB_braz_1[,2]

for(i in 1:length(AB_ENSEMBLE_thai_braz))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_thai_braz[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_thai_braz[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_thai_braz[[i]][[j]][,1] )
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


	AB_ENSEMBLE_thai_braz_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 3                          ##
##  Thailand predicting Solomons    ## 
##                                  ##                                  
######################################
######################################

AB_ENSEMBLE_thai_sol = list()

for(i in 1:7)
{
	AB_ENSEMBLE_thai_sol[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_thai_sol[[i]]))
	{
		AB_ENSEMBLE_thai_sol[[i]][[count]] = AB_LOGIT_thai_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_thai_sol[[i]]))
	{
		AB_ENSEMBLE_thai_sol[[i]][[count]] = AB_LDA_thai_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_thai_sol[[i]]))
	{
		AB_ENSEMBLE_thai_sol[[i]][[count]] = AB_QDA_thai_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_thai_sol[[i]]))
	{
		AB_ENSEMBLE_thai_sol[[i]][[count]] = AB_Dtree_thai_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_thai_sol[[i]]))
	{
		AB_ENSEMBLE_thai_sol[[i]][[count]] = AB_Rforest_thai_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_thai_sol[[i]]))
	{
		AB_ENSEMBLE_thai_sol[[i]][[count]] = AB_DYN1_thai_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_thai_sol[[i]]))
	{
		AB_ENSEMBLE_thai_sol[[i]][[count]] = AB_DYN2_thai_sol[[i]][[j]][[1]]

		count = count + 1
	}	
}


AB_ENSEMBLE_thai_sol_chull = list()

AB_ENSEMBLE_thai_sol_chull[[1]] = AB_sol_1

XX_VEC = AB_sol_1[,1]
YY_VEC = AB_sol_1[,2]

for(i in 1:length(AB_ENSEMBLE_thai_sol))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_thai_sol[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_thai_sol[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_thai_sol[[i]][[j]][,1] )
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


	AB_ENSEMBLE_thai_sol_chull[[1+i]] = XY_hull_plot 
}





######################################
######################################
##                                  ##
##  PART 4                          ##
##  Thailand predicting All         ## 
##                                  ##                                  
######################################
######################################

AB_ENSEMBLE_thai_all = list()

for(i in 1:7)
{
	AB_ENSEMBLE_thai_all[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_thai_all[[i]]))
	{
		AB_ENSEMBLE_thai_all[[i]][[count]] = AB_LOGIT_thai_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_thai_all[[i]]))
	{
		AB_ENSEMBLE_thai_all[[i]][[count]] = AB_LDA_thai_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_thai_all[[i]]))
	{
		AB_ENSEMBLE_thai_all[[i]][[count]] = AB_QDA_thai_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_thai_all[[i]]))
	{
		AB_ENSEMBLE_thai_all[[i]][[count]] = AB_Dtree_thai_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_thai_all[[i]]))
	{
		AB_ENSEMBLE_thai_all[[i]][[count]] = AB_Rforest_thai_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_thai_all[[i]]))
	{
		AB_ENSEMBLE_thai_all[[i]][[count]] = AB_DYN1_thai_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_thai_all[[i]]))
	{
		AB_ENSEMBLE_thai_all[[i]][[count]] = AB_DYN2_thai_all[[i]][[j]]

		count = count + 1
	}	
}


AB_ENSEMBLE_thai_all_chull = list()

AB_ENSEMBLE_thai_all_chull[[1]] = AB_all_1

XX_VEC = AB_all_1[,1]
YY_VEC = AB_all_1[,2]


for(i in 1:length(AB_ENSEMBLE_thai_all))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_thai_all[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_thai_all[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_thai_all[[i]][[j]][,1] )
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


	AB_ENSEMBLE_thai_all_chull[[1+i]] = XY_hull_plot 
}


######################################
######################################
##                                  ##
##  PART 5                          ##
##  Brazil predicting Thailand      ## 
##                                  ##                                  
######################################
######################################


AB_ENSEMBLE_braz_thai = list()

for(i in 1:7)
{
	AB_ENSEMBLE_braz_thai[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_braz_thai[[i]]))
	{
		AB_ENSEMBLE_braz_thai[[i]][[count]] = AB_LOGIT_braz_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_braz_thai[[i]]))
	{
		AB_ENSEMBLE_braz_thai[[i]][[count]] = AB_LDA_braz_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_braz_thai[[i]]))
	{
		AB_ENSEMBLE_braz_thai[[i]][[count]] = AB_QDA_braz_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_braz_thai[[i]]))
	{
		AB_ENSEMBLE_braz_thai[[i]][[count]] = AB_Dtree_braz_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_braz_thai[[i]]))
	{
		AB_ENSEMBLE_braz_thai[[i]][[count]] = AB_Rforest_braz_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_braz_thai[[i]]))
	{
		AB_ENSEMBLE_braz_thai[[i]][[count]] = AB_DYN1_braz_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_braz_thai[[i]]))
	{
		AB_ENSEMBLE_braz_thai[[i]][[count]] = AB_DYN2_braz_thai[[i]][[j]][[1]]

		count = count + 1
	}	
}


AB_ENSEMBLE_braz_thai_chull = list()


AB_ENSEMBLE_braz_thai_chull[[1]] = AB_thai_1

XX_VEC = AB_thai_1[,1]
YY_VEC = AB_thai_1[,2]

for(i in 1:length(AB_ENSEMBLE_braz_thai))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_braz_thai[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_braz_thai[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_braz_thai[[i]][[j]][,1] )
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


	AB_ENSEMBLE_braz_thai_chull[[1+i]] = XY_hull_plot 
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
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_braz_braz[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_braz_braz[[i]][[j]][,1] )
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
##  PART 7                          ##
##  Brazil predicting Solomons      ## 
##                                  ##                                  
######################################
######################################



AB_ENSEMBLE_braz_sol = list()

for(i in 1:7)
{
	AB_ENSEMBLE_braz_sol[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_braz_sol[[i]]))
	{
		AB_ENSEMBLE_braz_sol[[i]][[count]] = AB_LOGIT_braz_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_braz_sol[[i]]))
	{
		AB_ENSEMBLE_braz_sol[[i]][[count]] = AB_LDA_braz_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_braz_sol[[i]]))
	{
		AB_ENSEMBLE_braz_sol[[i]][[count]] = AB_QDA_braz_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_braz_sol[[i]]))
	{
		AB_ENSEMBLE_braz_sol[[i]][[count]] = AB_Dtree_braz_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_braz_sol[[i]]))
	{
		AB_ENSEMBLE_braz_sol[[i]][[count]] = AB_Rforest_braz_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_braz_sol[[i]]))
	{
		AB_ENSEMBLE_braz_sol[[i]][[count]] = AB_DYN1_braz_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_braz_sol[[i]]))
	{
		AB_ENSEMBLE_braz_sol[[i]][[count]] = AB_DYN2_braz_sol[[i]][[j]][[1]]

		count = count + 1
	}	
}


AB_ENSEMBLE_braz_sol_chull = list()

AB_ENSEMBLE_braz_sol_chull[[1]] = AB_sol_1

XX_VEC = AB_sol_1[,1]
YY_VEC = AB_sol_1[,2]

for(i in 1:length(AB_ENSEMBLE_braz_sol))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_braz_sol[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_braz_sol[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_braz_sol[[i]][[j]][,1] )
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


	AB_ENSEMBLE_braz_sol_chull[[1+i]] = XY_hull_plot 
}




######################################
######################################
##                                  ##
##  PART 8                          ##
##  Brazil predicting All           ## 
##                                  ##                                  
######################################
######################################


AB_ENSEMBLE_braz_all = list()

for(i in 1:7)
{
	AB_ENSEMBLE_braz_all[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_braz_all[[i]]))
	{
		AB_ENSEMBLE_braz_all[[i]][[count]] = AB_LOGIT_braz_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_braz_all[[i]]))
	{
		AB_ENSEMBLE_braz_all[[i]][[count]] = AB_LDA_braz_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_braz_all[[i]]))
	{
		AB_ENSEMBLE_braz_all[[i]][[count]] = AB_QDA_braz_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_braz_all[[i]]))
	{
		AB_ENSEMBLE_braz_all[[i]][[count]] = AB_Dtree_braz_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_braz_all[[i]]))
	{
		AB_ENSEMBLE_braz_all[[i]][[count]] = AB_Rforest_braz_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_braz_all[[i]]))
	{
		AB_ENSEMBLE_braz_all[[i]][[count]] = AB_DYN1_braz_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_braz_all[[i]]))
	{
		AB_ENSEMBLE_braz_all[[i]][[count]] = AB_DYN2_braz_all[[i]][[j]]

		count = count + 1
	}	
}


AB_ENSEMBLE_braz_all_chull = list()



AB_ENSEMBLE_braz_all_chull[[1]] = AB_all_1

XX_VEC = AB_all_1[,1]
YY_VEC = AB_all_1[,2]


for(i in 1:length(AB_ENSEMBLE_braz_all))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_braz_all[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_braz_all[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_braz_all[[i]][[j]][,1] )
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


	AB_ENSEMBLE_braz_all_chull[[1+i]] = XY_hull_plot 
}




######################################
######################################
##                                  ##
##  PART 9                          ##
##  Solomons predicting Thailand    ## 
##                                  ##                                  
######################################
######################################



AB_ENSEMBLE_sol_thai = list()

for(i in 1:7)
{
	AB_ENSEMBLE_sol_thai[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_sol_thai[[i]]))
	{
		AB_ENSEMBLE_sol_thai[[i]][[count]] = AB_LOGIT_sol_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_sol_thai[[i]]))
	{
		AB_ENSEMBLE_sol_thai[[i]][[count]] = AB_LDA_sol_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_sol_thai[[i]]))
	{
		AB_ENSEMBLE_sol_thai[[i]][[count]] = AB_QDA_sol_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_sol_thai[[i]]))
	{
		AB_ENSEMBLE_sol_thai[[i]][[count]] = AB_Dtree_sol_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_sol_thai[[i]]))
	{
		AB_ENSEMBLE_sol_thai[[i]][[count]] = AB_Rforest_sol_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_sol_thai[[i]]))
	{
		AB_ENSEMBLE_sol_thai[[i]][[count]] = AB_DYN1_sol_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_sol_thai[[i]]))
	{
		AB_ENSEMBLE_sol_thai[[i]][[count]] = AB_DYN2_sol_thai[[i]][[j]][[1]]

		count = count + 1
	}	
}


AB_ENSEMBLE_sol_thai_chull = list()

AB_ENSEMBLE_sol_thai_chull[[1]] = AB_thai_1

XX_VEC = AB_thai_1[,1]
YY_VEC = AB_thai_1[,2]

for(i in 1:length(AB_ENSEMBLE_sol_thai))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_sol_thai[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_sol_thai[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_sol_thai[[i]][[j]][,1] )
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


	AB_ENSEMBLE_sol_thai_chull[[1+i]] = XY_hull_plot 
}




######################################
######################################
##                                  ##
##  PART 10                         ##
##  Solomons predicting Brazil      ## 
##                                  ##                                  
######################################
######################################


AB_ENSEMBLE_sol_braz = list()

for(i in 1:7)
{
	AB_ENSEMBLE_sol_braz[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_sol_braz[[i]]))
	{
		AB_ENSEMBLE_sol_braz[[i]][[count]] = AB_LOGIT_sol_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_sol_braz[[i]]))
	{
		AB_ENSEMBLE_sol_braz[[i]][[count]] = AB_LDA_sol_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_sol_braz[[i]]))
	{
		AB_ENSEMBLE_sol_braz[[i]][[count]] = AB_QDA_sol_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_sol_braz[[i]]))
	{
		AB_ENSEMBLE_sol_braz[[i]][[count]] = AB_Dtree_sol_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_sol_braz[[i]]))
	{
		AB_ENSEMBLE_sol_braz[[i]][[count]] = AB_Rforest_sol_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_sol_braz[[i]]))
	{
		AB_ENSEMBLE_sol_braz[[i]][[count]] = AB_DYN1_sol_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_sol_braz[[i]]))
	{
		AB_ENSEMBLE_sol_braz[[i]][[count]] = AB_DYN2_sol_braz[[i]][[j]][[1]]

		count = count + 1
	}	
}


AB_ENSEMBLE_sol_braz_chull = list()

AB_ENSEMBLE_sol_braz_chull[[1]] = AB_braz_1

XX_VEC = AB_braz_1[,1]
YY_VEC = AB_braz_1[,2]


for(i in 1:length(AB_ENSEMBLE_sol_braz))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_sol_braz[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_sol_braz[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_sol_braz[[i]][[j]][,1] )
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


	AB_ENSEMBLE_sol_braz_chull[[1+i]] = XY_hull_plot 
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
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_sol_sol[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_sol_sol[[i]][[j]][,1] )
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
##  PART 12                         ##
##  Solomons predicting All         ## 
##                                  ##                                  
######################################
######################################


AB_ENSEMBLE_sol_all = list()

for(i in 1:7)
{
	AB_ENSEMBLE_sol_all[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_sol_all[[i]]))
	{
		AB_ENSEMBLE_sol_all[[i]][[count]] = AB_LOGIT_sol_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_sol_all[[i]]))
	{
		AB_ENSEMBLE_sol_all[[i]][[count]] = AB_LDA_sol_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_sol_all[[i]]))
	{
		AB_ENSEMBLE_sol_all[[i]][[count]] = AB_QDA_sol_all[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_sol_all[[i]]))
	{
		AB_ENSEMBLE_sol_all[[i]][[count]] = AB_Dtree_sol_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_sol_all[[i]]))
	{
		AB_ENSEMBLE_sol_all[[i]][[count]] = AB_Rforest_sol_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_sol_all[[i]]))
	{
		AB_ENSEMBLE_sol_all[[i]][[count]] = AB_DYN1_sol_all[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_sol_all[[i]]))
	{
		AB_ENSEMBLE_sol_all[[i]][[count]] = AB_DYN2_sol_all[[i]][[j]]

		count = count + 1
	}	
}


AB_ENSEMBLE_sol_all_chull = list()

AB_ENSEMBLE_sol_all_chull[[1]] = AB_all_1

XX_VEC = AB_all_1[,1]
YY_VEC = AB_all_1[,2]

for(i in 1:length(AB_ENSEMBLE_sol_all))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_sol_all[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_sol_all[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_sol_all[[i]][[j]][,1] )
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


	AB_ENSEMBLE_sol_all_chull[[1+i]] = XY_hull_plot 
}




######################################
######################################
##                                  ##
##  PART 13                         ##
##  All predicting Thailand         ## 
##                                  ##                                  
######################################
######################################


AB_ENSEMBLE_all_thai = list()

for(i in 1:7)
{
	AB_ENSEMBLE_all_thai[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_all_thai[[i]]))
	{
		AB_ENSEMBLE_all_thai[[i]][[count]] = AB_LOGIT_all_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_all_thai[[i]]))
	{
		AB_ENSEMBLE_all_thai[[i]][[count]] = AB_LDA_all_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_all_thai[[i]]))
	{
		AB_ENSEMBLE_all_thai[[i]][[count]] = AB_QDA_all_thai[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_all_thai[[i]]))
	{
		AB_ENSEMBLE_all_thai[[i]][[count]] = AB_Dtree_all_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_all_thai[[i]]))
	{
		AB_ENSEMBLE_all_thai[[i]][[count]] = AB_Rforest_all_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_all_thai[[i]]))
	{
		AB_ENSEMBLE_all_thai[[i]][[count]] = AB_DYN1_all_thai[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_all_thai[[i]]))
	{
		AB_ENSEMBLE_all_thai[[i]][[count]] = AB_DYN2_all_thai[[i]][[j]]

		count = count + 1
	}	
}


AB_ENSEMBLE_all_thai_chull = list()

AB_ENSEMBLE_all_thai_chull[[1]] = AB_thai_1

XX_VEC = AB_thai_1[,1]
YY_VEC = AB_thai_1[,2]

for(i in 1:length(AB_ENSEMBLE_all_thai))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_all_thai[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_all_thai[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_all_thai[[i]][[j]][,1] )
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

	AB_ENSEMBLE_all_thai_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 14                         ##
##  All predicting Brazil           ## 
##                                  ##                                  
######################################
######################################


AB_ENSEMBLE_all_braz = list()

for(i in 1:7)
{
	AB_ENSEMBLE_all_braz[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_all_braz[[i]]))
	{
		AB_ENSEMBLE_all_braz[[i]][[count]] = AB_LOGIT_all_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_all_braz[[i]]))
	{
		AB_ENSEMBLE_all_braz[[i]][[count]] = AB_LDA_all_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_all_braz[[i]]))
	{
		AB_ENSEMBLE_all_braz[[i]][[count]] = AB_QDA_all_braz[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_all_braz[[i]]))
	{
		AB_ENSEMBLE_all_braz[[i]][[count]] = AB_Dtree_all_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_all_braz[[i]]))
	{
		AB_ENSEMBLE_all_braz[[i]][[count]] = AB_Rforest_all_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_all_braz[[i]]))
	{
		AB_ENSEMBLE_all_braz[[i]][[count]] = AB_DYN1_all_braz[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_all_braz[[i]]))
	{
		AB_ENSEMBLE_all_braz[[i]][[count]] = AB_DYN2_all_braz[[i]][[j]]

		count = count + 1
	}	
}


AB_ENSEMBLE_all_braz_chull = list()

AB_ENSEMBLE_all_braz_chull[[1]] = AB_braz_1

XX_VEC = AB_braz_1[,1]
YY_VEC = AB_braz_1[,2]

for(i in 1:length(AB_ENSEMBLE_all_braz))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_all_braz[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_all_braz[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_all_braz[[i]][[j]][,1] )
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


	AB_ENSEMBLE_all_braz_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 15                         ##
##  All predicting Solomons         ## 
##                                  ##                                  
######################################
######################################

AB_ENSEMBLE_all_sol = list()

for(i in 1:7)
{
	AB_ENSEMBLE_all_sol[[i]] = list()

	count  = 1

	for(j in 1:length(AB_LOGIT_all_sol[[i]]))
	{
		AB_ENSEMBLE_all_sol[[i]][[count]] = AB_LOGIT_all_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_LDA_all_sol[[i]]))
	{
		AB_ENSEMBLE_all_sol[[i]][[count]] = AB_LDA_all_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_QDA_all_sol[[i]]))
	{
		AB_ENSEMBLE_all_sol[[i]][[count]] = AB_QDA_all_sol[[i]][[j]]

		count = count + 1
	}

	for(j in 1:length(AB_Dtree_all_sol[[i]]))
	{
		AB_ENSEMBLE_all_sol[[i]][[count]] = AB_Dtree_all_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_Rforest_all_sol[[i]]))
	{
		AB_ENSEMBLE_all_sol[[i]][[count]] = AB_Rforest_all_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN1_all_sol[[i]]))
	{
		AB_ENSEMBLE_all_sol[[i]][[count]] = AB_DYN1_all_sol[[i]][[j]]

		count = count + 1
	}	

	for(j in 1:length(AB_DYN2_all_sol[[i]]))
	{
		AB_ENSEMBLE_all_sol[[i]][[count]] = AB_DYN2_all_sol[[i]][[j]]

		count = count + 1
	}	
}


AB_ENSEMBLE_all_sol_chull = list()

AB_ENSEMBLE_all_sol_chull[[1]] = AB_sol_1

XX_VEC = AB_sol_1[,1]
YY_VEC = AB_sol_1[,2]

for(i in 1:length(AB_ENSEMBLE_all_sol))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_ENSEMBLE_all_sol[[i]]))
	{
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_all_sol[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_all_sol[[i]][[j]][,1] )
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


	AB_ENSEMBLE_all_sol_chull[[1+i]] = XY_hull_plot 
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
		XX_VEC = c( XX_VEC, 1 - AB_ENSEMBLE_all_all[[i]][[j]][,2] )
		YY_VEC = c( YY_VEC, AB_ENSEMBLE_all_all[[i]][[j]][,1] )
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



tiff(file="FigSX_AB_ENSEMBLE_ROC_curves.tif", width=32, height=26, units="cm", res=500)


lay.mat <- rbind( c( 1, 2, 3, 4 ),
                  c( 5, 6, 7, 8 ), 
                  c( 9,10,11,12 ),
                  c(13,14,15,16 ),
                  c(17,17,17,17 ) )
layout(lay.mat, heights=c(10,10,10,10,1))
layout.show(17)


par(mar=c(2.5,2.5,2,1))
par(mgp=c(1.5, 0.55,0))


line_seq <- c(0.2, 0.4, 0.6, 0.8)


line.size = 2
lab.size  = 1.25
axis.size = 1
main.size = 1.25


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
main="training = Thailand; testing = Thailand",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_thai_thai_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_thai_thai_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 



######################################
######################################
##                                  ##
##  PANEL 2                         ##
##  Thailand predicting Brazil      ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Thailand; testing = Brazil",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_thai_braz_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_thai_braz_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 






######################################
######################################
##                                  ##
##  PANEL 3                         ##
##  Thailand predicting Solomons    ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Thailand; testing = Solomons",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_thai_sol_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_thai_sol_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 




######################################
######################################
##                                  ##
##  PANEL 4                         ##
##  Thailand predicting All         ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Thailand; testing = all regions",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_thai_all_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_thai_all_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 



######################################
######################################
##                                  ##
##  PANEL 5                         ##
##  Brazil predicting Thailand      ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Brazil; testing = Thailand",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_braz_thai_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_braz_thai_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 





######################################
######################################
##                                  ##
##  PANEL 6                         ##
##  Brazil predicting Brazil        ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Brazil; testing = Brazil",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_braz_braz_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_braz_braz_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 




######################################
######################################
##                                  ##
##  PANEL 7                         ##
##  Brazil predicting Solomons      ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Brazil; testing = Solomons",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_braz_sol_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_braz_sol_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 



######################################
######################################
##                                  ##
##  PANEL 8                         ##
##  Brazil predicting All           ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Brazil; testing = all regions",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_braz_all_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_braz_all_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 





######################################
######################################
##                                  ##
##  PANEL 9                         ##
##  Solomons predicting Thailand    ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Solomons; testing = Thailand",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_sol_thai_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_sol_thai_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 






######################################
######################################
##                                  ##
##  PANEL 10                        ##
##  Solomons predicting Brazil      ## 
##                                  ##                                  
######################################
######################################

plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Solomons; testing = Brazil",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_sol_braz_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_sol_braz_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 






######################################
######################################
##                                  ##
##  PANEL 11                        ##
##  Solomons predicting Solomons    ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Solomons; testing = Solomons",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_sol_sol_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_sol_sol_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 



######################################
######################################
##                                  ##
##  PANEL 12                        ##
##  Solomons predicting All         ## 
##                                  ##                                  
######################################
######################################

plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = Solomons; testing = all regions",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_sol_all_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_sol_all_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 




######################################
######################################
##                                  ##
##  PANEL 13                        ##
##  All predicting Thailand         ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = all regions; testing = Thailand",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_all_thai_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_all_thai_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 





######################################
######################################
##                                  ##
##  PANEL 14                        ##
##  All predicting Brazil           ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = all regions; testing = Brazil",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_all_braz_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_all_braz_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 






######################################
######################################
##                                  ##
##  PANEL 15                        ##
##  All predicting Solomons         ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = all regions; testing = Solomons",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_all_sol_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_all_sol_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}
axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 




######################################
######################################
##                                  ##
##  PANEL 16                        ##
##  All predicting All              ## 
##                                  ##                                  
######################################
######################################


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="training = all regions; testing = all regions",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(i in 8:1)
{
	points( x=1-AB_ENSEMBLE_all_all_chull[[i]][,2], 
      	  y=1-AB_ENSEMBLE_all_all_chull[[i]][,1], 
      	  type='s', lwd=line.size, col=ant_cols[i])
}

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 ) 




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
       legend = c("1 antigen", "2 antigens", "3 antigens", "4 antigens", "5 antigens", "6 antigens", "7 antigens", "8 antigens"), 
       fill = ant_cols, 
       border = ant_cols, 
       ncol=8, cex=2, bty="n" )

dev.off()




