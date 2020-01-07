library(MASS)
library(ROCR)
library(randomForest)


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

thailand_data = read.csv("C:\\U\\GHIT\\NatMed_Paper\\Data\\proc\\thailand_ab_epi_data.csv")

brazil_data   = read.csv("C:\\U\\GHIT\\NatMed_Paper\\Data\\proc\\brazil_ab_epi_data.csv")

solomon_data  = read.csv("C:\\U\\GHIT\\NatMed_Paper\\Data\\proc\\solomon_ab_epi_data.csv")

control_data  = read.csv("C:\\U\\GHIT\\NatMed_Paper\\Data\\proc\\control_ab_epi_data.csv")


###################################
###################################
##                               ##
##  CATEGORISATION               ##
##                               ##
###################################
###################################

###################################
## Controls

control_cat = rep("VBDR", nrow(control_data) )

for(i in 1:length(control_cat))
{
	if( substr( control_data[i,1], 1, 3 ) == "TRC" )
	{
		control_cat[i] = "TRC"
	}
	if( substr( control_data[i,1], 1, 3 ) == "ARC" )
	{
		control_cat[i] = "ARC"
	}
	if( substr( control_data[i,1], 1, 3 ) == "BRC" )
	{
		control_cat[i] = "BRC"
	}
}

###################################
###################################
## Thailand

thailand_cat = rep("thai_never", nrow(thailand_data) )

thailand_cat[which( thailand_data[,10] == 0 )] <- "thai_current"

thailand_cat[intersect( which(thailand_data[,10]>0), which(thailand_data[,10]<=9*30) )] <- "thai_recent"

thailand_cat[which(thailand_data[,10]>9*30)] <- "thai_old"


####################################
## Put together Thai and control data

AB_thai = rbind( thailand_data[,17:81], control_data[,17:81] )

AB_thai = log(AB_thai)

N_part_thai = nrow(AB_thai)


inf_cat_thai = c( thailand_cat, control_cat )

bin_cat_thai = rep("old", N_part_thai)
bin_cat_thai[which(inf_cat_thai=="thai_current")] = "new"
bin_cat_thai[which(inf_cat_thai=="thai_recent")]  = "new"



###################################
###################################
## Brazil

brazil_cat = rep("braz_never", nrow(brazil_data) )

brazil_cat[which( brazil_data[,10] == 0 )] <- "braz_current"

brazil_cat[intersect( which(brazil_data[,10]>0), which(brazil_data[,10]<=9*30) )] <- "braz_recent"

brazil_cat[which(brazil_data[,10]>9*30)] <- "braz_old"


####################################
## Put together Brazilian and control data

AB_braz = rbind( brazil_data[,17:81], control_data[,17:81] )

AB_braz = log(AB_braz)

N_part_braz = nrow(AB_braz)


inf_cat_braz = c( brazil_cat, control_cat )

bin_cat_braz = rep("old", N_part_braz)
bin_cat_braz[which(inf_cat_braz=="braz_current")] = "new"
bin_cat_braz[which(inf_cat_braz=="braz_recent")]  = "new"


###################################
###################################
## Solomon Islands

solomon_cat = rep("sol_never", nrow(solomon_data) )

solomon_cat[which( solomon_data[,10] == 0 )] <- "sol_current"

solomon_cat[intersect( which(solomon_data[,10]>0), which(solomon_data[,10]<=9*30) )] <- "sol_recent"

solomon_cat[which(solomon_data[,10]>9*30)] <- "sol_old"


####################################
## Put together Solomon Islands and control data

AB_sol = rbind( solomon_data[,17:81], control_data[,17:81] )

AB_sol = log(AB_sol)

N_part_sol = nrow(AB_sol)


inf_cat_sol = c( solomon_cat, control_cat )

bin_cat_sol = rep("old", N_part_sol)
bin_cat_sol[which(inf_cat_sol=="sol_current")] = "new"
bin_cat_sol[which(inf_cat_sol=="sol_recent")]  = "new"


###################################
###################################
## All data

AB_all = rbind( thailand_data[,17:81], brazil_data[,17:81], solomon_data[,17:81], control_data[,17:81] )

AB_all = log(AB_all)

N_part_all = nrow(AB_all)


inf_cat_all = c( thailand_cat, brazil_cat, solomon_cat, control_cat )

bin_cat_all = rep("old", N_part_all)
bin_cat_all[which(inf_cat_all%in%c("thai_current", "braz_current", "sol_current"))] = "new"
bin_cat_all[which(inf_cat_all%in%c("thai_recent", "braz_recent", "sol_recent"))]    = "new"



###################################
###################################
## Trim out antibodies with too much 
## missing data


ant_drop <- c()

for(j in 1:ncol(AB_all))
{
	if( length(which(is.na(AB_all[,j]))) > 250 )
	{
		ant_drop = c(ant_drop, j)
	}
}




AB_all  = AB_all[,-ant_drop]
AB_thai = AB_thai[,-ant_drop]
AB_braz = AB_braz[,-ant_drop]
AB_sol  = AB_sol[,-ant_drop]


N_ant     = ncol(AB_all)
ant_names = colnames(AB_all)




#############################################################################
#############################################################################
##                                                                         ## 
##   ####  ##     ####   ####   ####  #### ##### #### ##### #####   ####   ## 
##  ##  ## ##    ##  ## ##     ##      ##  ##     ##  ##    ##  ## ##      ##
##  ##     ##    ######  ####   ####   ##  ####   ##  ####  #####   ####   ##
##  ##  ## ##    ##  ##     ##     ##  ##  ##     ##  ##    ## ##      ##  ##
##   ####  ##### ##  ##  ####   ####  #### ##    #### ##### ##  ##  ####   ## 
##                                                                         ##
#############################################################################
#############################################################################

top_8 = c(50, 1, 47, 2, 30, 8, 58, 39)



N_rep = 1000

N_plot = 100

N_prob = 500

prob_new_cut = seq(from=0, to=1, length=N_prob)



AB_QDA = function( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list )
{
	AB_train = as.matrix(AB_train[,antigen_list])

	if( length(which(is.na(AB_train))) > 0 )
	{
		bin_cat_train = bin_cat_train[-which(is.na(AB_train), arr.ind=TRUE)[,1]]

		AB_train = AB_train[-which(is.na(AB_train), arr.ind=TRUE)[,1],]
	}


	AB_test = as.matrix(AB_test[,antigen_list])

	if( length(which(is.na(AB_test))) > 0 )
	{
		bin_cat_test = bin_cat_test[-which(is.na(AB_test), arr.ind=TRUE)[,1]]

		AB_test = AB_test[-which(is.na(AB_test), arr.ind=TRUE)[,1],]
	}

	############################################
	## Fit QDA and create prediction object

	tryCatch(
	{
		QDA_2 = qda(AB_train, grouping=bin_cat_train)

		QDA_2_pred_obj = predict( QDA_2, AB_test)

		QDA_2_pred = prediction(QDA_2_pred_obj$posterior[,2], bin_cat_test)

		QDA_2_perf_xy = performance(QDA_2_pred, "spec", "sens")

	}, error=function(e){ NULL }
	)

	

	SS_mat = cbind( rev(1 - QDA_2_perf_xy@x.values[[1]]), rev(QDA_2_perf_xy@y.values[[1]]) )
	colnames(SS_mat) = c("x_plot", "y_plot")

	SS_mat
}


AB_QDA_Xval = function( AB_data, bin_cat_data, antigen_list, train_prop, N_rep )
{
	x_vals = c()	
	y_vals = c()

	for(n in 1:N_rep)
	{
		############################################
		## Prepare testing and training data

		index_train = sample( nrow(AB_data), train_prop*nrow(AB_data) )
		index_test  = setdiff( 1:nrow(AB_data), index_train )

		bin_cat_train = bin_cat_data[index_train]
		bin_cat_test  = bin_cat_data[index_test]

		AB_train = AB_data[index_train,]
		AB_test  = AB_data[index_test,]

		
		AB_class = AB_QDA( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list )

		x_vals = c( x_vals, AB_class[,1] )
		y_vals = c( y_vals, AB_class[,2] )
	}


	lowess_xy = lowess( x_vals, y_vals, f=0.2 )

	SS_plot = matrix(NA, nrow=N_plot, ncol=2)
	colnames(SS_plot) = c("x_plot", "y_plot")

	SS_plot[,1] = lowess_xy$x[seq(from=1, to=length(lowess_xy$x), length=N_plot)]
	SS_plot[,2] = lowess_xy$y[seq(from=1, to=length(lowess_xy$y), length=N_plot)]

	SS_plot = rbind( c(0,0), SS_plot, c(1,1) )

	SS_plot
}



AB_QDA_Xval_train_subset = function( AB_data, bin_cat_data, inf_cat_data, train_subset, 
                                              antigen_list, train_prop, N_rep )
{
	x_vals = c()	
	y_vals = c()

	for(n in 1:N_rep)
	{
		############################################
		## Prepare testing and training data

		index = which( inf_cat_data %in% train_subset )

		index_train = sample( index, train_prop*length(index) )
		index_test  = setdiff( 1:nrow(AB_data), index_train )

		bin_cat_train = bin_cat_data[index_train]
		bin_cat_test  = bin_cat_data[index_test]

		AB_train = AB_data[index_train,]
		AB_test  = AB_data[index_test,]

		
		AB_class = AB_QDA( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list )

		x_vals = c( x_vals, AB_class[,1] )
		y_vals = c( y_vals, AB_class[,2] )
	}


	lowess_xy = lowess( x_vals, y_vals, f=0.2 )

	SS_plot = matrix(NA, nrow=N_plot, ncol=2)
	colnames(SS_plot) = c("x_plot", "y_plot")

	SS_plot[,1] = lowess_xy$x[seq(from=1, to=length(lowess_xy$x), length=N_plot)]
	SS_plot[,2] = lowess_xy$y[seq(from=1, to=length(lowess_xy$y), length=N_plot)]

	SS_plot = rbind( c(0,0), SS_plot, c(1,1) )

	SS_plot
}



AB_QDA_Xval_test_subset = function( AB_data, bin_cat_data, inf_cat_data, test_subset, 
                                             antigen_list, train_prop, N_rep )
{
	x_vals = c()	
	y_vals = c()

	for(n in 1:N_rep)
	{
		############################################
		## Prepare testing and training data

		index = which( inf_cat_data %in% test_subset )
	
		index_test  = sample( index, (1-train_prop)*length(index) )
		index_train = setdiff( 1:nrow(AB_data), index_test )

		bin_cat_train = bin_cat_data[index_train]
		bin_cat_test  = bin_cat_data[index_test]

		AB_train = AB_data[index_train,]
		AB_test  = AB_data[index_test,]

		
		AB_class = AB_QDA( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list )

		x_vals = c( x_vals, AB_class[,1] )
		y_vals = c( y_vals, AB_class[,2] )
	}


	lowess_xy = lowess( x_vals, y_vals, f=0.2 )

	SS_plot = matrix(NA, nrow=N_plot, ncol=2)
	colnames(SS_plot) = c("x_plot", "y_plot")

	SS_plot[,1] = lowess_xy$x[seq(from=1, to=length(lowess_xy$x), length=N_plot)]
	SS_plot[,2] = lowess_xy$y[seq(from=1, to=length(lowess_xy$y), length=N_plot)]

	SS_plot = rbind( c(0,0), SS_plot, c(1,1) )

	SS_plot
}



#######################################################
#######################################################
##                                                   ## 
##   ####  ##  ## #####   ####  ##### ######  ####   ## 
##  ##     ##  ## ##  ## ##     ##      ##   ##      ##
##   ####  ##  ## #####   ####  ####    ##    ####   ##
##      ## ##  ## ##  ##     ## ##      ##       ##  ##
##   ####   ####  #####   ####  #####   ##    ####   ## 
##                                                   ##
#######################################################
#######################################################


AB_QDA_allant = function( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list )
{
	OUTPUT <- list()


	######################
	##                  ##
	## 2 antigens       ##
	##                  ##
	######################

	AB_QDA_2 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			if( i > j )
			{
				AB_QDA_2[[count]] = AB_QDA( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list[c(i,j)] )
	
				count = count + 1
			}
		}
	}

	OUTPUT[[1]] = AB_QDA_2


	######################
	##                  ##
	## 3 antigens       ##
	##                  ##
	######################

	AB_QDA_3 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				if( (i > j) && (j > k) )
				{
					AB_QDA_3[[count]] = AB_QDA( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list[c(i,j,k)] )
	
					count = count + 1
				}
			}
		}
	}

	OUTPUT[[2]] = AB_QDA_3


	######################
	##                  ##
	## 4 antigens       ##
	##                  ##
	######################

	AB_QDA_4 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					if( (i > j) && (j > k) && (k > l) )
					{
						AB_QDA_4[[count]] = AB_QDA( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list[c(i,j,k,l)] )
	
						count = count + 1
					}
				}
			}
		}
	}

	OUTPUT[[3]] = AB_QDA_4


	######################
	##                  ##
	## 5 antigens       ##
	##                  ##
	######################

	AB_QDA_5 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						if( (i > j) && (j > k) && (k > l) && (l > m) )
						{
							AB_QDA_5[[count]] = AB_QDA( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list[c(i,j,k,l,m)] )
	
							count = count + 1
						}
					}
				}
			}
		}
	}

	OUTPUT[[4]] = AB_QDA_5


	######################
	##                  ##
	## 6 antigens       ##
	##                  ##
	######################

	AB_QDA_6 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						for(n in 1:8)
						{
							if( (i > j) && (j > k) && (k > l) && (l > m) && (m > n) )
							{
								AB_QDA_6[[count]] = AB_QDA( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list[c(i,j,k,l,m,n)] )
	
								count = count + 1
							}
						}
					}
				}
			}
		}
	}

	OUTPUT[[5]] = AB_QDA_6


	######################
	##                  ##
	## 7 antigens       ##
	##                  ##
	######################

	AB_QDA_7 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						for(n in 1:8)
						{
							for(o in 1:8)
							{
								if( (i > j) && (j > k) && (k > l) && (l > m) && (m > n) && (n > o) )
								{
									AB_QDA_7[[count]] = AB_QDA( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list[c(i,j,k,l,m,n,o)] )
		
									count = count + 1
								}
							}
						}
					}
				}
			}
		}
	}

	OUTPUT[[6]] = AB_QDA_7



	######################
	##                  ##
	## 8 antigens       ##
	##                  ##
	######################

	AB_QDA_8 <- list()

	AB_QDA_8[[1]] = AB_QDA( AB_train, bin_cat_train, AB_test, bin_cat_test, antigen_list[1:8] )

	OUTPUT[[7]] = AB_QDA_8



	######################
	##                  ##
	## Output           ##
	##                  ##
	######################
	
	OUTPUT
}





AB_QDA_Xval_allant = function( AB_data, bin_cat_data, antigen_list, train_prop, N_rep )
{
	OUTPUT <- list()


	######################
	##                  ##
	## 2 antigens       ##
	##                  ##
	######################

	AB_QDA_2 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			if( i > j )
			{
				AB_QDA_2[[count]] = AB_QDA_Xval( AB_data, bin_cat_data, antigen_list[c(i,j)], train_prop, N_rep )
	
				count = count + 1
			}
		}
	}

	OUTPUT[[1]] = AB_QDA_2


	######################
	##                  ##
	## 3 antigens       ##
	##                  ##
	######################

	AB_QDA_3 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				if( (i > j) && (j > k) )
				{
					AB_QDA_3[[count]] = AB_QDA_Xval( AB_data, bin_cat_data, antigen_list[c(i,j,k)], train_prop, N_rep )
	
					count = count + 1
				}
			}
		}
	}

	OUTPUT[[2]] = AB_QDA_3


	######################
	##                  ##
	## 4 antigens       ##
	##                  ##
	######################

	AB_QDA_4 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					if( (i > j) && (j > k) && (k > l) )
					{
						AB_QDA_4[[count]] = AB_QDA_Xval( AB_data, bin_cat_data, antigen_list[c(i,j,k,l)], train_prop, N_rep )
	
						count = count + 1
					}
				}
			}
		}
	}

	OUTPUT[[3]] = AB_QDA_4


	######################
	##                  ##
	## 5 antigens       ##
	##                  ##
	######################

	AB_QDA_5 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						if( (i > j) && (j > k) && (k > l) && (l > m) )
						{
							AB_QDA_5[[count]] = AB_QDA_Xval( AB_data, bin_cat_data, antigen_list[c(i,j,k,l,m)], train_prop, N_rep )
	
							count = count + 1
						}
					}
				}
			}
		}
	}

	OUTPUT[[4]] = AB_QDA_5


	######################
	##                  ##
	## 6 antigens       ##
	##                  ##
	######################

	AB_QDA_6 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						for(n in 1:8)
						{
							if( (i > j) && (j > k) && (k > l) && (l > m) && (m > n) )
							{
								AB_QDA_6[[count]] = AB_QDA_Xval( AB_data, bin_cat_data, antigen_list[c(i,j,k,l,m,n)], train_prop, N_rep )
	
								count = count + 1
							}
						}
					}
				}
			}
		}
	}

	OUTPUT[[5]] = AB_QDA_6


	######################
	##                  ##
	## 7 antigens       ##
	##                  ##
	######################

	AB_QDA_7 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						for(n in 1:8)
						{
							for(o in 1:8)
							{
								if( (i > j) && (j > k) && (k > l) && (l > m) && (m > n) && (n > o) )
								{
									AB_QDA_7[[count]] = AB_QDA_Xval( AB_data, bin_cat_data, antigen_list[c(i,j,k,l,m,n,o)], train_prop, N_rep )
		
									count = count + 1
								}
							}
						}
					}
				}
			}
		}
	}

	OUTPUT[[6]] = AB_QDA_7



	######################
	##                  ##
	## 8 antigens       ##
	##                  ##
	######################

	AB_QDA_8 <- list()

	AB_QDA_8[[1]] = AB_QDA_Xval( AB_data, bin_cat_data, antigen_list[1:8], train_prop, N_rep )

	OUTPUT[[7]] = AB_QDA_8



	######################
	##                  ##
	## Output           ##
	##                  ##
	######################
	
	OUTPUT
}





AB_QDA_Xval_train_subset_allant = function( AB_data, bin_cat_data, inf_cat_data, train_subset, 
                                             antigen_list, train_prop, N_rep )
{
	OUTPUT <- list()


	######################
	##                  ##
	## 2 antigens       ##
	##                  ##
	######################

	AB_QDA_2 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			if( i > j )
			{
				AB_QDA_2[[count]] = AB_QDA_Xval_train_subset( AB_data, bin_cat_data, inf_cat_data, train_subset, 
                                                                        antigen_list[c(i,j)], train_prop, N_rep )

				count = count + 1
			}
		}
	}

	OUTPUT[[1]] = AB_QDA_2


	######################
	##                  ##
	## 3 antigens       ##
	##                  ##
	######################

	AB_QDA_3 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				if( (i > j) && (j > k) )
				{
					AB_QDA_3[[count]] = AB_QDA_Xval_train_subset( AB_data, bin_cat_data, inf_cat_data, train_subset, 
                                                                              antigen_list[c(i,j,k)], train_prop, N_rep )
	
					count = count + 1
				}
			}
		}
	}

	OUTPUT[[2]] = AB_QDA_3


	######################
	##                  ##
	## 4 antigens       ##
	##                  ##
	######################

	AB_QDA_4 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					if( (i > j) && (j > k) && (k > l) )
					{
						AB_QDA_4[[count]] = AB_QDA_Xval_train_subset( AB_data, bin_cat_data, inf_cat_data, train_subset, 
                                                                                    antigen_list[c(i,j,k,l)], train_prop, N_rep )
	
						count = count + 1
					}
				}
			}
		}
	}

	OUTPUT[[3]] = AB_QDA_4


	######################
	##                  ##
	## 5 antigens       ##
	##                  ##
	######################

	AB_QDA_5 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						if( (i > j) && (j > k) && (k > l) && (l > m) )
						{
							AB_QDA_5[[count]] = AB_QDA_Xval_train_subset( AB_data, bin_cat_data, inf_cat_data, train_subset, 
                                                                                          antigen_list[c(i,j,k,l,m)], train_prop, N_rep )
	
							count = count + 1
						}
					}
				}
			}
		}
	}

	OUTPUT[[4]] = AB_QDA_5


	######################
	##                  ##
	## 6 antigens       ##
	##                  ##
	######################

	AB_QDA_6 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						for(n in 1:8)
						{
							if( (i > j) && (j > k) && (k > l) && (l > m) && (m > n) )
							{
								AB_QDA_6[[count]] = AB_QDA_Xval_train_subset( AB_data, bin_cat_data, inf_cat_data, train_subset, 
                                                                                                antigen_list[c(i,j,k,l,m,n)], train_prop, N_rep )
	
								count = count + 1
							}
						}
					}
				}
			}
		}
	}

	OUTPUT[[5]] = AB_QDA_6


	######################
	##                  ##
	## 7 antigens       ##
	##                  ##
	######################

	AB_QDA_7 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						for(n in 1:8)
						{
							for(o in 1:8)
							{
								if( (i > j) && (j > k) && (k > l) && (l > m) && (m > n) && (n > o) )
								{
									AB_QDA_7[[count]] = AB_QDA_Xval_train_subset( AB_data, bin_cat_data, inf_cat_data, train_subset, 
                                                                                                      antigen_list[c(i,j,k,l,m,n,o)], train_prop, N_rep )
		
									count = count + 1
								}
							}
						}
					}
				}
			}
		}
	}

	OUTPUT[[6]] = AB_QDA_7



	######################
	##                  ##
	## 8 antigens       ##
	##                  ##
	######################

	AB_QDA_8 <- list()

	AB_QDA_8[[1]] = AB_QDA_Xval_train_subset( AB_data, bin_cat_data, inf_cat_data, train_subset, 
                                                  antigen_list[1:8], train_prop, N_rep )

	OUTPUT[[7]] = AB_QDA_8



	######################
	##                  ##
	## Output           ##
	##                  ##
	######################
	
	OUTPUT
}









AB_QDA_Xval_test_subset_allant = function( AB_data, bin_cat_data, inf_cat_data, test_subset, 
                                             antigen_list, train_prop, N_rep )
{
	OUTPUT <- list()


	######################
	##                  ##
	## 2 antigens       ##
	##                  ##
	######################

	AB_QDA_2 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			if( i > j )
			{
				AB_QDA_2[[count]] = AB_QDA_Xval_test_subset( AB_data, bin_cat_data, inf_cat_data, test_subset, 
                                                                        antigen_list[c(i,j)], train_prop, N_rep )

				count = count + 1
			}
		}
	}

	OUTPUT[[1]] = AB_QDA_2


	######################
	##                  ##
	## 3 antigens       ##
	##                  ##
	######################

	AB_QDA_3 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				if( (i > j) && (j > k) )
				{
					AB_QDA_3[[count]] = AB_QDA_Xval_test_subset( AB_data, bin_cat_data, inf_cat_data, test_subset, 
                                                                              antigen_list[c(i,j,k)], train_prop, N_rep )
	
					count = count + 1
				}
			}
		}
	}

	OUTPUT[[2]] = AB_QDA_3


	######################
	##                  ##
	## 4 antigens       ##
	##                  ##
	######################

	AB_QDA_4 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					if( (i > j) && (j > k) && (k > l) )
					{
						AB_QDA_4[[count]] = AB_QDA_Xval_test_subset( AB_data, bin_cat_data, inf_cat_data, test_subset, 
                                                                                    antigen_list[c(i,j,k,l)], train_prop, N_rep )
	
						count = count + 1
					}
				}
			}
		}
	}

	OUTPUT[[3]] = AB_QDA_4


	######################
	##                  ##
	## 5 antigens       ##
	##                  ##
	######################

	AB_QDA_5 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						if( (i > j) && (j > k) && (k > l) && (l > m) )
						{
							AB_QDA_5[[count]] = AB_QDA_Xval_test_subset( AB_data, bin_cat_data, inf_cat_data, test_subset, 
                                                                                          antigen_list[c(i,j,k,l,m)], train_prop, N_rep )
	
							count = count + 1
						}
					}
				}
			}
		}
	}

	OUTPUT[[4]] = AB_QDA_5


	######################
	##                  ##
	## 6 antigens       ##
	##                  ##
	######################

	AB_QDA_6 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						for(n in 1:8)
						{
							if( (i > j) && (j > k) && (k > l) && (l > m) && (m > n) )
							{
								AB_QDA_6[[count]] = AB_QDA_Xval_test_subset( AB_data, bin_cat_data, inf_cat_data, test_subset, 
                                                                                                antigen_list[c(i,j,k,l,m,n)], train_prop, N_rep )
	
								count = count + 1
							}
						}
					}
				}
			}
		}
	}

	OUTPUT[[5]] = AB_QDA_6


	######################
	##                  ##
	## 7 antigens       ##
	##                  ##
	######################

	AB_QDA_7 <- list()

	count = 1

	for(i in 1:8)
	{
		for(j in 1:8)
		{
			for(k in 1:8)
			{
				for(l in 1:8)
				{
					for(m in 1:8)
					{
						for(n in 1:8)
						{
							for(o in 1:8)
							{
								if( (i > j) && (j > k) && (k > l) && (l > m) && (m > n) && (n > o) )
								{
									AB_QDA_7[[count]] = AB_QDA_Xval_test_subset( AB_data, bin_cat_data, inf_cat_data, test_subset, 
                                                                                                      antigen_list[c(i,j,k,l,m,n,o)], train_prop, N_rep )
		
									count = count + 1
								}
							}
						}
					}
				}
			}
		}
	}

	OUTPUT[[6]] = AB_QDA_7



	######################
	##                  ##
	## 8 antigens       ##
	##                  ##
	######################

	AB_QDA_8 <- list()

	AB_QDA_8[[1]] = AB_QDA_Xval_test_subset( AB_data, bin_cat_data, inf_cat_data, test_subset, 
                                                  antigen_list[1:8], train_prop, N_rep )

	OUTPUT[[7]] = AB_QDA_8



	######################
	##                  ##
	## Output           ##
	##                  ##
	######################
	
	OUTPUT
}









############################
############################
##                        ##
##  #####  ##  ## #   ##  ##
##  ##  ## ##  ## ##  ##  ##
##  #####  ##  ## ### ##  ## 
##  ## ##  ##  ## ## ###  ##
##  ##  ##  ####  ##  ##  ##
##                        ##
############################ 
############################


######################################
######################################
##                                  ##
##  PART 1                          ##
##  Thailand predicting Thailand    ## 
##                                  ##                                  
######################################
######################################



AB_QDA_thai_thai = AB_QDA_Xval_allant( AB_thai, bin_cat_thai, top_8, 0.666, N_rep )



######################################
######################################
##                                  ##
##  PART 2                          ##
##  Thailand predicting Brazil      ## 
##                                  ##                                  
######################################
######################################



AB_QDA_thai_braz = AB_QDA_allant( AB_thai, bin_cat_thai, AB_braz, bin_cat_braz, top_8 )




######################################
######################################
##                                  ##
##  PART 3                          ##
##  Thailand predicting Solomons    ## 
##                                  ##                                  
######################################
######################################


AB_QDA_thai_sol = AB_QDA_allant( AB_thai, bin_cat_thai, AB_sol, bin_cat_sol, top_8 )





######################################
######################################
##                                  ##
##  PART 4                          ##
##  Thailand predicting All         ## 
##                                  ##                                  
######################################
######################################


AB_QDA_thai_all = AB_QDA_Xval_train_subset_allant( AB_all, bin_cat_all, inf_cat_all, 
                                                     c("thai_current", "thai_recent", "thai_old", "thai_never", "TRC", "ARC", "VBDR"), 
                                                     top_8, 0.666, N_rep )




######################################
######################################
##                                  ##
##  PART 5                          ##
##  Brazil predicting Thailand      ## 
##                                  ##                                  
######################################
######################################


AB_QDA_braz_thai = AB_QDA_allant( AB_braz, bin_cat_braz, AB_thai, bin_cat_thai, top_8 )



######################################
######################################
##                                  ##
##  PART 6                          ##
##  Brazil predicting Brazil        ## 
##                                  ##                                  
######################################
######################################

AB_QDA_braz_braz = AB_QDA_Xval_allant( AB_braz, bin_cat_braz, top_8, 0.666, N_rep )




######################################
######################################
##                                  ##
##  PART 7                          ##
##  Brazil predicting Solomons      ## 
##                                  ##                                  
######################################
######################################



AB_QDA_braz_sol = AB_QDA_allant( AB_braz, bin_cat_braz, AB_sol, bin_cat_sol, top_8 )




######################################
######################################
##                                  ##
##  PART 8                          ##
##  Brazil predicting All           ## 
##                                  ##                                  
######################################
######################################


AB_QDA_braz_all = AB_QDA_Xval_train_subset_allant( AB_all, bin_cat_all, inf_cat_all, 
                                                     c("braz_current", "braz_recent", "braz_old", "braz_never", "TRC", "ARC", "VBDR"), 
                                                     top_8, 0.666, N_rep )




######################################
######################################
##                                  ##
##  PART 9                          ##
##  Solomons predicting Thailand    ## 
##                                  ##                                  
######################################
######################################


AB_QDA_sol_thai = AB_QDA_allant( AB_sol, bin_cat_sol, AB_thai, bin_cat_thai, top_8 )




######################################
######################################
##                                  ##
##  PART 10                         ##
##  Solomons predicting Brazil      ## 
##                                  ##                                  
######################################
######################################


AB_QDA_sol_braz = AB_QDA_allant( AB_sol, bin_cat_sol, AB_braz, bin_cat_braz, top_8 )





######################################
######################################
##                                  ##
##  PART 11                         ##
##  Solomons predicting Solomons    ## 
##                                  ##                                  
######################################
######################################

AB_QDA_sol_sol = AB_QDA_Xval_allant( AB_sol, bin_cat_sol, top_8, 0.666, N_rep )






######################################
######################################
##                                  ##
##  PART 12                         ##
##  Solomons predicting All         ## 
##                                  ##                                  
######################################
######################################


AB_QDA_sol_all = AB_QDA_Xval_train_subset_allant( AB_all, bin_cat_all, inf_cat_all, 
                                                    c("sol_current", "sol_recent", "sol_old", "sol_never", "TRC", "ARC", "VBDR"), 
                                                    top_8, 0.666, N_rep )





######################################
######################################
##                                  ##
##  PART 13                         ##
##  All predicting Thailand         ## 
##                                  ##                                  
######################################
######################################


AB_QDA_all_thai = AB_QDA_Xval_test_subset_allant( AB_all, bin_cat_all, inf_cat_all, 
                                                    c("thai_current", "thai_recent", "thai_old", "thai_never", "TRC", "ARC", "VBDR"), 
                                                    top_8, 0.666, N_rep )



######################################
######################################
##                                  ##
##  PART 14                         ##
##  All predicting Brazil           ## 
##                                  ##                                  
######################################
######################################


AB_QDA_all_braz = AB_QDA_Xval_test_subset_allant( AB_all, bin_cat_all, inf_cat_all, 
                                                    c("braz_current", "braz_recent", "braz_old", "braz_never", "TRC", "ARC", "VBDR"), 
                                                    top_8, 0.666, N_rep )



######################################
######################################
##                                  ##
##  PART 15                         ##
##  All predicting Solomons         ## 
##                                  ##                                  
######################################
######################################


AB_QDA_all_sol = AB_QDA_Xval_test_subset_allant( AB_all, bin_cat_all, inf_cat_all, 
                                                    c("sol_current", "sol_recent", "sol_old", "sol_never", "TRC", "ARC", "VBDR"), 
                                                    top_8, 0.666, N_rep )



######################################
######################################
##                                  ##
##  PART 16                         ##
##  All predicting All              ## 
##                                  ##                                  
######################################
######################################


AB_QDA_all_all = AB_QDA_Xval_allant( AB_all, bin_cat_all, top_8, 0.666, N_rep )




######################################
######################################
##        ##                        ##
##   ##   ##   ####  #   ## ######  ##
##  ###   ##  ##  ## ##  ##   ##    ##
##   ##   ##  ###### ### ##   ##    ##
##   ##   ##  ##  ## ## ###   ##    ##
##  ####  ##  ##  ## ##  ##   ##    ##
##        ##                        ## 
######################################
######################################

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
##  Thailand                        ## 
##                                  ##                                  
######################################
######################################

if( ncol(AB_thai)==60 )
{
	AB_thai = AB_thai[,top_8]
}

N_SS <- 1000

SS_cut <- seq(from=-11, to=-3.5, length=N_SS)

sens_mat <- matrix(NA, nrow=8, ncol=N_SS)
spec_mat <- matrix(NA, nrow=8, ncol=N_SS)


for(i in 1:8)
{
	AB_i  = AB_thai[which(is.na(AB_thai[,i])==FALSE),i]

	bin_cat_i = bin_cat_thai[which(is.na(AB_thai[,i])==FALSE)]


	for(j in 1:N_SS)
	{
		bin_pred <- rep( "old", length(bin_cat_i) )

		bin_pred[ which(AB_i > SS_cut[j]) ] <- "new"

		Cmat <- table( bin_cat_i, bin_pred  )

		if( ncol(Cmat)==1 )
		{
			if( colnames(Cmat)=="new" )
			{
				Cmat <- cbind( Cmat, c(0,0) )
			}else{
				if( colnames(Cmat)=="old" )
				{
					Cmat = cbind( c(0,0), Cmat )
				}
			}
		}

		sens_mat[i,j] <- Cmat[1,1]/(Cmat[1,1]+Cmat[1,2])
	
		spec_mat[i,j] <- Cmat[2,2]/(Cmat[2,1]+Cmat[2,2])	
	}
}




XX_VEC = c(0, 1)
YY_VEC = c(0, 1)


for(i in 1:8)
{
	XX_VEC = c( XX_VEC, 1 - spec_mat[i,] )
	YY_VEC = c( YY_VEC, sens_mat[i,] )
}

AB_thai_1 = STEP_CHULL( XX_VEC, YY_VEC )


######################################
######################################
##                                  ##
##  Brazil                          ## 
##                                  ##                                  
######################################
######################################

if( ncol(AB_braz)==60 )
{
	AB_braz = AB_braz[,top_8]
}

N_SS <- 1000

SS_cut <- seq(from=-11, to=-3.5, length=N_SS)

sens_mat <- matrix(NA, nrow=8, ncol=N_SS)
spec_mat <- matrix(NA, nrow=8, ncol=N_SS)


for(i in 1:8)
{
	AB_i  = AB_braz[which(is.na(AB_braz[,i])==FALSE),i]

	bin_cat_i = bin_cat_braz[which(is.na(AB_braz[,i])==FALSE)]


	for(j in 1:N_SS)
	{
		bin_pred <- rep( "old", length(bin_cat_i) )

		bin_pred[ which(AB_i > SS_cut[j]) ] <- "new"

		Cmat <- table( bin_cat_i, bin_pred  )

		if( ncol(Cmat)==1 )
		{
			if( colnames(Cmat)=="new" )
			{
				Cmat <- cbind( Cmat, c(0,0) )
			}else{
				if( colnames(Cmat)=="old" )
				{
					Cmat = cbind( c(0,0), Cmat )
				}
			}
		}

		sens_mat[i,j] <- Cmat[1,1]/(Cmat[1,1]+Cmat[1,2])
	
		spec_mat[i,j] <- Cmat[2,2]/(Cmat[2,1]+Cmat[2,2])	
	}
}



XX_VEC = c(0, 1)
YY_VEC = c(0, 1)


for(i in 1:8)
{
	XX_VEC = c( XX_VEC, 1 - spec_mat[i,] )
	YY_VEC = c( YY_VEC, sens_mat[i,] )
}

AB_braz_1 = STEP_CHULL( XX_VEC, YY_VEC )



######################################
######################################
##                                  ##
##  Solomons                        ## 
##                                  ##                                  
######################################
######################################

if( ncol(AB_sol)==60 )
{
	AB_sol = AB_sol[,top_8]
}

N_SS <- 1000

SS_cut <- seq(from=-11, to=-3.5, length=N_SS)

sens_mat <- matrix(NA, nrow=8, ncol=N_SS)
spec_mat <- matrix(NA, nrow=8, ncol=N_SS)


for(i in 1:8)
{
	AB_i  = AB_sol[which(is.na(AB_sol[,i])==FALSE),i]

	bin_cat_i = bin_cat_sol[which(is.na(AB_sol[,i])==FALSE)]


	for(j in 1:N_SS)
	{
		bin_pred <- rep( "old", length(bin_cat_i) )

		bin_pred[ which(AB_i > SS_cut[j]) ] <- "new"

		Cmat <- table( bin_cat_i, bin_pred  )

		if( ncol(Cmat)==1 )
		{
			if( colnames(Cmat)=="new" )
			{
				Cmat <- cbind( Cmat, c(0,0) )
			}else{
				if( colnames(Cmat)=="old" )
				{
					Cmat = cbind( c(0,0), Cmat )
				}
			}
		}

		sens_mat[i,j] <- Cmat[1,1]/(Cmat[1,1]+Cmat[1,2])
	
		spec_mat[i,j] <- Cmat[2,2]/(Cmat[2,1]+Cmat[2,2])	
	}
}




XX_VEC = c(0, 1)
YY_VEC = c(0, 1)


for(i in 1:8)
{
	XX_VEC = c( XX_VEC, 1 - spec_mat[i,] )
	YY_VEC = c( YY_VEC, sens_mat[i,] )
}

AB_sol_1 = STEP_CHULL( XX_VEC, YY_VEC )




######################################
######################################
##                                  ##
##  All regions                     ## 
##                                  ##                                  
######################################
######################################

if( ncol(AB_all)==60 )
{
	AB_all = AB_all[,top_8]
}

N_SS <- 1000

SS_cut <- seq(from=-11, to=-3.5, length=N_SS)

sens_mat <- matrix(NA, nrow=8, ncol=N_SS)
spec_mat <- matrix(NA, nrow=8, ncol=N_SS)


for(i in 1:8)
{
	AB_i  = AB_all[which(is.na(AB_all[,i])==FALSE),i]

	bin_cat_i = bin_cat_all[which(is.na(AB_all[,i])==FALSE)]


	for(j in 1:N_SS)
	{
		bin_pred <- rep( "old", length(bin_cat_i) )

		bin_pred[ which(AB_i > SS_cut[j]) ] <- "new"

		Cmat <- table( bin_cat_i, bin_pred  )

		if( ncol(Cmat)==1 )
		{
			if( colnames(Cmat)=="new" )
			{
				Cmat <- cbind( Cmat, c(0,0) )
			}else{
				if( colnames(Cmat)=="old" )
				{
					Cmat = cbind( c(0,0), Cmat )
				}
			}
		}

		sens_mat[i,j] <- Cmat[1,1]/(Cmat[1,1]+Cmat[1,2])
	
		spec_mat[i,j] <- Cmat[2,2]/(Cmat[2,1]+Cmat[2,2])	
	}
}




XX_VEC = c(0, 1)
YY_VEC = c(0, 1)


for(i in 1:8)
{
	XX_VEC = c( XX_VEC, 1 - spec_mat[i,] )
	YY_VEC = c( YY_VEC, sens_mat[i,] )
}

AB_all_1 = STEP_CHULL( XX_VEC, YY_VEC )



#######################################################
#######################################################
##                                                   ## 
##  #####  #####   ####   ####  #####  ####   ####   ## 
##  ##  ## ##  ## ##  ## ##  ## ##    ##     ##      ##
##  #####  #####  ##  ## ##     ####   ####   ####   ## 
##  ##     ## ##  ##  ## ##  ## ##        ##     ##  ##
##  ##     ##  ##  ####   ####  #####  ####   ####   ## 
##                                                   ##
#######################################################
#######################################################


######################################
######################################
##                                  ##
##  PART 1                          ##
##  Thailand predicting Thailand    ## 
##                                  ##                                  
######################################
######################################


AB_QDA_thai_thai_chull = list()


AB_QDA_thai_thai_chull[[1]] = AB_thai_1

XX_VEC = AB_thai_1[,1]
YY_VEC = AB_thai_1[,2]

for(i in 1:length(AB_QDA_thai_thai))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_thai_thai[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_thai_thai[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_thai_thai[[i]][[j]][,2] )
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


	AB_QDA_thai_thai_chull[[1+i]] = XY_hull_plot 
}




######################################
######################################
##                                  ##
##  PART 2                          ##
##  Thailand predicting Brazil      ## 
##                                  ##                                  
######################################
######################################


AB_QDA_thai_braz_chull = list()


AB_QDA_thai_braz_chull[[1]] = AB_braz_1

XX_VEC = AB_braz_1[,1]
YY_VEC = AB_braz_1[,2]

for(i in 1:length(AB_QDA_thai_braz))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_thai_braz[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_thai_braz[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_thai_braz[[i]][[j]][,2] )
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


	AB_QDA_thai_braz_chull[[1+i]] = XY_hull_plot 
}





######################################
######################################
##                                  ##
##  PART 3                          ##
##  Thailand predicting Solomons    ## 
##                                  ##                                  
######################################
######################################


AB_QDA_thai_sol_chull = list()

AB_QDA_thai_sol_chull[[1]] = AB_sol_1

XX_VEC = AB_sol_1[,1]
YY_VEC = AB_sol_1[,2]

for(i in 1:length(AB_QDA_thai_sol))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_thai_sol[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_thai_sol[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_thai_sol[[i]][[j]][,2] )
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


	AB_QDA_thai_sol_chull[[1+i]] = XY_hull_plot 
}






######################################
######################################
##                                  ##
##  PART 4                          ##
##  Thailand predicting All         ## 
##                                  ##                                  
######################################
######################################


AB_QDA_thai_all_chull = list()

AB_QDA_thai_all_chull[[1]] = AB_all_1

XX_VEC = AB_all_1[,1]
YY_VEC = AB_all_1[,2]


for(i in 1:length(AB_QDA_thai_all))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_thai_all[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_thai_all[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_thai_all[[i]][[j]][,2] )
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


	AB_QDA_thai_all_chull[[1+i]] = XY_hull_plot 
}




######################################
######################################
##                                  ##
##  PART 5                          ##
##  Brazil predicting Thailand      ## 
##                                  ##                                  
######################################
######################################


AB_QDA_braz_thai_chull = list()

AB_QDA_braz_thai_chull[[1]] = AB_thai_1

XX_VEC = AB_thai_1[,1]
YY_VEC = AB_thai_1[,2]

for(i in 1:length(AB_QDA_braz_thai))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_braz_thai[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_braz_thai[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_braz_thai[[i]][[j]][,2] )
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


	AB_QDA_braz_thai_chull[[1+i]] = XY_hull_plot 
}




######################################
######################################
##                                  ##
##  PART 6                          ##
##  Brazil predicting Brazil        ## 
##                                  ##                                  
######################################
######################################


AB_QDA_braz_braz_chull = list()


AB_QDA_braz_braz_chull[[1]] = AB_braz_1

XX_VEC = AB_braz_1[,1]
YY_VEC = AB_braz_1[,2]

for(i in 1:length(AB_QDA_braz_braz))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_braz_braz[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_braz_braz[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_braz_braz[[i]][[j]][,2] )
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


	AB_QDA_braz_braz_chull[[1+i]] = XY_hull_plot 
}





######################################
######################################
##                                  ##
##  PART 7                          ##
##  Brazil predicting Solomnos      ## 
##                                  ##                                  
######################################
######################################


AB_QDA_braz_sol_chull = list()

AB_QDA_braz_sol_chull[[1]] = AB_sol_1

XX_VEC = AB_sol_1[,1]
YY_VEC = AB_sol_1[,2]

for(i in 1:length(AB_QDA_braz_sol))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_braz_sol[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_braz_sol[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_braz_sol[[i]][[j]][,2] )
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


	AB_QDA_braz_sol_chull[[1+i]] = XY_hull_plot 
}





######################################
######################################
##                                  ##
##  PART 8                          ##
##  Brazil predicting All           ## 
##                                  ##                                  
######################################
######################################


AB_QDA_braz_all_chull = list()

AB_QDA_braz_all_chull[[1]] = AB_all_1

XX_VEC = AB_all_1[,1]
YY_VEC = AB_all_1[,2]

for(i in 1:length(AB_QDA_braz_all))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_braz_all[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_braz_all[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_braz_all[[i]][[j]][,2] )
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


	AB_QDA_braz_all_chull[[1+i]] = XY_hull_plot 
}







######################################
######################################
##                                  ##
##  PART 9                          ##
##  Solomons predicting Thailand    ## 
##                                  ##                                  
######################################
######################################


AB_QDA_sol_thai_chull = list()

AB_QDA_sol_thai_chull[[1]] = AB_thai_1

XX_VEC = AB_thai_1[,1]
YY_VEC = AB_thai_1[,2]

for(i in 1:length(AB_QDA_sol_thai))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_sol_thai[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_sol_thai[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_sol_thai[[i]][[j]][,2] )
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


	AB_QDA_sol_thai_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 10                         ##
##  Solomons predicting Brazil      ## 
##                                  ##                                  
######################################
######################################


AB_QDA_sol_braz_chull = list()

AB_QDA_sol_braz_chull[[1]] = AB_braz_1

XX_VEC = AB_braz_1[,1]
YY_VEC = AB_braz_1[,2]

for(i in 1:length(AB_QDA_sol_braz))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_sol_braz[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_sol_braz[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_sol_braz[[i]][[j]][,2] )
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


	AB_QDA_sol_braz_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 11                         ##
##  Solomons predicting Solomons    ## 
##                                  ##                                  
######################################
######################################


AB_QDA_sol_sol_chull = list()


AB_QDA_sol_sol_chull[[1]] = AB_sol_1

XX_VEC = AB_sol_1[,1]
YY_VEC = AB_sol_1[,2]

for(i in 1:length(AB_QDA_sol_sol))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_sol_sol[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_sol_sol[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_sol_sol[[i]][[j]][,2] )
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


	AB_QDA_sol_sol_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 12                         ##
##  Solomons predicting All         ## 
##                                  ##                                  
######################################
######################################


AB_QDA_sol_all_chull = list()

AB_QDA_sol_all_chull[[1]] = AB_all_1

XX_VEC = AB_all_1[,1]
YY_VEC = AB_all_1[,2]

for(i in 1:length(AB_QDA_sol_all))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_sol_all[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_sol_all[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_sol_all[[i]][[j]][,2] )
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


	AB_QDA_sol_all_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 13                         ##
##  All predicting Thai             ## 
##                                  ##                                  
######################################
######################################


AB_QDA_all_thai_chull = list()

AB_QDA_all_thai_chull[[1]] = AB_thai_1

XX_VEC = AB_thai_1[,1]
YY_VEC = AB_thai_1[,2]

for(i in 1:length(AB_QDA_all_thai))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_all_thai[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_all_thai[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_all_thai[[i]][[j]][,2] )
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


	AB_QDA_all_thai_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 14                         ##
##  All predicting Brazil           ## 
##                                  ##                                  
######################################
######################################


AB_QDA_all_braz_chull = list()

AB_QDA_all_braz_chull[[1]] = AB_braz_1

XX_VEC = AB_braz_1[,1]
YY_VEC = AB_braz_1[,2]

for(i in 1:length(AB_QDA_all_braz))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_all_braz[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_all_braz[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_all_braz[[i]][[j]][,2] )
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


	AB_QDA_all_braz_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 15                         ##
##  All predicting Solomons         ## 
##                                  ##                                  
######################################
######################################


AB_QDA_all_sol_chull = list()

AB_QDA_all_sol_chull[[1]] = AB_sol_1

XX_VEC = AB_sol_1[,1]
YY_VEC = AB_sol_1[,2]

for(i in 1:length(AB_QDA_all_sol))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_all_sol[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_all_sol[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_all_sol[[i]][[j]][,2] )
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


	AB_QDA_all_sol_chull[[1+i]] = XY_hull_plot 
}



######################################
######################################
##                                  ##
##  PART 16                         ##
##  All predicting All              ## 
##                                  ##                                  
######################################
######################################


AB_QDA_all_all_chull = list()

AB_QDA_all_all_chull[[1]] = AB_all_1

XX_VEC = AB_all_1[,1]
YY_VEC = AB_all_1[,2]

for(i in 1:length(AB_QDA_all_all))
{
	##################
	##  i antigens  ##
	##################

	for(j in 1:length(AB_QDA_all_all[[i]]))
	{
		XX_VEC = c( XX_VEC, AB_QDA_all_all[[i]][[j]][,1] )
		YY_VEC = c( YY_VEC, AB_QDA_all_all[[i]][[j]][,2] )
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


	AB_QDA_all_all_chull[[1+i]] = XY_hull_plot 
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




tiff(file="FigSX_AB_QDA_ROC_curves.tif", width=32, height=26, units="cm", res=500)


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
	points( x= AB_QDA_thai_thai_chull[[i]][,1], 
      	  y= AB_QDA_thai_thai_chull[[i]][,2], 
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
	points( x= AB_QDA_thai_braz_chull[[i]][,1], 
      	  y= AB_QDA_thai_braz_chull[[i]][,2], 
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
	points( x= AB_QDA_thai_sol_chull[[i]][,1], 
      	  y= AB_QDA_thai_sol_chull[[i]][,2], 
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
	points( x= AB_QDA_thai_all_chull[[i]][,1], 
      	  y= AB_QDA_thai_all_chull[[i]][,2], 
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
	points( x= AB_QDA_braz_thai_chull[[i]][,1], 
      	  y= AB_QDA_braz_thai_chull[[i]][,2], 
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
	points( x= AB_QDA_braz_braz_chull[[i]][,1], 
      	  y= AB_QDA_braz_braz_chull[[i]][,2], 
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
	points( x= AB_QDA_braz_sol_chull[[i]][,1], 
      	  y= AB_QDA_braz_sol_chull[[i]][,2], 
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
	points( x= AB_QDA_braz_all_chull[[i]][,1], 
      	  y= AB_QDA_braz_all_chull[[i]][,2], 
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
	points( x= AB_QDA_sol_thai_chull[[i]][,1], 
      	  y= AB_QDA_sol_thai_chull[[i]][,2], 
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
	points( x= AB_QDA_sol_braz_chull[[i]][,1], 
      	  y= AB_QDA_sol_braz_chull[[i]][,2], 
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
	points( x= AB_QDA_sol_sol_chull[[i]][,1], 
      	  y= AB_QDA_sol_sol_chull[[i]][,2], 
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
	points( x= AB_QDA_sol_all_chull[[i]][,1], 
      	  y= AB_QDA_sol_all_chull[[i]][,2], 
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
	points( x= AB_QDA_all_thai_chull[[i]][,1], 
      	  y= AB_QDA_all_thai_chull[[i]][,2], 
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
	points( x= AB_QDA_all_braz_chull[[i]][,1], 
      	  y= AB_QDA_all_braz_chull[[i]][,2], 
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
	points( x= AB_QDA_all_sol_chull[[i]][,1], 
      	  y= AB_QDA_all_sol_chull[[i]][,2], 
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
	points( x= AB_QDA_all_all_chull[[i]][,1], 
      	  y= AB_QDA_all_all_chull[[i]][,2], 
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





save.image("AB_QDA.RData")



