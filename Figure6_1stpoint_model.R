library(MASS)
library(ROCR)
library(randomForest)
library(binom)
library(survival)



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

brazil_data = read.csv("C:\\U\\GHIT\\NatMed_Paper\\Data\\proc\\brazil_ab_epi_data.csv")

solomon_data = read.csv("C:\\U\\GHIT\\NatMed_Paper\\Data\\proc\\solomon_ab_epi_data.csv")

control_data = read.csv("C:\\U\\GHIT\\NatMed_Paper\\Data\\proc\\control_ab_epi_data.csv")



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
		control_cat[i] <- "BRC"
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

par(mfrow=c(1,1))

N_tree = 2500


top_8 = c(50, 1, 47, 2, 30, 8, 58, 39)


N_rep = 10

train_prop = 0.666

TRC         = rep(NA, N_rep)
TRC_correct = rep(NA, N_rep)

BRC         = rep(NA, N_rep)
BRC_correct = rep(NA, N_rep)

ARC         = rep(NA, N_rep)
ARC_correct = rep(NA, N_rep)

VBDR         = rep(NA, N_rep)
VBDR_correct = rep(NA, N_rep)


Thai_current         = rep(NA, N_rep)
Thai_current_correct = rep(NA, N_rep)

Thai_recent         = rep(NA, N_rep)
Thai_recent_correct = rep(NA, N_rep)

Thai_old         = rep(NA, N_rep)
Thai_old_correct = rep(NA, N_rep)

Thai_never         = rep(NA, N_rep)
Thai_never_correct = rep(NA, N_rep)


Braz_current         = rep(NA, N_rep)
Braz_current_correct = rep(NA, N_rep)

Braz_recent         = rep(NA, N_rep)
Braz_recent_correct = rep(NA, N_rep)

Braz_old         = rep(NA, N_rep)
Braz_old_correct = rep(NA, N_rep)

Braz_never         = rep(NA, N_rep)
Braz_never_correct = rep(NA, N_rep)


Sol_current         = rep(NA, N_rep)
Sol_current_correct = rep(NA, N_rep)

Sol_recent         = rep(NA, N_rep)
Sol_recent_correct = rep(NA, N_rep)

Sol_old         = rep(NA, N_rep)
Sol_old_correct = rep(NA, N_rep)

Sol_never         = rep(NA, N_rep)
Sol_never_correct = rep(NA, N_rep)

alpha_top_8 = rep(NA, N_rep)



par(mar=c(3,3,2,1))

for(n in 1:N_rep)
{
	############################################
	## Prepare testing and training data

	index_train = sample( nrow(AB_all), train_prop*nrow(AB_all) )
	index_test  = setdiff( 1:nrow(AB_all), index_train )

	bin_cat_train = bin_cat_all[index_train]
	bin_cat_test  = bin_cat_all[index_test]

	inf_cat_train = inf_cat_all[index_train]
	inf_cat_test  = inf_cat_all[index_test]

	AB_train = AB_all[index_train,top_8]
	AB_test  = AB_all[index_test,top_8]


	if( length(which(is.na(AB_train))) > 0 )
	{
		bin_cat_train = bin_cat_train[-which(is.na(AB_train), arr.ind=TRUE)[,1]]
		inf_cat_train = inf_cat_train[-which(is.na(AB_train), arr.ind=TRUE)[,1]]

		AB_train = AB_train[-which(is.na(AB_train), arr.ind=TRUE)[,1],]
	}



	if( length(which(is.na(AB_test))) > 0 )
	{
		bin_cat_test = bin_cat_test[-which(is.na(AB_test), arr.ind=TRUE)[,1]]
		inf_cat_test = inf_cat_test[-which(is.na(AB_test), arr.ind=TRUE)[,1]]
	
		AB_test = AB_test[-which(is.na(AB_test), arr.ind=TRUE)[,1],]
	}



	Rforest_top_8 = randomForest( as.factor(bin_cat_train) ~ ., data=as.data.frame(AB_train), 
                                      importance=TRUE, ntree=N_tree )

	Rforest_top_8_pred_obj = predict( Rforest_top_8, newdata=AB_test, predict.all=TRUE, type="prob")

	Rforest_top_8_votes = rowSums(Rforest_top_8_pred_obj$individual=="old")/N_tree

	Rforest_top_8_pred = prediction(Rforest_top_8_votes, bin_cat_test)

	Rforest_top_8_perf_xy = performance(Rforest_top_8_pred, "spec", "sens")


	plot( x = 1 - Rforest_top_8_perf_xy@x.values[[1]], y = Rforest_top_8_perf_xy@y.values[[1]], type='s')




	alpha_top_8[n] = Rforest_top_8_perf_xy@alpha.values[[1]][which.min(abs(Rforest_top_8_perf_xy@x.values[[1]] - Rforest_top_8_perf_xy@y.values[[1]]))]



	TRC[n]         = length(which(inf_cat_test=="TRC"))
	TRC_correct[n] = TRC[n] - length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="TRC"),2] < alpha_top_8[n]))

	BRC[n]         = length(which(inf_cat_test=="BRC"))
	BRC_correct[n] = BRC[n] - length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="BRC"),2] < alpha_top_8[n]))

	ARC[n]         = length(which(inf_cat_test=="ARC"))
	ARC_correct[n] = ARC[n] - length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="ARC"),2] < alpha_top_8[n]))

	VBDR[n]         = length(which(inf_cat_test=="VBDR"))
	VBDR_correct[n] = VBDR[n] - length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="VBDR"),2] < alpha_top_8[n]))



	Thai_current[n]         = length(which(inf_cat_test=="thai_current"))
	Thai_current_correct[n] = length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="thai_current"),2] < alpha_top_8[n]))

	Thai_recent[n]         = length(which(inf_cat_test=="thai_recent"))
	Thai_recent_correct[n] = length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="thai_recent"),2] < alpha_top_8[n]))

	Thai_old[n]         = length(which(inf_cat_test=="thai_old"))
	Thai_old_correct[n] = Thai_old[n] - length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="thai_old"),2] < alpha_top_8[n]))
	
	Thai_never[n]         = length(which(inf_cat_test=="thai_never"))
	Thai_never_correct[n] = Thai_never[n] - length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="thai_never"),2] < alpha_top_8[n]))



	Braz_current[n]         = length(which(inf_cat_test=="braz_current"))
	Braz_current_correct[n] = length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="braz_current"),2] < alpha_top_8[n]))

	Braz_recent[n]         = length(which(inf_cat_test=="braz_recent"))
	Braz_recent_correct[n] = length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="braz_recent"),2] < alpha_top_8[n]))

	Braz_old[n]         = length(which(inf_cat_test=="braz_old"))
	Braz_old_correct[n] = Braz_old[n] - length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="braz_old"),2] < alpha_top_8[n]))

	Braz_never[n]         = length(which(inf_cat_test=="braz_never"))
	Braz_never_correct[n] = Braz_never[n] - length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="braz_never"),2] < alpha_top_8[n]))


	Sol_current[n]         = length(which(inf_cat_test=="sol_current"))
	Sol_current_correct [n]= length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="sol_current"),2] < alpha_top_8[n]))

	Sol_recent[n]         = length(which(inf_cat_test=="sol_recent"))
	Sol_recent_correct[n] = length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="sol_recent"),2] < alpha_top_8[n]))

	Sol_old[n]         = length(which(inf_cat_test=="sol_old"))
	Sol_old_correct[n] = Sol_old[n] - length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="sol_old"),2] < alpha_top_8[n]))

	Sol_never[n]         = length(which(inf_cat_test=="sol_never"))
	Sol_never_correct[n] = Sol_never[n] - length(which(Rforest_top_8_pred_obj[[1]][which(inf_cat_test=="sol_never"),2] < alpha_top_8[n]))
}


##median(

(TRC_correct + BRC_correct + ARC_correct + VBDR_correct + 
Thai_current_correct + Thai_recent_correct + Thai_old_correct + Thai_never_correct +
Braz_current_correct + Braz_recent_correct + Braz_old_correct + Braz_never_correct +
Sol_current_correct + Sol_recent_correct + Sol_old_correct + Sol_never_correct) /
(TRC + BRC + ARC + VBDR + 
Thai_current + Thai_recent + Thai_old + Thai_never +
Braz_current + Braz_recent + Braz_old + Braz_never +
Sol_current + Sol_recent + Sol_old + Sol_never )

##)






alpha_top8_med <- median(alpha_top_8)




AB_top8 <- AB_all[,top_8]

bin_cat_all = bin_cat_all[-which(is.na(AB_top8), arr.ind=TRUE)[,1]]
inf_cat_all = inf_cat_all[-which(is.na(AB_top8), arr.ind=TRUE)[,1]]

AB_top8 = AB_top8[-which(is.na(AB_top8), arr.ind=TRUE)[,1],]

Rforest_top_8_all = randomForest( as.factor(bin_cat_all) ~ ., data=as.data.frame(AB_top8), 
                                      importance=TRUE, ntree=N_tree )


#######################################
#######################################
##                                   ##
##  ##### #### #####   ####  ######  ##
##  ##     ##  ##  ## ##       ##    ##
##  ####   ##  #####   ####    ##    ##
##  ##     ##  ## ##      ##   ##    ##
##  ##    #### ##  ##  ####    ##    ##
##                                   ##
#######################################
#######################################

#######################################
##                                   ##
##  Thailand                         ##
##                                   ##
#######################################


thai_firstdata = read.csv("C:\\U\\GHIT\\NatMed_Paper\\First_time_point\\Thailand\\thailand_firstab_epi_data.csv")

thai_FD <- thai_firstdata[-unique(which( is.na(thai_firstdata[,19:26]), arr.ind=TRUE )[,1]),]
colnames(thai_FD)[19:26] <- colnames(AB_top8)
thai_FD[19:26] <- log(thai_FD[19:26])

thai_FD <- cbind( thai_FD, matrix(NA, nrow=nrow(thai_FD), ncol=3) )
colnames(thai_FD)[27:29] <- c("prediction", "T_censor", "Pv_event")


thai_top8_pred_obj = predict( Rforest_top_8_all, newdata=thai_FD[,19:26], predict.all=TRUE, type="prob")

thai_top8_votes = rowSums(thai_top8_pred_obj$individual=="old")/N_tree


thai_FD$prediction <- as.vector("nohpz")
thai_FD$prediction[which(thai_top8_votes < alpha_top8_med)] <- as.vector("hpz")


thai_FD$T_censor <- thai_FD$T_follow
thai_FD$T_censor[which(is.na(thai_FD$T_first_Pv) == FALSE)] <- thai_FD$T_first_Pv[which(is.na(thai_FD$T_first_Pv) == FALSE)]

thai_FD$Pv_event <- 0
thai_FD$Pv_event[ which(thai_FD$N_Pv_pos > 0) ] <- 1



thai_FD_hpz_surv     = Surv( time=thai_FD$T_censor[which(thai_FD$prediction=="hpz")], event=thai_FD$Pv_event[which(thai_FD$prediction=="hpz")] )
thai_FD_hpz_surv_fit = survfit( thai_FD_hpz_surv ~ 1 )


thai_FD_nohpz_surv     = Surv( time=thai_FD$T_censor[which(thai_FD$prediction=="nohpz")], event=thai_FD$Pv_event[which(thai_FD$prediction=="nohpz")] )
thai_FD_nohpz_surv_fit = survfit( thai_FD_nohpz_surv ~ 1 )


par(mfrow=c(1,2))

plot( thai_FD_hpz_surv_fit )

plot( thai_FD_nohpz_surv_fit )



thai_FD_trim <- thai_FD[ which(thai_FD$T_censor > 0), ]


thai_FD_trim_hpz_surv     = Surv( time=thai_FD_trim$T_censor[which(thai_FD_trim$prediction=="hpz")], event=thai_FD_trim$Pv_event[which(thai_FD_trim$prediction=="hpz")] )
thai_FD_trim_hpz_surv_fit = survfit( thai_FD_trim_hpz_surv ~ 1 )


thai_FD_trim_nohpz_surv     = Surv( time=thai_FD_trim$T_censor[which(thai_FD_trim$prediction=="nohpz")], event=thai_FD_trim$Pv_event[which(thai_FD_trim$prediction=="nohpz")] )
thai_FD_trim_nohpz_surv_fit = survfit( thai_FD_trim_nohpz_surv ~ 1 )



par(mfrow=c(1,2))

plot( thai_FD_trim_hpz_surv_fit )

plot( thai_FD_trim_nohpz_surv_fit )




thai_FD$prediction <- relevel( x=as.factor(thai_FD$prediction), ref="nohpz" )

thai_FD_coxph = coxph( Surv( time=thai_FD$T_censor, event=thai_FD$Pv_event ) ~ thai_FD$prediction )

summary( thai_FD_coxph )$coef

summary( thai_FD_coxph )$conf.int




thai_FD_trim_coxph = coxph( Surv( time=thai_FD_trim$T_censor, event=thai_FD_trim$Pv_event ) ~ thai_FD_trim$prediction )




#######################################
##                                   ##
##  Brazil                           ##
##                                   ##
#######################################


braz_firstdata = read.csv("C:\\U\\GHIT\\NatMed_Paper\\First_time_point\\Brazil\\brazil_firstab_epi_data.csv")

braz_FD <- braz_firstdata[-unique(which( is.na(braz_firstdata[,19:26]), arr.ind=TRUE )[,1]),]
colnames(braz_FD)[19:26] <- colnames(AB_top8)
braz_FD[19:26] <- log(braz_FD[19:26])

braz_FD <- cbind( braz_FD, matrix(NA, nrow=nrow(braz_FD), ncol=3) )
colnames(braz_FD)[27:29] <- c("prediction", "T_censor", "Pv_event")


braz_top8_pred_obj = predict( Rforest_top_8_all, newdata=braz_FD[,19:26], predict.all=TRUE, type="prob")

braz_top8_votes = rowSums(braz_top8_pred_obj$individual=="old")/N_tree


braz_FD$prediction <- as.vector("nohpz")
braz_FD$prediction[which(braz_top8_votes < alpha_top8_med)] <- as.vector("hpz")


braz_FD$T_censor <- braz_FD$T_follow
braz_FD$T_censor[which(is.na(braz_FD$T_first_Pv) == FALSE)] <- braz_FD$T_first_Pv[which(is.na(braz_FD$T_first_Pv) == FALSE)]

braz_FD$Pv_event <- 0
braz_FD$Pv_event[ which(braz_FD$N_Pv_pos > 0) ] <- 1



braz_FD_hpz_surv     = Surv( time=braz_FD$T_censor[which(braz_FD$prediction=="hpz")], event=braz_FD$Pv_event[which(braz_FD$prediction=="hpz")] )
braz_FD_hpz_surv_fit = survfit( braz_FD_hpz_surv ~ 1 )


braz_FD_nohpz_surv     = Surv( time=braz_FD$T_censor[which(braz_FD$prediction=="nohpz")], event=braz_FD$Pv_event[which(braz_FD$prediction=="nohpz")] )
braz_FD_nohpz_surv_fit = survfit( braz_FD_nohpz_surv ~ 1 )


par(mfrow=c(1,2))

plot( braz_FD_hpz_surv_fit )

plot( braz_FD_nohpz_surv_fit )




braz_FD_trim <- braz_FD[ which(braz_FD$T_censor > 0), ]


braz_FD_trim_hpz_surv     = Surv( time=braz_FD_trim$T_censor[which(braz_FD_trim$prediction=="hpz")], event=braz_FD_trim$Pv_event[which(braz_FD_trim$prediction=="hpz")] )
braz_FD_trim_hpz_surv_fit = survfit( braz_FD_trim_hpz_surv ~ 1 )


braz_FD_trim_nohpz_surv     = Surv( time=braz_FD_trim$T_censor[which(braz_FD_trim$prediction=="nohpz")], event=braz_FD_trim$Pv_event[which(braz_FD_trim$prediction=="nohpz")] )
braz_FD_trim_nohpz_surv_fit = survfit( braz_FD_trim_nohpz_surv ~ 1 )



par(mfrow=c(1,2))

plot( braz_FD_trim_hpz_surv_fit )

plot( braz_FD_trim_nohpz_surv_fit )






braz_FD$prediction <- relevel( x=as.factor(braz_FD$prediction), ref="nohpz" )

braz_FD_coxph = coxph( Surv( time=braz_FD$T_censor, event=braz_FD$Pv_event ) ~ braz_FD$prediction )

summary( braz_FD_coxph )$coef

summary( braz_FD_coxph )$conf.int


#######################################
##                                   ##
##  Solomons                         ##
##                                   ##
#######################################


sol_firstdata = read.csv("C:\\U\\GHIT\\NatMed_Paper\\First_time_point\\Solomons\\solomon_firstab_epi_data.csv")

sol_FD <- sol_firstdata[-unique(which( is.na(sol_firstdata[,18:25]), arr.ind=TRUE )[,1]),]
colnames(sol_FD)[18:25] <- colnames(AB_top8)
sol_FD[18:25] <- log(sol_FD[18:25])

sol_FD <- cbind( sol_FD, matrix(NA, nrow=nrow(sol_FD), ncol=3) )
colnames(sol_FD)[26:28] <- c("prediction", "T_censor", "Pv_event")


sol_top8_pred_obj = predict( Rforest_top_8_all, newdata=sol_FD[,18:25], predict.all=TRUE, type="prob")

sol_top8_votes = rowSums(sol_top8_pred_obj$individual=="old")/N_tree


sol_FD$prediction <- as.vector("nohpz")
sol_FD$prediction[which(sol_top8_votes < alpha_top8_med)] <- as.vector("hpz")


sol_FD$T_censor <- sol_FD$T_follow
sol_FD$T_censor[which(is.na(sol_FD$T_first_Pv) == FALSE)] <- sol_FD$T_first_Pv[which(is.na(sol_FD$T_first_Pv) == FALSE)]

sol_FD$Pv_event <- 0
sol_FD$Pv_event[ which(sol_FD$N_Pv_pos > 0) ] <- 1



sol_FD_hpz_surv     = Surv( time=sol_FD$T_censor[which(sol_FD$prediction=="hpz")], event=sol_FD$Pv_event[which(sol_FD$prediction=="hpz")] )
sol_FD_hpz_surv_fit = survfit( sol_FD_hpz_surv ~ 1 )


sol_FD_nohpz_surv     = Surv( time=sol_FD$T_censor[which(sol_FD$prediction=="nohpz")], event=sol_FD$Pv_event[which(sol_FD$prediction=="nohpz")] )
sol_FD_nohpz_surv_fit = survfit( sol_FD_nohpz_surv ~ 1 )


par(mfrow=c(1,2))

plot( sol_FD_hpz_surv_fit )

plot( sol_FD_nohpz_surv_fit )







sol_FD_trim <- sol_FD[ which(sol_FD$T_censor > 0), ]


sol_FD_trim_hpz_surv     = Surv( time=sol_FD_trim$T_censor[which(sol_FD_trim$prediction=="hpz")], event=sol_FD_trim$Pv_event[which(sol_FD_trim$prediction=="hpz")] )
sol_FD_trim_hpz_surv_fit = survfit( sol_FD_trim_hpz_surv ~ 1 )


sol_FD_trim_nohpz_surv     = Surv( time=sol_FD_trim$T_censor[which(sol_FD_trim$prediction=="nohpz")], event=sol_FD_trim$Pv_event[which(sol_FD_trim$prediction=="nohpz")] )
sol_FD_trim_nohpz_surv_fit = survfit( sol_FD_trim_nohpz_surv ~ 1 )



par(mfrow=c(1,2))

plot( sol_FD_trim_hpz_surv_fit )

plot( sol_FD_trim_nohpz_surv_fit )






sol_FD$prediction <- relevel( x=as.factor(sol_FD$prediction), ref="nohpz" )

sol_FD_coxph = coxph( Surv( time=sol_FD$T_censor, event=sol_FD$Pv_event ) ~ sol_FD$prediction )

summary( sol_FD_coxph )$coef

summary( sol_FD_coxph )$conf.int




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



line_seq_x <- c(0, 50, 100, 150, 200, 250, 300, 350, 400)
line_seq_y <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)



tiff(file="First6_1stpoint_model.tif", width=30, height=25, units="cm", res=500)

lay.mat <- rbind( c( 1, 2, 3 ),
                  c( 4, 4, 4 ),
                  c( 5, 6, 7 ),
                  c( 8, 8, 8 ) )
layout(lay.mat, heights=c(10,1,10,1))
layout.show(8)

par(mar=c(4.5,5,2,1.5))
par(mgp=c(3.1, 0.6, 0))


line.size = 2
lab.size  = 2
axis.size = 1.35
main.size = 1.8	

#######################################
##                                   ##
##  PANEL 1: Thailand                ##
##                                   ##
#######################################


plot(x=1000, y=1000, 
xlim=c(0, 400), ylim=c(0.0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n', bty='n',
xlab="time (days)", ylab="Proportion uninfected by PCR", 
main="(A) Thailand: testing first time point",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(-10,10), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=c(0,3000), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}



polygon( x=c( summary(thai_FD_hpz_surv_fit)$time[-length(summary(thai_FD_hpz_surv_fit)$time)], rev(summary(thai_FD_hpz_surv_fit)$time[-length(summary(thai_FD_hpz_surv_fit)$time)]) ), 
	   y=c( summary(thai_FD_hpz_surv_fit)$lower[-length(summary(thai_FD_hpz_surv_fit)$time)], rev(summary(thai_FD_hpz_surv_fit)$upper[-length(summary(thai_FD_hpz_surv_fit)$time)]) ),
	   col=rgb(99/256,184/256,255/256,0.25), border=NA)

polygon( x=c( summary(thai_FD_nohpz_surv_fit)$time[-length(summary(thai_FD_nohpz_surv_fit)$time)], rev(summary(thai_FD_nohpz_surv_fit)$time[-length(summary(thai_FD_nohpz_surv_fit)$time)]) ), 
	   y=c( summary(thai_FD_nohpz_surv_fit)$lower[-length(summary(thai_FD_nohpz_surv_fit)$time)], rev(summary(thai_FD_nohpz_surv_fit)$upper[-length(summary(thai_FD_nohpz_surv_fit)$time)]) ),
	   col=rgb(205/256,133/256,63/256,0.25), border=NA)

points( x=summary(thai_FD_hpz_surv_fit)$time[-length(summary(thai_FD_hpz_surv_fit)$time)], y=summary(thai_FD_hpz_surv_fit)$surv[-length(summary(thai_FD_hpz_surv_fit)$time)],
        col="steelblue1", type='s', lwd=line.size)

points( x=summary(thai_FD_nohpz_surv_fit)$time[-length(summary(thai_FD_nohpz_surv_fit)$time)], y=summary(thai_FD_nohpz_surv_fit)$surv[-length(summary(thai_FD_nohpz_surv_fit)$time)],
        col="peru", type='s', lwd=line.size)


axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=axis.size ) 

axis(1, at = c(0, 50, 100, 150, 200, 250, 300, 350, 400),
        labels = c(0, 50, 100, 150, 200, 250, 300, 350, 400), 
        cex.axis=axis.size )


text( x = 130, y=0.05, 
      labels="HR = 5.88 (3.80, 9.10)",
      cex=1.75 )

#######################################
##                                   ##
##  PANEL 2: Brazil                  ##
##                                   ##
#######################################

plot(x=1000, y=1000, 
xlim=c(0, 400), ylim=c(0.0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n', bty='n',
xlab="time (days)", ylab="Proportion uninfected by PCR", 
main="(B) Brazil: testing first time point",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(-10,10), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=c(0,3000), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}


polygon( x=c( summary(braz_FD_hpz_surv_fit)$time[-length(summary(braz_FD_hpz_surv_fit)$time)], rev(summary(braz_FD_hpz_surv_fit)$time[-length(summary(braz_FD_hpz_surv_fit)$time)]) ), 
	   y=c( summary(braz_FD_hpz_surv_fit)$lower[-length(summary(braz_FD_hpz_surv_fit)$time)], rev(summary(braz_FD_hpz_surv_fit)$upper[-length(summary(braz_FD_hpz_surv_fit)$time)]) ),
	   col=rgb(99/256,184/256,255/256,0.25), border=NA)

polygon( x=c( summary(braz_FD_nohpz_surv_fit)$time[-length(summary(braz_FD_nohpz_surv_fit)$time)], rev(summary(braz_FD_nohpz_surv_fit)$time[-length(summary(braz_FD_nohpz_surv_fit)$time)]) ), 
	   y=c( summary(braz_FD_nohpz_surv_fit)$lower[-length(summary(braz_FD_nohpz_surv_fit)$time)], rev(summary(braz_FD_nohpz_surv_fit)$upper[-length(summary(braz_FD_nohpz_surv_fit)$time)]) ),
	   col=rgb(205/256,133/256,63/256,0.25), border=NA)

points( x=summary(braz_FD_hpz_surv_fit)$time[-length(summary(braz_FD_hpz_surv_fit)$time)], y=summary(braz_FD_hpz_surv_fit)$surv[-length(summary(braz_FD_hpz_surv_fit)$time)],
        col="steelblue1", type='s', lwd=line.size)

points( x=summary(braz_FD_nohpz_surv_fit)$time[-length(summary(braz_FD_nohpz_surv_fit)$time)], y=summary(braz_FD_nohpz_surv_fit)$surv[-length(summary(braz_FD_nohpz_surv_fit)$time)],
        col="peru", type='s', lwd=line.size)


axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=axis.size ) 

axis(1, at = c(0, 50, 100, 150, 200, 250, 300, 350, 400),
        labels = c(0, 50, 100, 150, 200, 250, 300, 350, 400), 
        cex.axis=axis.size )



text( x = 130, y=0.05, 
      labels="HR = 3.23 (2.52, 4.13)",
      cex=1.75 )




#######################################
##                                   ##
##  PANEL 3: Solomons                ##
##                                   ##
#######################################

plot(x=1000, y=1000, 
xlim=c(0, 400), ylim=c(0.0,1.01),
xaxs='i', yaxs='i', xaxt='n', yaxt='n', bty='n',
xlab="time (days)", ylab="Proportion uninfected by PCR", 
main="(C) Solomons: testing first time point",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)



for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(-10,10), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=c(0,3000), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}


polygon( x=c( summary(sol_FD_hpz_surv_fit)$time[-length(summary(sol_FD_hpz_surv_fit)$time)], rev(summary(sol_FD_hpz_surv_fit)$time[-length(summary(sol_FD_hpz_surv_fit)$time)]) ), 
	   y=c( summary(sol_FD_hpz_surv_fit)$lower[-length(summary(sol_FD_hpz_surv_fit)$time)], rev(summary(sol_FD_hpz_surv_fit)$upper[-length(summary(sol_FD_hpz_surv_fit)$time)]) ),
	   col=rgb(99/256,184/256,255/256,0.25), border=NA)

polygon( x=c( summary(sol_FD_nohpz_surv_fit)$time[-length(summary(sol_FD_nohpz_surv_fit)$time)], rev(summary(sol_FD_nohpz_surv_fit)$time[-length(summary(sol_FD_nohpz_surv_fit)$time)]) ), 
	   y=c( summary(sol_FD_nohpz_surv_fit)$lower[-length(summary(sol_FD_nohpz_surv_fit)$time)], rev(summary(sol_FD_nohpz_surv_fit)$upper[-length(summary(sol_FD_nohpz_surv_fit)$time)]) ),
	   col=rgb(205/256,133/256,63/256,0.25), border=NA)

points( x=summary(sol_FD_hpz_surv_fit)$time[-length(summary(sol_FD_hpz_surv_fit)$time)], y=summary(sol_FD_hpz_surv_fit)$surv[-length(summary(sol_FD_hpz_surv_fit)$time)],
        col="steelblue1", type='s', lwd=line.size)

points( x=summary(sol_FD_nohpz_surv_fit)$time[-length(summary(sol_FD_nohpz_surv_fit)$time)], y=summary(sol_FD_nohpz_surv_fit)$surv[-length(summary(sol_FD_nohpz_surv_fit)$time)],
        col="peru", type='s', lwd=line.size)


axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=axis.size ) 

axis(1, at = c(0, 50, 100, 150, 200, 250, 300, 350, 400),
        labels = c(0, 50, 100, 150, 200, 250, 300, 350, 400), 
        cex.axis=axis.size )


text( x = 135, y=0.05, 
      labels="HR = 8.55 (6.52, 11.22)",
      cex=1.75 )


#######################################
##                                   ##
##  LEGEND                           ##
##                                   ##
#######################################


par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("first time point classified positive", "first time point classified negative"), 
 	 col = c("steelblue1", "peru"), 
 	 lty = c("solid", "solid"),	
	 lwd = c(3,3),
       ncol=2, cex=2, bty="n" )


par(mar=c(4.5,5,2,1.5))

######################################
######################################
##                                  ## 
##  PANEL 4                         ##
##  Thailand model prediction       ##
##                                  ##
######################################   
######################################

load("C:\\U\\GHIT\\NatMed_Paper\\Pv_trans_model\\Thailand\\Processed_Thailand.RData")


time_plot <- THAILAND$time/365



line_seq_x <- c(2017, 2018, 2019, 2020, 2021, 2022)
line_seq_y <- c(0, 0.02, 0.04, 0.06, 0.08, 0.1)

plot(x=10000, y=10000, 
xlim=c(2018, 2022), ylim=c(0,0.1),
xaxs='i', yaxs='i', xaxt='n', yaxt='n', bty='n',
xlab="time (years)", ylab="PCR prevalence",
main="(D) Thailand: modelling treatment",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(-10,10), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=c(0,3000), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}

arrows(x0=2020, y0=0.1*0.95, x1=2020, y1=0.08, length=0.0, lwd=2)
arrows(x0=2021, y0=0.1*0.95, x1=2021, y1=0.08, length=0.0, lwd=2)

points(x=2020, y=0.5*(0.1*0.95 + 0.08), pch=25, cex=2)
points(x=2021, y=0.5*(0.1*0.95 + 0.08), pch=25, cex=2)


#################
## Baseline 
##
##points( x=time_plot, y=baseline_PCR_quant_smooth[1,], 
##col="black", type='l', lwd=line.size)
##
##polygon(x=c(time_plot, rev(time_plot)), 
##	  y=c( baseline_PCR_quant_smooth[2,], rev(baseline_PCR_quant_smooth[3,]) ),
##	  col=rgb(30/256,144/256,255/256,0.1), border=NA)


#################
## MDA 

points( x=time_plot, y=PQ_MDA_PCR_quant_smooth[1,], 
col="grey39", type='l', lwd=line.size)

polygon(x=c(time_plot, rev(time_plot)), 
	  y=c( PQ_MDA_PCR_quant_smooth[2,], rev(PQ_MDA_PCR_quant_smooth[3,]) ),
	  col=rgb(99/256,99/256,99/256,0.25), border=NA)


#################
## MSAT

points( x=time_plot, y=PQ_MSAT_PCR_quant_smooth[1,], 
col="sienna1", type='l', lwd=line.size)

polygon(x=c(time_plot, rev(time_plot)), 
	  y=c( PQ_MSAT_PCR_quant_smooth[2,], rev(PQ_MSAT_PCR_quant_smooth[3,]) ),
	  col=rgb(255/256,130/256,71/256,0.25), border=NA)


#################
## STAT_80_80

points( x=time_plot, y=PQ_STAT_80_80_PCR_quant_smooth[1,], 
col="forestgreen", type='l', lwd=line.size)

polygon(x=c(time_plot, rev(time_plot)), 
	  y=c( PQ_STAT_80_80_PCR_quant_smooth[2,], rev(PQ_STAT_80_80_PCR_quant_smooth[3,]) ),
	  col=rgb(34/256,139/256,34/256,0.25), border=NA)


#################
## Data

points( x=Thai_PvPCR_times, y=Thai_PvPCR_bins[1,], 
col="red", pch=19, cex=1.5)

for(j in 1:N_thai_bins)
{
	arrows(x0=Thai_PvPCR_times[j], y0=Thai_PvPCR_bins[2,j], 
             x1=Thai_PvPCR_times[j], y1=Thai_PvPCR_bins[3,j], 
             length=0.03, angle=90, code=3, col="red", lwd=1)	
}



axis(2, at=c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1), 
        labels=c("0%", "2%", "4%", "6%", "8%", "10%"), las=2, cex.axis=axis.size ) 

axis(1, at = c(2018, 2019, 2020, 2021, 2022),
labels=c("2018", "2019", "2020", "2021", "2022"), cex.axis=axis.size )




1 - PQ_MDA_PCR_quant_smooth[1,which.min( abs(time_plot - 2021.5) )]/baseline_PCR_quant_smooth[1,which.min( abs(time_plot - 2019.5) )]

1 - PQ_MSAT_PCR_quant_smooth[1,which.min( abs(time_plot - 2021.5) )]/baseline_PCR_quant_smooth[1,which.min( abs(time_plot - 2019.5) )]

1 - PQ_STAT_80_80_PCR_quant_smooth[1,which.min( abs(time_plot - 2021.5) )]/baseline_PCR_quant_smooth[1,which.min( abs(time_plot - 2019.5) )]




######################################
######################################
##                                  ## 
##  PANEL 5                         ##
##  Brazil model prediction         ##
##                                  ##
######################################   
######################################

load("C:\\U\\GHIT\\NatMed_Paper\\Pv_trans_model\\Brazil\\Processed_Brazil.RData")



time_plot <- BRAZIL$time/365


line_seq_x <- c(2017, 2018, 2019, 2020, 2021, 2022)
line_seq_y <- c(0, 0.02, 0.04, 0.06, 0.08, 0.1)

plot(x=10000, y=10000, 
xlim=c(2018, 2022), ylim=c(0,0.1),
xaxs='i', yaxs='i', xaxt='n', yaxt='n', bty='n',
xlab="time (years)", ylab="PCR prevalence",
main="(E) Brazil: modelling treatment",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(-10,10), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=c(0,3000), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}


arrows(x0=2020, y0=0.1*0.95, x1=2020, y1=0.08, length=0.0, lwd=2)
arrows(x0=2021, y0=0.1*0.95, x1=2021, y1=0.08, length=0.0, lwd=2)

points(x=2020, y=0.5*(0.1*0.95 + 0.08), pch=25, cex=2)
points(x=2021, y=0.5*(0.1*0.95 + 0.08), pch=25, cex=2)


#################
## Baseline 
##
##points( x=time_plot, y=baseline_PCR_quant_smooth[1,], 
##col="black", type='l', lwd=line.size)
##
##polygon(x=c(time_plot, rev(time_plot)), 
##	  y=c( baseline_PCR_quant_smooth[2,], rev(baseline_PCR_quant_smooth[3,]) ),
##	  col=rgb(30/256,144/256,255/256,0.1), border=NA)


#################
## MDA 

points( x=time_plot, y=PQ_MDA_PCR_quant_smooth[1,], 
col="grey39", type='l', lwd=line.size)

polygon(x=c(time_plot, rev(time_plot)), 
	  y=c( PQ_MDA_PCR_quant_smooth[2,], rev(PQ_MDA_PCR_quant_smooth[3,]) ),
	  col=rgb(99/256,99/256,99/256,0.25), border=NA)


#################
## MSAT

points( x=time_plot, y=PQ_MSAT_PCR_quant_smooth[1,], 
col="sienna1", type='l', lwd=line.size)

polygon(x=c(time_plot, rev(time_plot)), 
	  y=c( PQ_MSAT_PCR_quant_smooth[2,], rev(PQ_MSAT_PCR_quant_smooth[3,]) ),
	  col=rgb(255/256,130/256,71/256,0.25), border=NA)


#################
## STAT_80_80

points( x=time_plot, y=PQ_STAT_80_80_PCR_quant_smooth[1,], 
col="forestgreen", type='l', lwd=line.size)

polygon(x=c(time_plot, rev(time_plot)), 
	  y=c( PQ_STAT_80_80_PCR_quant_smooth[2,], rev(PQ_STAT_80_80_PCR_quant_smooth[3,]) ),
	  col=rgb(34/256,139/256,34/256,0.25), border=NA)


#################
## Data

points( x=Braz_PvPCR_times, y=Braz_PvPCR_bins[1,], 
col="red", pch=19, cex=1.5)

for(j in 1:N_braz_bins)
{
	arrows(x0=Braz_PvPCR_times[j], y0=Braz_PvPCR_bins[2,j], 
             x1=Braz_PvPCR_times[j], y1=Braz_PvPCR_bins[3,j], 
             length=0.03, angle=90, code=3, col="red", lwd=1)	
}



axis(2, at=c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1), 
        labels=c("0%", "2%", "4%", "6%", "8%", "10%"), las=2, cex.axis=axis.size ) 

axis(1, at = c(2018, 2019, 2020, 2021, 2022),
labels=c("2018", "2019", "2020", "2021", "2022"), cex.axis=axis.size )



1 - PQ_MDA_PCR_quant_smooth[1,which.min( abs(time_plot - 2021.5) )]/baseline_PCR_quant_smooth[1,which.min( abs(time_plot - 2019.5) )]

1 - PQ_MSAT_PCR_quant_smooth[1,which.min( abs(time_plot - 2021.5) )]/baseline_PCR_quant_smooth[1,which.min( abs(time_plot - 2019.5) )]

1 - PQ_STAT_80_80_PCR_quant_smooth[1,which.min( abs(time_plot - 2021.5) )]/baseline_PCR_quant_smooth[1,which.min( abs(time_plot - 2019.5) )]





######################################
######################################
##                                  ## 
##  PANEL 6                         ##
##  Solomon model prediction        ##
##                                  ##
######################################   
######################################

load("C:\\U\\GHIT\\NatMed_Paper\\Pv_trans_model\\Solomons\\Processed_Solomons.RData")


time_plot <- SOLOMON$time/365



line_seq_x <- c(2017, 2018, 2019, 2020, 2021, 2022)
line_seq_y <- c(0, 0.05, 0.1, 0.15, 0.2)

plot(x=10000, y=10000, 
xlim=c(2018, 2022), ylim=c(0,0.2),
xaxs='i', yaxs='i', xaxt='n', yaxt='n', bty='n',
xlab="time (years)", ylab="PCR prevalence",
main="(F) Solomons: modelling treatment",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(-10,10), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=c(0,3000), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}


arrows(x0=2020, y0=0.2*0.95, x1=2020, y1=0.16, length=0.0, lwd=2)
arrows(x0=2021, y0=0.2*0.95, x1=2021, y1=0.16, length=0.0, lwd=2)

points(x=2020, y=0.5*(0.2*0.95 + 0.16), pch=25, cex=2)
points(x=2021, y=0.5*(0.2*0.95 + 0.16), pch=25, cex=2)



#################
## Baseline 
##
##points( x=time_plot, y=baseline_PCR_quant_smooth[1,], 
##col="black", type='l', lwd=line.size)
##
##polygon(x=c(time_plot, rev(time_plot)), 
##	  y=c( baseline_PCR_quant_smooth[2,], rev(baseline_PCR_quant_smooth[3,]) ),
##	  col=rgb(30/256,144/256,255/256,0.1), border=NA)


#################
## MDA 

points( x=time_plot, y=PQ_MDA_PCR_quant_smooth[1,], 
col="grey39", type='l', lwd=line.size)

polygon(x=c(time_plot, rev(time_plot)), 
	  y=c( PQ_MDA_PCR_quant_smooth[2,], rev(PQ_MDA_PCR_quant_smooth[3,]) ),
	  col=rgb(99/256,99/256,99/256,0.25), border=NA)


#################
## MSAT

points( x=time_plot, y=PQ_MSAT_PCR_quant_smooth[1,], 
col="sienna1", type='l', lwd=line.size)

polygon(x=c(time_plot, rev(time_plot)), 
	  y=c( PQ_MSAT_PCR_quant_smooth[2,], rev(PQ_MSAT_PCR_quant_smooth[3,]) ),
	  col=rgb(255/256,130/256,71/256,0.25), border=NA)


#################
## STAT_80_80

points( x=time_plot, y=PQ_STAT_80_80_PCR_quant_smooth[1,], 
col="forestgreen", type='l', lwd=line.size)

polygon(x=c(time_plot, rev(time_plot)), 
	  y=c( PQ_STAT_80_80_PCR_quant_smooth[2,], rev(PQ_STAT_80_80_PCR_quant_smooth[3,]) ),
	  col=rgb(34/256,139/256,34/256,0.25), border=NA)


#################
## Data

points( x=Sol_PvPCR_times, y=Sol_PvPCR_bins[1,], 
col="red", pch=19, cex=1.5)

for(j in 1:N_sol_bins)
{
	arrows(x0=Sol_PvPCR_times[j], y0=Sol_PvPCR_bins[2,j], 
             x1=Sol_PvPCR_times[j], y1=Sol_PvPCR_bins[3,j], 
             length=0.03, angle=90, code=3, col="red", lwd=1)	
}





axis(2, at=c(0.0, 0.05, 0.1, 0.15, 0.2), 
        labels=c("0%", "5%", "10%", "15%", "20%"), las=2, cex.axis=axis.size ) 

axis(1, at = c(2018, 2019, 2020, 2021, 2022),
labels=c("2018", "2019", "2020", "2021", "2022"), cex.axis=axis.size )


1 - PQ_MDA_PCR_quant_smooth[1,which.min( abs(time_plot - 2021.5) )]/baseline_PCR_quant_smooth[1,which.min( abs(time_plot - 2019.5) )]

1 - PQ_MSAT_PCR_quant_smooth[1,which.min( abs(time_plot - 2021.5) )]/baseline_PCR_quant_smooth[1,which.min( abs(time_plot - 2019.5) )]

1 - PQ_STAT_80_80_PCR_quant_smooth[1,which.min( abs(time_plot - 2021.5) )]/baseline_PCR_quant_smooth[1,which.min( abs(time_plot - 2019.5) )]



######################################
######################################
##                                  ## 
##  Legend                          ##
##                                  ##
######################################   
######################################


par(mar = c(0,0,0,0))
plot.new()


legend(x='center', 
       legend = c("data: PCR prevalence", "population treatment", "MDA", "MSAT: LM", "SeroTAT: 80%, 80%"),
       col = c( "red", "black", "grey39", "sienna1", "forestgreen"  ), 
	 lty = c("solid", "solid", "solid", "solid", "solid"),	
	 pch = c(19, 25, NA, NA, NA ),
	 pt.cex = c(2, 2, 2, 2, 2),
	 lwd = c(3, 3, 3, 3, 3),
       ncol=5, cex=1.5, bty="n" )




dev.off()










