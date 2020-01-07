N_rep_sann = 100

N_ant = 104

N_panel_max = 40
panel_seq <- 1:N_panel_max

LL_panel <- matrix(NA, nrow=N_panel_max, ncol=N_rep_sann)

N_ant_panel <- matrix(NA, nrow=N_panel_max, ncol=N_rep_sann)

ant_prob <- matrix(0, nrow=N_ant, ncol=N_panel_max)




############################################
############################################
##       ##                               ##
##   ##  ##   ####   ####  #   ## #   ##  ##
##  ###  ##  ##     ##  ## ##  ## ##  ##  ##
##   ##  ##   ####  ###### ### ## ### ##  ##
##   ##  ##      ## ##  ## ## ### ## ###  ##
##  #### ##   ####  ##  ## ##  ## ##  ##  ##
##       ##                               ##
############################################
############################################

for(n in 1:N_rep_sann)
{
	sann_name <- paste("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//SANN//Alg_A//Output//sann_fit_", n, ".txt", sep="")

	sann <- read.table(sann_name)

	LL_panel[1:nrow(sann),n] <- sann[,2] 

	N_ant_panel[1:nrow(sann),n] <- rowSums( sann[,3:(2+N_ant)] )

	ant_prob[,1:nrow(sann)] <- ant_prob[,1:nrow(sann)] + t(sann[,3:(2+N_ant)])

	if( nrow(sann)<N_panel_max )
	{ 
		print( n )
		print( nrow(sann) ) 
	}
}


#########################################
##                                     ##
##  PANEL 1 (OLD)                      ## 
##  Diminishing returns of likelihood  ##
##                                     ##
#########################################

LL_panel_quant <- matrix(NA, nrow=5, ncol=N_panel_max)

for(j in 1:N_panel_max)
{
	LL_panel_quant[,j] <- quantile( LL_panel[j,], prob=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE )
}


#########################################
##                                     ##
##  PANEL 2 (OLD)                      ## 
##  Number of antigens                 ##
##                                     ##
#########################################


Thai_correct <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round2//Figure_gen//Thai_correct.csv")

Braz_correct <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round2//Figure_gen//Braz_correct.csv")




###################################
## Top 104 antigens

top_104_read <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//ANTIGENS_to_include.csv")

antig_ID <- as.vector(top_104_read[,1])


###################################
## SIGMA_A

thai_SIGMA_A <- read.table("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//SANN//Input_Data//thai_SIGMA_A_DS.txt")

thai_SIGMA_A <- as.matrix(thai_SIGMA_A)
rownames(thai_SIGMA_A) <- antig_ID
colnames(thai_SIGMA_A) <- antig_ID



braz_SIGMA_A <- read.table("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//SANN//Input_Data//braz_SIGMA_A_DS.txt")

braz_SIGMA_A <- as.matrix(braz_SIGMA_A)
rownames(braz_SIGMA_A) <- antig_ID
colnames(braz_SIGMA_A) <- antig_ID




###################################
## SIGMA_r


thai_SIGMA_r <- read.table("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//SANN//Input_Data//thai_SIGMA_r_DS.txt")

thai_SIGMA_r <- as.matrix(thai_SIGMA_r)
rownames(thai_SIGMA_r) <- antig_ID
colnames(thai_SIGMA_r) <- antig_ID


braz_SIGMA_r <- read.table("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//SANN//Input_Data//braz_SIGMA_r_DS.txt")

braz_SIGMA_r <- as.matrix(braz_SIGMA_r)
rownames(braz_SIGMA_r) <- antig_ID
colnames(braz_SIGMA_r) <- antig_ID



###################################
## SIGMA_m

thai_SIGMA_m <- read.table("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//SANN//Input_Data//thai_SIGMA_m_DS.txt")

thai_SIGMA_m <- as.matrix(thai_SIGMA_m)
rownames(thai_SIGMA_m) <- antig_ID
colnames(thai_SIGMA_m) <- antig_ID


braz_SIGMA_m <- read.table("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//SANN//Input_Data//braz_SIGMA_m_DS.txt")

braz_SIGMA_m <- as.matrix(braz_SIGMA_m)
rownames(braz_SIGMA_m) <- antig_ID
colnames(braz_SIGMA_m) <- antig_ID



###################################
## A_GMT_r_mean

thai_A_GMT_r_mean_read <- read.table("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//SANN//Input_Data//thai_A_GMT_r_mean_DS.txt")
rownames(thai_A_GMT_r_mean_read) <- antig_ID

thai_r_mean <- thai_A_GMT_r_mean_read[,2]

thai_d_half <- -log(2)/thai_r_mean


thai_noise <- diag( thai_SIGMA_A ) + 12*12*diag( thai_SIGMA_r ) + diag( thai_SIGMA_m )




braz_A_GMT_r_mean_read <- read.table("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//SANN//Input_Data//braz_A_GMT_r_mean_DS.txt")
rownames(braz_A_GMT_r_mean_read) <- antig_ID

braz_r_mean <- braz_A_GMT_r_mean_read[,2]

braz_d_half <- -log(2)/braz_r_mean


braz_noise <- diag( braz_SIGMA_A ) + 12*12*diag( braz_SIGMA_r ) + diag( braz_SIGMA_m )





d_half_order <- order(thai_d_half + braz_d_half)

ant_prob <- ant_prob[d_half_order,]

thai_d_half <- thai_d_half[d_half_order]
braz_d_half <- braz_d_half[d_half_order]

thai_noise <- thai_noise[d_half_order]
braz_noise <- braz_noise[d_half_order]



d_half_plot <- rep(NA, 2*N_ant)

d_half_plot[seq(from=1,by=2,length=N_ant)] <- thai_d_half
d_half_plot[seq(from=2,by=2,length=N_ant)] <- braz_d_half



noise_plot <- rep(NA, 2*N_ant)

noise_plot[seq(from=1,by=2,length=N_ant)] <- thai_noise
noise_plot[seq(from=2,by=2,length=N_ant)] <- braz_noise



#############################################################
#############################################################
##          ##                                             ##
##   ####   ##  ##    #### #   ##   #     #  ####  ####    ##
##  ##  ##  ##  ##     ##  ##  ##   ##   ## ##  ## ## ##   ##
##     ##   ##  ##     ##  ### ##   ####### ##  ## ##  ##  ##
##    ##    ##  ##     ##  ## ###   ## # ## ##  ## ## ##   ##
##   #####  ##  ##### #### ##  ##   ##   ##  ####  ####    ##
##          ##                                             ##
#############################################################
#############################################################

library(lme4)


####################
####################
####################
###              ### 
###   THAILAND   ###
###              ###
####################
####################
####################

N_thai_part = 32
N_sam  = 4
N_rep  = 3


############################################################
## Data - visit week 0

Thai_w0_read  <- read.csv("C://U//GHIT//1_AlphaScreen_Data//2_PLOS_NTDS_manuscript//Analysis//Thai_Alpha_w0_adj_reps.csv")

Thai_w0_read[384+(1:384),1]   <- Thai_w0_read[1:384,1]
Thai_w0_read[2*384+(1:384),1] <- Thai_w0_read[1:384,1]

Thai_w0 <- cbind( Thai_w0_read, rep("w0", nrow(Thai_w0_read)) )
colnames(Thai_w0)[ncol(Thai_w0)] <- "visit"


thai_part_ID <- colnames(Thai_w0)[2:(1+N_thai_part)]


############################################################
## Data - visit week 12

Thai_w12_read  <- read.csv("C://U//GHIT//1_AlphaScreen_Data//2_PLOS_NTDS_manuscript//Analysis//Thai_Alpha_w12_adj_reps.csv")

Thai_w12_read[384+(1:384),1]   <- Thai_w12_read[1:384,1]
Thai_w12_read[2*384+(1:384),1] <- Thai_w12_read[1:384,1]

Thai_w12 <- Thai_w12_read



############################################################
## Data - visit week 24

Thai_w24_read  <- read.csv("C://U//GHIT//1_AlphaScreen_Data//2_PLOS_NTDS_manuscript//Analysis//Thai_Alpha_w24_adj_reps.csv")

Thai_w24_read[384+(1:384),1]   <- Thai_w24_read[1:384,1]
Thai_w24_read[2*384+(1:384),1] <- Thai_w24_read[1:384,1]

Thai_w24 <- Thai_w24_read


############################################################
## Data - visit week 36

Thai_w36_read  <- read.csv("C://U//GHIT//1_AlphaScreen_Data//2_PLOS_NTDS_manuscript//Analysis//Thai_Alpha_w36_adj_reps.csv")

Thai_w36_read[384+(1:384),1]   <- Thai_w36_read[1:384,1]
Thai_w36_read[2*384+(1:384),1] <- Thai_w36_read[1:384,1]

Thai_w36 <- Thai_w36_read




Thai_all <- rbind( Thai_w0, Thai_w12, Thai_w24, Thai_w36 )


Thai_antig_ID = matrix(NA, nrow=12*length(antig_ID), ncol=ncol(Thai_all))
colnames(Thai_antig_ID) = colnames(Thai_all)

Thai_antig_ID = as.data.frame(Thai_antig_ID)

Thai_antig_ID[1:(4*length(antig_ID)),34]                       = 1
Thai_antig_ID[(4*length(antig_ID)+1):(8*length(antig_ID)),34]  = 2
Thai_antig_ID[(8*length(antig_ID)+1):(12*length(antig_ID)),34] = 3

Thai_antig_ID[1:(4*length(antig_ID)),35]                       = c( rep("w0", length(antig_ID)), rep("w12", length(antig_ID)), rep("w24", length(antig_ID)), rep("w36", length(antig_ID)) )
Thai_antig_ID[(4*length(antig_ID)+1):(8*length(antig_ID)),35]  = c( rep("w0", length(antig_ID)), rep("w12", length(antig_ID)), rep("w24", length(antig_ID)), rep("w36", length(antig_ID)) )
Thai_antig_ID[(8*length(antig_ID)+1):(12*length(antig_ID)),35] = c( rep("w0", length(antig_ID)), rep("w12", length(antig_ID)), rep("w24", length(antig_ID)), rep("w36", length(antig_ID)) )


Thai_antig_ID[1:(4*length(antig_ID)),1]                       = c( antig_ID, antig_ID, antig_ID, antig_ID )
Thai_antig_ID[(4*length(antig_ID)+1):(8*length(antig_ID)),1]  = c( antig_ID, antig_ID, antig_ID, antig_ID )
Thai_antig_ID[(8*length(antig_ID)+1):(12*length(antig_ID)),1] = c( antig_ID, antig_ID, antig_ID, antig_ID )


for(i in 1:nrow(Thai_antig_ID))
{
	index = intersect( which(Thai_all[,1] == Thai_antig_ID[i,1]), intersect( which(Thai_all[,34] == Thai_antig_ID[i,34]), which(Thai_all[,35] == Thai_antig_ID[i,35]) ) ) 	

	if( length(index) > 0 )
	{
		Thai_antig_ID[i,2:33] = Thai_all[index,2:33]
	}
}




t_sample <- 7*c(0, 12, 24, 36)
names(t_sample) <- c("w0", "w12", "w24", "w36")


thai_AB_data <- list()

for(n in 1:N_ant)
{
	temp <- Thai_antig_ID[ which(Thai_antig_ID[,1]==antig_ID[n]), ]

	N_rep <- length(unique(temp[,34]))


	AB_mat <- matrix(NA, nrow=N_thai_part, ncol=0)
	rownames(AB_mat) <- thai_part_ID


	for(k in 1:N_rep)
	{
		temp_rep <- temp[which(temp[,34]==k),]
	
		temp_add <- matrix(NA, nrow=N_thai_part, ncol=4)
		rownames(temp_add) <- thai_part_ID
		colnames(temp_add) <- names(t_sample)

		for(i in 1:N_thai_part)
		{
			for(j in 1:4)
			{
				temp_add[i,j] <- temp_rep[ which(temp_rep[,35]==colnames(temp_add)[j]),
                                                   which(colnames(temp_rep)==rownames(temp_add)[i]) ]
			}
		}
	
		AB_mat <- cbind(AB_mat, temp_add)
	}



	thai_AB_data[[n]] <- AB_mat
}



Prop_zero_first <- rep(NA, N_ant)

for(n in 1:N_ant)
{
	N_rep_n <- length(which(colnames(thai_AB_data[[n]])=="w0"))

	Prop_zero_first[n] <- length(which(thai_AB_data[[n]][,which(colnames(thai_AB_data[[n]])=="w0")]<1e-6))/( N_thai_part*N_rep_n )
}



Prop_zero <- rep(NA, N_ant)

for(n in 1:N_ant)
{
	N_rep_n <- length(which(colnames(thai_AB_data[[n]])=="w0"))

	Prop_zero[n] <- length(which(thai_AB_data[[n]]<1e-6))/( 4*N_thai_part*N_rep_n )

	if( Prop_zero[n] > 0 )
	{
		min_AB <- min( thai_AB_data[[n]][ which(thai_AB_data[[n]]>1e-6,arr.ind=TRUE) ], na.rm=TRUE )		
	
		thai_AB_data[[n]][ which(thai_AB_data[[n]]<1e-6,arr.ind=TRUE) ] <- 0.5*min_AB
	}
}	





##################################################
## Mixed effects model assuming antibody titres for each
## participant are described by local parameters:
##
## AB_0: initial antibody titre
## d_AB: half-life of antibody
##
## These local parameters then inform global parameters
## describing the average behaviour of the population




thai_AB_pars <- matrix(NA, nrow=N_ant, ncol=8)
rownames(thai_AB_pars) <- antig_ID
colnames(thai_AB_pars) <- c("AB_GMT", "sigma_A", "r_mean", "sigma_r", "sigma_m", "d_AB_est", "r_mean_low", "r_mean_high")



thai_AB_starts <- matrix(NA, nrow=N_ant, ncol=N_thai_part)
rownames(thai_AB_starts) <- antig_ID
colnames(thai_AB_starts) <- thai_part_ID


thai_AB_rates <- matrix(NA, nrow=N_ant, ncol=N_thai_part)
rownames(thai_AB_rates) <- antig_ID
colnames(thai_AB_rates) <- thai_part_ID

FLAG <- rep(NA, length=N_ant)

for(n in 1:N_ant)
{

	if( Prop_zero[n]<1 )
	{
		######################################
		## Put together dataframe of all data

		N_rep_n <- length(which(colnames(thai_AB_data[[n]])=="w0"))

		thai_AB_data_all <- as.vector( t(thai_AB_data[[n]]) )

		t_sample_all <- rep( t_sample, N_thai_part*N_rep_n )
		t_sample_all <- as.vector(t_sample_all)
	
		ID_all <- rep( 1, N_rep_n*length(t_sample))
		for(i in 2:N_thai_part)
		{
			ID_all <- c( ID_all, rep(i, N_rep_n*length(t_sample) ) )
		}


		thai_AB_data_frame <- cbind( thai_AB_data_all, t_sample_all, ID_all )
		thai_AB_data_frame <- as.data.frame( thai_AB_data_frame )


		######################################
		## Fit mixed effects linear model and 
		## extract parameters using lmer

		tryCatch(
		{
			mod_3_lmer <- lmer( log(thai_AB_data_all) ~ t_sample_all + (1+t_sample_all|ID_all),
            			     	  data=thai_AB_data_frame )

			thai_AB_pars[n,1] <- summary(mod_3_lmer)[[10]][1,1]
			thai_AB_pars[n,3] <- summary(mod_3_lmer)[[10]][2,1]

			thai_AB_pars[n,2] <- sqrt(summary( mod_3_lmer )[[13]]$ID_all[1,1])	
			thai_AB_pars[n,4] <- sqrt(summary( mod_3_lmer )[[13]]$ID_all[2,2])

			thai_AB_pars[n,5] <- summary(mod_3_lmer)[[11]]

			thai_AB_pars[n,6] <- - log(2)/thai_AB_pars[n,3]

			thai_AB_pars[n,7] <- summary(mod_3_lmer)[[10]][2,1] - 1.96*summary(mod_3_lmer)[[10]][2,2]
			thai_AB_pars[n,8] <- summary(mod_3_lmer)[[10]][2,1] + 1.96*summary(mod_3_lmer)[[10]][2,2]

			thai_AB_starts[n,] <- summary(mod_3_lmer)[[10]][1,1] + ranef(mod_3_lmer)$ID_all[,1]
			thai_AB_rates[n,]  <- summary(mod_3_lmer)[[10]][2,1] + ranef(mod_3_lmer)$ID_all[,2]


		}, error=function(e){ NULL }
		)

	}
}



####################
####################
####################
###              ### 
###   BRAZIL     ###
###              ###
####################
####################
####################


N_braz_part = 33
N_sam  = 4
N_rep  = 3


############################################################
## Data

Braz_read  <- read.csv("C://U//GHIT//1_AlphaScreen_Data//2_PLOS_NTDS_manuscript//Analysis//Brazil_Alpha_adj_reps.csv")

Braz_read[384+(1:384),1]   <- Braz_read[1:384,1]
Braz_read[2*384+(1:384),1] <- Braz_read[1:384,1]

##Braz_read <- Braz_read[which( Braz_read[,1] %in% proteins[,6] ),]


sample <- rep(NA, ncol(Braz_read)-2)

for(j in 1:length(sample))
{
	if( substr(colnames(Braz_read)[1+j], 1, 3) == "V01" )
	{
		sample[j] = "w0"
	}

	if( substr(colnames(Braz_read)[1+j], 1, 3) == "V08" )
	{
		sample[j] = "w12"
	}

	if( substr(colnames(Braz_read)[1+j], 1, 3) == "V14" )
	{
		sample[j] = "w24"
	}

	if( substr(colnames(Braz_read)[1+j], 1, 3) %in% c("V16","V17") )
	{
		sample[j] = "w36"
	}
}



braz_part_ID <- rep(NA, ncol(Braz_read)-2)

for(j in 1:length(braz_part_ID))
{
	braz_part_ID[j] <- substr(colnames(Braz_read)[1+j], 4, nchar(colnames(Braz_read)[1+j]))
}

braz_part_ID <- braz_part_ID[seq(from=1, by=4, length=33)]




Braz_w0 <- Braz_read[,c(1,1+which(sample=="w0"),ncol(Braz_read))]
Braz_w0 <- cbind( Braz_w0, rep("w0", nrow(Braz_w0)) )
colnames(Braz_w0) <- c("ant_names", braz_part_ID, "replicate", "visit")


Braz_w12 <- Braz_read[,c(1,1+which(sample=="w12"),ncol(Braz_read))]
Braz_w12 <- cbind( Braz_w12, rep("w12", nrow(Braz_w12)) )
colnames(Braz_w12) <- c("ant_names", braz_part_ID, "replicate", "visit")


Braz_w24 <- Braz_read[,c(1,1+which(sample=="w24"),ncol(Braz_read))]
Braz_w24 <- cbind( Braz_w24, rep("w24", nrow(Braz_w24)) )
colnames(Braz_w24) <- c("ant_names", braz_part_ID, "replicate", "visit")


Braz_w36 <- Braz_read[,c(1,1+which(sample=="w36"),ncol(Braz_read))]
Braz_w36 <- cbind( Braz_w36, rep("w36", nrow(Braz_w36)) )
colnames(Braz_w36) <- c("ant_names", braz_part_ID, "replicate", "visit")


Braz_all <- rbind( Braz_w0, Braz_w12, Braz_w24, Braz_w36 )




Braz_antig_ID = matrix(NA, nrow=12*length(antig_ID), ncol=ncol(Braz_all))
colnames(Braz_antig_ID) = colnames(Braz_all)

Braz_antig_ID = as.data.frame(Braz_antig_ID)

Braz_antig_ID[1:(4*length(antig_ID)),35]                       = 1
Braz_antig_ID[(4*length(antig_ID)+1):(8*length(antig_ID)),35]  = 2
Braz_antig_ID[(8*length(antig_ID)+1):(12*length(antig_ID)),35] = 3

Braz_antig_ID[1:(4*length(antig_ID)),36]                       = c( rep("w0", length(antig_ID)), rep("w12", length(antig_ID)), rep("w24", length(antig_ID)), rep("w36", length(antig_ID)) )
Braz_antig_ID[(4*length(antig_ID)+1):(8*length(antig_ID)),36]  = c( rep("w0", length(antig_ID)), rep("w12", length(antig_ID)), rep("w24", length(antig_ID)), rep("w36", length(antig_ID)) )
Braz_antig_ID[(8*length(antig_ID)+1):(12*length(antig_ID)),36] = c( rep("w0", length(antig_ID)), rep("w12", length(antig_ID)), rep("w24", length(antig_ID)), rep("w36", length(antig_ID)) )


Braz_antig_ID[1:(4*length(antig_ID)),1]                       = c( antig_ID, antig_ID, antig_ID, antig_ID )
Braz_antig_ID[(4*length(antig_ID)+1):(8*length(antig_ID)),1]  = c( antig_ID, antig_ID, antig_ID, antig_ID )
Braz_antig_ID[(8*length(antig_ID)+1):(12*length(antig_ID)),1] = c( antig_ID, antig_ID, antig_ID, antig_ID )


for(i in 1:nrow(Braz_antig_ID))
{
	index = intersect( which(Braz_all[,1] == Braz_antig_ID[i,1]), intersect( which(Braz_all[,35] == Braz_antig_ID[i,35]), which(Braz_all[,36] == Braz_antig_ID[i,36]) ) ) 	

	if( length(index) > 0 )
	{
		Braz_antig_ID[i,2:34] = Braz_all[index,2:34]
	}
}




t_sample <- 7*c(0, 12, 24, 36)
names(t_sample) <- c("w0", "w12", "w24", "w36")


braz_AB_data <- list()

for(n in 1:N_ant)
{
	temp <- Braz_antig_ID[ which(Braz_antig_ID[,1]==antig_ID[n]), ]

	N_rep <- length(unique(temp[,35]))


	AB_mat <- matrix(NA, nrow=N_braz_part, ncol=0)
	rownames(AB_mat) <- braz_part_ID


	for(k in 1:N_rep)
	{
		temp_rep <- temp[which(temp[,35]==k),]
	
		temp_add <- matrix(NA, nrow=N_braz_part, ncol=4)
		rownames(temp_add) <- braz_part_ID
		colnames(temp_add) <- names(t_sample)

		for(i in 1:N_braz_part)
		{
			for(j in 1:4)
			{
				temp_add[i,j] <- temp_rep[ which(temp_rep[,36]==colnames(temp_add)[j]),
                                                   which(colnames(temp_rep)==rownames(temp_add)[i]) ]
			}
		}
	
		AB_mat <- cbind(AB_mat, temp_add)
	}



	braz_AB_data[[n]] <- AB_mat
}



Prop_zero_first <- rep(NA, N_ant)

for(n in 1:N_ant)
{
	N_rep_n <- length(which(colnames(braz_AB_data[[n]])=="w0"))

	Prop_zero_first[n] <- length(which(braz_AB_data[[n]][,which(colnames(braz_AB_data[[n]])=="w0")]<1e-6))/( N_braz_part*N_rep_n )
}



Prop_zero <- rep(NA, N_ant)

for(n in 1:N_ant)
{
	N_rep_n <- length(which(colnames(braz_AB_data[[n]])=="w0"))

	Prop_zero[n] <- length(which(braz_AB_data[[n]]<1e-6))/( 4*N_braz_part*N_rep_n )

	if( Prop_zero[n] > 0 )
	{
		min_AB <- min( braz_AB_data[[n]][ which(braz_AB_data[[n]]>1e-6,arr.ind=TRUE) ], na.rm=TRUE )		
	
		braz_AB_data[[n]][ which(braz_AB_data[[n]]<1e-6,arr.ind=TRUE) ] <- 0.5*min_AB
	}
}	








##################################################
## Mixed effects model assuming antibody titres for each
## participant are described by local parameters:
##
## AB_0: initial antibody titre
## d_AB: half-life of antibody
##
## These local parameters then inform global parameters
## describing the average behaviour of the population


braz_AB_pars <- matrix(NA, nrow=N_ant, ncol=8)
rownames(braz_AB_pars) <- antig_ID
colnames(braz_AB_pars) <- c("AB_GMT", "sigma_A", "r_mean", "sigma_r", "sigma_m", "d_AB_est", "r_mean_low", "r_mean_high")



braz_AB_starts <- matrix(NA, nrow=N_ant, ncol=N_braz_part)
rownames(braz_AB_starts) <- antig_ID
colnames(braz_AB_starts) <- braz_part_ID


braz_AB_rates <- matrix(NA, nrow=N_ant, ncol=N_braz_part)
rownames(braz_AB_rates) <- antig_ID
colnames(braz_AB_rates) <- braz_part_ID

FLAG <- rep(NA, length=N_ant)

for(n in 1:N_ant)
{

	if( Prop_zero[n]<1 )
	{
		######################################
		## Put together dataframe of all data

		N_rep_n <- length(which(colnames(braz_AB_data[[n]])=="w0"))

		braz_AB_data_all <- as.vector( t(braz_AB_data[[n]]) )

		t_sample_all <- rep( t_sample, N_braz_part*N_rep_n )
		t_sample_all <- as.vector(t_sample_all)
	
		ID_all <- rep( 1, N_rep_n*length(t_sample))
		for(i in 2:N_braz_part)
		{
			ID_all <- c( ID_all, rep(i, N_rep_n*length(t_sample) ) )
		}


		braz_AB_data_frame <- cbind( braz_AB_data_all, t_sample_all, ID_all )
		braz_AB_data_frame <- as.data.frame( braz_AB_data_frame )


		######################################
		## Fit mixed effects linear model and 
		## extract parameters using lmer

		tryCatch(
		{
			mod_3_lmer <- lmer( log(braz_AB_data_all) ~ t_sample_all + (1+t_sample_all|ID_all),
            			     	  data=braz_AB_data_frame )

			braz_AB_pars[n,1] <- summary(mod_3_lmer)[[10]][1,1]
			braz_AB_pars[n,3] <- summary(mod_3_lmer)[[10]][2,1]

			braz_AB_pars[n,2] <- sqrt(summary( mod_3_lmer )[[13]]$ID_all[1,1])	
			braz_AB_pars[n,4] <- sqrt(summary( mod_3_lmer )[[13]]$ID_all[2,2])

			braz_AB_pars[n,5] <- summary(mod_3_lmer)[[11]]

			braz_AB_pars[n,6] <- - log(2)/braz_AB_pars[n,3]

			braz_AB_pars[n,7] <- summary(mod_3_lmer)[[10]][2,1] - 1.96*summary(mod_3_lmer)[[10]][2,2]
			braz_AB_pars[n,8] <- summary(mod_3_lmer)[[10]][2,1] + 1.96*summary(mod_3_lmer)[[10]][2,2]

			braz_AB_starts[n,] <- summary(mod_3_lmer)[[10]][1,1] + ranef(mod_3_lmer)$ID_all[,1]
			braz_AB_rates[n,]  <- summary(mod_3_lmer)[[10]][2,1] + ranef(mod_3_lmer)$ID_all[,2]


		}, error=function(e){ NULL }
		)

	}
}






#########################################
#########################################
##                                     ##
##  #     #  ####  ####   ##### ##     ##
##  ##   ## ##  ## ## ##  ##    ##     ##
##  ####### ##  ## ##  ## ####  ##     ##
##  ## # ## ##  ## ## ##  ##    ##     ##
##  ##   ##  ####  ####   ##### #####  ##
##                                     ##
######################################### 
#########################################

bb_5_index = order(ant_prob[,5], decreasing=TRUE)[1:5]

bb_5 = rep(0, N_ant)
bb_5[bb_5_index] = 1



####################
####################
####################
###              ###
###   Thailand   ###
###              ###
####################
####################

N_thai_part <- 32


thai_A_meas_0_rep1 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//thai_A_meas_0_rep1_DS.csv")

thai_A_meas_0_rep2 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//thai_A_meas_0_rep2_DS.csv")

thai_A_meas_12_rep1 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//thai_A_meas_12_rep1_DS.csv")

thai_A_meas_12_rep2 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//thai_A_meas_12_rep2_DS.csv")

thai_A_meas_24_rep1 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//thai_A_meas_24_rep1_DS.csv")

thai_A_meas_24_rep2 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//thai_A_meas_24_rep2_DS.csv")

thai_A_meas_36_rep1 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//thai_A_meas_36_rep1_DS.csv")

thai_A_meas_36_rep2 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//thai_A_meas_36_rep2_DS.csv")

thai_A_meas_neg <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//thai_A_meas_neg_DS.csv")

thai_A_meas_2y  <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//thai_A_meas_2y_DS.csv")


thai_A_meas_0_rep1 <- t(log(thai_A_meas_0_rep1[,2:(N_thai_part+1)]))
thai_A_meas_0_rep2 <- t(log(thai_A_meas_0_rep2[,2:(N_thai_part+1)]))

thai_A_meas_12_rep1 <- t(log(thai_A_meas_12_rep1[,2:(N_thai_part+1)]))
thai_A_meas_12_rep2 <- t(log(thai_A_meas_12_rep2[,2:(N_thai_part+1)]))

thai_A_meas_24_rep1 <- t(log(thai_A_meas_24_rep1[,2:(N_thai_part+1)]))
thai_A_meas_24_rep2 <- t(log(thai_A_meas_24_rep2[,2:(N_thai_part+1)]))

thai_A_meas_36_rep1 <- t(log(thai_A_meas_36_rep1[,2:(N_thai_part+1)]))
thai_A_meas_36_rep2 <- t(log(thai_A_meas_36_rep2[,2:(N_thai_part+1)]))

thai_A_meas_neg <- t(log(thai_A_meas_neg[,2:(N_thai_part+1)]))

thai_A_meas_2y <- t(log(thai_A_meas_2y[,2:(N_thai_part+1)]))



t_LL_seq <- seq(from=1, to=5000, by=1)
N_t <- length(t_LL_seq)



Thai_loglike <- function( bb, repx )
{
	#######################################
	## Create index of selected antigens
	## and create objects of interest

	index <- bb*(1:N_ant)
	index <- index[ which(index>0.5) ]
	index <- as.numeric(index)

	N_bb <- length(index)
	log_2pi_bb <- - 0.5*N_bb*log(2*pi)

	A_GMT_bb <- thai_A_GMT_r_mean_read[index,1]
	rr_bb    <- thai_A_GMT_r_mean_read[index,2]

	SIGMA_A_bb <- thai_SIGMA_A[index,index]
	SIGMA_r_bb <- thai_SIGMA_r[index,index]
	SIGMA_m_bb <- thai_SIGMA_m[index,index]

	if( repx==1 )
	{
		A_meas_12_bb  <- thai_A_meas_12_rep1[,index] 
		A_meas_24_bb  <- thai_A_meas_24_rep1[,index] 
		A_meas_36_bb  <- thai_A_meas_36_rep1[,index] 
	}	

	if( repx==2 )
	{
		A_meas_12_bb  <- thai_A_meas_12_rep2[,index] 
		A_meas_24_bb  <- thai_A_meas_24_rep2[,index] 
		A_meas_36_bb  <- thai_A_meas_36_rep2[,index] 
	}
	
	A_meas_2y_bb  <- thai_A_meas_2y[,index] 
	A_meas_neg_bb <- thai_A_meas_neg[,index] 

	loglike_t_12  <- matrix(NA, nrow=N_thai_part, ncol=N_t)
	loglike_t_24  <- matrix(NA, nrow=N_thai_part, ncol=N_t)
	loglike_t_36  <- matrix(NA, nrow=N_thai_part, ncol=N_t)
	loglike_t_2y  <- matrix(NA, nrow=N_thai_part, ncol=N_t)
	loglike_t_neg <- matrix(NA, nrow=N_thai_part, ncol=N_t)

	for(j in 1:N_t)
	{
		SIGMA_total <- SIGMA_A_bb + SIGMA_r_bb*t_LL_seq[j]^2 + SIGMA_m_bb
			
		log_det_SIGMA_total <- log( det(SIGMA_total) )
		SIGMA_total_inv     <- solve( SIGMA_total )

		A_pred <- A_GMT_bb + rr_bb*t_LL_seq[j]

		##############################
		## Loop through each of the participants

		for(i in 1:N_thai_part)
		{
			xx_12  <- A_pred - A_meas_12_bb[i,]
			xx_24  <- A_pred - A_meas_24_bb[i,]
			xx_36  <- A_pred - A_meas_36_bb[i,]
			xx_2y  <- A_pred - A_meas_2y_bb[i,]
			xx_neg <- A_pred - A_meas_neg_bb[i,]

			loglike_t_12[i,j]  = log_2pi_bb - 0.5*log_det_SIGMA_total - 0.5*xx_12%*%SIGMA_total_inv%*%xx_12
			loglike_t_24[i,j]  = log_2pi_bb - 0.5*log_det_SIGMA_total - 0.5*xx_24%*%SIGMA_total_inv%*%xx_24
			loglike_t_36[i,j]  = log_2pi_bb - 0.5*log_det_SIGMA_total - 0.5*xx_36%*%SIGMA_total_inv%*%xx_36
			loglike_t_2y[i,j]  = log_2pi_bb - 0.5*log_det_SIGMA_total - 0.5*xx_2y%*%SIGMA_total_inv%*%xx_2y
			loglike_t_neg[i,j] = log_2pi_bb - 0.5*log_det_SIGMA_total - 0.5*xx_neg%*%SIGMA_total_inv%*%xx_neg
		}
	}

	OUTPUT = list()

	OUTPUT[[1]] = loglike_t_12
	OUTPUT[[2]] = loglike_t_24
	OUTPUT[[3]] = loglike_t_36
	OUTPUT[[4]] = loglike_t_2y
	OUTPUT[[5]] = loglike_t_neg
	
	OUTPUT
}



####################
####################
####################
###              ###
###   Brazil     ###
###              ###
####################
####################



N_braz_part <- 33


braz_A_meas_0_rep1 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//braz_A_meas_0_rep1_DS.csv")

braz_A_meas_0_rep2 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//braz_A_meas_0_rep2_DS.csv")

braz_A_meas_12_rep1 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//braz_A_meas_12_rep1_DS.csv")

braz_A_meas_12_rep2 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//braz_A_meas_12_rep2_DS.csv")

braz_A_meas_24_rep1 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//braz_A_meas_24_rep1_DS.csv")

braz_A_meas_24_rep2 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//braz_A_meas_24_rep2_DS.csv")

braz_A_meas_36_rep1 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//braz_A_meas_36_rep1_DS.csv")

braz_A_meas_36_rep2 <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//braz_A_meas_36_rep2_DS.csv")

braz_A_meas_neg <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//braz_A_meas_neg_DS.csv")

braz_A_meas_2y  <- read.csv("C://U//GHIT//1_AlphaScreen_Data//1_AlphaScreen_Down_Selection//Round3//Processed_Data//braz_A_meas_2y_DS.csv")


braz_A_meas_0_rep1 <- t(log(braz_A_meas_0_rep1[,2:(N_braz_part+1)]))
braz_A_meas_0_rep2 <- t(log(braz_A_meas_0_rep2[,2:(N_braz_part)]))

braz_A_meas_12_rep1 <- t(log(braz_A_meas_12_rep1[,2:(N_braz_part+1)]))
braz_A_meas_12_rep2 <- t(log(braz_A_meas_12_rep2[,2:(N_braz_part+1)]))

braz_A_meas_24_rep1 <- t(log(braz_A_meas_24_rep1[,2:(N_braz_part+1)]))
braz_A_meas_24_rep2 <- t(log(braz_A_meas_24_rep2[,2:(N_braz_part+1)]))

braz_A_meas_36_rep1 <- t(log(braz_A_meas_36_rep1[,2:(N_braz_part+1)]))
braz_A_meas_36_rep2 <- t(log(braz_A_meas_36_rep2[,2:(N_braz_part+1)]))

braz_A_meas_neg <- t(log(braz_A_meas_neg[,2:(N_braz_part+1)]))

braz_A_meas_2y <- t(log(braz_A_meas_2y[,2:(N_braz_part+1)]))





Braz_loglike <- function( bb, repx )
{
	#######################################
	## Create index of selected antigens
	## and create objects of interest

	index <- bb*(1:N_ant)
	index <- index[ which(index>0.5) ]
	index <- as.numeric(index)

	N_bb <- length(index)
	log_2pi_bb <- - 0.5*N_bb*log(2*pi)

	A_GMT_bb <- braz_A_GMT_r_mean_read[index,1]
	rr_bb    <- braz_A_GMT_r_mean_read[index,2]

	SIGMA_A_bb <- braz_SIGMA_A[index,index]
	SIGMA_r_bb <- braz_SIGMA_r[index,index]
	SIGMA_m_bb <- braz_SIGMA_m[index,index]

	if( repx==1 )
	{
		A_meas_12_bb  <- braz_A_meas_12_rep1[,index] 
		A_meas_24_bb  <- braz_A_meas_24_rep1[,index] 
		A_meas_36_bb  <- braz_A_meas_36_rep1[,index] 
	}	

	if( repx==2 )
	{
		A_meas_12_bb  <- braz_A_meas_12_rep2[,index] 
		A_meas_24_bb  <- braz_A_meas_24_rep2[,index] 
		A_meas_36_bb  <- braz_A_meas_36_rep2[,index] 
	}
	
	A_meas_2y_bb  <- braz_A_meas_2y[,index] 
	A_meas_neg_bb <- braz_A_meas_neg[,index] 

	loglike_t_12  <- matrix(NA, nrow=N_braz_part, ncol=N_t)
	loglike_t_24  <- matrix(NA, nrow=N_braz_part, ncol=N_t)
	loglike_t_36  <- matrix(NA, nrow=N_braz_part, ncol=N_t)
	loglike_t_2y  <- matrix(NA, nrow=N_braz_part, ncol=N_t)
	loglike_t_neg <- matrix(NA, nrow=N_braz_part, ncol=N_t)

	for(j in 1:N_t)
	{
		SIGMA_total <- SIGMA_A_bb + SIGMA_r_bb*t_LL_seq[j]^2 + SIGMA_m_bb
			
		log_det_SIGMA_total <- log( det(SIGMA_total) )
		SIGMA_total_inv     <- solve( SIGMA_total )

		A_pred <- A_GMT_bb + rr_bb*t_LL_seq[j]

		##############################
		## Loop through each of the participants

		for(i in 1:N_thai_part)
		{
			xx_12  <- A_pred - A_meas_12_bb[i,]
			xx_24  <- A_pred - A_meas_24_bb[i,]
			xx_36  <- A_pred - A_meas_36_bb[i,]
			xx_2y  <- A_pred - A_meas_2y_bb[i,]
			xx_neg <- A_pred - A_meas_neg_bb[i,]

			loglike_t_12[i,j]  = log_2pi_bb - 0.5*log_det_SIGMA_total - 0.5*xx_12%*%SIGMA_total_inv%*%xx_12
			loglike_t_24[i,j]  = log_2pi_bb - 0.5*log_det_SIGMA_total - 0.5*xx_24%*%SIGMA_total_inv%*%xx_24
			loglike_t_36[i,j]  = log_2pi_bb - 0.5*log_det_SIGMA_total - 0.5*xx_36%*%SIGMA_total_inv%*%xx_36
			loglike_t_2y[i,j]  = log_2pi_bb - 0.5*log_det_SIGMA_total - 0.5*xx_2y%*%SIGMA_total_inv%*%xx_2y
			loglike_t_neg[i,j] = log_2pi_bb - 0.5*log_det_SIGMA_total - 0.5*xx_neg%*%SIGMA_total_inv%*%xx_neg
		}
	}

	OUTPUT = list()

	OUTPUT[[1]] = loglike_t_12
	OUTPUT[[2]] = loglike_t_24
	OUTPUT[[3]] = loglike_t_36
	OUTPUT[[4]] = loglike_t_2y
	OUTPUT[[5]] = loglike_t_neg
	
	OUTPUT
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

NAMES = read.csv("C://U//GHIT//2_All_Data//2_Data//Antigen_DownSelection_Files//1_PLOS_NTDS_combined_results.csv")

NAMES = NAMES[,c(2,5)]

well_names = rownames(thai_A_GMT_r_mean_read)

PVX_names = rep(NA, length(well_names))

for(i in 1:length(PVX_names))
{
	if( well_names[i] %in% NAMES[,2] )
	{
		PVX_names[i] = as.vector(NAMES[which( NAMES[,2] == well_names[i] ),1]) 
	}else{
		PVX_names[i] = "PVX_xxxxx"
	}
}

PVX_names[which( well_names == "342" )] = "PVX_002550"
PVX_names[which( well_names == "314" )] = "PVX_081560"
PVX_names[which( well_names == "312" )] = "PVX_118040"
PVX_names[which( well_names == "313" )] = "PVX_089585"
PVX_names[which( well_names == "337" )] = "PVX_001000"
PVX_names[which( well_names == "325" )] = "PVX_112625"


Thai_ID = 4
Braz_ID = 5




ant_5_cols = rainbow(5)

ant_prob_cols = rep("black", length(bb_5))
for(k in 1:5)
{
	ant_prob_cols[bb_5_index[k]] = ant_5_cols[k]
}


tiff(file="Figure3_Antigen_Discovery.tif", width=20, height=22, units="cm", res=500)
	


par(ask=FALSE)
main.size = 1.1
lab.size = 0.8
axis.size = 0.65
point.size = 0.8
line.size = 1
A.size = 1.25


#lay.mat <- rbind( c(1,1,1,1,2,2,2,2),
#                  c(3,3,3,3,4,4,4,4),
#			c(5,5,5,5,5,5,5,6),
#			c(7,7,7,7,7,7,7,8),
#			c(9,9,9,9,9,9,9,10) )
#layout(lay.mat, heights=c(1,1,2.5,0.75,0.75))
#layout.show(10)

lay.mat <- rbind( c(1,1,1,1,1,2,2,2,2,2),
                  c(3,3,3,3,3,4,4,4,4,4),
			c(5,5,5,5,5,5,5,5,5,6) )
layout(lay.mat, heights=c(1,1,2.5))
layout.show(6)

par(mar=c(3,3,2,1))
par(mgp=c(1.5,0.75,0))


###########################################
##                                       ##
##  PANEL 1                              ## 
##  Individual-level dynamics: Thailand  ##
##                                       ##
###########################################

tt_data = rep( 7*c(0, 12, 24, 36), 3)

tt_seq = seq(from=0, to=365, by=1)


par(mar=c(3,3,2,1))
par(mgp=c(1.5,0.75,0))

plot(x=1e10, y=1e10, log="y",
pch=19, 
xlim=c(0,365), ylim=c(0.005,100),
main="(A) Kinetics of 5 antigens: individual Thai_1",   
xlab="time (months)", ylab="antibody titre",
xaxt='n', yaxt='n', bty='n',
#xaxs='i',
yaxs='i',
cex.main=1.25, cex.lab=1.25, cex.axis=1.25)



for(k in 1:5)
{
	points( x=tt_data, y=thai_AB_data[[bb_5_index[k]]][Thai_ID,], 
              pch=19, col=ant_5_cols[k] )
}

for(k in 1:5)
{
	points(x=tt_seq, y=exp(thai_AB_starts[bb_5_index[k],Thai_ID] + thai_AB_rates[bb_5_index[k],Thai_ID]*tt_seq), 
 		 type='l', col=ant_5_cols[k] )
}


axis(1, at=c(0,84,168,252,365), label=c("0","3","6","9","12"), cex.axis=1.5*axis.size)
axis(2, at=c(0.01,0.1,1,10,100), label=c(0.01,0.1,1,10,100), cex.axis=axis.size, las=2)





###########################################
##                                       ##
##  PANEL 2                              ## 
##  Individual-level dynamics: Brazil    ##
##                                       ##
###########################################

tt_data = rep( 7*c(0, 12, 24, 36), 3)

tt_seq = seq(from=0, to=365, by=1)

plot(x=1e10, y=1e10, log="y",
pch=19, 
xlim=c(0,365), ylim=c(0.005,100),
main="(B) Kinetics of 5 antigens: individual Braz_1",  
xlab="time (months)", ylab="antibody titre",
xaxt='n', yaxt='n', bty='n',
#xaxs='i',
yaxs='i',
cex.main=1.25, cex.lab=1.25, cex.axis=1.25)




for(k in 1:5)
{
	points( x=tt_data, y=braz_AB_data[[bb_5_index[k]]][Braz_ID,], 
              pch=19, col=rainbow(5)[k] )
}

for(k in 1:5)
{
	points(x=tt_seq, y=exp(braz_AB_starts[bb_5_index[k],Braz_ID] + braz_AB_rates[bb_5_index[k],Braz_ID]*tt_seq), 
 		 type='l', col=rainbow(5)[k] )
}


axis(1, at=c(0,84,168,252,365), label=c("0","3","6","9","12"), cex.axis=1.5*axis.size)
axis(2, at=c(0.01,0.1,1,10,100), label=c(0.01,0.1,1,10,100), cex.axis=axis.size, las=2)



#################################################
##                                             ##
##  PANEL 3                                    ## 
##  Predicting time since infection: Thailand  ##
##                                             ##
#################################################


Thai_pred = Thai_loglike( bb_5, 1 )

Thai_LL_12 = exp(Thai_pred[[1]][Thai_ID,])
Thai_LL_12 = Thai_LL_12/sum(Thai_LL_12)

Thai_LL_24 = exp(Thai_pred[[2]][Thai_ID,])
Thai_LL_24 = Thai_LL_24/sum(Thai_LL_24)

Thai_LL_36 = exp(Thai_pred[[3]][Thai_ID,])
Thai_LL_36 = Thai_LL_36/sum(Thai_LL_36)





plot(x=t_LL_seq, y=Thai_LL_24, 
type='l', lwd=2,
xlim=c(0,18*30), ylim=c(0,1.1*max(Thai_LL_24)),
main="(C) Estimating infection time: individual Thai_1", 
xlab="time (months)", ylab="probability density",
xaxt='n', yaxt='n', bty='n',
#xaxs='i',
yaxs='i',
cex.main=1.25, cex.lab=1.25, cex.axis=1.25)


points(x=t_LL_seq[which.max(Thai_LL_24)], y=max(Thai_LL_24), pch=15, cex=2, col="black")

cut_95 = exp( log(max(Thai_LL_24)) - 1.92 )

low_95 = t_LL_seq[which.min(abs(Thai_LL_24[1:which.max(Thai_LL_24)] - cut_95))]
high_95 = t_LL_seq[which.min(abs(Thai_LL_24[which.max(Thai_LL_24):length(Thai_LL_24)] - cut_95))]


points(x=rep( low_95, 2), y=c(0,1000), type='l', lty="dashed", col="black" )

points(x=rep( high_95, 2), y=c(0,1000), type='l', lty="dashed", col="black" )



axis(1, at=c(0,90,180,270,360,450,540), label=c("0","3","6","9","12","15","18"), cex.axis=1.5*axis.size)
axis(2, at=c(0,1.1*max(Thai_LL_24)), label=c("", ""), cex.axis=axis.size)




#################################################
##                                             ##
##  PANEL 4                                    ## 
##  Predicting time since infection: Brazil    ##
##                                             ##
#################################################


Braz_pred = Braz_loglike( bb_5, 1 )

Braz_LL_12 = exp(Braz_pred[[1]][Braz_ID,])
Braz_LL_12 = Braz_LL_12/sum(Braz_LL_12)

Braz_LL_24 = exp(Braz_pred[[2]][Braz_ID,])
Braz_LL_24 = Braz_LL_24/sum(Braz_LL_24)

Braz_LL_36 = exp(Braz_pred[[3]][Braz_ID,])
Braz_LL_36 = Braz_LL_36/sum(Braz_LL_36)




par(mar=c(3,3,2,1))
par(mgp=c(1.5,0.75,0))

plot(x=t_LL_seq, y=Braz_LL_24, 
type='l', lwd=2,
xlim=c(0,18*30), ylim=c(0,1.1*max(Braz_LL_24)),
main="(D) Estimating infection time: individual Braz_1", 
xlab="time (months)", ylab="probability density",
xaxt='n', yaxt='n', bty='n',
#xaxs='i',
yaxs='i',
cex.main=1.25, cex.lab=1.25, cex.axis=1.25)



points(x=t_LL_seq[which.max(Braz_LL_24)], y=max(Braz_LL_24), pch=15, cex=2, col="black")

cut_95 = exp( log(max(Braz_LL_24)) - 1.92 )

low_95 = t_LL_seq[which.min(abs(Braz_LL_24[1:which.max(Braz_LL_24)] - cut_95))]
high_95 = t_LL_seq[which.min(abs(Braz_LL_24[which.max(Braz_LL_24):length(Braz_LL_24)] - cut_95))]


points(x=rep( low_95, 2), y=c(0,1000), type='l', lty="dashed", col="black" )

points(x=rep( high_95, 2), y=c(0,1000), type='l', lty="dashed", col="black" )



axis(1, at=c(0,90,180,270,360,450,540), label=c("0","3","6","9","12","15","18"), cex.axis=1.5*axis.size)
axis(2, at=c(0,1.1*max(Braz_LL_24)), label=c("", ""), cex.axis=axis.size)







#########################################
##                                     ##
##  PANEL 5                            ## 
##  Heatmap of antigen inclusion       ##
##  probability                        ##
##                                     ##
#########################################


par(mar=c(4,3,2,1))
par(mgp=c(1.7,0.75,0))

axis.size = 1
lab.size = 1


rain_cols <- heat.colors(N_rep_sann)


rain_cols <- c("white", rev(rain_cols))


plot(x=100, y=100,
xlim=c(0,N_ant), ylim=c(0,N_panel_max),
xlab="", ylab="number of antigens in combination",
main="(E) Probability of antigen inclusion",
xaxt='n', yaxt='n', xaxs='i', yaxs='i',
cex.main=1.5, cex.lab=1.5, cex.axis=1.5)


for(i in 1:N_ant)
{
	for(j in 1:N_panel_max)
	{
		polygon(x=c(i-1,i,i,i-1), y=c(j-1,j-1,j,j),
                    border=NA, col=rain_cols[ant_prob[i,j]+1])
	}
}

axis(2, at=seq(from=0, to=N_panel_max, by=5), label=seq(from=0, to=N_panel_max, by=5), cex.axis=axis.size, las=2)

axis(1, at=(1:N_ant) - 0.5, label=PVX_names, col.axis="black", las=2, cex.axis=0.5)

for(k in 1:5)
{
	axis(1, at=((1:N_ant) - 0.5)[bb_5_index[k]], label=PVX_names[bb_5_index[k]], col.axis=ant_5_cols[k], las=2, cex.axis=0.5)
}



#############
## LEGEND

par(mar=c(3,1,2,3))
par(mgp=c(1.5,0.75,0))

plot(x=100, y=100,
xlim=c(0,1), ylim=c(0,100),
xlab="", ylab="",
main="",
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n')

for(i in 1:length(rain_cols))
{
	polygon(x=c(0,1,1,0), y=c(i-1,i-1,i,i),
              border=NA, col=rain_cols[i])	
}


axis(4, at=100*c(0,0.25,0.5,0.75,1), label=c("0%", "25%", "50%", "75%", "100%"), las=2, cex.axis=axis.size)


dev.off()











