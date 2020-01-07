

###########################
###########################
##                       ##
##  #   ## ##### ######  ##
##  ##  ## ##      ##    ##
##  ### ## ####    ##    ##
##  ## ### ##      ##    ##
##  ##  ## #####   ##    ##
##                       ##
###########################
###########################


library(igraph)

load("C:\\U\\GHIT\\NatMed_Paper_OrigBrazData\\Figure4_analysis\\LDA_search.RData")


top_4ant = LDA_4ant[order(LDA_4ant[,5], decreasing=TRUE),]

VV_4 = rep(0, N_ant)
MM_4 = matrix(0, nrow=N_ant, ncol=N_ant)

for(i in 1:nrow(top_4ant))
{
	VV_4[which(ant_names_short == top_4ant[i,1])] = VV_4[which(ant_names_short == top_4ant[i,1])] + i	
	VV_4[which(ant_names_short == top_4ant[i,2])] = VV_4[which(ant_names_short == top_4ant[i,2])] + i
	VV_4[which(ant_names_short == top_4ant[i,3])] = VV_4[which(ant_names_short == top_4ant[i,3])] + i
	VV_4[which(ant_names_short == top_4ant[i,4])] = VV_4[which(ant_names_short == top_4ant[i,4])] + i

	MM_4[which(ant_names_short == top_4ant[i,1]),which(ant_names_short == top_4ant[i,2])] = MM_4[which(ant_names_short == top_4ant[i,1]),which(ant_names_short == top_4ant[i,2])] + i
	MM_4[which(ant_names_short == top_4ant[i,1]),which(ant_names_short == top_4ant[i,3])] = MM_4[which(ant_names_short == top_4ant[i,1]),which(ant_names_short == top_4ant[i,3])] + i
	MM_4[which(ant_names_short == top_4ant[i,2]),which(ant_names_short == top_4ant[i,3])] = MM_4[which(ant_names_short == top_4ant[i,2]),which(ant_names_short == top_4ant[i,3])] + i
	MM_4[which(ant_names_short == top_4ant[i,1]),which(ant_names_short == top_4ant[i,4])] = MM_4[which(ant_names_short == top_4ant[i,1]),which(ant_names_short == top_4ant[i,4])] + i
	MM_4[which(ant_names_short == top_4ant[i,2]),which(ant_names_short == top_4ant[i,4])] = MM_4[which(ant_names_short == top_4ant[i,2]),which(ant_names_short == top_4ant[i,4])] + i
	MM_4[which(ant_names_short == top_4ant[i,3]),which(ant_names_short == top_4ant[i,4])] = MM_4[which(ant_names_short == top_4ant[i,3]),which(ant_names_short == top_4ant[i,4])] + i
}

MM_4 = MM_4 + t(MM_4)


MM_4 = exp(-MM_4/quantile(MM_4, prob=0.3))

diag(MM_4) = 0

MM_4 = MM_4/mean(MM_4)




VV_4 = exp(-VV_4/quantile(VV_4, prob=0.5))

VV_4 = VV_4/mean(VV_4)




NET_4 = graph.adjacency(MM_4, mode="undirected", weighted=TRUE, diag=FALSE)

summary(NET_4)




color_scale = (E(NET_4)$weight - min(E(NET_4)$weight))/(max(E(NET_4)$weight) - min(E(NET_4)$weight))





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
## Thailand

thailand_cat <- rep("thai_never", nrow(thailand_data) )

thailand_cat[which( thailand_data[,10] == 0 )] <- "thai_current"

thailand_cat[intersect( which(thailand_data[,10]>0), which(thailand_data[,10]<=9*30) )] <- "thai_recent"

thailand_cat[which(thailand_data[,10]>9*30)] <- "thai_old"


###################################
## Brazil

brazil_cat <- rep("braz_never", nrow(brazil_data) )

brazil_cat[which( brazil_data[,10] == 0 )] <- "braz_current"

brazil_cat[intersect( which(brazil_data[,10]>0), which(brazil_data[,10]<=9*30) )] <- "braz_recent"

brazil_cat[which(brazil_data[,10]>9*30)] <- "braz_old"


###################################
## Solomon

solomon_cat <- rep("sol_never", nrow(solomon_data) )

solomon_cat[which( solomon_data[,10] == 0 )] <- "sol_current"

solomon_cat[intersect( which(solomon_data[,10]>0), which(solomon_data[,10]<=9*30) )] <- "sol_recent"

solomon_cat[which(solomon_data[,10]>9*30)] <- "sol_old"



###################################
## Controls

control_cat <- rep("VBDR", nrow(control_data) )

for(i in 1:length(control_cat))
{
	if( substr( control_data[i,1], 1, 3 ) == "TRC" )
	{
		control_cat[i] <- "TRC"
	}
	if( substr( control_data[i,1], 1, 3 ) == "ARC" )
	{
		control_cat[i] <- "ARC"
	}
	if( substr( control_data[i,1], 1, 3 ) == "BRC" )
	{
		control_cat[i] <- "BRC"
	}
}

####################################
## Put together Thai and control data

AB <- rbind( thailand_data[,17:81], brazil_data[,17:81], solomon_data[,17:81], control_data[,17:81] )

AB <- log(AB)



ant_drop <- c()

for(j in 1:ncol(AB))
{
	if( length(which(is.na(AB[,j]))) > 250 )
	{
		ant_drop = c(ant_drop, j)
	}
}

AB = AB[,-ant_drop]

N_ant <- ncol(AB)
ant_names <- colnames(AB)



############################################
## Create shortened antibody names

ant_names_short = c("W01", "W02", "W03", "W04", "W05", "W06", "W07", "W08", "W09", "W10",
                    "W11", "W12", "W13", "W14", "W15", "W16", "W17", "W18", "W19", "W20",
                    "W21", "W22", "W23", "W24", "W25", "W26", "W27", "W28", "W29", "W30",
                    "W31", "W32", "W33", "W34", "W35", "W36", "W37", "W38", "W39", "W40",
                    "W41", "W42", "W43", "W44", "W45", "W46", "W47", "W48", "W49", "W50",
                    "W51", "W52", "W53", "W54", "W55", "W56", "W57", "W58", "W59", "W60") 


inf_cat <- c( thailand_cat, brazil_cat, solomon_cat, control_cat )

N_part <- length(inf_cat)





###################################
## Binary category

bin_cat = rep("old", N_part)
bin_cat[which(inf_cat=="thai_current")] = "new"
bin_cat[which(inf_cat=="thai_recent")]  = "new"
bin_cat[which(inf_cat=="braz_current")] = "new"
bin_cat[which(inf_cat=="braz_recent")]  = "new"
bin_cat[which(inf_cat=="sol_current")] = "new"
bin_cat[which(inf_cat=="sol_recent")]  = "new"



##############################
##############################
##                          ## 
##  #####  ##   ## ##   ##  ##
##  ##  ## ##   ##  ## ##   ##
##  #####   ## ##    ###    ##
##  ##       ###    ## ##   ##
##  ##        #    ##   ##  ##
##                          ## 
##############################
##############################

PVX_read = read.csv("C:\\U\\GHIT\\NatMed_Paper\\Data\\protein_info\\Table2_protein.csv")

PVX_names <- rep(NA, N_ant)

for(i in 1:N_ant)
{
	PVX_names[i] <- as.vector(PVX_read[which(PVX_read[,1] == ant_names[i]),3])
}





################################################################################
################################################################################
##                                                                            ##
##   ####   ####  #####  #####  ##### ##     ####  ###### ####  ####  #   ##  ##
##  ##  ## ##  ## ##  ## ##  ## ##    ##    ##  ##   ##    ##  ##  ## ##  ##  ##
##  ##     ##  ## #####  #####  ####  ##    ######   ##    ##  ##  ## ### ##  ##
##  ##  ## ##  ## ## ##  ## ##  ##    ##    ##  ##   ##    ##  ##  ## ## ###  ##
##   ####   ####  ##  ## ##  ## ##### ##### ##  ##   ##   ####  ####  ##  ##  ##
##                                                                            ##
################################################################################
################################################################################

library(fields)

AB_cor <- matrix(NA, nrow=N_ant, ncol=N_ant)

for(i in 1:N_ant)
{
	for(j in 1:N_ant)
	{
		index_ij = intersect( which( is.na(AB[,i])==FALSE ), which( is.na(AB[,j])==FALSE ) )

		AB_cor[i,j] = cor( AB[index_ij,i], AB[index_ij,j], method="spearman" )
	}
}







####################################
####################################
##                                ##
##  ONE AT A TIME                 ##
##                                ##
####################################
####################################

N_SS <- 1000

SS_cut <- seq(from=-11, to=-3.5, length=N_SS)

sens_mat <- matrix(NA, nrow=N_ant, ncol=N_SS)
spec_mat <- matrix(NA, nrow=N_ant, ncol=N_SS)


for(i in 1:N_ant)
{
	AB_i  = AB[which(is.na(AB[,i])==FALSE),i]

	bin_cat_i = bin_cat[which(is.na(AB[,i])==FALSE)]


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

##sens_mat = cbind( rep(1,N_SS), sens_mat, rep(0,N_SS) )
##spec_mat = cbind( rep(0,N_SS), spec_mat, rep(1,N_SS) )



######################################
## Calculate Area Under Curve (AUC)


AUC_one_ant <- rep(NA, N_ant)

for(i in 1:N_ant)
{
	AUC_one_ant[i] <- sum( (sens_mat[i,1:(N_SS-1)] - sens_mat[i,2:N_SS])*
                             0.5*(spec_mat[i,1:(N_SS-1)] + spec_mat[i,2:N_SS]) )
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


top_8_final = order( VV_4, decreasing=TRUE )[c(1,3,4,7,9,5,8,14)]





top_8_names = ant_names[top_8_final] ##  ant_names_PVX[top_8_final]

top_8_PVX = PVX_names[top_8_final] ##  ant_names_PVX[top_8_final]


top_8_cols = rainbow(8)

alpha_seq = c("A", "B", "C", "D", "E", "F", "G", "H")

region_cols = c("magenta", "green3", "yellow2", "dodgerblue")



top_8_names_long <- rep(NA, 8)

for(i in 1:8)
{
	top_8_names_long[i] <- paste( ant_names[top_8_final[i]], ": ", PVX_names[top_8_final[i]], sep="" )
}



tiff(file="Figure4_Validation_final8_big_labels.tif", width=40, height=30, units="cm", res=500)

lay.mat <- rbind( c( 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 ),
                  c( 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8 ), 
                  c( 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 ),
                  c(10,10,10,10,11,11,11,11,13,13,13,13 ),
                  c(10,10,10,10,12,12,12,12,13,13,13,13 )  )
layout(lay.mat, heights=c(15, 15, 3, (4/3)*13, (4/3)*2))
layout.show(13)

par(mar=c(3,5,2.5,1))
par(mgp=c(3, 0.65, 0))

point.size = 0.75
lab.size   = 2
axis.size  = 2
main.size  = 2

#####################################
#####################################
##                                 ## 
##  PANELS 1-8                     ##
##  Boxplots for Ab distributions  ##
##                                 ##
#####################################
#####################################

for(n in 1:8)
{

	boxplot( AB[which(inf_cat=="thai_current"),top_8_final[n]],
      	   AB[which(inf_cat=="thai_recent"),top_8_final[n]],
      	   AB[which(inf_cat=="thai_old"),top_8_final[n]],
      	   AB[which(inf_cat=="thai_never"),top_8_final[n]],
		   AB[which(inf_cat=="braz_current"),top_8_final[n]],
      	   AB[which(inf_cat=="braz_recent"),top_8_final[n]],
      	   AB[which(inf_cat=="braz_old"),top_8_final[n]],
      	   AB[which(inf_cat=="braz_never"),top_8_final[n]],
		   AB[which(inf_cat=="sol_current"),top_8_final[n]],
      	   AB[which(inf_cat=="sol_recent"),top_8_final[n]],
      	   AB[which(inf_cat=="sol_old"),top_8_final[n]],
      	   AB[which(inf_cat=="sol_never"),top_8_final[n]],
      	   AB[which(inf_cat=="TRC"),top_8_final[n]],
      	   AB[which(inf_cat=="BRC"),top_8_final[n]],
      	   AB[which(inf_cat=="ARC"),top_8_final[n]],
      	   AB[which(inf_cat=="VBDR"),top_8_final[n]],
	pch=19, yaxt='n', xaxt='n',
	ylim=log(c(1e-5, 0.03)),
	col=c("darkgrey", "red", "orange" ,"green", 
	      "darkgrey", "red", "orange" ,"green", 
	      "darkgrey", "red", "orange" ,"green", 
	      "royalblue", "cornflowerblue", "dodgerblue", "cyan"),
	ylab="relative antibody unit", 
	main=paste("(", alpha_seq[n], ") anti-", top_8_PVX[n], " antibodies", sep=""),
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


	points(x=c(-1e6, 1e6), y=rep(log(1.95e-5),2), type='l', lty="dashed")
	points(x=c(-1e6, 1e6), y=rep(log(0.02),2), type='l', lty="dashed")
	
	axis(1, at=c(2.5, 7, 10.5, 14.5), label=c("Thailand", "Brazil", "Solomons", "Controls"), cex.axis=1.3 )
	axis(2, at=log(c(0.00001, 0.0001, 0.001, 0.01)), label=c("0.00001", "0.0001", "0.001", "0.01"), las=2, cex.axis=0.9 )
}

#####################################
#####################################
##                                 ## 
##  PANELS 9                       ##
##  Legend                         ##
##                                 ##
#####################################
#####################################

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("infected",              "negative controls: Thai RC", 
                  "infected 1-9 months",   "negative controls: Brazil RC", 
                  "infected 9-12 months",  "negative controls: Aus RC",
                  "no detected infection", "negative controls: VBDR"), 
       fill = c("darkgrey", "royalblue", 
                "red",      "cornflowerblue",  
                "orange",   "dodgerblue", 
                "green",    "cyan"), 
       border = c("darkgrey", "royalblue", 
                "red",        "cornflowerblue",  
                "orange",     "dodgerblue", 
                "green",      "cyan"), 
       ncol=4, cex=2.1, bty="n" )





#########################################
#########################################
##                                     ## 
##  PANEL 10                           ##
##  Binary classification (all sites)  ##
##                                     ##
#########################################
#########################################


line_seq <- c(0.2, 0.4, 0.6, 0.8)

par(mar=c(4.5,5,2.5,1.5))
par(mgp=c(3.25, 0.6,0))


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
xaxs='i', yaxs='i', xaxt='n', yaxt='n',
xlab="1 - specificity", ylab="sensitivity",
main="(I) Classification of recent infections",
cex.lab=1.2*lab.size, cex.axis=axis.size, cex.main=1.1*main.size)


for(i in 1:4)
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}



for(i in 1:N_ant)
{
	points( x=1-spec_mat[i,], y=sens_mat[i,], 
     	        type='S', lwd=1, col="grey" )
}

for(i in 1:N_ant)
{	
	if( i %in% top_8_final )
	{
		points( x=1-spec_mat[i,], y=sens_mat[i,], 
      	        type='s', lwd=2, col=top_8_cols[which(top_8_final==i)] )
	}

}


##points(x=0.2, y=0.8, pch=17, cex=2, col="red")
##points(x=0.02, y=0.5, pch=17, cex=2, col="green")
##points(x=0.5, y=0.98, pch=17, cex=2, col="blue")


legend(x="bottomright", 
cex=1.1, 
bg="white", box.col="white",
fill = top_8_cols,
border = top_8_cols, 
legend = top_8_names_long )


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.5 ) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.5, las=2 ) 



#########################################
#########################################
##                                     ## 
##  PANEL 11                           ##
##  Correlation                        ##
##                                     ##
#########################################
#########################################



par(mar=c(2,3,2.5,1))
par(mgp=c(1.5,0.75,0))


##N_cor_steps <- 100
##
##cor_cols <- rev(heat.colors(N_cor_steps))


cor_cols = c("springgreen4", "springgreen", "palegreen", "yellowgreen", "yellow", 
             "gold", "orange", "orangered", "firebrick1", "red3" )

N_cor_steps = length(cor_cols)

plot(x=100, y=100,
xlim=c(0,N_ant), ylim=c(0,N_ant),
xlab="", ylab="",
main="(J) Correlation between antibody titres",
xaxt='n', yaxt='n', xaxs='i', yaxs='i',
cex.lab=lab.size, cex.axis=axis.size, cex.main=1.1*main.size)


for(i in 1:N_ant)
{
	for(j in 1:N_ant)
	{
		polygon(x=c(i-1,i,i,i-1), y=c(j-1,j-1,j,j),
                    border=NA, col=cor_cols[N_cor_steps*AB_cor[i,j]])
	}
}


axis(1, at=seq(from=0.5, by=1, length=N_ant), label=ant_names_short, las=2, cex.axis=0.5)
axis(2, at=seq(from=0.5, by=1, length=N_ant), label=ant_names_short, las=2, cex.axis=0.5)



#############
## LEGEND

par(mar=c(2,1,1,1))


plot(x=100, y=100,
xlim=c(0,100), ylim=c(0,1),
xlab="", ylab="",
main="",
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n')

for(i in 1:length(cor_cols))
{
	polygon(y=c(0,1,1,0), x=c(i-1,i-1,i,i)*100/N_cor_steps,
              border=NA, col=cor_cols[i])	
}


axis(1, at=100*c(0,0.25,0.5,0.75,1), label=c("0%", "25%", "50%", "75%", "100%"), cex.axis=1.25)





#####################################
#####################################
##                                 ## 
##  PANELS 12                      ##
##  Multi-variate data view        ##
##                                 ##
#####################################
#####################################


par(mar=c(4.5,5,2.5,1))
par(mgp=c(3.25, 0.65, 0))


#####################################
## Colouring by infection status

plot(x=1e10, y=1e10, 
pch=19, yaxt='n', xaxt='n',
xlim=log(c(1e-5, 0.03)), ylim=log(c(1e-5, 0.03)),
xlab=paste("anti-", PVX_names[top_8_final[1]], " antibody titre", sep=""),
ylab=paste("anti-", PVX_names[top_8_final[2]], " antibody titre", sep=""),
main="(K) Distribution of antibody titres",
cex.lab=lab.size, cex.axis=axis.size, cex.main=1.1*main.size)

points( x=AB[ which(inf_cat%in%c("thai_never", "thai_old", "thai_recent", "thai_current")), top_8_final[1]],
        y=AB[ which(inf_cat%in%c("thai_never", "thai_old", "thai_recent", "thai_current")), top_8_final[2]],
        pch=19, cex=point.size, col=region_cols[1] )

points( x=AB[ which(inf_cat%in%c("braz_never", "braz_old", "braz_recent", "braz_current")), top_8_final[1]],
        y=AB[ which(inf_cat%in%c("braz_never", "braz_old", "braz_recent", "braz_current")), top_8_final[2]],
        pch=19, cex=point.size, col=region_cols[2] )

points( x=AB[ which(inf_cat%in%c("sol_never", "sol_old", "sol_recent", "sol_current")), top_8_final[1]],
        y=AB[ which(inf_cat%in%c("sol_never", "sol_old", "sol_recent", "sol_current")), top_8_final[2]],
        pch=19, cex=point.size, col=region_cols[3] )


points( x=AB[ which(inf_cat%in%c("TRC")), top_8_final[1]],
        y=AB[ which(inf_cat%in%c("TRC")), top_8_final[2]],
        pch=19, cex=point.size, col=region_cols[4] )

points( x=AB[ which(inf_cat%in%c("BRC")), top_8_final[1]],
        y=AB[ which(inf_cat%in%c("BRC")), top_8_final[2]],
        pch=19, cex=point.size, col=region_cols[4] )

points( x=AB[ which(inf_cat%in%c("ARC")), top_8_final[1]],
        y=AB[ which(inf_cat%in%c("ARC")), top_8_final[2]],
        pch=19, cex=point.size, col=region_cols[4] )

points( x=AB[ which(inf_cat%in%c("VBDR")), top_8_final[1]],
        y=AB[ which(inf_cat%in%c("VBDR")), top_8_final[2]],
        pch=19, cex=point.size, col=region_cols[4] )




points(x=c(-1e6, 1e6), y=rep(log(1.95e-5),2), type='l', lty="dashed")
points(x=c(-1e6, 1e6), y=rep(log(0.02),2), type='l', lty="dashed")

points(y=c(-1e6, 1e6), x=rep(log(1.95e-5),2), type='l', lty="dashed")
points(y=c(-1e6, 1e6), x=rep(log(0.02),2), type='l', lty="dashed")

legend(x="bottomright",
fill = region_cols,
border = region_cols,
bg="white", box.col="white",
legend = c("Thailand", "Brazil", "Solomons", "Controls"),
cex=1.25 )


axis(1, at=log(c(0.00001, 0.0001, 0.001, 0.01)), label=c("0.00001", "0.0001", "0.001", "0.01") )
axis(2, at=log(c(0.00001, 0.0001, 0.001, 0.01)), label=c("0.00001", "0.0001", "0.001", "0.01"), las=2, cex.axis=0.9 )




dev.off()







