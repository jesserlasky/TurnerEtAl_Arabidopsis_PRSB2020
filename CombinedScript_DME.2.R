

fieldD <- read.csv('CombinedFieldData_24Sep2020.csv', as.is = T)



#get mean per plant biomass in monos
biomC <- tapply(fieldD$FH_Wt[fieldD$plotType == 'mono' & fieldD$trt == 'C' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged']/fieldD$MaxPlantNum[fieldD$plotType == 'mono' & fieldD$trt == 'C'& fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], fieldD$lines1[fieldD$plotType == 'mono' & fieldD$trt == 'C'& fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)

biomS <- tapply(fieldD$FH_Wt[fieldD$plotType == 'mono' & fieldD$trt == 'S'& fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged']/fieldD$MaxPlantNum[fieldD$plotType == 'mono' & fieldD$trt == 'S' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], fieldD$lines1[fieldD$plotType == 'mono' & fieldD$trt == 'S'& fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)

#get mean SLA in monos
slaC <- tapply(fieldD$SLA[fieldD$plotType == 'mono' & fieldD$trt == 'C'], fieldD$lines1[fieldD$plotType == 'mono' & fieldD$trt == 'C'], mean, na.rm = T)

slaS <- tapply(fieldD$SLA[fieldD$plotType == 'mono' & fieldD$trt == 'S'], fieldD$lines1[fieldD$plotType == 'mono' & fieldD$trt == 'S'], mean, na.rm = T)



#get 90% and max monoculture for comp species, also figure what type of prediction plot it is
fieldD$biom90 <- NA
fieldD$biomMax <- NA
fieldD$biom50 <- NA
fieldD$slaPlot <- NA




for(i in 1:nrow(fieldD)){
	
	if(fieldD$trt[i] == 'C'){
		 fieldD$biom50[i] <- quantile(biomC[unlist(fieldD[i,paste0('lines', 1:20)])], 0.5, na.rm = T)
		 fieldD$biom90[i] <- quantile(biomC[unlist(fieldD[i,paste0('lines', 1:20)])], 0.9, na.rm = T)
		 fieldD$biomMax[i] <- max(biomC[unlist(fieldD[i,paste0('lines', 1:20)])],  na.rm = T)
		 
		 fieldD$slaPlot[i] <- mean(slaC[unlist(fieldD[i,paste0('lines', 1:20)])],  na.rm = T)


		}
	
		if(fieldD$trt[i] == 'S'){
		 fieldD$biom90[i] <- quantile(biomS[unlist(fieldD[i,paste0('lines', 1:20)])], 0.9, na.rm = T)
		 fieldD$biom50[i] <- quantile(biomS[unlist(fieldD[i,paste0('lines', 1:20)])], 0.5, na.rm = T)
		 fieldD$biomMax[i] <- max(biomS[unlist(fieldD[i,paste0('lines', 1:20)])],  na.rm = T)
		 
		 fieldD$slaPlot[i] <- mean(slaS[unlist(fieldD[i,paste0('lines', 1:20)])],  na.rm = T)


		}
		
		}







#get expected biomass for each plot
fieldD$ExpBiomPerPlant <- NA

for(i in 1:nrow(fieldD)){
	if(fieldD$plotType[i] != 'mono'){
		if(fieldD$trt[i] == 'C') fieldD$ExpBiomPerPlant[i] <- mean(biomC[unique(unlist(fieldD[i,paste0('lines', 1:20)]))],na.rm = T)
		if(fieldD$trt[i] == 'S') fieldD$ExpBiomPerPlant[i] <- mean(biomS[unique(unlist(fieldD[i,paste0('lines', 1:20)]))],na.rm = T)

	}
}

fieldD$DivEffect <- fieldD$FH_Wt - (fieldD$ExpBiomPerPlant * fieldD$MaxPlantNum)



cor.test(fieldD$FH_Wt[fieldD$MaxPlantNum > 10 & fieldD$trt == 'S' & fieldD$Damaged_pot != 'damaged'], (fieldD$ExpBiomPerPlant[fieldD$MaxPlantNum > 10 & fieldD$trt == 'S' & fieldD$Damaged_pot != 'damaged'] * fieldD$MaxPlantNum[fieldD$MaxPlantNum > 10 & fieldD$trt == 'S' & fieldD$Damaged_pot != 'damaged']))


cor.test(fieldD$FH_Wt[fieldD$MaxPlantNum > 10 & fieldD$trt == 'C' & fieldD$Damaged_pot != 'damaged'], (fieldD$ExpBiomPerPlant[fieldD$MaxPlantNum > 10 & fieldD$trt == 'C' & fieldD$Damaged_pot != 'damaged'] * fieldD$MaxPlantNum[fieldD$MaxPlantNum > 10 & fieldD$trt == 'C' & fieldD$Damaged_pot != 'damaged']))




#generate null expected compositions of diversity plots
#get survival rate of each line
isurvC <- tapply(fieldD$MaxPlantNum[fieldD$plotType == 'mono' & fieldD$trt == 'C'], fieldD$lines1[fieldD$plotType == 'mono' & fieldD$trt == 'C'], mean, na.rm = T)

isurvS <- tapply(fieldD$MaxPlantNum[fieldD$plotType == 'mono' & fieldD$trt == 'S'], fieldD$lines1[fieldD$plotType == 'mono' & fieldD$trt == 'S'], mean, na.rm = T)



#get expected biomass for each plot, accounting for different survival rates
RandExpBiomPerPlant <- matrix(NA, nrow = nrow(fieldD), ncol = 100)


for(i in 1:nrow(fieldD)){
	for(j in 1:100){

		if(fieldD$plotType[i] != 'mono'){
			
			if(fieldD$trt[i] == 'C'){
			expg <- 	names(isurvC)[names(isurvC) %in% fieldD[i,paste0('lines', 1:fieldD$divLevel[i])]]
			expg <- rep(expg, 20/fieldD$divLevel[i])
			
			#make random comm.
			if(fieldD$MaxPlantNum[i] < 20) {
				newcomm <- sample(expg, fieldD$MaxPlantNum[i], prob = isurvC[expg])
				}else{
					newcomm <- expg
				}
			RandExpBiomPerPlant[i,j] <- mean(biomC[newcomm],na.rm = T)
				
				 }
				 
			if(fieldD$trt[i] == 'S'){
			expg <- 	names(isurvS)[names(isurvS) %in% fieldD[i,paste0('lines', 1:fieldD$divLevel[i])]]
			expg <- rep(expg, 20/fieldD$divLevel[i])
			
			#make random comm.
			if(fieldD$MaxPlantNum[i] < 20) {
				newcomm <- sample(expg, fieldD$MaxPlantNum[i], prob = isurvS[expg])
				}else{
					newcomm <- expg
				}

			RandExpBiomPerPlant[i,j] <- mean(biomS[newcomm],na.rm = T)
				 }

			}
		}
	}



##calculate a new div effect (DME), based on resampled potential stand composition

fieldD$DivEffectR <- fieldD$FH_Wt - (rowMeans(RandExpBiomPerPlant) * fieldD$MaxPlantNum)




wilcox.test(fieldD$DivEffect[fieldD$trt == 'S' &  fieldD$MaxPlantNum  >  10 & fieldD$Damaged_pot != 'damaged'])
wilcox.test(fieldD$DivEffect[fieldD$trt == 'C' &  fieldD$MaxPlantNum  >  10 & fieldD$Damaged_pot != 'damaged'])

median(fieldD$DivEffect[fieldD$trt == 'S' &  fieldD$MaxPlantNum  >  10 & fieldD$Damaged_pot != 'damaged'], na.rm = T)
median(fieldD$DivEffect[fieldD$trt == 'C' &  fieldD$MaxPlantNum  >  10 & fieldD$Damaged_pot != 'damaged'], na.rm = T)


wilcox.test(fieldD$DivEffectR[fieldD$trt == 'S' &  fieldD$MaxPlantNum  >  10  & fieldD$Damaged_pot != 'damaged'])
wilcox.test(fieldD$DivEffectR[fieldD$trt == 'C' &  fieldD$MaxPlantNum  >  10  & fieldD$Damaged_pot != 'damaged'])

median(fieldD$DivEffectR[fieldD$trt == 'S' &  fieldD$MaxPlantNum  >  10], na.rm = T)
median(fieldD$DivEffectR[fieldD$trt == 'C' &  fieldD$MaxPlantNum  >  10], na.rm = T)


#####



#comparing different prediction stands

fieldD$typePlot2 <- fieldD$typePlot
fieldD$typePlot2[fieldD$plotType == 'div'] <- 'div'


a1 <- aov(DivEffect ~ typePlot2 + trt, data = fieldD[fieldD$MaxPlantNum > 10   & fieldD$Damaged_pot != 'damaged',])
summary(a1)
posthoc <- TukeyHSD(x=a1, 'typePlot2', conf.level=0.95)
posthoc

#
a1 <- aov(DivEffect ~ typePlot + trt, data = fieldD[fieldD$MaxPlantNum > 10   & fieldD$Damaged_pot != 'damaged',])
summary(a1)
posthoc <- TukeyHSD(x=a1, 'typePlot2', conf.level=0.95)
posthoc
#########


library(RColorBrewer)
pdf('~/Dropbox/jesse/Arabidopsis/BEF/Results/DivEffect_v_biomFT90_predComp_24Sep2020.pdf', pointsize = 12, height = 8)

par(mfrow = c(2,2))
par(mar = c(6,4,3,1))

plot(fieldD$biom90[fieldD$MaxPlantNum > 10    & fieldD$Damaged_pot != 'damaged'], fieldD$DivEffect[fieldD$MaxPlantNum > 10    & fieldD$Damaged_pot != 'damaged'], col = brewer.pal(9, 'RdBu')[c(1,9)][1 + (fieldD$trt[fieldD$MaxPlantNum > 10    & fieldD$Damaged_pot != 'damaged'] == 'S')], bty = 'n', xlab = '90th percentile of genotype biomass (g)', ylab = 'DME (g)', xlim = c(0, 0.032), pch = c(16,17)[1 + (fieldD$trt[fieldD$MaxPlantNum > 10    & fieldD$Damaged_pot != 'damaged'] == 'S')])
abline(h = 0, lty = 2)



plot(fieldD$FT1001_mean[fieldD$MaxPlantNum > 10    & fieldD$Damaged_pot != 'damaged'], fieldD$DivEffect[fieldD$MaxPlantNum > 10   & fieldD$Damaged_pot != 'damaged'], col = brewer.pal(9, 'RdBu')[c(1,9)][1 + (fieldD$trt[fieldD$MaxPlantNum > 10   & fieldD$Damaged_pot != 'damaged'] == 'S')], bty = 'n', xlab = '90th percentile of Flowering Time (d)', ylab = 'DME (g)', pch = c(16,17)[1 + (fieldD$trt[fieldD$MaxPlantNum > 10   & fieldD$Damaged_pot != 'damaged'] == 'S')])#, xlim = c(0, 0.032))

legend('topright', col = brewer.pal(9, 'RdBu')[c(1,9)], pch = 16:17, legend = c('High', 'Low'), cex = 0.9, title = 'Resource treatment')

par(xpd = T)

par(xpd = F)


abline(h = 0, lty = 2)

fieldD$typePlot <- factor(fieldD$typePlot, levels = c('EF', 'LF', 'Gdiv', 'SimClimGdiv', 'SimClimFdiv'))

fieldD$typePlot2a <- as.character(fieldD$typePlot2)
fieldD$typePlot2a[fieldD$typePlot2a == 'div'] <- 'Random'
fieldD$typePlot2a <- factor(fieldD$typePlot2a, levels = c('Random', 'EF', 'LF', 'Gdiv', 'SimClimGdiv', 'SimClimFdiv'))



boxplot(DivEffect ~ typePlot2a, data = fieldD[fieldD$trt == 'C' & fieldD$MaxPlantNum > 10    & fieldD$Damaged_pot != 'damaged',], range = 0, ylab = 'DME (g)', xlab = '', las = 2, names = c('Random','EF', 'LF', 'Gdiv', 'CliGdiv', 'CliFdiv'))

mtext('High resource', line = 0.5)
mtext('Stand assembly criterion', side = 1, line = 4.5)
abline(h = 0, lty = 2)



boxplot(DivEffect ~ typePlot2a, data = fieldD[fieldD$trt == 'S' & fieldD$MaxPlantNum > 10    & fieldD$Damaged_pot != 'damaged',], range = 0, ylab = 'DME (g)', xlab = '', las = 2, names = c('Random','EF', 'LF', 'Gdiv', 'CliGdiv', 'CliFdiv'))
abline(h = 0, lty = 2)

mtext('Stand assembly criterion', side = 1, line = 4.5)
mtext('Low resource', line = 0.5)


par(xpd = T)

par(xpd = F)



dev.off()





###regression of stand characteristics versus DME (DivEffect)
#tables to save results
tmp <- c()
tmp2 <- c() #w FT



tmod <- lm(DivEffect  ~ trt + trt: FT1001_q90 , data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )

tmp <- rbind(tmp, cbind(summary(tmod)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod)[[9]]))


tmod <- lm(DivEffect  ~ trt + trt: FT1001_mean , data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp <- rbind(tmp, cbind(summary(tmod)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod)[[9]]))


tmod <- lm(DivEffect  ~ trt + trt: biom90 , data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp <- rbind(tmp, cbind(summary(tmod)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod)[[9]]))


tmod <- lm(DivEffect  ~ trt + trt: biom50 , data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp <- rbind(tmp, cbind(summary(tmod)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod)[[9]]))

tmod2 <- lm(DivEffect  ~ trt + trt: biom50 + trt: FT1001_mean, data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp2 <- rbind(tmp2, cbind(summary(tmod2)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod2)[[9]]))


tmod <- lm(DivEffect  ~ trt + trt: slaPlot , data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp <- rbind(tmp, cbind(summary(tmod)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod)[[9]]))

tmod2 <- lm(DivEffect  ~ trt + trt: slaPlot + trt: FT1001_mean, data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp2 <- rbind(tmp2, cbind(summary(tmod2)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod2)[[9]]))





tmod <- lm(DivEffect  ~ trt + trt: meanKin , data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp <- rbind(tmp, cbind(summary(tmod)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod)[[9]]))

tmod2 <- lm(DivEffect  ~ trt + trt: meanKin + trt: FT1001_mean, data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp2 <- rbind(tmp2, cbind(summary(tmod2)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod2)[[9]]))






tmod <- lm(DivEffect  ~ trt + trt: PC5dist_mean , data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp <- rbind(tmp, cbind(summary(tmod)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod)[[9]]))

tmod2 <- lm(DivEffect  ~ trt + trt: PC5dist_mean + trt: FT1001_mean, data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp2 <- rbind(tmp2, cbind(summary(tmod2)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod2)[[9]]))



tmod <- lm(DivEffect  ~ trt + trt: divLevel , data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp <- rbind(tmp, cbind(summary(tmod)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod)[[9]]))

tmod2 <- lm(DivEffect  ~ trt + trt: divLevel + trt: FT1001_mean, data = fieldD[fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged',] )
tmp2 <- rbind(tmp2, cbind(summary(tmod2)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod2)[[9]]))





tmod <- lm(DivEffect  ~ trt + trt: treeDiv , data = fieldD[fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged',] )
tmp <- rbind(tmp, cbind(summary(tmod)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod)[[9]]))

tmod2 <- lm(DivEffect  ~ trt + trt: treeDiv + trt: FT1001_mean, data = fieldD[fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged',] )
tmp2 <- rbind(tmp2, cbind(summary(tmod2)[[4]][3:4,c('Estimate', 'Pr(>|t|)')], r2adj = summary(tmod2)[[9]]))

tmp <- data.frame(tmp)
tmp$trt <- NA
tmp$trt[grep('C', rownames(tmp))] <- 'High'
tmp$trt[grep('S', rownames(tmp))] <- 'Low'

#table S6
#write.csv(tmp, 'lmCovarsDivEffectTable_13Apr2020.csv')

tmp2 <- data.frame(tmp2)
tmp2$trt <- NA
tmp2$trt[grep('C', rownames(tmp2))] <- 'High'
tmp2$trt[grep('S', rownames(tmp2))] <- 'Low'

#write.csv(tmp2, 'lmCovars_wFT_DivEffectTable_13Apr2020.csv')





