
library(VCA)

fieldD <- read.csv('CombinedFieldData_24Sep2020.csv', as.is = T)

dats <- fieldD

####
dats$SLA[is.na(dats$SLA_Area)] <- NA
dats$logSLA <- log(dats$SLA)
dats$biomPplant <- dats$FH_Wt / dats$MaxPlantNum

dats$biomPplant[dats$PotID %in% fieldD$PotID[fieldD$Damaged_pot == 'damaged'] | dats$MaxPlantNum <= 10] <- NA



dats$potF2[dats$PotID %in% fieldD$PotID[fieldD$Damaged_pot == 'damaged'] | dats$MaxPlantNum <= 10] <- NA

dats$fecPplant <- dats$potF2 


##convert FT to numeric
dats$FT1 <- NA
dats$FT50 <- NA
for (i in 1:nrow(dats)){
	
	dats$FT1[i] <- which(dats[i,grep('FT\\.', colnames(dats))] > 0)[1]
	dats$FT50[i] <- which(dats[i,grep('FT\\.', colnames(dats))] > 49)[1]
	
	
}
dats$FT1d <- as.numeric(as.Date(c('5/29/18', '6/4/18', '6/8/18', '6/12/18', '6/15/18', '6/19/18', '6/21/18', '6/26/18', '6/29/18', '7/6/18', '7/10/18')[dats$FT1], format = '%m/%d/%y'))
dats$FT50d <- as.numeric(as.Date(c('5/29/18', '6/4/18', '6/8/18', '6/12/18', '6/15/18', '6/19/18', '6/21/18', '6/26/18', '6/29/18', '7/6/18', '7/10/18')[dats$FT50], format = '%m/%d/%y'))


max(dats$FT1d, na.rm = T)

dats$FT1d_InclNoFlow <- dats$FT1d
dats$FT1d_InclNoFlow[is.na(dats$FT1d_InclNoFlow)] <- max(dats$FT1d, na.rm = T) + 2



#heritability
hertable <- c()

hertable <- rbind(hertable, c('logSLA',
	summary(lm(logSLA ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'C',]))[[8]],
	summary(lm(logSLA ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'S',]))[[8]]
		))

hertable <- rbind(hertable, c('FT1d',
	summary(lm(FT1d ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'C',]))[[8]],
	summary(lm(FT1d ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'S',]))[[8]]
	))


hertable <- rbind(hertable, c('FT1d_InclNoFlow',

summary(lm(FT1d_InclNoFlow ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'C',]))[[8]],
summary(lm(FT1d_InclNoFlow ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'S',]))[[8]]
	))


hertable <- rbind(hertable, c('FT50d',

summary(lm(FT50d ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'C',]))[[8]],
summary(lm(FT50d ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'S',]))[[8]]
	))




hertable <- rbind(hertable, c('biomPplant',

	summary(lm(biomPplant ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'C',]))[[8]],
	summary(lm(biomPplant ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'S',]))[[8]]
	))


hertable <- rbind(hertable, c('fecPplant',

	summary(lm(fecPplant ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'C',]))[[8]],
	summary(lm(fecPplant ~ lines1, data = dats[dats$plotType == 'mono' & dats$trt == 'S',]))[[8]]
	))



colnames(hertable) <- c('trait', 'control', 'stress')

#write.csv(hertable, 'heritabilitiesTable.csv', row.names = F)


#get treatment effect - make sure this is just monocultures
dats$BinNoF <- paste0('B', dats$BinNo)
dats$X <- c(1:4)[match(dats$stripNo, c('A', 'B', 'C', 'D'))]
dats$Y <- floor((dats$PotID-1) / 4)
dats$Y[dats$Y >= 30] <- dats$Y[dats$Y >= 30] - 30
dats$Y[dats$Y >= 30] <- dats$Y[dats$Y >= 30] - 30
dats$Y[dats$Y >= 30] <- dats$Y[dats$Y >= 30] - 30
	
newtab <- c()
#SLA
		m1a <- anovaMM(logSLA ~ trt * (lines1) + (BinNoF) + X + Y, Data = na.omit(dats[dats$plotType == 'mono',c('logSLA', 'trt', 'lines1', 'BinNoF', 'X', 'Y')]))
		
		inf <- VCAinference(m1a, VarVC=TRUE)

		lc.mat <- VCA::getL(m1a, c("trtC-trtS"))
		
		test.fixef(m1a, lc.mat)

newtab <- rbind(newtab, c('SLA', test.fixef(m1a, lc.mat)[1,]))
	
	
#FT1d
		m1a <- anovaMM(FT1d ~ trt * (lines1) + (BinNoF) + X + Y, Data = na.omit(dats[which(dats$plotType == 'mono' ),c('FT1d', 'trt', 'lines1', 'BinNoF', 'X', 'Y')]))
		
		inf <- VCAinference(m1a, VarVC=TRUE)

		lc.mat <- VCA::getL(m1a, c("trtC-trtS"))
		
		test.fixef(m1a, lc.mat)

newtab <- rbind(newtab, c('FT1d', test.fixef(m1a, lc.mat)[1,]))
				
#FT1d_InclNoFlow
		m1a <- anovaMM(FT1d_InclNoFlow ~ trt * (lines1) + (BinNoF) + X + Y, Data = na.omit(dats[which(dats$plotType == 'mono' ),c('FT1d_InclNoFlow', 'trt', 'lines1', 'BinNoF', 'X', 'Y')]))
		
		inf <- VCAinference(m1a, VarVC=TRUE)

		lc.mat <- VCA::getL(m1a, c("trtC-trtS"))
		
		test.fixef(m1a, lc.mat)
								
newtab <- rbind(newtab, c('FT1d_InclNoFlow', test.fixef(m1a, lc.mat)[1,]))
				
								
##Biomass
		m1a <- anovaMM(biomPplant ~ trt * (lines1) + (BinNoF) + X + Y, Data = na.omit(dats[dats$plotType == 'mono',c('biomPplant', 'trt', 'lines1', 'BinNoF', 'X', 'Y')]))
		
		inf <- VCAinference(m1a, VarVC=TRUE)

		lc.mat <- VCA::getL(m1a, c("trtC-trtS"))
		
		test.fixef(m1a, lc.mat)
		
newtab <- rbind(newtab, c('Biomass', test.fixef(m1a, lc.mat)[1,]))
		
		
		
									
##fecudity average for measured plants
		m1a <- anovaMM(potF2 ~ trt * (lines1) + (BinNoF) + X + Y, Data = na.omit(dats[dats$plotType == 'mono',c('potF2', 'trt', 'lines1', 'BinNoF', 'X', 'Y')]))
		
		inf <- VCAinference(m1a, VarVC=TRUE)

		lc.mat <- VCA::getL(m1a, c("trtC-trtS"))
		
		test.fixef(m1a, lc.mat)
			

newtab <- rbind(newtab, c('avg of meas plants', test.fixef(m1a, lc.mat)[1,]))
		
		
		
		
#write.csv(newtab, 'PlasticitiesTable.csv', row.names = F)

		
		
		

