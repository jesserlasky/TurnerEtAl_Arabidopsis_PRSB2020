

fieldD <- read.csv('CombinedFieldData_24Sep2020.csv', as.is = T)


#get mean per plant biomass in monos
biomC <- tapply(fieldD$FH_Wt[fieldD$plotType == 'mono' & fieldD$trt == 'C' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged']/fieldD$MaxPlantNum[fieldD$plotType == 'mono' & fieldD$trt == 'C'& fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], fieldD$lines1[fieldD$plotType == 'mono' & fieldD$trt == 'C'& fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)

biomS <- tapply(fieldD$FH_Wt[fieldD$plotType == 'mono' & fieldD$trt == 'S'& fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged']/fieldD$MaxPlantNum[fieldD$plotType == 'mono' & fieldD$trt == 'S' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], fieldD$lines1[fieldD$plotType == 'mono' & fieldD$trt == 'S'& fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)



#get expected biomass for each plot
fieldD$ExpBiomPerPlant <- NA

for(i in 1:nrow(fieldD)){
	if(fieldD$plotType[i] != 'mono'){
		if(fieldD$trt[i] == 'C') fieldD$ExpBiomPerPlant[i] <- mean(biomC[unique(unlist(fieldD[i,paste0('lines', 1:20)]))],na.rm = T)
		if(fieldD$trt[i] == 'S') fieldD$ExpBiomPerPlant[i] <- mean(biomS[unique(unlist(fieldD[i,paste0('lines', 1:20)]))],na.rm = T)

	}
}

fieldD$DivEffect <- fieldD$FH_Wt - (fieldD$ExpBiomPerPlant * fieldD$MaxPlantNum)



combC2b <- unique(fieldD[fieldD$plotType == 'mono',c('ecotypeid', 'lines1')])
combC2b$CS_number <- combC2b$lines1





g1001 <- read.csv('1001genomesaccessionsInfo.csv')


#####



g60snps <- read.csv('60gSNPs.csv', header = T) #1001 Genomes binary SNPs for the 60 ecotype ids in fieldD
#https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX.tar.gz
chrStarts <- c(1, 1004590, 1700066, 2510573, 3181800)

g60snps$chr <- c(rep(1, chrStarts[2] - 1), rep(2, chrStarts[3] - chrStarts[2]), rep(3, chrStarts[4] - chrStarts[3]), rep(4, chrStarts[5] - chrStarts[4]), rep(5, nrow(g60snps) - chrStarts[5] + 1))

sumz <- rowSums(g60snps[,grep('X', colnames(g60snps))])
g60snps_05 <- g60snps[!(sumz < 3 | sumz >  57),]

#go through all the stand compositions, and calc allele freq for each


plotcom <- unique(fieldD[,paste0('lines', 1:20)])


csz <- combC2b$CS_number[match(sub('X', '', colnames(g60snps_05[,-1])), combC2b$ecotypeid)]


af <- c()
for(i in 1:nrow(plotcom)){
	if(plotcom[i,'lines2'] != ''){
		 af <- cbind(af, rowMeans(g60snps_05[,na.omit(match(plotcom[i,], csz) + 1)]))
		 }else{
		 	 af <- cbind(af, (g60snps_05[,na.omit(match(plotcom[i,], csz) + 1)]))

		 }
	}


#test vs af and var = af*(1-af) vs diversity effect and yield
avar <- af*(1-af)


#map rows of fieldD to columns of af/rows of plotcom
fieldD$plotcomRow <- NA
for(i in 1:nrow(fieldD)) for(j in 1:nrow(plotcom)) if(all.equal(fieldD[i,paste0('lines', 1:20)], plotcom[j,], check.attributes = F)[1] == T) fieldD$plotcomRow[i] <- j

fieldD$plotcomRow <- factor(fieldD$plotcomRow)

meanBiomC <- tapply((fieldD$FH_Wt/fieldD$MaxPlantNum)[fieldD$trt == 'C' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], fieldD$plotcomRow[fieldD$trt == 'C' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)

meanBiomS <- tapply((fieldD$FH_Wt/fieldD$MaxPlantNum)[fieldD$trt == 'S' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], fieldD$plotcomRow[fieldD$trt == 'S' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)



meanBiom <- tapply((fieldD$FH_Wt/fieldD$MaxPlantNum)[  fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], fieldD$plotcomRow[  fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)
meanBiom[is.na(meanBiomC)] <- NA
meanBiom[is.na(meanBiomS)] <- NA


meanDivEffect_C <- tapply(fieldD$DivEffect[fieldD$trt == 'C' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], fieldD$plotcomRow[fieldD$trt == 'C' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)

meanDivEffect_S <- tapply(fieldD$DivEffect[fieldD$trt == 'S' & fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged'], fieldD$plotcomRow[fieldD$trt == 'S' & fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)




#now make kinship matrices based on af and avar

af_kmat <- matrix(NA, nrow = ncol(af), ncol = ncol(af))

for(i in 1:(ncol(af)-1)){
	for(j in (i+1):ncol(af)){
		af_kmat[i,j] <- mean(1-abs(af[,i] - af[,j]))
		}
	cat(i, ' ')
	}

diag(af_kmat) <- 1
af_kmat[lower.tri(af_kmat)] <- t(af_kmat)[lower.tri(af_kmat)]



source('emma.v1_fixed.R') #EMMA scripts
source("emma.REML.t_fixed.R")

rs<-emma.REML.t(ys = t(meanBiomC), xs = af, K = af_kmat)

save(rs, file = 'BiomControl_vs_AF_inclMono_EMMA_14Apr2020.RDA')

rs<-emma.REML.t(ys = t(meanBiomS), xs = af, K = af_kmat)

save(rs, file = 'BiomStress_vs_AF_inclMono_EMMA_14Apr2020.RDA')




#find low maf
mafz <- rowMeans(af[,!is.na(meanDivEffect_C)])
tmpG <- mafz >= 0.05 & mafz <= 0.95
tmpG[tmpG == T][1275858] <- F #exclude SNP throwing error


rs<-emma.REML.t(ys = t(meanDivEffect_C[!is.na(meanDivEffect_C)]), xs = af[tmpG, !is.na(meanDivEffect_C)], K = af_kmat[!is.na(meanDivEffect_C), !is.na(meanDivEffect_C)])

save(rs, file = paste0('DivEffect_Control_vs_AF_EMMA_14Apr2020.RDA'))


	


rs<-emma.REML.t(ys = t(meanDivEffect_S), xs = af, K = af_kmat)

save(rs, file = 'DivEffect_Stress_vs_AF_EMMA_14Apr2020.RDA')




#find low maf
tmp <- rowMeans(af[,!is.na(meanDivEffect_C)])


source("ManhattanPlotFunctions.R")

rs2 <- rs
rs2$ps <- rs$ps[tmp >= 0.05 & tmp <= 0.95]



load('DivEffect_Control_vs_AF_EMMA_14Apr2020.RDA')

snp<- which(tmpG == T)[which(rs$ps == min(rs$ps))]



########### check association with flowering time
ft10 <- read.csv('FT10values.csv')
ft16 <- read.csv('FT16values.csv')

t.test(ft10$phenotype_value[paste0('X', ft10$accession_id) %in% colnames(g60snps_05)] ~ unlist(g60snps_05[snp, paste0('X', ft10$accession_id)[paste0('X', ft10$accession_id) %in% colnames(g60snps_05)]]))

t.test(ft16$phenotype_value[paste0('X', ft16$accession_id) %in% colnames(g60snps_05)] ~ unlist(g60snps_05[snp, paste0('X', ft16$accession_id)[paste0('X', ft16$accession_id) %in% colnames(g60snps_05)]]))



##convert FT to numeric
fieldD$FT1 <- NA
fieldD$FT50 <- NA
for (i in 1:nrow(fieldD)){
	
	fieldD$FT1[i] <- which(fieldD[i,grep('FT\\.', colnames(fieldD))] > 0)[1]
	fieldD$FT50[i] <- which(fieldD[i,grep('FT\\.', colnames(fieldD))] > 49)[1]
	
}

fieldD$FT1d <- as.numeric(as.Date(c('5/29/18', '6/4/18', '6/8/18', '6/12/18', '6/15/18', '6/19/18', '6/21/18', '6/26/18', '6/29/18', '7/6/18', '7/10/18')[fieldD$FT1], format = '%m/%d/%y'))
fieldD$FT50d <- as.numeric(as.Date(c('5/29/18', '6/4/18', '6/8/18', '6/12/18', '6/15/18', '6/19/18', '6/21/18', '6/26/18', '6/29/18', '7/6/18', '7/10/18')[fieldD$FT50], format = '%m/%d/%y'))

max(fieldD$FT1d, na.rm = T)

fieldD$FT1d_InclNoFlow <- fieldD$FT1d
fieldD$FT1d_InclNoFlow[is.na(fieldD$FT1d_InclNoFlow)] <- max(fieldD$FT1d, na.rm = T) + 2

meanFT_C<- tapply(fieldD$FT1d_InclNoFlow[fieldD$trt == 'C' ], fieldD$plotcomRow[fieldD$trt == 'C'], mean, na.rm = T)

summary(lm(meanFT_C ~ af[snp, ]))





### now do variance ##

avar_kmat <- matrix(NA, nrow = ncol(avar), ncol = ncol(avar))

for(i in 1:(ncol(avar)-1)){
	for(j in (i+1):ncol(avar)){
		avar_kmat[i,j] <- mean(1-abs(avar[,i] - avar[,j]))
		}
	cat(i, ' ')
	}

diag(avar_kmat) <- 1
avar_kmat[lower.tri(avar_kmat)] <- t(avar_kmat)[lower.tri(avar_kmat)]




rs<-emma.REML.t(ys = t(meanBiomC), xs = avar, K = avar_kmat)

save(rs, file = 'BiomControl_vs_AVar_inclMono_EMMA_16Apr2020.RDA')

rm('rs')



mafz <- rowMeans(avar[,!is.na(meanBiomS)])
tmpG <- mafz >= 0.05 & mafz <= 0.95
tmpG[tmpG == T][790347] <- F #exclude SNP throwing error


rs<-emma.REML.t(ys = t(meanBiomS), xs = avar[tmpG,], K = avar_kmat)


save(rs, file = 'BiomStress_vs_AVar_inclMono_EMMA_16Apr2020.RDA')

rm('rs')



rs<-emma.REML.t(ys = t(meanDivEffect_C), xs = avar, K = avar_kmat)

save(rs, file = 'DivEffect_Control_vs_AVar_EMMA_16Apr2020.RDA')


rm('rs')


rs<-emma.REML.t(ys = t(meanDivEffect_S), xs = avar, K = avar_kmat)

save(rs, file = 'DivEffect_Stress_vs_AVar_EMMA_16Apr2020.RDA')




####Plot main focal result ####

#recalc the MAF
mzC <- rowMeans(af[,!is.na(meanBiomC)]) #doing it like this includes the monocultures
mafC <- mzC
mafC[mafC > 0.5] <- 1-mafC[mafC > 0.5]

mzCdiv <- rowMeans(af[,!is.na(meanDivEffect_C)])
mafCdiv <- mzCdiv
mafCdiv[mafCdiv > 0.5] <- 1-mafCdiv[mafCdiv > 0.5]

mzS <- rowMeans(af[,!is.na(meanBiomS)])
mafS <- mzS
mafS[mafS > 0.5] <- 1-mafS[mafS > 0.5]

mzSdiv <- rowMeans(af[,!is.na(meanDivEffect_S)])
mafSdiv <- mzSdiv
mafSdiv[mafSdiv > 0.5] <- 1-mafSdiv[mafSdiv > 0.5]




source("ManhattanPlotFunctions.R")




load( 'DivEffect_Control_vs_AVar_EMMA_16Apr2020.RDA')





rs2 <- rs
rs2$ps <- rs$ps[mafCdiv >= 0.05]

man.h(rs = rs2, SNP_pos = g60snps_05[mafCdiv >= 0.05,c('chr', 'snp')], plot.file.name = 'DivEffect_Control_vs_AVar_EMMA_ManhattanPlot_21Jul2020.pdf',  pz = 10, add = F, high.loc = c(2, 9584630, 9584630), pvalthresh = max(rs2$ps[p.adjust(rs2$ps, method = 'fdr') < 0.1]))#, xlimz = NULL, ylimz = NULL, hlines = NULL, hlty = 2, bz = 10)




snp<- which(rs$ps == min(rs$ps))

pdf('DivEffect_Control_vs_AVar_inclMono_EMMA_topSNPboxplot_18Apr2020.pdf', pointsize = 21)
boxplot(meanDivEffect_C ~ avar[snp, ], varwidth = T, range = 0, xlab = 'Stand var(allele frequency)', ylab = 'DME (g)')
title(main = 'Chr. 2, 9584630 bp', col.main = brewer.pal(9, 'Reds')[7])
dev.off()

########### make plot for supplement
#

#af, a var, vs biom, dme in control and stress. 8 panels. 2 columns

###### FT candidates
ftcand <- read.csv('FTgenes.csv')
coords <- read.csv('TAIR10_loci_coords.csv')

coordsFT <- coords[match(unique(ftcand$GeneID), coords$locus),]

ftcanSNP <- c()


for(i in 1:5){
	tmpcoords <- coordsFT[coordsFT$chr == i,]
	tmpcoords <- tmpcoords[order(tmpcoords$start),]
	tmpcoords$start <- tmpcoords$start - 3e3
	tmpcoords$stop <- tmpcoords$start + 3e3

	tmp <- cut(g60snps_05$snp[g60snps_05$chr == i], c(t(tmpcoords[,c('start','stop')])), labels = F, include.lowest = T) #odds are genic #added 2 so that 2
	tmp <- (tmp %% 2) #== 0
	tmp[is.na(tmp)] <- 0
	
	ftcanSNP <- c(ftcanSNP, as.character(g60snpid[g60snps_05$chr == i][tmp == 1]))
	cat(length(ftcanSNP), '\n')

	}


ftcanTbl <- c()	




load( 'BiomControl_vs_AF_inclMono_EMMA_14Apr2020.RDA')

rs2 <- rs
rs2$ps <- rs$ps[mafC >= 0.05]


man.h(rs = rs2, SNP_pos = g60snps_05[mafC >= 0.05,c('chr', 'snp')], plot.file.name = 'BiomControl_vs_AFreq_inclMono_EMMA_ManhattanPlot_ForSupp_18Apr2020.pdf')#, xlimz = NULL, ylimz = NULL, hlines = NULL, hlty = 2, bz = 10)




tmpft <- data.frame(
			resp = 'Biomass', 
			trt = 'high resource', 
			cov = 'Allele frequency',
			g60snps_05[mafC >= 0.05,][g60snpid[mafC >= 0.05] %in% ftcanSNP,c('chr', 'snp')],
			p = rs2$ps[g60snpid[mafC >= 0.05] %in% ftcanSNP],
			q = p.adjust(rs2$ps[g60snpid[mafC >= 0.05] %in% ftcanSNP], method = 'fdr')
			)

ftcanTbl <- rbind(ftcanTbl, 
					tmpft
				)	





load( 'BiomStress_vs_AF_inclMono_EMMA_14Apr2020.RDA')


rs2 <- rs
rs2$ps <- rs$ps[mafS >= 0.05]

man.h(rs = rs2, SNP_pos = g60snps_05[mafS >= 0.05,c('chr', 'snp')], plot.file.name = 'BiomStress_vs_AFreq_inclMono_EMMA_ManhattanPlot_ForSupp_18Apr2020.pdf')#, xlimz = NULL, ylimz = NULL, hlines = NULL, hlty = 2, bz = 10)






tmpft <- data.frame(
			resp = 'Biomass', 
			trt = 'low resource', 
			cov = 'Allele frequency',
			g60snps_05[mafS >= 0.05,][g60snpid[mafS >= 0.05] %in% ftcanSNP,c('chr', 'snp')],
			p = rs2$ps[g60snpid[mafS >= 0.05] %in% ftcanSNP],
			q = p.adjust(rs2$ps[g60snpid[mafS >= 0.05] %in% ftcanSNP], method = 'fdr')
			)

ftcanTbl <- rbind(ftcanTbl, 
					tmpft
				)	







load( 'DivEffect_Control_vs_AF_EMMA_14Apr2020.RDA')

mafz <- rowMeans(af[,!is.na(meanDivEffect_C)])
tmpG <- mafz >= 0.05 & mafz <= 0.95
tmpG[tmpG == T][1275858] <- F #exclude SNP throwing error



man.h(rs = rs, SNP_pos = g60snps_05[tmpG,c('chr', 'snp')], plot.file.name = 'DivEffectControl_vs_AFreq_inclMono_EMMA_ManhattanPlot_ForSupp_18Apr2020.pdf')#, xlimz = NULL, ylimz = NULL, hlines = NULL, hlty = 2, bz = 10)







tmpft <- data.frame(
			resp = 'DME', 
			trt = 'high resource', 
			cov = 'Allele frequency',
			g60snps_05[tmpG,][g60snpid[tmpG] %in% ftcanSNP,c('chr', 'snp')],
			p = rs$ps[g60snpid[tmpG] %in% ftcanSNP],
			q = p.adjust(rs$ps[g60snpid[tmpG] %in% ftcanSNP], method = 'fdr')
			)

ftcanTbl <- rbind(ftcanTbl, 
					tmpft
				)	







load( 'DivEffect_Stress_vs_AF_EMMA_14Apr2020.RDA')


rs2 <- rs
rs2$ps <- rs$ps[mafSdiv >= 0.05]

man.h(rs = rs2, SNP_pos = g60snps_05[mafSdiv >= 0.05,c('chr', 'snp')], plot.file.name = 'DivEffectStress_vs_AFreq_inclMono_EMMA_ManhattanPlot_ForSupp_18Apr2020.pdf')#, xlimz = NULL, ylimz = NULL, hlines = NULL, hlty = 2, bz = 10)





tmpft <- data.frame(
			resp = 'DME', 
			trt = 'low resource', 
			cov = 'Allele frequency',
			g60snps_05[mafSdiv >= 0.05,][g60snpid[mafSdiv >= 0.05] %in% ftcanSNP,c('chr', 'snp')],
			p = rs2$ps[g60snpid[mafSdiv >= 0.05] %in% ftcanSNP],
			q = p.adjust(rs2$ps[g60snpid[mafSdiv >= 0.05] %in% ftcanSNP], method = 'fdr')
			)

ftcanTbl <- rbind(ftcanTbl, 
					tmpft
				)	







load( 'BiomControl_vs_AVar_inclMono_EMMA_16Apr2020.RDA')


rs2 <- rs
rs2$ps <- rs$ps[mafC > 0.05] #

man.h(rs = rs2, SNP_pos = g60snps_05[mafC > 0.05,c('chr', 'snp')], plot.file.name = 'BiomControl_vs_AVar_inclMono_EMMA_ManhattanPlot_ForSupp_27Jul2020.pdf')#, xlimz = NULL, ylimz = NULL, hlines = NULL, hlty = 2, bz = 10)





tmpft <- data.frame(
			resp = 'biomass', 
			trt = 'high resource', 
			cov = 'Allele frequency variance',
			g60snps_05[mafC >= 0.05,][g60snpid[mafC >= 0.05] %in% ftcanSNP,c('chr', 'snp')],
			p = rs2$ps[g60snpid[mafC >= 0.05] %in% ftcanSNP],
			q = p.adjust(rs2$ps[g60snpid[mafC >= 0.05] %in% ftcanSNP], method = 'fdr')
			)

ftcanTbl <- rbind(ftcanTbl, 
					tmpft
				)	








load( 'BiomStress_vs_AVar_inclMono_EMMA_16Apr2020.RDA')

mafz <- rowMeans(avar[,!is.na(meanBiomS)])
tmpG <- mafz >= 0.05 & mafz <= 0.95
tmpG[tmpG == T][790347] <- F #exclude SNP throwing error


man.h(rs = rs, SNP_pos = g60snps_05[tmpG,c('chr', 'snp')], plot.file.name = 'BiomStress_vs_AVar_inclMono_EMMA_ManhattanPlot_ForSupp_19Apr2020.pdf')#, xlimz = NULL, ylimz = NULL, hlines = NULL, hlty = 2, bz = 10)




tmpft <- data.frame(
			resp = 'biomass', 
			trt = 'low resource', 
			cov = 'Allele frequency variance',
			g60snps_05[tmpG,][g60snpid[tmpG] %in% ftcanSNP,c('chr', 'snp')],
			p = rs$ps[g60snpid[tmpG] %in% ftcanSNP],
			q = p.adjust(rs$ps[g60snpid[tmpG] %in% ftcanSNP], method = 'fdr')
			)

ftcanTbl <- rbind(ftcanTbl, 
					tmpft
				)	







load( 'DivEffect_Control_vs_AVar_EMMA_16Apr2020.RDA')


rs2 <- rs
rs2$ps <- rs$ps[mafCdiv >= 0.05]

man.h(rs = rs2, SNP_pos = g60snps_05[mafCdiv >= 0.05,c('chr', 'snp')], plot.file.name = 'DivEffectControl_vs_AVar_inclMono_EMMA_ManhattanPlot_ForSupp_19Apr2020.pdf')#, xlimz = NULL, ylimz = NULL, hlines = NULL, hlty = 2, bz = 10)





tmpft <- data.frame(
			resp = 'DME', 
			trt = 'high resource', 
			cov = 'Allele frequency variance',
			g60snps_05[mafCdiv >= 0.05,][g60snpid[mafCdiv >= 0.05] %in% ftcanSNP,c('chr', 'snp')],
			p = rs2$ps[g60snpid[mafCdiv >= 0.05] %in% ftcanSNP],
			q = p.adjust(rs2$ps[g60snpid[mafCdiv >= 0.05] %in% ftcanSNP], method = 'fdr')
			)

ftcanTbl <- rbind(ftcanTbl, 
					tmpft
				)	







load( 'DivEffect_Stress_vs_AVar_EMMA_16Apr2020.RDA')



rs2 <- rs
rs2$ps <- rs$ps[mafSdiv >= 0.05]

man.h(rs = rs2, SNP_pos = g60snps_05[mafSdiv >= 0.05,c('chr', 'snp')], plot.file.name = 'DivEffectStress_vs_AVar_inclMono_EMMA_ManhattanPlot_ForSupp_19Apr2020.pdf')#, xlimz = NULL, ylimz = NULL, hlines = NULL, hlty = 2, bz = 10)



tmpft <- data.frame(
			resp = 'DME', 
			trt = 'low resource', 
			cov = 'Allele frequency variance',
			g60snps_05[mafSdiv >= 0.05,][g60snpid[mafSdiv >= 0.05] %in% ftcanSNP,c('chr', 'snp')],
			p = rs2$ps[g60snpid[mafSdiv >= 0.05] %in% ftcanSNP],
			q = p.adjust(rs2$ps[g60snpid[mafSdiv >= 0.05] %in% ftcanSNP], method = 'fdr')
			)

ftcanTbl <- rbind(ftcanTbl, 
					tmpft
				)	


write.csv(ftcanTbl, 'FT_candSNPasscnTbl.csv', row.names =F)





### annotate genes back in to ftcanTbl


ftcanTbl$gene <- NA


for(i in 1:nrow(ftcanTbl)){

	tmpc <- coordsFT[coordsFT$chr == ftcanTbl$chr[i],]

	ftcanTbl$gene[i] <- as.character(tmpc$locus[which(abs(ftcanTbl$snp[i] - tmpc$stop) == min(c(abs(ftcanTbl$snp[i] - tmpc$stop), abs(ftcanTbl$snp[i] - tmpc$start))) | abs(ftcanTbl$snp[i] - tmpc$start) == min(c(abs(ftcanTbl$snp[i] - tmpc$stop), abs(ftcanTbl$snp[i] - tmpc$start))))])

	}






write.csv(ftcanTbl, 'FT_candSNPasscnTbl2.csv', row.names =F)




