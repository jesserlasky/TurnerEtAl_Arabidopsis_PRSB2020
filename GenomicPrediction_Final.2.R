
fieldD <- read.csv('CombinedFieldData_24Sep2020.csv', as.is = T)




#get mean per plant biomass in monos
biomC <- tapply(fieldD$FH_Wt[fieldD$plotType == 'mono' & fieldD$trt == 'C' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged']/fieldD$MaxPlantNum[fieldD$plotType == 'mono' & fieldD$trt == 'C'& fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], fieldD$ecotypeid[fieldD$plotType == 'mono' & fieldD$trt == 'C'& fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)

biomS <- tapply(fieldD$FH_Wt[fieldD$plotType == 'mono' & fieldD$trt == 'S'& fieldD$MaxPlantNum > 10  & fieldD$Damaged_pot != 'damaged']/fieldD$MaxPlantNum[fieldD$plotType == 'mono' & fieldD$trt == 'S' & fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], fieldD$ecotypeid[fieldD$plotType == 'mono' & fieldD$trt == 'S'& fieldD$MaxPlantNum > 10 & fieldD$Damaged_pot != 'damaged'], mean, na.rm = T)


library(rrBLUP)
library(rhdf5)  # load library for rhdf5 files


kmatrix <- H5Fopen("kinship_ibs_binary_mac5.h5py") # kinship matrix

all.equal(names(biomC), names(biomS))


kmat <- kmatrix$kinship[match(names(biomC), kmatrix$accessions), match(names(biomC), kmatrix$accessions)]
#rownames of    K must match the genotype labels in the data frame for phenotyped lines; missing phenotypes (NA) are simply omitted.
rownames(kmat) <- kmatrix$accessions[match(names(biomC), kmatrix$accessions)]

tmp_acc <- data.frame(ecotype_id = names(biomC), biomC = biomC, biomS = biomS)


#cross validation
set.seed(100)


foldz <- 10

repz <- 1000

corz <- matrix(nrow = repz, ncol = 2)

for(j in 1:repz){

	conK <- rep(NA, length(biomC))
	strK <- conK

	reco <- sample(1:length(biomC))

	for(i in 1:foldz){
		tmp_choose <- na.omit(reco[(1):(ceiling(length(biomC)/foldz)) + (i-1)* ceiling(length(biomC)/foldz)])

		
		tmp_accNA <- tmp_acc
		
		tmp_accNA[tmp_choose, 'biomC'] <- NA
		tmp_accNA[tmp_choose, 'biomS'] <- NA
		
	#	tmp_accNA <- tmp_accNA[c((1:length(biomC))[!(1:length(biomC)) %in% tmp_choose] , tmp_choose),]

		tmpC <- kin.blup(tmp_accNA, geno = 'ecotype_id', pheno = 'biomC', K = kmat)
		tmpS <- kin.blup(tmp_accNA, geno = 'ecotype_id', pheno = 'biomS', K = kmat)

		conK[tmp_choose] <- tmpC$pred[tmp_choose]
		strK[tmp_choose] <- tmpS$pred[tmp_choose]

	#cat(i)
	}
	cat(j)
	corz[j,1] <- cor(conK, biomC, use = 'pa')
	corz[j,2] <- cor(strK, biomS, use = 'pa')
}

colMeans(corz) 






#fecundity (avg per plant)

#get mean per plant fec in monos

fecC <- (tapply(fieldD$potF2[fieldD$plotType == 'mono' & fieldD$trt == 'C'], fieldD$ecotypeid[fieldD$plotType == 'mono' & fieldD$trt == 'C'], mean, na.rm = T))

fecS <- (tapply(fieldD$potF2[fieldD$plotType == 'mono' & fieldD$trt == 'S'], fieldD$ecotypeid[fieldD$plotType == 'mono' & fieldD$trt == 'S'], mean, na.rm = T))


all.equal(names(fecC), names(fecS))

kmat <- kmatrix$kinship[match(names(fecC), kmatrix$accessions), match(names(fecC), kmatrix$accessions)]

rownames(kmat) <- kmatrix$accessions[match(names(fecC), kmatrix$accessions)]


#cross validation

tmp_acc <- data.frame(ecotype_id = names(fecC), fecC = fecC, fecS = fecS)



foldz <- 10

repz <- 1000

corz <- matrix(nrow = repz, ncol = 2)

for(j in 1:repz){


conK <- rep(NA, length(fecC))
strK <- conK

reco <- sample(1:length(fecC))


for(i in 1:foldz){
	tmp_choose <- na.omit(reco[(1):(ceiling(length(fecC)/foldz)) + (i-1)* ceiling(length(fecC)/foldz)])

	
	tmp_accNA <- tmp_acc
	
	tmp_accNA[tmp_choose, 'fecC'] <- NA
	tmp_accNA[tmp_choose, 'fecS'] <- NA
	

	tmpC <- kin.blup(tmp_accNA, geno = 'ecotype_id', pheno = 'fecC', K = kmat)
	tmpS <- kin.blup(tmp_accNA, geno = 'ecotype_id', pheno = 'fecS', K = kmat)

	conK[tmp_choose] <- tmpC$pred[tmp_choose]
	strK[tmp_choose] <- tmpS$pred[tmp_choose]

}
cat(j)

corz[j,1] <- cor(conK, fecC, use = 'pa')
corz[j,2] <- cor(strK, fecS, use = 'pa')
}

colMeans(corz) 



########## plotting ##########

load( 'PairwiseDistMatrix_SNPrelate_6Sep2018.RDA') #Pairwise genetic distance matrix from 1001 genomes SNPs

#unrooted plot

g1001 <- read.csv('1001genomesaccessionsInfo.csv')

library(RColorBrewer)
mepal <- colorRampPalette(brewer.pal(9, 'YlGn')[1:9])

library(SNPRelate)
library(phyclust)





pdf('NJ_Tree_ReacNorm_Biom_28Sep2020.pdf', height = 2.7)


par(plt = c(0.1, 0.29, 0.28, 0.8))

matplot(x = t(matrix(1:2, byrow = T, nrow = length(biomC), ncol = 2)), y = t(cbind(biomC, biomS)) * 1000, type = 'l', col = gray(0.5), lty = 1, xaxt = 'n', xlab = '', ylab = '', lwd = 0.5)

points(rep(1, length(biomC)), biomC * 1000, pch = 19, col = (mepal(55))[round((biomC)*1000) + 5], cex = 0.2)

points(rep(2, length(biomS)), biomS * 1000, pch = 19, col = (mepal(55))[round((biomS)*1000) + 5], cex = 0.2)


axis(1, at = 1:2, labels = c('High', 'Low'))
mtext('Resources', side = 1, line = 2.5)
mtext('AGB per plant (mg)', side = 2, line = 2.5)

lines(1:2, c(mean(biomC, na.rm = T), mean(biomS, na.rm = T)) * 1000, lwd = 2, lty = 2)




par(new = T)

par(plt = c(0.32, 0.65, 0.27, 0.9))
plotnj(nj(gDis$diss[!gDis$sample.id %in% c('9600', '9817'), !gDis$sample.id %in% c('9600', '9817')]), show.tip.label =F, main = '') #tip.color = terrain.colors(200)[(combP$year-1866)]
mtext('High resource', line = 0.2)

#9600 and #9817 are not here - very poor germination for both
tiplabels(pch = 19, col =  (mepal(55))[round((biomC[as.character(gDis$sample.id[!gDis$sample.id %in% c('9600', '9817')])])*1000) + 5], bg = 'white', lwd = 2, cex = 1.5)

tiplabels(text = g1001$country[match(gDis$sample.id[!gDis$sample.id %in% c('9600', '9817')], g1001$tg_ecotypeid)], frame = 'none', bg = NA, cex = 0.3, col = gray(0.2))



par(new = T)
par(plt = c(0.67, 1, 0.27, 0.9))

plotnj(nj(gDis$diss[!gDis$sample.id %in% c('9600', '9817'), !gDis$sample.id %in% c('9600', '9817')]), show.tip.label =F, main = '') #tip.color = terrain.colors(200)[(combP$year-1866)]
mtext('Low resource', line = 0.2)
mtext('AGB/plant (mg)', line = 1.5, side = 1, at = 0.4)


tiplabels(pch = 19, col =  (mepal(55))[round((biomS[as.character(gDis$sample.id[!gDis$sample.id %in% c('9600', '9817')])])*1000) +5 ], bg = 'white', lwd = 2, cex = 1.5)

tiplabels(text = g1001$country[match(gDis$sample.id[!gDis$sample.id %in% c('9600', '9817')], g1001$tg_ecotypeid)], frame = 'none', bg = NA, cex = 0.3, col = gray(0.2))

range(c(biomC, biomS),na.rm = T) 

min(sqrt(c(biomC, biomS)), na.rm = T) + diff(range(sqrt(c(biomC, biomS)),na.rm = T))/2#halfway on the sqrt scale
(min(sqrt(c(biomC, biomS)), na.rm = T) + diff(range(sqrt(c(biomC, biomS)),na.rm = T))/2)^2


par(new = T)
par(plt = c(0.55, 0.67, 0.15, 0.2))
image(t(t(1:100)), col = (mepal(220))[c(20:192)], xaxt = 'n', yaxt = 'n', bty = 'n')
axis(1, at = c(0, 0.5, 1), las = 1, labels = c(0.5, 20, 40.0))




dev.off()








