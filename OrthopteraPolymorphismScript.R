##################################################################################
# Data preparation
##################################################################################

rm(list=ls())

# Settings
# Since the different sections are independent, this needs to be copy-pasted if required
path = "MYPATH\\OrthopteraPolymorphism\\"
version = "200326"
outfolder = "Figures\\"
filename = "Data\\OrthopteraPolymorphism"
cols=c(Brown=rgb(126,99,80, maxColorValue=255), Green=rgb(103,175,102, maxColorValue=255))
outToFile = TRUE

# Load data
mdnames=read.table(paste0(path, filename, "_Names_", version, ".txt"), header=TRUE, sep="\t")
mdmorphs=read.table(paste0(path, filename, "_Morphs_", version, ".txt"), header=TRUE, sep="\t")
mdhabitat = read.table(paste0(path, filename, "_Habitat_", version, ".txt"), header=TRUE, sep="\t", na.strings = ".")
mddist = read.table(paste0(path, filename, "_Distribution_", version, ".txt"), header=TRUE, sep="\t", na.strings = ".")

# Remove unnecesary columns
mdnames = mdnames[,-which(colnames(mdnames) %in% c("Synonyms", "AlternativeGenus", "Subspecies", "RemarksTaxonomy"))]
mdmorphs = mdmorphs[,-which(colnames(mdmorphs) %in% c("Species", "OtherVariants", "RemarksMorphs", "SourceMorphs"))]
mdhabitat = mdhabitat[,-which(colnames(mdhabitat) %in% c("Species", "RemarksHabitat"))]
mddist = mddist[,-which(colnames(mddist) %in% c("Species", "RemarksDistribution"))]

# Merging data
md = merge(mdnames, mdmorphs)
md = merge(md, mddist)
md = merge(md, mdhabitat)

# Additions of some useful columns
md$Polymorphic = factor(md$Polymorphic)
md$CricketsYN = as.numeric(md$Superfamily=="Grylloidea" | md$Superfamily=="Gryllotalpoidea" | md$Superfamily=="Rhaphidophoroidea")
md$Abundance = md$Brown
md$Abundance[md$Green == "absent" ] = "absent"
md$Abundance[md$Green == "exceptional"] = "exceptional"
md$Abundance[md$Green == "rare"] = "rare"
md$Abundance[md$Green == "regular"] = "regular"
md$Abundance[md$Green == "common"] = "common"
md$Abundance = factor(md$Abundance, levels=c("absent", "exceptional", "rare", "regular", "common"))
md$Polymorphic2 = md$Polymorphic
md$Polymorphic2[md$Green == "exceptional"] = "brown"
md$Polymorphic2[md$Brown == "exceptional"] = "green"

write.table(md, file = paste0(path, filename, "_Data_", version, ".txt"), row.names=FALSE)


##################################################################################
# Overall summaries
##################################################################################

rm(list=ls())

# Settings
path = "MYPATH\\OrthopteraPolymorphism\\"
version = "200326"
outfolder = "Figures\\"
filename = "Data\\OrthopteraPolymorphism"
cols = c("sandybrown", "darkolivegreen3", "purple")
outToFile = TRUE

# Load data
md = read.table(file = paste0(path, filename, "_Data_", version, ".txt"), header=TRUE)

# Some preparation
md$PolymorphismData = !is.na(md$Combined)
md$HabitatData = !is.na(md$H01_Marshes_and_reeds)
md$HabitatUseData = !is.na(md$U01_Bushes_and_trees)
md$ActivityData = !is.na(md$U06_Diurnal)
md$SeasonalityData = !is.na(md$U09_Spring)
md$DistributionData = !is.na(md$R01_Britain)

# Available data (Number of species)
nrow(md) # Number of species
sum(md$PolymorphismData)   # Number of species with polymorphism data
sum(md$DistributionData) # Number of species with distribution data
sum(md$HabitatData) # Number of species with habitat data
sum(md$HabitatUseData) # Number of species with habitat use data
sum(md$ActivityData) # Number of species with activity data
sum(md$SeasonalityData) # Number of species with seasonality data

# Available data (Proportions)
mean(md$PolymorphismData)  # Proportion of species with polymorphism data
mean(md$DistributionData) # Proprtion of species with distribution data
mean(md$HabitatData) # Proportion of species with habitat data
mean(md$HabitatUseData) # Proportion of species with habitat use data
mean(md$ActivityData) # Proportion of species with activity data
mean(md$SeasonalityData) # Proportion of species with seasonality data

# Overall proportion of polymorphic states (all species)
polymorph = table(md$Polymorphic[md$PolymorphismData])
polymorph # Number of species
polymorph/sum(polymorph) # Proportion of species

# Overall proportion of polymorphic states (without crickets)
sum(md$CricketsYN) # Number of cricket species
polymorph = table(md$Polymorphic[md$PolymorphismData & md$CricketsYN==0])
polymorph # Number of species
polymorph/sum(polymorph) # Proportion of species

# Proportion of polymorphic in Ensifera and Caelifera
polymorph = table(list(md$Suborder[md$PolymorphismData], md$Polymorphic[md$PolymorphismData]))
polymorph # Number of species
polymorph/rowSums(polymorph) # Proportion of species

# Polymorphism across subfamilies
polymorph = as.data.frame.matrix(table(list(factor(md$Subfamily[md$PolymorphismData & md$CricketsYN==0]), md$Polymorphic[md$PolymorphismData & md$CricketsYN==0])))
polymorph # Numbers by subfamily (besides crickets)
polymorph[polymorph[,"green"]==0 & polymorph[,"polymorphic"]==0,] # Exclusively brown subfamilies besides crickets

# Overview of relative abundances in polymorphic species
abundance = table(factor(md$Abundance[md$Abundance!="absent"]))
abundance # Number of species
abundance / sum(abundance) # Proportion of species

# Overview of relative abundances in polymorphic species
overall = table(md$Polymorphic[!is.na(md$Polymorphic2)])

# Overview of cases where one of the morphs is rare
# Abundance of the green morph is named first
rare = table(factor(md$Combined[md$Abundance=="rare"]))
sum(rare)
rare # Number of species
rare/sum(rare) # Proportion of species
chisq.test(rare, p=c(0.5,0.5))
# Overview of cases where one of the morphs is exception
# Abundance of the green morph is named first
exceptional = table(factor(md$Combined[md$Abundance=="exceptional"]))
sum(exceptional)
exceptional # Number of species
exceptional/sum(exceptional) # Proportion of species
chisq.test(exceptional, p=c(0.5,0.5))

# Overview of cases where there is only one morph
# Abundance of the green morph is named first
absent = table(factor(md$Combined[md$Abundance=="absent"]))
sum(absent)
absent # Number of species
absent/sum(absent) # Proportion of species
absent["absent-mono"]/absent["mono-absent"] # as a ratio

# Proportions within particular habitats
# Marshes and moist grasslands
md$TargetHabitatYN = md$H01_Marshes_and_reeds=="x" | md$H02_Moist_grassland=="x"
polymorph = table(factor(md$Polymorphic[md$TargetHabitatYN]))
polymorph/sum(polymorph)
# Grasslands overall
md$TargetHabitatYN = md$H01_Marshes_and_reeds=="x" | md$H02_Moist_grassland=="x" | md$H03_Mesic_grassland=="x" | md$H04_Dry_grassland=="x" | md$H05_Herbal_grassland=="x"
polymorph = table(factor(md$Polymorphic[md$TargetHabitatYN]))
polymorph/sum(polymorph)
# Alpine grasslands and shrubland
md$TargetHabitatYN = md$H13_Subalpine_shrubs=="x" | md$H14_Subalpine_grassland=="x"
polymorph = table(factor(md$Polymorphic[md$TargetHabitatYN]))
polymorph/sum(polymorph)

# Proportions of pied morphs
pied = as.data.frame.matrix(table(list(md$Subfamily, md$Pied)))
pied # Numbers by subfamily
pied/rowSums(pied) # Proportion by subfamily


##################################################################################
# Phylogenetic distribution (Figure 2)
##################################################################################

rm(list=ls())

# Settings
path = "MYPATH\\OrthopteraPolymorphism\\"
version = "200326"
outfolder = "Figures\\"
filename = "Data\\OrthopteraPolymorphism"
cols = c("sandybrown", "darkolivegreen3", "purple")
outToFile = TRUE

# Load data
md = read.table(file = paste0(path, filename, "_Data_", version, ".txt"), header=TRUE)

# Load packages
require(ape)

# Function for aggregation at different levels
aggFunct = function(md, level="Tribe", type="Polymorphic") {
	md$GreenOnly = ifelse(!is.na(md[,type]), ifelse(md[,type]=="green",1,0), 0)
	md$BrownOnly = ifelse(!is.na(md[,type]), ifelse(md[,type]=="brown",1,0), 0)
	md$PolymorphicOnly = ifelse(!is.na(md[,type]), ifelse(md[,type]=="polymorphic",1,0), 0)
	md$MissingData = as.numeric(is.na(md[,type]))
	md$AvailableData = as.numeric(!is.na(md[,type]))
	if(level=="Species") cols = 1:7
	if(level=="Genus") cols = 1:6
	if(level=="Tribe") cols = 1:5
	if(level=="Subfamily") cols = 1:4
	if(level=="Family") cols = 1:3
	if(level=="Superfamily") cols = 1:2
	if(level=="Suborder") cols = 1:1
	taxonlevels = c("Suborder", "Superfamily", "Family", "Subfamily", "Tribe", "Genus", "Species")
	if(level=="Species") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$AvailableData, md$MissingData), by=list(md$Suborder, md$Superfamily, md$Family, md$Subfamily, md$Tribe, md$Genus, md$Species), FUN=sum)
	if(level=="Genus") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$AvailableData, md$MissingData), by=list(md$Suborder, md$Superfamily, md$Family, md$Subfamily, md$Tribe, md$Genus), FUN=sum)
	if(level=="Tribe") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$AvailableData, md$MissingData), by=list(md$Suborder, md$Superfamily, md$Family, md$Subfamily, md$Tribe), FUN=sum)
	if(level=="Subfamily") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$AvailableData, md$MissingData), by=list(md$Suborder, md$Superfamily, md$Family, md$Subfamily), FUN=sum)
	if(level=="Family") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$AvailableData, md$MissingData), by=list(md$Suborder, md$Superfamily, md$Family), FUN=sum)
	if(level=="Superfamily") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$AvailableData, md$MissingData), by=list(md$Suborder, md$Superfamily), FUN=sum)
	if(level=="Suborder") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$AvailableData, md$MissingData), by=list(md$Suborder), FUN=sum)
	colnames(res) = c(taxonlevels[cols], "Green", "Brown", "Polymorphic", "AvailableData", "MissingData")
	if(level=="Species") index = order(res$Suborder, res$Superfamily, res$Family, res$Subfamily, res$Tribe, res$Genus, res$Species, decreasing = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), method="radix")
	if(level=="Genus") index = order(res$Suborder, res$Superfamily, res$Family, res$Subfamily, res$Tribe, res$Genus, decreasing = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), method="radix")
	if(level=="Tribe") index = order(res$Suborder, res$Superfamily, res$Family, res$Subfamily, res$Tribe, decreasing = c(TRUE, TRUE, TRUE, TRUE, TRUE), method="radix")
	if(level=="Subfamily") index = order(res$Suborder, res$Superfamily, res$Family, res$Subfamily, decreasing = c(TRUE, TRUE, TRUE, TRUE), method="radix")
	if(level=="Family") index = order(res$Suborder, res$Superfamily, res$Family, decreasing = c(TRUE, TRUE, TRUE), method="radix")
	if(level=="Superfamily") index = order(res$Suborder, res$Superfamily, decreasing = c(TRUE, TRUE), method="radix")
	if(level=="Suborder") index = order(res$Suborder, decreasing = c(TRUE), method="radix")
	res = res[index,]
	res$Name = res[,level]
	return(res)
}

# Aggregation at various levels
resSpec = aggFunct(md, "Species")
resGenus = aggFunct(md, "Genus")
resTribe = aggFunct(md, "Tribe")
resSubfam = aggFunct(md, "Subfamily")
resFam = aggFunct(md, "Family")
resSuperfam = aggFunct(md, "Superfamily")
resSuborder = aggFunct(md, "Suborder")

# Function for calculation of pie chart values
pieValues = function(resTribe, resGenus, resSubfam, resFam, resSubord, suborder=c("Ensifera", "Caelifera"), percent=TRUE) {
	resPie = rbind(resTribe[resTribe$Suborder==suborder,c("Name", "Green", "Brown", "Polymorphic")],
					resSubfam[resSubfam$Suborder==suborder,c("Name", "Green", "Brown", "Polymorphic")],
					resFam[resFam$Suborder==suborder,c("Name", "Green", "Brown", "Polymorphic")],
					resSuperfam[resSuperfam$Suborder==suborder,c("Name", "Green", "Brown", "Polymorphic")],
					resSuborder[resSuborder$Suborder==suborder,c("Name", "Green", "Brown", "Polymorphic")])
	resPie[,c("Green", "Brown", "Polymorphic")] / rowSums(resPie[,c("Green", "Brown", "Polymorphic")])
	resPie$Missing = 0
	resPie$Missing[resPie$Green==0 & resPie$Brown==0 & resPie$Polymorphic==0] = 1
	resPiePerc = round(resPie[,c("Green", "Brown", "Polymorphic", "Missing")] / rowSums(resPie[,c("Green", "Brown", "Polymorphic", "Missing")]),3)
	resPiePerc = cbind(resPiePerc[,1], resPiePerc[,2], resPiePerc[,3], resPiePerc[,4])
	ifelse(percent, return(resPiePerc), return(resPie))
}

# Calculation of pie chart values
resPieCaelifPerc = pieValues(resTribe, resGenus, resSubfam, resFam, resSubord, "Caelifera")
resPieEnsifPerc = pieValues(resTribe, resGenus, resSubfam, resFam, resSubord, "Ensifera")
resPieCaelif = pieValues(resTribe, resGenus, resSubfam, resFam, resSubord, "Caelifera", percent=FALSE)
resPieEnsif = pieValues(resTribe, resGenus, resSubfam, resFam, resSubord, "Ensifera", percent=FALSE)

# Function for writing the tree in Newick format
makeTree = function(res, suborder="all") {
	if(suborder!="all") res = subset(res, Suborder==suborder)
	tree = "(((("
	for(i in 2:nrow(res)) 
		if(res$Subfamily[i]==res$Subfamily[i-1]) tree=paste0(tree, ",")
		else
			if(res$Family[i]==res$Family[i-1]) tree=paste0(tree,"),(")
				else 
					if(res$Superfamily[i]==res$Superfamily[i-1]) tree=paste0(tree, ")),((")
						else 
							if(res$Suborder[i]==res$Suborder[i-1]) tree=paste0(tree,"))),(((")
								else tree=paste0(tree,")))),((((")
	tree = paste0(tree, "))));")
	return(tree)
}

# Function for node values
treeNodes = function(res, suborder="all") {
	if(suborder!="all") res = subset(res, Suborder==suborder)
	nodes = c(as.character(res$Suborder[1]), as.character(res$Superfamily[1]), as.character(res$Family[1]), as.character(res$Subfamily[1]))
	for(i in 2:nrow(res)) 
		if(res$Superfamily[i]!=res$Superfamily[i-1] | res$Family[i]!=res$Family[i-1] | res$Subfamily[i]!=res$Subfamily[i-1])
			if(res$Superfamily[i]!=res$Superfamily[i-1]) nodes=c(nodes, as.character(res$Superfamily[i]), as.character(res$Family[i]), as.character(res$Subfamily[i]))
				else
					if(res$Family[i]!=res$Family[i-1]) nodes=c(nodes, as.character(res$Family[i]), as.character(res$Subfamily[i]))
						else
							if(res$Subfamily[i]!=res$Subfamily[i-1]) nodes=c(nodes, as.character(res$Subfamily[i]))
	return(nodes)
}

# Function for tip values
treeTips = function(res, suborder="all") {
	if(suborder!="all") res = subset(res, Suborder==suborder)
	tips = rep("-",nrow(res)) 
	for(i in 1:nrow(res)) tips[i]=paste0(as.character(res$Tribe[i]), " (", res$AvailableData[i], "/", res$AvailableData+res$MissingData[i], ")")
	return(tips)
}

# Calculations for phylogenetic tree
treestrEnsif = makeTree(resTribe, suborder="Ensifera")
trEnsif = read.tree(text=treestrEnsif)
tipsEnsif = treeTips(resTribe, suborder="Ensifera")
nodesEnsif = treeNodes(resTribe, suborder="Ensifera")
nnodesEnsif = length(nodesEnsif)
ntipsEnsif = length(tipsEnsif) # 37
treestrCaelif = makeTree(resTribe, suborder="Caelifera")
trCaelif = read.tree(text=treestrCaelif)
tipsCaelif = treeTips(resTribe, suborder="Caelifera")
nodesCaelif = treeNodes(resTribe, suborder="Caelifera")
nnodesCaelif = length(nodesCaelif)
ntipsCaelif = length(tipsCaelif) #39

# Function for plotting the tree
plotTree = function(tr, resPie, resPiePerc, nodes, tips, ntips, label="") {
	plot(tr, node.depth=2, x.lim=8.1) # , no.margin=TRUE
	tiplabels(pie=resPiePerc[1:ntips,], cex=1.3, piecol=c("darkolivegreen3", "sandybrown", "purple", "white"))
	indices = rep(NA, length(nodes))
	for(i in 1:length(nodes)) indices[i] = which(resPie$Name==nodes[i])
	nodelabels(pie=resPiePerc[indices,], cex=1.3, piecol=c("darkolivegreen3", "sandybrown", "purple", "white"))
	tiplabels(tips, bg=NULL, adj=0, frame="none", cex=0.7, offset=0.4)
	#nodelabels(c(1:ntipsEnsif,1:nnodesEnsif), 1:66)
	text(-0.2, ntips, label, cex=1.2, pos=4)
}

# Figure 2
if(outToFile) png(file=paste0(path, outfolder, "Figure2.png"), width = 7, height = 7, units = "in", res=300) 
if(!outToFile) windows(7,7)
par(mfrow=c(1,2), mar=c(0,0,0,0))
plotTree(trEnsif, resPieEnsif, resPieEnsifPerc, nodesEnsif, tipsEnsif, ntipsEnsif, "Ensifera")
plotTree(trCaelif, resPieCaelif, resPieCaelifPerc, nodesCaelif, tipsCaelif, ntipsCaelif, "Caelifera")
if(outToFile)  dev.off()


##################################################################################
# Polymorphisms in taxa with green, brown and both morphs (Figure 3)
##################################################################################

rm(list=ls())

# Settings
path = "MYPATH\\OrthopteraPolymorphism\\"
version = "200326"
outfolder = "Figures\\"
filename = "Data\\OrthopteraPolymorphism"
cols = c("sandybrown", "darkolivegreen3", "purple")
outToFile = TRUE

# Load data
md = read.table(file = paste0(path, filename, "_Data_", version, ".txt"), header=TRUE)

# Function for aggregation at different levels
aggFunct = function(md, suborder="all", level="Tribe", select="all", minSpec=1) {
	if(suborder!="all" & all(select != "all")) mdsub = subset(md, !is.na(md$Polymorphic) & select & md$Suborder==suborder)
	if(suborder!="all" & all(select == "all")) mdsub = subset(md, !is.na(md$Polymorphic)          & md$Suborder==suborder)
	if(suborder=="all" & all(select != "all")) mdsub = subset(md, !is.na(md$Polymorphic) & select                        )
	mdagg = table(list(factor(mdsub[,level]), mdsub$Polymorphic))
	nSpecies = rowSums(mdagg)
	mdagg = mdagg[nSpecies>=minSpec,]
	mdagg = ceiling(mdagg/1000)
	taxonsum = data.frame(Aggregate=suborder, TaxonLevel=level, Taxon=rownames(mdagg), Morphs=paste(ifelse(mdagg[,1]>0,"brown",""),ifelse(mdagg[,2]>0,"green",""),ifelse(mdagg[,3]>0,"polymorphic",""), sep="-"), SingleTypeYN=as.numeric(rowSums(mdagg)==1))
	rownames(taxonsum) = NULL
	return(taxonsum)
}

# Function for calculation of suborder sums
subordsumFunct = function(taxonsum) {
	subordersum = table(taxonsum$Morphs, taxonsum$Aggregate)
	Morphs=factor(rownames(subordersum), levels=c("-green-", "brown--", "brown-green-", "brown--polymorphic", "-green-polymorphic", "brown-green-polymorphic", "--polymorphic"))
	subordersum = data.frame(Morphs=Morphs, Ensifera=subordersum[,1], Caelifera=subordersum[,2], All=subordersum[,3], row.names=NULL)
}

# Function
taxonClassFnc = function(minSpec = 1, select, filterByTribe=TRUE) {
	tribesumEnsif = aggFunct(md, "Ensifera", level="Tribe", minSpec=minSpec, select=select)
	tribesumCaelif = aggFunct(md, "Caelifera", level="Tribe", minSpec=minSpec, select=select)
	tribesumAll = aggFunct(md, "all", level="Tribe", minSpec=minSpec, select=select)
	subordersumTribes = subordsumFunct(rbind(tribesumEnsif, tribesumCaelif, tribesumAll))
	if(filterByTribe) selectEnsif = select & as.character(md$Tribe) %in% tribesumEnsif$Taxon[which(tribesumEnsif$SingleTypeYN==0)]
	if(filterByTribe) selectCaelif = select & as.character(md$Tribe) %in% tribesumCaelif$Taxon[which(tribesumCaelif$SingleTypeYN==0)]
	if(filterByTribe) selectAll = select & as.character(md$Tribe) %in% tribesumAll$Taxon[which(tribesumAll$SingleTypeYN==0)]
	if(!filterByTribe) selectEnsif = select
	if(!filterByTribe) selectCaelif = select
	if(!filterByTribe) selectAll = select
	generasumEnsif = aggFunct(md, "Ensifera", level="Genus", select=selectEnsif, minSpec=minSpec)
	generasumCaelif = aggFunct(md, "Caelifera", level="Genus", select=selectCaelif, minSpec=minSpec)
	generasumAll = aggFunct(md, "all", level="Genus", select=selectAll, minSpec=minSpec)
	subordersumGenera = subordsumFunct(rbind(generasumEnsif, generasumCaelif, generasumAll))
	res = merge(subordersumTribes, subordersumGenera, by="Morphs", all=TRUE)
	for(i in 2:7) res[is.na(res[,i]),i] = 0
	colnames(res) = c("Morphs", "EnsiferaTribes", "CaeliferaTribes", "AllTribes", "EnsiferaGenera", "CaeliferaGenera", "AllGenera")
	res$EnsiferaTribesProp = round(res$EnsiferaTribes/sum(res$EnsiferaTribes),2)
	res$EnsiferaGeneraProp = round(res$EnsiferaGenera/sum(res$EnsiferaGenera),2)
	res$CaeliferaTribesProp = round(res$CaeliferaTribes/sum(res$CaeliferaTribes),2)
	res$CaeliferaGeneraProp = round(res$CaeliferaGenera/sum(res$CaeliferaGenera),2)
	res$AllTribesProp = round(res$AllTribes/sum(res$AllTribes),2)
	res$AllGeneraProp = round(res$AllGenera/sum(res$AllGenera),2)
	res$Morphs = factor(res$Morphs, levels=c(levels(res$Morphs), "---")) 
	res = rbind(res, c("---", rep(0,8)))
	for(i in 2:13) res[,i] = as.numeric(res[,i])
	return(res)
}

# Function for plotting
barplotFnc = function(dat1, dat2, group,  ylim, rows=1:2, type=1, labels="", col=c("darkolivegreen3","purple"), col2=NULL, ylab="Number of groups") {
	dat1 = as.numeric(dat1[rows,group])
	dat2 = as.numeric(dat2[rows,group])
	rgb1 = as.vector(col2rgb(col[1]))/255
	rgb2 = as.vector(col2rgb(col[2]))/255
	barplot(dat1, col=c(rgb(rgb1[1], rgb1[2], rgb1[3], alpha=0.2), rgb(rgb2[1], rgb2[2], rgb2[3], alpha=0.2)), ylim=ylim, las=2)
	#barplot(dat1, col="white", ylim=ylim, las=2)
	if(type==1) title(ylab=ylab, cex.lab=1.6, font.lab=2)
	barplot(dat2, col=col, ylim=ylim, add=TRUE, las=2)
	if(type==3) barplot(dat2, col=col2, ylim=ylim, las=2, density=40, add=TRUE)
	axis(1, at=c(0.7,1.9), labels=labels, tick=FALSE)
}

# Data preparation
res1 = taxonClassFnc(minSpec=1, select = md$CricketsYN==0, filterByTribe=FALSE)
res2 = taxonClassFnc(minSpec=2, select = md$CricketsYN==0, filterByTribe=FALSE)
colGreenPoly = which(res1$Morphs=="-green-" | res1$Morphs=="-green-polymorphic")
colBrownPoly = which(res1$Morphs=="brown--" | res1$Morphs=="brown--polymorphic")
colBrGrPoly = which(res1$Morphs=="brown-green-" | res1$Morphs=="brown-green-polymorphic")
colPoly = which(res1$Morphs=="--polymorphic"  | res1$Morphs=="---")
nmaxtriEnsif=max(res1$EnsiferaTribes)
nmaxtriCaelif=max(res1$CaeliferaTribes)
nmaxtriAll=max(res1$AllTribes)
nmaxgenEnsif=max(res1$EnsiferaGenera)
nmaxgenCaelif=max(res1$CaeliferaGenera)
nmaxgenAll=max(res1$AllGenera)

# Figure 3
if(outToFile) png(file=paste0(path, outfolder, "Figure3.png"), width = 7, height = 5, units = "in", res=300) 
if(!outToFile) windows(7,5)
par(mfrow=c(2,4), mar=c(3,5,5,0.5))
ylim = c(0,nmaxtriAll)
group="AllTribes"
barplotFnc(res1, res2, group, ylim, rows=colGreenPoly, labels=c("non-poly", "poly"), col=c("darkolivegreen3","purple"), type=1, ylab="Number of tribes    ")
	title("Green and\nno brown", cex.main=1.6, line=1.8)
barplotFnc(res1, res2, group, ylim, rows=colBrownPoly, labels=c("non-poly", "poly"), col=c("sandybrown","purple"), type=2)
	title("Brown and\nno green", cex.main=1.6, line=1.8)
barplotFnc(res1, res2, group, ylim, rows=colBrGrPoly, labels=c("non-poly", "poly"), col=c("darkolivegreen3","purple"), col2=c("sandybrown","purple"), type=3)
	title("Green and\nbrown", cex.main=1.6, line=1.8)
barplotFnc(res1, res2, group, ylim, rows=rev(colPoly), labels=c("", "poly"), col=c("purple","purple"), type=4)
	title("Only\npolymorphic", cex.main=1.6, line=1.8)
ylim = c(0,nmaxgenAll)
par(mar=c(8,5,1,0.5))
group="AllGenera"
barplotFnc(res1, res2, group, ylim, rows=colGreenPoly, labels=c("non-poly", "poly"), col=c("darkolivegreen3","purple"), type=1, ylab="Number of genera  ")
barplotFnc(res1, res2, group, ylim, rows=colBrownPoly, labels=c("non-poly", "poly"), col=c("sandybrown","purple"), type=2)
barplotFnc(res1, res2, group, ylim, rows=colBrGrPoly, labels=c("non-poly", "poly"), col=c("darkolivegreen3","purple"), col2=c("sandybrown","purple"), type=3)
barplotFnc(res1, res2, group, ylim, rows=rev(colPoly), labels=c("", "poly"), col=c("purple","purple"), type=4)
ylim = c(0,nmaxgenCaelif)
if(outToFile) dev.off()

# Significance tests for differences between group "Green and no brown" and "Brown and no green" and "Green and brown"
# Tribes
mat = matrix(c(res2[res2$Morphs=="-green-","AllTribes"], res2[res2$Morphs=="-green-polymorphic", "AllTribes"],
				res2[res2$Morphs=="brown--", "AllTribes"], res2[res2$Morphs=="brown--polymorphic", "AllTribes"],
				res2[res2$Morphs=="brown-green-", "AllTribes"], res2[res2$Morphs=="brown-green-polymorphic", "AllTribes"]), ncol=3)
rownames(mat) = c("monomorphic", "polymorphic")
colnames(mat) = c("Green", "Brown", "Both")
chisq.test(mat)
# Genera
mat = matrix(c(res2[res2$Morphs=="-green-","AllGenera"], res2[res2$Morphs=="-green-polymorphic", "AllGenera"],
				res2[res2$Morphs=="brown--", "AllGenera"], res2[res2$Morphs=="brown--polymorphic", "AllGenera"],
				res2[res2$Morphs=="brown-green-", "AllTribes"], res2[res2$Morphs=="brown-green-polymorphic", "AllGenera"]), ncol=3)
rownames(mat) = c("monomorphic", "polymorphic")
colnames(mat) = c("Green", "Brown", "Both")
chisq.test(mat)


##################################################################################
# Rareness (Figure 4)
##################################################################################

rm(list=ls())

# Settings
path = "MYPATH\\OrthopteraPolymorphism\\"
version = "200326"
outfolder = "Figures\\"
filename = "Data\\OrthopteraPolymorphism"
cols = c("sandybrown", "darkolivegreen3", "purple")
outToFile = TRUE

# Load data
md = read.table(file = paste0(path, filename, "_Data_", version, ".txt"), header=TRUE)

# Data preparation
md$Abundance = md$Brown
md$Abundance[md$Green == "absent" ] = "absent"
md$Abundance[md$Green == "exceptional"] = "exceptional"
md$Abundance[md$Green == "rare"] = "rare"
md$Abundance[md$Green == "regular"] = "regular"
md$Abundance[md$Green == "common"] = "common"
md$Abundance = factor(md$Abundance, levels=c("absent", "exceptional", "rare", "regular", "common"))

# Data preparation
md$Combined[md$Combined=="NA-NA"] = NA
md$Combined = factor(md$Combined,
			  levels=c("absent-mono", "exceptional-dominant", "rare-dominant", 
						"regular-dominant", "common-dominant", 
						"dominant-common", "dominant-regular", "dominant-rare", "dominant-exceptional", 
						"mono-absent"))
levels(md$Combined)[which(levels(md$Combined)=="absent-mono")] = "brown"
levels(md$Combined)[which(levels(md$Combined)=="mono-absent")] = "green"

# Figure 4
if(outToFile) png(file=paste0(path, outfolder, "Figure4.png"), width = 7, height = 7, units = "in", res=300) 
if(!outToFile) windows(7,7)
par(mfrow=c(1,1), mar=c(10,4,1,1))
barplot(table(md$Combined[md$Combined != "brown" & md$Combined != "green"]),las=2, ylab="Number of species", col="darkolivegreen3", xaxt="n")
barplot(table(md$Combined[md$Combined != "brown" & md$Combined != "green"]),las=2, ylab="Number of species", col="sandybrown", 
	density=seq(100, 0, length.out=10), xaxt="n", add=TRUE)
axis(1, at=seq(0.7,11.5, length.out=10), labels=c("Monomorphic brown", "Green exceptional", "Green rare", "Green regular", "Green common",
	"Brown common", "Brown regular", "Brown rare", "Brown exceptional", "Monomorphic green"), las=2)
barplot(table(md$Combined[md$Combined == "brown" | md$Combined == "green"]),las=2, ylab="", col=c("sandybrown", "darkolivegreen3"), xaxt="n", add=TRUE)
text(0.75, 2, paste0(sum(md$Combined == "brown", na.rm=TRUE), " species"), adj=0, srt=90)
text(11.5, 2, paste0(sum(md$Combined == "green", na.rm=TRUE), " species"), adj=0, srt=90)
if(outToFile) dev.off()


##################################################################################
# Distribution (Figure 5)
##################################################################################

rm(list=ls())

# Settings
path = "MYPATH\\OrthopteraPolymorphism\\"
version = "200326"
outfolder = "Figures\\"
filename = "Data\\OrthopteraPolymorphism"
cols = c("sandybrown", "darkolivegreen3", "purple")
outToFile = TRUE

# Load data
md = read.table(file = paste0(path, filename, "_Data_", version, ".txt"), header=TRUE)

# Find columns that contain data on distribution
distcolumns = which(substr(colnames(md),1,3) %in% c(paste0("R0", 1:9), paste0("R", 10:13)))

# Function for aggregation at different taxonimic levels
aggDistribution = function(md, level="Species", columns=18:30) {
	res=data.frame(Distribution=colnames(md)[columns], Brown=NA, Green=NA, Polymorphic=NA, nTaxa=NA, Level=level)
	for(i in 1:length(columns)) {
		mdsub = md[md[,columns[1]-1+i]=="x" & !is.na(md$Combined),]
		mdagg = data.frame(table(list(mdsub$Polymorphic, mdsub[,level])))
		colnames(mdagg)  = c("Combined", "Taxon", "Count")
		res[i,2:4] = tapply(I(mdagg$Count>0),mdagg$Combined,sum)
		res$nTaxa[i] = sum(table(mdsub[,level])>0)
	}
	res$BrownProp = res$Brown/res$nTaxa #rowSums(res[,2:4])
	res$GreenProp = res$Green/res$nTaxa #rowSums(res[,2:4])
	res$PolymorphicProp = res$Polymorphic/res$nTaxa #rowSums(res[,2:4])
	res = cbind(res, BrownAvg = NA, GreenAvg = NA, PolymorphicAvg = NA, BrownP = NA, GreenP = NA, PolymorphicP = NA)
	mdsub = md[!is.na(md$Combined),]
	mdagg = data.frame(table(list(mdsub$Polymorphic, mdsub[,level])))
	colnames(mdagg)  = c("Combined", "Taxon", "Count")
	avg = tapply(I(mdagg$Count>0),mdagg$Combined,sum) / sum(table(mdsub[,level])>0)
	res$BrownAvg = avg[1]
	res$GreenAvg = avg[2]
	res$PolymorphicAvg = avg[3]
	for(i in 1:nrow(res)) {
		res$BrownP[i] = chisq.test(c(res$Brown[i], res$nTaxa[i]-res$Brown[i]), p=c(res$BrownAvg[i], 1-res$BrownAvg[i]))$p.value
		res$GreenP[i] = chisq.test(c(res$Green[i], res$nTaxa[i]-res$Green[i]), p=c(res$GreenAvg[i], 1-res$GreenAvg[i]))$p.value
		res$PolymorphicP[i] = chisq.test(c(res$Polymorphic[i], res$nTaxa[i]-res$Polymorphic[i]), p=c(res$PolymorphicAvg[i], 1-res$PolymorphicAvg[i]))$p.value
	}
	return(res)
}

# Function for plotting
plotDistribution = function(res, regions, main="", ymax=1, add=FALSE, alpha=1, ylab="Proportion", col=c("sandybrown", "darkolivegreen3", "purple")) {
	rgb1 = as.vector(col2rgb(col[1]))/255
	rgb2 = as.vector(col2rgb(col[2]))/255
	rgb3 = as.vector(col2rgb(col[3]))/255
	col1 = rgb(rgb1[1], rgb1[2], rgb1[3], alpha=alpha)
	col2 = rgb(rgb2[1], rgb2[2], rgb2[3], alpha=alpha)
	col3 = rgb(rgb3[1], rgb3[2], rgb3[3], alpha=alpha)	
	res$Distribution = factor(res$Distribution, levels=regions) 
	levels(res$Distribution) = substr(levels(res$Distribution),5,100)
	if(!add) plot(res$Distribution, rep(-2,13), ylim=c(0,ymax), las=2, ylab=ylab, type="n", xaxt="n", xlab="", main=main)
	axis(1, at=1:13, labels = gsub("SwedenNorway", "Sweden/Norway", gsub("FinlandBaltic", "Finland/Baltic", gsub("BalkanBulgariaRomania", "Balkan/\nBulgaria/Romania", gsub("_", " ", as.character(res$Distribution)), fixed=TRUE))), las=2)
	pointsFnc = function(res, cols) {
		points(res$Distribution[cols], res$BrownProp[cols], type="b", col=col1, lwd=2, cex=c(2,1)[add+1], pch=c(19,1)[add+1])	
		points(res$Distribution[cols], res$GreenProp[cols], type="b", col=col2, lwd=2, cex=c(2,1)[add+1], pch=c(19,1)[add+1])	
		points(res$Distribution[cols], res$PolymorphicProp[cols], type="b", col=col3, lwd=2, cex=c(2,1)[add+1], pch=c(19,1)[add+1])	
	}
	pointsFnc(res, cols=c(1:3))
	pointsFnc(res, cols=c(4:6))
	pointsFnc(res, cols=c(7:13))
	segments(x0=c(1,4,7), x1=c(3,6,13), y0=ymax*0.95, y1=ymax*0.95)
	text(x=c(2,5,10), y=ymax, labels=c("North", "Central", "South"))
	abline(v=c(0.5,3.5, 6.5,13.5))
}

# Figure 5
if(outToFile) png(file=paste0(path, outfolder, "Figure5.png"), width = 7, height = 7, units = "in", res=300) 
if(!outToFile) windows(7,7)
par(mfrow=c(1,1), mar=c(10,4,1,1))
regions = colnames(md)[distcolumns]
res = aggDistribution(md[md$CricketsYN==0,], "Species", columns=distcolumns)
plotDistribution(res, regions, main="", ymax=1.15, alpha=1, ylab="Proportion of species (without crickets)", col=cols)
if(outToFile) dev.off() 


##################################################################################
# Habitats (Figure 6)
##################################################################################

rm(list=ls())

# Settings
path = "MYPATH\\OrthopteraPolymorphism\\"
version = "200326"
outfolder = "Figures\\"
filename = "Data\\OrthopteraPolymorphism"
cols = c("sandybrown", "darkolivegreen3", "purple")
outToFile = TRUE

# Load data
md = read.table(file = paste0(path, filename, "_Data_", version, ".txt"), header=TRUE)

# Find columns that contain data on distribution
habitatcolumns = which(substr(colnames(md),1,3) %in% c(paste0("H0", 1:9), paste0("H", 10:21)))

# Function for aggregation at different taxonimic levels
aggHabitat = function(md, level="Species", columns=33:53) {
	res=data.frame(Habitat=colnames(md)[columns], Brown=NA, Green=NA, Polymorphic=NA, nTaxa=NA, Level=level)
	for(i in 1:length(columns)) {
		mdsub = md[md[,columns[1]-1+i]=="x" & !is.na(md$Combined),]
		mdagg = data.frame(table(list(mdsub$Polymorphic, mdsub[,level])))
		colnames(mdagg)  = c("Combined", "Taxon", "Count")
		res[i,2:4] = tapply(I(mdagg$Count>0),mdagg$Combined,sum)
		res$nTaxa[i] = sum(table(mdsub[,level])>0)
	}
	res$BrownProp = res$Brown/res$nTaxa #rowSums(res[,2:4])
	res$GreenProp = res$Green/res$nTaxa #rowSums(res[,2:4])
	res$PolymorphicProp = res$Polymorphic/res$nTaxa #rowSums(res[,2:4])
	mdsub = md[!is.na(md$Combined),]
	mdagg = data.frame(table(list(mdsub$Polymorphic, mdsub[,level])))
	colnames(mdagg)  = c("Combined", "Taxon", "Count")
	return(res)
}

# Function for plotting
plotHabitat = function(res, main="", ymax=1, add=FALSE, alpha=1, order=FALSE, columns=30:50, ylab="Proportion", col=c("sandybrown", "darkolivegreen3", "purple")) {
	rgb1 = as.vector(col2rgb(col[1]))/255
	rgb2 = as.vector(col2rgb(col[2]))/255
	rgb3 = as.vector(col2rgb(col[3]))/255
	col1 = rgb(rgb1[1], rgb1[2], rgb1[3], alpha=alpha)
	col2 = rgb(rgb2[1], rgb2[2], rgb2[3], alpha=alpha)
	col3 = rgb(rgb3[1], rgb3[2], rgb3[3], alpha=alpha)	
	#res$Habitat = factor(res$Habitat, levels=colnames(md)[30:50])
	if(order) index = sort(res$PolymorphicProp, index.return=TRUE)$ix
	if(!order) index = 1:nrow(res)
	lev = colnames(md)[columns][index]
	lev = c(lev[1:7], "_", lev[8:12], "__", lev[13:15], "___", lev[16:21])
	res$Habitat = factor(res$Habitat, levels=lev)
	index = 1:length(lev)
	if(!add) plot(res$Habitat[index], rep(-2,length(lev)), ylim=c(0,ymax), las=2, ylab=ylab, xaxt="n", xlab="", type="n", main=main)
	axis(1, at=1:length(lev), labels = gsub("_", " ", substr(as.character(lev),5,100)), las=2)
	pointsFnc = function(res, index, cols) {
		points(res$Habitat[index][cols], res$BrownProp[index][cols], type="b", col=col1, lwd=3, cex=c(2,2)[as.numeric(res$BrownP<0.5)+1], pch=c(19,19)[as.numeric(res$BrownP<0.5)+1])	
		points(res$Habitat[index][cols], res$GreenProp[index][cols], type="b", col=col2, lwd=3, cex=c(2,2)[as.numeric(res$GreenP<0.5)+1], pch=c(19,19)[as.numeric(res$GreenP<0.5)+1])	
		points(res$Habitat[index][cols], res$PolymorphicProp[index][cols], type="b", col=col3, lwd=3, cex=c(2,2)[as.numeric(res$PolymorphicP<0.5)+1], pch=c(19,19)[as.numeric(res$PolymorphicP<0.5)+1])	
	}
	pointsFnc(res, index, cols=c(1:7))
	pointsFnc(res, index, cols=c(8:12))
	pointsFnc(res, index, cols=c(13:15))
	pointsFnc(res, index, cols=c(16:21))
	segments(x0=c(0.5,8.5,14.5,18.5), x1=c(7.5,13.5,17.5,24.5), y0=ymax*0.95, y1=ymax*0.95)
	text(x=c(4,11,16,21.5), y=ymax, labels=c("Grasslands", "Bushland", "Alpine", "Other habitats"))
	#abline(v=c(0.5,2.5,9.5,14.5,17.5,21.5))
	abline(v=c(0,8,14,18,25))
}

# Figure 6
if(outToFile) png(file=paste0(path, outfolder, "Figure6.png"), width = 7, height = 7, units = "in", res=300) 
if(!outToFile) windows(7,7)
par(mfrow=c(1,1), mar=c(10,4,1,1))
res = aggHabitat(md, level="Species", columns=habitatcolumns)
plotHabitat(res, main="", ymax=1.15, add=FALSE, alpha=1, order=FALSE, ylab="Proportion of species", columns=habitatcolumns, col=cols)
if(outToFile) dev.off() 


##################################################################################
# Habitat use (Figure 7)
##################################################################################

rm(list=ls())

# Settings
path = "MYPATH\\OrthopteraPolymorphism\\"
version = "200326"
outfolder = "Figures\\"
filename = "Data\\OrthopteraPolymorphism"
cols = c("sandybrown", "darkolivegreen3", "purple")
outToFile = TRUE

# Load data
md = read.table(file = paste0(path, filename, "_Data_", version, ".txt"), header=TRUE)

# Find columns that contain data on distribution
habitatusecolumns = which(substr(colnames(md),1,3) %in% c(paste0("U0", 1:9), paste0("U", 10:11)))

# Function for aggregation at different taxonimic levels
aggHabitatUse = function(md, level="Species", columns=52:63) {
	res=data.frame(HabitatUse=colnames(md)[columns], Brown=NA, Green=NA, Polymorphic=NA, nTaxa=NA, Level=level)
	for(i in 1:length(columns)) {
		mdsub = md[md[,columns[1]-1+i]=="x" & !is.na(md$Combined),]
		mdagg = data.frame(table(list(mdsub$Polymorphic, mdsub[,level])))
		colnames(mdagg)  = c("Combined", "Taxon", "Count")
		res[i,2:4] = tapply(I(mdagg$Count>0),mdagg$Combined,sum)
		res$nTaxa[i] = sum(table(mdsub[,level])>0)
	}
	res$BrownProp = res$Brown/res$nTaxa #rowSums(res[,2:4])
	res$GreenProp = res$Green/res$nTaxa #rowSums(res[,2:4])
	res$PolymorphicProp = res$Polymorphic/res$nTaxa #rowSums(res[,2:4])
	res$HabitatUse = factor(res$HabitatUse, levels=colnames(md)[columns])
	mdsub = md[!is.na(md$Combined),]
	mdagg = data.frame(table(list(mdsub$Polymorphic, mdsub[,level])))
	colnames(mdagg)  = c("Combined", "Taxon", "Count")
	avg = tapply(I(mdagg$Count>0),mdagg$Combined,sum) / sum(table(mdsub[,level])>0)
	return(res)
}

# Function for plotting
plotHabitatuse = function(res, main="", ymax=1, add=FALSE, alpha=1, ylab="Proportion", col=c("sandybrown", "darkolivegreen3", "purple")) {
	rgb1 = as.vector(col2rgb(col[1]))/255
	rgb2 = as.vector(col2rgb(col[2]))/255
	rgb3 = as.vector(col2rgb(col[3]))/255
	col1 = rgb(rgb1[1], rgb1[2], rgb1[3], alpha=alpha)
	col2 = rgb(rgb2[1], rgb2[2], rgb2[3], alpha=alpha)
	col3 = rgb(rgb3[1], rgb3[2], rgb3[3], alpha=alpha)	
	#res$Habitat = factor(res$Habitat, levels=colnames(md)[30:50])
	if(!add) plot(res$HabitatUse, rep(-2,11), ylim=c(0,ymax), las=2, ylab=ylab, type="n", xaxt="n", xlab="", main=main)
	axis(1, at=1:11, labels = gsub(".", "\n", gsub("_", " ", substr(as.character(res$HabitatUse),5,100)), fixed=TRUE), las=2)
	pointsFnc = function(res, cols) {
		points(res$HabitatUse[cols], res$BrownProp[cols], type="b", col=col1, lwd=3, cex=c(2,2)[as.numeric(res$BrownP<0.5)+1], pch=c(19,19)[as.numeric(res$BrownP<0.5)+1])	
		points(res$HabitatUse[cols], res$GreenProp[cols], type="b", col=col2, lwd=3, cex=c(2,2)[as.numeric(res$GreenP<0.5)+1], pch=c(19,19)[as.numeric(res$GreenP<0.5)+1])	
		points(res$HabitatUse[cols], res$PolymorphicProp[cols], type="b", col=col3, lwd=3, cex=c(2,2)[as.numeric(res$PolymorphicP<0.5)+1], pch=c(19,19)[as.numeric(res$PolymorphicP<0.5)+1])	
	}
	pointsFnc(res, cols=c(1:5))
	pointsFnc(res, cols=c(6:8))
	pointsFnc(res, cols=c(9:11))
	segments(x0=c(1,6,9), x1=c(5,8,11), y0=ymax*0.95, y1=ymax*0.95)
	text(x=c(3,7,10), y=ymax, labels=c("Habitat use", "Activity", "Seasonality"))
	abline(v=c(0.5,5.5,8.5,11.5))
}
	
# Figure 7
if(outToFile) png(file=paste0(path, outfolder, "Figure7.png"), width = 7, height = 7, units = "in", res=300) 
if(!outToFile) windows(7,7)
par(mfrow=c(1,1), mar=c(10,4,1,1))
res = aggHabitatUse(md, level="Species", columns=habitatusecolumns)
plotHabitatuse(res, main="", ymax=1.15, add=FALSE, alpha=1, ylab="Proportion of species", col=cols)
if(outToFile) dev.off() 


##################################################################################
# Species richness, range size and habitat specificity (Figure 8)
##################################################################################

rm(list=ls())

# Settings
path = "MYPATH\\OrthopteraPolymorphism\\"
version = "200326"
outfolder = "Figures\\"
filename = "Data\\OrthopteraPolymorphism"
cols = c("sandybrown", "darkolivegreen3", "purple")
outToFile = TRUE

# Load data
md = read.table(file = paste0(path, filename, "_Data_", version, ".txt"), header=TRUE)

# Find columns that contain data on distribution
habitatcolumns = which(substr(colnames(md),1,3) %in% c(paste0("H0", 1:9), paste0("H", 10:21)))
distcolumns = which(substr(colnames(md),1,3) %in% c(paste0("R0", 1:9), paste0("R", 10:13)))

# Function for aggregation at different taxonimic levels
aggFunct = function(md, level="Tribe", type="Polymorphic") {
	md$GreenOnly = ifelse(!is.na(md[,type]), ifelse(md[,type]=="green",1,0), 0)
	md$BrownOnly = ifelse(!is.na(md[,type]), ifelse(md[,type]=="brown",1,0), 0)
	md$PolymorphicOnly = ifelse(!is.na(md[,type]), ifelse(md[,type]=="polymorphic",1,0), 0)
	md$MissingOnly = ifelse(is.na(md[,type]), 1,0)
	if(level=="Species") cols = 1:7
	if(level=="Genus") cols = 1:6
	if(level=="Tribe") cols = 1:5
	if(level=="Subfamily") cols = 1:4
	if(level=="Family") cols = 1:3
	if(level=="Superfamily") cols = 1:2
	if(level=="Suborder") cols = 1:1
	taxonlevels = c("Suborder", "Superfamily", "Family", "Subfamily", "Tribe", "Genus", "Species")
	if(level=="Species") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$MissingOnly), by=list(md$Suborder, md$Superfamily, md$Family, md$Subfamily, md$Tribe, md$Genus, md$Species), FUN=sum)
	if(level=="Genus") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$MissingOnly), by=list(md$Suborder, md$Superfamily, md$Family, md$Subfamily, md$Tribe, md$Genus), FUN=sum)
	if(level=="Tribe") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$MissingOnly), by=list(md$Suborder, md$Superfamily, md$Family, md$Subfamily, md$Tribe), FUN=sum)
	if(level=="Subfamily") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$MissingOnly), by=list(md$Suborder, md$Superfamily, md$Family, md$Subfamily), FUN=sum)
	if(level=="Family") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$MissingOnly), by=list(md$Suborder, md$Superfamily, md$Family), FUN=sum)
	if(level=="Superfamily") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$MissingOnly), by=list(md$Suborder, md$Superfamily), FUN=sum)
	if(level=="Suborder") res = aggregate(list(md$GreenOnly, md$BrownOnly, md$PolymorphicOnly, md$MissingOnly), by=list(md$Suborder), FUN=sum)
	colnames(res) = c(taxonlevels[cols], "Green", "Brown", "Polymorphic", "Missing")
	res$nSpecies = rowSums(res[,c("Green", "Brown", "Polymorphic", "Missing")])
	res$PolymorphicPerc = res$Polymorphic/rowSums(res[,c("Green", "Brown", "Polymorphic")])
	res$Name = res[,level]
	return(res)
}

# Habitat use breadth
md$HabitatUseBreadth = NA
for(i in 1:nrow(md))
	md$HabitatUseBreadth[i] = ifelse(is.na(md[i,habitatcolumns[1]]) | md[i,habitatcolumns[1]]==".", NA, sum(md[i,habitatcolumns]=="x"))
md$DistributionBreadth = NA
for(i in 1:nrow(md))
	md$DistributionBreadth[i] = ifelse(is.na(md[i,distcolumns[1]]),NA, sum(md[i,distcolumns]=="x"))
res = aggFunct(md, level="Genus", type="Polymorphic")

# Figure 8
if(outToFile) png(file=paste0(path, outfolder, "Figure8.png"), width = 7, height = 7, units = "in", res=300) 
if(!outToFile) windows(7,7)
par(mfrow=c(2,2), mar=c(5,4,3,1))
# Subplot a
Suborder="Ensifera"
jitx = rnorm(length(res$PolymorphicPerc[res$Suborder==Suborder]),0,0.01)
jity = rnorm(length(res$PolymorphicPerc[res$Suborder==Suborder]),0,0.1)
plot(res$PolymorphicPerc[res$Suborder==Suborder]+jitx, log10(res$nSpecies[res$Suborder==Suborder]+jity), yaxt="n", xaxt="n", col=rgb(0.5,0.5,0.5,alpha=0.4), ylim=log10(c(0.9,68)), pch=19, cex=1.2, type="n", ylab="Number of species", xlab="Proportion polymorphic species")
	title(main="(a) Species diversity: Ensifera")
	axis(1, at=c(0,0.25,0.5,0.75,1), labels=as.character(c(0,0.25,0.5,0.75,1)))
	axis(2, at=log10(c(1,2,4,8,16,32,64)), labels=c(1,2,4,8,16,32,64), las=2)
	abline(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>0]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>0]), col="grey20")
	abline(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>1]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>1]), col="grey30")
	abline(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>2]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>2]), col="grey50")
	abline(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>4]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>4]), col="grey60")
	abline(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>8]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>8]), col="grey70")
	points(res$PolymorphicPerc[res$Suborder==Suborder]+jitx, log10(res$nSpecies[res$Suborder==Suborder]+jity), col=rgb(0.5,0.5,0.5,alpha=0.4), pch=19, cex=1.2)
	# Pvalues for regression analyses
	coef(summary(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>0]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>0])))[2,"Pr(>|t|)"]
	coef(summary(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>1]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>1])))[2,"Pr(>|t|)"]
	coef(summary(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>2]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>2])))[2,"Pr(>|t|)"]
	coef(summary(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>4]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>4])))[2,"Pr(>|t|)"]
	coef(summary(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>8]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>8])))[2,"Pr(>|t|)"]
# Subplot b
Suborder="Caelifera"
jitx = rnorm(length(res$PolymorphicPerc[res$Suborder==Suborder]),0,0.01)
jity = rnorm(length(res$PolymorphicPerc[res$Suborder==Suborder]),0,0.1)
plot(res$PolymorphicPerc[res$Suborder==Suborder]+jitx, log10(res$nSpecies[res$Suborder==Suborder]+jity), xaxt="n", yaxt="n", col=rgb(0.5,0.5,0.5,alpha=0.4), ylim=log10(c(0.9,68)), pch=19, cex=1.2, type="n", ylab="Number of species", xlab="Proportion polymorphic species")
	title(main="(b) Species diversity: Caelifera")
	axis(1, at=c(0,0.25,0.5,0.75,1), labels=as.character(c(0,0.25,0.5,0.75,1)))
	axis(2, at=log10(c(1,2,4,8,16,32,64)), labels=c(1,2,4,8,16,32,64), las=2)
	abline(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>0]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>0]), col="grey20")
	abline(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>1]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>1]), col="grey30")
	abline(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>2]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>2]), col="grey50")
	abline(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>4]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>4]), col="grey60")
	abline(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>8]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>8]), col="grey70")
	points(res$PolymorphicPerc[res$Suborder==Suborder]+jitx, log10(res$nSpecies[res$Suborder==Suborder]+jity), col=rgb(0.5,0.5,0.5,alpha=0.4), pch=19, cex=1.2)
	coef(summary(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>0]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>0])))[2,"Pr(>|t|)"]
	coef(summary(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>1]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>1])))[2,"Pr(>|t|)"]
	coef(summary(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>2]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>2])))[2,"Pr(>|t|)"]
	coef(summary(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>4]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>4])))[2,"Pr(>|t|)"]
	coef(summary(lm(log10(res$nSpecies[res$Suborder==Suborder & res$nSpecies>8]) ~ res$PolymorphicPerc[res$Suborder==Suborder & res$nSpecies>8])))[2,"Pr(>|t|)"]
# Subplot c
res = data.frame(tapply(md$Polymorphic, list(md$DistributionBreadth, md$Polymorphic), length))
res$BrownProp = ifelse(is.na(res$brown), 0, res$brown) / rowSums(res[,1:3], na.rm=TRUE)
res$GreenProp = ifelse(is.na(res$green), 0, res$green) / rowSums(res[,1:3], na.rm=TRUE)
res$PolymorphicProp = ifelse(is.na(res$polymorphic), 0, res$polymorphic) / rowSums(res[,1:3], na.rm=TRUE)
plot(as.numeric(rownames(res)), res$PolymorphicProp, ylim=c(0,1), col="purple",  ylab="Proportion of species", xlab="Geographic regions occupied", type="b", lwd=3, cex=2, pch=19, las=1)
title(main="(c) Range size")
points(as.numeric(rownames(res)), res$BrownProp, col="sandybrown",  type="b", lwd=3, cex=2, pch=19)
points(as.numeric(rownames(res)), res$GreenProp, col="darkolivegreen3",  type="b", lwd=3, cex=2, pch=19)
abline(lm(res$PolymorphicProp ~ as.numeric(rownames(res))), col="purple")
coef(summary(lm(res$PolymorphicProp ~ as.numeric(rownames(res)))))[2,"Pr(>|t|)"]
# Subplot d
res = with(subset(md, !is.na(HabitatUseBreadth)), data.frame(tapply(Polymorphic, list(HabitatUseBreadth, Polymorphic), length)))
res$BrownProp = ifelse(is.na(res$brown), 0, res$brown) / rowSums(res[,1:3], na.rm=TRUE)
res$GreenProp = ifelse(is.na(res$green), 0, res$green) / rowSums(res[,1:3], na.rm=TRUE)
res$PolymorphicProp = ifelse(is.na(res$polymorphic), 0, res$polymorphic) / rowSums(res[,1:3], na.rm=TRUE)
plot(as.numeric(rownames(res)), res$PolymorphicProp, ylim=c(0,1), col="purple",  ylab="Proportion of species", xlab="Number of habitat categories", type="b", lwd=3, cex=2, pch=19, las=1)
title(main="(d) Habitat specificity")
points(as.numeric(rownames(res)), res$BrownProp, col="sandybrown",  type="b", lwd=3, cex=2, pch=19)
points(as.numeric(rownames(res)), res$GreenProp, col="darkolivegreen3",  type="b", lwd=3, cex=2, pch=19)
abline(lm(res$PolymorphicProp ~ as.numeric(rownames(res))), col="purple")
coef(summary(lm(res$PolymorphicProp ~ as.numeric(rownames(res)))))[2,"Pr(>|t|)"]
if(outToFile) dev.off()
