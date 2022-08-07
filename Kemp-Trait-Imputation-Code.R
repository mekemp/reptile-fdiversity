#Imputation

require(picante)
require(phytools)
require(TreeTools)
require(missForest)
require(PVR)
#require(Rphylopars)
phy<-read.tree("squamata.txt")
#plot(phy)


#write.csv(labels, "tiplabels282022.csv")
#go through and scan the file, edit tips if necessary
#edit necessary tip labels
phy$tip.label[2913]<-c("Mitophis_asbolepis")
phy$tip.label[2914]<-c("Mitophis_leptepileptus")
phy$tip.label[2915]<-c("Mitophis_pyrites")
phy$tip.label[1493]<-c("Pholidoscelis_auberi")
phy$tip.label[1497]<-c("Pholidoscelis_chrysolaemus")
phy$tip.label[1496]<-c("Pholidoscelis_exsul")
phy$tip.label[1506]<-c("Pholidoscelis_griswoldi")
phy$tip.label[1501]<-c("Pholidoscelis_lineolatus")
phy$tip.label[1503]<-c("Pholidoscelis_plei")
phy$tip.label[1499]<-c("Pholidoscelis_taeniurus")
phy$tip.label[1492]<-c("Pholidoscelis_dorsalis")
phy$tip.label[3920]<-c("Nerodia_clarkii")


labels<-cbind(TipLabels(phy))

spdb<-read.csv("CARIB-REP-TRAITS.csv")
species<-cbind(spdb$Binomial)

labels<-gsub("_", " ", labels)
phy$tip.label[which(labels%in%species)]

#match.phylo.data(phy,species)
drop<-subset(labels, !labels%in%species)
drop<-gsub(" ","_",drop)
drop
phy2<-drop.tip(phy, drop)
phy2
#plot(phy2)
length(species)
length(phy2$tip.label)
?add.species.to.genus

write.csv(drop, "dropped_tips.csv")
tips<-TipLabels(phy2)
rem<-subset(species, !species%in%gsub("_"," ", tips))
rem<-sort(rem)
length(rem)
#write.csv(rem,"species-to-add-to-tree-.csv")
#go through previous file and edit genus names as needed
spadd<-read.csv("species-to-add-to-tree-wgenus.csv")

phy3<-force.ultrametric(phy2, method = "extend")
spadd$species<-sub(" ","_", spadd$species)
length(phy3$tip.label)
spadd
#phy4<-add.species.to.genus(phy3,spadd$x,genus=spadd$genus,where="root")
added<-c()
for(i in 1:length(spadd[,1])){
  phy3<-add.species.to.genus(phy3,spadd$species[i], genus=spadd$genus[i],where="root")
  added<-append(added,spadd$species[i])
}
#write.csv(phy3$tip.label, "phy3tips02102022.csv")
length(phy3$tip.label)
subset(spadd$species, !spadd$species%in%phy3$tip.label)
phy3<-add.species.to.genus(phy3, "Clelia_cf._clelia", genus="Liophis", where="root")
phy3<-add.species.to.genus(phy3, "Typhlops_cf._silus", genus="Typhlops", where="root")
length(phy3$tip.label)
plot.phylo(phy, type="radial", cex=.5)

phylDiss <- sqrt(cophenetic(phy3))
hist(phylDiss)


hist(phylDiss)
pcoaPhyl <- cmdscale(phylDiss, k=10) 
colnames(pcoaPhyl) <- paste0("Eigen.", 1:ncol(pcoaPhyl))
rownames(pcoaPhyl)<-sub("_", " ", rownames(pcoaPhyl))
write.csv(pcoaPhyl, "PCOA-squam.csv")

spdb<-read.csv("CARIB-REP-TRAITS.csv")
spdb<-as.data.frame(spdb)
rownames(spdb)<-spdb[,1]

rownames(pcoaPhyl)<-gsub("_", " ", rownames(pcoaPhyl))
traitsInPhyl <- rownames(spdb)[which(rownames(spdb) %in% rownames(pcoaPhyl))]
phylInTraits <- rownames(pcoaPhyl)[which(rownames(pcoaPhyl) %in% rownames(spdb))]
commonSpecies <- unique(c(traitsInPhyl, phylInTraits))

pcoaPhylWithTraits <- pcoaPhyl[commonSpecies, ]
traitsWithPhyl <- spdb[commonSpecies, ]
traitsAndPCOA <- cbind(traitsWithPhyl, pcoaPhylWithTraits)
write.csv(traitsAndPCOA, "CARIB-REP-TRAITS-pcoa-squamates.csv")

reptileData<-read.csv("CARIB-REP-TRAITS-pcoa-squamates.csv")

rownames(reptileData)<-reptileData$Binomial
reptileDatacopy<-reptileData
reptileData<-reptileData[,c(-1)]
colnames(reptileData)
reptileData<-reptileData[,c(1:7,15:24)]
colnames(reptileData)
reptileData<-reptileData[,c(-1)]
reptileData$mass<-log10(reptileData$mass)
reptileData$activity<-as.factor(reptileData$activity)
reptileData$habitatex<-as.factor(reptileData$habitatex)
reptileData$diet<-as.factor(reptileData$diet)
reptileData$foraging<-as.factor(reptileData$foraging)
reptileData$reprod<-as.factor(reptileData$reprod)

#cutoffs<-list(1, rep(.05, 20), rep(.25, 4), rep(.5,.2), rep(1/3,3), rep(1/3,3), 1,1,1,1,1,1,1,1,1,1 )
cutoffs<-list(1, c(.75,.20, rep(.05/18,18)), c(.75,.20, rep(.05/2,2)), c(.75,.25), c(.75,.20,.05), c(.75,.2,.05), 1,1,1,1,1,1,1,1,1,1 )
wts<-list("NULL", rep(.05, 20), rep(.25, 4), rep(.5, 2), rep(1/3,3), rep(1/3,3), "NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL" )
reptileTraitslogImputed<-as.matrix(missForest(xmis=reptileData, classwt=wts,cutoff=cutoffs)$ximp)
write.csv(reptileTraitslogImputed, "Carib-Traits-imputed.csv")

  