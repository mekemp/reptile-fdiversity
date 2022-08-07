#set working directory

library(mFD)
library(dplyr)

###### island categorization and color code ######
islands<-c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3),rep(7,3),rep(8,3),rep(9,3),rep(10,3),rep(11,3),rep(12,3),rep(13,3),rep(14,3),rep(15,3),rep(16,3),rep(17,3))
isnames<-c("Anguilla",
           "Antigua",
           "Barbuda",
           "Cuba",
           "Culebra",
           "Grand Terre",
           "Grand Turk",
           "Great Abaco",
           "Hispaniola",
           "Jamaica",
           "La Desirade",
           "Marie Galante",
           "Middle Caicos",
           "Mona",
           "Navassa",
           "Puerto Rico",
           "Sombrero")

ancient_islands<-seq(1,51, by = 3)
extantnative_islands<-seq(2,51,by=3)
modern_islands<-seq(3,51,by=3)

time<-rep(1:3,17)

LA_ancient<-c(1,4,7,16,31,34,49)
LA_ext<-c(2,5,8,17,32,35,50)
LA_mod<-c(3,6,9,18,33,36,51)

GA_ancient<-c(10,13,25,28,40,43,46)
GA_ext<-c(11,14,26,29,41,44,47)
GA_mod<-c(12,15,27,30,42,45,48)

BA_ancient<-c(19,22,37)
BA_ext<-c(20,23,38)
BA_mod<-c(21,24,39)

orenji<-rgb(255,176,0, 128, maxColorValue = 255)
sherblue<-rgb(100,143,255, 128, maxColorValue = 255)
brightmag<-rgb(220,38,127, 128, maxColorValue = 255)
colors<-c(orenji, sherblue, brightmag)

#same colors but darker
orenji<-rgb(255,176,0, 220, maxColorValue = 255)
sherblue<-rgb(100,143,255, 220, maxColorValue = 255)
brightmag<-rgb(220,38,127, 220, maxColorValue = 255)

colors_op<-c(orenji, sherblue, brightmag)
###### end island categorization and color code ######

#read in presence absense data and trait data, along with trait description data
Carib_PA<-read.csv("Dataset_S2_Caribbean_Reptile_PA_Matrix.csv")

#read in caribbean trait category data
Carib_trait_cat_simplehab<-read.csv("Carib_trait_cat.csv")


# read in Caribbean trait data
#the data needs to be subset depending on which habitat is being uses
#Carib_trait_2 corresponds to Carib_trait_cat_simplehab (uses the habitat column)
#Carib_trait corresponds to Carib_trait_cat (uses the habitatex column)

Carib_trait_2<-read.csv("Dataset_S1_Caribbean_Reptile_Traits.csv")

#reformat PA data, take first column and set it as rownames then remove it, turn it into matrix
rownames(Carib_PA)<-Carib_PA[,1]
Carib_PA<-Carib_PA[,-1]
ncol(Carib_PA)
Carib_PA<-as.matrix(Carib_PA)

#format Carib trait data, give them rownames based on the first column, remove first column, and turn colnames of Carib_PA into rownames of Carib_trait (species names)
rownames(Carib_trait_2)<-Carib_trait_2[,1]
colnames(Carib_PA)<-rownames(Carib_trait_2)
Carib_trait_2<-Carib_trait_2[,-1]
Carib_trait_2<-Carib_trait_2[,-1]
Carib_trait_2<-Carib_trait_2[,-1]
Carib_trait_2<-Carib_trait_2[,c(1:6)]
colnames(Carib_trait_2)
Carib_trait_2$sizeclass<-as.factor(Carib_trait_2$sizeclass)
Carib_trait_2$sizeclass<-ordered(Carib_trait_2$sizeclass, levels=c("xxsmall","xsmall","small","medium","large","mega"))
Carib_trait_2$habitat<-as.factor(Carib_trait_2$habitat)
Carib_trait_2$diet<-as.factor(Carib_trait_2$diet)
Carib_trait_2$reprod<-as.factor(Carib_trait_2$reprod)
Carib_trait_2$foraging<-as.factor(Carib_trait_2$foraging)
Carib_trait_2$activity<-as.factor(Carib_trait_2$activity)
colnames(Carib_PA)<-rownames(Carib_trait_2)

#turn species into functional entities using the FE function. this collapses species with similar traits into functional entities
Carib_FE2<-mFD::sp.to.fe(sp_tr = Carib_trait_2,
                         tr_cat = Carib_trait_cat_simplehab,
                         fe_nm_type = "fe_rank",
                         check_input = TRUE)

#create a document that shows which trait combinations correspond to which FEs
write.csv(Carib_FE2$details_fe$fe_codes, "functional entity codes.csv")
#write a document that shows which species belong to which FEs
write.csv(Carib_FE2$sp_fe, "functional entity db.csv")

#compute functional distance between species/FEs 
Carib_fe_dist2 <- mFD::funct.dist(
  sp_tr         = Carib_FE2$fe_tr,
  tr_cat        = Carib_trait_cat_simplehab,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

#write a document that saves the functional distances between species 
write.csv(as.matrix(Carib_fe_dist2), "fe_dist.csv")
#write a document that shows which species belong to which FEs
write.csv(Carib_FE2$sp_fe, "species in each fe.csv")
#write a document that shows the number of species in each FE
write.csv(Carib_FE2$fe_nb_sp, "number of sp per fe.csv")

#compute the functional diversity indices for the assemblages in the dataset 
Carib_FE_div2<-mFD::alpha.fd.fe(
  asb_sp_occ       = Carib_PA, 
  sp_to_fe         = Carib_FE2,
  ind_nm           = c("fred", "fored", "fvuln"),
  check_input      = TRUE,
  details_returned = TRUE) 

#write a document that contains the FE metrics 
write.csv(Carib_FE_div2$asb_fdfe, "Carib-FE-metrics.csv") 
#write a document that shows which FEs are present/absent in an assemblage
write.csv(Carib_FE_div2$details_fdfe$asb_fe_nbsp,"Carib-PA-FE.csv")

fred_ancient<-c()
ancient_fred_index<-c() 

for (i in 1:length(ancient_islands)){
  if(Carib_FE_div2$asb_fdfe[ancient_islands[i],3] > 1){
  fred_ancient<-c(fred_ancient, Carib_FE_div2$asb_fdfe[ancient_islands[i],3])
  ancient_fred_index<-c(ancient_fred_index, i)
  } else {}
}



fred_native_extant<-Carib_FE_div2$asb_fdfe[extantnative_islands[ancient_fred_index],3]
hist(fred_native_extant)
fred_modern<-Carib_FE_div2$asb_fdfe[modern_islands[ancient_fred_index],3]
fred_modern
mean(fred_ancient)
mean(fred_native_extant)
mean(fred_modern)
wilcox.test(fred_ancient, fred_native_extant)
wilcox.test(fred_ancient, fred_modern)


#a matrix that shows us the difference in FE PA between ancient and extant native fauna, this is basically a list of the extinct FEs (they are the "1"s in the datasheet)
write.csv(Carib_FE_div2$details_fdfe$asb_fe_nbsp[ancient_islands,] - Carib_FE_div2$details_fdfe$asb_fe_nbsp[extantnative_islands,], "differences in PA for ancient v extant.csv")

###EXTINCT SPECIES: create two documents, the first one (extinctions-final) tells us the FEs that were lost from an assemblage and how many species within the FE were lost. 
extinct2<-as.matrix(Carib_FE_div2$details_fdfe$asb_fe_nbsp[ancient_islands,] - Carib_FE_div2$details_fdfe$asb_fe_nbsp[extantnative_islands,])
rownames(extinct2)<-isnames
lostfes<-c()
lostfenum<-c()
lostfeisland<-c()

for (i in 1:nrow(extinct2)){
  for (j in 1:ncol(extinct2)){
    if (extinct2[i,j] > 0){
      lostfes<-append(lostfes,colnames(extinct2)[j])
      lostfenum<-append(lostfenum,extinct2[i,j]) 
      lostfeisland<-append(lostfeisland,isnames[i])
    }
  }
}

write.csv(cbind(lostfes,lostfenum,lostfeisland), "extinctions-final.csv")

#calculate the extinct species that correspond to the fes above
extinct_sp2<-Carib_PA[ancient_islands,]-Carib_PA[extantnative_islands,]
rownames(extinct_sp2)<-isnames
lostfes<-c()
lostfenum<-c()
lostfeisland<-c()
extinct_sp<-c()
for (i in 1:nrow(extinct_sp2)){
  for (j in 1:ncol(extinct_sp2)){
    if (extinct_sp2[i,j] > 0){
      lostfes<-append(lostfes,colnames(extinct_sp2)[j])
      lostfenum<-append(lostfenum,extinct_sp2[i,j]) 
      lostfeisland<-append(lostfeisland,isnames[i])
    }
  }
}
extinctsp_island<-cbind(lostfes,lostfenum,lostfeisland)
#list of extinct species per island, ignore fenum
write.csv(extinctsp_island, "extinct-extirpated-species.csv")

###Introduced Species- determine which FEs have introduced species and what islands those FEs are on
introduced2<-as.matrix(Carib_FE_div2$details_fdfe$asb_fe_nbsp[modern_islands,] - Carib_FE_div2$details_fdfe$asb_fe_nbsp[extantnative_islands,])
introduced2
rownames(introduced2)<-isnames
intfes<-c()
intfenum<-c()
intfeisland<-c()

for (i in 1:nrow(introduced2)){
  for (j in 1:ncol(introduced2)){
    if (introduced2[i,j] > 0){
      intfes<-append(intfes,colnames(introduced2)[j])
      intfenum<-append(intfenum,introduced2[i,j]) 
      intfeisland<-append(intfeisland,isnames[i])
    }
  }
}

introduced_FEs2<-cbind(intfes,intfenum,intfeisland)
introduced_FEs2
write.csv(introduced_FEs2, "introduced_FEs.csv")

fspaces_quality_Carib2<-mFD::quality.fspaces(
  sp_dist = Carib_fe_dist2,
  maxdim_pcoa = 10,
  deviation_weighting = "absolute",
  fdendro = "average"
)
fspaces_quality_Carib2$quality_fspaces #4 spaces for simple habitat

#PCo coordinates for each FE
sp_faxes_coord2<-fspaces_quality_Carib2$"details_fspaces"$"sp_pc_coord"
sp_faxes_coord2[,1]#PC coordinates for each fe
plot(sp_faxes_coord2[,1], sp_faxes_coord2[,2])

#functional richness
alpha_fd_indices_Carib2<-mFD::alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord2[,c("PC1","PC2","PC3","PC4")],
  asb_sp_w =  Carib_FE_div2$details_fdfe$asb_fe_nbsp[c(1:12,16:19,21:34,36:40,43,46:48),],
  ind_vect = c("fric", "feve", "fdiv"),
  scaling = TRUE,
  check_input = TRUE,
  details_returned = TRUE) #compute alpha functional diversity indices

write.csv(alpha_fd_indices_Carib2$functional_diversity_indices, "FE_diversityindices.csv")

Carib2_tr_faxes<-mFD::traits.faxes.cor(
  sp_tr = Carib_FE2$fe_tr,
  sp_faxes_coord = sp_faxes_coord2[,c("PC1","PC2","PC3","PC4")],
  plot = TRUE,
)

write.csv(Carib2_tr_faxes$tr_faxes_stat, "Trait and faxes correlation.csv")

alpha_fd_indices_Carib2$functional_diversity_indices
div_FE2<-read.csv("4PC-FE_diversityindices.csv") 

#see if there is a significant difference between richness for different time bins
wilcox.test(div_FE2$fric[ancient_islands],div_FE2$fric[extantnative_islands])
wilcox.test(div_FE2$fric[ancient_islands],div_FE2$fric[modern_islands])
wilcox.test(div_FE2$fric[extantnative_islands],div_FE2$fric[modern_islands])

#see if there is a significant differnce between time bins taking into account archipelago
wilcox.test(div_FE2$fric[GA_ancient],div_FE2$fric[GA_ext])
wilcox.test(div_FE2$fric[GA_ancient],div_FE2$fric[GA_mod])
wilcox.test(div_FE2$fric[GA_ext],div_FE2$fric[GA_mod])

wilcox.test(div_FE2$fric[LA_ancient],div_FE2$fric[LA_ext])  #.05468
wilcox.test(div_FE2$fric[LA_ancient],div_FE2$fric[LA_mod]) #.06494
wilcox.test(div_FE2$fric[LA_ext],div_FE2$fric[LA_mod])

wilcox.test(div_FE2$fric[BA_ancient],div_FE2$fric[BA_ext])
wilcox.test(div_FE2$fric[BA_ancient],div_FE2$fric[BA_mod])
wilcox.test(div_FE2$fric[BA_ext],div_FE2$fric[BA_mod])

##divide island based on size 2000

anc_big<-which(div_FE2$island.area[ancient_islands] > 2000)
div_FE2$island.area[ancient_islands[anc_big]]
anc_sm<-which(div_FE2$island.area[ancient_islands] < 2000)
ext_big<-which(div_FE2$island.area[extantnative_islands] > 2000)
ext_sm<-which(div_FE2$island.area[extantnative_islands] < 2000)
mod_big<-which(div_FE2$island.area[modern_islands] > 2000)
mod_sm<-which(div_FE2$island.area[modern_islands] < 2000)


### sig difference between richness in small ancient v native and small ancient v. modern
boxplot(div_FE2$fric[ancient_islands[anc_sm]],div_FE2$fric[extantnative_islands[ext_sm]], div_FE2$fric[modern_islands[mod_sm]], col=colors, ylab="functional richness", names=c("ancient", "native extant", "modern")) ####FIGURE

boxplot(div_FE2$fric[ancient_islands[anc_sm]],div_FE2$fric[extantnative_islands[ext_sm]], div_FE2$fric[modern_islands[mod_sm]], col=colors, ylab="functional richness", names=c("ancient", "native extant", "modern")) ####FIGURE

wilcox.test(div_FE2$fric[ancient_islands[anc_sm]],div_FE2$fric[extantnative_islands[ext_sm]]) #.01118
wilcox.test(div_FE2$fric[ancient_islands[anc_sm]],div_FE2$fric[modern_islands[mod_sm]]) #.06744
wilcox.test(div_FE2$fric[extantnative_islands[ext_sm]],div_FE2$fric[modern_islands[mod_sm]]) 

wilcox.test(div_FE2$fric[ancient_islands[anc_big]],div_FE2$fric[extantnative_islands[ext_big]]) 
wilcox.test(div_FE2$fric[ancient_islands[anc_big]],div_FE2$fric[modern_islands[mod_big]]) 
wilcox.test(div_FE2$fric[extantnative_islands[ext_big]],div_FE2$fric[modern_islands[mod_big]]) 

wilcox.test(div_FE2$fric[ancient_islands[anc_sm]],div_FE2$fric[ancient_islands[anc_big]]) #.001
wilcox.test(div_FE2$fric[modern_islands[mod_sm]],div_FE2$fric[modern_islands[mod_big]]) #.0027
wilcox.test(div_FE2$fric[extantnative_islands[ext_sm]],div_FE2$fric[extantnative_islands[ext_big]]) #.01056
#on average 5.6% of functional diversity is loss for simple habitat
difference<-div_FE2$fric[ancient_islands[anc_sm]] - div_FE2$fric[extantnative_islands[ext_sm]]
difference
mean(difference, na.rm=TRUE)

boxplot(div_FE2$fric[ancient_islands[anc_big]],div_FE2$fric[extantnative_islands[ext_big]], div_FE2$fric[modern_islands[mod_big]], col=colors, ylab="functional richness", names=c("ancient", "native extant", "modern"))
wilcox.test(div_FE2$fric[ancient_islands[anc_big]],div_FE2$fric[extantnative_islands[ext_big]])
wilcox.test(div_FE2$fric[ancient_islands[anc_big]],div_FE2$fric[modern_islands[mod_big]]) 
##FIG 3 
boxplot(div_FE2$fric[ancient_islands[anc_sm]],div_FE2$fric[extantnative_islands[ext_sm]], div_FE2$fric[modern_islands[mod_sm]], div_FE2$fric[ancient_islands[anc_big]],div_FE2$fric[extantnative_islands[ext_big]], div_FE2$fric[modern_islands[mod_big]], col=colors, ylab="functional richness", names=c("Ancient small islands", "Native extant small islands", "Modern small islands","Ancient large islands", "Native extant large islands", "Modern large islands"))
abline(v =3.5, col="black")
title("Reptile Functional Richness of Caribbean Assemblages", cex=2)


############# simple habitat fig 1 a and b

FE_data2<-read.csv("Carib-FE-metrics.csv")
frich<-read.csv("4PC-FE_diversityindices.csv")
area<-log10(frich$island.area)
nbsp<-log10(FE_data2$nb_sp)
nbfe<-log10(FE_data2$nb_fe)
ric<-log10(frich$fric)
frich
par(mfrow=c(1,2))
#rank order of FES
Carib_FE2$fe_nb_sp
plot(c(1:123),Carib_FE2$fe_nb_sp, type="s", xlim=c(1,125), ylab="# species per FE", col="darkslategray4", lwd=2, xlab="FE number")

#island area with species richness and fe richness
plot(area,nbsp, bg=colors, col=colors_op, pch=24, cex=1, ylim=c(0,2), xlab="log Island Area", ylab="log Richness")
points(area,nbfe, bg=colors, col = colors_op, pch=25, cex=1)
#abline(v=log10(2000), col="gray")

abline(lm(nbsp[ancient_islands]~area[ancient_islands]), col=colors_op[1], lwd=4)
abline(lm(nbsp[extantnative_islands]~area[extantnative_islands]), col=colors_op[2], lwd=4)
abline(lm(nbsp[modern_islands]~area[modern_islands]), col=colors_op[3], lwd=4)
abline(lm(nbfe[ancient_islands]~area[ancient_islands]), col=colors_op[1], lty=3, lwd=4)
abline(lm(nbfe[extantnative_islands]~area[extantnative_islands]), col=colors_op[2], lty=3, lwd=4)
abline(lm(nbfe[modern_islands]~area[modern_islands]), col=colors_op[3], lty=3, lwd=4)

legend(2,.4, c("log species richness", "log FE richness"), pch=c(24,25), lty = c(1,3), bty="n", lwd=2, col="light gray", pt.cex=1, cex=.75)
legend(2.05,.25, c("ancient", "native extant", "modern"), pch=15, col=colors_op, bty="n", pt.cex = 1.5, cex=.75)
### end fig 1 a/b


#construct linear models that look at the relationship between island area and each richness metric 
ancient_lm_nbsp<-lm(nbsp[ancient_islands]~area[ancient_islands])
extant_lm_nbsp<-lm(nbsp[extantnative_islands]~area[extantnative_islands])
modern_lm_nbsp<-lm(nbsp[modern_islands]~area[modern_islands])

ancient_lm_nbfe<-lm(nbfe[ancient_islands]~area[ancient_islands])
extant_lm_nbfe<-lm(nbfe[extantnative_islands]~area[extantnative_islands])
modern_lm_nbfe<-lm(nbfe[modern_islands]~area[modern_islands])

lmparts<-rbind(ancient_lm_nbsp$coefficients,extant_lm_nbsp$coefficients,modern_lm_nbsp$coefficients,ancient_lm_nbfe$coefficients,extant_lm_nbfe$coefficients,modern_lm_nbfe$coefficients)


r.sq<-rbind(summary(ancient_lm_nbsp)$adj.r.squared,summary(extant_lm_nbsp)$adj.r.squared,summary(modern_lm_nbsp)$adj.r.squared,summary(ancient_lm_nbfe)$adj.r.squared,summary(extant_lm_nbfe)$adj.r.squared,summary(modern_lm_nbfe)$adj.r.squared)
pval<-rbind(anova(ancient_lm_nbsp)$'Pr(>F)'[1],anova(extant_lm_nbsp)$'Pr(>F)'[1],anova(modern_lm_nbsp)$'Pr(>F)'[1],anova(ancient_lm_nbfe)$'Pr(>F)'[1],anova(extant_lm_nbfe)$'Pr(>F)'[1],anova(modern_lm_nbfe)$'Pr(>F)'[1])

lmnames<-rbind("ancient_nsbp", "extant_nsbp","modern_nsbp","ancient_nbfe","extant_nbfe","modern_nbfe")
lmoutput<-cbind(lmnames,lmparts,r.sq,pval)
colnames(lmoutput)<-c("","intercept","slope","adjusted r squared","pval")
lmoutput
write.csv(lmoutput,"Linear Models for Area-SPRic and Area-FE.csv")

######### end fig 1a/b


##supplemental frich figure S1
ancient_lm_rich<-lm(frich$fric[ancient_islands]~area[ancient_islands])
extant_lm_rich<-lm(frich$fric[extantnative_islands]~area[extantnative_islands])
modern_lm_rich<-lm(frich$fric[modern_islands]~area[modern_islands])

par(mfrow=c(1,1))
plot(area, frich$fric, col=colors, bg=colors_op,pch=16, xlab="log Island Area", ylab="Functional Richness")
abline(ancient_lm_rich, col=colors_op[1], lwd=4)
abline(extant_lm_rich, col=colors_op[2], lwd=4)
abline(modern_lm_rich, col=colors_op[3], lwd=4)

legend(-.6,.57, c("ancient", "native extant", "modern"), pch=15, col=colors_op, bty="n", pt.cex = 1.5, cex=.75)
text(1, .6, expression("Functional Richness (ancient) = .12 x log (Island Area) - .17, R"^2*" = .76"), cex=.75)
text(1.12, .585, expression("Functional Richness (native extant) = .13 x log (Island Area) - .27, R"^2*" = .81"), cex=.75)
text(1, .57, expression("Functional Richness (modern) = .15 x log (Island Area) - .26, R"^2*" = .78"), cex=.75)


anc_ric_lm<-cbind(rbind(ancient_lm_rich$coefficients), summary(ancient_lm_rich)$adj.r.squared, anova(ancient_lm_rich)$'Pr(>F)'[1])
ext_ric_lm<-cbind(rbind(extant_lm_rich$coefficients), summary(extant_lm_rich)$adj.r.squared, anova(extant_lm_rich)$'Pr(>F)'[1])
mod_ric_lm<-cbind(rbind(modern_lm_rich$coefficients), summary(modern_lm_rich)$adj.r.squared, anova(modern_lm_rich)$'Pr(>F)'[1])

frichness_lm<-rbind(anc_ric_lm,ext_ric_lm,mod_ric_lm)
rownames(frichness_lm)<-c("ancient","extant", "modern")
colnames(frichness_lm)<-c("intercept","slope","adjusted r-squared", "pval")
write.csv(frichness_lm, "linear model info for frichness and island area.csv")

##end supplemental fig 2

###begin figure 2 ###
FE_int<-read.csv("introduced FE database.csv")
FE_loss<-read.csv("lostFEdb.csv")  

##FIGURE 2

orenji2<-rgb(255,176,0, 175, maxColorValue = 255)
sherblue2<-rgb(100,143,255, 175, maxColorValue = 255)
brightmag2<-rgb(220,38,127, 175, maxColorValue = 255)
par(mfrow=c(1,2))
#1 &2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], col="white", xlim=c(-.6,.6), ylim=c(-.6,.6), pch=18,cex=2, xlab="PCo 1", ylab="PCo 2")
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")
points(sp_faxes_coord2[c(FE_loss$fe_num),1] ,sp_faxes_coord2[c(FE_loss$fe_num),2], pch=19, col=orenji2, cex=2) #scaled based on %FE extinct
polygon(sp_faxes_coord2[FE_loss$fe_num[chull(sp_faxes_coord2[c(FE_loss$fe_num),1:2])],1:2], col=colors[1], border=NA)
#text(sp_faxes_coord2[c(FE_loss$fe_num),1] ,sp_faxes_coord2[c(FE_loss$fe_num),2], labels = FE_loss$fe_num, cex=.75)

points(sp_faxes_coord2[c(FE_int$fe),1] ,sp_faxes_coord2[c(FE_int$fe),2], pch=19,  col=brightmag2, cex=2)
polygon(sp_faxes_coord2[FE_int$fe[chull(sp_faxes_coord2[c(FE_int$fe),1:2])],1:2], col=colors[3], border=NA)
#text(sp_faxes_coord2[c(FE_int$fe),1] ,sp_faxes_coord2[c(FE_int$fe),2], labels = FE_int$fe, cex=.75)

#3&4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], col="white", xlim=c(-.6,.6), ylim=c(-.6,.6), pch=18,cex=2, xlab="PCo 3", ylab="PCo 4")
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")
points(sp_faxes_coord2[c(FE_loss$fe_num),3] ,sp_faxes_coord2[c(FE_loss$fe_num),4], pch=19, col=orenji2, cex=2) 
polygon(sp_faxes_coord2[FE_loss$fe_num[chull(sp_faxes_coord2[c(FE_loss$fe_num),3:4])],3:4], col=colors[1], border=NA)
#text(sp_faxes_coord2[c(FE_loss$fe_num),3] ,sp_faxes_coord2[c(FE_loss$fe_num),4], labels = FE_loss$fe_num, cex=.75)
points(sp_faxes_coord2[c(FE_int$fe),3] ,sp_faxes_coord2[c(FE_int$fe),4], pch=19,  col=brightmag2, cex=2)
polygon(sp_faxes_coord2[FE_int$fe[chull(sp_faxes_coord2[c(FE_int$fe),3:4])],3:4], col=colors[3], border=NA)
#text(sp_faxes_coord2[c(FE_int$fe),3] ,sp_faxes_coord2[c(FE_int$fe),4], labels = FE_int$fe, cex=.75)

legend(0,-.4,c("Extinct/Extirpated FEs","Introduced FEs"), pch=19,col= colors[c(1,3)], bty="n", cex=1.5)

## data for FIGURE 4

prop<-(FE_data2$nb_fe[ancient_islands] - FE_data2$nb_fe[extantnative_islands])/FE_data2$nb_fe[ancient_islands]
names(prop)<-isnames
sort(prop)
1 - prop
#could plot circles as 1 and then 1- proportion inside to show change PART OF FIGURE
par(mfrow=c(1,1))
plot(seq(1,425, by=25), FE_data2$nb_fe[ancient_islands], ylim=c(-1,80),cex=10, bg="white",axes=F, xlab="",ylab="", pch=1)
points(seq(1,425, by=25), FE_data2$nb_fe[ancient_islands], cex=10*(1-prop), pch=21, bg=colors[2], col=colors[2])
text(seq(1,425, by=25), -.5, isnames)
points(rep(400,4), rep(65,4), cex=c(1,5,9,10))

t.test(sp_faxes_coord2[c(FE_loss$fe_num),3], sp_faxes_coord2[c(FE_int$fe),3])
wilcox.test(sp_faxes_coord2[c(FE_loss$fe_num),3], sp_faxes_coord2[c(FE_int$fe),3])

#####  functional richness of introduced and extinct
test<-read.csv("int_ext_matrices.csv")

rownames(test)<-test[,1]
test<-test[,-c(1)]
test<-as.matrix(test)
test
testrich<-mFD::alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord2[,c("PC1","PC2","PC3","PC4")],
  asb_sp_w =  test,
  ind_vect = c("fric", "feve", "fdiv"),
  scaling = TRUE,
  check_input = TRUE,
  details_returned = TRUE) #compute alpha functional diversity indices


testrich$functional_diversity_indices

write.csv(testrich$functional_diversity_indices, "functional diversity indices for extinct and introduced.csv")

###### Supplemental figure S2 ######

par(mfrow=c(5,6))
###fix indices to match the archipelagos that are actually in the data...ancient_islands etc don't work b/c of missing archipelagos
isnamfin<-isnames[c(1:4,6:16)]
isnames


orenji<-rgb(255,176,0, 150, maxColorValue = 255)
sherblue<-rgb(100,143,255, 0, maxColorValue = 255)
brightmag<-rgb(220,38,127, 70, maxColorValue = 255)

colors<-c(orenji, sherblue, brightmag)

orenji<-rgb(255,176,0, 0, maxColorValue = 255)
sherblue<-rgb(100,143,255, 255, maxColorValue = 255)
brightmag<-rgb(220,38,127, 0, maxColorValue = 255)

colors_op<-c(orenji, sherblue, brightmag)
#Anguilla - Grand terre
for (i in 1:5){
  
  #1 and 2
  plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[i], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
  hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
  polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")
  
  ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[ancient_islands[i]]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[ancient_islands[i]]),2])
  points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[ancient_islands[i]]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[ancient_islands[i]]),2], pch=".", col=colors[1])
  is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[ancient_islands[i]])
  polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])
  
  ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),2])
  points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),2], pch=".", col=colors[2])
  is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]])
  polygon(sp_faxes_coord2[is2[ang2],1:2],col=colors[2])
  
  ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[modern_islands[i]]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[modern_islands[i]]),2])
  points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[modern_islands[i]]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[modern_islands[i]]),2], pch=".", col=colors[3])
  is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[modern_islands[i]])
  polygon(sp_faxes_coord2[is3[ang3],1:2],col=colors[3], border=colors[3])
  
  #3 and 4
  plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], pch=".",col="lightgray", axes=T, main=isnamfin[i], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
  hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
  polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")
  
  ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[ancient_islands[i]]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[ancient_islands[i]]),4])
  points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[ancient_islands[i]]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[ancient_islands[i]]),4], pch=".", col=colors[1])
  is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[ancient_islands[i]])
  polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])
  
  ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),4])
  points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),4], pch=".", col=colors[])
  is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]])
  polygon(sp_faxes_coord2[is2[ang2],3:4],col=colors[2])
  
  ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[modern_islands[i]]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[modern_islands[i]]),4])
  points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[modern_islands[i]]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[modern_islands[i]]),4], pch=".", col=colors[3])
  is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[modern_islands[i]])
  polygon(sp_faxes_coord2[is3[ang3],3:4],col=colors[3], border=colors[3])
  
  
}
#Grand Turk ancient and modern

#1 and 2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[6], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[16]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[16]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[16]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[16]),2], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[16])
polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])

#ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),2])
#points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),2], pch=".", col=colors[2])
#is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]])
#polygon(sp_faxes_coord2[is2[ang2],1:2],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[17]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[17]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[17]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[17]),2], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[17])
polygon(sp_faxes_coord2[is3[ang3],1:2],col=colors[3], border=colors[3])

#3 and 4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], pch=".",col="lightgray", axes=T, main=isnamfin[6], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[16]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[16]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[16]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[16]),4], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[16])
polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])

#ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),4])
#points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]]),4], pch=".", col=colors[])
#is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[extantnative_islands[i]])
#polygon(sp_faxes_coord2[is2[ang2],3:4],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[17]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[17]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[17]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[17]),4], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[17])
polygon(sp_faxes_coord2[is3[ang3],3:4],col=colors[3], border=colors[3])

#Great Abaco

#1 and 2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[7], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[18]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[18]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[18]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[18]),2], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[18])
polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[19]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[19]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[19]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[19]),2], pch=".", col=colors[2])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[19])
polygon(sp_faxes_coord2[is2[ang2],1:2],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[20]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[20]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[20]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[20]),2], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[20])
polygon(sp_faxes_coord2[is3[ang3],1:2],col=colors[3], border=colors[3])

#3 and 4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], pch=".",col="lightgray", axes=T, main=isnamfin[7], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[18]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[18]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[18]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[18]),4], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[18])
polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[19]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[19]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[19]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[19]),4], pch=".", col=colors[])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[19])
polygon(sp_faxes_coord2[is2[ang2],3:4],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[20]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[20]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[20]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[20]),4], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[20])
polygon(sp_faxes_coord2[is3[ang3],3:4],col=colors[3], border=colors[3])


#hispaniola 

#1 and 2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[8], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[21]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[21]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[21]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[21]),2], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[21])
polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[22]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[22]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[22]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[22]),2], pch=".", col=colors[2])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[22])
polygon(sp_faxes_coord2[is2[ang2],1:2],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[23]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[23]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[23]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[23]),2], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[23])
polygon(sp_faxes_coord2[is3[ang3],1:2],col=colors[3], border=colors[3])

#3 and 4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], pch=".",col="lightgray", axes=T, main=isnamfin[8], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[21]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[21]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[21]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[21]),4], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[21])
polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[22]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[22]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[22]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[22]),4], pch=".", col=colors[])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[22])
polygon(sp_faxes_coord2[is2[ang2],3:4],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[23]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[23]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[23]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[23]),4], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[23])
polygon(sp_faxes_coord2[is3[ang3],3:4],col=colors[3], border=colors[3])


#jamaica

#1 and 2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[9], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[24]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[24]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[24]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[24]),2], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[24])
polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[25]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[25]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[25]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[25]),2], pch=".", col=colors[2])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[25])
polygon(sp_faxes_coord2[is2[ang2],1:2],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[26]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[26]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[26]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[26]),2], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[26])
polygon(sp_faxes_coord2[is3[ang3],1:2],col=colors[3], border=colors[3])

#3 and 4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], pch=".",col="lightgray", axes=T, main=isnamfin[9], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[24]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[24]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[24]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[24]),4], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[24])
polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[25]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[25]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[25]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[25]),4], pch=".", col=colors[])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[25])
polygon(sp_faxes_coord2[is2[ang2],3:4],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[26]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[26]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[26]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[26]),4], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[26])
polygon(sp_faxes_coord2[is3[ang3],3:4],col=colors[3], border=colors[3])

#la desirade

#1 and 2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[10], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[27]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[27]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[27]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[27]),2], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[27])
polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[28]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[28]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[28]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[28]),2], pch=".", col=colors[2])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[28])
polygon(sp_faxes_coord2[is2[ang2],1:2],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[29]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[29]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[29]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[29]),2], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[29])
polygon(sp_faxes_coord2[is3[ang3],1:2],col=colors[3], border=colors[3])

#3 and 4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], pch=".",col="lightgray", axes=T, main=isnamfin[10], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[27]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[27]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[27]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[27]),4], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[27])
polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[28]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[28]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[28]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[28]),4], pch=".", col=colors[])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[28])
polygon(sp_faxes_coord2[is2[ang2],3:4],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[29]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[29]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[29]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[29]),4], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[29])
polygon(sp_faxes_coord2[is3[ang3],3:4],col=colors[3], border=colors[3])

#marie galante

#1 and 2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[11], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[30]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[30]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[30]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[30]),2], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[30])
polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])

#ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[null]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[null]),2])
#points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[null]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[null]),2], pch=".", col=colors[2])
#is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[null])
#polygon(sp_faxes_coord2[is2[ang2],1:2],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[31]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[31]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[31]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[31]),2], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[31])
polygon(sp_faxes_coord2[is3[ang3],1:2],col=colors[3], border=colors[3])

#3 and 4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], pch=".",col="lightgray", axes=T, main=isnamfin[11], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[30]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[30]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[30]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[30]),4], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[30])
polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])

#ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[null]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[null]),4])
#points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[null]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[null]),4], pch=".", col=colors[])
#is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[null])
#polygon(sp_faxes_coord2[is2[ang2],3:4],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[31]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[31]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[31]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[31]),4], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[31])
polygon(sp_faxes_coord2[is3[ang3],3:4],col=colors[3], border=colors[3])

#middle caicos

#1 and 2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[12], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[32]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[32]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[32]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[32]),2], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[32])
polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[33]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[33]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[33]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[33]),2], pch=".", col=colors[2])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[33])
polygon(sp_faxes_coord2[is2[ang2],1:2],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),2], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34])
polygon(sp_faxes_coord2[is3[ang3],1:2],col=colors[3], border=colors[3])

#3 and 4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], pch=".",col="lightgray", axes=T, main=isnamfin[12], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[32]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[32]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[32]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[32]),4], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[32])
polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[33]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[33]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[33]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[33]),4], pch=".", col=colors[])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[33])
polygon(sp_faxes_coord2[is2[ang2],3:4],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),4], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34])
polygon(sp_faxes_coord2[is3[ang3],3:4],col=colors[3], border=colors[3])

#mona

#1 and 2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[13], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[35]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[35]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[35]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[35]),2], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[35])
polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])


#3 and 4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4],pch=".",col="lightgray", axes=T, main=isnamfin[13], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[35]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[35]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[35]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[35]),4], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[35])
polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])

#navassa

#1 and 2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[14], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[36]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[36]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[36]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[36]),2], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[36])
polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])


#3 and 4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], pch=".",col="lightgray", axes=T, main=isnamfin[14], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[36]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[36]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[36]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[36]),4], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[36])
polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])


#puerto rico 

#1 and 2
plot(sp_faxes_coord2[,1],sp_faxes_coord2[,2], pch=".",col="lightgray", axes=T, main=isnamfin[15], xlab="PCoA 1", ylab="PCoA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,1], sp_faxes_coord2[,2])
polygon(sp_faxes_coord2[hpts,1:2],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[37]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[37]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[37]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[37]),2], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[37])
polygon(sp_faxes_coord2[is1[ang1],1:2],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[38]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[38]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[38]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[38]),2], pch=".", col=colors[2])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[38])
polygon(sp_faxes_coord2[is2[ang2],1:2],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),2])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),1],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),2], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34])
polygon(sp_faxes_coord2[is3[ang3],1:2],col=colors[3], border=colors[3])

#3 and 4
plot(sp_faxes_coord2[,3],sp_faxes_coord2[,4], pch=".",col="lightgray", axes=T, main=isnamfin[15], xlab="PCoA 3", ylab="PCoA 4", xlim=c(-.6,.6), ylim=c(-.6,.6))
hpts<-chull(sp_faxes_coord2[,3], sp_faxes_coord2[,4])
polygon(sp_faxes_coord2[hpts,3:4],border="lightgray")

ang1<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[37]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[37]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[37]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[37]),4], pch=".", col=colors[1])
is1<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[37])
polygon(sp_faxes_coord2[is1[ang1],3:4],col=colors[1], border=colors[1])

ang2<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[38]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[38]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[38]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[38]),4], pch=".", col=colors[])
is2<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[38])
polygon(sp_faxes_coord2[is2[ang2],3:4],col=colors[2])

ang3<-chull(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),4])
points(sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),3],sp_faxes_coord2[unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34]),4], pch=".", col=colors[3])
is3<-unlist(alpha_fd_indices_Carib2$details$asb_vert_nm[34])
polygon(sp_faxes_coord2[is3[ang3],3:4],col=colors[3], border=colors[3])


##legend for figure
par(mfrow=c(1,1))
#plot(c(10, 20, 30, 40),rep(1,4), pch=c(15, 22,15,22), col=c(colors[1],"black",colors[3],"lightgray"), bg=c(colors,NA), xlim =c(0,60), cex=3)
#text(c(11.8,22.5,32,43.4),rep(1,4),labels=c("Ancient","Native Extant", "Modern","Global Trait Space"))

### END Supplemental Figure 2 ######
