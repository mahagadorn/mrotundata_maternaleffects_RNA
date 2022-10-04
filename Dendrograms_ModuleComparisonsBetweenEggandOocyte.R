####MAKE MANUSCRIPT DENDROGRAM FIGURE#####
## M.A. HAGADORN
## last modified - 5/12/2022
## Manuscript: Maternal RNA M. rotundata

#The script imports the WGCNA data (only the first time) and then saves what we need for generating the figure
#Reimports only that data to save space
#Then we order the names of the original dendrogram leaves so we can appropriately rename them
#draw dendrograms and labels + connecting lines






#OOCYTE
# #Module Eigengenes list = MEs
# load("Mrotundata_OOCYTE_step2_ModuleAssignment&Merge.RData")
# 
# MEs_oocytes <- MEs
# 
# MEs_oocytes_MEDiss <- 1 - cor(MEs_oocytes, use = "pairwise.complete.obs")
# MEs_oocytes_MEDiss_noNA <- subset(MEs_oocytes_MEDiss, select = -c(MEgrey)) #MEgrey has NaN--divided by zero
# MEs_oocytes_MEDiss_noNA <- MEs_oocytes_MEDiss_noNA[!(row.names(MEs_oocytes_MEDiss_noNA) %in% "MEgrey"),]
# 
# MEOocyte_Tree <- hclust(as.dist(MEs_oocytes_MEDiss_noNA), method="average")
# plot(MEOocyte_Tree)
# 
# moduleLabels_oocytes <- moduleLabels
# moduleColors_oocytes <- moduleColors
# 
# 
# # Save module colors and labels for use in subsequent parts
# save(MEs_oocytes, moduleLabels_oocytes, moduleColors_oocytes, MEOocyte_Tree, MEs_oocytes_MEDiss_noNA,
#      file = "Mrotundata_OOCYTE_dendrograminfo")
# 
# 
# 
#EGG
# #Module Eigengenes list = MEs
# load("Mrotundata_EGG_step2_ModuleAssignment&Merge.RData")
# 
# MEs_eggs <- MEs
# 
# MEs_eggs_MEDiss <- 1 - cor(MEs_eggs, use = "pairwise.complete.obs")
# MEs_eggs_MEDiss_noNA <- subset(MEs_eggs_MEDiss, select = -c(MEgrey)) #MEgrey has NaN--divided by zero
# MEs_eggs_MEDiss_noNA <- MEs_eggs_MEDiss_noNA[!(row.names(MEs_eggs_MEDiss_noNA) %in% "MEgrey"),]
# 
# MEEgg_Tree <- hclust(as.dist(MEs_eggs_MEDiss_noNA), method="average")
# plot(MEEgg_Tree)
# 
# moduleLabels_eggs <- moduleLabels
# moduleColors_eggs <- moduleColors
# 
# 
# # Save module colors and labels for use in subsequent parts
# save(MEs_eggs, moduleLabels_eggs, moduleColors_eggs, MEEgg_Tree, MEs_eggs_MEDiss_noNA,
#      file = "Mrotundata_EGG_dendrograminfo")



##Load in just the data we want:
load("Mrotundata_OOCYTE_dendrograminfo")
load("Mrotundata_EGG_dendrograminfo")


#important table with module assignment for each gene
oldnew <- read.csv2("ModuleOverlap/ModuleAssignment_OldNewNaming.csv", header=TRUE, sep=",")
rownames(oldnew) <- oldnew$X
oldnew <- oldnew[,-1]


#Look up table for old new name matching for egg modules
new.modname <- unique(oldnew$new)

ModLookup <- data.frame(new.modname, old.modname = rep(NA, length(new.modname)))

for(i in new.modname){
  ModLookup[i,2] <- head(oldnew[which(oldnew$new == i),2], 1)
}  
  

ModLookup <- data.frame(cbind(new.modname = ModLookup[1:27,1], old.modname = ModLookup[28:54,2]))
#Verified




##########################################
##Assign Simplified Names for Manuscript##
##########################################
Oocyte.modules <- unique(oldnew$Oocyte)
list.oocytemodules <- cbind(rep("oocyteMod", length(Oocyte.modules)), seq_along(Oocyte.modules))
Oocyte.modules_simplifiednames <- cbind(originalnames=Oocyte.modules, simplifiednames=paste0(list.oocytemodules[,1],list.oocytemodules[,2]))

Egg.modules <- unique(oldnew$new)
list.eggmodules <- cbind(rep("eggMod", length(Egg.modules)), seq_along(Egg.modules))
Egg.modules_simplifiednames <- cbind(newmatchednames=Egg.modules, simplifiednames=paste0(list.eggmodules[,1],list.eggmodules[,2]))

NewModNames.formanuscript <- rbind(Oocyte.modules_simplifiednames, Egg.modules_simplifiednames)
colnames(NewModNames.formanuscript) <- c("WGCNA_ColorNamesMatched", "Simplified_NumericNames")

write.table(NewModNames.formanuscript, "SimplifiedNumericModuleNamesForManuscriptUse_newmodname.csv", sep = ",", row.names = FALSE)






#Make Dendrograms
l <- length(MEOocyte_Tree$order) + 1  #this leaves one extra space at the end so that the extra 6% from the y axis can be removed and that lables stay in their appropriate location

# # The matrix to draw the arrows:
# cbind((1:l)[order(hc1$order)],(1:l)[order(hc2$order)]) -> ord_arrow

# The two vectors of ordered leave labels:
#get order for oocyte
leavesoocyte <- MEOocyte_Tree$labels[MEOocyte_Tree$order]
leavesoocyte.simplified <- gsub("ME","", leavesoocyte)
leavesoocyte.simplified2 <- NewModNames.formanuscript[1:30,2]


#get order for egg
leavesegg <- MEEgg_Tree$labels[MEEgg_Tree$order]
leavesegg.simplified <- gsub("ME","", leavesegg)
leavesegg.simplified2 <- NewModNames.formanuscript[31:57,2]


#match the order of the dendrogram eggs
ModLookup2 <- cbind(ModLookup, leavesegg.simplified2)
leavesegg.newposition <- ModLookup2[match(leavesegg.simplified, ModLookup$old.modname),]
#subset out the order for new labels
leavesegg.newnames <- leavesegg.newposition$leavesegg.simplified2

#match the order of the dendrogram oocyte
ModLookup.oocyte <- as.data.frame(cbind(Oocyte.modules, leavesoocyte.simplified2))
leavesoocyte.newposition <- ModLookup.oocyte[match(leavesoocyte.simplified, ModLookup.oocyte$Oocyte.modules),]
#subset out the order for new labels
leavesoocyte.newnames <- leavesoocyte.newposition$leavesoocyte.simplified2
leavesoocyte.newnames <- c(leavesoocyte.newnames, NA) #this leaves one extra space at the end so that the extra 6% from the y axis can be removed and that lables stay in their appropriate location


NewLabels.MEs_eggs_MEDiss_noNA <- MEs_eggs_MEDiss_noNA
rownames(NewLabels.MEs_eggs_MEDiss_noNA) <- leavesegg.newnames
##NEW EGG TREE WITH PROPER LABELS
MEEgg_Tree.newLabels <- hclust(as.dist(NewLabels.MEs_eggs_MEDiss_noNA), method="average")
plot(MEEgg_Tree.newLabels)


NoMe.MEs_oocytes_MEDiss_noNA <- MEs_oocytes_MEDiss_noNA
rownames(NoMe.MEs_oocytes_MEDiss_noNA) <- leavesoocyte.newnames
##NEW OOCYTE TREE WITH PROPER LABELS
MEOocyte_Tree.NoMELabels <- hclust(as.dist(NoMe.MEs_oocytes_MEDiss_noNA), method="average")
plot(MEOocyte_Tree.NoMELabels)

########################################
##make a table to see where genes went##
########################################

library(dplyr)
#Grey60
Grey60Distribution <- as.data.frame(table(oldnew[oldnew$Oocyte=="grey60", 3]))
write.table(Grey60Distribution, "ModuleOverlap/Grey60DistributionAcross.csv", sep=",")

#violet
VioletDistribution <- as.data.frame(table(oldnew[oldnew$Oocyte=="violet", 3]))
write.table(VioletDistribution, "ModuleOverlap/VioletDistributionAcross.csv", sep=",")

#Midnightblue
MidnightblueDistribution <- as.data.frame(table(oldnew[oldnew$Oocyte=="midnightblue", 3]))
write.table(MidnightblueDistribution, "ModuleOverlap/MidnightblueDistributionAcross.csv", sep=",")

#DarkTurquoise
DarkTurquoiseDistribution <- as.data.frame(table(oldnew[oldnew$Oocyte=="darkturquoise", 3]))
write.table(DarkTurquoiseDistribution, "ModuleOverlap/DarkTurquoiseDistributionAcross.csv", sep=",")



###############
##Dendrograms##
###############

#file containing distribution information and color coding
alldist <- read.csv2(file = "ModuleOverlap/AllSigModulesDistributionAcross.csv", header=TRUE, sep=",")

library(scales)
# And the plot:
jpeg(file="../MRot_WGCNA/ModuleOverlap/ComparisonDendrograms_moduleredistribution.jpeg", width=7.5, height=9.5, units="in", res=1500)
  

layout(matrix(1:5,nrow=1),width=c(5,2.15,4,1.7,5.1))

Add.Line <- function(all.dist=all.dist, line.num, color){
  lines(y=c(alldist$y.1[line.num], alldist$y.2[line.num]), x=c(0,1), lwd = alldist$linewdt[line.num], col=alpha(alldist$Color[line.num],1.0))
}

Make.Box <- function(line, color){
  xx = c(0,1,1,0)
  yy = c(line-.45,line-.45,line+.45,line+.45)
  # When density=0, col refers to the line colour
  polygon(xx, yy, density=0, col=color, lwd=2)
}


###color option 1
# c1 <- c("#E69F00")
# c2 <- c("#0072B2")
# c3 <- c("#009E73")
# c4 <- c("#CC79A7")
# # 
# ###color option 2
# c1 <- c("#661100")
# c2 <- c("#009E73")
# c3 <- c("#CC79A7")
# c4 <- c("#F0E442")
# # 
# ###color option 3
# c1 <- c("#661100")
# c2 <- c("#009E73")
# c3 <- c("#CC79A7")
# c4 <- c("#E69F00")
# # 
# ###color option 4
# c1 <- c("#661100")
# c2 <- c("#D55E00")
# c3 <- c("#6699CC")
# c4 <- c("#999933")
# # 
# # ###color option 5
# c1 <- c("#661100")
# c2 <- c("#E69F00")
# c3 <- c("#191970")
# c4 <- c("#808080")

# # ###color option 6
# c1 <- c("#023663")
# c2 <- c("#E0AB09")
# c3 <- c("#D24982")
# c4 <- c("#0FB79A")

par <- c(0, 0, 0, 0)
# The first dendrogram:
par(mai = par, yaxs="i") #b,l,t,r; mar=c(3,3,3,0),
Oocyte <- as.dendrogram(MEOocyte_Tree.NoMELabels)
plot(Oocyte, horiz=TRUE, ylim=c(0,l), leaflab="none", xlab = "", axes=F)
title(main="Oocyte Modules", line=-1.5, cex.main=1.2, adj=.7)

# The first serie of labels (i draw them separately because, for the second serie, I didn't find a simple way to draw them nicely on the cluster):
par(mai = par) #mar=c(3,0,3,0)
plot(NA, bty="n",axes=FALSE,xlim=c(0,1), ylim=c(0,l),ylab="",xlab="")
#sapply(1:l,function(x)text(x=0,y=x,labels=leavesoocyte.newnames[x], pos=4, cex=1.6))
text(x=0,y=c(1:9,11:12,15,17:29),labels=leavesoocyte.newnames[c(1:9,11:12,15,17:29)], pos=4, cex=1, font=1)
text(x=0,y=c(10,13:14,16),labels=leavesoocyte.newnames[c(10,13:14,16)], pos=4, cex=1, font=2)

#draw boxes around bolded names
Make.Box(10, "#D24982") #oocyte29
Make.Box(13, "#E0AB09") #oocyteMod4
Make.Box(14, "#0FB79A") #oocyteMod20
Make.Box(16, "#023663") #oocyteMod16


# Space for lines
par(mai = par)
#linetoLSB1
plot(y=c(alldist$y.1[1], alldist$y.2[1]), x=c(0,1), type="l", lwd = alldist$linewdt[1], axes=FALSE, xlim=c(0,1), ylim=c(0,l),ylab="",xlab="", col=alldist$Color[1])
Add.Line(all.dist = all.dist, line.num = 2)
Add.Line(all.dist = all.dist, line.num = 3)
Add.Line(all.dist = all.dist, line.num = 4)
Add.Line(all.dist = all.dist, line.num = 5)
Add.Line(all.dist = all.dist, line.num = 6)
Add.Line(all.dist = all.dist, line.num = 7)
Add.Line(all.dist = all.dist, line.num = 8)
Add.Line(all.dist = all.dist, line.num = 9)
Add.Line(all.dist = all.dist, line.num = 11) #skip 10 b/c NA
Add.Line(all.dist = all.dist, line.num = 12)
Add.Line(all.dist = all.dist, line.num = 13)
Add.Line(all.dist = all.dist, line.num = 14)
Add.Line(all.dist = all.dist, line.num = 15)
Add.Line(all.dist = all.dist, line.num = 16)
Add.Line(all.dist = all.dist, line.num = 17)
Add.Line(all.dist = all.dist, line.num = 18)
Add.Line(all.dist = all.dist, line.num = 19)
Add.Line(all.dist = all.dist, line.num = 20)
Add.Line(all.dist = all.dist, line.num = 21)
Add.Line(all.dist = all.dist, line.num = 22)
Add.Line(all.dist = all.dist, line.num = 23)
Add.Line(all.dist = all.dist, line.num = 24)
Add.Line(all.dist = all.dist, line.num = 25)
Add.Line(all.dist = all.dist, line.num = 26)

Add.Line(all.dist = all.dist, line.num = 27)
Add.Line(all.dist = all.dist, line.num = 28)
Add.Line(all.dist = all.dist, line.num = 29)
Add.Line(all.dist = all.dist, line.num = 30)
Add.Line(all.dist = all.dist, line.num = 31)
Add.Line(all.dist = all.dist, line.num = 32)
Add.Line(all.dist = all.dist, line.num = 33)
Add.Line(all.dist = all.dist, line.num = 34)
Add.Line(all.dist = all.dist, line.num = 35)
Add.Line(all.dist = all.dist, line.num = 36)
Add.Line(all.dist = all.dist, line.num = 37)
Add.Line(all.dist = all.dist, line.num = 38)
Add.Line(all.dist = all.dist, line.num = 39)
Add.Line(all.dist = all.dist, line.num = 40)
Add.Line(all.dist = all.dist, line.num = 41)
Add.Line(all.dist = all.dist, line.num = 42)
Add.Line(all.dist = all.dist, line.num = 43)
Add.Line(all.dist = all.dist, line.num = 44)
Add.Line(all.dist = all.dist, line.num = 45)
Add.Line(all.dist = all.dist, line.num = 46)

Add.Line(all.dist = all.dist, line.num = 47)
Add.Line(all.dist = all.dist, line.num = 48)
Add.Line(all.dist = all.dist, line.num = 49)
Add.Line(all.dist = all.dist, line.num = 50)
Add.Line(all.dist = all.dist, line.num = 51)
Add.Line(all.dist = all.dist, line.num = 52)
Add.Line(all.dist = all.dist, line.num = 53)
Add.Line(all.dist = all.dist, line.num = 54)
Add.Line(all.dist = all.dist, line.num = 55)
Add.Line(all.dist = all.dist, line.num = 56)
Add.Line(all.dist = all.dist, line.num = 57)
Add.Line(all.dist = all.dist, line.num = 58)
Add.Line(all.dist = all.dist, line.num = 59)
Add.Line(all.dist = all.dist, line.num = 60)

Add.Line(all.dist = all.dist, line.num = 61)
Add.Line(all.dist = all.dist, line.num = 62)
Add.Line(all.dist = all.dist, line.num = 63)
Add.Line(all.dist = all.dist, line.num = 64)
Add.Line(all.dist = all.dist, line.num = 65)
Add.Line(all.dist = all.dist, line.num = 66)
Add.Line(all.dist = all.dist, line.num = 67)
Add.Line(all.dist = all.dist, line.num = 68)
Add.Line(all.dist = all.dist, line.num = 70) #skip 69 b/c NA
Add.Line(all.dist = all.dist, line.num = 71)
Add.Line(all.dist = all.dist, line.num = 72)
Add.Line(all.dist = all.dist, line.num = 73)
Add.Line(all.dist = all.dist, line.num = 74)
Add.Line(all.dist = all.dist, line.num = 75)
Add.Line(all.dist = all.dist, line.num = 76)
Add.Line(all.dist = all.dist, line.num = 77)
Add.Line(all.dist = all.dist, line.num = 78)
Add.Line(all.dist = all.dist, line.num = 79)
Add.Line(all.dist = all.dist, line.num = 80)
Add.Line(all.dist = all.dist, line.num = 81)
Add.Line(all.dist = all.dist, line.num = 82)

# The second serie of labels:
par(mai = par)
plot(NA, bty="n",axes=FALSE, xlim=c(0,1), ylim=c(0,l), ylab="",xlab="")
#sapply(1:l,function(x)text(x=1,y=x,labels=leavesegg.newnames[x], pos=2, cex=1.6, font=2))
text(x=1,y=c(1:2,9:11,13,16:26),labels=leavesegg.newnames[c(1:2,9:11,13,16:26)], pos=2, cex=1, font=1)
text(x=1,y=c(3:8,12,14:15),labels=leavesegg.newnames[c(3:8,12,14:15)], pos=2, cex=1, font=2)

#draw boxes around bolded names
Make.Box(3, "grey45") #eggMod20
Make.Box(4, "grey45") #eggMod15
Make.Box(5, "grey45") #eggMod4
Make.Box(6, "grey45") #eggMod24
Make.Box(7, "grey45") #eggMod13
Make.Box(8, "grey45") #eggMod25
Make.Box(8, "grey45") #eggMod25
Make.Box(12, "grey45") #eggMod17
Make.Box(14, "grey45") #eggMod23
Make.Box(15, "grey45") #eggMod5



# And the second dendrogram (to reverse it I reversed the xlim vector:
par(mai = par)
plot(as.dendrogram(MEEgg_Tree.newLabels), horiz=TRUE, xlim=c(0,1.07), leaflab="none", ylim=c(0,l), axes=F)
title(main="Egg Modules", line=-1.5, cex.main=1.2, adj=.3)


dev.off()


write.table(Oocyte.modules_simplifiednames, "Oocyte.modules_simplifiednamesformanuscript.csv", sep = ",")



