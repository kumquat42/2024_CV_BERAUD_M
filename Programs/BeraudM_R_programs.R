# 2020 - BeraudM_Lipofabrik-PCA_CAH_Données Transcriptomiques

rm (list = ls())

library(ade4TkGUI)
library(FactoMineR)
library (Factoshiny)
#https://www.youtube.com/watch?v=4T9tDX4aVS4
library(corrplot)
library(gplots)
library(RColorBrewer)
library(vegan)
library(rioja)

#choose Brewercolor
#display.brewer.all()
##si bp=boplot -> bp + scale_fill_brewer(palette = "Dark2")

#Data
ade4TkGUI()
decaALL
summary(decaALL)
v=length(decaALL)
# t(x) is creating a matrix...like everything in one line. So need to translate as a table
tdeca<-as.data.frame(t(decaALL))
class(tdeca)

summary (tdeca)
header(tdeca)
tdeca[0:1]



#correlation
cor(decaALL[3:16])
corrplot(cor(decaALL[3:16]),method="number",hclust.method="ward",main=levels(decaALL$Uniprot)[1])

#PCA
#if take deca§SHORT to avoid having the quali.sup line
decaSHORT<-as.data.frame(t(decaALL[2:4]))
res2<-PCA(decaALL,quali.sup=1,2)
res2<-PCA(decaSHORT)
plotellipses(res2)
plot(res2)
summary(res2,nbelements=Inf)


#fACTOSHINY-PCA
PCAshiny(decaSHORT)

#CAH
## For CAH, data need to be centered and scaled before being used for heatmap at least.
##trp.sc<-scale(trp ,center=T,scale=T);
##sc<-t(trp.sc);
## -> Done with log2FC, but need to be done for Data 

resALL<-PCA(decaALL,quali.sup=1)
HCPCshiny(decaALL)
# According to PCA, only 8 dimensions are sufficent to have the whole picture. SO I have to do the PCA again with only8 dimension
resALL<-PCA(decaALL,ncp=8,quali.sup=1)
ALL.hcpc<-HCPC(resALL,kk=Inf,min=3, max=10,consol=TRUE)
ALL.hcpc$call
ALL.hcpc$desc.var$quanti$"1"
ALL.hcpc$desc.var$quanti.var

#res.HCPC<-HCPC(decaALL,nb.clust=10,consol=FALSE,graph=FALSE,metric='euclidean')
#plot.HCPC(res.HCPC,choice='map',draw.tree=FALSE,title='Plan factoriel',axes=c(1,2))
#plot.HCPC(res.HCPC,choice='3D.map',ind.names=FALSE,centers.plot=FALSE,title='',angle=60,axes=c(1,2))
#plot.HCPC(res.HCPC,choice='tree',title='Arbre hiérarchique')

str(ALL.hcpc)
##uses euclidian and ward.D2 method

#heatmap


##colorLFK
mycol<-colorRampPalette(c("black","grey","white","darkgreen"))(100);
#mytopcolor<-c("gold1","gold1","gold1","green","green","green")

##heatmap - see later with CAH
#bou<-heatmap(as.matrix(decaALL[2:16]),col=mycol, xlab = "Ratio", ylab = "Probe")
#bou<-heatmap(as.matrix(hm[2:18]),col=mycol,ColSideColors=mytopcolor,density.info="none")

##clustering


#chclust/hclust : the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

##top cluster
mymat<-as.matrix(t(decaALL[2:16]));
eu.d<-vegdist(mymat,method="euclidean")
#reord<-chclust(eu.d,method="coniss")
reord<-hclust(eu.d,method="ward.D2");
top.clust<-as.dendrogram(reord)

##side cluster
mymat2<-as.matrix(decaALL[2:16]);
eu.d2<-vegdist(mymat2,method="euclidean");
#eu.d2<-vegdist(mymat2,method="manhattan"); #absolute value between 2 samples...
#eu.d2<-vegdist(mymat2,method="mahalanobis");#Mahalanobis distances are Euclidean distances of a matrix where columns are centred, have unit variance, and are uncorrelated
#reord2<-chclust(eu.d2,method="coniss");
reord2<-hclust(eu.d2,method="ward.D2");
#reord2<-hclust(eu.d2,method="centroid");
side.clust<-as.dendrogram(reord2,horiz=TRUE);



#hmp<-heatmap(as.matrix(decaALL[2:16]),col=mycol,margins=c(5,14),cexCol=1.0,cexRow=0.5,Colv=top.clust,Rowv=side.clust, xlab = "Ratio", ylab = "Probe");
#hmp<-heatmap(as.matrix(decaALL[2:16]),col=mycol,margins=c(5,14),cexCol=1.0,cexRow=0.5,Colv=top.clust,Rowv=side.clust)
hmp<-heatmap(as.matrix(decaALL[2:16]),col=mycol,cexCol=1.0,cexRow=0.5,Colv=top.clust,Rowv=side.clust)


#legend("right",col=mycol)
# Create a legend at (1, max_y) that is slightly smaller
# (cex) and uses the same line colors and points used by
# the actual plots
#legend(1, max_y, names(autos_data), cex=0.8, col=plot_colors, pch=21:23, lty=1:3);
#legend("topleft", c("Mon","Tue","Wed","Thu","Fri"), cex=0.6,bty="n", fill=rainbow(5));

#either put na as 0 for log2FC , either use:"na.color=par("bg")"

#recover the tree:
#hmpR<-as.hclust(hmp$rowDendrogram)
str(side.clust) # give the tree
str(side.clust[1:10])

#hc.rows <- hclust(dist(analyses_scaled), method="ward.D2")
#plot(hc.rows,sub="",xlab="",main="Clusters des individus")

#write.csv2(str(side.clust),"C:/Users/PCPORT-LIPOFABRIK-04/Dropbox (Lipofabrik)/MOLECULAR BIOLOGY/TRANSCRIPTOMICS/Results/Lipofabrik/Analyse Transcriptomique/2by2-comparisons-Table/Analysis R/dendro.Tree.csv")
# doesn't work

reord2$order->dendro.orderRow
write.csv2(dendro.orderRow,"C:/Users/PCPORT-LIPOFABRIK-04/Dropbox (Lipofabrik)/MOLECULAR BIOLOGY/TRANSCRIPTOMICS/Results/Lipofabrik/Analyse Transcriptomique/2by2-comparisons-Table/Analysis R/dendro.orderRow.csv")
reord$order->dendro.orderCol
write.csv2(dendro.orderCol,"C:/Users/PCPORT-LIPOFABRIK-04/Dropbox (Lipofabrik)/MOLECULAR BIOLOGY/TRANSCRIPTOMICS/Results/Lipofabrik/Analyse Transcriptomique/2by2-comparisons-Table/Analysis R/dendro.orderCol.csv")

##cut the heatmap, with k being the number of clusters where to cut and ==1, the first of this k cluster.
# Classe 1
heatmap(as.matrix(decaALL[2:16])[cutree(reord2,k=5)==4,], Colv=top.clust, scale='none', col=mycol,cexCol=1.0,cexRow=0.5)



#normal plot
plot.new
plot(decaALL$bLIP14,decaALL$M2.6h)
lines(decaALL$bLIP14,decaALL$M2.6h,col="forestgreen")
lines(decaALL$bLIP14,decaALL$M2.12h,col="blue")
lines(decaALL$bLIP14,decaALL$yya.f1.10h ,col="red")
#points(x) add points on an existing graph


# Start PNG device driver to save output to figure.png
#png(filename="C:/R/figure.png", height=295, width=300,bg="white")
# Graph autos using y axis that ranges from 0 to max_y.
# Turn off axes and annotations (axis labels) so we can
# specify them ourself
plot(autos_data$cars, type="o", col=plot_colors[1],
     ylim=c(0,max_y), axes=FALSE, ann=FALSE)
# Make x axis using Mon-Fri labels
axis(1, at=1:5, lab=c("Mon", "Tue", "Wed", "Thu", "Fri"))


# 201-2020 - BeraudM_Lipofabrik-Heatmap-Données Transcriptomiques

rm (list = ls())

library(ade4TkGUI)
library(FactoMineR)
library(gplots)
library(RColorBrewer)
library(vegan)
library(rioja)

# Ouvrir les dossiers:
#hm<-read.table("C:/Users/Utilisateur/Documents/R/MOOC_FUN_StatR/151105-Q-cm.csv",sep=";",row.names = 1,na.strings="NA",dec=".", strip.white=TRUE,fileEncoding="latin1",header=TRUE)
#hm2<-read.table("C:/Users/Utilisateur/Documents/R/MOOC_FUN_StatR/160215.csv",sep=";",row.names = 1,na.strings="NA",dec=".", strip.white=TRUE,fileEncoding="latin1",header=TRUE)
#hm<- read.table("C:\\Users\\User\\Documents\\UMONS\\heamap.txt", header=TRUE, sep="\t", row.names = 1,na.strings="NA", dec=".", strip.white=TRUE)


ade4TkGUI()

hm
summary(hm)

?heatmap

hm2<- hm[,1:12]
hm2
#Center-scale your data;For that, need to be in the other orientation
#for PSEPK = 116 prot
trp<-t(hm);
#trp<-hm
trp[2:16]
trp[2:4546]
trp.sc<-scale(trp ,center=T,scale=T);
sc<-t(trp.sc);
#sc<-trp.sc
#sc<-as.matrix(hm)


# color pattern with RColorBrewer
mytopcolor<-c("gold1","gold1","gold1","gray","gray","gray","deeppink4","deeppink4","deeppink4","darkcyan","darkcyan","darkcyan")
mytopcolor<-c("gold1","gold1","gold1","green","green","green")

mycol<-colorRampPalette(c("cyan","darkcyan","deeppink4","deeppink3"))(100);

#For PSEPK
#mysidecolor<-c("#0000FF","#1A11E4","#3522C9","#5034AE","#6B4593","#865678","#865678","#865678","#6B4593","#A1685D","#A1685D","#5034AE","gray","#6B4593","#F19C0D","gray","#BB7943","#6B4593","gray","#5034AE","gray","#5034AE","#FF5600","#1A11E4","#6B4593","#F19C0D","#5034AE","#5034AE","#1A11E4","gray","gray","#5034AE","#F19C0D","gray","#F19C0D","gray","#6B4593","gray","#F19C0D","#FF9C00","#6B4593","#865678","#BB7943","gray","gray","#6B4593","#6B4593","#FF8A00","#FF8A00","#BB7943","#FF3400","gray","#6B4593","#6B4593","gray","gray","#5034AE","#FF5600","#6B4593","gray","#1A11E4","#BB7943","#BB7943","#6B4593","#0000FF","gray","#6B4593","gray","#FF0000","#6B4593","#6B4593","gray","gray","gray","#0000FF","#FF5600","#6B4593","#BB7943","#5034AE","#5034AE","#FF8A00","#5034AE","gray","#FF5600","#F19C0D","#6B4593","#FF8A00","gray","#6B4593","#5034AE","#1A11E4","gray","gray","#5034AE","#5034AE","#5034AE","#5034AE","#5034AE","#5034AE","#5034AE","#5034AE","#5034AE","gray","gray","gray","#5034AE","gray","#5034AE","#6B4593","#1A11E4","#865678","#FF9C00","#5034AE","#FF0000","#FF0000")
#For CUPMC
#mysidecolor<-c("gray","#FF3900","gray","gray","white","white","#3F29BF","gray","white","gray","white","white","gray","white","gray","#F29C0C","white","white","white","gray","#FF4100","#FF7300","#FF7300","#FF7300","gray","#FF7300","#FF7300","#FF7300","#FF7300","#FF7300","#FF7300","#FF7300","#FF7300","#FF7300","#FF4100","#FF7300","white","white","#FF1000","white","white","white","#FF7300","#FF8C00","gray","white","white","white","white","gray","#E59419","white","gray","#FF7300","#FF4100","gray","#FF4100","white","#3F29BF","gray","#FF4100","gray","white","#2618D8","gray","#FF1000","white","gray","#FF8300","#FF5A00","white","#FF8300","#724A8C","#FF7300","white","white","gray","gray","white","white","gray","gray","gray","white","gray","white","#FF3100","gray","white","#FF2900","#FF6200","#FF6200","white","#3F29BF","gray","#FF7300","#FF7300","#1910E5","#FF2000","gray","gray","white","gray")
#mysidecolor<-c("#068488","#5A3764","white","#8B0A50","#5A3764","#8B0A50","white","#8B0A50","white","white","#5A3764","#841052","white","#5A3764","white")
mysidecolor<-Rcolor[,1]
side.col<-as.character(mysidecolor)
side.col

#ouvre de nouvelles fenêtre: x11(largeur, hauteur). largeur 2000 = un peu moins que la largeur de Zazie
x11(2000,12000)
dev.new()

bou<-heatmap(sc,col=mycol,ColSideColors=mytopcolor,RowSideColors = side.col,density.info="none",)
bou<-heatmap(sc,col=mycol,density.info="none",cexCol=1.0,cexRow=0.2)


#no side.col
bou<-heatmap(sc,col=mycol,ColSideColors=mytopcolor,density.info="none")
bou<-heatmap(sc,col=mycol,density.info="none")

bou<-heatmap(t(hm[2:16]),col=mycol,ColSideColors=mytopcolor,density.info="none")
bou<-heatmap(t(hm[2:16]),col=mycol,density.info="none")

summary (hm[2:16])


#colnames(sc)<-c("NoZ-t5","NoZ-t5","NoZ-t5","Z-t5","Z-t5","Z-t5","Z-t5","NoZ-t28","NoZ-t28","NoZ-t28","Z-t28","Z-t28","Z-t28","Z-t28")
#ColSideColors:(optional) character vector of length ncol(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x.
#RowSideColors:(optional) character vector of length nrow(x) containing the color names for a vertical side bar that may be used to annotate the rows of x.
#colRow, colCol	color of row/column labels, either a scalar to set the color of all labels the same, or a vector providing the colors of each label item
#labRow, labCol

#Prepare the color-key you like and value/color breaks if needed (colors=breaks-1);
##mycol<-colorRampPalette(c("dodgerblue3","dodgerblue","cornsilk2","orange","red"),space="rgb")(100);
##mysidecolor<-colorRampPalette(c("blue","orange","red"),space="rgb")(100)
##brk<-c(-3,-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3);
##col.brk<-colorRampPalette(c("dodgerblue3","cornsilk2","orange","red"),space="rgb")(18);
##heatmap.2(sc,col=mycol,margins=c(5,14),trace="none",density.info="none",lhei=c(2,8),lwid=c(2,8),cexCol=1.2,cexRow=1.6);

#Reordering according to euclidean matrices
#Re-order the samples according to entry design;
mymat<-as.matrix(trp);
eu.d<-vegdist(mymat,method="euclidean");
reord<-chclust(eu.d,method="coniss");
top.clust<-as.dendrogram(reord);
#plot(top.clust)

#not working for bou
mymat<-as.matrix(t(hm[2:16]);
eu.d <- vegdist(mymat, method = "euclidean")
                 
reord <- chclust(eu.d, method = "coniss")
top.clust <- as.dendrogram(reord)
    
mymat2 <- as.matrix(sc)
eu.d2 <- vegdist(mymat2, method = "euclidean")
                 
#eu.d2<-hclust(mymat2,method="euclidean")
reord2 <- chclust(eu.d2, method = "coniss")

side.clust <- as.dendrogram(reord2, horiz = TRUE)

#heatmap.2(trp,col=mycol,margins=c(5,14),trace="none",density.info="none",lhei=c(2,8),lwid=c(2,8),cexCol=0.2,cexRow=0.2);
#testvegdist<-euclidean, bray, gower,jaccard,canberra,manhattan,altGower,morisita,horn,mountford,raup,binomial, chao, cao,  mahalanobis
#print(eu.d2)
#x11()
#plot(side.clust)
str(side.clust) # tree
#sans top cluster
#heatmap.2(sc,col=mycol,margins=c(5,14),trace="none",density.info="none",lhei=c(2,8),lwid=c(2,8),cexCol=1.0,cexRow=0.2,Rowv=side.clust);

###Test matrix for proteins
#hclust<-"average" (= UPGMA),"ward.D", "ward.D2", "single", "complete", "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#average, complete, centroid, median
mymat2 <- as.matrix(sc)

eu.d2 <- vegdist(mymat2, method = "UPGMA")

reord2 <- hclust(eu.d2, method = "average")

side.clust <- as.dendrogram(reord2, horiz = TRUE)

heatmap.2(
  sc,
  col = mycol,
  margins = c(5, 14),
  trace = "none",
  density.info = "none",
  lhei = c(2, 8),
  lwid = c(2, 8),
  cexCol = 1.0,
  cexRow = 0.2,
  Colv = top.clust,
  Rowv = side.clust
)

heatmap.2(
  sc,
  col = mycol,
  margins = c(5, 14),
  trace = "none",
  density.info = "none",
  lhei = c(2, 8),
  lwid = c(2, 8),
  cexCol = 1.0,
  cexRow = 0.2,
  Rowv = side.clust
)


##dist<-"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
mymat2 <- as.matrix(sc)

eu.d2 <- dist(mymat2, method = "canberra")

reord2 <- hclust(eu.d2, method = "complete")

#reord2<-chclust(reord2,method="coniss");
side.clust <- as.dendrogram(reord2, horiz = TRUE)

heatmap.2(
  sc,
  col = mycol,
  margins = c(5, 14),
  trace = "none",
  density.info = "none",
  lhei = c(2, 8),
  lwid = c(2, 8),
  cexCol = 1.0,
  cexRow = 0.2,
  Colv = top.clust,
  Rowv = side.clust
)

#dist canberra + hclust complete/canberra + hclust binary



#heatmap.2(sc,col=mycol,margins=c(5,14),trace="none",density.info="none",lhei=c(2,8),lwid=c(2,8),cexCol=1.2,cexRow=1.6,Colv=top.clust);
#or
#heatmap.2(sc,margins=c(5,14),trace="none",density.info="none",lhei=c(2,8),lwid=c(2,8),cexCol=1.2,cexRow=1.6,Colv=top.clust,ColSideColors=mytopcolor);
heatmap.2(
  sc,
  col = mycol,
  margins = c(5, 14),
  trace = "none",
  density.info = "none",
  lhei = c(2, 8),
  lwid = c(2, 8),
  cexCol = 1.2,
  cexRow = 2.6,
  Colv = top.clust,
  RowSideColors = side.col,
  ColSideColors = mytopcolor
)

heatmap.2(
  sc,
  col = mycol,
  margins = c(5, 14),
  trace = "none",
  density.info = "none",
  lhei = c(2, 8),
  lwid = c(2, 8),
  cexCol = 1.2,
  cexRow = 2.6,
  Colv = top.clust,
  ColSideColors = mytopcolor
)

# with side.clust
heatmap.2(
  sc,
  col = mycol,
  margins = c(5, 14),
  trace = "none",
  density.info = "none",
  lhei = c(2, 8),
  lwid = c(2, 8),
  cexCol = 1.2,
  cexRow = 1.2,
  Colv = top.clust,
  Rowv = side.clust,
  ColSideColors = mytopcolor
)


heatmap.2(
  sc,
  col = mycol,
  margins = c(5, 14),
  trace = "none",
  density.info = "none",
  lhei = c(2, 8),
  lwid = c(2, 8),
  cexCol = 1.0,
  cexRow = 0.6,
  Colv = top.clust,
  ColSideColors = mytopcolor
)

heatmap.2(
  sc,
  col = mycol,
  margins = c(5, 14),
  trace = "none",
  density.info = "none",
  lhei = c(2, 8),
  lwid = c(2, 8),
  cexCol = 1.0,
  cexRow = 0.2,
  Colv = top.clust
)

heatmap.2(
  sc,
  col = mycol,
  margins = c(5, 14),
  trace = "none",
  density.info = "none",
  lhei = c(2, 8),
  lwid = c(2, 8),
  cexCol = 1.0,
  cexRow = 0.2,
  Colv = top.clust,
  Rowv = side.clust
)


#Add column and row side color key for samples and taxa;
##top.col<-as.character(mytopcolors$topcolors);
##side.col<-as.character(mysidecolors$sidecolors);
##heatmap.2(sc[,1:27],col=mycol,margins=c(5,14),trace="none",density.info="none",lhei=c(2,8),lwid=c(2,8),cexCol=1.2,cexRow=1.6,Colv=top.clust,RowSideColors=side.col,ColSideColors=top.col);


#creer une palette
mytopcolor <-
  colorRampPalette(c("blue", "orange", "red"), space = "rgb")(40)
mytopcolor <-
  colorRampPalette(c("darkcyan", "deeppink4"))(41)
mytopcolor <-
  colorRampPalette(c("blue", "green", "yellow", "orange", "red", "purple"))(12)
doss$COLOR <- as.character(doss$COLOR)
barplot(doss$VALUE, col = mytopcolor)#doesn't works ometimes because rgb not always here and might be problematic then use mytopcolor instead
mytopcolor

#PCA to select interesting cluster
resALL <- PCA(trp)
resALL$ind
#check the number of dimension necessary to have more than 70%, here = 5
resALL$eig
#redo the PCA with the right dimensions
resALL <- PCA(trp, ncp = 5)
#to have the contribution for each uniprot, to have them all
summary(resALL, nbelements = Inf)
write.csv2(
  summary(resALL, nbelements = Inf),
  "C:/Users/Utilisateur/Documents/R/MOOC_FUN_StatR/151105-Q-cmresults.csv"
)




hm2 <- heatmap.2(sc)
hc <- as.hclust(hm2$rowDendrogram)

#doesn't work
hm2 <- heatmap.2(sc, Rowv = side.clust)
vegdist(mymat2, method = "euclidean")
#Export the order of the taxa in the left dendrogram
hc$order -> dendro.order
write.csv2(
  dendro.order,
  "C:/Users/Utilisateur/Documents/R/MOOC_FUN_StatR/151105-Q-cmresults.csv"
)
#dendogram avec hclust = affiche un dendogram ?
#no

write.csv2(
  heatmap.2(sc, col = mycol)$rowDendogram,
  "C:/Users/Utilisateur/Documents/R/MOOC_FUN_StatR/Results.csv"
)
write.csv2(rownames(hm),
           "C:/Users/Utilisateur/Documents/R/MOOC_FUN_StatR/Results.csv")

cuttree(hm)
library(Factoshiny)
HCPC(hm)$data.clust
HCPC(hm)$call$t

# a tester

jpeg(filename = "C:/Users/Utilisateur/Documents/R/MOOC_FUN_StatR/hm2.jpg",
     width = 1000 ,
     height = 1000)
jpeg(filename = "hm2.jpg",
     width = 1000 ,
     height = 1000)

hm2
dev.off()
hm2$carpet #=> données inititales
bou <- t(hm2$carpet)
bou
as.dendrogram(reord2, horiz = TRUE) -> dendro.order


# Umons

# 2016_analyse efficacité qPCR en batch
library(ade4TkGUI)
ade4TkGUI()
library(qpcR)

data<-read.table("C:/Users/Utilisateur/Documents/R/MOOC_FUN_StatR/160215.csv",header= TRUE,sep=";",dec=".")
data


#fisher test entre samples

species=levels(data1$SPECIES)
# [1] "16281" "Bg"    "Bx"    "cm"    "Ec"    "Mv"    "Sf"    "SP902"
sample=levels(data1$SAMPLE)

tab2<-subset(data1,data1$SPECIES==species[2])
#tab21<-subset(tab2,tab2$SAMPLE=="B")
#tab21<-subset(tab2,tab2$SAMPLE==sample[1])
tab21<-rbind(subset(tab2,tab2$SAMPLE==sample[1]),subset(tab2,tab2$SAMPLE==sample[2]))


fisher.test(tab21$ASSAY,tab21$SAMPLE)
#chisq.test(tab21$ASSAY,tab21$SAMPLE,correct=FALSE)

# doesn't work, when you add several fluo, it maks only 1 efficiency for all...

tab<-table(c("Name","Eff.2","Rmsd"))
as.integer(bou)
as.integer(i)

bou=length(names(data))-1
i=2

while (i < bou+1){
  
  fit<-pcrfit(data, cyc = 1, fluo=i, model = l4)
  tab<-rbind(tab,cbind(names(data)[i],efficiency(fit, plot=FALSE)$eff,efficiency(fit, plot=FALSE)$Rsq))
  i=i+1
}

print (tab) 
write.csv2(file="C:/Users/Utilisateur/Documents/R/MOOC_FUN_StatR/results.csv",tab)



#fit<-pcrfit(data, cyc = 1, fluo=i, model = l4)
#tab<-rbind(tab,cbind(names(data)[i],efficiency(fit, plot=FALSE)$eff,efficiency(fit, plot=FALSE)$Rsq))




Recover column names
print(names(data))
length(names(data))
bou=length(names(data))-1
bou
names(data)[2]


fit<-pcrfit(data, cyc = 1, fluo=i, model = l4)

print(i,efficiency(fit, plot=TRUE)$eff,efficiency(fit, plot=TRUE)$Rsq )

efficiency(fit, plot=TRUE)$eff
efficiency(fit, plot=TRUE)$Rsq

fit
eff(fit)$effmax.y
efficiency(fit, plot=TRUE)$eff

remove(tab)

pcrGOF(fit,PRESS=FALSE)


