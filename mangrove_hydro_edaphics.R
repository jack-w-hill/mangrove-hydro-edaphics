# ############################################################## ####
# Title: Limited relationships between mangrove forest structure ####
#        and hydro-edaphic conditions over a latitudinal         ####  
#        gradient in sub-tropical Queensland, Australia          ####  
# Author: Timothy L Staples                                      ####
# Collaborators: Jack W Hill                                     #### 
#                Lachlan A Bourke                                ####
#                Ciara M Horton                                  ####
#                Catherine E Lovelock                            ####
# ############################################################## ####

setwd("LOCATION OF THIS SCRIPT")

rm(list=ls())
library(lme4) # mixed-effect modelling
library(shape) # good arrows
library(multcomp) # linear hypothesis testing
library(mgcv) # additive modelling
library(lsmeans) # mean comparisons from models (post-hoc testing)
library(merTools) # predict from lmer/glmer models with confidence intervals

# mapping libraries
library(maptools)
library(raster)
library(rworldmap)

# PCQ functions
source("http://math.hws.edu/pcqm/pcqm.txt")

# read-in custom functions
sapply(paste0("./functions/", list.files("./functions")), source)

####### DATA IMPORT ####

mang <- read.csv("./data/Hill_et_al_2020-mangrove-survey-data.csv")

# convert points along transect to landward, transition and seaward categories
mang$transect.cat <- mang$point
levels(mang$transect.cat) <- list(land = c("P1", "P2"),
                                 trans = c("P3", "P4"),
                                 sea = c("P5", "P6"))

# correct missing point-level data, remove duplicate point rows
mang$index <- paste0(mang$site, mang$transect, mang$point)
mang.point <- mang[!duplicated(paste0(mang$site, mang$transect, mang$point)),]

# variables we want to keep for the point-level data
point.vars <- c("cover", "salinity", "soil_bulk_density", 
                "soil_water_content", "soil_carbon_content")

# merge across point-level data to replicate abiotic vars for each tree
mang <- mang[, !colnames(mang) %in% point.vars]
mang <- merge(mang, mang.point[,c("index", point.vars)], all.x=TRUE, all.y=FALSE)

# create transect-site index for future reference
mang$transect.index <- paste0(mang$site, ":", mang$transect)

#         Data correction ####

# calculation for basal area from dbh.
ba.calc <- function(dbh){pi*(dbh^2)/4}

# several data points were transcribed incorrectly. Fixing here.
mang[mang$index == "SouthStradbroke8P5" & mang$quarter == 4, "dbh"] = 3.5
mang[mang$index == "SouthStradbroke8P5" & mang$quarter == 4, "ba"] = ba.calc(3.5)

mang[mang$index == "SouthStradbroke9P3" & mang$quarter == 2, "dbh"] = 9.5
mang[mang$index == "SouthStradbroke9P3" & mang$quarter == 2, "ba"] = ba.calc(9.5)

#         PCQ attributes ####

# per-transect assessment of species importance

mang.import <- do.call("rbind", lapply(unique(mang$transect.index), function(x){
  
  sub.df <- mang[mang$transect.index == x,]
  sub.imp <- importance.val.local(sub.df[,c("point", "quarter", "species", 
                                      "distance_from_point", "dbh")])
  colnames(sub.imp) = c("rel.density", "rel.cover", "rel.freq", "import", "rel.import",
                        "abs.density")
  sub.imp$transect = sub.df$transect[1]
  sub.imp$site = sub.df$site[1]
  sub.imp$species <- rownames(sub.imp)
  
  return(sub.imp)
  }))

mang.dens <- do.call("rbind", lapply(unique(mang$transect.index), function(x){
  
  sub.df <- mang[mang$transect.index == x,]  
  sub.table <- wide.form(group.by = sub.df$point, 
                         spread.by = sub.df$quarter, 
                         values=sub.df$distance_from_point)
  sub.dens <- density.est(z=sub.table)
  sub.dens <- as.data.frame(t(unlist(sapply(sub.dens, as.numeric))[1:4]))
  rownames(sub.dens) = 1
  colnames(sub.dens) = c("n.points", "lwr", "upr", "density")
  sub.dens$transect = sub.df$transect[1]
  sub.dens$site = sub.df$site[1]
  
  return(sub.dens)
}))

boxplot(log(mang.dens$density) ~ mang.dens$site)

###### ABIOTIC PCA ####

# run PCA on abiotic variables, using correlation matrix (un-scaled)
mang.pca <- princomp(mang.point[,c("salinity", "soil_bulk_density", 
                                   "soil_water_content", "soil_carbon_content")],
                     cor=TRUE)

# variance explained by each PC
pca.var.prop <- mang.pca$sdev^2 / sum(mang.pca$sdev^2)

# k-means clustering to test for abiotic clusters

# set seed so cluster 1 always equals the same group for standardized plotting/colors
set.seed(250220) 
mang.kmeans <- lapply(2:50, function(n){kmeans(mang.pca$scores, centers=n)})
kmeans.ssprop <- sapply(mang.kmeans, function(x){x$betweenss / x$totss})
plot(kmeans.ssprop)
# we only need 4 groups to separate points with ~80% of post-PCA variation in
# abiotic conditions

# get order of sites by latitude for plotting
site.lat <- as.character(unique(mang.point$site[order(mang.point$lat, decreasing=TRUE)]))

# CONSTRUCT PCA DF ####

# copy PCA scores and cluster numbers to main data-frame
mang.pca.df <- data.frame(index = mang.point[,c("index")])
mang.pca.df$PC1 <- mang.pca$scores[,1]
mang.pca.df$PC2 <- mang.pca$scores[,2]
mang.pca.df$PC.clust2 <- mang.kmeans[[1]]$cluster
mang.pca.df$PC.clust3 <- mang.kmeans[[2]]$cluster
mang.pca.df$PC.clust4 <- mang.kmeans[[3]]$cluster

mang <- merge(mang, mang.pca.df,
              by.x="index", by.y="index",
              all.x=TRUE, all.y=FALSE)
mang$PC.clust2 <- as.factor(mang$PC.clust2)
mang$PC.clust3 <- as.factor(mang$PC.clust3)
mang$PC.clust4 <- as.factor(mang$PC.clust4)

# FIGURE 2 - OVERALL PCA PLOT ####

# plot of pca clusters
cand.cols <- col2rgb(c("blue", "darkgreen", "orange", "purple")) / 255
cand.pch <- c(21,22,24,25)
pdf("./plots/Figure 2 (overall PCA).pdf", height=3.5, width=3.5)
par(mar=c(0,0,0,0), oma=c(2,2,1,1), las=1, ps=8, tcl=-0.25, mgp=c(3,0.5,0))

plot(mang.pca$scores[,1:2],
     xlim=range(mang.pca$scores[,1]) * 1.1, 
     ylim=range(mang.pca$scores[,2]) * c(1.1, 1.2),
     axes=FALSE, type="n", xlab="", ylab="")

axis(side=2, at=seq(-2,2,2))
axis(side=2, at=seq(-3,3,1), tcl=-0.125, labels=NA)

mtext(side=2, line=1.25, text=paste0("PC2 (", round(pca.var.prop[2]*100,2), "%)"), las=0)
axis(side=1, mgp=c(3,0,0))
axis(side=1, mgp=c(3,0,0), at=seq(-5,5,1), labels=NA, tcl=-0.125)
mtext(side=1, line=0.75, text=paste0("PC1 (", round(pca.var.prop[1]*100,2), "%)"))

sapply(1:4, function(n){
  
  # create convex hull
  cluster.scores <- mang.pca$scores[mang.kmeans[[3]]$cluster == n, 1:2]
  temp.hull <- chull(cluster.scores)
  
  polygon(x=cluster.scores[temp.hull,1], y=cluster.scores[temp.hull,2],
          col=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.1),
          border=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.5))
  
  
  text(x=c(-1.35, 3.1, -4.1, -4.5)[n],
       y=c(-2.1, 1.75, 2.25, -0.5)[n],
       col=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1),
       font=2, labels=c("RI", "MD", "LD", "GS")[n], cex=1.15)
  
  site.scores <- mang.pca$scores[mang.kmeans[[3]]$cluster == n, 1:2]
  
  if(is.null(dim(site.scores))){
    points(y=site.scores[2], x=site.scores[1], 
           pch=cand.pch[n], cex=0.6, lwd=0.5,
           bg=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1),
           col=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1))
  } else {
    points(y=site.scores[,2], x=site.scores[,1], 
           pch=cand.pch[n], cex=0.6, lwd=0.5,
           bg=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1),
           col=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1))
  }
  
  
})

loadings.mult = 4
sapply(1:nrow(mang.pca$loadings), function(n){
  
  Arrows(x0=0, y0=0,
         x1 = mang.pca$loadings[n,1] * loadings.mult,
         y1 = mang.pca$loadings[n,2] * loadings.mult,
         lwd=1.5, arr.length=0.075, arr.width=0.1,
         col="black",
         arr.type="triangle", lend=2)
})

text(x = mang.pca$loadings[,1] * loadings.mult,
     y = mang.pca$loadings[,2] * loadings.mult,
     labels=c("Salinity", "Bulk density", "Water", "Carbon"), col="black",
     pos=c(4, 1, 2,3), offset=0.25, font=2)

box()
    dev.off()
    
# SUP. PLOT - PCA SCORES WITH VARYING CLUSTER NUMBER ####
# unused in published version

pdf("./plots/Overall PCA (varying clusters).pdf", height=3.5, width=3.5)
par(mar=c(0,0,0,0), oma=c(3,3,1,1), las=1, ps=8, tcl=-0.25, mgp=c(3,0.5,0))

sapply(2:10, function(n1){

plot(mang.pca$scores[,1:2],
     xlim=range(mang.pca$scores[,1]), 
     ylim=range(mang.pca$scores[,2]),
     axes=FALSE, type="n", xlab="", ylab="")

axis(side=2)
mtext(side=2, line=1.5, text=paste0("PC2 (", round(pca.var.prop[2]*100,2), "%)"), las=0)

axis(side=1, mgp=c(3,0,0))
mtext(side=1, line=1, text=paste0("PC1 (", round(pca.var.prop[1]*100,2), "%)"))

cand.cols <- col2rgb(c("blue","orange","darkgreen","purple", "red",
                       "goldenrod", "pink", "green", "skyblue", "grey70",
                       "grey20")) /255

sapply(1:n1, function(n){
    
    cluster.scores <- mang.pca$scores[mang.kmeans[[n1-1]]$cluster == n, 1:2]   
    
    # create convex hull of group
    temp.hull <- chull(cluster.scores)
    
    polygon(x=cluster.scores[temp.hull,1], y=cluster.scores[temp.hull,2],
            col=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.1),
            border=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.5))
    
    points(cluster.scores, pch=21, cex=0.6, lwd=0.5,
           bg=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.3),
           col=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0))
    
    text(x=mean(cluster.scores[temp.hull,1]), y=mean(cluster.scores[temp.hull,2]),
         col=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1),
         font=2, labels=n, cex=2)
    
})

sapply(1:nrow(mang.pca$loadings), function(n){
    
    Arrows(x0=0, y0=0,
           x1 = mang.pca$loadings[n,1] * 3,
           y1 = mang.pca$loadings[n,2] * 3,
           lwd=3, arr.length=0.075, arr.width=0.1,
           col="black",
           arr.type="triangle", lend=2)
    
    Arrows(x0=0, y0=0,
           x1 = mang.pca$loadings[n,1] * 3,
           y1 = mang.pca$loadings[n,2] * 3,
           lwd=1.5, arr.length=0.075, arr.width=0.1,
           col="black",
           arr.type="triangle", lend=2)
})

text(x=mang.pca$loadings[,1] * 3,
     y = mang.pca$loadings[,2] * 3 + c(0,0,-0.1,0.1),
     labels=c("Salinity", "Density", "Water", "Carbon"), font=2, col="black",
     pos=c(4, 4,2,2), offset=0.35)

mtext(side=3, line=0.1,
      text=paste0(n1, " Clusters (", round(kmeans.ssprop*100,2)[n1], "% variance)"),
      font=2)

box()    

})
dev.off()

# FIGURE 4 - MAP w/ SITE ABIOTIC CLUSTERS ####

# Australia outline shape file obtainable at www.gadm.org.
map <- readShapeSpatial("AUS SHAPEFILE LOCATION")
data(countriesCoarse)
map.low <- countriesCoarse
map.low <- map.low[map.low@data$SOVEREIGNT == "Australia",]

# Queensland local government area shape file obtainable at https://www.data.qld.gov.au/.
bris.lga <- readShapeSpatial("QLD LGA SHAPEFILE LOCATION")
bris.lga <- bris.lga[bris.lga@data$LGA == "Brisbane City",]
bris.lga <- crop(bris.lga, extent(152,153.2,-28,-27))
qld <- crop(map, extent(150,155,-32,-22.5))

site.mang <- mang[!duplicated(mang$site),]
site.mang$site <- as.character(site.mang$site)
site.mang <- site.mang[order(site.mang$lat, decreasing=TRUE),]
site.mang$site[site.mang$site=="SouthStradbroke"] = "Currigee"
site.mang$site[site.mang$site=="NorthStradbroke"] = "Minjerribah"
site.mang$site[site.mang$site=="Fraser"] = "K'gari"
site.mang$site[site.mang$site=="Bribie"] = "Boorabee"
site.mang$site[site.mang$site=="Moreton"] = "Moorgumpin"
site.mang$site[site.mang$site=="Noosa"] = "Dauwadhan"

pca.ylims <- seq(0.05,0.97,len=8)
pca.ylims <- cbind(pca.ylims[-length(pca.ylims)], pca.ylims[-1])

cand.cols <- col2rgb(c("blue", "darkgreen", "orange", "purple")) / 255
cand.pch <- c(21,22,24,25)

pdf("./plots/Figure 4 (map and PCA).pdf", height=7, width=6.5)

split.screen(rbind(c(0.10,0.7,0.05,0.99),
                   cbind(0.5,0.7,rev(pca.ylims[,1]), rev(pca.ylims[,2])),
                   c(0.10,0.30,0.845,0.99),
                   c(0.7,0.99,min(as.vector(pca.ylims)), 0.99),
                   c(0.10,0.99,0.05,0.99)))

screen(1)
par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0.1,0))
plot(qld, xlim=c(152.6,154.85), ylim=c(-27,-26), yaxs="i", xaxs="i",
     col="grey75", border=NA, axes=FALSE)
map.pars <- par("usr")
plot(bris.lga, border="black", col=rgb(0,0,0,0.2), add=TRUE, lwd=1)
points(site.mang$lat ~ site.mang$long, pch=16, cex=1)

par(lheight=0.8)
text(x=site.mang$long, 
     y=site.mang$lat, 
     labels=paste0(site.mang$site, " (",letters[-1][1:7],")"), pos=4, offset=0.4)
#font=ifelse(site.mang$island, 3, 1)
par(lheight=1)

text(x=153, y=-27.7, labels="Brisbane City", adj=0.5, font=2)

axis(side=1, at=seq(151,153.5,0.5), labels=parse(text=paste(seq(151,153.5,0.5), "*degree~E", sep="")))
axis(side=1, at=seq(151,153.5,0.1), tcl=-0.125, labels=NA)
axis(side=2, at=seq(-29,-25,0.5), mgp=c(3,0.5,0),
     labels=parse(text=paste(seq(-29,-25,0.5), "*degree~S", sep="")))
axis(side=2, at=seq(-29,-25,0.1), tcl=-0.125, labels=NA, mgp=c(3,0.5,0))

par(xpd=NA)
scalebar(d=50, xy=c(relative.axis.point(0.1, "x"),
                    relative.axis.point(0.05, "y")), 
         lonlat=TRUE, label=c("","50 km", ""), lwd=1)
par(xpd=FALSE)

text(x=relative.axis.point(0.36, "x"),
     y=relative.axis.point(0.9875, "y"),
     labels="(a)", font=2)
close.screen(1)

sapply(1:length(site.lat), function(n1){
  
  screen(n1+1)
  par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0.5,0))
  
  site <- site.lat[n1]
  
  site.bin <- mang.point$site == site
  plot(mang.pca$scores[site.bin,1:2],
       xlim=range(mang.pca$scores[,1]) * 1.1, 
       ylim=range(mang.pca$scores[,2]) * c(1.1, 1.2),
       axes=FALSE, type="n", xlab="", ylab="")
  
  axis(side=2, at=seq(-2,2,2))
  axis(side=2, at=seq(-3,3,1), tcl=-0.125, labels=NA)
  
  if(n1 == 4){mtext(side=2, line=1.25, 
                    text=paste0("PC2 (", round(pca.var.prop[2]*100,2), "%)"), las=0)}
  
  if(n1 == 7){
    axis(side=1, mgp=c(3,0,0))
    axis(side=1, mgp=c(3,0,0), at=seq(-5,5,1), labels=NA, tcl=-0.125)
    mtext(side=1, line=0.75, text=paste0("PC1 (", round(pca.var.prop[1]*100,2), "%)"))
  } else {axis(side=1, labels=NA)}
  
  sapply(1:nrow(mang.pca$loadings), function(n){
    
    Arrows(x0=0, y0=0,
           x1 = mang.pca$loadings[n,1] * 3,
           y1 = mang.pca$loadings[n,2] * 3,
           lwd=1.5, arr.length=0.075, arr.width=0.1,
           col="grey65",
           arr.type="triangle", lend=2)
  })
  
  text(x = mang.pca$loadings[,1] * 3,
       y = mang.pca$loadings[,2] * 3,
       labels=c("Salinity", "Density", "Water", "Carbon"), font=1, col="grey65",
       pos=c(4, 1, 2,3), offset=0.25, cex=0.8)
  
  sapply(1:4, function(n){
    
    # create convex hull
    cluster.scores <- mang.pca$scores[mang.kmeans[[3]]$cluster == n, 1:2]
    temp.hull <- chull(cluster.scores)
    
    polygon(x=cluster.scores[temp.hull,1], y=cluster.scores[temp.hull,2],
            col=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.1),
            border=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.5))
    
    
    text(x=c(-1.35, 3.1, -4.1, -4.5)[n],
         y=c(-2.1, 1.75, 2.25, -0.5)[n],
         col=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1),
         font=2, labels=c("RI", "MD", "LD", "GS")[n], cex=1.15)
    
    site.scores <- mang.pca$scores[mang.point$site == site &
                                     mang.kmeans[[3]]$cluster == n, 1:2]
    
    if(is.null(dim(site.scores))){
      points(y=site.scores[2], x=site.scores[1], 
             pch=cand.pch[n], cex=0.6, lwd=0.5,
             bg=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1),
             col="black")
    } else {
      points(y=site.scores[,2], x=site.scores[,1], 
             pch=cand.pch[n], cex=0.6, lwd=0.5,
             bg=rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1),
             col="black")
    }
    
    
  })
  
  box()
  text(x=relative.axis.point(0.0025, "x"),
       y=relative.axis.point(0.9, "y"),
       labels=paste0(" (", letters[-1][n1], ")"), font=2, adj=0)
  
  close.screen(n1+1)
})

screen(9)
par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0.1,0))
plot(x=NULL, y=NULL, xlim=c(112,155), ylim=c(-45,-10), axes=FALSE, xlab="", ylab="")
rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
     col="white", border=NA)
plot(map.low, add=TRUE, border="grey75", col="grey75")
rect(xleft=map.pars[1], xright=map.pars[2], ybottom=map.pars[3], ytop=map.pars[4],
     border="black", lwd=2, col=rgb(0,0,0,0.2))
box()
close.screen(9)

screen(10)
par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0,0))
mang.point <- merge(mang.point, mang.pca.df,
                    all.x=TRUE, all.y=FALSE)
mang.point$transect.index = paste0(mang.point$site,":",mang.point$transect)
mang.point <- mang.point[order(mang.point$lat),]

plot(x=NULL, y=NULL, xlim=c(0,max(mang.point$distance_along_tape)), 
     ylim=c(0.75, 75), axes=FALSE, yaxs="i", xlab="", ylab="")
box()

axis(side=1)

mtext(side=1, line=0.75, text="Distance along transect (m)")

site.tran <- table(mang.point$site[!duplicated(mang.point$transect.index)])
site.tran <- site.tran[site.lat]

site.centers <- seq(relative.axis.point(0.07, "y"), 
                    relative.axis.point(0.91, "y"), len=length(site.tran))
site.rad <- 4
tran.locs <- sapply(site.centers, function(x){ seq(x-site.rad,x + site.rad, len=10)})

segments(x0=0, x1=0, 
         y0=tran.locs[1,], 
         y1=tran.locs[nrow(tran.locs),], lwd=0.75)

sapply(1:length(unique(mang.point$transect.index)), function(n){
  
  tran.ind <- unique(mang.point$transect.index)[n]  
  sub.df <- mang.point[mang.point$transect.index == tran.ind,]
  
  segments(x0=0, x1=max(sub.df$distance_along_tape),
           y0=tran.locs[n], y1=tran.locs[n], lwd=0.75)
  
  points(y=rep(tran.locs[n], nrow(sub.df)), x=sub.df$distance_along_tape,
         pch=cand.pch[sub.df$PC.clust4], 
         bg=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1)[sub.df$PC.clust4],
         cex=0.8, lwd=0.5)
})

text(x=relative.axis.point(0.045, "x"),
     y=relative.axis.point(0.9875, "y"),
     labels="(i)", font=2)

par(lheight=0.8)
legend(x=150, y=76,
       legend=paste0(c("Regularly\ninundated", "Mineral\ndrained",
                       "Limited\ndrainage", "Groundwater\nsupplied"),
                     " (", c("RI", "MD", "LD", "GS"), ")"),
       pt.bg=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,]), 
       pch = cand.pch, bty="n",
       pt.lwd=0.5, y.intersp=1.3, x.intersp=0.7)
par(lheight=1)
rect(xleft=150, xright=par("usr")[2], ybottom=62.75, ytop=par("usr")[4])
close.screen(10)

screen(11)
par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0,0))
plot.new()
box()
close.screen(11)

close.screen(all.screens=TRUE)
dev.off()

# FIGURE 1 - MAP w/ SIMPLE BIOTIC CHARACTERS ####
pdf("./plots/Figure 1 (map and tree structure).pdf", height=7, width=6.5)

pca.ylims <- seq(0.05,0.97,len=8)
pca.ylims <- cbind(pca.ylims[-length(pca.ylims)], pca.ylims[-1])

split.screen(rbind(c(0.10,0.7,0.05,0.99),
                   c(0.10,0.30,0.845,0.99),
                   c(0.51,0.67,0.05,0.99),
                   c(0.67,0.83,0.05,0.99),
                   c(0.83,0.99,0.05,0.99),
                   c(0.10,0.99,0.05,0.99)))

screen(1)
par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0.1,0))
plot(qld, xlim=c(152.6,154.85), ylim=c(-27,-26), yaxs="i", xaxs="i",
     col="grey75", border=NA, axes=FALSE)
map.pars <- par("usr")
plot(bris.lga, border="black", col=rgb(0,0,0,0.2), add=TRUE, lwd=1)
points(site.mang$lat ~ site.mang$long, pch=16, cex=1)

par(lheight=0.8)
text(x=site.mang$long, 
     y=site.mang$lat, 
     labels=paste0(site.mang$site, " (",1:7,")"), pos=4, offset=0.4)
#font=ifelse(site.mang$island, 3, 1)
par(lheight=1)

text(x=153, y=-27.7, labels="Brisbane City", adj=0.5, font=2)

axis(side=1, at=seq(151,153.5,0.5), labels=parse(text=paste(seq(151,153.5,0.5), "*degree~E", sep="")))
axis(side=1, at=seq(151,153.5,0.1), tcl=-0.125, labels=NA)
axis(side=2, at=seq(-29,-25,0.5), mgp=c(3,0.5,0),
     labels=parse(text=paste(seq(-29,-25,0.5), "*degree~S", sep="")))
axis(side=2, at=seq(-29,-25,0.1), tcl=-0.125, labels=NA, mgp=c(3,0.5,0))

par(xpd=NA)
scalebar(d=50, xy=c(relative.axis.point(0.1, "x"),
                    relative.axis.point(0.05, "y")), 
         lonlat=TRUE, label=c("","50 km", ""), lwd=1)
par(xpd=FALSE)

text(x=relative.axis.point(0.36, "x"),
     y=relative.axis.point(0.9875, "y"),
     labels="(a)", font=2)
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0.1,0))
plot(x=NULL, y=NULL, xlim=c(112,155), ylim=c(-45,-10), axes=FALSE, xlab="", ylab="")
rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
     col="white", border=NA)
plot(map.low, add=TRUE, border="grey75", col="grey75")
rect(xleft=map.pars[1], xright=map.pars[2], ybottom=map.pars[3], ytop=map.pars[4],
     border="black", lwd=2, col=rgb(0,0,0,0.2))
box()
close.screen(2)

screen(3)
par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0.1,0))
plot(x=NULL, y=NULL, xlim=c(0,8.5), 
     ylim=c(0.75, 75), axes=FALSE, yaxs="i", xlab="", ylab="")
box()
axis(side=1)
mtext(side=1, line=0.75, text="Tree height (m)")

# height model
mang.height <- lmer(height ~ site*point + (1|transect.index), data=mang)
pred.df <- expand.grid(site=factor(levels(mang$site), levels=levels(mang$site)),
                       point=factor(levels(mang$point), levels=levels(mang$point)),
                       transect.index="a")
pred.df <- cbind(pred.df, predictInterval(mang.height,
                                          newdata = pred.df,
                                          n.sims = 3999,
                                          which="fixed",
                                          level=0.95,
                                          include.resid.var=FALSE,
                                          type="linear.prediction"))

site.tran <- table(mang.point$site[!duplicated(mang.point$transect.index)])
site.tran <- site.tran[site.lat]

site.centers <- seq(relative.axis.point(0.07, "y"), 
                    relative.axis.point(0.93, "y"), len=length(site.tran))
site.rad <- 4
point.locs <- sapply(site.centers, function(x){ seq(x-site.rad,x + site.rad, len=6)})
point.locs <- point.locs[,ncol(point.locs):1]

rect.divs <- rowMeans(cbind(c(par("usr")[4], point.locs[1,]),
                            c(point.locs[nrow(point.locs),], par("usr")[3])))
rect.divs[1] = par("usr")[4]
rect.divs[length(rect.divs)] = par("usr")[3]

rect(xleft=par("usr")[1], xright=par("usr")[2],
     ybottom=rect.divs[-1],
     ytop=rect.divs[-length(rect.divs)],
     border=NA, col=c("white","grey90"))

sapply(1:length(site.lat), function(n){
  site <- site.lat[n]
  site.preds <- pred.df[pred.df$site == site,]
  site.preds <- site.preds[order(site.preds$point),]

  arrows(x0=site.preds$lwr, x1=site.preds$upr,
         y0=point.locs[,n], y1=point.locs[,n],
         angle=90, length=0.025, code=3)  
  segments(x0=site.preds$fit[-nrow(site.preds)],
           x1=site.preds$fit[-1],
           y0=point.locs[-nrow(point.locs),n],
           y1=point.locs[-1,n], col="grey70")
  points(x = site.preds$fit, y = point.locs[,n], pch=21, bg="grey70")
  
  axis(side=2, at=range(point.locs[,n]), 
       labels=c("Land", "Sea"), mgp=c(3,0.5,0))
  par(xpd=NA)
  Arrows(x0=relative.axis.point(-0.175, "x"),
         x1=relative.axis.point(-0.175, "x"),
         y0=range(point.locs[,n])[1] + 0.75,
         y1=range(point.locs[,n])[2] - 1,
         arr.type="triangle", arr.length=0.1, arr.width=0.1)
  par(xpd=FALSE)


})
box()
close.screen(3)

screen(4)
par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0.1,0))
plot(x=NULL, y=NULL, xlim=c(0,16), 
     ylim=c(0.75, 75), axes=FALSE, yaxs="i", xlab="", ylab="")
box()
axis(side=1)
mtext(side=1, line=0.75, text="DBH (cm)")

rect(xleft=par("usr")[1], xright=par("usr")[2],
     ybottom=rect.divs[-1],
     ytop=rect.divs[-length(rect.divs)],
     border=NA, col=c("white","grey90"))

# height model
mang.dbh <- lmer(dbh ~ site*point + (1|transect.index), data=mang)
pred.df <- expand.grid(site=factor(levels(mang$site), levels=levels(mang$site)),
                       point=factor(levels(mang$point), levels=levels(mang$point)),
                       transect.index="a")
pred.df <- cbind(pred.df, predictInterval(mang.dbh,
                                          newdata = pred.df,
                                          n.sims = 3999,
                                          which="fixed",
                                          level=0.95,
                                          include.resid.var=FALSE,
                                          type="linear.prediction"))

sapply(1:length(site.lat), function(n){
  site <- site.lat[n]
  site.preds <- pred.df[pred.df$site == site,]
  site.preds <- site.preds[order(site.preds$point),]
  
  arrows(x0=site.preds$lwr, x1=site.preds$upr,
         y0=point.locs[,n], y1=point.locs[,n],
         angle=90, length=0.025, code=3)  
  segments(x0=site.preds$fit[-nrow(site.preds)],
           x1=site.preds$fit[-1],
           y0=point.locs[-nrow(point.locs),n],
           y1=point.locs[-1,n], col="grey70")
  points(x = site.preds$fit, y = point.locs[,n], pch=21, bg="grey70")
  
  axis(side=2, at=range(point.locs[,n]), labels=NA,
       mgp=c(3,0.5,0))
  
})

box()
close.screen(4)

screen(5)
par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0.1,0))
plot(x=NULL, y=NULL, xlim=c(0,150), 
     ylim=c(0.75, 75), axes=FALSE, yaxs="i", xlab="", ylab="")
box()
axis(side=1)
mtext(side=1, line=0.9, text=expression("BA (m"^2*" ha"^-1*")"))

rect(xleft=par("usr")[1], xright=par("usr")[2],
     ybottom=rect.divs[-1],
     ytop=rect.divs[-length(rect.divs)],
     border=NA, col=c("white","grey90"))

text(x=relative.axis.point(0.85, "x"),
     y=site.centers + c(0,-3,0,0,0,0,0), 
     labels=paste0("(",7:1,")"), font=2, 
     col=rgb(0.5,0.5,0.5,0.8), cex=1.5)

# height model
mang.ba <- lmer(log(ba) ~ site*point + (1|transect.index), data=mang)
pred.df <- expand.grid(site=factor(levels(mang$site), levels=levels(mang$site)),
                       point=factor(levels(mang$point), levels=levels(mang$point)),
                       transect.index="a")
pred.df <- cbind(pred.df, sapply(predictInterval(mang.ba,
                                          newdata = pred.df,
                                          n.sims = 3999,
                                          which="fixed",
                                          level=0.95,
                                          include.resid.var=FALSE,
                                          type="linear.prediction"), exp))

sapply(1:length(site.lat), function(n){
  site <- site.lat[n]
  site.preds <- pred.df[pred.df$site == site,]
  site.preds <- site.preds[order(site.preds$point),]
  
  arrows(x0=site.preds$lwr, x1=site.preds$upr,
         y0=point.locs[,n], y1=point.locs[,n],
         angle=90, length=0.025, code=3)  
  segments(x0=site.preds$fit[-nrow(site.preds)],
           x1=site.preds$fit[-1],
           y0=point.locs[-nrow(point.locs),n],
           y1=point.locs[-1,n], col="grey70")
  points(x = site.preds$fit, y = point.locs[,n], pch=21, bg="grey70")
  
  axis(side=2, at=range(point.locs[,n]), labels=NA,
       mgp=c(3,0.5,0))
  
})

box()
close.screen(5)

screen(6)
par(mar=c(0,0,0,0), las=1, ps=8, tcl=-0.25, mgp=c(3,0,0))
plot.new()
box()
close.screen(6)

close.screen(all.screens=TRUE)
dev.off()

###### BIOTIC PCA ####

mang$log.height <- log(mang$height)
mang$log.dbh <- log(mang$dbh)
mang$log.ba <- log(mang$ba)

cor(mang[,c("log.height", "log.dbh")])
# lots of shared information across the three biotic variables

# quick PCA of biotic variables
mang.bio.pca <- princomp(mang[,c("log.height", "log.dbh")],
                         cor = TRUE)
bio.pca <- as.data.frame(mang.bio.pca$scores)
colnames(bio.pca) = paste0("bioPC", 1:ncol(bio.pca))
mang <- cbind(mang, bio.pca)
biplot(mang.bio.pca)
bio.pca.scores <- mang.bio.pca$sdev^2 / sum(mang.bio.pca$sdev^2)

# variation loaded almost completly onto two axes
# axis 1 = TREE SIZE big trees vs small trees
# axis 2 = TREE SHAPE tall, thin trees vs short, broad trees

# MEAN EFFECT MODELS ####

# tree size model
bioPC1.m <- lmer(bioPC1 ~ -1 + PC.clust4 + (1|site/transect/point) + (1|species), data=mang)
bioPCA1.glht <- glht(bioPC1.m, linfct = mcp(PC.clust4 = "Tukey"))
summary(bioPCA1.glht)
bioPC1.coef <- summary(bioPC1.m)$coefficients

# tree height/width model
bioPC2.m <- lmer(bioPC2 ~ -1 + PC.clust4 + (1|site/transect/point) + (1|species), data=mang)
bioPCA2.glht <- glht(bioPC2.m, linfct = mcp(PC.clust4 = "Tukey"))
summary(bioPCA2.glht)
bioPC2.coef <- summary(bioPC2.m)$coefficients

plot(x=NULL, y=NULL, xlim=c(-2.5,1), ylim=c(-1.5,1))
points(bioPC2.coef[,1] ~ bioPC1.coef[,1], pch=cand.pch)
segments(x0 = bioPC1.coef[,1] - 1.96 * bioPC1.coef[,2],
         x1 = bioPC1.coef[,1] + 1.96 * bioPC1.coef[,2],
         y0 = bioPC2.coef[,1], y1 = bioPC2.coef[,1])
segments(y0 = bioPC2.coef[,1] - 1.96 * bioPC2.coef[,2],
         y1 = bioPC2.coef[,1] + 1.96 * bioPC2.coef[,2],
         x0 = bioPC1.coef[,1], x1 = bioPC1.coef[,1])

# not super useful, but there is a bit of info there. Too much variation in
# tree size/shape within each cluster. It may be that it's the distribution of tree
# size/shape that differs among clusters rather than mean differences

# SIZE DISTRIBUTION MODELS  ####
# DENSITY KERNEL ESTIMATES  ####

# test density functions

bioPC1.tran.means <- as.data.frame(tapply(mang$bioPC1, list(rownames(mang),
                                                            mang$PC.clust4), 
                                          mean))
bioPC2.tran.means <- as.data.frame(tapply(mang$bioPC2, list(rownames(mang),
                                                            mang$PC.clust4), 
                                          mean))

PC.clust.dens <- lapply(1:4, function(n){
  
  temp1 <- bioPC1.tran.means[,n]
  temp1 <-  temp1[!is.na(temp1)]
  temp2 <- bioPC2.tran.means[,n]
  temp2 <-  temp2[!is.na(temp2)]
  
  return(list(PC1 = density(temp1, na.rm=TRUE,
                            from=min(temp1), to=max(temp1)),
              PC2 = density(temp2, na.rm=TRUE,
                            from=min(temp2), to=max(temp2))))
  
})

# unused in published version
pdf("./plots/Overall biotic PCA density functions.pdf", height=3.5, width=8, useDingbats = FALSE)
par(mfrow=c(1,2), mar=c(0,2,0,0), oma=c(3,1,1,1), ps=8, tcl=-0.25, mgp=c(3,0,0), las=1)
plot(x=NULL, y=NULL, xlim=range(bioPC1.tran.means, na.rm=TRUE), ylim=c(0,0.45), 
     xlab="", ylab="", yaxs="i", axes=FALSE)
axis(side=2, mgp=c(3,0.5,0))

mtext(side=2, line=1.75, text="Density", las=0)
axis(side=1)

mtext(side=1, line=1, text=paste0("Tree size (PC1: ", sprintf("%.2f", 
                                                              bio.pca.scores[1]*100),"%)"))
sapply(1:4, function(n){
  
  polygon(x = c(min(PC.clust.dens[[n]]$PC1$x),
                 PC.clust.dens[[n]]$PC1$x,
                 max(PC.clust.dens[[n]]$PC1$x)),
           y = c(0, PC.clust.dens[[n]]$PC1$y, 0),
           border = rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1), 
           col = rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.1), 
           lwd=1.5)
  
})
box()
text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.96,"y"),
     labels="(a)", font=2)

plot(x=NULL, y=NULL, xlim=range(bioPC2.tran.means, na.rm=TRUE), ylim=c(0,0.85), 
     xlab="", ylab="", yaxs="i", axes=FALSE)
axis(side=2, mgp=c(3,0.5,0))

axis(side=1)
mtext(side=1, line=1, text=paste0("Tree shape (PC2: ", 
                                  sprintf("%.2f", bio.pca.scores[2]*100),"%)"))
sapply(1:4, function(n){
  
  polygon(x = c(min(PC.clust.dens[[n]]$PC2$x),
                PC.clust.dens[[n]]$PC2$x,
                max(PC.clust.dens[[n]]$PC2$x)),
          y = c(0, PC.clust.dens[[n]]$PC2$y, 0),
          border = rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1), 
          col = rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.1), 
          lwd=1.5)
})
box()
text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.96,"y"),
     labels="(b)", font=2)

dev.off()

# Figure 6 - SITE-LEVEL DENSITY KERNEL PLOT ####

pdf("./plots/Figure 6 (site-level density functions).pdf", height=7.25, width=4.5, useDingbats = FALSE)
par(mfcol=c(7,2),mar=c(0,2,0,0), oma=c(4.5,2,2,1), tcl=-0.25, mgp=c(3,0.2,0), las=1)

sapply(1:length(site.lat), function(n){
  
  s <- site.lat[n]
  plot(x=NULL, y=NULL, xlim=c(-3,4), ylim=c(0,0.7), yaxs="i", axes=FALSE,
       xlab="", ylab="")
  axis(side=2, at=c(0,0.2,0.4), mgp=c(3,0.5,0))
  abline(v=0, col="grey80", lty="31")
  
  if(n==length(site.lat)){
    axis(side=1)
    mtext(side=1, line=2.75, text=paste0("Tree size\nPC1 (", 
                                      round(bio.pca.scores[1]*100,2),"%)"),
          cex=0.8)
  } else {
    axis(side=1, labels=NA)
  }
  if(n==4){mtext(side=2, line=2.25, text="Density", las=0, cex=0.8)}
  
  if(n==1){
    mtext(side=3, at=relative.axis.point(0.1, "x"), 
          line=0.1, text="Small", cex=0.65, font=3)
    mtext(side=3, at=relative.axis.point(0.9, "x"), 
          line=0.1, text="Large", cex=0.65, font=3)
    par(xpd=NA)
    Arrows(x0=relative.axis.point(0.25, "x"),
           x1=relative.axis.point(0.75, "x"),
           y0=relative.axis.point(1.1, "y"),
           y1=relative.axis.point(1.1, "y"),
           arr.type="triangle", arr.width=0.2, arr.length=0.2, lwd=1.5)
   par(xpd=FALSE) 
  }
  
  box()
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(ifelse(n==1, 0.90, 0.875), "y"),
       labels = paste0("(", letters[seq(1,2*length(site.lat),2)][n], ") ", 
                       site.mang$site[n]), font=2, adj=0)
  
sapply(1:4, function(n){
  
    sub.tree <- mang[mang$site == s & mang$PC.clust4 ==n,]
    
    if(nrow(sub.tree) < 5){return(NA)}
    
    sub.PC1 <- density(sub.tree$bioPC1, from=min(sub.tree$bioPC1), to=max(sub.tree$bioPC1))  
    
    polygon(x = c(min(sub.PC1$x),sub.PC1$x, max(sub.PC1$x)),
            y = c(0, sub.PC1$y, 0),
            border = rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1), 
            col = rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.1), 
            lwd=c(1.5,1,1.5,1.5)[n],
            lty=c("solid", "solid", "31", "11")[n])

  })
  
})

sapply(1:length(site.lat), function(n){
  
  s <- site.lat[n]
  plot(x=NULL, y=NULL, xlim=c(-3,1.5), ylim=c(0,1.3), yaxs="i", axes=FALSE,
       xlab="", ylab="")
  axis(side=2, at=c(0,0.5,1), mgp=c(3,0.5,0))
  abline(v=0, col="grey80", lty="31")
  
  if(n==length(site.lat)){
    axis(side=1)
    mtext(side=1, line=2.75, text=paste0("Tree shape\nPC2 (", 
                                         round(bio.pca.scores[2]*100,2),"%)"),
          cex=0.8)
  } else {
    axis(side=1, labels=NA)
  }
  if(n==1){
    mtext(side=3, at=relative.axis.point(0.15, "x"), 
          line=0.1, text="Short, thick", cex=0.65, font=3)
    mtext(side=3, at=relative.axis.point(0.875, "x"), 
          line=0.1, text="Tall, thin", cex=0.65, font=3)
    par(xpd=NA)
    Arrows(x0=relative.axis.point(0.335, "x"),
           x1=relative.axis.point(0.71, "x"),
           y0=relative.axis.point(1.1, "y"),
           y1=relative.axis.point(1.1, "y"),
           arr.type="triangle", arr.width=0.2, arr.length=0.2, lwd=1.5)
    par(xpd=FALSE) 
  }
  
  box()
  
  sapply(1:4, function(n){
    
    sub.tree <- mang[mang$site == s & mang$PC.clust4 ==n,]
    
    if(nrow(sub.tree) < 5){return(NA)}
    
    sub.PC <- density(sub.tree$bioPC2, from=min(sub.tree$bioPC2), to=max(sub.tree$bioPC2))  
    
    polygon(x = c(min(sub.PC$x),sub.PC$x, max(sub.PC$x)),
            y = c(0, sub.PC$y, 0),
            border = rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 1), 
            col = rgb(cand.cols[1,n], cand.cols[2,n], cand.cols[3,n], 0.1), 
            lwd=c(1.5,1,1.5,1.5)[n],
            lty=c("solid", "solid", "31", "11")[n])
    })
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(ifelse(n==1, 0.90, 0.875), "y"),
       labels = paste0("(", letters[seq(2,2*length(site.lat),2)][n], ") ", 
                       site.mang$site[n]), font=2, adj=0)
 
})

dev.off()

###### SPECIES PROBABILITY MODELS ####

# model as multinomial, treating AM as reference level
# use mgcv as it has a simple implementation of multinomial regressions
library(mgcv)

# re-level factor so AM is the reference level
mang.multi <- mang
mang$species <- relevel(mang$species, ref="AM")

# remove AA that is only present twice in the dataset. We highlight this
# result in text.
mang.multi <- droplevels(mang.multi[mang$species != "AA",])
  
# mgcv multinomial models need numeric factors starting at "0" for reference level
mang.multi$species.mn <- as.numeric(mang.multi$species)-1

# with site-level random effects implemented as random effect splines
# rather than true mixed-effect models (gamm4 and gamm do not like "extended"
# distributions like multinomial)
mang.multi$dummy=1
sp.gam <- gam(list(species.mn ~ -1 + PC.clust4 + s(site, bs="re", by=dummy),
                   ~ -1 + PC.clust4  + s(site, bs="re", by=dummy), 
                   ~ -1 + PC.clust4  + s(site, bs="re", by=dummy),
                   ~ -1 + PC.clust4  + s(site, bs="re", by=dummy)),
              family=multinom(K=4), data=mang.multi, control=list(maxit=500))

# Simple logistic regression to independently estimate probability of occurrence
# of reference level (AM).
mang.multi$is.AM <- ifelse(mang.multi$species == "AM", 1, 0)
AM.point <- do.call("rbind", lapply(split(mang.multi, f=mang.multi$index), function(x){
  temp <- x[1,]
  temp$success = sum(x$is.AM)
  temp$failure = 4-temp$success
  return(temp)
  }))

am.glm <- glmer(cbind(success,failure) ~ PC.clust4 + (1|site), family=binomial,
                data=AM.point)

# Predict estimates and CIs for each PC cluster from both models
pred.df <- expand.grid(site=levels(mang$site),
                       PC.clust4 = factor(1:4, levels=levels(mang$PC.clust4)))
pred.df$dummy = 1

sp.pred <- predict(sp.gam, newdata=pred.df, se.fit=TRUE)
sp.pred$fit[sp.pred$fit < log(0.01)] = NA
sp.pred$upper <- sp.pred$fit + 1.96 * sp.pred$se.fit
sp.pred$lower <- sp.pred$fit - 1.96 * sp.pred$se.fit
sp.pred$raw.fit <- sp.pred$fit

am.pred <-  predictInterval(am.glm,
                            newdata = pred.df,
                            n.sims = 3999,
                            which="full",
                            level=0.95,
                            include.resid.var=FALSE,
                            type="linear.prediction")
am.pred <- sapply(am.pred, plogis)

# site-level predictions without random effects
pred.df$dummy=0
pred.df$site="a"
pred.df <- pred.df[!duplicated(pred.df$PC.clust4),]
site.pred <- predict(sp.gam, newdata=pred.df, se.fit=TRUE)
site.pred$fit[site.pred$fit < -10] = NA
site.pred$upper <- site.pred$fit + 1.96 * site.pred$se.fit
site.pred$lower <- site.pred$fit - 1.96 * site.pred$se.fit
site.pred$raw.fit <- site.pred$fit
# site.pred$upper <- exp(site.pred$fit + 1.96 * site.pred$se.fit)
# site.pred$lower <- exp(site.pred$fit - 1.96 * site.pred$se.fit)
# site.pred$raw.fit <- exp(site.pred$fit)
#site.pred$raw.fit[is.na(site.pred$raw.fit)] = 0

am.pred.site <-  predictInterval(am.glm,
                            newdata = pred.df,
                            n.sims = 3999,
                            which="fixed",
                            level=0.95,
                            include.resid.var=FALSE,
                            type="linear.prediction")
am.pred.site <- sapply(am.pred.site, plogis)

# Figure 5 - species occurence probability ####

pdf("./plots/Figure 5 (species occurence probability).pdf", height=3.5, width=7, useDingbats = FALSE)

split.screen(rbind(c(0.1,0.25,0.1,0.99),
                   c(0.35,0.99,0.1,0.99)))

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(0.5, 4.5), ylim=c(-9,1.5),
     xlab="", ylab="", axes=FALSE)

abline(h=log(1), lty="31", col="grey60", lwd=1.5)
text(x=relative.axis.point(0.2, "x"), y=log(1)+0.2,
     labels=expression("Equivalent to "*italic("Avicennia marina")),
     adj=0.5, col="black")

# other species results

sapply(1:ncol(sp.pred$raw.fit), function(sp.n){

  sp.pos <- c(1:5)[sp.n]
  
sapply(1:4, function(pos.n){

  pos.pos <- sp.pos + seq(-0.3,0.3,len=4)[pos.n]
  pos.pos.off <- pos.pos + seq(-0.06, 0.06, len=7)

  temp.upper <- sp.pred$upper[((pos.n-1)*7+1):(pos.n*7),sp.n]
  temp.lower <- sp.pred$lower[((pos.n-1)*7+1):(pos.n*7),sp.n]
  temp.fit <- sp.pred$raw.fit[((pos.n-1)*7+1):(pos.n*7),sp.n]
  
segments(x0=pos.pos.off, x1=pos.pos.off, 
         y0=temp.upper, y1=temp.lower,
         col=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1)[pos.n])

points(y=temp.fit, x=pos.pos.off, 
       col=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1)[pos.n],
       bg=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1)[pos.n],
       pch=c(16,15,17,25)[pos.n], cex=0.7)
})
  
})
rect(xleft=par("usr")[1], xright=par("usr")[2],
     ybottom=par("usr")[3], ytop=par("usr")[4],
     col=rgb(1,1,1,0.8), border=NA)

# site-level preds
sapply(1:ncol(site.pred$raw.fit), function(sp.n){
  
  sp.pos <- c(1:5)[sp.n]
  
  pos.pos <- sp.pos + seq(-0.3,0.3,len=4)
  
  temp.upper <- site.pred$upper[,sp.n]
  temp.lower <- site.pred$lower[,sp.n]
  temp.fit <- site.pred$raw.fit[,sp.n]
    
  arrows(x0=pos.pos, x1=pos.pos,
           y0=temp.upper, y1=temp.lower, lwd=2,
           col=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1),
         angle=90, length=0.05, code=3)
    
    points(y=temp.fit, x=pos.pos, 
           bg=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1),
           pch=cand.pch)
  })
  
axis(side=1, at=1:4, labels=c("Aegiceras\ncorniculatum",
                            "Bruguieria\n gymnohiza",
                            "Ceriops\naustralis",
                            "Rhizophora\nstylosa"), font=3, mgp=c(3,0.65,0))

axis(side=2, at=log(c(0.001,0.01,0.1,1, 5)), labels=c(0.001,0.01,0.1,1,5))
axis(side=2, at=log(c(seq(1e-4, 0.001, 1e-4),
                      seq(0.001,0.01,0.001),
                      seq(0.01,0.1,0.01),
                      seq(0.1,1,0.1),
                      seq(1,5,1))), tcl=-0.125, labels=NA)

mtext(side=2, line=1.75, text="Relative probability ratio", las=0)

par(lheight=0.8)
legend(x="bottomright",
       legend=paste0(c("Regularly\ninundated", "Mineral\ndrained", 
                       "Limited\ndrainage", "Groundwater\nsupplied"),
                     " (", c("RI", "MD", "LD", "GS"), ")"),
       pt.bg=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,]), 
       pch = cand.pch, bty="n",
       pt.lwd=0.5, y.intersp=1.3, x.intersp=0.7)
par(lheight=1)
rect(xleft=relative.axis.point(0.78,"x"),
     xright=par("usr")[2], 
     ybottom=par("usr")[3], 
     ytop=relative.axis.point(0.38,"y"))

text(x=relative.axis.point(0.0125, "x"),
     y=relative.axis.point(0.97, "y"),
     labels="(B)", font=2, adj=0)
box()
close.screen(2)

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(0.5,4.5), ylim=c(0,0.9),
     xlab="", ylab="", xaxt="n")

# AM results

sapply(1:4, function(pos.n){
    
    am.pos <- (1:4)[pos.n]
    pos.pos.off <- am.pos + seq(-0.25, 0.25, len=7)
    
    temp.upper <- am.pred[((pos.n-1)*7+1):(pos.n*7),2]
    temp.lower <- am.pred[((pos.n-1)*7+1):(pos.n*7),3]
    temp.fit <- am.pred[((pos.n-1)*7+1):(pos.n*7),1]
    
    segments(x0=pos.pos.off, x1=pos.pos.off, 
             y0=temp.upper, y1=temp.lower,
             col=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1)[pos.n])
    
    points(y=temp.fit, 
           x=pos.pos.off, 
           col=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1)[pos.n],
           bg=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1)[pos.n],
           pch=c(16,15,17,25)[pos.n], cex=0.7)
  })
  rect(xleft=par("usr")[1], xright=par("usr")[2],
     ybottom=par("usr")[3], ytop=par("usr")[4],
     col=rgb(1,1,1,0.8), border=NA)

arrows(x0=1:4, 
         x1=1:4, 
         y0=am.pred.site[,2], y1=am.pred.site[,3],
         col=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1), lwd=2,
       angle=90, length=0.05, code=3)

points(y=am.pred.site[,1], x=1:4, 
       bg=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,], 1),
       pch=cand.pch)

axis(side=1, at=relative.axis.point(0.5,"x"), labels="Avicennia\nmarina", font=3, mgp=c(3,0.65,0))
mtext(side=2, line=1.5, text="Occurrence probability", las=0)
text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.97, "y"),
     labels="(A)", font=2, adj=0)

close.screen(1)

dev.off()

# STACKED BARCHART ####
sp.order <- table(mang$species)
sp.order <- names(sp.order)[order(sp.order, decreasing=TRUE)]

sp.table <- table(mang$PC.clust4, mang$site, mang$species)
sp.clust.levels <- expand.grid(dimnames(sp.table)[[1]], dimnames(sp.table)[[3]])
sp.clust.levels <- paste0(sp.clust.levels[,1], ":", sp.clust.levels[,2])

sp.table <- do.call("rbind", lapply(1:dim(sp.table)[3], function(n){
  
  x <- sp.table[,,n]
  
  data.frame(count = as.vector(x),
             PC.clust4 = rep(rownames(x), ncol(x)),
             site = rep(colnames(x), each=nrow(x)),
             species=dimnames(sp.table)[[3]][n])
}))
sp.table$species.clust <- factor(paste0(sp.table$PC.clust4,":", sp.table$species),
                                 levels=sp.clust.levels)

sp.wide <- as.data.frame(wide.form(sp.table$species.clust, 
                                   sp.table$site, 
                                   sp.table$count, prune=FALSE, prop=FALSE))
sp.wide[is.na(sp.wide)]=0

sp.wide <- split(sp.wide, f=substr(rownames(sp.wide),1,1))

sp.col <- col2rgb(c("brown", "aquamarine4", "mediumorchid4", 
                    "royalblue4", "steelblue4", "seashell3"))/255

# SUP. PLOT - species relative abundance ####

pdf("./plots/Supplementary figure 1 (species relative abundance).pdf", height=3.5, width=7, useDingbats = FALSE)
par(mar=c(4,3,0.5,0.5), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(0, 40), ylim=c(0,1), yaxs="i",
     axes=FALSE, xlab="", ylab="")

x.pos <- seq(relative.axis.point(0.125, "x"),
             relative.axis.point(0.875, "x"), len=4)
bar.pos <-  seq(-3.5,3.5,len=ncol(sp.wide[[1]]))
width=0.5 * diff(bar.pos)[1]
y.line.buffer=0.002
x.line.buffer= 0.035

sapply(1:length(x.pos), function(n){
  x <- x.pos[n]
axis(side=1, at=x + bar.pos,
     labels = NA, las=0)
  par(xpd=NA)
text(x=x+bar.pos, y=relative.axis.point(-0.035, "y"),
     labels=site.mang$site, srt=45, adj=1, cex=0.85,
     col=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,])[n])
axis(side=1, line=2.75, at=range(x+bar.pos) - c(1.5,0), 
     labels=NA, tcl=0.25,
     col=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,])[n],
     col.ticks=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,])[n])
mtext(side=1, line=2.65, at=mean(range(x+bar.pos) - c(1.5,0)),
      text=c("Regularly inundated", "Mineral drained", "Limited drainage",
             "Groundwater supplied")[n],
      col=rgb(cand.cols[1,], cand.cols[2,], cand.cols[3,])[n])
par(xpd=FALSE)
})

axis(side=2)
mtext(side=2, line=1.75, text="Proportion", las=0)

lapply(1:length(sp.wide), function(n){

  data <- prop.table(as.matrix(sp.wide[[n]]),2)
  temp.clust = substr(rownames(data),1,1)
  temp.sp = substr(rownames(data),3,4)
    
  bar.pos <- x.pos[n] + bar.pos

  
  sapply(1:ncol(data), function(n1){
    
    raw.data <- data[match(sp.order, temp.sp),n1]
    temp.data <- cumsum(raw.data)
    temp.limits <- cbind(c(0, temp.data[-length(temp.data)]),
                         temp.data)
    
    
    rect(xleft=bar.pos[n1]-width + x.line.buffer, 
         xright=bar.pos[n1]+width - x.line.buffer,
         ybottom=c(0,temp.data[-length(temp.data)]) + y.line.buffer,
         ytop=temp.data - y.line.buffer, 
         border=rgb(sp.col[1,],sp.col[2,],sp.col[3,],1),
         col=rgb(sp.col[1,],sp.col[2,],sp.col[3,],0.25),
         lend=2, ljoin=1)
         
    sp.labs <- which(raw.data>0.05)
    
    text(x=bar.pos[n1],
         y=temp.limits[sp.labs,1] + 0.5 * raw.data[sp.labs],
         labels=sp.order[raw.data>0.05], cex=0.7, font=2,
         col=rgb(sp.col[1,],sp.col[2,],sp.col[3,],1)[which(raw.data>0.05)])
    
  })
  
  })
box()
dev.off()