# Author Elisabeth Huellbusch

rm(list = ls())

##############################-
library(data.table)
library(raster)
library(readxl)
library(rworldmap)
library(sp)

library(RColorBrewer)
library(colorRamps)
library(maptools)
library(classInt)


##############################-
# load data ----

# go through states with a loop and process each state
# in the end do plotting
test <- read.csv("Data/Alabama_TREE_DATA.csv")

# read ECM file
myco <- read_xlsx("Data/ECM_and_references.xlsx", "Combined")
colnames(myco)[colnames(myco)=="Mycorrhizal status"] <- "Mycorrhizal_status"

# get data from worldclim
# TODO: specify lat and lon
r <- getData("worldclim",var="bio",res=10)

# read species reference file
spec <- read.csv("Data/REF_SPECIES.csv")

# for loop for reading tree files
treeFiles <- list.files("Data", pattern = "*TREE_DATA.csv")
treeFiles <- treeFiles[!grepl("part", treeFiles)]
trees <- list()
year <- 2014 # years to look at

now <- proc.time()
for(i in treeFiles){
  temp <- fread(paste0("Data/",i))
  temp <- temp[temp$INV_YR %in% year]
  state <- unlist(strsplit(i,"_"))[1]
  temp$STATE <- state
  eval(parse(text=paste0("trees$`",state, "` <- temp")))
  rm(temp)
}

trees <- rbindlist(trees)
trees$PLOT_ID <- paste(trees$STETE_CD, trees$PLOT, sep="_")

later <- proc.time()
later-now

#############################-
# check coordinates #########

newmap <- getMap(resolution = "low")

plot(newmap, xlim=c(min(trees$LON, na.rm=T), max(trees$LON, na.rm=T)), ylim=c(min(trees$LAT, na.rm=T), max(trees$LAT, na.rm=T)), asp=1)
points(trees$LON, trees$LAT, col=2, pch=16, cex=1)

trees$LAT[trees$LAT==0] <- NA
trees$LON[trees$LON==0] <- NA

plot(newmap, xlim=c(min(trees$LON, na.rm=T), max(trees$LON, na.rm=T)), ylim=c(min(trees$LAT, na.rm=T), max(trees$LAT, na.rm=T)), asp=1)
points(trees$LON, trees$LAT, col=2, pch=16, cex=1)


df1 <- data.frame(aggregate(LON ~ PLOT_ID, trees, mean, na.rm=T))
df2 <- data.frame(aggregate(LON ~ PLOT_ID, trees, sd, na.rm=T))
df3 <- data.frame(aggregate(LAT ~ PLOT_ID, trees, mean, na.rm=T))
df4 <- data.frame(aggregate(LAT ~ PLOT_ID, trees, sd, na.rm=T))

colnames(df1)[2] <- "LON.mean"
colnames(df2)[2] <- "LON.sd"
colnames(df3)[2] <- "LAT.mean"
colnames(df4)[2] <- "LAT.sd"

df <- merge(df1, df2, by="PLOT_ID")
df <- merge(df, df3, by="PLOT_ID")
df <- merge(df, df4, by="PLOT_ID")

head(df)
tail(df)

df <- df[order(df$LON.sd, decreasing = T),]


# in plots with NAs, any which are not NA?
bla <- unique(trees$PLOT_ID[is.na(trees$LAT)])

for(i in bla)
{
  checksum <- sum(abs(trees[trees$PLOT_ID==i, c("LON", "LAT")]))
  if(!is.na(checksum)) print(i)
}


#############################-
# mycorrhizal types #########
str(myco)
sort(unique(myco$Mycorrhizal_status))

myco$AM <- 0
myco$AM[setdiff(grep("AM",myco$Mycorrhizal_status), grep("AM-like", myco$Mycorrhizal_status))] <- 1

myco$NONE <- 0
myco$NONE[myco$Mycorrhizal_status=="None"] <- 1

# unique(myco[, c("ECM", "AM", "NONE")])

##############################-
# match species names #####
trees <- merge(trees, spec[,c("SPCD","GENUS","SPECIES")], by=c("SPCD"), all.x=T)

# match mycorrhizal data
trees <- merge(trees, myco, by.x=c("GENUS", "SPECIES"), by.y=c("Genus","Species"))

##############################-
# get coordinates ######

coords <- data.frame(PLOT_ID=unique(trees$PLOT_ID))
coords <- merge(coords, aggregate(LON ~ PLOT_ID, trees, mean, na.rm=T), by="PLOT_ID")
coords <- merge(coords, aggregate(LAT ~ PLOT_ID, trees, mean, na.rm=T), by="PLOT_ID")

# check if coords make sense
newmap <- getMap(resolution = "low")

plot(newmap, xlim=c(min(trees$LON, na.rm=T), max(trees$LON, na.rm=T)), ylim=c(min(trees$LAT, na.rm=T), max(trees$LAT, na.rm=T)), asp=1)
points(trees$LON, trees$LAT, col=2, cex=1, pch=16)
points(coords$LON, coords$LAT, col=3, cex=1, pch=16)

#############################-
# count how many species with ECM, AM, no mycorrhiza (NONE) #############

plots.number <- aggregate(SPCD ~ PLOT_ID, trees, function(x) {length(unique(x))})
colnames(plots.number)[2] <- "TOTAL"

temp <- aggregate(ECM ~ PLOT_ID+SPCD, trees, function(x){ifelse(sum(x)==0,0,1)})
plots.number <- merge(plots.number, aggregate(ECM ~ PLOT_ID, temp, sum), by="PLOT_ID")

temp <- aggregate(AM ~ PLOT_ID+SPCD, trees, function(x){ifelse(sum(x)==0,0,1)})
plots.number <- merge(plots.number, aggregate(AM ~ PLOT_ID, temp, sum), by="PLOT_ID")

temp <- aggregate(NONE ~ PLOT_ID+SPCD, trees, function(x){ifelse(sum(x)==0,0,1)})
plots.number <- merge(plots.number, aggregate(NONE ~ PLOT_ID, temp, sum), by="PLOT_ID")
rm(temp)

# calculate fractions
plots.number$ECM.frac <- plots.number$ECM/plots.number$TOTAL
plots.number$AM.frac <- plots.number$AM/plots.number$TOTAL
plots.number$NONE.frac <- plots.number$NONE/plots.number$TOTAL

head(plots.number)

# add coordinates
plots.number <- merge(plots.number, coords[, c("PLOT_ID", "LON", "LAT")], by="PLOT_ID")

head(plots.number);tail(plots.number)



#############################-
# aggregate by abundance #####

plots.abund <- aggregate(SPCD ~ PLOT_ID, trees, length)
colnames(plots.abund)[2] <- "TOTAL"

plots.abund <- merge(plots.abund, aggregate(ECM ~ PLOT_ID, trees, sum))
plots.abund <- merge(plots.abund, aggregate(AM ~ PLOT_ID, trees, sum))
plots.abund <- merge(plots.abund, aggregate(NONE ~ PLOT_ID, trees, sum))

plots.abund$ECM.frac <- plots.abund$ECM/plots.abund$TOTAL
plots.abund$AM.frac <- plots.abund$AM/plots.abund$TOTAL
plots.abund$NONE.frac <- plots.abund$NONE/plots.abund$TOTAL

plots.abund <- merge(plots.abund, coords, by="PLOT_ID")

head(plots.abund); tail(plots.abund)



#############################-
# aggregate by biomass #######

plots.bio <- aggregate(`ROOT(lb)` ~ PLOT_ID, trees, sum, na.rm=T)
colnames(plots.bio)[2] <- "TOTAL"

trees$ECM.bio <- trees$ECM * trees$`ROOT(lb)`
trees$AM.bio <- trees$AM * trees$`ROOT(lb)`
trees$NONE.bio <- trees$NONE * trees$`ROOT(lb)`

plots.bio <- merge(plots.bio, aggregate(ECM.bio ~ PLOT_ID, trees, sum))
plots.bio <- merge(plots.bio, aggregate(AM.bio ~ PLOT_ID, trees, sum))
plots.bio <- merge(plots.bio, aggregate(NONE.bio ~ PLOT_ID, trees, sum))

trees <- trees[, -(c("ECM.bio", "AM.bio", "NONE.bio"))]

plots.bio$ECM.frac <- with(plots.bio, ECM.bio/TOTAL)
plots.bio$AM.frac <- with(plots.bio, AM.bio/TOTAL)
plots.bio$NONE.frac <- with(plots.bio, NONE.bio/TOTAL)

plots.bio <- merge(plots.bio, coords, by="PLOT_ID")

head(plots.bio); tail(plots.bio)


############################################-
# plot ########

breaks <- 7
style <- "pretty"

todo <- c("plots.number", "plots.abund", "plots.bio")
todo2 <- c("ECM", "AM")
colors <- colorRampPalette(c("blue","red"))(breaks)

newmap <- getMap(resolution = "low")

xlim <- c(min(coords$LON), max(coords$LON))
ylim <- c(min(coords$LAT), max(coords$LAT))

x <- (par("usr"))[1]
y <- par("usr")[4]-(par("usr")[4]-par("usr")[3])*(1.05)

pdf(file=paste0("Figures/myco_dist_", min(year), "-", max(year), ".pdf") ,onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54) 

par(mfrow=c(length(todo2),1))

for(i in todo){
  temp <- get(i)
  
  for(type in todo2){
    
    brks <- classIntervals(temp[, paste0(type, ".frac")], n=breaks, style=style)
    colcode <- findColours(brks, colors)
    
    plot(newmap, xlim=xlim, ylim=ylim, asp=1, main=i)
    points(temp$LON, temp$LAT, col=colcode, pch=1, cex=0.5)
    
    legend("topleft", legend=c(names(attr(colcode, 'table'))), pch=1, col=c(attr(colcode, 'palette'))
           , title = paste("Fraction of Species with", type), horiz=F, cex=0.8, xpd=T, bg="white")
  }

}

par(mfrow=c(1,1))

dev.off()



###########################-
########## temperature

plot(r[[10]]/10)
newext <- extent(xlim, ylim)
map <- crop(r, newext)

plot(map[[1]]/10, asp=1)


# plot(r$bio1/10)
