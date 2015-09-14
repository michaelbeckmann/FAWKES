setwd("~/Dropbox/FAWKES/Reading Material/data")
#~/Dropbox/FAWKES/Reading Material/data/Europe_wide_analysis
#library("xlsx") # actually, this seems to be rubbish. extremely slow. its better to create csv files from xlsx in libreoffice
library(plyr)
require(sp)
require(rgdal)
require(maptools)
require(spatstat)
require (raster)
require(rgeos)
library(pryr)

### load all files and prepare data (step 1), 
## test
# OSM data extracted through overpass-turbo and converted into csv
slipway<-read.csv("OSM_data/slipway.csv", header=TRUE, sep=";" )
marina<-(read.csv("OSM_data/marina.csv", header=TRUE, sep=";" ))
watermills<-(read.csv("OSM_data/watermills.csv", header=TRUE, sep=";" ))
fishing<-read.csv("OSM_data/fishing.csv", header=TRUE, sep=";" )
water_works<-read.csv("OSM_data/water_works.csv", header=TRUE, sep=";" )

# the following section may be skipped and the fullycombined.csv can be loaded directly
{
# # load Waterbase extracted tables (coming from ms access) for water quality.
# # files still contain alle measurement timesteps
# EQR_Phytobenthos_G<-read.csv("Waterbase/EQR_Phytobenthos_G.csv", header=TRUE, sep=";" )
# EQR_Invertebrate<-read.csv("Waterbase/EQR_Invertebrate.csv", header=TRUE, sep=";" )
# EQR_Phytobenthos_E<-read.csv("Waterbase/EQR_Phytobenthos_E.csv", header=TRUE, sep=";" )
# Nutrients_Nitrate<-read.csv("Waterbase/Nutrients_Nitrate.csv", header=TRUE, sep=";" )
# Nutrients_Nitrite<-read.csv("Waterbase/Nutrients_Nitrite.csv", header=TRUE, sep=";" )
# Nutrients_PH<-read.csv("Waterbase/Nutrients_PH.csv", header=TRUE, sep=";" )
# Nutrients_Temperature<-read.csv("Waterbase/Nutrients_Temperature.csv", header=TRUE, sep=";" )
# Nutrients_TOC<-read.csv("Waterbase/Nutrients_TOC.csv", header=TRUE, sep=";" )
# Nutrients_Total_Nitrogen<-read.csv("Waterbase/Nutrients_Total_Nitrogen.csv", header=TRUE, sep=";" )
# Nutrients_Total_Phosphorous<-read.csv("Waterbase/Nutrients_Total_Phosphorous.csv", header=TRUE, sep=";" )
# 
# # Waterbase extracted tablefor station data.
# Waterbase_rivers_v14_Stations<-read.csv("Waterbase/Waterbase_rivers_v14_Stations.csv", header=TRUE, sep=";" )
# 
# # make dataframes for further handling
# Nutrients_Nitrate_MaxYear<-data.frame()
# Nutrients_Nitrite_MaxYear<-data.frame()
# Nutrients_PH_MaxYear<-data.frame()
# Nutrients_Temperature_MaxYear<-data.frame()
# Nutrients_TOC_MaxYear<-data.frame()
# Nutrients_Total_Nitrogen_MaxYear<-data.frame()
# Nutrients_Total_Phosphorous_MaxYear<-data.frame()
# 
# # subset all dataset to contain only one measurement per station (that of the last mentioned year)
# 
# for(i in 1:length(unique(Nutrients_Nitrate$WaterbaseID))){ 
#   df<-subset(Nutrients_Nitrate,Nutrients_Nitrate[,"WaterbaseID"]==(paste(unique(Nutrients_Nitrate$WaterbaseID)[i])))
#   subdf<-(df[df$Year==(max(df$Year)), ])
#   Nutrients_Nitrate_MaxYear<-rbind(Nutrients_Nitrate_MaxYear, subdf)
#   #print(i)
# } 
# 
# for(i in 1:length(unique(Nutrients_Nitrite$WaterbaseID))){ 
#   df<-subset(Nutrients_Nitrite,Nutrients_Nitrite[,"WaterbaseID"]==(paste(unique(Nutrients_Nitrite$WaterbaseID)[i])))
#   subdf<-(df[df$Year==(max(df$Year)), ])
#   Nutrients_Nitrite_MaxYear<-rbind(Nutrients_Nitrite_MaxYear, subdf)
#   #print(i)
# } 
# 
# for(i in 1:length(unique(Nutrients_PH$WaterbaseID))){ 
#   df<-subset(Nutrients_PH,Nutrients_PH[,"WaterbaseID"]==(paste(unique(Nutrients_PH$WaterbaseID)[i])))
#   subdf<-(df[df$Year==(max(df$Year)), ])
#   Nutrients_PH_MaxYear<-rbind(Nutrients_PH_MaxYear, subdf)
#   #print(i)
# } 
# 
# for(i in 1:length(unique(Nutrients_Temperature$WaterbaseID))){ 
#   df<-subset(Nutrients_Temperature,Nutrients_Temperature[,"WaterbaseID"]==(paste(unique(Nutrients_Temperature$WaterbaseID)[i])))
#   subdf<-(df[df$Year==(max(df$Year)), ])
#   Nutrients_Temperature_MaxYear<-rbind(Nutrients_Temperature_MaxYear, subdf)
#   #print(i)
# } 
# 
# for(i in 1:length(unique(Nutrients_TOC$WaterbaseID))){ 
#   df<-subset(Nutrients_TOC,Nutrients_TOC[,"WaterbaseID"]==(paste(unique(Nutrients_TOC$WaterbaseID)[i])))
#   subdf<-(df[df$Year==(max(df$Year)), ])
#   Nutrients_TOC_MaxYear<-rbind(Nutrients_TOC_MaxYear, subdf)
#   #print(i)
# } 
# 
# for(i in 1:length(unique(Nutrients_Total_Nitrogen$WaterbaseID))){ 
#   df<-subset(Nutrients_Total_Nitrogen,Nutrients_Total_Nitrogen[,"WaterbaseID"]==(paste(unique(Nutrients_Total_Nitrogen$WaterbaseID)[i])))
#   subdf<-(df[df$Year==(max(df$Year)), ])
#   Nutrients_Total_Nitrogen_MaxYear<-rbind(Nutrients_Total_Nitrogen_MaxYear, subdf)
#   #print(i)
# } 
# 
# for(i in 1:length(unique(Nutrients_Total_Phosphorous$WaterbaseID))){ 
#   df<-subset(Nutrients_Total_Phosphorous,Nutrients_Total_Phosphorous[,"WaterbaseID"]==(paste(unique(Nutrients_Total_Phosphorous$WaterbaseID)[i])))
#   subdf<-(df[df$Year==(max(df$Year)), ])
#   Nutrients_Total_Phosphorous_MaxYear<-rbind(Nutrients_Total_Phosphorous_MaxYear, subdf)
#   #print(i)
# } 
# 
# # create list of all sub-setted datasets
# datalist<-c("EQR_Phytobenthos_G","EQR_Phytobenthos_E", "EQR_Invertebrate", "Nutrients_Nitrate_MaxYear","Nutrients_Nitrite_MaxYear","Nutrients_PH_MaxYear","Nutrients_Temperature_MaxYear","Nutrients_TOC_MaxYear","Nutrients_Total_Nitrogen_MaxYear","Nutrients_Total_Phosphorous_MaxYear" )
# 
# # add meaningful column names 
# for (j in 1:length(datalist)){
#   tmp<-get(datalist[j])
#   colnames(tmp) <- paste(datalist[j], colnames(tmp), sep = "_")
#   colnames(tmp)[grep("WaterbaseID", colnames(tmp))]<- "WaterbaseID"
#   assign(datalist[j], tmp)
# }
# 
# # join everything using the waterbaseID, needs to be done in steps for memory reasons
# total1<-join_all(list(EQR_Phytobenthos_G,EQR_Phytobenthos_E, EQR_Invertebrate), by="WaterbaseID", type = 'full')
# total1<-cbind(total1[c(11,7,9,18,20,28,30)])
# total2<-join_all(list(Nutrients_Nitrate_MaxYear,Nutrients_Nitrite_MaxYear,Nutrients_PH_MaxYear,Nutrients_Temperature_MaxYear), by="WaterbaseID", type = 'full')
# total2<-cbind(total2[c(1,6,14,22,30)])
# total3<-join_all(list(Nutrients_TOC_MaxYear,Nutrients_Total_Nitrogen_MaxYear,Nutrients_Total_Phosphorous_MaxYear), by="WaterbaseID", type = 'full')
# total3<-cbind(total3[c(1,6,14,22)])
# total<-join_all(list(total1,total2,total3), by="WaterbaseID", type = 'full', match="first")
# fullycombined<-join_all(list(Waterbase_rivers_v14_Stations, total), by="WaterbaseID", type = 'full', match="first")
# 
# # change column names for Lat and Lon (needed later for the loop to work)
# colnames(fullycombined)[21:22]<-c("Lon","Lat")
# 
# # write combined waterbase dataframe as csv
# write.csv(fullycombined, file="Waterbase_maxyear_fullycombined.csv")
}

### prepare data (step 2), snap all points to rivers

# read prepared waterbase dataframe, remove all incomplety cases based on coordinates
waterbase_maxyear<-read.csv("Waterbase_maxyear_fullycombined.csv", header=TRUE )
waterbase_maxyear<-waterbase_maxyear[complete.cases(waterbase_maxyear[,21:22]),]

# lines (rivers)
# read lines (CCM2 river data), its a shapefile. Has been combined using ArcGIS from all sub-datasets available for Europe
lines_df<-readOGR("CCM2/","Riversegments_All")

# list of point dataframes as character. contains all point datasets we would like to analyse as ES proxies + the station data (i.e. water quality)
#points_df_list<-c("watermills","marina","fishing","water_works","slipway", "waterbase_maxyear")
points_df_list<-c("marina")

# outer loop start
for(k in 1:length(points_df_list)){

  # get data, keep as tmp and order points by Lat  
  tmp<-get(points_df_list[[k]])
  tmp<-tmp[order(tmp$Lat),,drop=FALSE]
  
  #elevate to spatialpointsdataframe
  coordinates(tmp) <- c("Lon", "Lat")      
  proj4string(tmp)=CRS("+proj=longlat +datum=WGS84")
  
  #correct projection
  tmp<-spTransform(tmp, CRS( "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs" ))
  
  #assign unsnapped points back to input data for comparison
  #assign(points_df_list[[k]], tmp)

  # split points into chunks of 100 or 50 and define splitpointslist
  split_points_df_list<-split(tmp, (0:nrow(tmp)%/%400))
  split_points_df_list_snapped<-split_points_df_list

  # inner loop start
  for(i in 1:length(split_points_df_list)){ 
    # create bounding box of points
    bb<-bbox(split_points_df_list[[i]])
    b_poly <- as(extent(bb), "SpatialPolygons")
    lines_df_clipped<-gIntersection(lines_df, b_poly, byid = T)
    
### ATTENTION: bb needs to be adjusted to include at least maxDist around it!
    
    # clip rivers according to bounding box
    #lines_df_clipped<-gClip(lines_df,bb)
    # snap points to rivers, maximum distance is 500 meters
    snapped_points<-snapPointsToLines(split_points_df_list[[i]],lines_df_clipped,maxDist=500)
    # merge snapped points back to the rest
    split_points_df_list_snapped[[i]]<-snapped_points
    # show progress
    print(paste("data-subset",i,"of",length(split_points_df_list),sep=" "))
    rm(lines_df_clipped,snapped_points)
    gc()
  } 
  
  # create snapped points dataframe
  tmp_snapped_merged<-do.call("rbind", split_points_df_list_snapped)
  #assign(paste(points_df_list[[k]],"_snapped",sep=""),tmp_snapped_merged)
  # write csv
  write.csv(tmp_snapped_merged, file=(paste(points_df_list[[k]],"_snapped.csv",sep="")))
  writeOGR(tmp_snapped_merged, dsn = ".",layer=paste(points_df_list[k],"_snapped",sep=""),driver = "ESRI Shapefile")
  rm(tmp_snapped_merged,split_points_df_list_snapped,tmp,split_points_df_list)
  gc()
} 

### analyse distribution of point features (ES proxies) in relation to water quality
# match points with water quality values from stations 
# at the moment this is a "quick and dirty" approach that will not take into account whether or not 
# a point is actually located at the same river as the station it is linked to.
# in the future we should explore cost path analysis or similar approache to do this more thoroughly

watermills_snapped<-readOGR(".","watermills_snapped")
marina_snapped<-readOGR(".","marina_snapped")
fishing_snapped<-readOGR(".","fishing_snapped")
water_works_snapped<-readOGR(".","water_works_snapped")
slipway_snapped<-readOGR(".","slipway_snapped")

waterbase_maxyear_snapped<-read.csv("waterbase_maxyear_snapped.csv", header=TRUE )
coordinates(waterbase_maxyear_snapped) <- c("X.1", "Y")      
#waterbase_maxyear_snapped<-readOGR(".","waterbase_maxyear_snapped")

#points_df_list_snapped<-c("watermills_snapped","marina_snapped","fishing_snapped","water_works_snapped","slipway_snapped")
points_df_list_snapped<-c("marina_snapped")

#water_variables_list<-c("EQR_Phytobenthos_G_DeterminandStatusClass","EQR_Phytobenthos_G_MeanValueEQR","EQR_Phytobenthos_E_DeterminandStatusClass", "EQR_Phytobenthos_E_MeanValueEQR","EQR_Invertebrate_DeterminandStatusClass","EQR_Invertebrate_MeanValueEQR","Nutrients_Nitrate_MaxYear_Mean","Nutrients_Nitrite_MaxYear_Mean","Nutrients_PH_MaxYear_Mean","Nutrients_Temperature_MaxYear_Mean","Nutrients_TOC_MaxYear_Mean","Nutrients_Total_Nitrogen_MaxYear_Mean","Nutrients_Total_Phosphorous_MaxYear_Mean")
water_variables_list<-c("EQR_Phytobenthos_G_MeanValueEQR", "EQR_Phytobenthos_E_MeanValueEQR","EQR_Invertebrate_MeanValueEQR","Nutrients_Nitrate_MaxYear_Mean","Nutrients_Nitrite_MaxYear_Mean","Nutrients_PH_MaxYear_Mean","Nutrients_Temperature_MaxYear_Mean","Nutrients_TOC_MaxYear_Mean","Nutrients_Total_Nitrogen_MaxYear_Mean","Nutrients_Total_Phosphorous_MaxYear_Mean")
#points_df_list_snapped<-c("slipway_snapped")


    for (l in 1 : length(points_df_list_snapped)){
      tmp<-get(points_df_list_snapped[[l]])
      closestSiteVec <- vector(mode = "numeric",length = nrow(tmp))
      minDistVec     <- vector(mode = "numeric",length = nrow(tmp))
    
        for (m in 1 : nrow(tmp))
        {
          distVec <- spDistsN1(waterbase_maxyear_snapped,tmp[m,],longlat = TRUE)
          minDistVec[m] <- min(distVec)
          closestSiteVec[m] <- which.min(distVec)
        }
      
      
      for (n in 1 : length(water_variables_list)){
      PointAssignStations <- as(as.data.frame(waterbase_maxyear_snapped)[closestSiteVec,][,paste(water_variables_list[n])],"numeric")
      tmp_FinalTable = data.frame(coordinates(tmp),tmp$Name,closestSiteVec,minDistVec,PointAssignStations)
      assign(paste(points_df_list_snapped[[l]],water_variables_list[n],"_matched",sep=""),tmp_FinalTable)
      write.csv(tmp_FinalTable, file=(paste(points_df_list_snapped[[l]],"_",water_variables_list[[n]] ,"_matched.csv",sep="")))
      
      png(filename=paste("~/UFZ/FAWKES/",points_df_list_snapped[[l]],"_",water_variables_list[[n]],sep=""))
      hist(as.data.frame(waterbase_maxyear_snapped)[,(paste(water_variables_list[n]))], breaks=200, main=paste(points_df_list_snapped[[l]],"_matched",sep=""))
      #hist(as.data.frame(waterbase_maxyear_snapped)[,(paste(water_variables_list[n]))], breaks=200,xlim=c(quantile(as.data.frame(waterbase_maxyear_snapped)[,(paste(water_variables_list[n]))],na.rm=TRUE)[[2]],quantile(as.data.frame(waterbase_maxyear_snapped)[,(paste(water_variables_list[n]))],na.rm=TRUE)[[4]]), main=paste(points_df_list_snapped[[l]],"_matched",sep=""))
      
      hist(tmp_FinalTable$PointAssignStations, xlim=c(0,30),add=TRUE,breaks=200,col="red")
      dev.off()
      }
    }

### next steps undertaken in ArcGis

### load data from ArcGIS outpu

## lets start with the waterbase all
## because ArcGIS is stupid and exports badly formatted txt files we have to do the following little dance

# read tables function
setAs("character", "num.with.commas", function(from) as.numeric(gsub(",", "", from) ) )
read.tables <- function(file.names, ...) {
    require(plyr)
  tmp<-read.table(file.names[1], header=TRUE, sep=";", quote="\"")
  classes <- sapply(tmp, class)
  classes["FacilityID"]<-"num.with.commas" # more columns will be needed to be changed
  classes["IncidentID"]<-"num.with.commas" # more columns will be needed to be changed
  classes["ObjectID"]<-"num.with.commas" # more columns will be needed to be changed
  classes["Total_Length"]<-"num.with.commas" # more columns will be needed to be changed
  ldply(file.names, function(fn) data.frame(Filename=fn, read.table(fn,header=TRUE, sep=";",quote="\"", colClasses=classes, ...)))
}

# read waterbase data

#define new class for importing numbers containing commas that will be automatically removed
setAs("character", "num.with.commas", function(from) as.numeric(gsub(",", "", from) ) )

tmp<-read.table("Europe_wide_analysis/waterbase_all.txt", header=TRUE, sep=";", quote="\"")
classes <- sapply(tmp, class)
classes["FID"]<-"num.with.commas" # more columns will be needed to be changed
waterbase_all<-read.table("Europe_wide_analysis/waterbase_all.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
### ATTENTION! # for some weired reason, we need to add 1 to the FIDs otherwise its a complete mismatch - *!*&%$ YOU ArcGIS!
waterbase_all$FID<-(waterbase_all$FID+1) 
names(waterbase_all)[3:57]<-names(waterbase_maxyear_snapped)[3:57]


# read non-Strahler based tables
setwd("Europe_wide_analysis/")
slipway_all<-read.tables(c("slipway_all.txt"))
colnames(slipway_all)[3]<-"FID"
waterworks_all<-read.tables(c("waterworks_all.txt"))
colnames(waterworks_all)[3]<-"FID"
watermills_all<-read.tables(c("watermills_all.txt"))
colnames(watermills_all)[3]<-"FID"
marinas_all<-read.tables(c("marinas_all.txt"))
colnames(marinas_all)[3]<-"FID"
fishing_all<-read.tables(c("fishing_all.txt"))
colnames(fishing_all)[3]<-"FID"
setwd("..")

# merge everything according to FID
slipway_waterbase_all<-merge(slipway_all,waterbase_all,by="FID")
watermills_waterbase_all<-merge(watermills_all,waterbase_all,by="FID")
marinas_waterbase_all<-merge(marinas_all,waterbase_all,by="FID")
fishing_waterbase_all<-merge(fishing_all,waterbase_all,by="FID")

# read Strahler based tables
setwd("Europe_wide_analysis/")
slipway_strahler<-read.tables(c("slipway_s1.txt","slipway_s2.txt","slipway_s3.txt","slipway_s4.txt","slipway_s5.txt","slipway_s6.txt"))
colnames(slipway_strahler)[3]<-"FID"
#waterworks_strahler<-read.tables(c("waterworks_s1.txt","waterworks_s2.txt","waterworks_s3.txt","waterworks_s4.txt","waterworks_s5.txt","waterworks_s6.txt"))
watermills_strahler<-read.tables(c("watermills_s1.txt","watermills_s2.txt","watermills_s3.txt","watermills_s4.txt","watermills_s5.txt","watermills_s6.txt"))
colnames(watermills_strahler)[3]<-"FID"
marinas_strahler<-read.tables(c("marinas_s1.txt","marinas_s2.txt","marinas_s3.txt","marinas_s4.txt","marinas_s5.txt","marinas_s6.txt"))
colnames(marinas_strahler)[3]<-"FID"
fishing_strahler<-read.tables(c("fishing_s1.txt","fishing_s2.txt","fishing_s3.txt","fishing_s4.txt","fishing_s5.txt","fishing_s6.txt"))
colnames(fishing_strahler)[3]<-"FID"
setwd("..")

# merge everything according to FID
slipway_waterbase_strahler<-merge(slipway_strahler,waterbase_all,by="FID")
watermills_waterbase_strahler<-merge(watermills_strahler,waterbase_all,by="FID")
marinas_waterbase_strahler<-merge(marinas_strahler,waterbase_all,by="FID")
fishing_waterbase_strahler<-merge(fishing_strahler,waterbase_all,by="FID")

# histograms for Strahler based data

strahler_list<-c("slipway_waterbase_strahler","watermills_waterbase_strahler", "marinas_waterbase_strahler","fishing_waterbase_strahler")
strahler_list_names<-c("OSM: Slipways","OSM: Watermills", "OSM: Marinas","OSM: Fishing")

for (l in 1 : length(strahler_list)){
  tmp<-get(strahler_list[[l]])
hist(tmp$Nutrients_Total_Nitrogen_MaxYear_Mean[tmp$Nutrients_Total_Nitrogen_MaxYear_Mean>=0][tmp$Total_Length<=10000], breaks=20, main=paste(strahler_list_names[l]))
abline(v=3, col="red")

mtext(sum(tmp$Nutrients_Total_Nitrogen_MaxYear_Mean>=0&tmp$Nutrients_Total_Nitrogen_MaxYear_Mean<=3&tmp$Total_Length<=10000), side=3,line=-1.5, at=par("usr")[1]+0.7*diff(par("usr")[1:2]), cex=1.2)
mtext(sum(tmp$Nutrients_Total_Nitrogen_MaxYear_Mean>=3&tmp$Total_Length<=10000), side=3,line=-1.5, at=par("usr")[1]+0.8*diff(par("usr")[1:2]), cex=1.2)
mtext("10000m")
}

for (l in 1 : length(strahler_list)){
  tmp<-get(strahler_list[[l]])
  hist(tmp$Nutrients_Total_Nitrogen_MaxYear_Mean[tmp$Nutrients_Total_Nitrogen_MaxYear_Mean>=0][tmp$Total_Length<=1000], breaks=20, main=paste(strahler_list_names[l]))
  abline(v=1, col="red")
  all<-length(tmp$Nutrients_Total_Nitrogen_MaxYear_Mean&tmp$Total_Length<=1000)
  percent_good<-round((sum(tmp$Nutrients_Total_Nitrogen_MaxYear_Mean>=0&tmp$Nutrients_Total_Nitrogen_MaxYear_Mean<=3&tmp$Total_Length<=1000))/all*100)
  percent_bad<-round((sum(tmp$Nutrients_Total_Nitrogen_MaxYear_Mean>=3&tmp$Total_Length<=1000))/all*100)
  mtext(sum(tmp$Nutrients_Total_Nitrogen_MaxYear_Mean>=0&tmp$Nutrients_Total_Nitrogen_MaxYear_Mean<=1&tmp$Total_Length<=1000), side=3,line=-1.5, at=par("usr")[1]+0.7*diff(par("usr")[1:2]), cex=1.2)
  mtext(sum(tmp$Nutrients_Total_Nitrogen_MaxYear_Mean>=1&tmp$Total_Length<=1000), side=3,line=-1.5, at=par("usr")[1]+0.8*diff(par("usr")[1:2]), cex=1.2)
  mtext(percent_good, side=3,line=-2.5, at=par("usr")[1]+0.7*diff(par("usr")[1:2]), cex=1.2)
  mtext(percent_bad, side=3,line=-2.5, at=par("usr")[1]+0.8*diff(par("usr")[1:2]), cex=1.2)
  
  mtext("1000m")
}




sum(waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean>=0&waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean<=3)
sum(waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean>=3)




hist(slipway_waterbase_all$EQR_Phytobenthos_G_MeanValueEQR[slipway_waterbase_all$EQR_Phytobenthos_G_MeanValueEQR>0][slipway_waterbase_all$Total_Length<10000], breaks=20)

par(mar = c(5, 4, 4, 4) + 0.3)

hist(slipway_waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean[slipway_waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean>0][slipway_waterbase_all$Total_Length<1000], breaks=20,xlim=c(0,8))
abline(v=3, col="red")
mtext(sum(slipway_waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean>0&slipway_waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean<3&slipway_waterbase_all$Total_Length<1000), line=-2, adj=0.1)
mtext(sum(slipway_waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean>3&slipway_waterbase_all$Total_Length<1000), line=-2, adj=0.2)
mtext("1000m")

par(new = TRUE)
hist(waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean[waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean>0], breaks=200, ylim=c(0,1200),axes=FALSE,xlim=c(0,8),border="grey", xlab="",ylab="",main="")
abline(v=3, col="red")
mtext(sum(waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean>0&waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean<3), line=-2, adj=0.7)
mtext(sum(waterbase_all$Nutrients_Total_Nitrogen_MaxYear_Mean>3), line=-2, adj=0.9)


### Resterampe
# 
# # promote the input lists to SpatialPointsDataFrames
# 
# fullycombined<-fullycombined[complete.cases(fullycombined[,21:22]),]
# 
# coordinates(fullycombined) <- c("Longitude", "Latitude")
# coordinates(watermills_half) <- c("Lon", "Lat")      
# proj4string(watermills_half)=CRS("+proj=longlat +datum=WGS84")
# 
# pts_buffer=gBuffer(watermills_half_ETRS,width=500,byid=T) # This is the points with a buffer 
# point_buffer_intersect_river=over(pts_buffer,CCM2_riversegments_clipped) 
# 
# snapped_watermills<-snapPointsToLines(watermills_half_ETRS,CCM2_riversegments_clipped,maxDist=500)
# snapped_watermills<-snapPointsToLines(watermills_half_ETRS,point_buffer_intersect_river,maxDist=500)
# 
# closestSiteVec <- vector(mode = "numeric",length = nrow(watermills))
# minDistVec     <- vector(mode = "numeric",length = nrow(watermills))
# 
# # Get the vector index of the temperature station closest to each field station.
# # Use the spDistsN1 function to compute the distance vector between each
# # field station site and all of the temperature stations. Then, find and
# # retain the actual temperature, and the index of the closest temperature
# # to each transect station.
# #
# # spDistsN1 usage: spDistsN1(pointList, pointToMatch, longlat)
# #
# # where:
# #         pointList   : List of candidate points.
# #         pointToMatch: Single point for which we seek the closest point in pointList.
# #         longlat     : TRUE  computes Great Circle distance in km,
# #                       FALSE computes Euclidean distance in units of input geographic coordinates
# #
# # We use Great Circle distance to increase distance calculation accuracy at high latitudes
# # See the discussion of distance units in the header portion of this file
# #
# # minDistVec stores distance from the closest temperature station to each density measurement point.
# # closestSiteVec stores the index of the closest temperature station to each density measurement point.
# #
# 
# for (i in 1 : nrow(watermills))
# {
#   distVec <- spDistsN1(fullycombined,watermills[i,],longlat = TRUE)
#   minDistVec[i] <- min(distVec)
#   closestSiteVec[i] <- which.min(distVec)
# }
# 
# #
# # Create the station Assignment table: merge the station point list with the transect point list
# # into a five-column table by merging the temperature point and transect point lists.
# #
# PointAssignStations <- as(fullycombined[closestSiteVec,]$Nutrients_Total_Nitrogen_MaxYear_Mean,"numeric")
# FinalTable = data.frame(coordinates(watermills),watermills$Name,closestSiteVec,minDistVec,PointAssignStations)
# #
# # Update the Temperature Assignment table column names 
# #
# names(FinalTable) <- c("Lon","Lat","watermills.name","CloseStationIndex","Distance","Nutrients_Total_Nitrogen_MaxYear_Mean")
# #
# # And, at the end, write the point assignment file.
# #
# # message("Write temperature/density assignment table to disk in .csv format")
# # write.csv(FinalTable,file="FinalTempAssignments.csv")
# # #
# 
# hist(fullycombined$Nutrients_Total_Nitrogen_MaxYear_Mean, xlim=c(0,30), breaks=200)
# hist(FinalTable$Nutrients_Total_Nitrogen_MaxYear_Mean, xlim=c(0,30),add=TRUE,breaks=200)
# 
# # total<-merge(Waterbase_rivers_v14_Stations, EQR_Invertebrate, by="WaterbaseID", suffixes=c('_Waterbase_rivers_v14_Stations', '_EQR_Invertebrate'))
# # total<-merge(total, EQR_Phytobenthos_E, by="WaterbaseID", suffixes=c('2', 'EQR_Phytobenthos_E'))
# # total<-merge(total, EQR_Phytobenthos_G, by="WaterbaseID", suffixes=c('3', 'EQR_Phytobenthos_G'))

