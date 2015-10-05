setwd("~/Dropbox/FAWKES/Reading Material/data")
# ~/Dropbox/FAWKES/Reading Material/data/Europe_wide_analysis
#library("xlsx") # actually, this seems to be rubbish. extremely slow. its better to create csv files from xlsx in libreoffice
library(plyr)
require(sp)
require(rgdal)
require(maptools)
require(spatstat)
require (raster)
require(rgeos)
library(pryr)

### ecoregion based data


### waterbase_s1_ecor$FID_waterb + 1 = fishing_s1$FacilityID
waterbase_maxyear_snapped<-read.csv("waterbase_maxyear_snapped.csv", header=TRUE )


#define new class for importing numbers containing commas that will be automatically removed
setAs("character", "num.with.commas", function(from) as.numeric(gsub(",", "", from) ) )

tmp<-read.table("Europe_wide_analysis/waterbase_all_ecor.txt", header=TRUE, sep=";", quote="\"")
classes <- sapply(tmp, class)
classes["FID_waterb"]<-"num.with.commas" # more columns will be needed to be changed
waterbase_all_ecor<-read.table("Europe_wide_analysis/waterbase_all_ecor.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
waterbase_s1_ecor<-read.table("Europe_wide_analysis/waterbase_s1_ecor.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
waterbase_s2_ecor<-read.table("Europe_wide_analysis/waterbase_s2_ecor.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
waterbase_s3_ecor<-read.table("Europe_wide_analysis/waterbase_s3_ecor.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
waterbase_s4_ecor<-read.table("Europe_wide_analysis/waterbase_s4_ecor.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
waterbase_s5_ecor<-read.table("Europe_wide_analysis/waterbase_s5_ecor.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
waterbase_s6_ecor<-read.table("Europe_wide_analysis/waterbase_s6_ecor.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)

### ATTENTION! # for some weired reason, we need to add 1 to the FID_waterb otherwise its a complete mismatch - *!*&%$ YOU ArcGIS!
waterbase_all_ecor$FID_waterb<-(waterbase_all_ecor$FID_waterb+1) 
names(waterbase_all_ecor)[4:58]<-names(waterbase_maxyear_snapped)[3:57]
waterbase_s1_ecor$FID_waterb<-(waterbase_s1_ecor$FID_waterb+1) 
names(waterbase_s1_ecor)[4:58]<-names(waterbase_maxyear_snapped)[3:57]
waterbase_s2_ecor$FID_waterb<-(waterbase_s2_ecor$FID_waterb+1) 
names(waterbase_s2_ecor)[4:58]<-names(waterbase_maxyear_snapped)[3:57]
waterbase_s3_ecor$FID_waterb<-(waterbase_s3_ecor$FID_waterb+1) 
names(waterbase_s3_ecor)[4:58]<-names(waterbase_maxyear_snapped)[3:57]
waterbase_s4_ecor$FID_waterb<-(waterbase_s4_ecor$FID_waterb+1) 
names(waterbase_s4_ecor)[4:58]<-names(waterbase_maxyear_snapped)[3:57]
waterbase_s5_ecor$FID_waterb<-(waterbase_s5_ecor$FID_waterb+1) 
names(waterbase_s5_ecor)[4:58]<-names(waterbase_maxyear_snapped)[3:57]
waterbase_s6_ecor$FID_waterb<-(waterbase_s6_ecor$FID_waterb+1) 
names(waterbase_s6_ecor)[4:58]<-names(waterbase_maxyear_snapped)[3:57]

waterbase_all_strahler_ecor<-rbind(waterbase_s1_ecor,waterbase_s2_ecor,waterbase_s3_ecor,waterbase_s4_ecor,waterbase_s5_ecor,waterbase_s6_ecor)


### read and merge slipway data


tmp<-read.table("Europe_wide_analysis/slipway_all.txt", header=TRUE, sep=";", quote="\"")
classes <- sapply(tmp, class)
classes["FacilityID"]<-"num.with.commas" # more columns will be needed to be changed
classes["Total_Length"]<-"num.with.commas" # more columns will be needed to be changed

slipway_all<-read.table("Europe_wide_analysis/slipway_all.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
slipway_s1<-read.table("Europe_wide_analysis/slipway_s1.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
slipway_s2<-read.table("Europe_wide_analysis/slipway_s2.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
slipway_s3<-read.table("Europe_wide_analysis/slipway_s3.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
slipway_s4<-read.table("Europe_wide_analysis/slipway_s4.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
slipway_s5<-read.table("Europe_wide_analysis/slipway_s5.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
slipway_s6<-read.table("Europe_wide_analysis/slipway_s6.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)

colnames(slipway_all)[2]<-"FID_waterb"
colnames(slipway_s1)[2]<-"FID_waterb"
colnames(slipway_s2)[2]<-"FID_waterb"
colnames(slipway_s3)[2]<-"FID_waterb"
colnames(slipway_s4)[2]<-"FID_waterb"
colnames(slipway_s5)[2]<-"FID_waterb" 
colnames(slipway_s6)[2]<-"FID_waterb"


# merge everything according to FID
slipway_waterbase_all_ecor<-merge(slipway_all,waterbase_all_ecor,by="FID_waterb")
slipway_waterbase_s1_ecor<-merge(slipway_s1,waterbase_s1_ecor,by="FID_waterb")
slipway_waterbase_s2_ecor<-merge(slipway_s2,waterbase_s2_ecor,by="FID_waterb")
slipway_waterbase_s3_ecor<-merge(slipway_s3,waterbase_s3_ecor,by="FID_waterb")
slipway_waterbase_s4_ecor<-merge(slipway_s4,waterbase_s4_ecor,by="FID_waterb")
slipway_waterbase_s5_ecor<-merge(slipway_s5,waterbase_s5_ecor,by="FID_waterb")
slipway_waterbase_s6_ecor<-merge(slipway_s6,waterbase_s6_ecor,by="FID_waterb")

slipway_waterbase_strahler_ecor<-rbind(slipway_waterbase_s1_ecor,slipway_waterbase_s2_ecor,slipway_waterbase_s3_ecor,slipway_waterbase_s4_ecor,slipway_waterbase_s5_ecor,slipway_waterbase_s6_ecor)



### read and merge fishing data


tmp<-read.table("Europe_wide_analysis/fishing_all.txt", header=TRUE, sep=";", quote="\"")
classes <- sapply(tmp, class)
classes["FacilityID"]<-"num.with.commas" # more columns will be needed to be changed
classes["Total_Length"]<-"num.with.commas" # more columns will be needed to be changed

fishing_all<-read.table("Europe_wide_analysis/fishing_all.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
fishing_s1<-read.table("Europe_wide_analysis/fishing_s1.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
fishing_s2<-read.table("Europe_wide_analysis/fishing_s2.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
fishing_s3<-read.table("Europe_wide_analysis/fishing_s3.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
fishing_s4<-read.table("Europe_wide_analysis/fishing_s4.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
fishing_s5<-read.table("Europe_wide_analysis/fishing_s5.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
fishing_s6<-read.table("Europe_wide_analysis/fishing_s6.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)

colnames(fishing_all)[2]<-"FID_waterb"
colnames(fishing_s1)[2]<-"FID_waterb"
colnames(fishing_s2)[2]<-"FID_waterb"
colnames(fishing_s3)[2]<-"FID_waterb"
colnames(fishing_s4)[2]<-"FID_waterb"
colnames(fishing_s5)[2]<-"FID_waterb" 
colnames(fishing_s6)[2]<-"FID_waterb"


# merge everything according to FID
fishing_waterbase_all_ecor<-merge(fishing_all,waterbase_all_ecor,by="FID_waterb")
fishing_waterbase_s1_ecor<-merge(fishing_s1,waterbase_s1_ecor,by="FID_waterb")
fishing_waterbase_s2_ecor<-merge(fishing_s2,waterbase_s2_ecor,by="FID_waterb")
fishing_waterbase_s3_ecor<-merge(fishing_s3,waterbase_s3_ecor,by="FID_waterb")
fishing_waterbase_s4_ecor<-merge(fishing_s4,waterbase_s4_ecor,by="FID_waterb")
fishing_waterbase_s5_ecor<-merge(fishing_s5,waterbase_s5_ecor,by="FID_waterb")
fishing_waterbase_s6_ecor<-merge(fishing_s6,waterbase_s6_ecor,by="FID_waterb")

fishing_waterbase_strahler_ecor<-rbind(fishing_waterbase_s1_ecor,fishing_waterbase_s2_ecor,fishing_waterbase_s3_ecor,fishing_waterbase_s4_ecor,fishing_waterbase_s5_ecor,fishing_waterbase_s6_ecor)


### read and merge watermills data


tmp<-read.table("Europe_wide_analysis/watermills_all.txt", header=TRUE, sep=";", quote="\"")
classes <- sapply(tmp, class)
classes["FacilityID"]<-"num.with.commas" # more columns will be needed to be changed
classes["Total_Length"]<-"num.with.commas" # more columns will be needed to be changed

watermills_all<-read.table("Europe_wide_analysis/watermills_all.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
watermills_s1<-read.table("Europe_wide_analysis/watermills_s1.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
watermills_s2<-read.table("Europe_wide_analysis/watermills_s2.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
watermills_s3<-read.table("Europe_wide_analysis/watermills_s3.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
watermills_s4<-read.table("Europe_wide_analysis/watermills_s4.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
watermills_s5<-read.table("Europe_wide_analysis/watermills_s5.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
watermills_s6<-read.table("Europe_wide_analysis/watermills_s6.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)

colnames(watermills_all)[2]<-"FID_waterb"
colnames(watermills_s1)[2]<-"FID_waterb"
colnames(watermills_s2)[2]<-"FID_waterb"
colnames(watermills_s3)[2]<-"FID_waterb"
colnames(watermills_s4)[2]<-"FID_waterb"
colnames(watermills_s5)[2]<-"FID_waterb" 
colnames(watermills_s6)[2]<-"FID_waterb"


# merge everything according to FID
watermills_waterbase_all_ecor<-merge(watermills_all,waterbase_all_ecor,by="FID_waterb")
watermills_waterbase_s1_ecor<-merge(watermills_s1,waterbase_s1_ecor,by="FID_waterb")
watermills_waterbase_s2_ecor<-merge(watermills_s2,waterbase_s2_ecor,by="FID_waterb")
watermills_waterbase_s3_ecor<-merge(watermills_s3,waterbase_s3_ecor,by="FID_waterb")
watermills_waterbase_s4_ecor<-merge(watermills_s4,waterbase_s4_ecor,by="FID_waterb")
watermills_waterbase_s5_ecor<-merge(watermills_s5,waterbase_s5_ecor,by="FID_waterb")
watermills_waterbase_s6_ecor<-merge(watermills_s6,waterbase_s6_ecor,by="FID_waterb")

watermills_waterbase_strahler_ecor<-rbind(watermills_waterbase_s1_ecor,watermills_waterbase_s2_ecor,watermills_waterbase_s3_ecor,watermills_waterbase_s4_ecor,watermills_waterbase_s5_ecor,watermills_waterbase_s6_ecor)


### read and merge marinas data


tmp<-read.table("Europe_wide_analysis/marinas_all.txt", header=TRUE, sep=";", quote="\"")
classes <- sapply(tmp, class)
classes["FacilityID"]<-"num.with.commas" # more columns will be needed to be changed
classes["Total_Length"]<-"num.with.commas" # more columns will be needed to be changed

marinas_all<-read.table("Europe_wide_analysis/marinas_all.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
marinas_s1<-read.table("Europe_wide_analysis/marinas_s1.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
marinas_s2<-read.table("Europe_wide_analysis/marinas_s2.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
marinas_s3<-read.table("Europe_wide_analysis/marinas_s3.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
marinas_s4<-read.table("Europe_wide_analysis/marinas_s4.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
marinas_s5<-read.table("Europe_wide_analysis/marinas_s5.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)
marinas_s6<-read.table("Europe_wide_analysis/marinas_s6.txt", header=TRUE, sep=";", quote="\"", colClasses=classes)

colnames(marinas_all)[2]<-"FID_waterb"
colnames(marinas_s1)[2]<-"FID_waterb"
colnames(marinas_s2)[2]<-"FID_waterb"
colnames(marinas_s3)[2]<-"FID_waterb"
colnames(marinas_s4)[2]<-"FID_waterb"
colnames(marinas_s5)[2]<-"FID_waterb" 
colnames(marinas_s6)[2]<-"FID_waterb"


# merge everything according to FID
marinas_waterbase_all_ecor<-merge(marinas_all,waterbase_all_ecor,by="FID_waterb")
marinas_waterbase_s1_ecor<-merge(marinas_s1,waterbase_s1_ecor,by="FID_waterb")
marinas_waterbase_s2_ecor<-merge(marinas_s2,waterbase_s2_ecor,by="FID_waterb")
marinas_waterbase_s3_ecor<-merge(marinas_s3,waterbase_s3_ecor,by="FID_waterb")
marinas_waterbase_s4_ecor<-merge(marinas_s4,waterbase_s4_ecor,by="FID_waterb")
marinas_waterbase_s5_ecor<-merge(marinas_s5,waterbase_s5_ecor,by="FID_waterb")
marinas_waterbase_s6_ecor<-merge(marinas_s6,waterbase_s6_ecor,by="FID_waterb")

marinas_waterbase_strahler_ecor<-rbind(marinas_waterbase_s1_ecor,marinas_waterbase_s2_ecor,marinas_waterbase_s3_ecor,marinas_waterbase_s4_ecor,marinas_waterbase_s5_ecor,marinas_waterbase_s6_ecor)


### Mann Whitney U test

ecoregions<-unique(waterbase_all_ecor$NAME)
ecoregions<-ecoregions[c(1,16,17,8,9,14,15)]
#ecoregions<-ecoregions[17] # just for GB
strahler_ecor_list<-c("slipway_waterbase_strahler_ecor", "marinas_waterbase_strahler_ecor")#,"watermills_waterbase_strahler_ecor","fishing_waterbase_strahler_ecor")
variable_list<-c("EQR_Phytobenthos_G_MeanValueEQR","EQR_Phytobenthos_E_MeanValueEQR","EQR_Invertebrate_MeanValueEQR","Nutrients_Total_Nitrogen_MaxYear_Mean","Nutrients_Total_Phosphorous_MaxYear_Mean")
all_ecor_list<-c("slipway_waterbase_all_ecor", "marinas_waterbase_all_ecor")#,"watermills_waterbase_all_ecor","fishing_waterbase_all_ecor")

variable_names<-c("EQR_Phytobenthos_G","EQR_Phytobenthos_E","EQR_Invertebrate","Total_Nitrogen","Total_Phosphorous")
all_ecor_names<-c("OSM: Slipways", "OSM: Marinas", "OSM: Fishing")#,"watermills_waterbase_all_ecor","fishing_waterbase_all_ecor")


### Europe without starhler

for (x in 1: length(all_ecor_list) ){
  tmp_poi<-all_ecor_list[x]
  tmp_poi_data<-get(tmp_poi)
  tmp_poi_data<-tmp_poi_data[tmp_poi_data$Total_Length<=1000,]
  
  for (y in 1: length(variable_list) ){
    tmp_var<-variable_list[y] 
    
    waterbase_sub<-subset(waterbase_all_ecor)
    tmp_poi_sub<-subset(tmp_poi_data)
    
    for (z in 1: nrow(tmp_poi_sub)){
      waterbase_sub<-waterbase_sub[waterbase_sub$WaterbaseID!= paste(tmp_poi_sub$WaterbaseID[z]),]
    }
    
    if (nrow(tmp_poi_sub)>10){
      tmp_poi_sub_tmp_var<-tmp_poi_sub[,tmp_var]
      tmp_poi_sub_tmp_var<-data.frame(tmp_poi_sub_tmp_var,rep("tmp_poi"))
      names(tmp_poi_sub_tmp_var)<-c("value","group")
      
      waterbase_sub_tmp_var<-waterbase_sub[,tmp_var]
      waterbase_sub_tmp_var<-data.frame(waterbase_sub_tmp_var,rep("all"))
      names(waterbase_sub_tmp_var)<-c("value","group")
      
      tmp_poi_all_sub_tmp_var<-rbind.data.frame(tmp_poi_sub_tmp_var,waterbase_sub_tmp_var)
     
      png(paste("/tmp/","EU without Strahler",all_ecor_names[x],variable_names[y],".png",sep="_"),width = 1760, height = 960, units = "px",pointsize=40)
      par(mfrow = c(1, 4), pty = "s")
      hist(as.numeric(tmp_poi_sub_tmp_var[,1][as.numeric(tmp_poi_sub_tmp_var[,1])>0]), main=paste(all_ecor_names[x]),xlab=variable_names[y]) 
      mtext("EU without Strahler")
      boxplot(tmp_poi_sub_tmp_var[,1][(tmp_poi_sub_tmp_var[,1])>0],outline=FALSE)
      mtext(paste("n=", nrow(tmp_poi_sub)))
      hist(as.numeric(waterbase_sub_tmp_var[,1][as.numeric(waterbase_sub_tmp_var[,1])>0]), main=paste("all waterbase stations"),xlab=variable_names[y]) 
      mtext("EU without Strahler")
      boxplot(waterbase_sub_tmp_var[,1][as.numeric(waterbase_sub_tmp_var[,1])>0],outline=FALSE)
      dev.off()
      
      result<-(wilcox.test(tmp_poi_all_sub_tmp_var$value[tmp_poi_all_sub_tmp_var$value>0]~tmp_poi_all_sub_tmp_var$group[tmp_poi_all_sub_tmp_var$value>0]) )
      
      
      if (is.na(result$p.value)){ 
        print(paste("NA", strahler_ecor_list[x],variable_list[y],ecoregions[l]))
      } else {
        if (result$p.value<0.1){
          #           print((ecoregions[l]))
          #           print(strahler_ecor_list[x])
          #           print(tmp_var<-variable_list[y])
          #           print(result)
          if (median(tmp_poi_sub_tmp_var[,1][(tmp_poi_sub_tmp_var[,1])>0])>median(waterbase_sub_tmp_var[,1][as.numeric(waterbase_sub_tmp_var[,1])>0])){
            print(paste(result$p.value, strahler_ecor_list[x],variable_list[y], "n=", nrow(tmp_poi_sub)))
          } else {
            print(paste("-",result$p.value, strahler_ecor_list[x],variable_list[y], "n=", nrow(tmp_poi_sub)))
          }
        } else {
          print(paste("NS", strahler_ecor_list[x],variable_list[y]))
        }
      }
      
      
    }
  }  
}

### - strahler + ecoregions

all_ecor_list<-c("slipway_waterbase_all_ecor", "marinas_waterbase_all_ecor","fishing_waterbase_all_ecor")

for (x in 1: length(all_ecor_list) ){
  tmp_poi<-all_ecor_list[x]
  tmp_poi_data<-get(tmp_poi)
  tmp_poi_data<-tmp_poi_data[tmp_poi_data$Total_Length<=1000,]
  
  for (y in 1: length(variable_list) ){
    tmp_var<-variable_list[y] 
    
    for (l in 1 : length(ecoregions)){
      
      waterbase_sub<-subset(waterbase_all_ecor, waterbase_all_ecor$NAME==paste(ecoregions[l]))
      tmp_poi_sub<-subset(tmp_poi_data, tmp_poi_data$NAME==paste(ecoregions[l]))
      
      for (z in 1: nrow(tmp_poi_sub)){
        waterbase_sub<-waterbase_sub[waterbase_sub$WaterbaseID!= paste(tmp_poi_sub$WaterbaseID[z]),]
      }
      
      if (nrow(tmp_poi_sub)>10){
        tmp_poi_sub_tmp_var<-tmp_poi_sub[,tmp_var]
        tmp_poi_sub_tmp_var<-data.frame(tmp_poi_sub_tmp_var,rep("tmp_poi"))
        names(tmp_poi_sub_tmp_var)<-c("value","group")
        
        waterbase_sub_tmp_var<-waterbase_sub[,tmp_var]
        waterbase_sub_tmp_var<-data.frame(waterbase_sub_tmp_var,rep("all"))
        names(waterbase_sub_tmp_var)<-c("value","group")
        
        tmp_poi_all_sub_tmp_var<-rbind.data.frame(tmp_poi_sub_tmp_var,waterbase_sub_tmp_var)
        
        
        
        

        
        if (length(unique(tmp_poi_all_sub_tmp_var$group[tmp_poi_all_sub_tmp_var$value>0]))==2){
          result<-(wilcox.test(tmp_poi_all_sub_tmp_var$value[tmp_poi_all_sub_tmp_var$value>0]~tmp_poi_all_sub_tmp_var$group[tmp_poi_all_sub_tmp_var$value>0]) )
        } else {
          print(paste ("all zero",ecoregions[l]))
          {next}
        }  
        
        png(paste("/tmp/","Ecoregions without Strahler",ecoregions[l],all_ecor_names[x],variable_names[y],".png",sep="_"),width = 1760, height = 960, units = "px",pointsize=40)
        par(mfrow = c(1, 4), pty = "s")
        hist(as.numeric(tmp_poi_sub_tmp_var[,1]), main=paste(all_ecor_names[x]),xlab=variable_names[y]) 
        mtext(paste(ecoregions[l],"without Strahler"))
        boxplot(tmp_poi_sub_tmp_var[,1][(tmp_poi_sub_tmp_var[,1])>0],outline=FALSE)
        mtext(paste("n=", nrow(tmp_poi_sub)))
        hist(as.numeric(waterbase_sub_tmp_var[,1]), main=paste("all waterbase stations"),xlab=variable_names[y]) 
        mtext(paste(ecoregions[l],"without Strahler"))
        boxplot(waterbase_sub_tmp_var[,1][as.numeric(waterbase_sub_tmp_var[,1])>0],outline=FALSE)
        dev.off()
        
        if (is.na(result$p.value)){ 
          print(paste("NA", all_ecor_names[x],variable_list[y],ecoregions[l]))
        } else {
          if (result$p.value<0.1){
            #           print((ecoregions[l]))
            #           print(strahler_ecor_list[x])
            #           print(tmp_var<-variable_list[y])
            #           print(result)
            if (median(tmp_poi_sub_tmp_var[,1][(tmp_poi_sub_tmp_var[,1])>0])>median(waterbase_sub_tmp_var[,1][as.numeric(waterbase_sub_tmp_var[,1])>0])&(length(unique(tmp_poi_all_sub_tmp_var$group[tmp_poi_all_sub_tmp_var$value>0]))==2)){
              print(paste(result$p.value,ecoregions[l], strahler_ecor_list[x],variable_list[y], "n=", nrow(tmp_poi_sub)))
            } else {
              print(paste("-",result$p.value, ecoregions[l], all_ecor_names[x],variable_list[y], "n=", nrow(tmp_poi_sub)))
            }
            
          } else {
            print(paste("NS", ecoregions[l],all_ecor_names[x],variable_list[y]))
          }
        }
        
      } else {
        print(paste("less than 10 points", ecoregions[l],all_ecor_names[x],variable_list[y]))
      }
    }  
  }
}

#### Europe and Strahler

for (x in 1: length(strahler_ecor_list) ){
  tmp_poi<-strahler_ecor_list[x]
  tmp_poi_data<-get(tmp_poi)
  tmp_poi_data<-tmp_poi_data[tmp_poi_data$Total_Length<=1000,]
  
  for (y in 1: length(variable_list) ){
    tmp_var<-variable_list[y] 
    
    
    
    waterbase_sub<-subset(waterbase_all_strahler_ecor)
    tmp_poi_sub<-subset(tmp_poi_data)
    
    for (z in 1: nrow(tmp_poi_sub)){
      waterbase_sub<-waterbase_sub[waterbase_sub$WaterbaseID!= paste(tmp_poi_sub$WaterbaseID[z]),]
    }
    
    if (nrow(tmp_poi_sub)>10){
      tmp_poi_sub_tmp_var<-tmp_poi_sub[,tmp_var]
      tmp_poi_sub_tmp_var<-data.frame(tmp_poi_sub_tmp_var,rep("tmp_poi"))
      names(tmp_poi_sub_tmp_var)<-c("value","group")
      
      waterbase_sub_tmp_var<-waterbase_sub[,tmp_var]
      waterbase_sub_tmp_var<-data.frame(waterbase_sub_tmp_var,rep("all"))
      names(waterbase_sub_tmp_var)<-c("value","group")
      
      tmp_poi_all_sub_tmp_var<-rbind.data.frame(tmp_poi_sub_tmp_var,waterbase_sub_tmp_var)
      #tmp_poi_all_sub_tmp_var<-as.data.frame(tmp_poi_all_sub_tmp_var)
      #names(tmp_poi_all_sub_tmp_var)<-c("value","group")
      
      png(paste("/tmp/","EU with Strahler",all_ecor_names[x],variable_names[y],".png",sep="_"),width = 1760, height = 960, units = "px",pointsize=40)
      par(mfrow = c(1, 4), pty = "s")
      hist((tmp_poi_sub_tmp_var[,1][(tmp_poi_sub_tmp_var[,1])>0]), main=paste(all_ecor_names[x]),xlab=variable_names[y])
      mtext(paste("Strahler"))
      boxplot(tmp_poi_sub_tmp_var[,1][(tmp_poi_sub_tmp_var[,1])>0],outline=FALSE)
      mtext(paste("n=", nrow(tmp_poi_sub)))
      hist(as.numeric(waterbase_sub_tmp_var[,1][as.numeric(waterbase_sub_tmp_var[,1])>0]), main=paste("all waterbase stations"),xlab=variable_names[y]) 
      mtext(paste("Strahler"))
      boxplot(waterbase_sub_tmp_var[,1][as.numeric(waterbase_sub_tmp_var[,1])>0],outline=FALSE)
      dev.off()
      
      result<-(wilcox.test(tmp_poi_all_sub_tmp_var$value[tmp_poi_all_sub_tmp_var$value>0]~tmp_poi_all_sub_tmp_var$group[tmp_poi_all_sub_tmp_var$value>0]) )
      #result<-(wilcox.test(as.numeric(as.character(tmp_poi_all_sub_tmp_var$value))[as.numeric(as.character(tmp_poi_all_sub_tmp_var$value))>0]~tmp_poi_all_sub_tmp_var$group[as.numeric(as.character(tmp_poi_all_sub_tmp_var$value))>0]) )
      
      if (is.na(result$p.value)){ 
        print(paste("NA", strahler_ecor_list[x],variable_list[y],ecoregions[l]))
      } else {
        if (result$p.value<0.1){
          #           print((ecoregions[l]))
          #           print(strahler_ecor_list[x])
          #           print(tmp_var<-variable_list[y])
          #           print(result)
          if (median(tmp_poi_sub_tmp_var[,1][(tmp_poi_sub_tmp_var[,1])>0])>median(waterbase_sub_tmp_var[,1][as.numeric(waterbase_sub_tmp_var[,1])>0])){
          print(paste(result$p.value, strahler_ecor_list[x],variable_list[y], "n=", nrow(tmp_poi_sub)))
          } else {
            print(paste("-",result$p.value, strahler_ecor_list[x],variable_list[y], "n=", nrow(tmp_poi_sub)))
          }
        } else {
          print(paste("NS", strahler_ecor_list[x],variable_list[y]))
        }
      }
      
      
    }
  }  
}

### strahler and ecoregions:


for (x in 1: length(strahler_ecor_list) ){
 tmp_poi<-strahler_ecor_list[x]
 tmp_poi_data<-get(tmp_poi)
 tmp_poi_data<-tmp_poi_data[tmp_poi_data$Total_Length<=1000,]
 
for (y in 1: length(variable_list) ){
 tmp_var<-variable_list[y] 

for (l in 1 : length(ecoregions)){
  
  
  waterbase_sub<-subset(waterbase_all_strahler_ecor, waterbase_all_strahler_ecor$NAME==paste(ecoregions[l]))
  tmp_poi_sub<-subset(tmp_poi_data, tmp_poi_data$NAME==paste(ecoregions[l]))
  
   for (z in 1: nrow(tmp_poi_sub)){
       waterbase_sub<-waterbase_sub[waterbase_sub$WaterbaseID!= paste(tmp_poi_sub$WaterbaseID[z]),]
   }
  
  if (nrow(tmp_poi_sub)>10){
    tmp_poi_sub_tmp_var<-tmp_poi_sub[,tmp_var]
    tmp_poi_sub_tmp_var<-data.frame(tmp_poi_sub_tmp_var,rep("tmp_poi"))
    names(tmp_poi_sub_tmp_var)<-c("value","group")
    
    waterbase_sub_tmp_var<-waterbase_sub[,tmp_var]
    waterbase_sub_tmp_var<-data.frame(waterbase_sub_tmp_var,rep("all"))
    names(waterbase_sub_tmp_var)<-c("value","group")
    
    tmp_poi_all_sub_tmp_var<-rbind.data.frame(tmp_poi_sub_tmp_var,waterbase_sub_tmp_var)
    
    

    
      if (length(unique(tmp_poi_all_sub_tmp_var$group[tmp_poi_all_sub_tmp_var$value>0]))==2){
              result<-(wilcox.test(tmp_poi_all_sub_tmp_var$value[tmp_poi_all_sub_tmp_var$value>0]~tmp_poi_all_sub_tmp_var$group[tmp_poi_all_sub_tmp_var$value>0]) )
      } else {
        print(paste ("all zero",ecoregions[l]))
        {next}
      }
        
    png(paste("/tmp/","Ecoregions with Strahler",ecoregions[l],all_ecor_names[x],variable_names[y],".png",sep="_"),width = 1760, height = 960, units = "px",pointsize=40)
    par(mfrow = c(1, 4), pty = "s")
    hist(as.numeric(tmp_poi_sub_tmp_var[,1]), main=paste(all_ecor_names[x]),xlab=variable_names[y]) 
    mtext(paste(ecoregions[l],"Strahler"))
    boxplot(tmp_poi_sub_tmp_var[,1][(tmp_poi_sub_tmp_var[,1])>0],outline=FALSE)
    mtext(paste("n=", nrow(tmp_poi_sub)))
    hist(as.numeric(waterbase_sub_tmp_var[,1]), main=paste("all waterbase stations"),xlab=variable_names[y]) 
    mtext(paste(ecoregions[l],"Strahler"))
    boxplot(waterbase_sub_tmp_var[,1][as.numeric(waterbase_sub_tmp_var[,1])>0],outline=FALSE)
    dev.off()
    
    if (is.na(result$p.value)){ 
      print(paste("NA", strahler_ecor_list[x],variable_list[y],ecoregions[l]))
    } else {
      if (result$p.value<0.1){
        #           print((ecoregions[l]))
        #           print(strahler_ecor_list[x])
        #           print(tmp_var<-variable_list[y])
        #           print(result)
        if (median(tmp_poi_sub_tmp_var[,1][(tmp_poi_sub_tmp_var[,1])>0])>median(waterbase_sub_tmp_var[,1][as.numeric(waterbase_sub_tmp_var[,1])>0])&(length(unique(tmp_poi_all_sub_tmp_var$group[tmp_poi_all_sub_tmp_var$value>0]))==2)){
          print(paste(result$p.value,ecoregions[l], strahler_ecor_list[x],variable_list[y], "n=", nrow(tmp_poi_sub)))
        } else {
          print(paste("-",result$p.value, ecoregions[l], strahler_ecor_list[x],variable_list[y], "n=", nrow(tmp_poi_sub)))
        }
      } else {
        print(paste("NS", ecoregions[l],strahler_ecor_list[x],variable_list[y]))
      }
    }
    }
    
  }
}  
}







    
