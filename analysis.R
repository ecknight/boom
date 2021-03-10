#title: Analysis of common nighthawk territoriality and significance of the wingboom display
#author: Elly C. Knight

options(scipen = 999)
set.seed(1234)

library(adehabitatHR) #For home range mapping
library(sp) #For other mapping stuffs
library(devtools)
library(readxl)
library(tidyverse)
library(lubridate)
library(maptools)#For writing shapefiles
library(geosphere)
library(sf)
library(lme4)
library(MuMIn)
library(ResourceSelection)
library(AICcmodavg)
library(data.table)
library(gridExtra)
library(rgdal)
library(patchwork)
library(rlist)
library(merTools)
library(nlme)

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=16),
        plot.title=element_text(size=12))

load("BoomWorkingSpace.Rdata")

#1. WRANGLING####
nest1 <- read.csv("NestsForTerritoryMapping.csv") %>% 
  rename(BirdID=MaleID, NestX=LocationX, NestY=LocationY) %>% 
  dplyr::filter(BirdID!=1) %>% 
  dplyr::mutate(ID=paste0(BirdID,"-",Year))

boom3 <- read.csv("Booms3.csv") %>% 
  mutate(DateTime = ymd_hms(DateTime)) %>% 
  mutate(ID=paste0(BirdID, "-", Year)) %>% 
  left_join(nest1) %>% 
  mutate(nest=ifelse(is.na(NestID), 0, 1)) %>% 
  dplyr::select(-NestID, -NestX, -NestY)

#Create dataframe to identify sites for each bird
table(boom3$Year, boom3$BirdID)
sites <- data.frame(BirdID = unique(boom3$BirdID),
                    site=c(2,2,3,3,3,1,1,2,4,4,4,4,4,5,5,5,5,5,5,1,4,4,2,2)) %>% 
  mutate(site=case_when(site==3 ~ 1,
                        site==5 ~ 2,
                        site==2 ~ 3,
                        site==1 ~ 4,
                        site==4 ~ 5)) %>% 
  arrange(site, BirdID)

write.csv(sites, "BirdIDSites.csv", row.names = FALSE)

#Add sample size
boom4 <- boom3 %>% 
  group_by(ID) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  left_join(boom3)

#2. SAMPLE SIZE####

#Setup----
table(boom4$Year,boom4$n)
#maxn <- 70
maxn <- max(boom4$n)
boot <- 100
samplesize.kde.list <- list()

#Bootstraps----
set.seed(999)
for(h in 1:boot){
  
  samplesize.kde.dat <- data.frame()
  
  for(i in seq(5, maxn, 5)){
    
    boom.i <- boom4 %>% 
      dplyr::filter(n >= i) %>%
      group_by(ID) %>% 
      sample_n(i, replace=FALSE) %>% 
      ungroup()
    
    boom.sp <- SpatialPointsDataFrame(coords=cbind(boom.i$BoomX, boom.i$BoomY), 
                                      data=data.frame(ID=boom.i$ID),
                                      proj4string = CRS("+proj=utm +zone=12 +datum=WGS84"))
    
    kd <- kernelUD(boom.sp[,1], grid = 1000, h="href", extent=2, same4all=FALSE)
    
    samplesize.kde.dat <- kernel.area(kd, percent=100,
                                      unin="m", unout="ha",
                                      standardize = FALSE) %>% 
      data.frame() %>% 
      transpose(keep.names="ID") %>% 
      rename(hr95=V1) %>% 
      mutate(n=i,
             nbirds=n(),
             boot=h) %>% 
      rbind(samplesize.kde.dat)
    
    print(paste0("FINISHED SAMPLE SIZE ", i))
    
  }
  
  samplesize.kde.list[[h]] <- samplesize.kde.dat
  print(paste0("FINISHED BOOTSTRAP ", h))
  rm(samplesize.kde.dat)
  
}

#Collapse results----
samplesize.kde <- rbindlist(samplesize.kde.list)
write.csv(samplesize.kde, "KDESampleSizeBootstrapResults.csv", row.names = FALSE)
samplesize.kde <- read.csv("KDESampleSizeBootstrapResults.csv")

#NLS
table(samplesize.kde$ID)

IDs <- unique(samplesize.kde$ID)
ls.sum <- data.frame()
ls.fit <- data.frame()
Asym <- data.frame()

for(i in 1:length(IDs)){
  samplesize.kde.i <- samplesize.kde %>% 
    dplyr::filter(ID==IDs[i])
  
  samplesize.nls <- try(nls(hr95 ~ SSlogis(n, Asym, xmid, scal), data = samplesize.kde.i))
  
  if(class(samplesize.nls)=="try-error"){
    next
  }
  
  #Determine asymptote of recall
  Asym.i <-
    round(environment(samplesize.nls[["m"]][["fitted"]])[["env"]][["Asym"]], 3)
  scale <- environment(samplesize.nls[["m"]][["fitted"]])[["env"]][["scal"]]
  xmid <- environment(samplesize.nls[["m"]][["fitted"]])[["env"]][["xmid"]]
  
  Asym <- data.frame(ID=IDs[i],
                     Asym=Asym.i) %>% 
    rbind(Asym)
  
  #Fit growth curve to new data
  new_frac <- seq(min(samplesize.kde.i$n), max(samplesize.kde.i$n),
                  length.out = 1000)
  
  ls.fit.i <- data.frame(
    r=predict(newdata = data.frame(n = new_frac), object = samplesize.nls),
    n = new_frac) %>% 
    mutate(ID = IDs[i])
  
  ls.fit <- rbind(ls.fit, ls.fit.i)
  
  #Find sample size for 99% of asymptote
  ls.sum.i <- ls.fit.i %>%
    mutate(r = round(r, digits = 3)) %>%
    filter(r >= 0.90 * Asym.i) %>% 
    mutate(ID = IDs[i])
  ls.sum <- rbind(ls.sum, ls.sum.i)
  
  print(paste0("Finished bird ", i, " of ", length(IDs)))
  
}

write.csv(ls.fit, "KDESampleSizeNLSPredictions.csv", row.names = FALSE)
write.csv(Asym, "KDESampleSizeNLSAsymptotes.csv", row.names = FALSE)

#Visualize----
ggplot(ls.fit) +
  geom_point(data=samplesize.kde, aes(x=n, y=hr95)) +
  geom_line(aes(x=n, y=r), colour="blue") +
  geom_hline(data=Asym, aes(yintercept=Asym), colour="red") +
  facet_wrap(~ID, scales="free")

#Summarize into selected IDs----
ls.sum.summary <- ls.sum %>% 
  group_by(ID) %>% 
  summarize(n=min(n)) %>% 
  ungroup() %>% 
  mutate(ID=str_sub(ID, 2, 10),
         ID=gsub(ID, pattern=".201", replacement="-201"))
ls.sum.summary

#Subset data according to results----
boom5 <- boom4 %>% 
  dplyr::mutate(select=ifelse(ID %in% ls.sum.summary$ID, 1, 0))

boom5 %>% 
  dplyr::select(ID, n, select) %>% 
  unique() %>% 
  arrange(n) %>% 
  data.frame()

boom6 <- boom5 %>% 
  dplyr::filter(select==1)

#3. CALCULATE KDE FOR ALL####

#Calculate KDE----
boom.all <- SpatialPointsDataFrame(coords=cbind(boom6$BoomX, boom6$BoomY), 
                                   data=data.frame(ID=boom6$ID),
                                   proj4string = CRS("+proj=utm +zone=12 +datum=WGS84"))

kd.all <- kernelUD(boom.all, grid = 1000, extent=2, h="href", same4all=FALSE)

#Calculate area----
kd.area.all <- kernel.area(kd.all, percent=seq(5, 95, by=5),
                           unin="m", unout="ha",
                           standardize = FALSE) %>% 
  data.frame() %>% 
  transpose(keep.names="ID") %>% 
  rename(hr05=V1,
         hr10=V2,
         hr15=V3,
         hr20=V4,
         hr25=V5,
         hr30=V6,
         hr35=V7,
         hr40=V8,
         hr45=V9,
         hr50=V10,
         hr55=V11,
         hr60=V12,
         hr65=V13,
         hr70=V14,
         hr75=V15,
         hr80=V16,
         hr85=V17,
         hr90=V18,
         hr95=V19) %>% 
  mutate(ID=str_sub(ID, 2, 20)) %>% 
  separate(ID, into=c("BirdID", "Year")) %>% 
  mutate(ID=paste0(BirdID, "-", Year)) %>% 
  dplyr::select(-BirdID, -Year) %>% 
  left_join(boom6)

write.csv(kd.area.all, "KDEAreaSummaryAllBirds.csv", row.names=FALSE)

#Make shapefiles
#Vector home range calculation
hr.95 <- getverticeshr(kd.all, percent=95)
hr.90 <- getverticeshr(kd.all, percent=90)
hr.85 <- getverticeshr(kd.all, percent=85)
hr.80 <- getverticeshr(kd.all, percent=80)
hr.75 <- getverticeshr(kd.all, percent=75)
hr.70 <- getverticeshr(kd.all, percent=70)
hr.65 <- getverticeshr(kd.all, percent=65)
hr.60 <- getverticeshr(kd.all, percent=60)
hr.55 <- getverticeshr(kd.all, percent=55)
hr.50 <- getverticeshr(kd.all, percent=50)
hr.45 <- getverticeshr(kd.all, percent=45)
hr.40 <- getverticeshr(kd.all, percent=40)
hr.35 <- getverticeshr(kd.all, percent=35)
hr.30 <- getverticeshr(kd.all, percent=30)
hr.25 <- getverticeshr(kd.all, percent=25)
hr.20 <- getverticeshr(kd.all, percent=20)
hr.15 <- getverticeshr(kd.all, percent=15)
hr.10 <- getverticeshr(kd.all, percent=10)
hr.05 <- getverticeshr(kd.all, percent=5)

#Export kernel density shapefiles
writeOGR(obj=hr.95, dsn="shapefiles", layer="HR95", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.90, dsn="shapefiles", layer="HR90", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.85, dsn="shapefiles", layer="HR85", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.80, dsn="shapefiles", layer="HR80", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.75, dsn="shapefiles", layer="HR75", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.70, dsn="shapefiles", layer="HR70", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.65, dsn="shapefiles", layer="HR65", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.60, dsn="shapefiles", layer="HR60", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.55, dsn="shapefiles", layer="HR55", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.50, dsn="shapefiles", layer="HR50", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.45, dsn="shapefiles", layer="HR45", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.40, dsn="shapefiles", layer="HR40", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.35, dsn="shapefiles", layer="HR35", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.30, dsn="shapefiles", layer="HR30", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.25, dsn="shapefiles", layer="HR25", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.20, dsn="shapefiles", layer="HR20", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.15, dsn="shapefiles", layer="HR15", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.10, dsn="shapefiles", layer="HR10", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=hr.05, dsn="shapefiles", layer="HR05", driver="ESRI Shapefile", overwrite_layer = TRUE)

#Summarize area----
kd.area.all <- read.csv("KDEAreaSummaryAllBirds.csv")
summary(kd.area.all)
sd(kd.area.all$hr95)

#Without the one huge one
kd.area.all.sub <- kd.area.all %>% 
  dplyr::filter(hr95 < 40)

ggplot(kd.area.all.sub) +
  geom_violin(aes(x=factor(Year), y=hr95))

t.test(hr95 ~ Year, data=kd.area.all.sub)

#4. ADJACENT OVERLAP####
#Identify adjacent pairs---
adjacent <- data.frame(BirdPr = c("4-5",
                                  "4-6",
                                  "5-6",
                                  "32-36",
                                  "29-32",
                                  "15-17",
                                  "15-44",
                                  "19-22",
                                  "20-21",
                                  "21-76",
                                  "76-52"))

#2016----
birds.2016 <- boom6 %>% 
  dplyr::filter(Year==2016) %>% 
  dplyr::select(BirdID) %>% 
  unique()

boom.2016 <- boom6 %>% 
  dplyr::filter(Year==2016)
boom.2016.sp <- SpatialPointsDataFrame(coords=cbind(boom.2016$BoomX, boom.2016$BoomY), 
                                       data=data.frame(ID=boom.2016$ID),
                                       proj4string = CRS("+proj=utm +zone=12 +datum=WGS84"))

overlap1 <- kerneloverlap(boom.2016.sp, method="PHR", grid=2000, percent=95)

overlap.2016 <- overlap1 %>% 
  as.data.frame() %>% 
  mutate(ID1=row.names(overlap1)) %>% 
  pivot_longer(cols='1-2016':'6-2016',
               names_to="ID2",
               values_to="Overlap") %>% 
  separate(ID1, into=c("BirdID1", "Year1"), remove=FALSE) %>% 
  separate(ID2, into=c("BirdID2", "Year2"), remove=FALSE) %>% 
  mutate(BirdID1 = as.numeric(BirdID1),
         BirdID2 = as.numeric(BirdID2)) %>% 
  left_join(sites %>% 
              dplyr::rename(BirdID1=BirdID, site1=site)) %>% 
  left_join(sites %>% 
              dplyr::rename(BirdID2=BirdID, site2=site)) %>% 
  dplyr::filter(site1==site2,
                BirdID1!=BirdID2) %>% 
  unique()

#2017----
birds.2017 <- boom6 %>% 
  dplyr::filter(Year==2017) %>% 
  dplyr::select(BirdID) %>% 
  unique()

boom.2017 <- boom6 %>% 
  dplyr::filter(Year==2017)
boom.2017.sp <- SpatialPointsDataFrame(coords=cbind(boom.2017$BoomX, boom.2017$BoomY), 
                                       data=data.frame(ID=boom.2017$ID),
                                       proj4string = CRS("+proj=utm +zone=12 +datum=WGS84"))

overlap2 <- kerneloverlap(boom.2017.sp, method="PHR", grid=2000, percent=95)

overlap.2017 <- overlap2 %>% 
  as.data.frame() %>% 
  mutate(ID1=row.names(overlap2)) %>% 
  pivot_longer(cols='15-2017':'93-2017',
               names_to="ID2",
               values_to="Overlap") %>% 
  separate(ID1, into=c("BirdID1", "Year1"), remove=FALSE) %>% 
  separate(ID2, into=c("BirdID2", "Year2"), remove=FALSE) %>% 
  mutate(BirdID1 = as.numeric(BirdID1),
         BirdID2 = as.numeric(BirdID2)) %>% 
  left_join(sites %>% 
              dplyr::rename(BirdID1=BirdID, site1=site)) %>% 
  left_join(sites %>% 
              dplyr::rename(BirdID2=BirdID, site2=site)) %>% 
  dplyr::filter(site1==site2,
                BirdID1!=BirdID2) %>% 
  unique()

#Put together----
overlap.neighbours <- rbind(overlap.2016, overlap.2017) %>% 
  rowwise() %>% 
  mutate(minbird = min(BirdID1, BirdID2),
         maxbird = max(BirdID1, BirdID2),
         BirdPr = paste0(minbird, "-", maxbird)) %>% 
  group_by(BirdPr, Year1) %>% 
  summarize(Overlap.mn=round(mean(Overlap), 3),
            Overlap.sd=round(sd(Overlap), 3)) %>% 
  ungroup() %>% 
  inner_join(adjacent %>% 
               mutate(BirdPr=as.character(BirdPr))) %>% 
  arrange(-Overlap.mn)
write.csv(overlap.neighbours, "PHRForNeighbours.csv", row.names = FALSE) 

#T-test----
t.test(overlap.neighbours$Overlap.mn, alternative="greater", mu=0)

#5. INTERANNUAL OVERLAP####
birds.year <- boom6 %>% 
  dplyr::select(BirdID, Year) %>% 
  unique() %>% 
  group_by(BirdID) %>% 
  summarize(years=n()) %>% 
  ungroup() %>% 
  dplyr::filter(years==2)

#Calculate overlap----
boom.year <- boom6 %>% 
  dplyr::filter(BirdID %in% birds.year$BirdID)
boom.year.sp <- SpatialPointsDataFrame(coords=cbind(boom.year$BoomX, boom.year$BoomY), 
                                       data=data.frame(ID=boom.year$ID),
                                       proj4string = CRS("+proj=utm +zone=12 +datum=WGS84"))

overlap3 <- kerneloverlap(boom.year.sp, method="PHR", grid=2000, percent=95)

overlap.years <- overlap3 %>% 
  as.data.frame() %>% 
  mutate(ID1=row.names(overlap3)) %>% 
  pivot_longer(cols='15-2016':'6-2017',
               names_to="ID2",
               values_to="Overlap") %>% 
  separate(ID1, into=c("BirdID1", "Year1"), remove=FALSE) %>% 
  separate(ID2, into=c("BirdID2", "Year2"), remove=FALSE) %>% 
  mutate(BirdID1 = as.numeric(BirdID1),
         BirdID2 = as.numeric(BirdID2)) %>% 
  left_join(sites %>% 
              dplyr::rename(BirdID1=BirdID)) %>% 
  dplyr::filter(BirdID1==BirdID2,
                Year1 != Year2) %>% 
  unique() %>% 
  group_by(BirdID1, site) %>% 
  summarize(Overlap.mn=mean(Overlap),
            Overlap.sd=sd(Overlap)) %>% 
  ungroup() %>% 
  rename(BirdID = BirdID1)
overlap.years

write.csv(overlap.years, "PHRBetweenYears.csv", row.names = FALSE)  

#T-test of overlap----
overlap.years <- read.csv("PHRBetweenYears.csv")

t.test(overlap.years$Overlap.mn, mu=0)

summary(overlap.years)
sd(overlap.years$Overlap.mn)

#T-test between years----
kd.area.pr <- kd.area.all %>% 
  dplyr::filter(BirdID %in% overlap.years$BirdID) %>% 
  dplyr::select(BirdID, hr95, Year) %>% 
  unique()

t.test(kd.area.pr$hr95 ~ kd.area.pr$Year, alternative="two.sided", paired=TRUE)

#Distance between nests----
birds.year.nest <- expand.grid(BirdID = birds.year$BirdID, Year=c("2016", "2017")) %>% 
  left_join(nest1) %>% 
  dplyr::filter(!is.na(NestID)) %>% 
  group_by(BirdID) %>% 
  mutate(n=n()) %>% 
  ungroup() %>% 
  dplyr::filter(n > 1)

nest.2016.a <- birds.year.nest %>% 
  dplyr::filter(Year==2016,
                NestID !=27) %>% 
  st_as_sf(coords=c("NestX", "NestY"), crs="+proj=utm +zone=12 +datum=WGS84") %>% 
  dplyr::select(BirdID)

nest.2016.b <- birds.year.nest %>% 
  dplyr::filter(Year==2016,
                NestID !=12) %>% 
  st_as_sf(coords=c("NestX", "NestY"), crs="+proj=utm +zone=12 +datum=WGS84") %>% 
  dplyr::select(BirdID)

nest.2017 <- birds.year.nest %>% 
  dplyr::filter(Year==2017,
                NestID !=12) %>% 
  st_as_sf(coords=c("NestX", "NestY"), crs="+proj=utm +zone=12 +datum=WGS84") %>% 
  dplyr::select(BirdID)

nest.distances <- data.frame(distance = as.numeric(st_distance(nest.2016.a,
                                                               nest.2017,
                                                               by_element=TRUE)),
                             BirdID = nest.2017$BirdID) %>% 
  rbind(data.frame(distance = as.numeric(st_distance(nest.2016.b,
                                                     nest.2017,
                                                     by_element=TRUE)),
                   BirdID = nest.2017$BirdID)) %>% 
  unique()

summary(nest.distances)
sd(nest.distances$distance)

#5. NEST RSF####
#Wrangle----
#add nest data
nest1 <- read.csv("NestsForTerritoryMapping.csv") %>% 
  dplyr::filter(BirdID!=1) %>% 
  dplyr::mutate(ID=paste0(BirdID,"-",Year))

boom6 <- boom6 %>% 
  dplyr::left_join(nest1, by=c("BirdID", "Year")) %>% 
  dplyr::mutate(ID=paste0(BirdID,"-",Year))


#Clean up year with two nests for same male
boom7 <- boom6 %>% 
  dplyr::filter(!(NestID==12 & DateTime > "2016-07-02 00:00:00"),
                !(NestID==27 & DateTime < "2016-07-10 00:00:00"))



#Calculate KDE----
boom.sp <- SpatialPointsDataFrame(coords=cbind(boom6$BoomX, boom6$BoomY), 
                                  data=data.frame(ID=boom7$ID),
                                  proj4string = CRS("+proj=utm +zone=12 +datum=WGS84"))

kd <- kernelUD(boom.sp, grid = 1000, h="href", same4all=FALSE)

#Calculate area----
kd.area <- kernel.area(kd, percent=seq(5, 95, by=5),
                       unin="m", unout="ha",
                       standardize = FALSE) %>% 
  data.frame() %>% 
  transpose(keep.names="ID") %>% 
  rename(hr05=V1,
         hr10=V2,
         hr15=V3,
         hr20=V4,
         hr25=V5,
         hr30=V6,
         hr35=V7,
         hr40=V8,
         hr45=V9,
         hr50=V10,
         hr55=V11,
         hr60=V12,
         hr65=V13,
         hr70=V14,
         hr75=V15,
         hr80=V16,
         hr85=V17,
         hr90=V18,
         hr95=V19) %>% 
  mutate(ID=str_sub(ID, 2, 20)) %>% 
  separate(ID, into=c("BirdID", "Year")) %>% 
  mutate(ID=paste0(BirdID, "-", Year)) %>% 
  dplyr::select(-BirdID, -Year)


#4. CREATE AVAILABLE POINTS####
#Create polygons
kd.shp <- getverticeshr(kd, 95)

#Create 100 random available points in each 95% KDE
#https://ecosystems.psu.edu/research/labs/walter-lab/manual/chapter-8-resource-selection/8.2-preparing-additional-covariates
random <- sapply(slot(kd.shp, "polygons"), 
                 function(i) spsample(i, n=100, type="random", offset=c(0,0),
                                      proj4string = CRS("+proj=utm +zone=12 +datum=WGS84")))
random.merged <- do.call("rbind", random) 
ids <- sapply(slot(kd.shp, "polygons"), function(i) slot(i, "ID"))
newpts <- sapply(random, function(i) nrow(i@coords))
newpts
pt_id <- rep(ids, newpts)
random.final <- random.merged %>% 
  st_as_sf() %>% 
  mutate(ID=pt_id,
         type=0) 

#Plot to check
ggplot() +
  geom_sf(aes(colour=ID), data=random.final) + 
  geom_point(aes(x=BoomX, y=BoomY, colour=ID), data=boom7, shape=4)

#Put used and available together
rand.sf <- st_as_sf(random.final) %>% 
  mutate(Type=0)

all.sf <- st_as_sf(boom.sp) %>% 
  mutate(type=1) %>% 
  rbind(random.final) %>% 
  left_join(kd.area, by=c("ID")) %>% 
  left_join(nest1, by=c("ID")) 

#5. CALCULATE DISTANCE TO NEST####
boom.coord <- SpatialPoints(coords=st_coordinates(all.sf), proj4string = CRS("+proj=utm +zone=12 +datum=WGS84")) %>% 
  sp::spTransform("+proj=longlat +datum=WGS84")
nest.coord <- SpatialPoints(coords=cbind(all.sf$NestX, all.sf$NestY), proj4string = CRS("+proj=utm +zone=12 +datum=WGS84")) %>% 
  sp::spTransform(CRS("+proj=longlat +datum=WGS84"))
Distance <- distHaversine(boom.coord, nest.coord)

final <- all.sf %>% 
  st_coordinates() %>% 
  cbind(all.sf) %>% 
  cbind(Distance) %>% 
  mutate(Distance.st = scale(Distance)[,1],
         area=hr95,
         area.st=scale(area)[,1]) %>% 
  data.frame() %>% 
  dplyr::select(-geometry)

table(final$BirdID, final$Year, final$type)

write.csv(final, "UsedAvailableData.csv", row.names=FALSE)
final <- read.csv("UsedAvailableData.csv")

#Visualize
ggplot(final) +
  geom_violin(aes(x=factor(type), y=Distance.st, colour=factor(type)))

#6. RSF####
#https://ecosystems.psu.edu/research/labs/walter-lab/manual/chapter-8-resource-selection/8.4-resource-selection-functions
#Model
fit1 <- glmer(type ~ Distance.st + (1|Year/BirdID), data=final,
              family=binomial(link="logit")) #Distance model
fit2 <- glmer(type ~ 1 + (1|Year/BirdID), data=final,
              family=binomial(link="logit")) #Null model
fit3 <- glmer(type ~ Distance.st*area + (1|Year/BirdID), data=final,
              family=binomial(link="logit")) #Model with MCP
aictab(list(fit1, fit2, fit3), sort=F)

#Taking out the year RE to try and deal with singular fit
fit1 <- glmer(type ~ Distance + (1|BirdID), data=final,
              family=binomial(link="logit")) #Distance model
fit2 <- glmer(type ~ 1 + (1|BirdID), data=final,
              family=binomial(link="logit")) #Null model
fit3 <- glmer(type ~ Distance*area + (1|BirdID), data=final,
              family=binomial(link="logit")) #Model with MCP
aictab(list(fit1, fit2, fit3), sort=F)

#Use scaled distance to deal with large eigenvalue
fit1 <- glmer(type ~ Distance.st + (1|BirdID), data=final,
              family=binomial(link="logit")) #Distance model
fit2 <- glmer(type ~ 1 + (1|BirdID), data=final,
              family=binomial(link="logit")) #Null model
fit3 <- glmer(type ~ Distance.st*area.st + (1|BirdID), data=final,
              family=binomial(link="logit")) #Model with MCP
aictab(list(fit1, fit2, fit3), sort=F)

#7. PLOT RSF####
#Backtransform
att.dist <- attributes(scale(final$Distance))
att.area <- attributes(scale(final$area))

newdat <- expand.grid(Distance.st=seq(min(final$Distance.st), max(final$Distance.st), 0.01),
                      area.st=c(-1.0065, 0.072, 1.342, 1.976),
                      BirdID=unique(final$BirdID)) %>% 
  mutate(Distance=Distance.st*att.dist$`scaled:scale` + att.dist$`scaled:center`,
         Area=area.st*att.area$`scaled:scale` + att.area$`scaled:center`,
         Area.rd=round(Area, 1))
table(newdat$Area)
  
pred <- predictInterval(fit3, newdata = newdat, which="fixed", n.sims = 1000, type="probability", level=0.95) %>% 
  cbind(newdat)

#get mean across the bird IDs
pred.mn <- pred %>% 
  group_by(Distance, Area.rd) %>% 
  summarize(fit=mean(fit),
            upr=mean(upr),
            lwr=mean(lwr)) %>% 
  ungroup()

#visualize
ggplot(pred.mn) +
  geom_line(aes(x=Distance, y=fit, colour=factor(Area.rd))) +
  geom_ribbon(aes(x=Distance, ymin=lwr, ymax=upr, group=factor(Area.rd)), alpha=0.3)

write.csv(pred.mn, "RSFPredictions.csv", row.names=FALSE)

save.image("BoomWorkingSpace.Rdata")
