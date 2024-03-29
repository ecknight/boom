#title: Figures for analysis of common nighthawk territoriality and significance of the wingboom display
#author: Elly C. Knight

library(tidyverse)
library(lubridate)
library(ggmap)
library(sf)
library(gridExtra)
library(rosm)
library(prettymapr)
library(ggspatial)
library(patchwork)
library(RColorBrewer)
library(ggsn)
library(seewave)
library(tuneR)
library(grid)
library(nord)
library(cowplot)
library(colorspace)
library(pBrackets) #for brackets on spectrogram
library(Cairo)

my.theme <- theme_classic() +
  theme(text=element_text(size=12, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(size=12, hjust = 0.5))

#WRANGLING####

boom3 <- read.csv("Booms3.csv") %>% 
  mutate(DateTime = ymd_hms(DateTime)) %>% 
  mutate(ID=paste0(BirdID, "-", Year))

sites <- read.csv("BirdIDSites.csv")

nest <- read.csv("NestsForTerritoryMapping.csv") %>% 
  dplyr::filter(BirdID!=1) %>% 
  dplyr::mutate(ID=paste0(BirdID,"-",Year)) %>% 
  dplyr::select(BirdID, Year, NestID)

#Put together
boom4 <- boom3 %>% 
  group_by(ID) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  left_join(boom3) %>% 
  left_join(sites) %>% 
  dplyr::filter(n >= 5) %>% 
  left_join(nest) %>% 
  mutate(nest = ifelse(is.na(NestID), 0, 1)) %>% 
  dplyr::select(-NestID)

#FIGURE 1. SPECTROGRAM####

bracketsGrob <- function(...){
  l <- list(...)
  e <- new.env()
  e$l <- l
  grid:::recordGrob(  {
    do.call(grid.brackets, l)
  }, e)
}

wav <- readWave("/Users/ellyknight/Documents/UoA/Presentations/Materials/Recordings/CONI Common nighthawk - boom.wav")
wav.f <- bwfilter(wav, from=200, to=5400, output="Wave")

Cairo(width = 6, height = 8, file="Figures/Fig1Spectrogram.jpeg", type="jpeg", pointsize=12,
      bg = "transparent", canvas = "white", units = "in", dpi = 300)

spectro <- ggspectro(wav.f, flim=c(0,6), wl=1024, ovlp=60) +
  stat_contour(geom="polygon", aes(fill=..level..), bins=30, show.legend=FALSE) +
  scale_fill_continuous(name="Amplitude\n(dB)\n", limits=c(-40,0),
                        na.value="transparent", low="white", high="black") + 
  xlim(c(0, 0.9)) +
  xlab("Time (s)") +
  ylab("Frequency (kHz)") +
  geom_text(x=0, y=1, label="Wingboom", angle=90, size=5, family="Arial") + 
  geom_text(x=0, y=4, label="Call", angle=90, size=5, family="Arial") +
  my.theme +
  theme(axis.title = element_text(size=14))
spectro
grid.brackets(x1=400, x2=400, y1=1480, y2=100, lwd=2, col="black")
grid.brackets(x1=400, x2=400, y1=2100, y2=1480, lwd=2, col="black")

dev.off()

 
#FIGURE 2. STUDY AREA####
#Part 1. Study area----
nam <- map_data("world", region=c("Canada", 
                                  "USA", 
                                  "Mexico",
                                  "Guatemala", 
                                  "Belize", 
                                  "El Salvador",
                                  "Honduras", 
                                  "Nicaragua", 
                                  "Costa Rica",
                                  "Panama", 
                                  "Jamaica", 
                                  "Cuba", 
                                  "The Bahamas",
                                  "Haiti", 
                                  "Dominican Republic", 
                                  "Antigua and Barbuda",
                                  "Dominica", 
                                  "Saint Lucia", 
                                  "Saint Vincent and the Grenadines", 
                                  "Barbados",
                                  "Grenada",
                                  "Trinidad and Tobago")) %>% 
  dplyr::filter(!group%in%c(258:264))

nam.eq <- nam %>% 
  st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
  st_transform(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") %>%
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(nam)

area.eq.center <- boom4 %>% 
  st_as_sf(coords=c("BoomX", "BoomY"), crs="+proj=utm +zone=12 +datum=WGS84") %>% 
  st_transform(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") %>%
  st_coordinates() %>% 
  as.data.frame() %>% 
  dplyr::summarize(X=mean(X),
                   Y=mean(Y))
map.nam <- ggplot() +
  geom_polygon(data=nam.eq, aes(x=X, y=Y, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_text(label="★", aes(x=X, y=Y), size=10, family = "HiraKakuPro-W3", data=area.eq.center, colour="black") +
  xlim(c(-4000000, 3000000)) +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  my.theme +
  theme(plot.margin = unit(c(0,0,-1,-1), "cm"))
map.nam

#Part 2. Study sites----
#Wrangle study sites
pts.sf <- boom4 %>% 
  st_as_sf(coords=c("BoomX", "BoomY"), crs="+proj=utm +zone=12 +datum=WGS84") %>% 
  st_transform(crs=4326)

pts.sf.2016 <- pts.sf %>% 
  dplyr::filter(Year==2016)

pts.sf.2017 <- pts.sf %>% 
  dplyr::filter(Year==2017)

pts.wgs <- pts.sf %>% 
  st_coordinates() %>% 
  cbind(boom4) %>% 
  as.data.frame()

pts.wgs.2016 <- pts.wgs %>% 
  dplyr::filter(Year==2016)

pts.wgs.2017 <- pts.wgs %>% 
  dplyr::filter(Year==2017)

area.wgs.center <- pts.wgs %>% 
  dplyr::summarize(X=mean(X),
                   Y=mean(Y))

sites.wgs.center <- pts.wgs %>% 
  group_by(site) %>% 
  dplyr::summarize(X=mean(X),
                   Y=mean(Y)) %>% 
  ungroup()

sites.wgs.sf <- sites.wgs.center %>% 
  st_as_sf(coords=c("X", "Y"), crs=4326)

#Get background data
register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

map <- get_map(area.wgs.center, zoom=9, force=TRUE, maptype="satellite", color="color")

map_attributes <- attributes(map)

map_transparent <- matrix(adjustcolor(map, 
                                      alpha.f = 0.8), 
                          nrow = nrow(map))
attributes(map_transparent) <- map_attributes

#Set colours
clrs.sa <- c("#D35000","#9C3E5D","#8175CB","#5087C1","#00767A")

#map
map.site <- ggmap(map_transparent) +
  geom_spatial_point(aes(x = X, y = Y,
                         colour=factor(site)),
                     data = sites.wgs.center, 
                     crs=4326,
                     alpha = 0.7,
                     size=10,
                     show.legend = FALSE) +
  geom_text(aes(x=X, y=Y, label=site),
            data = sites.wgs.center,
            size=6,
            colour="grey20") +
  geom_text(label="McLelland\nLake", aes(x=-111.315, y=57.485), size=4, colour="grey80") +
  geom_text(label="Athabasca River", aes(x=-111.508, y=57.57), size=4, colour="grey80", angle = 85) +
  scale_colour_manual(values=clrs.sa) +
  ggspatial::annotation_north_arrow(location = "tr",
                                    style = ggspatial::north_arrow_orienteering(fill = c("grey80", "grey20"), line_col = "grey20")) +
  ggsn::scalebar(x.min = -111.75, x.max = -111.55, 
                 y.min = 57.38, y.max = 57.7, 
                 transform=TRUE, model="WGS84",
                 dist=5, dist_unit="km",
                 box.fill=c("grey80", "grey20"),
                 box.color="grey20") +
  coord_sf(crs=4326) +
  my.theme +
  xlab("") +
  ylab("") +
  xlim(c(min(area.wgs.center$X)-0.35, max(area.wgs.center$X)+0.25)) +
  ylim(c(min(area.wgs.center$Y)-0.15, max(area.wgs.center$Y)+0.15)) +
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"),
        legend.position = "right",
        axis.text.x.bottom=element_text(size=10),
        axis.text.y.left=element_text(size=10)) +
  guides(colour = guide_legend(override.aes = list(size=5),
                               ncol=2))
#map.site

#Put it together####
plot.sa <- map.site +
  inset_element(map.nam,
                right=0.35,
                bottom=0.6,
                left=0.02,
                top=0.98)
#plot.sa

ggsave("Figures/Fig2StudyArea.jpeg", device="jpeg", width=8, height=7, units="in", plot=plot.sa, dpi=600)


#FIGURE 3. KDE RESULTS####

#Read stuff in----
area <- read.csv("KDEAreaSummaryAllBirds.csv") %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE) %>% 
  dplyr::select(BirdID, Year, hr95, hr50) %>% 
  unique()

hr.95 <- read_sf("shapefiles/HR95.shp") %>% 
  dplyr::rename(ID=id) %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE) %>% 
  mutate(BirdID = as.numeric(BirdID)) %>% 
  left_join(sites) %>% 
  st_transform(crs=4326)

hr.50 <- read_sf("shapefiles/HR50.shp") %>% 
  dplyr::rename(ID=id) %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE) %>% 
  mutate(BirdID = as.numeric(BirdID)) %>% 
  left_join(sites) %>% 
  st_transform(crs=4326)

nest1 <- read.csv("NestsForTerritoryMapping.csv") %>% 
  dplyr::filter(BirdID!=1) %>% 
  mutate(Year = as.character(Year),
         nest="Nest") %>% 
  dplyr::filter(ID %in% c(hr.95$ID))

#Set up colours----
birds <- hr.95 %>% 
  as.data.frame() %>% 
  dplyr::select(site, BirdID) %>% 
  unique() %>% 
  arrange(site, BirdID)
birds

clrs1 <- data.frame(BirdID=c(NA,4,5,6,NA),
                   col=sequential_hcl(5, "Oranges")) %>% 
  dplyr::filter(!is.na(BirdID))
pie(rep(1,3), col=clrs1$col)

clrs2 <- data.frame(BirdID=c(NA,32,NA,36,NA),
                    col=sequential_hcl(5, "Burg")) %>% 
  dplyr::filter(!is.na(BirdID))
pie(rep(1,2), col=clrs2$col)

clrs3 <- data.frame(BirdID=c(1,2,93,NA),
                    col=sequential_hcl(4, "Purples 3")) %>% 
  dplyr::filter(!is.na(BirdID))
pie(rep(1,3), col=clrs3$col)

clrs4 <- data.frame(BirdID=c(15,17,44,NA),
                    col=sequential_hcl(4, "Blues")) %>% 
  dplyr::filter(!is.na(BirdID))
pie(rep(1,3), col=clrs4$col)

clrs5 <- data.frame(BirdID=c(19,20,21,22,54,52),
                    col=sequential_hcl(6, "ag_Grnyl")) %>% 
  dplyr::filter(!is.na(BirdID))
pie(rep(1,6), col=clrs5$col)

clrs <- rbind(clrs1, clrs2, clrs3, clrs4, clrs5)
pie(rep(1,17), col=clrs$col)

#Wrangling----
all <- rbind(hr.95 %>% 
               mutate(iso = 95),
             hr.50 %>% 
               mutate(iso=50)) %>% 
  left_join(nest1) %>% 
  dplyr::filter(!is.na(area)) %>% 
  left_join(clrs)

#Plot KDE----
sizes <- data.frame(site=c(1:5),
                    size=c(200,300,400,300,400))

plots <- data.frame(all) %>% 
  dplyr::select(site, Year) %>% 
  rename(year=Year) %>% 
  unique() %>% 
  arrange(site, year)

plots.list <- list()
for(i in 1:nrow(plots)){
  
  site.i <- plots$site[i]
  year.i <- plots$year[i]
  
  iso.i <- all %>% 
    dplyr::filter(site==site.i,
                  Year==year.i) %>% 
    dplyr::select(-NestID, -NestX, -NestY, -nest, -X, -Y) %>% 
    unique()
  
  iso.site.i <- all %>% 
    dplyr::filter(site==site.i) %>% 
    dplyr::select(-NestID, -NestX, -NestY, -nest, -X, -Y) %>% 
    unique()
  
  size.i <- sizes %>% 
    dplyr::filter(site==site.i)
  
  nest.i <- all %>% 
    dplyr::filter(site==site.i,
                  Year==year.i) 
  
  clrs.i <- nest.i %>% 
    as.data.frame() %>% 
    dplyr::select(BirdID, Year, ID, col) %>% 
    unique() %>% 
    arrange(Year, BirdID)
  
  new_bb = c(st_bbox(iso.site.i)[['xmin']], st_bbox(iso.site.i)[['ymin']]-size.i$size*0.000005,
             st_bbox(iso.site.i)[['xmax']], st_bbox(iso.site.i)[['ymax']])
  names(new_bb) = c("xmin", "ymin", "xmax", "ymax")
  attr(new_bb, "class") = "bbox"
  
  attr(st_geometry(iso.i), "bbox") = new_bb
  
  xmin <- st_bbox(iso.site.i)[['xmin']]-size.i$size*0.00001
  xmax <- st_bbox(iso.site.i)[['xmax']]+size.i$size*0.00001
  ymin <- st_bbox(iso.site.i)[['ymin']]-size.i$size*0.000005
  ymax <- st_bbox(iso.site.i)[['ymax']]
  
  plot.i <- ggplot() +
    geom_sf(data=iso.i, aes(colour=factor(BirdID), alpha=factor(iso)), lwd=1, inherit.aes = FALSE, fill="black", show.legend = FALSE) +
    scale_colour_manual(values=clrs.i$col) +
    scale_alpha_manual(values=c(0.4, 0.1), name="Isopleth") +
    xlab("") +
    ylab("") +
    xlim(c(xmin, xmax)) +
    ylim(c(ymin, ymax)) +
    ggtitle(paste0("Site ", site.i, "\n", "Year ", year.i)) +
    my.theme +
    theme(plot.margin = unit(c(0,-1,0,-1), "cm"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size=12)) +
    coord_sf(datum=NA) +
    ggsn::scalebar(data=iso.i,
                   transform=TRUE, model="WGS84",
                   dist=size.i$size, dist_unit="m",
                   box.fill=c("grey80", "grey20"),
                   box.color="grey20",
                   st.bottom=FALSE,
                   st.size=3,
                   st.dist=0.04,
                   border.size = 0.5,
                   height = 0.05,
                   location="bottomleft")
  plot.i
  
  plots.list[[i]] <- plot.i + 
    geom_point(data=nest.i, aes(x=X, y=Y, fill=factor(BirdID)), shape=23, size=3, show.legend = FALSE, colour="black") +
    scale_fill_manual(values=clrs.i$col)
  
}

#Make legends----
plot.legend.iso <- ggplot() +
  geom_sf(data=all, aes(alpha=factor(iso)), fill="black", inherit.aes = FALSE, colour="grey30", show.legend = TRUE) +
  scale_alpha_manual(values=c(0.4, 0.1), name="Isopleth", labels=c("50%", "95%")) + 
  theme(legend.position="bottom")
legend.iso <- get_legend(plot.legend.iso)

plot.legend.nest <- ggplot() +
  geom_point(data=nest1, aes(x=X, y=Y, shape=nest), colour="black", fill="white", size=4) +
  scale_shape_manual(values=c(23), name="Nest", labels=c("", "")) + 
  theme(legend.position="bottom")
legend.nest <- get_legend(plot.legend.nest)

#Plot distribution of area----
plot.area <- ggplot(area) + 
  geom_histogram(aes(x=hr95, fill=Year), colour="black", show.legend=FALSE) +
  scale_fill_manual(values=c("grey50", "grey20")) +
  xlab("95% isopleth area (ha)") +
  ylab("Number of UDs") +
  xlim(c(0,60)) +
  my.theme + 
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  ggtitle("")
#plot.area

plot.legend.year <- ggplot(area) + 
  geom_histogram(aes(x=hr95, fill=Year), colour="black", show.legend=TRUE) +
  scale_fill_manual(values=c("grey50", "grey20")) +
  xlab("95% isopleth area (ha)") +
  ylab("Number of UDs") +
  xlim(c(0,60)) +
  my.theme + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        legend.position = "bottom") + 
  ggtitle("")
legend.year <- get_legend(plot.legend.year)

#Put it together----
ggsave("Figures/Fig3UDArea.jpeg", device="jpeg", width=8, height=8, units="in", dpi=600,
       plot=grid.arrange(plots.list[[1]], plots.list[[2]], plots.list[[3]], plots.list[[4]], plots.list[[5]], plots.list[[6]], plots.list[[7]], plots.list[[8]], plots.list[[9]], plots.list[[10]],
                         plot.area,
                         legend.iso, legend.nest, legend.year,
                         widths=c(4,4,4,4),
                         heights=c(4,4,4,1),
                         layout_matrix=rbind(c(1,2,3,4),
                                             c(5,6,7,8),
                                             c(9,10,11,11),
                                             c(13,12,14,14))))

#FIGURE 4: INTERANNUAL OVERLAP####

overlap.years <- read.csv("PHRBetweenYears.csv") %>%  
  arrange(BirdID)

hr.95 <- read_sf("shapefiles/HR95.shp") %>% 
  dplyr::rename(ID=id) %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE) %>% 
  mutate(BirdID = as.numeric(BirdID)) %>% 
  left_join(sites) %>% 
  st_transform(crs=4326)

hr.50 <- read_sf("shapefiles/HR50.shp") %>% 
  dplyr::rename(ID=id) %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE) %>% 
  mutate(BirdID = as.numeric(BirdID)) %>% 
  left_join(sites) %>% 
  st_transform(crs=4326)

hr.10 <- read_sf("shapefiles/HR10.shp") %>% 
  dplyr::rename(ID=id) %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE) %>% 
  mutate(BirdID = as.numeric(BirdID)) %>% 
  left_join(sites) %>% 
  st_transform(crs=4326)

nest1 <- read.csv("NestsForTerritoryMapping.csv") %>% 
  dplyr::filter(BirdID!=1) %>% 
  mutate(Year = as.character(Year),
         nest="Nest") %>% 
  dplyr::filter(ID %in% c(hr.95$ID))

#Set scale bar sizes----
sizes <- data.frame(BirdID=unique(overlap.years$BirdID)) %>% 
  mutate(size=case_when(BirdID %in% c(21, 32) ~ 200,
                        BirdID %in% c(22) ~ 200,
                        BirdID %in% c(52) ~ 400,
                        !is.na(BirdID) ~ 100))

#Wrangling----
all <- rbind(hr.95 %>% 
               mutate(iso = 95),
             hr.50 %>% 
               mutate(iso = 50)) %>% 
  left_join(nest1) %>% 
  dplyr::filter(!is.na(area),
                BirdID %in% overlap.years$BirdID) 


#Set plot dimensions----
plot.year.list <- list()
for(i in 1:nrow(overlap.years)){
  
  bird.i <- overlap.years$BirdID[i]
  
  iso.i <- all %>% 
    dplyr::filter(BirdID==bird.i) %>% 
    dplyr::select(-NestID, -NestX, -NestY, -nest, -X, -Y) %>% 
    unique()
  
  nest.i <- all %>% 
    dplyr::filter(BirdID==bird.i) 
  
  size.i <- sizes %>% 
    dplyr::filter(BirdID==bird.i)
  
  new_bb = c(st_bbox(iso.i)[['xmin']], st_bbox(iso.i)[['ymin']]-size.i$size*0.00001,
             st_bbox(iso.i)[['xmax']], st_bbox(iso.i)[['ymax']])
  names(new_bb) = c("xmin", "ymin", "xmax", "ymax")
  attr(new_bb, "class") = "bbox"
  
  attr(st_geometry(iso.i), "bbox") = new_bb
  
  plot.i <- ggplot() +
    geom_sf(data=iso.i, aes(colour=Year, alpha=factor(iso)), lwd=1, inherit.aes = FALSE, fill="black", show.legend = FALSE) +
    scale_alpha_manual(values=c(0.4, 0.1), name="Isopleth") +
    scale_colour_manual(values=c("grey50", "grey20")) +
#    scale_colour_manual(values=diverging_hcl(2, "Blue-Yellow 3")) +
    xlab("") +
    ylab("") +
    ggtitle(paste0("Bird ", bird.i, "\n", round(overlap.years$Overlap.mn[i],3)*100, " % overlap")) +
    my.theme +
    theme(plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = ifelse(bird.i==17, "right", "none"),
          plot.title = element_text(hjust = 0.5, size=12)) +
    coord_sf(datum=NA) +
    ggsn::scalebar(data=iso.i,
                   transform=TRUE, model="WGS84",
                   dist=size.i$size, dist_unit="m",
                   box.fill=c("grey80", "grey20"),
                   box.color="grey20",
                   st.bottom=FALSE,
                   st.size=3,
                   st.dist=0.04,
                   border.size = 0.5,
                   height = 0.05,
                   facet.var="Year",
                   facet.lev="2017",
                   location="bottomleft")
  
  plot.year.list[[i]] <- plot.i + 
    geom_point(data=nest.i, aes(x=X, y=Y, fill=Year), shape=23, size=3, show.legend = FALSE, colour="black") +
    scale_fill_manual(values=c("grey50", "grey20"))
  
}

#Make legends----
plot.legend.year <- ggplot() +
  geom_sf(data=all, aes(colour=Year), lwd=2, inherit.aes = FALSE, fill="grey30", show.legend = TRUE, alpha=0.1) +
  scale_colour_manual(values=c("grey50", "grey20")) +
  theme(legend.position="bottom")
legend.year <- get_legend(plot.legend.year)

plot.legend.iso <- ggplot() +
  geom_sf(data=all, aes(alpha=factor(iso)), fill="black", inherit.aes = FALSE, colour="grey30", show.legend = TRUE) +
  scale_alpha_manual(values=c(0.4, 0.1), name="Isopleth", labels=c("50%", "95%")) +
  theme(legend.position="bottom")
legend.iso <- get_legend(plot.legend.iso)

plot.legend.nest <- ggplot() +
  geom_point(data=nest1, aes(x=X, y=Y, shape=nest), colour="black", fill="white", size=4) +
  scale_shape_manual(values=c(23), name="") +
  theme(legend.position="bottom")
legend.nest <- get_legend(plot.legend.nest)

ggsave("Figures/Fig4InterannualOverlap.jpeg", device="jpeg", width=8, height=8, units="in", dpi=600,
       plot=grid.arrange(plot.year.list[[1]], plot.year.list[[2]], plot.year.list[[3]], plot.year.list[[4]], plot.year.list[[5]], 
                         plot.year.list[[6]], plot.year.list[[7]], plot.year.list[[8]], plot.year.list[[9]], plot.year.list[[10]],
                         legend.nest,
                         legend.iso,
                         legend.year,
                         heights=c(10,10,2,2,2,2,2),
                         widths=c(10,10,10,10),
                         layout_matrix=rbind(c(1,2,3,4),
                                             c(5,6,7,8),
                                             c(9,10,NA,NA),
                                             c(9,10,11,11),
                                             c(9,10,12,12),
                                             c(9,10,13,13),
                                             c(9,10,NA,NA))))


#FIGURE 5. RSF####
pred <- read.csv("RSFPredictions.csv")

plot.rsf <- ggplot(pred) +
  geom_ribbon(aes(x=Distance, ymin=lwr, ymax=upr, group=factor(Area.rd)), alpha=0.4) +
  geom_line(aes(x=Distance, y=fit, colour=factor(Area.rd))) +
  scale_colour_manual(values=diverging_hcl(4, "Blue-Yellow 3"), name="95% isopleth\narea (ha)") +
  my.theme +
  xlab("Distance from nest (m)") +
  ylab("Relative selection probability\nfor wing-boom display") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  scale_y_continuous(limits = c(0,1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
plot.rsf

ggsave(plot.rsf, file="Figures/Fig5RSF.jpeg", device = "jpeg", height=4, width=6, units="in", dpi=600)

#APPENDICES####

#Appendix 1. Sample size----
samplesize.kde <- read.csv("KDESampleSizeBootstrapResults.csv") %>% 
  unique() %>% 
  mutate(ID=str_sub(ID, 2, 10),
         ID=gsub(ID, pattern=".201", replacement="-201")) %>% 
  left_join(boom4 %>% 
              dplyr::select(ID, n) %>% 
              rename(truen = n)) %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE) %>% 
  dplyr::filter(truen >= 5)

id.levels <- samplesize.kde %>% 
  dplyr::select(ID, truen) %>% 
  unique() %>% 
  arrange(truen) %>% 
  mutate(order = row_number())

samplesize.kde.levels <- samplesize.kde %>% 
  left_join(id.levels)

ls.fit <- read.csv("KDESampleSizeNLSPredictions.csv") %>% 
  mutate(ID=str_sub(ID, 2, 10),
         ID=gsub(ID, pattern=".201", replacement="-201")) %>% 
  left_join(id.levels) %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE)

Asym <- read.csv("KDESampleSizeNLSAsymptotes.csv") %>% 
  mutate(ID=str_sub(ID, 2, 10),
         ID=gsub(ID, pattern=".201", replacement="-201")) %>% 
  left_join(id.levels) %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE)

labs <- id.levels$ID
names(labs) <- id.levels$order

scaleFUN <- function(x) round(x, 0)

plot.n <- ggplot() +
  geom_point(data=samplesize.kde.levels, aes(x=n, y=hr95, colour=BirdID)) +
  geom_line(data=ls.fit, aes(x=n, y=r), colour="black", size=1) +
  geom_hline(data=Asym, aes(yintercept=Asym), colour="black", linetype="dashed", size=1) +
  guides(col = guide_legend(nrow=1), byrow = TRUE) +
  scale_x_continuous(labels=scaleFUN) +
  scale_y_continuous(labels=scaleFUN) +
  facet_wrap(Year~BirdID, scales="free", ncol=5) + 
  my.theme + 
  xlab("Number of wingboom points") +
  ylab("100% isopleth area (ha)") +
  theme(legend.position = "none")

ggsave(plot.n, file="Figures/Appendix2SampleSize.jpeg", device = "jpeg", height=15, width=12, units="in", dpi=600)

#DATA FOR DRYAD#####

boom3 <- read.csv("Booms3.csv") %>% 
  mutate(DateTime = ymd_hms(DateTime)) %>% 
  mutate(ID=paste0(BirdID, "-", Year))

app1 <- boom3 %>% 
  group_by(ID) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  left_join(boom3 %>%
              st_as_sf(coords=c("BoomX", "BoomY"), crs="+proj=utm +zone=12 +datum=WGS84") %>% 
              st_transform(crs=4326) %>% 
              st_coordinates() %>% 
              cbind(boom3)) %>% 
  dplyr::filter(n >= 5) %>% 
  mutate(Date=as.Date(DateTime),
         Time=str_sub(DateTime, -8, -1)) %>% 
  dplyr::select(BirdID, Date, Time, X, Y) %>% 
  rename(Longitude=X, Latitude=Y) 

write.csv(app1, "ORNITH-21-054_WingboomLocations.csv", row.names=FALSE)

nest1 <- read.csv("NestsForTerritoryMapping.csv") %>% 
  dplyr::filter(BirdID %in% app1$BirdID) %>% 
  rename(Longitude=X, Latitude=Y) %>% 
  dplyr::select(BirdID, Year, Longitude, Latitude)

write.csv(app1, "ORNITH-21-054_NestLocations.csv", row.names=FALSE)

#SUMMARY STATS####
birdscaught <- boom4 %>% 
  select(Year, BirdID, site, n, nest) %>% 
  unique() %>% 
  pivot_wider(id_cols = c(BirdID, site, nest), names_from=Year, values_from=n, names_prefix="Y", values_fill=0)

#Number of UDS estimated----
#Total birds caught
nrow(birdscaught)

#Total bird-years
nrow(birdscaught %>% dplyr::filter(Y2016 > 0)) + nrow(birdscaught %>% dplyr::filter(Y2017 > 0))

#Total bird-years with sufficient sample size
nrow(birdscaught %>% dplyr::filter(Y2016 >= 30)) + nrow(birdscaught %>% dplyr::filter(Y2017 >=30))

#Birds tracked in 2016
birdscaught %>% 
  dplyr::filter(Y2016 > 0) %>% 
  nrow()

#Birds tracked in 2016 with sufficient sample size
birdscaught %>% 
  dplyr::filter(Y2016 >= 30) %>% 
  nrow()

#Birds tracked in 2016 with sufficient sample size
birdscaught %>% 
  dplyr::filter(Y2017 >= 30) %>% 
  nrow()

#Birds recaptured and tracked in 2017
birdscaught %>% 
  dplyr::filter(Y2016 > 5, Y2017 >0) %>% 
  nrow()

#Birds recaptured and tracked in 2017 with sufficient sample size
birdscaught %>% 
  dplyr::filter(Y2016 >= 30, Y2017 >= 30) %>% 
  nrow()

#New birds caught and tracked in 2017
birdscaught %>% 
  dplyr::filter(Y2016==0) %>% 
  nrow()

birdscaught %>% 
  dplyr::filter(Y2016==0, Y2017 >= 30) %>% 
  nrow()

#Number of wingbooms per UD----
summary(boom4 %>% 
          select(Year, BirdID, site, n) %>% 
          unique() %>% 
          dplyr::filter(n >= 30))
boom4 %>% 
  select(Year, BirdID, site, n) %>% 
  unique() %>% 
  dplyr::filter(n >= 30) %>%
  summarize(sd(n))

#Number of UDs with nest----
boom4 %>% 
  select(Year, BirdID, site, n, nest) %>% 
  unique() %>% 
  dplyr::filter(n >= 30) %>%
  summarize(sum(nest))

#UD size----
#With 54 ha bird
area <- read.csv("KDEAreaSummaryAllBirds.csv") %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE)
summary(area$hr95)
sd(area$hr95)
ggplot(area) +
  geom_histogram(aes(x=hr95, colour=Year))

t.test(area$hr95 ~ area$Year)

#Without 54 ha bird
area2 <- read.csv("KDEAreaSummaryAllBirds.csv") %>% 
  separate(ID, into=c("BirdID", "Year"), remove=FALSE) %>% 
  dplyr::filter(hr95 < 40)
summary(area2$hr95)
sd(area2$hr95)

t.test(area2$hr95 ~ area2$Year)
