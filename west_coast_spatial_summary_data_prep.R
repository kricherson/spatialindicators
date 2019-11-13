library(broom)
library(tidyverse)
library(ggforce) # for plotting ellipses
library(gridExtra)
library(geosphere)
library(viridis)
library(janitor)

source("Spatial_indicators_functions_Woillez2009_modified.R")

source("~/observer/Input/load_data_2019-08-09.R") 
load_data(c("WCGOP","ASHOP"))

source("~/observer/Input/Richerson/Spatial_indicators_functions_Woillez2009_modified.r")

#Port coordinates
ports <- read_csv("~/observer/Input/Richerson/ports.csv") %>% 
  filter(IOPCID != "Whidbey Island ports") %>% #This one is wierd, just use other ONP lat/lon
  clean_names() %>% 
  dplyr::select(pcid, pcid_lat, pcid_long)
#Note: missing coords for SFA --  make same as SF? 

out_drive<-("~/observer/Output/Richerson other/Groundfish spatial patterns/")

#What species are we interested in? Including Lewis' suggestions plus a few more high rev species 
sel_species <- c("Arrowtooth Flounder", "Big Skate", "Dover Sole", "English Sole", "Lingcod" ,"Longnose Skate","Longspine Thornyhead", "Pacific Cod", "Pacific Ocean Perch", "Pacific Hake", "Petrale Sole", "Rex Sole", "Sablefish", "Shortspine Thornyhead", "Spiny Dogfish Shark")

#Get haul-level data for selected species. Only look at retained? Note: should we include MWT for hake/MWR? Including for now.
ob_bt_sel <- OBOrig_Proc %>% 
  clean_names() %>% 
  filter(sector %in% c("Catch Shares", "Limited Entry Trawl") & gear %in% c("Bottom Trawl", "Midwater Trawl")) %>% 
  filter(datatype == "Analysis Data") %>% 
  filter(species %in% sel_species & catch_disposition == "R") %>% 
  left_join(ports) %>% 
  group_by(year, drvid, species, haul_id, pcid, pcid_lat, pcid_long, avg_lat, avg_long) %>% 
  summarise(mlon = sum(avg_long*mt, na.rm=T)/sum(mt[which(!is.na(avg_long))],na.rm=T),
            mlat = sum(avg_lat*mt, na.rm=T)/sum(mt[which(!is.na(avg_lat))],na.rm=T),
            mt = sum(mt)) %>% 
  rowwise() %>% 
  mutate(dist = distm(c(avg_long, avg_lat), c(pcid_long, pcid_lat), fun=distHaversine)) %>% 
  dplyr::filter(mlon!=0 & mlat!=0)

#Let's also look at total retained, discarded, just for fun. Put in same form as above for easier plotting.
ob_bt_disret <- OBOrig_Proc %>% 
  clean_names() %>% 
  filter(sector %in% c("Catch Shares", "Limited Entry Trawl") & gear == "Bottom Trawl") %>% 
  filter(datatype == "Analysis Data") %>%
  left_join(ports) %>% 
  group_by(year, drvid,  species, haul_id, pcid, pcid_lat, pcid_long, avg_lat, avg_long) %>% 
  summarise(all_retained = sum(ret_mt),
            all_discard = sum(dis_mt)) %>% 
  gather(species, mt, all_retained:all_discard) %>% 
  group_by(year, species, haul_id, pcid, pcid_lat, pcid_long, avg_lat, avg_long) %>% 
  summarise(mlon = sum(avg_long*mt, na.rm=T)/sum(mt[which(!is.na(avg_long))],na.rm=T),
            mlat = sum(avg_lat*mt, na.rm=T)/sum(mt[which(!is.na(avg_lat))],na.rm=T),
            mt = sum(mt)) %>% 
  rowwise() %>% 
  mutate(dist = distm(c(avg_long, avg_lat), c(pcid_long, pcid_lat), fun=distHaversine)) %>% 
  dplyr::filter(mlon!=0 & mlat!=0)

bt_dat <- bind_rows(ob_bt_sel, ob_bt_disret) %>% 
  mutate(dist_km = dist/1000)

save(bt_dat, file = paste0(out_drive, "bt_dat_ret.rda"))

#Calculate CGI
bt_cgi <- bt_dat %>% 
  group_by(year, species) %>% 
  do((cgi(x=.$mlon,y=.$mlat) %>% data.frame)) %>% 
  mutate(disp_species = gsub(" ", "\n", species)) %>% 
  data.frame

save(bt_cgi, file = paste0(out_drive, "bt_cgi_ret.rda"))

#also do ashop
#First make haul-level DF of PWHT catches
ashop <- ASOrig_Proc %>% 
  clean_names() %>% 
  group_by(year, drvid, haul_id, avg_lat, avg_long) %>% 
  summarise(mlon = sum(avg_long*pwht_mt, na.rm=T)/sum(pwht_mt[which(!is.na(avg_long))],na.rm=T),
            mlat = sum(avg_lat*pwht_mt, na.rm=T)/sum(pwht_mt[which(!is.na(avg_lat))],na.rm=T),
            mt = sum(pwht_mt)) %>% 
  dplyr::filter(mlon!=0 & mlat!=0) 

ashop_cgi <- ashop %>% 
  group_by(year) %>% 
  do((cgi(x=.$mlon,y=.$mlat) %>% data.frame)) %>% 
  mutate(species = "Pacific hake",
         disp_species = "At-sea\nhake") %>% 
  data.frame

save(ashop_cgi, file = paste0(out_drive, "ashop_cgi.rda"))

#combine CGIs
wc_cgi <- bt_cgi %>% 
  bind_rows(ashop_cgi)

save(wc_cgi, file = paste0(out_drive, "wc_cgi.rda"))


####get GIC
#WCGOP GIC
bt_gic <- bt_dat %>% 
  as.data.frame() %>% 
  #filter(species == "Big Skate") %>% 
  #filter(year>2008) %>% 
  group_by(species, year) %>% 
  summarise(gic=gic(x1=avg_long,
                    y1=avg_lat,
                    z1=mt,
                    x2=bt_dat$avg_long[bt_dat$species==species],
                    y2=bt_dat$avg_lat[bt_dat$species==species],
                    z2=bt_dat$mt[bt_dat$species==species])) %>% 
  mutate(disp_species = gsub(" ", "\n", species))


save(bt_gic, file = paste0(out_drive,"bt_gic.rda"))

#ASHOP GIC
#Need a haul-level DF to start with 
ashop_gic <- ashop %>% 
  as.data.frame() %>% 
  group_by(year) %>% 
  summarise(gic=gic(x1=avg_long,
                    y1=avg_lat,
                    z1=mt,
                    x2=ashop$avg_long,
                    y2=ashop$avg_lat,
                    z2=ashop$mt)) %>% 
  mutate(species = "Pacific hake",
         disp_species = "At-sea\nhake") %>% 
  data.frame

save(ashop_gic, file = "ashop_gic.rda")

#combine
wc_gic <- bind_rows(bt_gic, ashop_gic)

save(wc_gic, file = paste0(out_drive, "wc_gic.rda"))