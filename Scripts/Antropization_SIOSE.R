rm(list=ls())
set.seed(10)
setwd("~")
buffdist = 250

# Librerias####
options(timeout=10000)
extrafont::loadfonts(device="win")
pacman::p_load(terra, exactextractr, dplyr, sf, readr, raster, vegan, tidyr) 


# Leer shapefile ####
SIOSE <- vect("~/GIS/España/Usos del suelo/36_Pontevedra.gpkg", layer="SAR_36_T_COMBINADA")
Leyenda_SIOSE <- read_csv2("~/GIS/España/Usos del suelo/Leyenda_SIOSE_2017_bulk.csv")
Leyenda_SIOSE$ID_COBERTURA_MAX <- Leyenda_SIOSE$ID
CRS <- crs(SIOSE)

Sites_data <- vect("D:/CloudDrive/OneDrive - Universitat de les Illes Balears/Investigación/Liquenes Prince/Zonas_Liquenes.gpkg")
Sites <- project(Sites_data, CRS)
SIOSE <- crop(SIOSE, ext(Sites)+1000)
r <- rast(ext(SIOSE))
res(r) <- 10
crs(r) <- CRS

SIOSEr <- rasterize(SIOSE, r, field="ID_COBERTURA_MAX", fun=max)
levels(SIOSEr) <- data.frame(id=Leyenda_SIOSE$ID, ANTR=Leyenda_SIOSE$Type)
names(SIOSEr)
buffers <- terra::buffer(Sites, width = buffdist)

terra::plot(SIOSEr)
terra::plot(buffers, add=T, col="blue")
terra::plot(Sites, add=T, col="red")

buffers <-  sf::st_as_sf(buffers)

# SIOSE_extracted <- terra::extract(SIOSE, buffers, fun="table")
# head(SIOSE_extracted, 20)

SIOSE <- merge(SIOSE, Leyenda_SIOSE, by = "ID_COBERTURA_MAX")
water_features <- SIOSE[SIOSE$Type == "Water"]

water_dist <- distance(Sites, water_features)
Sites$Min_Dist_to_Water <- apply(water_dist, 1, min) 

Sites$Min_Dist_to_Water
write.table(Sites$Min_Dist_to_Water, "clipboard", sep="\t", row.names=F)

# Crear un SpatRaster binario donde los valores "Water" son 1 y el resto 0

### Summarizing function
sum_cover <- function(x){
  list(x %>%
         group_by(value) %>%
         summarize(total_area = sum(coverage_area)) %>%
         mutate(proportion = total_area/sum(total_area)))
  
}

#extract the area of each raster cell covered by the plot and summarize
SIOSE_extracted <- exactextractr::exact_extract(SIOSEr, buffers, coverage_area = TRUE, 
                                                summarize_df = TRUE, fun = sum_cover, stack_apply = T)
names(SIOSE_extracted) <- buffers$Name

#merge the list elements into a df
Area <- bind_rows(SIOSE_extracted, .id = "Zonas")
Area$ID <- Area$value
Area2 <- Area %>%
  inner_join(Leyenda_SIOSE, by = join_by(ID)) %>%
  group_by(Zonas, Type) %>%
  summarize(Area = sum(proportion)) %>%
  arrange(Zonas)
Area2

write.table(Area2, "clipboard", sep="\t", row.names=F)

Diver <- Area %>% 
  group_by(Zonas) %>%
  count(ID) %>%
  pivot_wider(names_from="ID", values_from="n", values_fill=0)
DivH <- cbind(Diver[,1],"H"=diversity(Diver[,-1], index="shannon")) %>% arrange(Zonas)
DivH

write.table(DivH, "clipboard", sep="\t", row.names=F)
