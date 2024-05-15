library(aws.s3)
library(mapview)

root_outpath = "/Users/ahmadjan/Core/SimulationsData/preprocessing/CAMELS_2024"
gpkg_path = glue("{root_outpath}/basins_all")
print (gpkg_path)
dir.create(gpkg_path, recursive = TRUE, showWarnings = FALSE)

setwd(gpkg_path)
wbt_wd(getwd())

###########################################################################################
# Example: Download one CAMELS basin
cat_id <- "4115265"
df = get_bucket_df(bucket = "s3://lynker-spatial/v20.1", region = 'us-west-2', max = Inf) %>% 
  filter(grepl(cat_id, Key),!grepl("[.]_", Key))


print (df$Key)

aws.s3::save_object(object =df$Key,
                    bucket = df$Bucket,
                    file = basename(df$Key),
                    region = 'us-west-2')
###########################################################################################


###########################################################################################
# Example: Download all CAMELS basins

df_all = get_bucket_df(bucket = "s3://lynker-spatial/v20.1", region = 'us-west-2', max = Inf) %>% 
  filter(grepl("Gage_", Key), !grepl("[.]_", Key))
#df_all$Key
df_all$Bucket[1]

for (key in df_all$Key){
  print (key)
  aws.s3::save_object(object = key,
                      bucket = df_all$Bucket[1],
                      file = basename(key),
                      region = 'us-west-2')
}


usa_states <- ne_states(country = "united states of america", returnclass = "sf")

usa_states50 <- usa_states[!usa_states$name %in% c("Alaska","Rhode Island", "Hawaii"), ]

camels_plot <- ggplot() +
  geom_sf(data = usa_states50, fill = "#f9f9f9", color = "black", size = 0.001) +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
  ) 


file_names <- list.files(getwd(), full.names = TRUE)

sf_objects <- lapply(file_names, read_sf, 'divides')
combined_sf <- do.call(rbind, sf_objects)

st_write(combined_sf, dsn = "test.gpkg", driver = "GPKG")

combined_512 <- read_sf("all_combined_512.gpkg")
combined_512$geom
camels_plot + geom_sf(data = combined_512$geom)#, fill = "red", color = "red") 

#combined_sf <- do.call(st_union, sf_objects)
#combined_sf
#camels_plot +
#  geom_sf(data = combined_sf$geom, fill = "blue", color = "red") 

for (fname in file_names) {
  cat_path1 <- read_sf(fname, 'divides')
  #camels_plot <- camels_plot + geom_sf(data = cat_path1$geom, fill = "red", color = "red") 
  #camels_plot <- camels_plot + geom_sf(data = cat_path1, fill = "red", color = "black") 
  camels_plot <- camels_plot + geom_sf(data = cat_path1, fill = "red", color = "red") 
  #cat_path1 <- st_union(cat_path1)
}

print (camels_plot)


cat_path1 <- read_sf(glue("{gpkg_path}/Gage_1022500.gpkg"), 'divides')

cat_path2 <- read_sf(glue("{gpkg_path}/Gage_3456500.gpkg"), 'divides')

#geom_sf(data = cat_path$geom, fill = "blue", color = "red", cex=2.75) 

ggplot() +
  geom_sf(data = vec$geom, fill = "blue", color = "red", cex=2.75) 





#+++++++++++++++++++
usa_map <- maps::map("state", interior=FALSE, col='black')
states_bdry <- maps::map("state", boundary=FALSE, col="#f9f9f9", fill= TRUE, add=TRUE) 


ggplot() +
  geom_sf(data = cat_path$geom, fill = "blue", color = "red") 

+ 
  geom_polygon (data = states_bdry,  aes(x = long, y = lat, group = group))


plot(cat_path$geom[1], add = TRUE, color = 'red', cex=20.75)
cat_path$geom[1].xmin()

plot(cat_path$geom)
mapview(cat_path, col.regions='red')



mapview(cat_path, color = "red") + states_bdry
  
maps::map("state", interior=FALSE, col='black') +
states_bdry <- maps::map("state", boundary=FALSE, col="#f9f9f9", fill= TRUE, add=TRUE)

# Plot the map of the USA with state boundaries
usa_map <- maps::map("state", interior = FALSE, col = 'black')



ggplot() +
  geom_sf(data = cat_path$geom, fill = "red", color = "red") 

# Plot the map of the USA with state boundaries and fill
states_bdry <- maps::map("state", boundary = FALSE, col = "#f9f9f9", fill = TRUE, add = TRUE) %>%
  ggplot() +
  geom_sf(data = cat_path$geom, fill = "red", color = "red") 


ggplot() +
  geom_sf(data = usa_states50, fill = "#f9f9f9", color = "black", size = 0.001) +
  geom_sf(data = cat_path$geom, fill = "blue", color = "red", cex=2.75) 

# Create a ggplot object for the map
map_plot <- ggplot() +
  geom_sf(data = cat_path$geom, fill = "blue", color = "red")  + # Add the polygon
  geom_polygon(data = fortify(states_bdry), aes(x = long, y = lat),
               fill = 'green', color = "black") 
  


