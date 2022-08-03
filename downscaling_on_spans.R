source('scripts/routines_functions.R')

## load shape files 

#poly = readOGR("GIS_Polimi/projected/vento/projected_poly.shp")
poly = readOGR("GIS_Polimi/projected/vento/original_crs_poly.shp")

crs_set <- proj4string(poly)

## load the mapping pixels-spans
mapping = read_excel('input/ID_campate_TO_vento.xls', sheet = 1)
## mapping ha 32220 rows. In totale ci sono però 29439 campate. mapping non ha duplicati (ID,STRNO) -> alcune campate sono a cavallo di due pixels. elimino i duplicati -> devo per forza avere corrispondenza 1:1 tra campata e pixel.. soluzione ottimale: associare una campata al pixel cui appartiene il suo punto medio..TODO in futuro (maybe)
## delle 29439 campate noi però consideriamo solo le 28684 campate per cui sono stati calcolati i tempi di ritorno
length(unique(mapping$STRNO))
map_id_span <- merge(poly@data[,c('Id', 'PC1', 'PC2')],mapping[, c('Id','STRNO', 'STRNO_SUP')], by = 'Id')
map_id_span <- map_id_span[!duplicated(map_id_span$STRNO),] ## ora ho 29439 campate associate al relativo pixel Id
head(map_id_span)
map_id_span <- map_id_span[order(map_id_span$STRNO),]

### c'era inconsistenza tra i nomi delle campate nel file kml che mi ha mandato FF (lat_lon) rispetto al file 'ID_campate_TO_vento.xls'... risolto in lat_lon_spans.ipynb

spans = read.table('csv_files/spans_lat_lon.csv')
spans$lat1 <- as.numeric(spans$lat1)
spans$lat2 <- as.numeric(spans$lat2)
spans$lon1 <- as.numeric(spans$lon1)
spans$lon2 <- as.numeric(spans$lon2)

#head(spans)

## spans contiene i nomi delle campate in esame -> elimino da map_id_span le campate extra
map_id_span <- map_id_span[which(map_id_span$STRNO %in% unique(spans$campata)),]
# ok ora ho 28684 campate in map_id_span

## convert spans into a SpatialLinesDataFrame

st_segment = function(r){st_linestring(t(matrix(unlist(r), 2, 2)))}
df = spans[, c(1,8)]
df$geometry = st_sfc(sapply(1:nrow(spans), function(i){st_segment(spans[i,c(3,2,6,5)])},simplify=FALSE))
p_lines = st_sf(df)
sp_line = as(p_lines, "Spatial")
tmp_df <- p_lines
tmp_df$geometry <- NULL
tmp_df <- as.data.frame(tmp_df)

## create the SpatialLinesDataFrame
sp_lns_dfr <- sp::SpatialLinesDataFrame(sp_line, data = tmp_df)
sp_lns_dfr@proj4string <- CRS(crs_set)
## project in the correct crs
crs=CRS('+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs')
sp_lns_dfr_t <- spTransform(sp_lns_dfr, CRS('+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs'))
#plot(sp_lns_dfr_t)
tmp <- st_as_sf(sp_lns_dfr_t)

# ottengo i punti medi di tutte le campate 
spans_mid <- st_centroid(tmp)
to_sim <- spans_mid
to_sim <- as(spans_mid, "Spatial")

## controllo che non ci siano campate con nomi diversi ma stessa geometry (i.e stessa localizzazione)
check_coords <- coordinates(sp_lns_dfr_t )
doubles <- zerodist(to_sim, zero = 0.0)

## invece ci sono duplicati.. Ad es:
#check_coords[2]
#check_coords[86]

to_sim2 <- remove.duplicates(to_sim)  ## se ci sono più punti a distanza 0 (i.e duplicati) kriging non funziona : problemi di singolarità

# tengo traccia delle campate rimosse
removed = to_sim[which(!to_sim@data$campata %in% to_sim2@data$campata),]

## salvo le coppie di campate con nomi diversi ma identiche coordinate geografiche

save1 <- st_drop_geometry(tmp[doubles[,1],])$campata
save2 <- st_drop_geometry(tmp[doubles[,2],])$campata
save3 <- st_drop_geometry(tmp[doubles[,1],])$linea
save4 <- st_drop_geometry(tmp[doubles[,2],])$linea
final_save <- data.frame(cbind('campata1'=save1,'campata2' = save2, 'linea1' = save3, 'linea2' = save4))
write.csv(final_save, "extra/geometries_duplicates.csv")


## load the mapping pixels-spans

sp_lns_dfr@data = merge(sp_lns_dfr@data, map_id_span, by.x = 'campata', by.y = 'STRNO')

sp_lns_dfr_t@data = merge(sp_lns_dfr_t@data, map_id_span, by.x = 'campata', by.y = 'STRNO')

######################################################################################################################################################
########### First method for downscaling relying mainly on the atakrig package  ######################################################################
######################################################################################################################################################

####################################
### DOWNSCALING 4X4 Km -> campata ###
####################################


### 1. project on the correct CRS (UMT 32)
crs=CRS('+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs')

poly_t <- spTransform(poly, crs)


### 2. employ atakrig for variogram deconvolution 

# Compute empirical variograms (as a reference to choose the variogram model).
emp_vg1 = variogram(PC1 ~ 1, data = poly_t) 
emp_vg2 = variogram(PC2 ~ 1, data = poly_t) 

plot(emp_vg1)
plot(emp_vg2)


g_emp_vg1vg2 <- gstat(id = "PC1", formula = PC1~1, data = poly_t)
g_emp_vg1vg2 <- gstat(g_emp_vg1vg2, id = "PC2", formula = PC2~1, data = poly_t)
emp_vg1vg2=variogram(g_emp_vg1vg2)
quartz()
plot(emp_vg1vg2)


semi_variograms_plot<-data.frame(emp_vg1vg2$gamma[1:15],emp_vg1vg2$gamma[16:30],emp_vg1vg2$gamma[31:45])
semi_variograms_plot_x<-data.frame(emp_vg1vg2$dist[1:15],emp_vg1vg2$dist[16:30],emp_vg1vg2$dist[31:45])
quartz()
par(cex.lab=1.5,cex.main=1.5,cex.axis=1.5)
matplot(semi_variograms_plot_x,semi_variograms_plot,col=c("red","blue","green"),pch=20,ylab = "Semi-variograms of residuals",xlab = "Lag (m)")
legend("topleft", legend=c('Semi-cross-variogram','Semi-variogram of PC2', 'Semi-variogram of PC1'), 
       col = c("red","blue","green"), lty = c(2,1,1), cex = .8)
grid(col = 'cornsilk2')
quartz.save('gamma_1-2.png', type = "png")


v1 = fit.variogram(emp_vg1, vgm(model = "Exp"))
quartz()
plot(emp_vg1, v1, xlab=TeX("$h$"),ylab=TeX("$\\gamma(h)$")) # Exponential
v2 = fit.variogram(emp_vg2, vgm(model = "Exp"))
quartz()
plot(emp_vg2, v2) # Exponential

v12 = fit.variogram(emp_vg1vg2, vgm(model = "Sph"))
quartz()
plot(emp_vg1vg2, v12) # Exponential

# Perform variogram deconvolution.

discr_poly1 = discretizePolygon(poly_t,cellsize = 500,value = "PC1",showProgressBar = T)
discr_poly2 = discretizePolygon(poly_t,cellsize = 500,value = "PC2",showProgressBar = T)

discr_tot = list(PC1 = discr_poly1, PC2 = discr_poly2)

# independent deconvolutions

deconv1 = deconvPointVgm(discr_poly1,model = "Sph")

deconv2 = deconvPointVgm(discr_poly2,model = "Sph")

# Cross-deconvolution
deconv = deconvPointVgmForCoKriging(discr_tot, model = "Sph", fig = F)
plotDeconvVgm(deconv)


# 3. perform downscaling simulations in the centroids of each span
## check that there are no duplicates in to_sim

#PC1
vgm1 <- extractPointVgm(deconv1)
PC1_downscaled <- krige(PC1~1, locations = poly_t, newdata = to_sim2, model = vgm1, nmax = 10, nsim=500, beta = 0)

#PC2
vgm2 <- extractPointVgm(deconv2)
PC2_downscaled <- krige(PC2~1, locations = poly_t, newdata = to_sim2, model = vgm2, nmax = 10, nsim=500, beta = 0)

### plot PC1 and PC2 on original resolution 
p_poly = st_as_sf(sp_lns_dfr_t)
class(p_poly)

g1 = p_poly %>% ggplot() + geom_sf(aes(color = PC1), lwd = 1) + scale_color_viridis(option = "inferno") + 
  ggtitle("PC1 scores (original resolution)")  + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + theme(legend.position="bottom") + rremove("xlab")  + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")

g2 = p_poly %>% ggplot() + geom_sf(aes(color = PC2), lwd = 1) + scale_color_viridis(option = "inferno") + 
  ggtitle("PC2 scores (original resolution)")  + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + theme(legend.position="bottom") + rremove("xlab")  + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")

ggarrange(g1,g2,common.legend = TRUE, legend = 'bottom')


## arrange data for plots and save

tmp1 = cbind(PC1_downscaled@data, to_sim2@data)
PC1_df = merge(p_poly,tmp1, by = 'campata' )

tmp2 = cbind(PC2_downscaled@data, to_sim2@data)
PC2_df = merge(p_poly,tmp2, by = 'campata' )

final_save <- final_save[which(!duplicated(final_save$campata1)),]
merged_tmp <- merge(removed, final_save[,c(1,2)], by.x = 'campata', by.y = 'campata2', all.x = TRUE, all.y = FALSE)
tmp_PC1 <- PC1_df
tmp_PC1$geometry <- NULL
merged_tmp2 <- merge(merged_tmp,tmp_PC1, by.x = 'campata1', by.y = 'campata' )
merged_tmp2$campata1 <- NULL
merged_tmp2$linea.x <- merged_tmp2$linea
merged_tmp2$linea.y <- merged_tmp2$linea
merged_tmp2$linea <- NULL
merged_tmp2 <- st_as_sf(merged_tmp2)
PC1_df <- rbind(PC1_df, merged_tmp2)

final_save <- final_save[which(!duplicated(final_save$campata1)),]
merged_tmp <- merge(removed, final_save[,c(1,2)], by.x = 'campata', by.y = 'campata2', all.x = TRUE, all.y = FALSE)
tmp_PC2 <- PC2_df
tmp_PC2$geometry <- NULL
merged_tmp2 <- merge(merged_tmp,tmp_PC2, by.x = 'campata1', by.y = 'campata' )
merged_tmp2$campata1 <- NULL
merged_tmp2$linea.x <- merged_tmp2$linea
merged_tmp2$linea.y <- merged_tmp2$linea
merged_tmp2$linea <- NULL
merged_tmp2 <- st_as_sf(merged_tmp2)
PC2_df <- rbind(PC2_df, merged_tmp2)

## PC1 plots

g1_ds1 = PC1_df %>% ggplot() + geom_sf(aes(color = sim1), lwd = 1) + scale_color_viridis(option = "inferno") + ggtitle("PC1 scores on each span") + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")
g1_ds2 = PC1_df %>% ggplot() + geom_sf(aes(color = sim2), lwd = 1) + scale_color_viridis(option = "inferno") + ggtitle("PC1 scores on each span") + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")
g1_ds3 = PC1_df %>% ggplot() + geom_sf(aes(color = sim3), lwd = 1) + scale_color_viridis(option = "inferno") + ggtitle("PC1 scores on each span") + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")
g1_ds4 = PC1_df %>% ggplot() + geom_sf(aes(color = sim4), lwd = 1) + scale_color_viridis(option = "inferno") + ggtitle("PC1 scores on each span") + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")
g1_ds5 = PC1_df %>% ggplot() + geom_sf(aes(color = sim5), lwd = 1) + scale_color_viridis(option = "inferno") + ggtitle("PC1 scores on each span") + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")

ggarrange(g1, g1_ds1,g1_ds2,g1_ds3, ncol = 2, nrow = 2, common.legend = TRUE, legend = 'bottom')



## PC2 plots



g2_ds1 =  PC2_df %>% ggplot() + geom_sf(aes(color = sim1), lwd = 1) + scale_color_viridis(option = "inferno") + ggtitle("PC1 scores on each span") + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")
g2_ds2 =  PC2_df %>% ggplot() + geom_sf(aes(color = sim2), lwd = 1) + scale_color_viridis(option = "inferno") + ggtitle("PC1 scores on each span") + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")
g2_ds3 =  PC2_df %>% ggplot() + geom_sf(aes(color = sim3), lwd = 1) + scale_color_viridis(option = "inferno") + ggtitle("PC1 scores on each span") + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")
g2_ds4 =  PC2_df %>% ggplot() + geom_sf(aes(color = sim4), lwd = 1) + scale_color_viridis(option = "inferno") + ggtitle("PC1 scores on each span") + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")
g2_ds5 =  PC2_df %>% ggplot() + geom_sf(aes(color = sim5), lwd = 1) + scale_color_viridis(option = "inferno") + ggtitle("PC1 scores on each span") + clean_theme()+theme_pubclean(base_size = 10) + labs_pubr(base_size = 10) + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")

ggarrange(g2, g2_ds1,g2_ds2,g2_ds3, ncol = 2, nrow = 2, common.legend = TRUE, legend = 'bottom')

## save downscaling
save(PC1_df,PC2_df, poly, poly_t, file = "Rdata/downscales_predictions_spans.Rdata")









