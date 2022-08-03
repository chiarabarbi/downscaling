source('scripts/routines_functions.R')
load('RData/FPCA_eigen.RData')
load('RData/downscales_predictions_spans.RData')

####

linee <- unique(PC1_df$linea.x)


store_results <- list()
## level 1 -> linea
## level 2 -> campata 
## level 3 -> "scores_sim", "curve_sim", "avg_curve_std","original_scores", "original_curve"


for( l in 1:length(linee)){
  tmp0 <- list()
  campate <- unique( st_drop_geometry(PC1_df[which(PC1_df$linea.x == linee[l]),c(1)]))
  
  for( cc in 1:dim(campate)[1]){
    tmp = list() 
    
    tmp[['scores_sim']] <- data.frame(t(rbind('PC1' = st_drop_geometry(PC1_df[which(PC1_df$campata == campate[cc,] ),-c(1,2,3,4,5,6,507)]),'PC2' = st_drop_geometry(PC2_df[which(PC2_df$campata == campate[cc,] ),-c(1,2,3,4,5,6,507)])))) 
    scr <- t(rbind('PC1' = st_drop_geometry(PC1_df[which(PC1_df$campata == campate[cc,] ),-c(1,2,3,4,5,6,507)]),'PC2' = st_drop_geometry(PC2_df[which(PC2_df$campata == campate[cc,] ),-c(1,2,3,4,5,6,507)])))
    data_sim <- m.s+pc.fc$vec[,1:K]%*%t(scr)
    densitiesK_sim =clr2density.mv(data_sim, t_cond, t_step) 
    tmp[['curve_sim']] <- densitiesK_sim
    m.sim =apply(densitiesK_sim, 1, mean)
    sd.sim = apply(densitiesK_sim, 1, sd)
    tmp[['avg_curve_std']] <- data.frame(t(rbind("mean" = m.sim, "sd"=sd.sim)))
    tmp[['original_scores']] <- st_drop_geometry(PC1_df[which(PC1_df$campata == campate[cc,] ), c(4,5)])
    #tmp['original_curve'] <- pdf_originals[which(pdf_originals$STRNO == campate[cc,]),-c(1,2,3,4,5,6)]
    
    tmp0[[campate[cc,]]] <- tmp
  }
  store_results[[linee[l]]] <-tmp0
}

save(store_results, file = "RData/Downscaling_results.Rdata")



