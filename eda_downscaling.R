source('scripts/routines_functions.R')
load('RData/FPCA_eigen.RData')
load('RData/downscales_predictions_spans.RData')

## carico la mappa pixel-campata
K = 2
t = seq(40,140,1)
map_id <- data.frame(read_excel('input/ID_campate_TO_vento.xls', sheet = 1))

## 6 linee considerate 
linee <- c("22223A1", "22288A1", "23841D1", "23842B1", "23876B1")


## fisso una linea 
linea <- "22223A1"

## ottengo le campate in questa linea 
campate_df <- map_id[which(map_id$STRNO_SUP == linea),] 


## ottengo le pdf (smoothate) associate alle campate della linea 
## associo al ciascuna campata della linea l'id pixel associato (quindi associo a ciascuna campata della linea la pdf associata)

densities_cond <- read.csv('csv_files/smoothed_pdf_conditioned.csv', header = TRUE)
ids <- unlist(lapply(colnames(densities_cond), FUN = function(i) substr(i, 2, 100)))
colnames(densities_cond) <- ids
pdf_df <- data.frame(t(densities_cond))
pdf_df$Id <- ids
pdf_originals <- merge(campate_df, pdf_df, by = 'Id')

## plot delle pdf non downscaled
matplot(t,t(pdf_originals[-c(1,2,3,4,5,6)]),type="l",lty=1,xlab = "soglia vento [Km/h]", ylab="density",main = "PDFs originali associate alle campate della linea",cex.lab=1.2,  ylim = c(0,0.03))

### ottengo le curve medie per tutte le campate in questa linea
linea
linea_camp <- pdf_originals[which(pdf_originals$STRNO_SUP == linea),]

campate_in_pix <- unique(linea_camp$STRNO)
sim_df_PC1 <- st_drop_geometry(PC1_df[which(PC1_df$campata %in% campate_in_pix), -c(1,2,3,4,5,6, 507, 508)])
sim_df_PC2 <- st_drop_geometry(PC2_df[which(PC2_df$campata %in% campate_in_pix), -c(1,2,3,4,5,6, 507, 508)])

matplot(t,t(linea_camp[,-c(1,2,3,4,5,6)]),type="l",lty=1,lwd = 1,xlab = "soglia vento [Km/h]", ylab="", main = paste("PDFs associate alle", length(campate_in_pix), "campate della linea", linea, "(no downscaling)", sep = " "),cex.lab=1.2,  ylim = c(0,0.05))

mean_densities <- c()

for( i in 1:length(campate_in_pix)){
  sim_df <- rbind(sim_df_PC1[i,], sim_df_PC2[i,])
  sim_df <- t(sim_df)
  colnames(sim_df) <- c('PC1', 'PC2')
  data_sim <- m.s+pc.fc$vec[,1:K]%*%t(sim_df)
  densitiesK_sim =clr2density.mv(data_sim, t_cond, t_step) 
  m.sim =apply(densitiesK_sim, 1, mean)
  mean_densities <- cbind(mean_densities, m.sim)
  sd.sim = apply(densitiesK_sim, 1, sd)
  #matplot(t_cond,cbind(densitiesK_sim,m.sim, m.sim + 3*sd.sim, m.sim - 3*sd.sim ), col = c(rep('grey', dim(densitiesK_sim)[2]),'black', 'black', 'black'), type='l', lwd = c(rep(1,500), 2, 1, 1),  lty=c(rep(1, 501), 2, 2), ylim = c(0,0.05),xlab = "soglia vento [Km/h]", ylab=" ",main = campate_in_pix[i])
  
  
  #matplot(t_cond,cbind(t(pdf_originals[which(pdf_originals$STRNO == camp),-c(1,2,3,4,5,6)]),m.sim, m.sim + 3*sd.sim, m.sim - 3*sd.sim ), col = c('red','black', 'black', 'black'), type='l', lwd = c(2,2, 1, 1),  lty=c(1,1, 2, 2), ylim = c(0,0.05),xlab = "soglia vento [Km/h]", ylab=" ",main = campate_in_pix[i])
  
}

matplot(t_cond, mean_densities, type='l', lwd = 1,  lty=1, ylim = c(0,0.05),xlab = "soglia vento [Km/h]", ylab=" ", main = paste("PDFs associate alle", length(campate_in_pix), "campate della linea", linea, "(downscaling)", sep = " "))












### focus su una singola campata (plot delle 500 simulazioni )

camp <- "22223A1-077-----"
sim_df_PC1 <- st_drop_geometry(PC1_df[which(PC1_df$campata == camp), -c(1,2,3,4,5,6, 507, 508)])
sim_df_PC2 <- st_drop_geometry(PC2_df[which(PC2_df$campata == camp), -c(1,2,3,4,5,6, 507, 508)])
sim_df <- rbind(sim_df_PC1, sim_df_PC2)
sim_df <- t(sim_df)
colnames(sim_df) <- c('PC1', 'PC2')

matplot(t,t(pdf_originals[which(pdf_originals$STRNO == camp),-c(1,2,3,4,5,6)]),type="l",lty=1,lwd = 2, col = 'red',xlab = "soglia vento [Km/h]", ylab="",main = "PDFs associata alla campata 22223A1-077-----",cex.lab=1.2,  ylim = c(0,0.05))


data_sim <- m.s+pc.fc$vec[,1:K]%*%t(sim_df)
densitiesK_sim =clr2density.mv(data_sim, t_cond, t_step) 

m.sim =apply(densitiesK_sim, 1, mean)
sd.sim = apply(densitiesK_sim, 1, sd)

matplot(t_cond,cbind(densitiesK_sim,m.sim, m.sim + 3*sd.sim, m.sim - 3*sd.sim ), col = c(rep('grey', dim(densitiesK_sim)[2]),'black', 'black', 'black'), type='l', lwd = c(rep(1,500), 2, 1, 1),  lty=c(rep(1, 501), 2, 2), ylim = c(0,0.05),xlab = "soglia vento [Km/h]", ylab=" ",main = "Simulazione di 500 PDFs (downscaling) per la campata 22223A1-077-----")


matplot(t_cond,cbind(t(pdf_originals[which(pdf_originals$STRNO == camp),-c(1,2,3,4,5,6)]),m.sim, m.sim + 3*sd.sim, m.sim - 3*sd.sim ), col = c('red','black', 'black', 'black'), type='l', lwd = c(2,2, 1, 1),  lty=c(1,1, 2, 2), ylim = c(0,0.05),xlab = "soglia vento [Km/h]", ylab=" ",main = "PDF originale, media delle simulazioni ed incertezza per la campata 22223A1-077-----")


### ottengo le curve medie per tutte le campate in questo pixel 
id <- pdf_originals[which(pdf_originals$STRNO == camp),c(1)]
pixel_camp <- pdf_originals[which(pdf_originals$Id == id),]

matplot(t,t(pdf_originals[which(pdf_originals$Id == id),-c(1,2,3,4,5,6)]),type="l",lty=1,lwd = 2, col = 'red',xlab = "soglia vento [Km/h]", ylab="",main = "PDFs associate alle 12 campate nel pixel ID = 516",cex.lab=1.2,  ylim = c(0,0.05))

campate_in_pix <- unique(pixel_camp$STRNO)
sim_df_PC1 <- st_drop_geometry(PC1_df[which(PC1_df$campata %in% campate_in_pix), -c(1,2,3,4,5,6, 507, 508)])
sim_df_PC2 <- st_drop_geometry(PC2_df[which(PC2_df$campata %in% campate_in_pix), -c(1,2,3,4,5,6, 507, 508)])


mean_densities <- c()

for( i in 1:length(campate_in_pix)){
  sim_df <- rbind(sim_df_PC1[i,], sim_df_PC2[i,])
  sim_df <- t(sim_df)
  colnames(sim_df) <- c('PC1', 'PC2')
  data_sim <- m.s+pc.fc$vec[,1:K]%*%t(sim_df)
  densitiesK_sim =clr2density.mv(data_sim, t_cond, t_step) 
  m.sim =apply(densitiesK_sim, 1, mean)
  mean_densities <- cbind(mean_densities, m.sim)
  sd.sim = apply(densitiesK_sim, 1, sd)
  #matplot(t_cond,cbind(densitiesK_sim,m.sim, m.sim + 3*sd.sim, m.sim - 3*sd.sim ), col = c(rep('grey', dim(densitiesK_sim)[2]),'black', 'black', 'black'), type='l', lwd = c(rep(1,500), 2, 1, 1),  lty=c(rep(1, 501), 2, 2), ylim = c(0,0.05),xlab = "soglia vento [Km/h]", ylab=" ",main = campate_in_pix[i])
  
  
  matplot(t_cond,cbind(t(pdf_originals[which(pdf_originals$STRNO == camp),-c(1,2,3,4,5,6)]),m.sim, m.sim + 3*sd.sim, m.sim - 3*sd.sim ), col = c('red','black', 'black', 'black'), type='l', lwd = c(2,2, 1, 1),  lty=c(1,1, 2, 2), ylim = c(0,0.05),xlab = "soglia vento [Km/h]", ylab=" ",main = campate_in_pix[i])
  
  }

matplot(t_cond, mean_densities, col = c('grey60'), type='l', lwd = 1,  lty=1, ylim = c(0,0.05),xlab = "soglia vento [Km/h]", ylab=" ", main = "PDFs associate alle 12 campate del pixel ID=516 (downscaling)")




