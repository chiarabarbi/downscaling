
source('scripts/routines_functions.R')

#______________________________________________________________________________________
###  LOAD DATA  

densities <- read.csv('csv_files/smoothed_pdf.csv', header=TRUE)
densities_cond <- read.csv('csv_files/smoothed_pdf_conditioned.csv', header = TRUE)


if(work_with_cond == 1){
  start = 40
  stop = 140
}else{
  start = 0
  stop = 180  ## uso 180 come max per evitare problemi numerici (code quasi tutte nulle..)
}

t = seq(min_dom,max_dom,step)  ## define the whole domain (step=1, defined in routines_functions.R)
t_step=t[2]-t[1]
ids <- colnames(densities)

matplot(t,densities,type="l",lty=1,xlab="", ylab="density",cex.lab=1.2,col=grey(0.1+(1:21)/30))


# Clr-transforms 
## clr-transform è invariante al dominio considerato (sia che se lavoro sulle condizionate che sulle totali)

t_cond <- t[which(t >= start & t <=stop )]

data=(clr.mv(densities[which(t >= start & t <=stop),], t_cond, t_cond[2]-t_cond[1]))
colnames(data)<- ids

#dim(data)


drop <- which(is.na(data), arr.ind = TRUE)[,2]
if( length(drop) != 0) {
  data <- data[, -drop]
}

# Define data dimensions
N_samples=dim(data)[2]
n=dim(data)[1]

# Check zero-integral condition
int=rep(0,N_samples)
for(i in 1:N_samples){
  int[i]=trapzc(t[2]-t[1], data[,i])}
int

# Plot clr-transforms
matplot(t_cond,data, type='l',lty=1,
        ylab='clr transformed pdf', xlab='soglia [Km/h]', main = 'CLR transformed PDFs' )

#______________________________________________________________________________________

### SFPCA 
# SFPCA on data
lmod=data
lmod=lmod*sqrt(t_step)
lmod[1,]=lmod[1,]/sqrt(2)
lmod[n,]=lmod[n,]/sqrt(2)


SS=(N_samples-1)/N_samples*cov(t(lmod))
S=SS

pc.fc=eigen(S)
pc.fc$vec[1,]=pc.fc$vec[1,]*sqrt(2)
pc.fc$vec[n,]=pc.fc$vec[n,]*sqrt(2)
pc.fc$vec=pc.fc$vec/sqrt(t_step)

K=Nmax.harm=10

# scores
m.s=apply(data, 1, mean)

data.c=NULL
for(i in 1:N_samples)
  data.c=cbind(data.c, data[,i]-m.s)

matplot(t_cond,data.c, type='l', lty=1,main = 'centered clr densities', xlab = 'soglia [Km/h]',ylab='centered clr densities')

sc.a=matrix(NA, ncol=Nmax.harm, nrow=N_samples)
for(i in 1:N_samples)
{
  for(j in 1:Nmax.harm)
    sc.a[i,j]=trapzc(t_step, data.c[,i]*pc.fc$vec[,j])
}

boxplot(sc.a, las=1, col='gold', main='Principal Components')

# Scatterplot of scores
pairs(sc.a[, 1:4], pch=19)
plot(sc.a[,1], rep(0, N_samples), col='transparent', xlab='Score PC1', ylab='')
text(sc.a[,1], rep(0, N_samples), 1:N_samples)

# projected curves
K=2  ## 2 principal components considered

dataK=m.s+pc.fc$vec[,1:K]%*%t(sc.a[,1:K])
densitiesK=clr2density.mv(dataK, t_cond, t_step)

# saving results for further analysis
save(t_cond, K, t_step,m.s,pc.fc, file = "Rdata/FPCA_eigen.Rdata")


#---- Plots results ----
matplot(t_cond,data, type='l', lty=1, ylab='clr-density',col=grey(0.1+(0:21)/30))
matplot(t_cond,dataK, type='l', lty=1, ylab='clr-density (proxy)', 
        col=grey(0.1+(0:21)/30))
title(main=paste0('K = ',K))

matplot(t_cond, densities_cond, type='l', lty=1, ylab='density',col=grey(0.1+(0:21)/30), main = 'Original Data')
matplot(t_cond,densitiesK, type='l', lty=1, ylab='density (proxy)', col=grey(0.1+(0:21)/30))
title(main=paste0('K = ',K))

par(mfrow = c(1,K), cex = 1.5)

totvar = sum(pc.fc$values)

for(k in 1:K){
  sdev = sqrt(pc.fc$values[k])
  perc <- round(sdev^2/totvar*100,1)
  matplot(t_cond, clr2density(m.s, t_cond, t_step), type='l', lty=1, lwd = 2, ylab='density', ylim = c(0,0.05), col = 'red', xlab = paste0(perc,"% explained variance"))
  matlines(t_cond, clr2density(m.s+2*sqrt(pc.fc$values[k])*pc.fc$vectors[,k], t_cond, t_step),lwd = 2, col='gold1' )
  matlines(t_cond, clr2density(m.s-2*sqrt(pc.fc$values[k])*pc.fc$vectors[,k], 
                          t_cond, t_step), lwd = 2,col='purple4' )
  title(main=paste0('k = ',k))
}



scores <- sc.a[, 1:k]

scores <- data.frame(cbind('Id'=substring(colnames(data), 2),scores))
colnames(scores) <- c('Id', 'PC1', 'PC2')
scores$PC1=as.numeric(scores$PC1)
scores$PC2=as.numeric(scores$PC2)

if(work_with_cond == 1){
write.csv(scores,"csv_files/scores_conditioned.csv", row.names = FALSE)}else{
  write.csv(scores,"csv_files/scores_unconditioned.csv", row.names = FALSE)
}


## spatial plots

poly = readOGR("GIS_Polimi/vento/griglia_vento.shp")

poly = st_as_sf(poly)
tail(poly)
poly = merge( poly, scores, by = 'Id')
tail(poly)

poly_t = as_Spatial(poly)

plot(poly_t)

p_poly = st_as_sf(poly_t)

g1 = p_poly %>% ggplot() + geom_sf(aes(fill = PC1), lwd = 0) + scale_fill_viridis(option = "inferno") + 
  ggtitle("PC 1")  + clean_theme()+theme_pubclean(base_size = 20) + labs_pubr(base_size = 20) + theme(legend.position="bottom") + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")
g2 = p_poly %>% ggplot() + geom_sf(aes(fill = PC2), lwd = 0) + scale_fill_viridis(option = "inferno") + 
  ggtitle("PC 2")  + clean_theme()+theme_pubclean(base_size = 20) + labs_pubr(base_size = 20) + theme(legend.position="bottom") + rremove("xlab") + rremove("ylab") + rremove("xy.text") + rremove("ticks") + rremove("legend.title")

ggarrange(g1, g2, ncol = 2, nrow = 1, common.legend = TRUE, legend = 'bottom')
