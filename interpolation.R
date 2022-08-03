source('scripts/routines_functions.R')

## range per curve condizionate

if(work_with_cond == 1){
start = 40
stop = 140} else{
  start = 0
  stop = 200
}


## dataset griglia-valori superamento soglia

vento_df <- read_excel('input/ID_celle_vento_TO_MI.xls', sheet = 2)
vento_df <- as.data.frame(vento_df)
ids <- vento_df[,1]     ##id cells
N_samples <- length(ids)  ## number of pixels
#matplot(c(40,60,80,100,120,140), t(vento_df[,-1]), type="l", xlab = 'soglia vento', ylab = 'Probabilità', main = 'CCDF')  ## plot raw data

## modifico dominio [0,200] imponendo 1 in 0 e 0 in 200 
dom <- c(0,40,60,80,100,120,140,200)
mdf_vento <- cbind('0' = rep(1,dim(vento_df)[1]), vento_df[,-1],'200' = rep(0,dim(vento_df)[1]) )
matplot(dom, t(mdf_vento), type="l", xlab = 'soglia vento', ylab = 'Probabilità'  ,main = 'CCDF in [0-200]Km/h')
abline(v = 40, lty = 2)
abline(v = 140, lty = 2)

## calcolo la cdf facendo 1-ccdf
cdf <- 1 - mdf_vento

matplot(dom, t(cdf), type="l", xlab = 'soglia vento', ylab = 'probability'  ,main = 'CDF in [0-200]Km/h')
abline(v = 40, lty = 2)
abline(v = 140, lty = 2)

## smoothing con Bernstein Polynomials in the whole domain [0,200] 

# sqeeze domain between 0 and 1
## [0,200] -> [0,1]

x_eval <- seq(min_dom,max_dom,step)
x_sq <- squeeze(dom, dom)
x_eval_sq <- squeeze(x_eval, dom)

# 1. select order of polynomial minimizing SSE
m = seq(10,140,10)

## smoothing per tutte le cdf
#select the number of polynomials to employ 

SSE <- matrix(nrow = length(m), ncol = N_samples )
for(i in 1:length(m)){
  curr_m = m[i]
  
  for( j in 1:N_samples){
    tmp <- data.frame(cbind(t = x_sq, val =array(t(cdf)[,j])))
    total_data <- smoothing(x_sq,tmp, curr_m)
    SSE[i,j] <- sum((total_data-tmp$val)^2)
  }
}

rownames(SSE) <- m
rowMeans(SSE)
boxplot(as.data.frame(t(SSE)) , col='gold', main = 'SSE', xlab = 'number of basis', ylab = "SSE")
abline(h = 0.01, col = 'red', lty = 2)

# selected m
m = 40

# 2. obtain the smoothed cdfs
tmp <- data.frame(cbind(t = x_sq, val =array(t(cdf)[,1])))
total_data <- data.frame(smoothing(x_eval_sq,tmp, m))

for( i in 2:dim(t(cdf))[2]){
  tmp <- data.frame(cbind(t = x_sq, val =array(t(cdf)[,i])))
  total_data <- cbind(total_data, smoothing(x_eval_sq,tmp, m))
}

colnames(total_data) <- ids

matplot(x_eval, total_data, type="l", xlab = 'soglia vento [Km/h]', ylab = 'probability'  ,main = 'Smoothed CDF in [0-200]Km/h')

# save 
write.csv(total_data,"csv_files/smoothed_cdf.csv", row.names = FALSE)

## 3. obtain the smoothing for the pdfs

tmp <- data.frame(cbind(t = x_sq, val =array(t(cdf)[,1])))
total_data_pdf <- data.frame(pdf_from_smoothing(x_eval_sq,tmp, m))

for( i in 2:dim(t(cdf))[2]){
  tmp <- data.frame(cbind(t = x_sq, val =array(t(cdf)[,i])))
  total_data_pdf <- cbind(total_data_pdf, pdf_from_smoothing(x_eval_sq,tmp, m))
}

colnames(total_data_pdf) <- ids

# normalize the densities to 1 integral
prova <- normalize_density(total_data_pdf, step )
#matplot(x_eval, total_data_pdf, type="l", lty = 1, xlab = 'soglia vento', ylab = ''  ,main = 'PDF in [0-200]Km/h')
colnames(prova) <- ids
matplot(x_eval, prova, type="l", lty = 1, xlab = 'soglia vento', ylab = ''  ,main = 'PDF in [0-200]Km/h')

# save densities (IDs of pixels are column names)
write.csv(prova,"csv_files/smoothed_pdf.csv", row.names = FALSE)

################ conditoning pdfs on the domain [start,stop]

tmp_data <- data.frame(cbind(t = x_eval, val = prova[,1]))
df_conditioned <- data.frame(conditioning(tmp_data,start,stop))

for( i in 2:dim(prova)[2]){
  tmp_data <- data.frame(cbind(t = x_eval, val = prova[,i]))
  df_conditioned <- cbind(df_conditioned, conditioning(tmp_data,start,stop))
}

matplot(seq(start,stop,step), df_conditioned, type="l", lty = 1, xlab = 'soglia vento', ylab = ''  ,main = 'conditioned PDF in [40-140]Km/h')
colnames(df_conditioned) <- ids

# save conditioned densities
write.csv(df_conditioned,"csv_files/smoothed_pdf_conditioned.csv", row.names = FALSE)


