file_url <- 'http://users.stat.umn.edu/~sandy/alr3ed/website/data/UN2.txt'
datos <- read.table(file_url, T)
rownames(datos) <- datos$Locality
datos <- datos[,-4]

y <- datos[,'logFertility']
X <- cbind(1, datos[,c('logPPgdp', 'Purban')])
X <- as.matrix(X)

bh  <- solve(t(X) %*% X) %*% t(X) %*% y; round(bh, 5)

s2h <- t(y - X%*%bh) %*% (y - X%*%bh) / (193 - 3)
s2h <- as.numeric(s2h); round(s2h, 5)

vbh <- s2h * solve(t(X) %*% X); round(vbh, 6)

gamma1 <- qchisq(0.025, 193 - 3)
gamma2 <- qchisq(0.975, 193 - 3)
s2_LL <- (193 - 3) * s2h / gamma2; s2_LL

s2_UL <- (193 - 3) * s2h / gamma1; s2_UL

b0_LL <- bh[1] - qt(0.975, 93 - 3) * sqrt(vbh[1,1]); b0_LL

b0_UL <- bh[1] + qt(0.975, 93 - 3) * sqrt(vbh[1,1]); b0_UL

t0 <- abs(bh[1]/sqrt(vbh[1,1])); t0 # Estadístico

tt <- qt(0.975, 93 - 3); tt         # Cuantil de la distribución t  

b1_LL <- bh[2] - qt(0.975, 93 - 3) * sqrt(vbh[2,2]); b1_LL
b1_UL <- bh[2] + qt(0.975, 93 - 3) * sqrt(vbh[2,2]); b1_UL

t1 <- abs(bh[2]/sqrt(vbh[2,2])); t1 # Estadístico
tt <- qt(0.975, 93 - 3); tt         # Cuantil de la distribución t  
#para beta_2
b2_LL <- bh[3] - qt(0.975, 93 - 3) * sqrt(vbh[3,3]); b2_LL
b2_UL <- bh[3] + qt(0.975, 93 - 3) * sqrt(vbh[3,3]); b2_UL

t2 <- abs(bh[3]/sqrt(vbh[3,3])); t2 # Estadístico
tt <- qt(0.975, 93 - 3); tt         # Cuantil de la distribución t  

modelo1 <- lm(logFertility ~ logPPgdp, datos)
summary(modelo1)

modelo2 <- lm(Purban ~ logPPgdp, datos)
summary(modelo2)

res_logfert <- residuals(modelo1)
res_Purban <- residuals(modelo2)
modelo3 <- lm(res_logfert ~ res_Purban)
summary(modelo3)

plot(res_Purban, res_logfert, pch = 16, col = 'steelblue',
     xlab = 'Purban | logPPgdp', ylab = 'logFertility | logPPgdp')

#Análisis de la varianza
sc_tc <- sum((y - mean(y))^2);    sc_tc
sc_error <- sum((y - X%*%bh)^2);  sc_error
sc_reg <- sc_tc - sc_error;       sc_reg

#región de confianza para b1 y b2
region2 <- function(x, y, bh, X){
  b0 <- matrix(c(x, y), 2, 1)
  b <- bh[2:3]
  aux <- solve(t(X)%*%X)
  m <- solve(aux[2:3, 2:3])
  res <- t(b - b0)%*%m%*%(b-b0)
  return(res)
}

seqx <- seq(bh[2] - 0.06, bh[2] + 0.06, length.out = 250)
seqy <- seq(bh[3]-0.006, bh[3]+0.006, length.out = 250)

m_aux <- matrix(0, 250, 250)
for (i in 1:250){
  for (j in 1:250){
    m_aux[i, j] <- region2(seqx[i], seqy[j], bh, X)
  }
}

confianza <- c(2 * s2h * qf(0.90, 2, 190),
               2 * s2h * qf(0.95, 2, 190),
               2 * s2h * qf(0.99, 2, 190))
image(seqx, seqy, m_aux, breaks = c(0, confianza), xlim = c(-0.2, 0),
      col = c('steelblue3', 'steelblue2', 'steelblue1'),
      xlab = 'beta1', ylab = 'beta2')
contour(seqx, seqy, m_aux, levels = confianza, drawlabels = F, add = T)
legend('topright', legend = c('90%', '95%', '99%'), pch = 16,
       col = c('steelblue3', 'steelblue2', 'steelblue1'), cex = 1.5)

#pruebas simultáneas+
#beta0
b0_LL_B <- bh[1] - qt(1 - 0.05/6, 193 - 3) * sqrt(vbh[1,1]); b0_LL_B
b0_UL_B <- bh[1] + qt(1 - 0.05/6, 193 - 3) * sqrt(vbh[1,1]); b0_UL_B
#beta1
b1_LL_B <- bh[2] - qt(1 - 0.05/6, 193 - 3) * sqrt(vbh[2,2]); b1_LL_B
b1_UL_B <- bh[2] + qt(1 - 0.05/6, 193 - 3) * sqrt(vbh[2,2]); b1_UL_B
#beta2
b2_LL_B <- bh[3] - qt(1 - 0.05/6, 193 - 3) * sqrt(vbh[3,3]); b2_LL_B
b2_UL_B <- bh[3] + qt(1 - 0.05/6, 193 - 3) * sqrt(vbh[3,3]); b2_UL_B

