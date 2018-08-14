IMR_mtr = read.table('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/5_matrix_IMR90_Coverage.txt.100.1500',header = F)
IMR_band = read.table('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/5_matrix_IMR90_Coverage.txt.100.1500.band.txt')
NHEK_mtr = read.table('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/NHEK/5_matrix_NHEK_Coverage.txt.100.1500',header = F)
NHEK_band = read.table('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/NHEK/5_matrix_NHEK_Coverage.txt.100.1500.band.txt')
GM12878_mtr = read.table('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/GM12878/5_matrix_GM12878_Coverage.txt.100.1500',header = F)
GM12878_band = read.table('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/GM12878/5_matrix_GM12878_Coverage.txt.100.1500.band.txt')
HUVEC_mtr = read.table('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HUVEC/5_matrix_HUVEC_Coverage.txt.100.1500',header = F)
HUVEC_band = read.table('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HUVEC/5_matrix_HUVEC_Coverage.txt.100.1500.band.txt')
HMEC_mtr = read.table('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HMEC/optimaze/5_matrix_HMEC_Coverage.txt.100.1500',header = F)
HMEC_band = read.table('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HMEC/optimaze/5_matrix_HMEC_Coverage.txt.100.1500.band.txt',header = F)


reg_linear <- function(matrix, band_right, band_lef, band_up, band_down){
  require(plyr)
    y = c()
    x = c()
    n = 1
    for(i in band_up:band_down){
        right = max(band_right,i)
        for(j in right:band_lef){
          if(matrix[i,j]!=0 & (j - i) !=0){
            y[n] = log(matrix[i, j])
            x[n] = log(j - i)
            n = n+1
          }
        }
    }
    
    pnt = data.frame(x,y)
    pnt$group = as.numeric(cut(pnt[,1],breaks = 30))
    pnt2 = ddply(pnt, "group", numcolwise(mean))
    
    return(pnt2)
}

band = as.vector(GM12878_band[2,])
pnt2 = reg_linear(GM12878_mtr, as.numeric(band[2]), as.numeric(band[3]), as.numeric(band[2]), as.numeric(band[3]))
band = as.vector(GM12878_band[1,])
pnt = reg_linear(GM12878_mtr, as.numeric(band[2]), as.numeric(band[3]), as.numeric(band[2]), as.numeric(band[3]))
pnt3 = reg_linear(GM12878_mtr, 66, 106, 66, 106)
pnt4 = reg_linear(GM12878_mtr, 106, 441, 66, 106)


plot(pnt2[[1]],pnt2[[2]])
points(pnt[[1]],pnt[[2]],col='red')
points(pnt3[[1]],pnt3[[2]],col='blue')
points(pnt4[[1]],pnt4[[2]],col='green')



pnt = reg_linear(IMR_mtr, 110, 280, 110, 280)
pnt2 = reg_linear(IMR_mtr, 280, 400, 110, 280)
pnt3 = reg_linear(IMR_mtr, 66, 110, 66, 110)
pnt4 = reg_linear(IMR_mtr, 110, 400, 66, 110)


plot(pnt[[1]],pnt[[2]])
points(pnt2[[1]],pnt2[[2]],col='red')
points(pnt3[[1]],pnt3[[2]],col='blue')
points(pnt4[[1]],pnt4[[2]],col='green')


pnt_IMR = reg_linear(IMR_mtr, 110, 280, 110, 280)
pnt_HMEC = reg_linear(HMEC_mtr, 110, 280, 110, 280)
pnt2_IMR = reg_linear(IMR_mtr, 261, 388, 105, 261)
pnt2_HMEC = reg_linear(HMEC_mtr, 261, 388, 105, 261)
pnt3_HMEC = reg_linear(HMEC_mtr, 105, 261, 60, 105)



plot(pnt2_IMR[,2],pnt2_IMR[,3],ylim=c(0,6))
points(pnt_IMR[,2],pnt_IMR[,3],col='blue')
points(pnt_HMEC[,2],pnt_HMEC[,3],col='green')
points(pnt2_HMEC[,2],pnt2_HMEC[,3],col='red')

lmMod_IMR = lm(pnt_IMR[,3]~pnt_IMR[,2])
lmMod_IMR2 = lm(pnt2_IMR[,3]~pnt2_IMR[,2])
summary (lmMod_IMR) # 4.738460, -0.656097 
summary (lmMod_IMR2) # 4.470579, -0.644222
lmMod_IMR_pred = predict(lmMod_IMR)
lmMod_IMR_pred2 = predict(lmMod_IMR2)



lmMod_HMEC = lm(pnt_HMEC[pnt_HMEC[,2]>1,3]~pnt_HMEC[pnt_HMEC[,2]>1,2])
lmMod_HMEC2 = lm(pnt2_HMEC[pnt2_HMEC[,2]>1,3]~pnt2_HMEC[pnt2_HMEC[,2]>1,2])
lmMod_HMEC3 = lm(pnt3_HMEC[pnt3_HMEC[,2]>1,3]~pnt3_HMEC[pnt3_HMEC[,2]>1,2])
summary (lmMod_HMEC) # 3.46083, -0.62116 
summary (lmMod_HMEC2) # 1.45658, -0.20872
summary (lmMod_HMEC3) # 2.95210, -0.53374
lmMod_HMEC_pred = predict(lmMod_HMEC)
lmMod_HMEC_pred2 = predict(lmMod_HMEC2)
lmMod_HMEC_pred3 = predict(lmMod_HMEC3)


s = smooth.spline(pnt_HMEC[,2],lmMod_HMEC_pred)
p_s = predict(s)

plot(pnt2_IMR[,2],pnt2_IMR[,3],ylim=c(0,6))
points(pnt_IMR[,2],pnt_IMR[,3],col='blue')
points(pnt_HMEC[,2],pnt_HMEC[,3],col='green')
points(pnt2_HMEC[,2],pnt2_HMEC[,3],col='red')
points(pnt3_HMEC[,2],pnt3_HMEC[,3],col='blue')
lines(pnt_IMR[,2],lmMod_IMR_pred)
lines(pnt2_IMR[,2],lmMod_IMR_pred2)
lines(pnt_HMEC[pnt_HMEC[,2]>1,2],lmMod_HMEC_pred)
lines(pnt2_HMEC[pnt2_HMEC[,2]>1,2],lmMod_HMEC_pred2)
lines(pnt3_HMEC[pnt3_HMEC[,2]>1,2],lmMod_HMEC_pred3)

as.numeric(cut(pnt_IMR[,1],breaks = 30))

