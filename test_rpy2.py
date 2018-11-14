import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

r = robjects.r

r('''
library("changepoint")
edge_row <- function(l){
    m.binseg <- cpt.meanvar(l,test.stat = "Normal", method = "BinSeg")
    change = cpts(m.binseg)
    if(length(change)>2){
        change = change
    }
    return(change)
}
        
bound_v2 <- function(y,fold){
    require(bcp)
    if(fold <= 1){
    p_lim = 0.6
    }
    else{p_lim = 0.4}
    
    #y = band2
    #hist(y,breaks = 500)
    #plot(freq[,1],freq[,3],type = 'h')
    freq = data.frame(table(unlist(y)))
    freq[,1] <- as.numeric(as.character(freq[,1]))
    freq[,2] <- as.numeric(as.character(freq[,2]))
    freq[,3] <- freq[,2]/max(freq[,2])
    
    ff = floor(length(freq[,3])/5000)
    point = freq[1,]
    i = 1
    for(i in 1:(ff+1)){
    f0 = freq[(1+(i-1)*5000):min((i*5000),length(freq[,3])),]
    bcp.1 <- bcp(f0[,3],mcmc = 5000)
    #bcp.1 <- bcp(freq[(1+(i-1)*5000):(i*5000),3],mcmc = 5000)
    #plot(bcp.1,xlim = c(0,200), main="Univariate Change Point Example")
    pre.p = summary(bcp.1)
    point = as.matrix(f0[pre.p[,1]>0.2,])
    #point = rbind(as.matrix(point),as.matrix(f0[pre.p[,1]>p_lim,]))
    
    point = point[!is.na(point[,1]),]
    }
    if(!is.null(nrow(point))){
    return(point[-1,])
    }
    else{
    return(point)
    }
}

find_edge <- function(points,cut){
  require('cluster')
  points = as.data.frame(points)
  d=dist(points[,1],method = "euclidean")
  if(nrow(points)>0){
    cluster = hclust(d,method = "complete")
    #cut=50
    clusterCut <- cutree(cluster, h = cut)
    points[,4] = clusterCut
    pk1 = c()
    j=1
    for(i in unique(points[,4])){
      m = max(points[points[,4]==i,2])
      pk1[j]=points[points[,2]==m&points[,4]==i,1]
      j=j+1
    }
    return(pk1)
  }
  else{
    return(points)
  }
  
}

processFile_v2 = function(filepath,cover_limup,cover_limdown,fold) {
  if(fold<=1){
    range = 500
    range2 = 150
    range3 = 20
  }
  else{
    range = 500/fold
    range2 = 150/fold
    range3 = 20/fold
  }
  
  bps = data.frame(0,0,0,0,0,0)
  con = file(filepath, "r")
  line = readLines(con, n = 1)
  l = as.numeric(as.character(strsplit(line,'\t')[[1]]))
  end = length(l) - 50
  i = 1
  while ( i<end ) {
    line = readLines(con, n = 1)
    l = as.numeric(as.character(strsplit(line,'\t')[[1]]))
    l1=l[i:min((i+range),length(l))]
    #l1[l1>10]=10
    l1[l1>cover_limup]=cover_limup
    l1[l1<=cover_limdown]=0
    l2=l[i:max((i-range),1)]
    #l2[l2>10]=10
    l2[l2>cover_limup]=cover_limup
    l2[l2<=cover_limdown]=0
    if ( length(line) == 0 ) {
      break
    }
    else{
      if(i>50){
        band1 = edge_row(l1)
        #print(band1)
        band2 = edge_row(l2)
        #print(band2)
        if(length(band1)>0&length(band2)>0){
          bps[i,1] = i
          bps[i,2] = i + band1[2]
          bps[i,3] = i - band2[2]
          bps[i,4] = i + band1[3]
          bps[i,5] = i - band2[3]
        }
      }
    }
    i = i+1
  }
  bps[,1] = 1:length(bps[,1])
  bandary_col = bps[(bps[,2]-bps[,1])<range2&(bps[,2]-bps[,1])>range3,2]
  bandary_col = bandary_col[!is.na(bandary_col)]
  bandary_row = bps[(bps[,1]-bps[,3])<range2&(bps[,1]-bps[,3])>range3,3]
  bandary_row = bandary_row[!is.na(bandary_row)]
  close(con)
  lis = list(bps,bandary_col,bandary_row)
  return(lis)
}

find_bound<-function(band2,band3,band1,dist_lim,fold){
  require(MASS) 
  require(spatstat)
  #band2 = b1
  #band3 = b2
  #band1 = b5
  #dist_lim = 25
  
  prob_col = bound_v2(band2,fold)
  if(FALSE){
    return(data.frame(0,0,0))
  }
  else{
    if(is.null(nrow(prob_col)) |nrow(prob_col)==0){
      return(data.frame(0,0,0))
    }
    else{
      prob_col2 = prob_col[order(prob_col[,1]),]
      
      edg_col = find_edge(prob_col2,dist_lim)
      
      prob_row = bound_v2(band3,fold)
      if(!is.null(nrow(prob_row))){
        if(nrow(prob_row)>0){
          prob_row2 = prob_row[order(prob_row[,1]),]
          edg_row = find_edge(prob_row2,dist_lim)
          
          #hist(band3[1:700],breaks = 200)
          
          pnts = data.frame(0,0)
          t=1
          for(i in 1:length(edg_col)){
            for(j in 1:length(edg_row)){
              if(edg_col[i]-edg_row[j]>10&edg_col[i]-edg_row[j]<1000){
                pnts[t,1]=edg_col[i]
                pnts[t,2]=edg_row[j]
                t=t+1
              }
            }
          }
          
          cord = band1
          #cord = cord[1:lim,] # test
          cord = na.omit(cord)
          ex=data.frame()
          
          t=1
          for(t in 1:ceiling(length(cord[,1])/3000)){
            n1=1+3000*(t-1)
            n2=3000*t
            x = na.omit(cord[n1:min(n2,length(cord[,1])),])
            xy.kde <- kde2d(x[,1],x[,2], n=(n2-n1)/3,h=20)
            xy.im <- im(t(xy.kde$z), xcol=xy.kde$x, yrow=xy.kde$y) # Allows interpolation $
            #filled.contour(xy.kde$z,levels=c(0,0.00000123,0.00000726,1),col=c('#FFFFFF','#C0FFC0','#80FF80'))
            #image(xy.kde$z,col = topo.colors(12))
            #image(xy.kde,col = topo.colors(12))
            #points(pnts[,2],pnts[,1],col='red')
            #points(bandary[,2],bandary[,3],col='black')
            test_pnt = pnts[pnts[,2]<xy.im$xrange[2]&pnts[,2]>xy.im$xrange[1],]
            test_pnt = test_pnt[test_pnt[,1]<xy.im$yrange[2]&test_pnt[,1]>xy.im$yrange[1],]
            c.0 <- sum(xy.kde$z)
            ex_trmp = data.frame()
            if(length(test_pnt[,1])>0){
              #i=1
              for(i in 1:length(test_pnt[,1])){
                z <- interp.im(xy.im, test_pnt[i,2], test_pnt[i,1])                            # Densities at the probe points
                #z <- interp.im(xy.im, squr1, squr2)                            # Densities at the probe points
                #squr1 = rep((test_pnt[i,2]-20):(test_pnt[i,2]+20),41)
                #squr2 = rep((test_pnt[i,1]-20):(test_pnt[i,1]+20),each = 41)
                #squr = data.frame(squr2,squr1)
                #p <- sum(interp.im(xy.kde$z,squr1,squr2))
                #z <- interp.im(xy.im, (test_pnt[i,2]-10):(test_pnt[i,2]+10), rep(test_pnt[i,1],21))  #######
                p <- sapply(z, function(a) sum(xy.kde$z[xy.kde$z < a])) / c.0
                ex_trmp[i,1] = p
                ex_trmp[i,2] = test_pnt[i,2]
                ex_trmp[i,3] = test_pnt[i,1]
              }
              ex = rbind(ex,ex_trmp)
            }
          }
          if(fold<=1){
            #ex[,4] = ex[,1]/(ex[,3]-ex[,2])
            bandary = ex[ex[,1]>0.8,]
            #bandary = ex[ex[,4]>0.001,]
          }
          else(bandary = ex[ex[,1]>quantile(ex[,1],0.8),])
          return(bandary)
        }
        else{
          return(data.frame(0,0,0))
        }
      }
      else{
        return(data.frame(0,0,0))
      }
    }
  }
}

main_bound_v2 <- function(file, color_limup, color_limdown, dist_lim,fold){
  require(MASS)     # kde2d
  require(spatstat) 
  band = processFile_v2(file, color_limup,color_limdown, fold)
  band2 = band[[2]]
  band3 = band[[3]]
  if(length(band3)==0|length(band2)==0){
    return(t(matrix(c(0,0,0))))
  }
  else{
    band4 = band[[1]][,c(3,2)]
    band5 = band[[1]][,c(5,4)]
    bandary = t(matrix(c(1,1,1)))
    for(i in 1:max(1,floor(length(band2)/2000))){
      b1 = band2[(1+(i-1)*2000):min((i*2000),length(band2))]
      b2 = band3[(1+(i-1)*2000):min((i*2000),length(band2))]
      b1 = na.omit(b1)
      b2 = na.omit(b2)
      b3 = band4[min(b1,b2):max(b1,b2),]
      b4 = band5[min(b1,b2):max(b1,b2),]
      b5 = rbind(as.matrix(b3),as.matrix(b3),as.matrix(b4))
      b5 = b5[order(b5[,1]),]
      bandary_temp = find_bound(b1,b2,b5,dist_lim,fold)
      print(i)
      bandary = rbind(bandary,as.matrix(bandary_temp))
    }
    return(bandary[-1,])
  }
}
        ''')

TAD = r['main_bound_v2']

chr = '13'
down = '9000'
up = '10000'

input = r.paste('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HUVEC/', chr,
                '_matrix_HUVEC_Coverage.txt.', down, '.', up, sep='')

bb3 = TAD(input,10,1,25,1)

print(bb3)

