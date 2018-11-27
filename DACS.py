import rpy2.robjects as robjects
from itertools import combinations
import sys
import argparse
import cv2
import math
import os
import warnings
from rpy2.rinterface import RRuntimeWarning
import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri
import numpy as np
from scipy.sparse.linalg import eigsh
import imagehash
from PIL import Image


def printHelp():
    print('\nDACTAD version 1.0.0')
    print('For help information for each function, try:\npython3 DACTAD.py <function> -h')
    print('\nFunctions:')
    print('\tTAD_calculator:\n\t\tidentify the topological domain\n')
    print('\tcorner_split:\n\t\tcorner split algorithm for identifying split TAD\n')
    print('\tTAD_similarity:\n\t\tfour algorithms for calculating TAD similarity\n')
    #print('\nKaifu Chen, et al. chenkaifu@gmail.com, Li lab, Biostatistics department, Dan L. Duncan cancer center, Baylor College of Medicine.')
    print('')


def TAD_calculator(command='TAD_calculator'):
    '''
    Description:
        This function provide an entrance to identify TAD.
    parameters:
        none
    '''

    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython3 DACS.py TAD_calculator <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python DACS.py TAD_calculator -h\n")
        return 0

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython3 DACS.py TAD_calculator "
                                           "[optional arguments]\n\n", description='')
    parser.add_argument('command', default=None,
                        help="set as 'TAD_calculator' to identify TAD relative to Hi-C contact map")

    parser.add_argument('-c', '--contact_map', dest="contact_map", default=None,
                        help="path to Hi-C contact map")

    parser.add_argument('-u', '--up_limit', dest="up", default=None, type=float,
                        help="up limit for Hi-C contact map")

    parser.add_argument('-d', '--down_limit', dest="down", default=0, type=float,
                        help="down limit for Hi-C contact map")

    parser.add_argument('-o', '--TAD_output', dest="output", default=None,
                        help="path and filename for the output file of TADs")

    parser.add_argument('-p', '--TAD_plot', dest="plot", default=1, type=int,
                        help="Set to 1 to plot the contact map and TAD, else set to 0 to cancel this analysis.")

    args = parser.parse_args()

    head, tail = os.path.split(args.output)
    base = os.path.splitext(tail)
    if (not os.path.exists(head)) and head != '':
        print('path does not exist\n')
    else:
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        r = robjects.r
        r('''
        options(warn=-1)
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
            options(warn=-1)
            require(bcp)
            if(fold <= 1){
            p_lim = 0.6
            }
            else{p_lim = 0.4}
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
            pre.p = data.frame(bcp.1$posterior.prob, bcp.1$posterior.mean[,1])                       
            colnames(pre.p) = c('Probability', 'X1')
            point = as.matrix(f0[pre.p[,1]>0.2,])
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
          options(warn=-1)
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
            l1[l1>cover_limup]=cover_limup
            l1[l1<=cover_limdown]=0
            l2=l[i:max((i-range),1)]
            l2[l2>cover_limup]=cover_limup
            l2[l2<=cover_limdown]=0
            if ( length(line) == 0 ) {
              break
            }
            else{
              if(i>50){
                band1 = edge_row(l1)
                band2 = edge_row(l2)
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
          options(warn=-1)
          require(MASS) 
          require(spatstat)
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
                  cord = na.omit(cord)
                  ex=data.frame()
    
                  t=1
                  for(t in 1:ceiling(length(cord[,1])/3000)){
                    n1=1+3000*(t-1)
                    n2=3000*t
                    x = na.omit(cord[n1:min(n2,length(cord[,1])),])
                    xy.kde <- kde2d(x[,1],x[,2], n=(n2-n1)/3,h=20)
                    xy.im <- im(t(xy.kde$z), xcol=xy.kde$x, yrow=xy.kde$y) # Allows interpolation $
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
                    bandary = ex[ex[,1]>0.8,]
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
    
        main_bound_v2 <- function(file, color_limup, color_limdown, dist_lim = 25, fold = 1){
          options(warn=-1)
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
              bandary = rbind(bandary,as.matrix(bandary_temp))
            }
            return(bandary[-1,])
          }
        }
                ''')
        TAD_calling = r['main_bound_v2']

        TAD = np.array(TAD_calling(args.contact_map, args.up, args.down))
        np.savetxt(args.output, TAD, delimiter="\t")

        def TAD_plot(f1, down, up):
            '''
            :param f1: contact map
            :param down: down limit for color key
            :param up: up limit for color key
            :return:
            '''
            img = np.loadtxt(f1)
            img[img > up] = up
            img[img < down] = 0
            img = img * (255 / up)
            img = img.astype(np.uint8)
            color_img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
            color_img[:, :, 0] = 255 - color_img[:, :, 0]
            color_img[:, :, 1] = 255 - color_img[:, :, 1]
            color_img[:, :, 2] = 255
            for x in TAD:
                if float(x[0]) > 0.8:
                    color_img = cv2.rectangle(color_img, (math.ceil(float(x[1])), math.ceil(float(x[1]))),
                                              (math.ceil(float(x[2])), math.ceil(float(x[2]))), (0, 0, 0), 2,
                                              lineType=4)
            cv2.imwrite(os.path.join(head, base[0]+'.tiff'), color_img)
        if args.plot == 1:
            TAD_plot(args.contact_map, args.down, args.up)


def corner_split(command='corner_split'):
    '''
    corner split algorithm for identifying split TAD
    '''

    def region_detect(TAD1, TAD2):
        flag1 = []
        flag2 = []
        band_s = np.array([0, 0, 0])
        try:
            for i in range(0, TAD1.shape[0]):
                for j in range(0, TAD2.shape[0]):
                    if abs(TAD1[i, 2] - TAD2[j, 2]) <= 10 and abs(TAD1[i, 1] - TAD2[j, 1]) <= 10:
                        band_s = np.vstack([band_s, (TAD1[i, :] + TAD2[j, :]) / 2])
                        flag1.append(i)
                        flag2.append(j)
            TAD1_only = np.delete(TAD1, flag1, axis=0)
            TAD2_only = np.delete(TAD2, flag2, axis=0)
        except:
            try:
                TAD2 = TAD2.reshape([1, 3])
                for i in range(0, TAD1.shape[0]):
                    for j in range(0, TAD2.shape[0]):
                        if abs(TAD1[i, 2] - TAD2[j, 2]) <= 10 and abs(TAD1[i, 1] - TAD2[j, 1]) <= 10:
                            band_s = np.vstack([band_s, (TAD1[i, :] + TAD2[j, :]) / 2])
                            flag1.append(i)
                            flag2.append(j)
                TAD1_only = np.delete(TAD1, flag1, axis=0)
                TAD2_only = np.delete(TAD2, flag2, axis=0)
            except:
                TAD1 = TAD1.reshape([1, 3])
                for i in range(0, TAD1.shape[0]):
                    for j in range(0, TAD2.shape[0]):
                        if abs(TAD1[i, 2] - TAD2[j, 2]) <= 10 and abs(TAD1[i, 1] - TAD2[j, 1]) <= 10:
                            band_s = np.vstack([band_s, (TAD1[i, :] + TAD2[j, :]) / 2])
                            flag1.append(i)
                            flag2.append(j)
                TAD1_only = np.delete(TAD1, flag1, axis=0)
                TAD2_only = np.delete(TAD2, flag2, axis=0)
        try:
            return band_s[1:, :], TAD1_only, TAD2_only
        except:
            return band_s, TAD1_only, TAD2_only

    def RangInList(L,TAD):
        flag1 = 1
        for t in TAD:
            if abs(L[1] - t[1]) < 10 and abs(L[2] - t[2]) < 10:
                flag1 = 0
        return flag1

    def TAD_divid_V3(range1, range2, TAD_s):
        flag = 0
        for t in range2:
            if abs(t[1] - range1[1]) < 10:
                flag = flag + 1
            if abs(t[2] - range1[2]) < 10:
                flag = flag + 1
        if flag >= 2:
            dic = []
            for i in range(2, range2.shape[0] + 1):
                comb = list(combinations(range(range2.shape[0]), i))
                for j, t in enumerate(comb):
                    x = range2[t, :]
                    length = 0
                    flag3 = []  # all divide region are similar region
                    for s in x:
                        length = length + s[2] - s[1]
                        f2 = RangInList(s, TAD_s.tolist())
                        flag3.append(f2)
                    r = length / (range1[2] - range1[1])
                    if r > 0.8 and r < 1.1:  # change parameter
                        flag2 = 0  # find all real sep
                        for l in range(x.shape[0] - 1):
                            if x[l + 1, 1] - x[l, 2] > 20 or x[l, 2] - x[l + 1, 1] > 20:  # change parameter
                                flag2 = 1
                        if flag2 == 0:
                            dic.append([x, r])
            return (dic)

    def region_divid_v3(TAD1, TAD2, TAD_s):
        '''
        :param TAD1: whole region
        :param TAD2: divide region
        :param TAD_s: similar region
        :return: D1: TAD of unit region, D2: divide TAD with coverage rate, D3: divide TAD without coverage rate
        '''
        TAD_in = np.empty(shape=[0, 3])
        diff1 = []
        diff2 = []
        diff3 = []
        count = 0
        try:
            for i, t1 in enumerate(TAD1):
                for t2 in TAD2:
                    if t1[1] - t2[1] < 20 and t2[2] - t1[2] < 20:
                        TAD_in = np.append(TAD_in, [t2], axis=0)
                D = TAD_divid_V3(t1, TAD_in, TAD_s)

                if D is not None and len(D) != 0:
                    m = []
                    for j in range(0, len(D)):
                        m.append(D[j][0].shape[0])
                    for j in range(0, len(D)):
                        diff1.append([t1, count])
                        diff2.append(D[j])
                        diff3.append([D[j][0][:, 1:3], count])
                        count = count + 1
                TAD_in = np.empty(shape=[0, 3])
        except:
            try:
                TAD1 = TAD1.reshape([1, 3])
                for i, t1 in enumerate(TAD1):
                    for t2 in TAD2:
                        if t1[1] - t2[1] < 20 and t2[2] - t1[2] < 20:
                            TAD_in = np.append(TAD_in, [t2], axis=0)
                    D = TAD_divid_V3(t1, TAD_in, TAD_s)
                    # print(t1, D)
                    if D is not None and len(D) != 0:
                        # dif.append([t1, D])
                        m = []
                        for j in range(0, len(D)):
                            m.append(D[j][0].shape[0])
                        # print(t1, D, max(m))

                        for j in range(0, len(D)):
                            if D[j][0].shape[0] == max(m):
                                diff1.append([t1, count])
                                diff2.append(D[j])
                                diff3.append([D[j][0][:, 1:3], count])
                                count = count + 1
                                max_region = D[j]
                                # print(t1, max_region[0][:, 1:3])
                    TAD_in = np.empty(shape=[0, 3])
            except:
                TAD2 = TAD2.reshape([1, 3])
                for i, t1 in enumerate(TAD1):
                    for t2 in TAD2:
                        if t1[1] - t2[1] < 20 and t2[2] - t1[2] < 20:
                            TAD_in = np.append(TAD_in, [t2], axis=0)
                    D = TAD_divid_V3(t1, TAD_in, TAD_s)
                    # print(t1, D)
                    if D is not None and len(D) != 0:
                        # dif.append([t1, D])
                        m = []
                        for j in range(0, len(D)):
                            m.append(D[j][0].shape[0])
                        # print(t1, D, max(m))

                        for j in range(0, len(D)):
                            if D[j][0].shape[0] == max(m):
                                diff1.append([t1, count])
                                diff2.append(D[j])
                                diff3.append([D[j][0][:, 1:3], count])
                                count = count + 1
                                max_region = D[j]
                                # print(t1, max_region[0][:, 1:3])
                    TAD_in = np.empty(shape=[0, 3])
        return np.array(diff1), np.array(diff2), np.array(diff3)

    def calculate_region_mean(D3, map1, map2, p1, p2):
        '''
        :param D3: divided region
        :param map1: contact map of cell1
        :param map2: contact map of cell2
        :return: average interaction between different divided region
        '''
        avr_all1 = []
        avr_all2 = []
        avr_diff = []
        for count, d0 in enumerate(D3):
            d = d0[0]
            avr_all1_temp = np.empty(shape=[d.shape[0], d.shape[0]])
            avr_all2_temp = np.empty(shape=[d.shape[0], d.shape[0]])
            for i in range(d.shape[0]):
                for j in range(d.shape[0]):
                    m1 = map1[int(d[i][0]): int(d[i][1]), int(d[j][0]): int(d[j][1])]
                    m2 = map2[int(d[i][0]): int(d[i][1]), int(d[j][0]): int(d[j][1])]
                    # print(p1)
                    avr1 = np.mean(m1[m1 < p1])  # need to change!!!!
                    avr_all1_temp[i, j] = avr1
                    avr2 = np.mean(m2[m2 < p2])  # need to change!!!!
                    avr_all2_temp[i, j] = avr2
            avr_all1.append([avr_all1_temp, count])
            avr_all2.append([avr_all2_temp, count])
            avr_diff.append([avr_all1_temp - avr_all2_temp, count])
        return avr_all1, avr_all2, avr_diff

    def count_cross_prob_eachregion(intract):
        R = np.zeros((intract.shape[0] - 1, intract.shape[0] - 1))
        for i in range(intract.shape[0] - 1):
            for j in range(intract.shape[0] - 1):
                R[i, j] = (intract[i, j + 1] * 2) / (intract[i, j] + intract[i + 1, j + 1])
        return R

    def count_cross_prob(avr, fold):
        '''
        count the ratio of interaction in conner region/ interaction in diagonal region
        :param avr: mean of interaction in each region
        :return: matrix of ratio, flag == 1 is real merge region
        '''
        ratio = []
        for count, t in enumerate(avr):
            intract = t[0]
            R = count_cross_prob_eachregion(intract)
            x1 = R.diagonal()[R.diagonal() > 0.47].shape[0]  # high ratio region define by o.45 as cutoff
            x2 = R.diagonal()[R.diagonal() < 0.47].shape[0]  # low ratio region
            flag = 0
            if R.shape[0] == 1:
                if R[0, 0] > 0.45:
                    flag = 1
            else:
                if x2 < x1:
                    flag = 1
            ratio.append([R / fold, count, flag])
        return ratio

    def count_intrc(idx, avr):
        '''
        count average interaction
        :param idx: index of real region
        :param avr: interact of all region
        :return: average interaction
        '''
        inter = np.empty(shape=(0, 3))
        for t in avr:
            if t[1] in idx:
                T = t[0]  # interaction array
                for i in range(T.shape[0] - 1):
                    inter_tmp = [T[i, i], T[i, i + 1], T[i + 1, i + 1]]
                    inter = np.append(inter, [inter_tmp], axis=0)
        return inter

    def loc_divid(D3, idx):
        '''
        find the location of divide region
        :param D3: all divided region
        :param idx: real region
        :return: location of real region
        '''
        reg = np.empty(shape=(0, 2))
        for d in D3:
            if d[1] in idx:
                reg = np.append(reg, d[0], axis=0)
        return reg

    def loc_union(D1, idx):
        '''
        find the location of union region
        :param D3: all divided region
        :param idx: real region
        :return: location of real region
        '''
        reg = np.empty(shape=(0, 3))
        for d in D1:
            if d[1] in idx:
                reg = np.append(reg, [d[0]], axis=0)
        return reg[:, 1:3]

    def remove_diff(loc_u, loc_d, TAD_s):
        '''
        remove differential region from similar region
        :param loc_u: TAD of union region
        :param loc_d: TAD of divide region
        :param TAD_s: TAD of similar region
        :return: similar region without differential region
        '''
        diff = np.append(loc_u, loc_d, axis=0)
        dlt = np.empty(shape=(0, 3))
        idx = []
        for i, t in enumerate(TAD_s):
            if min(abs(t[1] - diff[:, 1])) < 10 and min(abs(t[2] - diff[:, 2])) < 10:
                dlt = np.append(dlt, [t], axis=0)
                idx.append(i)
        TAD_s = np.delete(TAD_s, idx, 0)
        return TAD_s

    def compar_prob(ratio1, ratio2):
        '''
        find real split site
        :param ratio1: split region
        :param ratio2: merge region
        :return: real split site
        '''
        idx = []
        for i in range(len(ratio1)):
            # if ratio2[i][2] == 1 and ratio1[i][2] != 1:       # not need to unit
            if ratio1[i][2] != 1:
                ratio_tmp = ratio2[i][0] - ratio1[i][0]
                if np.mean(ratio_tmp.diagonal()) > 0.09:  # real divide region by cutoff: -0.1
                    idx.append(ratio2[i][1])
        return idx

    def divid_region(f1, f2, D3, D1, TAD_s, p1, p2, fold):
        '''
        main function for divided region detection
        :param f1: cell line1 (divided)
        :param f2: cell line2
        :param D3: divided region in cell line 1
        :param D1: Unite region in cell line 2
        :param TAD_s: similar region
        :return: divid1: mean interaction in divide region, divid2: mean interaction in union region, loc_d: location of
        divided region, loc_u: location of union region, dlt: similar region after removing differential region
        '''
        map1 = np.loadtxt(f1)
        map2 = np.loadtxt(f2)
        avr1, avr2, avr_diff = calculate_region_mean(D3, map1, map2, p1, p2)  # count mean interaction in different region
        ratio1 = count_cross_prob(avr1, 1)  # count diff ratio
        ratio2 = count_cross_prob(avr2, fold)
        idx = compar_prob(ratio1, ratio2)  # find real divide region index
        if len(idx) > 0:
            divid1 = count_intrc(idx, avr1)  # count interaction in real divide region
            divid2 = count_intrc(idx, avr2)
            loc_u = loc_union(D1, idx)
            loc_u = np.insert(loc_u, 0, np.array(range(loc_u.shape[0])).transpose(), axis=1)
            loc_d = loc_divid(D3, idx)
            loc_d = np.insert(loc_d, 0, np.array(range(loc_d.shape[0])).transpose(), axis=1)
            dlt = remove_diff(loc_u, loc_d, TAD_s)
            return divid1, divid2, loc_d, loc_u, dlt
        else:
            return 0, 0, 0, 0, 0

    def find_t1(TAD):
        def is_t1(t, TAD):
            '''
            t is T1 TAD return 1, else return 0
            '''
            flag = 1
            for t0 in TAD:
                if not np.array_equal(t0, t):
                    # print([t0, t[1] - t0[1], t0[2] - t[2]])
                    if t[1] - t0[1] <= 0 and t0[2] - t[2] <= 0:
                        flag = 0
            return flag

        t1 = np.empty(shape=[0, 3])
        for t in TAD:
            if is_t1(t, TAD) == 1:
                t1 = np.append(t1, [t], axis=0)
        return t1

    def get_ratio(map, TAD, up):
        def neighbor_TAD(TAD):
            '''
            find each pair of nearby TAD
            '''
            pair = []
            for i in range(TAD.shape[0]):
                for j in range(TAD.shape[0]):
                    if abs(TAD[i, 2] - TAD[j, 1]) < 20:
                        pair.append([TAD[i], TAD[j]])
            return np.array(pair)

        def calculate_ratio(map, pair):
            ratio = []
            for T in pair:
                t1 = T[0]
                t2 = T[1]
                m1 = map[int(t1[1]):int(t1[2]), int(t1[1]):int(t1[2])]
                m2 = map[int(t2[1]):int(t2[2]), int(t2[1]):int(t2[2])]
                m3 = map[int(t1[1]):int(t1[2]), int(t2[1]):int(t2[2])]
                mean_m1 = np.mean(m1[m1 < up])
                mean_m2 = np.mean(m2[m2 < up])
                mean_m3 = np.mean(m3[m3 < up])
                r = mean_m3 / (mean_m1 + mean_m2)
                ratio.append(r)
            return ratio

        pair = neighbor_TAD(TAD)
        ratio = calculate_ratio(map, pair)
        return ratio

    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython3 DACTAD.py TAD_calculator <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python DACS.py TAD_calculator -h\n")
        return 0

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython3 DACS.py corner_split <contact_map_file_paths> "
                                           "[optional arguments]\n\n", description='')
    parser.add_argument('command', default=None,
                        help="set as 'corner_split' to identify split TAD")

    parser.add_argument('-c', '--contact_maps', dest="contact_map", default=None,
                        help="paths to two compared Hi-C contact maps. paths must be separated by the comma ','.")

    parser.add_argument('--contact_maps_aliases', dest="aliases", default=None,
                        help="A set of short aliases for the contact map. Paths must be separated by the comma ','.")

    parser.add_argument('-t', '--TAD', dest="TAD", default=None,
                        help="input files of TADs for two compared Hi-C contact maps. Paths must be separated by the"
                             " comma ','.")

    parser.add_argument('-u', '--up_limits', dest="up", default=None,
                        help="up limits for two compared Hi-C contact maps, paths must be separated by the comma ','.")

    parser.add_argument('-j', '--adjust_quality', dest="adjust_quality", default=0, type=int,
                        help="set as 1 to normalize sequence quality for two Hi-C contact maps, set as 0 not to "
                             "normalize sequence quality for two Hi-C contact maps")

    parser.add_argument('-o', '--output', dest="output", default=None,
                        help="path to output files")

    parser.add_argument('-d', '--split_direction', dest="direction", default=0, type=int,
                        help="set as 0: output TADs spliced in both two contact maps, set as 1: output TADs spliced in "
                             "contact map 1, set as 2: output TADs spliced in contact map 2")

    args = parser.parse_args()

    file = args.contact_map.split(',')
    TAD = args.TAD.split(',')
    up = args.up.split(',')
    aliases = args.aliases.split(',')

    TAD1 = np.loadtxt(TAD[0])
    TAD2 = np.loadtxt(TAD[1])

    TAD_s,TAD1_only,TAD2_only=region_detect(TAD2, TAD1)

    D1, D2, D3 = region_divid_v3(TAD1, TAD2, TAD_s)

    if args.adjust_quality == 0:
        divid1_1, divid2_1, loc_d, loc_u, dlt = divid_region(file[1], file[0], D3, D1, TAD_s, float(up[1]), float(up[0]),
                                                             1)      # contact_map1 -> contact_map2
    elif args.adjust_quality == 1:
        map1 = np.loadtxt(file[0])
        map2 = np.loadtxt(file[1])
        ratio1 = get_ratio(map1, TAD1, float(up[0]))
        ratio2 = get_ratio(map2, TAD2, float(up[1]))
        fold = np.mean(ratio1)/np.mean(ratio2)
        print(fold)
        split1_1, split2_1, loc_d, loc_u, dlt = divid_region(file[1], file[0], D3, D1, TAD_s, float(up[1]),
                                                             float(up[0]), fold)  # contact_map1 -> contact_map2
    else:
        print("\nusage:\npython3 DACS.py TAD_calculator <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python DACS.py TAD_calculator -h\n")

    try:
        if split1_1 == 0:
            print('No spliced region')
    except:
        if (not os.path.exists(args.output)):
            print('path does not exist\n')
        else:
            np.savetxt(os.path.join(args.output, aliases[0]+'->'+aliases[1]+'.split.txt'), loc_d)
            np.savetxt(os.path.join(args.output, aliases[0]+'->'+aliases[1]+'.merge.txt'), loc_u)


def TAD_similarity(command='TAD_similarity'):
    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython3 DACTAD.py TAD_similarity <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python DACTAD.py TAD_similarity -h\n")
        return 0

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython3 DACTAD.py TAD_similarity <contact_map_file_paths> "
                                           "[optional arguments]\n\n", description='')
    parser.add_argument('command', default=None,
                        help="set as 'TAD_similarity' to calculate similarity of two TADs")

    parser.add_argument('-c', '--contact_maps', dest="contact_map", default=None,
                        help="paths to two compared Hi-C contact maps. paths must be separated by the comma ','.")

    parser.add_argument('-t', '--TAD', dest="TAD", default=None,
                        help="input files of TADs for two compared Hi-C contact maps. Paths must be separated by the"
                             " comma ','.")

    parser.add_argument('-o', '--output', dest="output", default=None,
                        help="path to output files")

    args = parser.parse_args()

    def similarity_scc(map1, map2, TAD):
        r = robjects.r

        r('''
            library(hicrep)
            hicrep_similarity<-function(map1, map2, TAD){
              scc = c()
              for(i in 1:nrow(TAD)){
                start = TAD[i, 2]
                end = TAD[i, 3]
                hic1 = data.frame('chr', 1:(end - start+1), 2:(end - start+2), map1[start:end,start:end])
                hic2 = data.frame('chr', 1:(end - start+1), 2:(end - start+2), map2[start:end,start:end])
                processed <- prep(hic1, hic2, 1, 0, 5)
                scc.out = get.scc(processed, 1, 5)
                scc=c(scc,scc.out$scc)
              }
              TAD2 = data.frame(TAD,scc)
              return(TAD2)
            }
        ''')

        rpy2.robjects.numpy2ri.activate()
        pandas2ri.activate()

        nr, nc = TAD.shape
        TAD_r = r.matrix(TAD, nrow=nr, ncol=nc)
        r.assign("TAD1", TAD_r)

        hicrep_similarity = r('hicrep_similarity')
        scc = hicrep_similarity(map1, map2, TAD_r)
        scc = pandas2ri.ri2py(scc)
        return scc

    def similarity_Laplacian(map1, map2, TAD):
        def get_Laplacian(M):
            S = M.sum(1)
            i_nz = np.where(S > 0)[0]
            S = S[i_nz]
            M = (M[i_nz].T)[i_nz].T
            S = 1 / np.sqrt(S)
            M = S * M
            M = (S * M.T).T
            n = np.size(S)
            M = np.identity(n) - M
            M = (M + M.T) / 2
            return M

        def evec_distance(v1, v2):
            d1 = np.dot(v1 - v2, v1 - v2)
            d2 = np.dot(v1 + v2, v1 + v2)
            if d1 < d2:
                d = d1
            else:
                d = d2
            return np.sqrt(d)

        def get_ipr(evec):
            ipr = 1.0 / (evec * evec * evec * evec).sum()
            return ipr

        def get_reproducibility(M1, M2, num_evec=20):
            k1 = np.sign(M1).sum(1)
            d1 = np.diag(M1)
            kd1 = ~((k1 == 1) * (d1 > 0))
            k2 = np.sign(M2).sum(1)
            d2 = np.diag(M2)
            kd2 = ~((k2 == 1) * (d2 > 0))
            iz = np.nonzero((k1 + k2 > 0) * (kd1 > 0) * (kd2 > 0))[0]
            M1b = (M1[iz].T)[iz].T
            M2b = (M2[iz].T)[iz].T

            i_nz1 = np.where(M1b.sum(1) > 0)[0]
            i_nz2 = np.where(M2b.sum(1) > 0)[0]
            M1b_L = get_Laplacian(M1b)
            M2b_L = get_Laplacian(M2b)

            a1, b1 = eigsh(M1b_L, k=num_evec, which="SM")
            a2, b2 = eigsh(M2b_L, k=num_evec, which="SM")

            b1_extend = np.zeros((np.size(M1b, 0), num_evec))
            b2_extend = np.zeros((np.size(M2b, 0), num_evec))
            for i in range(num_evec):
                b1_extend[i_nz1, i] = b1[:, i]
                b2_extend[i_nz2, i] = b2[:, i]

            ipr_cut = 5
            ipr1 = np.zeros(num_evec)
            ipr2 = np.zeros(num_evec)
            for i in range(num_evec):
                ipr1[i] = get_ipr(b1_extend[:, i])
                ipr2[i] = get_ipr(b2_extend[:, i])

            b1_extend_eff = b1_extend[:, ipr1 > ipr_cut]
            b2_extend_eff = b2_extend[:, ipr2 > ipr_cut]
            num_evec_eff = min(np.size(b1_extend_eff, 1), np.size(b2_extend_eff, 1))

            evd = np.zeros(num_evec_eff)
            for i in range(num_evec_eff):
                evd[i] = evec_distance(b1_extend_eff[:, i], b2_extend_eff[:, i])

            Sd = evd.sum()
            l = np.sqrt(2)
            evs = abs(l - Sd / num_evec_eff) / l
            return evs

        def similarity(map1, map2, TAD_f):
            e = []

            for i in TAD_f:
                n1 = int(i[1])
                n2 = int(i[2])
                x1 = map1[n1:(n2 + 1), n1:(n2 + 1)]
                x2 = map2[n1:(n2 + 1), n1:(n2 + 1)]
                e.append(get_reproducibility(x1, x2))

            e = np.array(e)
            e = e.reshape(e.shape[0], 1)
            TAD_f = np.concatenate((TAD_f, e), 1)
            return TAD_f
        Laplacian = similarity(map1, map2, TAD)
        return Laplacian

    def hash_similarity(map1, map2, TAD):
        e_ahash = []
        for i in TAD:
            n1 = int(i[0])
            n2 = int(i[1])
            x1 = map1[n1:(n2 + 1), n1:(n2 + 1)]
            x2 = map2[n1:(n2 + 1), n1:(n2 + 1)]
            x1 = x1.astype(np.uint8)
            x2 = x2.astype(np.uint8)
            x2 = x2 / np.mean(x2)
            img1 = Image.fromarray(x1)
            img2 = Image.fromarray(x2)
            hash1 = imagehash.average_hash(img1, 16)
            hash2 = imagehash.average_hash(img2, 16)
            d = 1 - (hash2 - hash1) / 256
            e_ahash.append(d)
        e_ahash = np.array(e_ahash)
        e_ahash = e_ahash.reshape(e_ahash.shape[0], 1)
        TAD = np.concatenate((TAD, e_ahash), 1)
        return TAD

    file1 = args.contact_map.split(',')
    file2 = args.TAD.split(',')
    map1 = np.loadtxt(file1[0])
    map2 = np.loadtxt(file1[1])
    TAD1 = np.loadtxt(file2[0])
    TAD2 = np.loadtxt(file2[1])
    TAD = np.concatenate((TAD1, TAD2), axis=0)
    scc = similarity_scc(map1, map2, TAD)
    Laplacian = similarity_Laplacian(map1, map2, TAD)
    hash = hash_similarity(map1, map2, TAD)
    np.savetxt(os.path.join(args.output, 'similarity_scc.txt'), scc, delimiter='\t')
    np.savetxt(os.path.join(args.output, 'similarity_laplacian.txt'), Laplacian, delimiter='\t')
    np.savetxt(os.path.join(args.output, 'similarity_hash.txt'), hash, delimiter='\t')


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == 'TAD_calculator':
            TAD_calculator(command='TAD_calculator')
        elif sys.argv[1] == 'corner_split':
            corner_split(command='corner_split')
        elif sys.argv[1] == 'TAD_similarity':
            TAD_similarity(command='TAD_similarity')
        else:
            printHelp()
    else:
        print('\nDACTAD version 1.0.0')
        print('For a list of functions in DACTAD, please try:\npython DACTAD.py -h')
        print('')








