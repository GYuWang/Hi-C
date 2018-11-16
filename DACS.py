import rpy2.robjects as robjects
import numpy as np
from itertools import combinations
import sys
import argparse
import cv2
import math
import os
import warnings
from rpy2.rinterface import RRuntimeWarning



def printHelp():
    print('\nDACS version 1.0.0')
    print('For help information for each function, try:\npython3 DACS.py <function> -h')
    print('\nFunctions:')
    print('\tTAD_calculator:\n\t\tidentify the topological domain\n')
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
                                     usage="\n\npython3 DACS.py TAD_calculator <contact_map_file_paths> "
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
                        help="path to the output file of TADs for Hi-C contact map")

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


def sim_range(TAD):
    r = []
    for t in TAD:
        r1 = (t[8] - t[7] + t[5] - t[4]) / (t[2] - t[1])
        r.append(r1)

    r = np.array(r)
    r = r.reshape(r.shape[0], 1)
    TAD = np.concatenate((TAD, r), 1)
    return TAD


def RangInList(L, TAD):
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
        # f = RangInList(range1, TAD_s.tolist())  # remove similar large region
        # print(range1)
        # if f == 1:
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
                # print(x, r)
                # if r > 0.9 and r < 1.1 and sum(flag3) != 0:
                if r > 0.8 and r < 1.1:  # change parameter
                    flag2 = 0  # find all real sep
                    for l in range(x.shape[0] - 1):
                        # print(x[l,:], x.shape[0])
                        if x[l + 1, 1] - x[l, 2] > 20 or x[l, 2] - x[l + 1, 1] > 20:  # change parameter
                            flag2 = 1
                    if flag2 == 0:
                        dic.append([x, r])
                        # print(x, r)
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
            # print('-------')
            # print(t1, TAD_in)
            # print('+++++++')
            D = TAD_divid_V3(t1, TAD_in, TAD_s)

            if D is not None and len(D) != 0:
                # dif.append([t1, D])
                m = []
                for j in range(0, len(D)):
                    m.append(D[j][0].shape[0])

                # print(t1, D, max(m))

                for j in range(0, len(D)):
                    # if D[j][0].shape[0] == max(m):             # choose the most split (not need)
                    diff1.append([t1, count])
                    diff2.append(D[j])
                    diff3.append([D[j][0][:, 1:3], count])
                    count = count + 1
                    # max_region = D[j]
                    # print(t1, max_region[0][:, 1:3])
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


###### find divide regions

def linear_transform(map):
    '''
    use linear regression to transform matrix
    :param map: contact map
    :return: transformed contact map
    '''
    x = np.empty(shape=[1])
    y = np.empty(shape=[1])


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
            # print(i, j, R[i, j], intract[i, j + 1], intract[i, j], intract[i + 1, j + 1])
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


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == 'TAD_calculator':
            TAD_calculator(command='TAD_calculator')

        else:
            printHelp()
    else:
        print('\nDACS version 1.0.0')
        print('For a list of functions in DACS, please try:\npython DACS.py -h')
        print('')








