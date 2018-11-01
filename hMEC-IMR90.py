import numpy as np
from sklearn import linear_model
from itertools import combinations
import os
import pandas as pd


def region_detect(TAD1, TAD2):
    flag1 = []
    flag2 = []
    band_s = np.array([0,0,0])
    try:
        for i in range(0,TAD1.shape[0]):
            for j in range(0, TAD2.shape[0]):
                if abs(TAD1[i,2] - TAD2[j,2]) <= 10 and abs(TAD1[i,1] - TAD2[j,1]) <= 10:
                    band_s = np.vstack([band_s, (TAD1[i,:]+TAD2[j,:])/2])
                    flag1.append(i)
                    flag2.append(j)
        TAD1_only = np.delete(TAD1, flag1, axis=0)
        TAD2_only = np.delete(TAD2, flag2, axis=0)
    except:
        try:
            TAD2 = TAD2.reshape([1,3])
            for i in range(0,TAD1.shape[0]):
                for j in range(0, TAD2.shape[0]):
                    if abs(TAD1[i,2] - TAD2[j,2]) <= 10 and abs(TAD1[i,1] - TAD2[j,1]) <= 10:
                        band_s = np.vstack([band_s, (TAD1[i,:]+TAD2[j,:])/2])
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
        return band_s[1:,:], TAD1_only, TAD2_only
    except:
        return band_s, TAD1_only, TAD2_only


def sim_range(TAD):
    r = []
    for t in TAD:
        r1 = (t[8] - t[7] + t[5] - t[4])/(t[2] - t[1])
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
        #f = RangInList(range1, TAD_s.tolist())  # remove similar large region
        #print(range1)
        #if f == 1:
        dic = []
        for i in range(2, range2.shape[0] + 1):
            comb = list(combinations(range(range2.shape[0]), i))
            for j, t in enumerate(comb):
                x = range2[t, :]
                length = 0
                flag3 = []                                                  # all divide region are similar region
                for s in x:
                    length = length + s[2] - s[1]
                    f2 = RangInList(s, TAD_s.tolist())
                    flag3.append(f2)
                r = length / (range1[2] - range1[1])
                #print(x, r)
                #if r > 0.9 and r < 1.1 and sum(flag3) != 0:
                if r > 0.8 and r < 1.1:                                 # change parameter
                    flag2 = 0  # find all real sep
                    for l in range(x.shape[0] - 1):
                        # print(x[l,:], x.shape[0])
                        if x[l + 1, 1] - x[l, 2] > 20 or x[l, 2] - x[l + 1, 1] > 20:     # change parameter
                            flag2 = 1
                    if flag2 == 0:
                        dic.append([x, r])
                        # print(x, r)
        return(dic)




def region_divid_v3(TAD1, TAD2, TAD_s):
    '''
    :param TAD1: whole region
    :param TAD2: divide region
    :param TAD_s: similar region
    :return: D1: TAD of unit region, D2: divide TAD with coverage rate, D3: divide TAD without coverage rate
    '''
    TAD_in = np.empty(shape = [0, 3])
    diff1 = []
    diff2 = []
    diff3 = []
    count = 0
    try:
        for i, t1 in enumerate(TAD1):
            for t2 in TAD2:
                if t1[1] - t2[1] < 20 and t2[2] - t1[2] < 20:
                    TAD_in = np.append(TAD_in, [t2], axis=0)
            #print('-------')
            #print(t1, TAD_in)
            #print('+++++++')
            D = TAD_divid_V3(t1, TAD_in, TAD_s)

            if D is not None and len(D) != 0:
                #dif.append([t1, D])
                m = []
                for j in range(0, len(D)):
                    m.append(D[j][0].shape[0])

                #print(t1, D, max(m))

                for j in range(0, len(D)):
                    #if D[j][0].shape[0] == max(m):             # choose the most split (not need)
                        diff1.append([t1, count])
                        diff2.append(D[j])
                        diff3.append([D[j][0][:, 1:3], count])
                        count = count + 1
                        #max_region = D[j]
                #print(t1, max_region[0][:, 1:3])
            TAD_in = np.empty(shape=[0, 3])
    except:
        try:
            TAD1 = TAD1.reshape([1,3])
            for i, t1 in enumerate(TAD1):
                for t2 in TAD2:
                    if t1[1] - t2[1] < 20 and t2[2] - t1[2] < 20:
                        TAD_in = np.append(TAD_in, [t2], axis=0)
                D = TAD_divid_V3(t1, TAD_in, TAD_s)
                #print(t1, D)
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
                #print(t1, D)
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
                #print(p1)
                avr1 = np.mean(m1[m1 < p1])             # need to change!!!!
                avr_all1_temp[i, j] = avr1
                avr2 = np.mean(m2[m2 < p2])             # need to change!!!!
                avr_all2_temp[i, j] = avr2
        avr_all1.append([avr_all1_temp, count])
        avr_all2.append([avr_all2_temp, count])
        avr_diff.append([avr_all1_temp - avr_all2_temp, count])
    return avr_all1, avr_all2, avr_diff




def calculate_region_split100(D3, map, n):
    '''
    split each region into n small region
    :param D3: divided region
    :param map1: contact map of cell1
    :param n: region number
    :return: average interaction between different divided region
    '''
    avr_all = np.empty(shape=(0, 3, 30, 30))
    for count, d0 in enumerate(D3):
        d = d0[0]
        for i in range(d.shape[0]-1):
            s1 = np.array_split(range(int(d[i][0]), int(d[i][1])), n)
            s2 = np.array_split(range(int(d[i+1][0]), int(d[i+1][1])), n)
            avr_all2_temp_1 = np.empty(shape=[n,n])
            avr_all2_temp_2 = np.empty(shape=[n, n])
            avr_all2_temp_3 = np.empty(shape=[n, n])
            for t1 in range(n):
                for t2 in range(n):
                    m = map[s1[t1][0]:(s1[t1][-1]+1), s1[t2][0]:(s1[t2][-1]+1)]
                    m[m > 10] = 10
                    avr1 = np.mean(m)
                    avr_all2_temp_1[t1, t2] = avr1
            for t1 in range(n):
                for t2 in range(n):
                    m = map[s1[t1][0]:(s1[t1][-1]+1), s2[t2][0]:(s2[t2][-1]+1)]
                    m[m > 10] = 10
                    avr1 = np.mean(m)
                    avr_all2_temp_2[t1, t2] = avr1
            for t1 in range(n):
                for t2 in range(n):
                    m = map[s2[t1][0]:(s2[t1][-1]+1), s2[t2][0]:(s2[t2][-1]+1)]
                    m[m > 10] = 10
                    avr1 = np.mean(m)
                    avr_all2_temp_3[t1, t2] = avr1
            #avr_all0_temp.append(np.array([avr_all2_temp_1,avr_all2_temp_2,avr_all2_temp_3]))
            #print(np.array([[avr_all2_temp_1, avr_all2_temp_2, avr_all2_temp_3]]).shape)
            #print(avr_all.shape)
            avr_all = np.append(avr_all, [[avr_all2_temp_1, avr_all2_temp_2, avr_all2_temp_3]], axis=0)
        #avr_all.append(np.array(avr_all0_temp))
    return avr_all






def count_cross_prob_eachregion(intract):
    R = np.zeros((intract.shape[0] - 1, intract.shape[0] - 1))
    for i in range(intract.shape[0] - 1):
        for j in range(intract.shape[0] - 1):
            R[i, j] = (intract[i, j + 1] * 2) / (intract[i, j] + intract[i + 1, j + 1])
            #print(i, j, R[i, j], intract[i, j + 1], intract[i, j], intract[i + 1, j + 1])
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
        x1 = R.diagonal()[R.diagonal() > 0.47].shape[0]   # high ratio region define by o.45 as cutoff
        x2 = R.diagonal()[R.diagonal() < 0.47].shape[0]   # low ratio region
        flag = 0
        if R.shape[0] == 1:
            if R[0, 0] > 0.45:
                flag = 1
        else:
            if x2 < x1:
                flag = 1
        ratio.append([R/fold, count, flag])
    return ratio



def count_intrc(idx, avr):
    '''
    count average interaction
    :param idx: index of real region
    :param avr: interact of all region
    :return: average interaction
    '''
    inter = np.empty(shape=(0,3))
    for t in avr:
        if t[1] in idx:
            T = t[0]                    # interaction array
            for i in range(T.shape[0] - 1):
                inter_tmp = [T[i, i], T[i, i+1], T[i+1, i+1]]
                inter = np.append(inter, [inter_tmp], axis=0)
    return inter



def loc_divid(D3, idx):
    '''
    find the location of divide region
    :param D3: all divided region
    :param idx: real region
    :return: location of real region
    '''
    reg = np.empty(shape=(0,2))
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
    reg = np.empty(shape=(0,3))
    for d in D1:
        if d[1] in idx:
            reg = np.append(reg, [d[0]], axis=0)
    return reg[:,1:3]





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
        #if ratio2[i][2] == 1 and ratio1[i][2] != 1:       # not need to unit
        if ratio1[i][2] != 1:
            ratio_tmp = ratio2[i][0] - ratio1[i][0]
            if np.mean(ratio_tmp.diagonal()) > 0.09:    # real divide region by cutoff: -0.1
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
    avr1, avr2, avr_diff = calculate_region_mean(D3, map1, map2, p1, p2)      # count mean interaction in different region
    ratio1 = count_cross_prob(avr1, 1)                                   # count diff ratio
    ratio2 = count_cross_prob(avr2, fold)
    idx = compar_prob(ratio1, ratio2)                                 # find real divide region index
    if len(idx) > 0:
        divid1 = count_intrc(idx, avr1)                                 # count interaction in real divide region
        divid2 = count_intrc(idx, avr2)
        loc_u = loc_union(D1, idx)
        loc_u = np.insert(loc_u,0,  np.array(range(loc_u.shape[0])).transpose(), axis=1)
        loc_d = loc_divid(D3, idx)
        loc_d = np.insert(loc_d, 0, np.array(range(loc_d.shape[0])).transpose(), axis=1)
        dlt = remove_diff(loc_u, loc_d, TAD_s)
        return divid1, divid2, loc_d, loc_u, dlt
    else:
        return 0, 0, 0, 0, 0



# whole genome analysis

list_IRM90 = os.listdir('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file')
list_IMR90_band = []
list_IMR90_map = []
for l in list_IRM90:
    if('band' in l):
        list_IMR90_band.append(l)
    else:
        list_IMR90_map.append(l)

list_HUVEC = os.listdir('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HUVEC/new_file')
list_HUVEC_band = []
list_HUVEC_map = []
for l in list_HUVEC:
    if('band' in l):
        list_HUVEC_band.append(l)
    else:
        list_HUVEC_map.append(l)

chr = []
up = []
down = []
for l in list_IMR90_map:
    s = l.replace('.', '_').split('_')
    if(len(s)>5):
        chr.append(s[0])
        up.append(s[6])
        down.append(s[5])


lim = pd.read_csv('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/color_lim2_withHESC.txt', sep='\t')



'''
heatmap hMSC -> IMR90
'''

print(chr)

for i in range(1, len(chr)):
    #print(i)
    try:

        #if chr[i] == '13':

            TAD1 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file/'+chr[i]+'_matrix_IMR90_Coverage.txt.'+down[i] + '.'+up[i] +'.new.band.txt')
            # TAD1 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/band_hESC/' + chr[
            #     i] + '_matrix_IMR90_Coverage.txt.' + down[i] + '.' + up[i] + '.new.band.txt')

            TAD2 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/hESC/Contact_map/map/band/'+chr[i]+'_hESC_Coverage.txt.' + down[i] + '.'+up[i] +'.new.band.txt')

            TAD_s, TAD1_only, TAD2_only = region_detect(TAD2, TAD1)
            D1, D2, D3 = region_divid_v3(TAD1, TAD2, TAD_s)

            # print('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file2/'+chr[i]+'_matrix_IMR90_Coverage.txt.'+down[i] + '.'+up[i] +'.new2.band.txt')
            # print('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HUVEC/new_file2/'+chr[i]+'_matrix_HUVEC_Coverage.txt.' + down[i] + '.'+up[i] +'.new2.band.txt')
            f1 = '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file/'+chr[i]+'_matrix_IMR90_Coverage.txt.'+down[i] + '.' +up[i] + '.new'
            f2 = '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/hESC/Contact_map/map/band/'+chr[i]+'_hESC_Coverage.txt.' + down[i] + '.'+up[i] +'.new'
            lim1 = lim.iloc[int(chr[i])-1, 10]
            lim2 = lim.iloc[int(chr[i])-1, 8]
            #print(D3)
            #print(i)

            divid1_1, divid2_1, loc_d, loc_u, dlt = divid_region(f2, f1, D3, D1, TAD_s, lim2, lim1, 0.5)
            #print('YYYY')


            try:
                if divid1_1 == 0:
                    print('No')
            except:
                map1 = np.loadtxt(f1)
                map2 = np.loadtxt(f2)

                test1 = calculate_region_split100(D3, map2, 30)
                mean1_1 = sum(test1[:, 0]) / test1.shape[0]
                mean1_2 = sum(test1[:, 1]) / test1.shape[0]
                mean1_3 = sum(test1[:, 2]) / test1.shape[0]

                a1 = np.append(mean1_1, mean1_2, axis=1)
                a2 = np.append(mean1_2, mean1_3, axis=1)

                np.savetxt(
                    '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/hMSC-IMR90/' + chr[i] + '.' + down[
                        i] + '.' + up[i] + '.divid1_1', divid1_1)
                np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/hMSC-IMR90/'+chr[i]+'.' + down[i] + '.'+up[i] +'.divid2_1', divid2_1)
                np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/hMSC-IMR90/'+chr[i]+'.' + down[i] + '.'+up[i] +'.loc_d', loc_d)
                np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/hMSC-IMR90/'+chr[i]+'.' + down[i] + '.'+up[i] +'.loc_u', loc_u)
                np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/hMSC-IMR90/'+chr[i]+'.' + down[i] + '.'+up[i] +'.dlt', dlt)
                print(f1)
    except:
        #print(i)
        pass




# for i in range(1, len(chr)):
#     print(i)
#     try:
#
#         #if chr[i] == '13':
#             print(chr[i])
#             TAD2 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file/'+chr[i]+'_matrix_IMR90_Coverage.txt.'+down[i] + '.'+up[i] +'.new.band.txt')
#             # TAD2 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/band_hESC/' + chr[
#             #     i] + '_matrix_IMR90_Coverage.txt.' + down[i] + '.' + up[i] + '.new.band.txt')
#
#             TAD1 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/hESC/Contact_map/map/band/'+chr[i]+'_hESC_Coverage.txt.' + down[i] + '.'+up[i] +'.new.band.txt')
#             TAD_s, TAD1_only, TAD2_only = region_detect(TAD2, TAD1)
#             D1, D2, D3 = region_divid_v3(TAD1, TAD2, TAD_s)
#
#             # print('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file2/'+chr[i]+'_matrix_IMR90_Coverage.txt.'+down[i] + '.'+up[i] +'.new2.band.txt')
#             # print('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HUVEC/new_file2/'+chr[i]+'_matrix_HUVEC_Coverage.txt.' + down[i] + '.'+up[i] +'.new2.band.txt')
#             f2 = '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file/'+chr[i]+'_matrix_IMR90_Coverage.txt.'+down[i] + '.' +up[i] + '.new'
#             f1 = '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/hESC/Contact_map/map/band/'+chr[i]+'_hESC_Coverage.txt.' + down[i] + '.'+up[i] +'.new'
#             lim2 = lim.iloc[int(chr[i])-1, 10]
#             lim1 = lim.iloc[int(chr[i])-1, 8]
#             #print(D3)
#             #print(i)
#
#             divid1_1, divid2_1, loc_d, loc_u, dlt = divid_region(f2, f1, D3, D1, TAD_s, lim2, lim1, 2.2)
#             #print('YYYY')
#
#
#             try:
#                 if divid1_1 == 0:
#                     print('No')
#             except:
#                 map1 = np.loadtxt(f1)
#                 map2 = np.loadtxt(f2)
#
#                 test1 = calculate_region_split100(D3, map2, 30)
#                 mean1_1 = sum(test1[:, 0]) / test1.shape[0]
#                 mean1_2 = sum(test1[:, 1]) / test1.shape[0]
#                 mean1_3 = sum(test1[:, 2]) / test1.shape[0]
#
#                 a1 = np.append(mean1_1, mean1_2, axis=1)
#                 a2 = np.append(mean1_2, mean1_3, axis=1)
#                 #b = np.append(a1, a2, axis=0)
#
#                 np.savetxt(
#                     '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/IMR90-hMSC/' + chr[i] + '.' + down[
#                         i] + '.' + up[i] + '.divid1_1', divid1_1)
#                 np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/IMR90-hMSC/'+chr[i]+'.' + down[i] + '.'+up[i] +'.divid2_1', divid2_1)
#                 np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/IMR90-hMSC/'+chr[i]+'.' + down[i] + '.'+up[i] +'.loc_d', loc_d)
#                 np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/IMR90-hMSC/'+chr[i]+'.' + down[i] + '.'+up[i] +'.loc_u', loc_u)
#                 np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/IMR90-hMSC/'+chr[i]+'.' + down[i] + '.'+up[i] +'.dlt', dlt)
#                 print(f1)
#     except:
#         #print(i)
#         pass


# chr2 = '3'
# down2='17000'
# up2='18000'
#
#
# TAD2 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file/'+chr2+'_matrix_IMR90_Coverage.txt.'+down2 + '.'+up2 +'.new.band.txt')
# TAD1 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/hESC/Contact_map/map/band/'+chr2+'_hESC_Coverage.txt.' + down2 + '.'+up2+'.new.band.txt')
# TAD_s, TAD1_only, TAD2_only = region_detect(TAD2, TAD1)
# D1, D2, D3 = region_divid_v3(TAD1, TAD2, TAD_s)
#
#
#
# f2 = '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file/'+chr2+'_matrix_IMR90_Coverage.txt.'+down2 + '.' +up2 + '.new'
# f1 = '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/hESC/Contact_map/map/band/'+chr2+'_hESC_Coverage.txt.' + down2 + '.'+up2 +'.new'
#
#
# divid1_1, divid2_1, loc_d, loc_u, dlt = divid_region(f2, f1, D3, D1, TAD_s, 30, 5, 2.2)
#
# print(loc_u)
#
#
# for i in range(1, len(chr)):
#
#     TAD2 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file/'+chr[i]+'_matrix_IMR90_Coverage.txt.'+down[i] + '.'+up[i] +'.new.band.txt')
#     # TAD2 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/band_hESC/' + chr[
#     #     i] + '_matrix_IMR90_Coverage.txt.' + down[i] + '.' + up[i] + '.new.band.txt')
#
#     TAD1 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/hESC/Contact_map/map/band/'+chr[i]+'_hESC_Coverage.txt.' + down[i] + '.'+up[i] +'.new.band.txt')
#     TAD_s, TAD1_only, TAD2_only = region_detect(TAD2, TAD1)
#     D1, D2, D3 = region_divid_v3(TAD1, TAD2, TAD_s)
#
#     f4 = '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file/'+chr[i]+'_matrix_IMR90_Coverage.txt.'+down[i] + '.' +up[i] + '.new'
#     f3 = '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/hESC/Contact_map/map/band/'+chr[i]+'_hESC_Coverage.txt.' + down[i] + '.'+up[i] +'.new'
#     lim2 = 30
#     lim1 = 4
#     #print(D3)
#     #print(i)
#
#     divid1_12, divid2_12, loc_d2, loc_u, dlt2 = divid_region(f4, f3, D3, D1, TAD_s, lim2, lim1)
#     print('YYYY')
#
#
#             try:
#                 if divid1_1 == 0:
#                     print('No')
#             except:
#                 map1 = np.loadtxt(f1)
#                 map2 = np.loadtxt(f2)
#
#                 test1 = calculate_region_split100(D3, map2, 30)
#                 mean1_1 = sum(test1[:, 0]) / test1.shape[0]
#                 mean1_2 = sum(test1[:, 1]) / test1.shape[0]
#                 mean1_3 = sum(test1[:, 2]) / test1.shape[0]
#
#                 a1 = np.append(mean1_1, mean1_2, axis=1)
#                 a2 = np.append(mean1_2, mean1_3, axis=1)
#                 #b = np.append(a1, a2, axis=0)
#
#                 np.savetxt(
#                     '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/IMR90-hMSC/' + chr[i] + '.' + down[
#                         i] + '.' + up[i] + '.divid1_1', divid1_1)
#                 np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/IMR90-hMSC/'+chr[i]+'.' + down[i] + '.'+up[i] +'.divid2_1', divid2_1)
#                 np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/IMR90-hMSC/'+chr[i]+'.' + down[i] + '.'+up[i] +'.loc_d', loc_d)
#                 np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/IMR90-hMSC/'+chr[i]+'.' + down[i] + '.'+up[i] +'.loc_u', loc_u)
#                 np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/IMR90-hMSC/'+chr[i]+'.' + down[i] + '.'+up[i] +'.dlt', dlt)
#                 print(f1)
#     # except:
#     #     #print(i)
#     #     pass