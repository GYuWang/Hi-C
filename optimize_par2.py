# import the necessary packages
import numpy as np
import pandas as pd
import cv2
import math


def write_subset(f1, lim_down, lim_up):
    M = lim_down
    N = lim_up
    f = open(f1)
    l = f.readline().strip().split()
    A = np.array([l[M:(N - 1)]], dtype=np.float32)

    for i, line in enumerate(f):
        if i < N and i >= M:
            l = line.strip().split()
            a = np.array(l, dtype=np.float32)
            A = np.append(A, [a[M:(N - 1)]], axis=0)
        else:
            if i >= N:
                break
            else:
                pass
    f.close()
    img = np.uint8(A[0:N, 0:N])
    np.savetxt(f1 + '.' + str(lim_down) + '.' + str(lim_up) + '.new', img, delimiter='\t')

def plot_TAD(f1, f2, lim_down, lim_up):
    '''

    :param f1: contact map
    :param f2: band
    :param lim_down: color range down
    :param lim_up: color range up
    :return:
    '''
    img = np.loadtxt(f1)
    img[img > lim_up] = lim_up
    img[img < lim_down] = 0
    img = img * (255 / lim_up)
    img = img.astype(np.uint8)
    color_img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    color_img[:, :, 0] = 255 - color_img[:, :, 0]
    color_img[:, :, 1] = 255 - color_img[:, :, 1]
    color_img[:, :, 2] = 255
    band = open(f2)
    for line in band:
        x = line.strip().split('\t')

        if float(x[0]) > 0.8:
            color_img = cv2.rectangle(color_img, (math.ceil(float(x[1])), math.ceil(float(x[1]))),
                                      (math.ceil(float(x[2])), math.ceil(float(x[2]))), (0, 0, 0), 2,
                                      lineType=4)
            # color_img = cv2.rectangle(color_img, (math.ceil(float(x[1]))-4000, math.ceil(float(x[1]))-4000),
            #                           (math.ceil(float(x[2]))-4000, math.ceil(float(x[2]))-4000), (0, 0, 0), 2,
            #                           lineType=4)
    band.close()
    cv2.imwrite(f1 + '.tiff', color_img)


def main():
    length = pd.read_csv('/Users/guangyu/Work/Ref_data/hg19/hg19.fa.fai', sep='\t', header=None)
    for i in range(1,23):
        chr = 'chr' + str(i)
        s = length[length.iloc[:, 0] == chr].iloc[:, 1]
        s2 = int(math.ceil(s/10000000))
        for j in range(s2):
            down = j*1000
            up = (j+1)*1000
            write_subset('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HUVEC/' + str(i) + '_matrix_HUVEC_Coverage.txt', down, up)
            write_subset('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/' + str(i) + '_matrix_IMR90_Coverage.txt', down,
                     up)




if __name__== "__main__":
    main()




