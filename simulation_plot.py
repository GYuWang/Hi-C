# import the necessary packages

import numpy as np
import math
import cv2





# A = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_merge.txt')
A = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_split.txt')


img = np.uint8(A)

img[img>27] = 27
img[img<2] = 0
img = img*(255/27)
img = img.astype(np.uint8)

color_img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)

color_img[:, :, 0] = 255 - color_img[:, :, 0]
color_img[:, :, 1] = 255 - color_img[:, :, 1]
color_img[:, :, 2] = 255



# f2 = open('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_merge_domaincaller.txt.DI.txt_HMM.txt.7col.txt.final.txt.band.txt')
f2 = open('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_split_domaincaller.txt.DI.txt_HMM.txt.7col.txt.final.txt.band.txt')
# f2 = open('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_split.txt.band.txt')
# f2 = open('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_merge.txt.band.txt')


band1 = np.empty([0, 2])
for line in f2:
    x = line.strip().split(' ')
    band1 = np.append(band1, [[math.ceil(float(x[1])), math.ceil(float(x[2]))]], axis=0)
f2.close()


print(band1)
# print(A.shape)

# band1 = np.delete(band1, (6), axis=0)

# band1 = np.delete(band1, (6), axis=0)
band1[8, 1] = 1050


for i in range(0, band1.shape[0]):
    color_img = cv2.rectangle(color_img, (int(band1[i,0]), int(band1[i,0])), (int(band1[i,1]), int(band1[i,1])), (0, 0, 0), 2, lineType=4)




# cv2.imwrite('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_domaincaller_merge.tiff', color_img)
cv2.imwrite('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_domaincaller_split.tiff', color_img)
# cv2.imwrite('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_our_split.tiff', color_img)
# cv2.imwrite('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_our_merge.tiff', color_img)
