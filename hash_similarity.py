import cv2
import numpy as np
import imagehash
from PIL import Image


def dhash(image, hashSize=8):
    # resize the input image, adding a single column (width) so we
    # can compute the horizontal gradient
    resized = cv2.resize(image, (hashSize + 1, hashSize))

    # compute the (relative) horizontal gradient between adjacent
    # column pixels
    diff = resized[:, 1:] > resized[:, :-1]

    # convert the difference image to a hash
    return sum([2 ** i for (i, v) in enumerate(diff.flatten()) if v])

def read_img(file,M,N):

    f = open(file)
    #f = open('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/matrix/chr5_chr5_matrix_HUVEC_Coverage.txt')
    #f = open('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/5_matrix_IMR_Coverage.txt')

    l = f.readline().strip().split()
    A = np.array([l[M:(N-1)]],dtype=np.float32)


    for i,line in enumerate(f):
        if i < N and i >= M:
            l = line.strip().split()
            a = np.array(l,dtype=np.float32)
            #print(a[1:10])
            A = np.append(A, [a[M:(N-1)]], axis=0)
        else:
            if i >= N:
                break
            else:
                pass
    f.close()
    img = np.uint8(A[0:N,0:N])
    img[img>10] = 10
    img[img<1] = 0
    img = img*(255/10)
    img = img.astype(np.uint8)
    return img

# file1 = '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/matrix/chr5_chr5_matrix_HMEC_Coverage.txt'
# file2 = '/Users/guangyu/Work/Hi-C/Data/Contactmatrix/matrix/chr5_chr5_matrix_HUVEC_Coverage.txt'
# M = 677
# N = 785
#
# # M = 166
# # N = 488
#
# # M = 551
# # N = 659
#
# # M = 794
# # N = 1024
#
# # M = 1400
# # N = 1467
#
# # M = 1072
# # N = 1396
# #
# img1 = read_img(file1,M,N)
# img2 = read_img(file2,M,N)
#
# #
# # imageHash1 = dhash(img1)
# # imageHash2 = dhash(img2)


def hash_similarity(map1, map2, band):

    e_whash = []
    e_ahash = []
    e_phash = []
    e_dhash = []

    for i in band:
        n1 = int(i[1])
        n2 = int(i[2])
        x1 = map1[n1:(n2 + 1), n1:(n2 + 1)]
        x2 = map2[n1:(n2 + 1), n1:(n2 + 1)]
        x1[x1 > 10] = 10
        x1[x1 < 1] = 0
        x1 = x1 * (255 / 10)
        x1 = x1.astype(np.uint8)
        x2[x2 > 10] = 10
        x2[x2 < 1] = 0
        x2 = x2 * (255 / 10)
        x2 = x2.astype(np.uint8)
        img1 = Image.fromarray(x1)
        img2 = Image.fromarray(x2)
        hash1 = imagehash.whash(img1, 16)
        hash2 = imagehash.whash(img2, 16)
        d = 1 - (hash2 - hash1) / 256
        e_whash.append(d)
        hash1 = imagehash.average_hash(img1, 16)
        hash2 = imagehash.average_hash(img2, 16)
        d = 1 - (hash2 - hash1) / 256
        e_ahash.append(d)
        hash1 = imagehash.phash(img1, 16)
        hash2 = imagehash.phash(img2, 16)
        d = 1 - (hash2 - hash1) / 256
        e_phash.append(d)
        hash1 = imagehash.dhash(img1, 16)
        hash2 = imagehash.dhash(img2, 16)
        d = 1 - (hash2 - hash1) / 256
        e_dhash.append(d)

    e_whash = np.array(e_whash)
    e_whash = e_whash.reshape(e_whash.shape[0], 1)
    e_ahash = np.array(e_ahash)
    e_ahash = e_ahash.reshape(e_ahash.shape[0], 1)
    e_phash = np.array(e_phash)
    e_phash = e_phash.reshape(e_phash.shape[0], 1)
    e_dhash = np.array(e_dhash)
    e_dhash = e_dhash.reshape(e_dhash.shape[0], 1)
    band = np.concatenate((band, e_whash, e_ahash, e_phash, e_dhash), 1)
    return band


map1 = np.loadtxt("/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HMEC/5_matrix_HMEC_Coverage.txt")
map2 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HUVEC/5_matrix_HUVEC_Coverage.txt')
band = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/matrix/HUVEC_HMEC.txt')


hash_similarity(map1, map2, band)














band1 = np.loadtxt("/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HMEC/5_matrix_HMEC_Coverage.txt")
band2 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HUVEC/5_matrix_HUVEC_Coverage.txt')
TAD = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/matrix/HUVEC_HMEC.txt')


e = []

for i in TAD:
    n1 = int(i[1])
    n2 = int(i[2])
    x1 = band1[n1:(n2 + 1), n1:(n2 + 1)]
    x2 = band2[n1:(n2 + 1), n1:(n2 + 1)]
    x1[x1 > 10] = 10
    x1[x1 < 1] = 0
    x1 = x1 * (255 / 10)
    x1 = x1.astype(np.uint8)
    x2[x2 > 10] = 10
    x2[x2 < 1] = 0
    x2 = x2 * (255 / 10)
    x2 = x2.astype(np.uint8)
    img1 = Image.fromarray(x1)
    img2 = Image.fromarray(x2)
    hash1 = imagehash.whash(img1, 16)
    hash2 = imagehash.whash(img2, 16)
    d = 1 - (hash2-hash1)/256
    e.append(d)


e = np.array(e)
e = e.reshape(e.shape[0],1)
TAD = np.concatenate((TAD, e), 1)


np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/matrix/HUVEC_HMEC3.txt', TAD)




# hash1.hash
#
#
# #
# im1 = Image.fromarray(img1)
# im2 = Image.fromarray(img2)
#
# # hash1 = imagehash.average_hash(im1,16)
# # hash2 = imagehash.average_hash(im2,16)
#
# # hash1 = imagehash.dhash(im1)
# # hash2 = imagehash.dhash(im2)
#
# # hash1 = imagehash.phash(im1)
# # hash2 = imagehash.phash(im2)
#
# hash1 = imagehash.whash(im1,16)
# hash2 = imagehash.whash(im2,16)
#
#
# print(hash1,hash2)
# print(hash2 - hash1)
