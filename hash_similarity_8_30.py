import cv2
import numpy as np
import imagehash
from PIL import Image



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