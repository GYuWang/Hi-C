import numpy as np


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
        # print('---------------')
        # print([t, is_t1(t, TAD)])
        if is_t1(t, TAD) == 1:
            t1 = np.append(t1, [t], axis=0)
    return t1

def get_ratio(map, TAD):
    def neighbor_TAD(TAD):
        '''
        find each pair of nearby TAD
        '''
        pair = []
        for i in range(TAD.shape[0]):
            for j in range(TAD.shape[0]):
                if abs(TAD[i,2] - TAD[j, 1]) < 20:
                    pair.append([TAD[i], TAD[j]])
        return np.array(pair)

    def calculate_ratio(map, pair):
        ratio = []
        for T in pair:
            t1 = T[0]
            t2 = T[1]
            m1 = np.mean(map[int(t1[1]):int(t1[2]), int(t1[1]):int(t1[2])])
            m2 = np.mean(map[int(t2[1]):int(t2[2]), int(t2[1]):int(t2[2])])
            m3 = np.mean(map[int(t1[1]):int(t1[2]), int(t2[1]):int(t2[2])])
            r = m3/(m1 + m2)
            ratio.append(r)
        return ratio

    pair = neighbor_TAD(TAD)
    ratio = calculate_ratio(map, pair)

    return ratio




TAD = np.loadtxt("/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/pipline/output/test.HUVEC.txt")
# TAD = np.loadtxt("/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/pipline/output/test.IMR90.txt")

map = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/HUVEC/new_file/8_matrix_HUVEC_Coverage.txt.10000.11000.new')
# map = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file/8_matrix_IMR90_Coverage.txt.10000.11000.new')


TAD_t1 = find_t1(TAD)

ratio = get_ratio(map, TAD)
print(ratio)
print(np.mean(ratio))

TAD = np.loadtxt("/Users/guangyu/Work/Hi-C/Data/Contactmatrix/differential/pipline/output/test.IMR90.txt")

map = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/IRM90/matrix/new_file/8_matrix_IMR90_Coverage.txt.10000.11000.new')


ratio = get_ratio(map, TAD)
print(ratio)
print(np.mean(ratio))
