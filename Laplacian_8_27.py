import numpy as np
from scipy.sparse.linalg import eigsh



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
    i_z1 = np.where(M1b.sum(1) == 0)[0]
    i_z2 = np.where(M2b.sum(1) == 0)[0]

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

    N = float(M1.shape[1]);
    # if (np.sum(ipr1 > N / 100) <= 1) | (np.sum(ipr2 > N / 100) <= 1):
    #     print("at least one of the maps does not look like typical Hi-C maps")
    # else:
    #     print("size of maps: %d" % (np.size(M1, 0)))
    #     print("reproducibility score: %6.3f " % (evs))
    #     print("num_evec_eff: %d" % (num_evec_eff))
    return evs


def similarity(f1, f2, TAD_f):
    e = []
    m1 = []
    St1 = []
    m2 = []
    St2 = []
    band = np.loadtxt(f1)
    band2 = np.loadtxt(f2)

    for i in TAD_f:
        n1 = int(i[1])
        n2 = int(i[2])
        x1 = band[n1:(n2 + 1), n1:(n2 + 1)]
        x2 = band2[n1:(n2 + 1), n1:(n2 + 1)]
        # x1[x1 > 10] = 10
        # x1[x1 < 1] = 0
        # x1 = x1 * (255 / 10)
        # x1 = x1.astype(np.uint8)
        # x2[x2 > 10] = 10
        # x2[x2 < 1] = 0
        # x2 = x2 * (255 / 10)
        # x2 = x2.astype(np.uint8)

        e.append(get_reproducibility(x1, x2))

    for i in TAD_f:
        n1 = int(i[1])
        n2 = int(i[2])
        x1 = band[n1:(n2 + 1), n1:(n2 + 1)]
        x2 = band2[n1:(n2 + 1), n1:(n2 + 1)]
        x1[x1 > 10] = 10
        x1[x1 < 1] = 0
        x1 = x1 * (255 / 10)
        x1 = x1.astype(np.uint8)
        x2[x2 > 10] = 10
        x2[x2 < 1] = 0
        x2 = x2 * (255 / 10)
        x2 = x2.astype(np.uint8)
        m1.append(x1.mean())
        St1.append(x1.var())
        m2.append(x2.mean())
        St2.append(x2.var())


    e = np.array(e)
    e = e.reshape(e.shape[0],1)
    m1 = np.array(m1)
    m1 = m1.reshape(m1.shape[0], 1)
    St1 = np.array(St1)
    St1 = St1.reshape(St1.shape[0], 1)
    m2 = np.array(m2)
    m2 = m2.reshape(m2.shape[0], 1)
    St2 = np.array(St2)
    St2 = St2.reshape(St2.shape[0], 1)
    TAD_f = np.concatenate((TAD_f, e, m1, m2, St1, St2), 1)
    return TAD_f



def merge(TAD1, TAD2):
    TAD_merge = np.append(TAD1,TAD2, 0)
    return TAD_merge



# f1 = "/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_merge.txt"
# f2 = "/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_split.txt"
#
# TAD1 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_merge.txt.band.txt')
# TAD2 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_split.txt.band.txt')
# TAD_all = merge(TAD1, TAD2)
#
# TAD = similarity(f1,f2,TAD_all)
# np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation2_laplacian.txt', TAD)
#

for i in range(2,27):
    f1 = "/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation" + str(i) + "_merge.txt"
    f2 = "/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation" + str(i) + "_split.txt"

    TAD1 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation'
                      + str(i) + '_merge.txt.band.txt')
    TAD2 = np.loadtxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation'
                      + str(i) + '_split.txt.band.txt')
    TAD_all = merge(TAD1, TAD2)

    TAD = similarity(f1, f2, TAD_all)
    np.savetxt('/Users/guangyu/Work/Hi-C/Data/Contactmatrix/simulation/figure2/simulation'+
               str(i) + '_laplacian.txt', TAD)






