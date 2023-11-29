import numpy as np
from matplotlib import pyplot as plt
from numpy import transpose, meshgrid
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext

def pinjie(spp_matrix_small,dimen5):
    [hang, lie, num] = spp_matrix_small.shape

    print('num',num)
    print('dimen5',dimen5)
    group = int(num/dimen5)
    dalie = lie*dimen5
    spp_matrix_big = np.zeros((hang, dalie, group), dtype=np.float_)
    # spp_matrix_mid = 10*np.ones((hang, dalie), dtype=np.float_)
    for j in range(group):
        for i in range(dimen5):
            # spp_matrix_mid[:,i] = np.append(spp_matrix_big,spp_matrix_small[:,:,j*dimen5+i],axis=1)
            spp_matrix_big[:,i*lie:(i+1)*lie,j] = spp_matrix_small[:, :, j*dimen5+i]
    # mm = np.max(spp_matrix_big)
    # print('max= ',mm) #不等于3就对了
    # print('size of small', spp_matrix_small.shape)
    # print('size of big' , spp_matrix_big.shape)
    return spp_matrix_big




def spp_plot(spp_matrix, fs, jiequshijian):
    [flen_half, tlen] = spp_matrix.shape
    # t = np.linspace(start=jiequshijian / tlen, stop=jiequshijian, num=tlen)
    # t = transpose(t)
    #
    # f = np.linspace(start=fs / 2 / flen_half, stop=fs / 2, num=flen)
    # f = transpose(f)
    # mapmap = meshgrid(t,f/1000,spp_matrix)

    #
    # plt.imshow(spp_matrix,origin='lower',vmin=0, vmax=0.3)
    # plt.imshow(spp_matrix, origin='lower')

    plt.pcolormesh(spp_matrix, vmin=0.05, vmax=1, norm =LogNorm())


    # plt.xticks(np.linspace(0, tlen-1, tlen), np.linspace(0, jiequshijian, tlen))
    # plt.yticks(np.linspace(0, flen_half-1, flen_half), np.linspace(0, fs/2, flen_half))
    # plt.imshow(X=spp_matrix,cmap=plt.cm.winter,interpolation='nearest')
    plt.ion()

    plt.colorbar(use_gridspec=True);


    plt.show()
    #


# x = np.array(range(1,4))
# y = np.array(range(1,3))
# z = np.array(range(1,13)).reshape(4,3)
# print(z.shape)
# plt.imshow(z,origin='lower')
# plt.xticks(np.linspace(0, 2, 3), np.linspace(1, 20, 3))
# plt.yticks(np.linspace(0, 3, 4), np.linspace(0, 30, 4))
# plt.colorbar()
# plt.ion()
# plt.show()
# plt.imshow(X)
# plt.show()
