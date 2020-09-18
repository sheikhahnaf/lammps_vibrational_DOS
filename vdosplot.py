import numpy as np
import matplotlib.pyplot as plt

filename=input("name of VDOS file-\n")
freq,x,y,z,avg=np.genfromtxt(filename,unpack=True,skip_header=2,usecols=(0,1,2,3,4),delimiter=' ')
plt.plot(freq,avg)
plt.xlim((0,20.0))
xx,locs = plt.xticks()
ll = ['%2.0f' % a for a in xx]
plt.xticks(xx, ll)
plt.xlabel('THz')
plt.ylabel(r'VDOS [THz $^{-1}$]')
plt.show()
