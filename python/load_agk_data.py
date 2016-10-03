import numpy as np

data_dir = '/Volumes/yuguang_sd/data/agk/lin_response/lin_alf_1'
filename = data_dir + '/lin_alf_1.apar'

dtype = {'names': ('f1', 't', 'f2', 'apar2', 'f3', 'f4', 'h1', 'h2'),
         'formats': ('S2', 'f4', 'S10', 'f4', 'S10', 'S10', 'f4', 'f4')}
data = np.loadtxt(filename, dtype = dtype)

apar2 = np.array([elem[3] for elem in data])
apar = np.sqrt(apar2)
t = np.array([elem[1] for elem in data])

import matplotlib.pyplot as plt

plt.plot(t, apar)
plt.xlim([0, 5])
plt.show()

