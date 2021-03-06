import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns

sns.set()

# QDaedalus summer and autumn N-S vertical deflections
summer = np.array([-2.436, -2.298, -2.195, -2.430, -2.418, -2.361, -2.440, -2.271, -2.087, -2.407,
 -2.380, -2.313, -2.161, -2.362, -2.129, -2.328, -2.293, -2.341, -2.312, -1.973, -2.117, -2.219,
 -2.366, -2.219, -2.285, -2.297, -2.081, -2.123, -2.080, -2.386, -2.324, -2.373, -2.422, -2.324,
 -2.268, -2.356, -2.104, -2.373, -2.580, -2.389, -2.240, -2.255, -2.460, -2.169])

autumn = np.array([-2.790, -2.393, -3.156, -3.096, -2.230, -2.373, -2.946, -2.686, -2.189, -2.314,
 -2.805, -2.919, -2.630, -2.731, -2.771, -2.547, -2.310, -2.550, -2.757, -2.451, -2.572, -2.673,
 -2.866, -2.687, -2.582, -2.828, -2.802, -2.639, -2.390, -2.674, -2.473, -2.452, -2.725, -2.631,
 -2.921, -3.199, -2.594, -2.809, -2.530, -2.880, -3.090, -3.170, -2.650, -2.990, -2.842, -2.858,
 -2.670, -2.825, -2.690, -2.983, -2.851, -2.340])

plt.figure(figsize=(6,4))
su = np.sort(summer)
au = np.sort(autumn)

su = np.insert(su,0,np.amin(summer))
ndat1 = len(su)
ecx1 = np.linspace(0,1,ndat1)
# ECDF nyári
plt.step(su, ecx1, label="summer")
sns.set_style("darkgrid", {"axes.linewidth": ".5"})

au = np.insert(au,0,np.amin(autumn))
ndat2 = len(au)
ecx2 = np.linspace(0,1,ndat2)
# ECDF őszi
plt.step(au, ecx2, label="autumn")

plt.xlabel('N-S vertical deflection (")')
plt.legend()
plt.show()

#plt.savefig('xisa_ecdf_en.png',dpi=150, linewidth=0.3, color='r')
