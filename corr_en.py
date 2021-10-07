import numpy as np

data = np.loadtxt('values.txt', skiprows=1)
nd = data.shape[0]
nv = data.shape[1]
covmx = np.cov(data.T)
v = ['x','y']

print("Covariance matrix from Monte Carlo simulation")
print(" number of data: {:d}".format(nd))
print("covariance matrix:")
for i in range(nv):
    for j in range(nv):
        print("{0:2s}{1:2s}: {2:8.4e}".format(v[i],v[j],covmx[i,j]))



