import numpy as np

x = np.array([-12.5, -6.7, -2, -1.5, 0.1, 2.4, 6.8, 9.8, 15, 23.5, 30])
M1 = np.mean(x)
epsilon1 = 0.5 * np.sqrt(3) * (np.max(x) - np.min(x))
itermax = 30

M = np.zeros(itermax)
epsilon = np.zeros(itermax)

for j in range(itermax):
    counter = 0
    counter2 = 0
    name = 0
    name2 = 0
    seg = 0
    seg2 = 0

    if j == 0:
        epsilon[j] = epsilon1
        M[j] = M1
    else:
        for i in range(len(x)):
            seg = (x[i] - M[j - 1]) ** 2
            counter += 3 * (seg / ((epsilon[j - 1] ** 2) + seg) ** 2)
            name += 1 / ((epsilon[j - 1] ** 2) + seg) ** 2

        epsilon[j] = np.sqrt(counter / name)

        for i in range(len(x)):
            seg2 = (epsilon[j - 1] ** 2) / ((epsilon[j - 1] ** 2) + (x[i] - M[j - 1]) ** 2)
            counter2 += seg2 * x[i]
            name2 += seg2

        M[j] = counter2 / name2

print(M)
print(epsilon)

