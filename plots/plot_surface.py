import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

vec_exact = np.loadtxt("dataexact.txt")
vec_solve = np.loadtxt("datasolve.txt")

size_mat = int(np.sqrt(vec_exact.shape[0]))
x = np.linspace(0,1,size_mat)
X, Y = np.meshgrid(x,x)

def change(vector):
    size = int(np.sqrt(vector.shape[0]))
    x = np.empty([size, size])
    for i in range(0, size):
        for j in range(0,size):
            x[i,j] = vector[size*i+j]
    return x

Z = change(vec_solve)
#plt.figure()
fig = plt.figure()
ax = fig.gca(projection='3d')
#plt.axis('equal')
CS = ax.plot_surface(X, Y, Z)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('Calulated solution')
#plt.show()
plt.savefig("calc_sol_surface.png")
#plt.clabel(CS, inline=1, fontsize=10)
Z = change(vec_exact)
#plt.figure()
fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')
#plt.axis('equal')
CS= ax1.plot_surface(X, Y, Z)
#CS = plt.contour(X, Y, Z)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('Exact solution')
plt.savefig("exact_sol_surface.png")
#plt.show()
