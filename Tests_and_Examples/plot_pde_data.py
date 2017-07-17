

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
file=open('pde_data1')
X_MIN,X_MAX,X_STEPS=map(int,file.readline().strip().split(' '))
Y_MIN,Y_MAX,Y_STEPS=map(int,file.readline().strip().split(' '))
Z_MIN,Z_MAX=map(float,file.readline().strip().split(' '))
X_STEPSIZE=(X_MAX-X_MIN)/float(X_STEPS)
Y_STEPSIZE=(Y_MAX-Y_MIN)/float(Y_STEPS)
X = np.linspace(X_MIN,X_MAX,X_STEPS+1)
Y = np.linspace(Y_MIN,Y_MAX,Y_STEPS+1)
X, Y = np.meshgrid(X, Y)
Z=np.empty((Y_STEPS+1,X_STEPS+1))
for i in range(Y_STEPS+1):
    Z[i]=map(float,file.readline().strip().split(' '))

file.close()

# Plot the surface.
surf = ax.plot_wireframe(X, Y, Z)

# customize axes
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_zlim(Z_MIN, Z_MAX)


plt.show()
