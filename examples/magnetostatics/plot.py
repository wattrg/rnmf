import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt
import sys

if len(sys.argv) == 2:
    iter = sys.argv[1]
else:
    iter = "psi_final"


psi = pd.read_csv("./examples/magnetostatics/" + iter + ".csv", header=None)
psi = psi.to_numpy()
Hy, Hx = np.gradient(psi)
Hx *= -1
Hy *= -1

x = np.linspace(0,200,201)
y = np.linspace(0,200,201)
X,Y = np.meshgrid(x,y)

plt.imshow(np.sqrt(Hx**2 + Hy**2))
plt.colorbar()
#plt.streamplot(X,Y,Hx,Hy)
plt.show()