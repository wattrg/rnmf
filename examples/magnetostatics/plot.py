import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt

psi = pd.read_csv("./examples/magnetostatics/test.csv", header=None)
psi = psi.to_numpy()
Hy, Hx = np.gradient(psi)
Hx *= -1
Hy *= -1
print(Hy)
plt.contourf(np.sqrt(Hx**2 + Hy**2))
plt.quiver(Hx, Hy)
plt.show()