import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt
import sys

if len(sys.argv) == 2:
    iter = sys.argv[1]
else:
    iter = "psi_final"

mu_1 = 1
mu_2 = 3
Hy_far = 1.0
Hx_far = 0.0

psi = pd.read_csv("./examples/magnetostatics/" + iter + ".csv", header=None)
psi = psi.to_numpy()
Hy, Hx = np.gradient(psi)
Hx *= -1
Hy *= -1


H = np.sqrt(((3*mu_1)/(2*mu_1 + mu_2))**2*(Hy_far**2 + Hx_far**2))
print(f"H analytical = {H}, H computational = {np.sqrt(Hx[125,125]**2 + Hy[125,125]**2)}")

plt.figure()
plt.imshow(np.sqrt(Hx**2 + Hy**2))
plt.colorbar()
plt.show()