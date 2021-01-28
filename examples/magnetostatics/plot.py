import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt
import sys
import pyvista as pv

if len(sys.argv) == 2:
    iter = sys.argv[1]
else:
    iter = "_ferro_droplet_00000020"

mu_1 = 1
mu_2 = 3
Hy_far = 1.0
Hx_far = 0.0

data = pv.read("./examples/magnetostatics/" + iter + ".vti")
H = data.cell_arrays["h"].reshape((data.extent[1], data.extent[3],2))


H_analytical = np.sqrt(((3*mu_1)/(2*mu_1 + mu_2))**2*(Hy_far**2 + Hx_far**2))
print(f"H analytical = {H}, H computational = {np.sqrt(Hx[int(data.extent[1]/2),int(data.extent[1]/2)]**2 + Hy[int(data.extent[1]/2),int(data.extent[1]/2)]**2)}")

plt.figure()
plt.imshow(np.sqrt(Hx**2 + Hy**2))
plt.colorbar()
plt.show()