
# Bibliotheken einbinden
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import funktions
import LDcalc as LD
import lmfit
from lmfit import minimize, Parameters, printfuncs
from scipy import constants

npts = 300  # Anzahl der Punkte je Achse
r0 = np.linspace(10, 350, npts)  # Partikelradiusbereich in nm
lam = np.linspace(250, 1250, npts)  # Wellenlängenbereich in nm

# Importieren der Lorentz-Drude-Modellparameter
model_data = pd.read_excel('LDdata_au.xlsx', index_col=0)
# Filtern der Parameter für Silicium
model_params = np.array(model_data["Au2"])
model_params = model_params[~np.isnan(model_params)]
params = Parameters()

params.add('Omega_p', value =model_params[0] )


params.add('f0', value=model_params[1]) 
params.add('Gamma0', value=model_params[2] )

params.add('f1', value=model_params[3]) 
params.add('Gamma1', value=model_params[4] ) 
params.add('Omega1', value =model_params[5]) 

params.add('f2', value=model_params[6]) 
params.add('Gamma2', value= model_params[7])
params.add('Omega2', value =model_params[8]) 

params.add('f3', value=model_params[9])  
params.add('Gamma3', value=model_params[10])
params.add('Omega3', value =model_params[11] )

params.add('f4', value=model_params[12])  
params.add('Gamma4', value= model_params[13])
params.add('Omega4', value =model_params[14] )
params.add('f5', value=model_params[15]) 
params.add('Gamma5', value=model_params[16] )
params.add('Omega5', value = model_params[17])

# 2-dimensionales Arrays der Eingangsparameter erzeugen
R0, LAM = np.meshgrid(r0, lam)
omega = 2 * np.pi * constants.c / (LAM * 1e-9)
omega *= constants.hbar / constants.e  # eV
# Berechnung des komplexen Brechungsindexes
#epsilon = funktions.LD(LAM, model_params )
epsilon = LD.drude_lorentz_model(params,omega)
#print(epsilon)
n_s = epsilon ** 0.5  # kompl. Brechungsindex Partikel
n_m = 1  # kompl. Brechungsindex Umgebungsmedium
m = n_s / n_m  # Verhältnis der Brechungsindizes

# Berechnung von k und x
k = np.array(2 * np.pi * n_m / LAM)
x = np.array(k * R0)

# Berechnung der maximalen Ordnung
N = int(np.max(2 + x + 4 * x ** (1 / 3)))
n = np.arange(1, N + 1)

# Berechnung der Mie-Koeffizienten a_n und b_n
a_n, b_n = funktions.an_bn_3d(x, m, n)

# Berechnung der Streueffizienz Q_sca
Q_sca = 2 / (k ** 2 * R0 ** 2) * np.sum(((2 * n + 1) * (np.abs(a_n) ** 2 + np.abs(b_n) ** 2)), axis=2)

# Darstellung des Plots
fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)

extents = np.array([r0[0], r0[-1], lam[0], lam[-1]])  #
im = ax.imshow(Q_sca, cmap='hot', interpolation='spline16', origin='lower', extent=extents, aspect='auto', vmin=2)
cbar = fig.colorbar(im)
ax.grid(False)
ax.set_xlabel(r'$r_0\,/\,$nm')
ax.set_ylabel(r'$\lambda\,/\,$nm')

fig.savefig("Q_sca_Au.png") 