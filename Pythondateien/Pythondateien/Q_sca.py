"""
Bachelorthesis
Titel: Berechnung einer Metaoberfläche mit fokussierenden Eigenschaften auf Basis der Mie-Theorie
Autor: Torben Hespe
Matrikel-Nr.: 5010294
Erstprüfer: Prof. Dr. C. Reinhardt
Zweitprüfer: Dr. D. Hilbig
Hochschule Bremen
Bearbeitungszeitraum: 29.03.2021 - 31.05.2021

Datei zur Berechnung der Streueffizienz in Abhängigkeit des Partikelradius und der Wellenlänge
"""

# Bibliotheken einbinden
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import funktionen

npts = 300  # Anzahl der Punkte je Achse
r0 = np.linspace(10, 350, npts)  # Partikelradiusbereich in nm
lam = np.linspace(250, 1250, npts)  # Wellenlängenbereich in nm

# Importieren der Lorentz-Drude-Modellparameter
model_data = pd.read_excel('LDdata.xls', index_col=0)
# Filtern der Parameter für Silicium
model_params = np.array(model_data["Si"])
model_params = model_params[~np.isnan(model_params)]

# 2-dimensionales Arrays der Eingangsparameter erzeugen
R0, LAM = np.meshgrid(r0, lam)

# Berechnung des komplexen Brechungsindexes
epsilon = funktionen.LD(LAM, model_params)
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
a_n, b_n = funktionen.an_bn_3d(x, m, n)

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

fig.savefig("Q_sca_Si.png")
