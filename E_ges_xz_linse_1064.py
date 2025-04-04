
# Bibliotheken einbinden
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import funktionen

from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'qt')

lam = 780  # Wellenlänge in nm
r0 = 100  # Partikelradius in nm

npts = 101  # Anzahl der Punkte je Achse
xrange = np.linspace(-6500, 6500, npts)  # Bereich der x-Ebene
yrange = 0  # 0, da Betrachtung in der xz-Ebene
zrange = np.linspace(0, 20000, npts)  # Bereich der x-Ebene

delta_min = 1000  # minimaler Abstand zweier Partikel

# Importieren der Lorentz-Drude-Modellparameter
model_data = pd.read_excel('LDdata.xls', index_col=0)
# Filtern der Parameter für Silicium
model_params = np.array(model_data["Si"])
model_params = model_params[~np.isnan(model_params)]

# Berechnung des komplexen Brechungsindexes
epsilon = funktionen.LD(lam, model_params)
n_s = epsilon ** 0.5  # kompl. Brechungsindex Partikel
n_m = 1  # kompl. Brechungsindex Umgebungsmedium
m = n_s / n_m  # Verhältnis der Brechungsindizes

# Berechnung von k und x
k = np.array(2 * np.pi * n_m / lam)
x = np.array(k * r0)

# Berechnung der maximalen Ordnung
N = int(2 + x + 4 * x ** (1 / 3))
n = np.arange(1, N + 1)

# 2-dimensionales Arrays der x-y- und z-Arrays erzeugen
xrange_2d, zrange_2d = np.meshgrid(xrange, zrange)
yrange_2d = np.zeros_like(xrange_2d)

E0 = 1
E_n = 1j ** n * E0 * (2 * n + 1) / (n * (n + 1))

# Einfallende Welle
Ein_x = np.array(E0 * np.exp(1j * k * zrange_2d))
Ein_y = np.zeros_like(Ein_x)
Ein_z = np.zeros_like(Ein_x)
Eges = np.array([Ein_x, Ein_y, Ein_z])

# Berechnung der Mie-Koeffizienten
a_n, b_n = funktionen.an_bn_3d(x.reshape((1, 1)), m.reshape((1, 1)), n)
a_n = a_n.reshape(N)
b_n = b_n.reshape(N)

# Ringradien der Linse
Ri = np.array([ 1350,  2350,  4300,  6250,  8200,  9750, 11100, 12350, 13550,
       14650, 15750, 16800, 17800, 18800, 20750])

# Berechnung der Partikelpositionen der Linse
particles = funktionen.prtcls(Ri, delta_min)

# Schleife zur Berechnung des Feldes für jeden Partikel
for particle in particles:
    # Differenz zwischen zu betrachtenden Punkten und Partikelpositionen
    x_diff = xrange_2d - particle[0]
    y_diff = yrange_2d - particle[1]
    z_diff = zrange_2d - particle[2]

    # Radiale Differenz Bestimmen
    r_diff = np.linalg.norm((x_diff, y_diff, z_diff), axis=0)
    r_diff[r_diff < r0] = r0  # Eliminiere Werte < r0, da Gleichungen für gestreutes Feld sonst nicht gelten

    # Berechnung der Winkel
    phi = np.arctan2(y_diff, x_diff)
    theta = np.arccos(z_diff / r_diff)

    # Berechnung von rho
    rho = k * r_diff

    # Berechnung der Riccati-Bessel-Funktionen
    xi_n_rho = funktionen.xi_3d(rho, n)

    # Berechnung der Vektorsphärischen Harmonischen
    M_o1n, N_e1n = funktionen.Mo1n_Ne1n(phi, theta, rho, n)

    # Berechnung des gestreuten elektrischen Feldes für den Partikel
    E_sca = np.sum((E_n * (1j * a_n * N_e1n - b_n * M_o1n)), axis=3)

    # Umwandlung in kartesische Basisvektoren
    Ex = np.cos(theta) * np.cos(phi) * E_sca[1, :, :] - np.sin(phi) * E_sca[2, :, :]
    Ey = np.cos(theta) * np.sin(phi) * E_sca[1, :, :] + np.cos(phi) * E_sca[2, :, :]
    Ez = -np.sin(theta) * E_sca[1, :, :]
    E_sca_xyz = np.array([Ex, Ey, Ez])

    # Addieren des gestreuten Feldes des Partikels zum Gesamtfeld
    Eges += E_sca_xyz

# Betragsquadrat des elektrischen Feldes berechnen
E_square = np.linalg.norm(np.abs(Eges), axis=0) ** 2

# Darstellung des Plots
fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)

extents = np.array([zrange[0], zrange[-1], xrange[0], xrange[-1]])
im = ax.imshow(E_square.T, cmap='seismic', interpolation="spline16", extent=extents, origin="lower", aspect="auto")  # Transponieren, um x auf der y-Achse zu setzen und z auf die x-Achse
cbar = fig.colorbar(im)
ax.grid(False)
ax.set_xlabel(r'$z\,/\,$nm')
ax.set_ylabel(r'$x\,/\,$nm')
cbar.set_label(r'$|E|^2$')

#fig.savefig("E_lam{}_r{}_f14600.png".format(lam,r0),dpi=300)
