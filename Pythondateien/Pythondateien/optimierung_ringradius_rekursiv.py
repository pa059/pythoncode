"""
Bachelorthesis
Titel: Berechnung einer Metaoberfläche mit fokussierenden Eigenschaften auf Basis der Mie-Theorie
Autor: Torben Hespe
Matrikel-Nr.: 5010294
Erstprüfer: Prof. Dr. C. Reinhardt
Zweitprüfer: Dr. D. Hilbig
Hochschule Bremen
Bearbeitungszeitraum: 29.03.2021 - 31.05.2021

Datei zur rekursiven Optimierung der Ringradien
"""

# Bibliotheken einbinden
import numpy as np
import pandas as pd
import funktionen
import matplotlib.pylab as plt

lam = 780 # Wellenlänge in nm
r0 = 100 # Partikelradius in nm

# Zu optimierenden Fokuspunkt festlegen
xfocus = np.array(0)
yfocus = np.array(0)
zfocus = np.array(14600)

delta_min = 1000 # Minimaler Abstand zweier Partikel
delta_max = 2000 # Maximaler Abstand zweier Partikel
stepsize = 50 # Schrittweite der Erhöhung der Ringradien

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

E0 = 1
E_n = 1j ** n * E0 * (2 * n + 1) / (n * (n + 1))

# Einfallende Welle
Ein_x = np.array(E0 * np.exp(1j * k * zfocus))
Ein_y = np.zeros_like(Ein_x)
Ein_z = np.zeros_like(Ein_x)
Eges = np.array([Ein_x, Ein_y, Ein_z])

# Berechnugn der Mie-Koeffizienten a_n und b_n
a_n, b_n = funktionen.an_bn_3d(x.reshape((1, 1)), m.reshape((1, 1)), n)
a_n = a_n.reshape(N, )
b_n = b_n.reshape(N, )

# Erstellen von Arrays
Ri = np.zeros(15, dtype=int)
Ri_best = np.zeros(15, dtype=int)
Esca_ring_best = np.zeros((15, 3), dtype=complex)

# Schleife zur Berechnung erzeugen
for i in range(len(Ri)):
    if i == 0:
        Ri[i] = delta_min
    elif Ri[i] == 0:
        Ri[i] = Ri_best[i - 1] + delta_min

    E_square_best = 0
    
    # Erhöhe Ri so lange, bis delta_max erreicht ist
    while (Ri[i] < (Ri_best[i - 1] + delta_max)):
        # Berechnung der Partikelpositionen auf Basis des Ringradius
        particles = funktionen.prtcls(np.array([Ri[i]]), delta_min)
        
        # Erstellen eines leeren Arrays für das gestreute Feld des Rings
        Esca_ring = np.zeros((3), dtype=complex)
        
        # Berechnung des Feldes für jeden Partikel
        for particle in particles:
            # Differenz zwischen Partikelposition und Fokuspunkt
            x_diff = xfocus - particle[0]
            y_diff = yfocus - particle[1]
            z_diff = zfocus - particle[2]
            
            # Berechnung des radialen Abstands
            r_diff = np.linalg.norm((x_diff, y_diff, z_diff), axis=0)
            
            # Berechnung der Winkel
            phi = np.arctan2(y_diff, x_diff)
            theta = np.arccos(z_diff / r_diff)

            # Berechnung von Rho
            rho = k * r_diff

            # Berechnung der Riccati-Bessel-Funktion und ihrer Ableitung
            xi_n_rho, d_xi_n_rho = funktionen.xi_1d(n, rho)
            
            # Berechnung der Vektorsphärischen Harmonischen
            M_o1n, N_e1n = funktionen.Mo1n_Ne1n(phi.reshape((1, 1)), theta.reshape((1, 1)), rho.reshape((1, 1)), n)
            M_o1n = M_o1n.reshape(3, N)
            N_e1n = N_e1n.reshape(3, N)
            
            # Berechnung des gestreuten elektrischen Feldes für den Partikel
            Esca_particle = np.sum((E_n * (1j * a_n * N_e1n - b_n * M_o1n)), axis=1)
            
            # Umwandlung in kartesische Basisvektoren
            Ex = np.cos(theta) * np.cos(phi) * Esca_particle[1] - np.sin(phi) * Esca_particle[2]
            Ey = np.cos(theta) * np.sin(phi) * Esca_particle[1] + np.cos(phi) * Esca_particle[2]
            Ez = -np.sin(theta) * Esca_particle[1]
            Exyz_particle = np.array([Ex, Ey, Ez])

            # Addieren des Feldes des Partikels zu dem Feld des Rings
            Esca_ring += Exyz_particle
        
        # Betragsquadrat bilden
        E_square = np.linalg.norm(np.abs(Eges + Esca_ring), axis=0) ** 2
        
        # Wenn Feldüberhöhung höher als bisher, nehme diesen Wert als optimalen Parameter an
        if (E_square > E_square_best):
            Ri_best[i] = Ri[i]
            Esca_ring_best[i] = Esca_ring
            E_square_best = E_square
        
        # Erhöhe Ringradius um stepsize
        Ri[i] += stepsize
    
    # Füge das am besten gestreute Feld dem Gesamtfeld hinzu
    Eges += Esca_ring_best[i]

# Berechnung des Betragsquadrats
E_square_ges = np.linalg.norm(np.abs(Eges), axis=0) ** 2

# Ausgabe der Ergebnisse
print("|E_ges|^2 =", E_square_ges)
print("Ringradien:", Ri_best)

#%%


fig2 = plt.figure(tight_layout=True)
ax2 = fig2.add_subplot(111)

particles = funktionen.prtcls(Ri_best, delta_min)

for i in range(len(particles)):
    circle1 = plt.Circle((particles[i,0], particles[i,1]), r0, color='k', fill=True)
    ax2.add_patch(circle1)
#ax.plot(particles[:,0],particles[:,1],'ko')
ax2.axis('equal')
ax2.set_box_aspect(1)
ax2.set_xlabel('$x\,/\,$nm')
ax2.set_ylabel('$y\,/\,$nm')

fig2.savefig("Linse_lam{}_r0{}_f{}.png".format(lam,r0,zfocus),dpi=300)
np.savetxt("Ringradien_lam{}_r0{}_f{}.txt".format(lam,r0,zfocus),Ri_best)
