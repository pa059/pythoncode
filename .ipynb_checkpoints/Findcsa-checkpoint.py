import numpy as np
import matplotlib.pylab as plt
#import pandas as pd
import funktionen


npts = 300 
r0 = np.linspace (10,250, npts) 
lam = np.linspace (200,1250,npts)

# Importieren der Lorentz−Drude−Modellparameter

model_data = np.loadtxt('LD-Modellparameter_Silicium.txt')

# Filtern der Parameter für Silicium
model_params = np.array(model_data)

model_params = model_params[~np.isnan(model_params)]

# 2−dimensionales Arrays der Eingans parametererzeugen
R0,LAM = np.meshgrid(r0,lam)

# Berechnung de s komplexen B rechun g sinde xe s
epsilon = funktionen.LD(LAM, model_params )
n_s = epsilon ** 0.5 # kompl . B rechung sindex P a r t i k e l
n_m = 1 # kompl . B rechung sindex Umgebungsmedium
m = n_s / n_m # Verhä l t n i s de r B r e c h u n g si n di z e s

# Berechnung von k und x
k = np.array(2 * np . pi * n_m / LAM)
x = np.array(k * R0)

# Berechnung de r maximalen Ordnung
N = int (np.max( 2 + x + 4 * x ** ( 1 / 3 )))
n = np.arange (1,N + 1 )

# Berechnung de r Mie−K o e f f i z i e n t e n a_n und b_n
a_n , b_n = funktionen.an_bn_3d( x , m, n )

# Berechnung de r S t r e u e f f i z i e n z Q_sca
c_sca = 2 / ( k ** 2 * R0 **2 ) * np.sum ((( 2 * n + 1 ) * ( np . abs ( a_n ) ** 2 + np . abs (
b_n) ** 2 ) ) , axis =2)
print(c_sca)

# Darstellung des Plots
fig = plt.figure(tight_layout=True )
ax = fig.add_subplot(111)

extents = np.array([r0[0] , r0[-1] ,lam [ 0 ] , lam [-1]])
im = ax.imshow (c_sca , cmap='hot'  , interpolation= 'spline16'  , origin='lower'  , extent=
extents , aspect='auto' , vmin=2)
cbar = fig.colorbar(im )
ax.grid(False)
ax.set_xlabel (r'$r_0 \ , / \ ,$nm')
ax.set_ylabel(r'$\lambda \ , / \ ,$nm'  )
plt.savefig('csa.png')
#plt.show
