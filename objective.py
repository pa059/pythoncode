# Bi bli o t h e k e n ei n bi n d e n
import numpy as np
import scipy.optimize as optimize
from funktions import LD

import matplotlib.pyplot as plt


 # Messdaten importieren
ref_ind = np.loadtxt("silicon-ref.txt")
arr= (ref_ind) .T
lam_data = arr[0] *1e3 #Wellenlänge 250−1450 n
n_data = arr[1] # Brechungsindex−Daten
k_data = arr[2] # E x ti n k s ti o n s −Daten
# Zu optimierende Funktion
#print(type(ref_ind))
#print(lam_data
print(lam_data)

def objective( x ) :
    # Berechnung von n und k
    epsilon = LD(lam_data , x)
    n = (epsilon ** 0.5).real
    k = (epsilon ** 0.5).imag

    # Summe de r q u a d r a ti s c h e n D i f f e r e n z von n
    delta_n = np . abs ( n - n_data )
    delta_n = np . square ( delta_n )
    delta_n = np . sum( delta_n )

    # Summe de r q u a d r a ti s c h e n D i f f e r e n z von k
    delta_k = np . abs ( k - k_data )
    delta_k = np . square ( delta_k )
    delta_k = np . sum( delta_k )
    
    # Summe de r gesamten q u a d r a ti s c h e n D i f f e r e n z
    delta_ges = delta_n + delta_k
    return delta_ges

    # Startwerte sind die Modellparameter von Gold
xGuess = np.array((9.030 , 0.760 , 0.053 , 0.0 , 0.241 , 0.415,
                   0.010 , 0.345 , 0.830 , 0.071 , 0.870 , 2.969,
                   0.601 , 2.494 , 4.304 , 4.384 , 2.214 , 13.32))

#model_params = model_params[~np.isnan(model_params)]

xGuess=xGuess[~np.isnan(xGuess)]
print(xGuess)

while(True) :
    min_result = optimize.minimize(objective , xGuess , method= 'Nelder-Mead')
    if objective ( xGuess ) > objective(min_result.x) :
        xGuess = min_result.x
        #print(xGuess)
        print(" Differenz : " , objective(min_result.x))
    else :
        xGuess = min_result.x
        print( " Keinen besseren Wert gefunden " )
        break
print(min_result.x)
np.savetxt(str( "LD-Modellparameter_Silicium.txt") , xGuess )

""" plt.plot(lam_data, n_data)
plt.plot(lam_data, k_data)
plt.show() """