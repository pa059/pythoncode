import numpy as np
from scipy import special
from scipy import constants


def LD( lam,model_params ) :
    #print(type(model_params))
    omega_p = model_params[0]
    f0 = model_params[1]
    Gamma_0 = model_params[2]
    omega = 2 * np . pi * constants. c / (lam*1e-9)
    omega *= constants.hbar / constants.e # eV
    Omega_p = (f0 ** (0.5 *omega_p)) # added brackets  without brackets calculation error occur
     #Omega_p = f0 ** 0.5 *omega_p

    epsilon = 1-Omega_p ** 2 /(omega*(omega + 1j *Gamma_0))

    k = int ((len(model_params )-3)/3)
    for i in range (k) :
        f_i = model_params [3 + i*3]
        Gamma_i = model_params[4+i*3]
        omega_i = model_params[5+ i*3]

        epsilon+= f_i * (omega_p ** 2) / (((omega_i ** 2) - (omega ** 2) ) - (1j * omega *Gamma_i))

    return epsilon