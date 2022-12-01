import numpy as np
from scipy import special
from scipy import constants


def LD( lam,model_params ) :
    #print(type(model_params))
    omega_p = model_params[0]
    f0 = model_params[1]
    Gamma_0 = model_params[2]
    
    omega = (2* np.pi * constants.c) / (lam*1e-9)
    omega *= constants.hbar / constants.e # eV
    Omega_p = (f0 ** 0.5) * omega_p # added brackets  without bracket calculation error
     #Omega_p = f0 ** 0.5 *omega_p

    epsilon = 1-((Omega_p ** 2) /(omega*(omega + 1j *Gamma_0)))

    k = int ((len(model_params )-3)/3)
    for i in range (k) :
        f_i = model_params [3 + i*3]
        Gamma_i = model_params[4+i*3]
        omega_i = model_params[5+ i*3]

        epsilon+= f_i * (omega_p ** 2) / ((omega_i ** 2) - (omega ** 2)) - (1j * omega *Gamma_i)

    return epsilon

def psi_1d ( n , z ) :
    psi_n_z = special.spherical_jn (n ,z) * z
    d_psi_n_z = special.spherical_jn(n-1 , z ) * z - (n * special.spherical_jn (n , z ))
    return psi_n_z , d_psi_n_z

def xi_1d ( n , z ) :
    xi_n_z = (special . spherical_jn( n , z ) + 1j * special.spherical_yn( n , z ) ) * z
    d_xi_n_z = (special.spherical_jn (
    n - 1 , z ) + 1j * special.spherical_yn ( n-1 , z ) ) * z - n * (
    special.spherical_jn(n,z) + 1j * special.spherical_yn(n , z))
    return xi_n_z , d_xi_n_z



def psi_3d ( x , n ) :
    N = len (n)
    psi_n = np . zeros( ( x . shape [ 0 ] , x . shape [ 1 ] , N) , dtype=complex )
    d_psi_n = np . zeros ( ( x . shape [ 0 ] , x . shape [ 1 ] , N) , dtype=complex )
    for i in range(N) :
        psi_n [ : , : , i ] , d_psi_n [ : , : , i ] = psi_1d ( i + 1 , x )

    return psi_n , d_psi_n

def xi_3d ( x , n ) :
    N = len(n)
    xi_n = np . zeros ( ( x . shape [ 0 ] , x . shape [ 1 ] , N) , dtype=complex )
    d_xi_n = np . zeros (( x . shape [ 0 ] , x . shape [ 1 ] , N) , dtype=complex )

    for i in range(N) :
     xi_n [ : , : , i ] , d_xi_n [ : , : , i ] = xi_1d ( i + 1 , x )

    return xi_n , d_xi_n


def an_bn_3d(x,m,n):
    psi_n_x , d_psi_n_x = psi_3d ( x , n )
    psi_n_mx , d_psi_n_mx = psi_3d (m * x , n )
    xi_n_x , d_xi_n_x = xi_3d ( x , n )

    m = m.reshape ((m.shape[0] , m.shape[1] , 1))

    a_n = (m * psi_n_mx * d_psi_n_x - psi_n_x - d_psi_n_mx ) / (m * psi_n_mx *
    d_xi_n_x - xi_n_x * d_psi_n_mx )
    b_n = ( psi_n_mx * d_psi_n_x - m * psi_n_x * d_psi_n_mx ) / ( psi_n_mx * d_xi_n_x
-   m * xi_n_x * d_psi_n_mx )
    return a_n , b_n