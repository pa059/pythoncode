"""
Bachelorthesis
Titel: Berechnung einer Metaoberfläche mit fokussierenden Eigenschaften auf Basis der Mie-Theorie
Autor: Torben Hespe
Matrikel-Nr.: 5010294
Erstprüfer: Prof. Dr. C. Reinhardt
Zweitprüfer: Dr. D. Hilbig
Hochschule Bremen
Bearbeitungszeitraum: 29.03.2021 - 31.05.2021
"""

import numpy as np
from scipy import special
from scipy import constants
from scipy.special import logsumexp


def LD(lam, model_params):
    """
    Funktion zur Bestimmung der dielektrischen Funktion über das Lorentz-Drude Modell

    Parameters
    ----------
    lam : nd-array
        array der Wellenlänge in nm
    model_params : 1d-array
        Lorentz-Drude Modellparameter

    Returns
    -------
    epsilon : nd-array
        Dielektrische Funktion

    """
    omega_p = model_params[0]
    f0 = model_params[1]
    Gamma_0 = model_params[2]

    omega = 2 * np.pi * constants.c / (lam * 1e-9)
    omega *= constants.hbar / constants.e  # eV
    Omega_p = logsumexp(f0 ** 0.5) * logsumexp(omega_p)

    epsilon = 1 - Omega_p ** 2 / (omega * (omega + 1j * Gamma_0))

    k = int((len(model_params) - 3) / 3)
    for i in range(k):
        f_i = model_params[3 + i * 3]
        Gamma_i = model_params[4 + i * 3]
        omega_i = model_params[5 + i * 3]

        epsilon += f_i * omega_p ** 2 / ((omega_i ** 2 - omega ** 2) - 1j * omega * Gamma_i)

    return epsilon


def psi_1d(n, z):
    """
    Berechnung 2-dimensionaler Arrays für die Ricatti-Bessel-Funktion Psi und ihrer Ableitung

    Parameters
    ----------
    n : int
        Ordnung
    z : 2d-array
        Argument der Funktion

    Returns
    -------
    psi_n_z : 2d-array
        Riccati-Bessel-Funktion
    d_psi_n_z : 2d-array
        Ableitung der Riccati-Bessel-Funktion

    """
    psi_n_z = special.spherical_jn(n, z) * z
    d_psi_n_z = special.spherical_jn(n - 1, z) * z - n * special.spherical_jn(n, z)
    return psi_n_z, d_psi_n_z


def xi_1d(n, z):
    """
    Berechnung 2-dimensionaler Arrays für die Ricatti-Bessel-Funktion Xi und ihrer Ableitung

    Parameters
    ----------
    n : int
        Ordnung
    z : 2d-array
        Argument der Funktion

    Returns
    -------
    xi_n_z : 2d-array
        Riccati-Bessel-Funktion
    d_xi_n_z : 2d-array
        Ableitung der Riccati-Bessel-Funktion

    """
    xi_n_z = (special.spherical_jn(n, z) + 1j * special.spherical_yn(n, z)) * z
    d_xi_n_z = (special.spherical_jn(
        n - 1, z) + 1j * special.spherical_yn(n - 1, z)) * z - n * (
                       special.spherical_jn(n, z) + 1j * special.spherical_yn(n, z))
    return xi_n_z, d_xi_n_z


def psi_3d(x, n):
    """
    Berechnung 3-dimensionaler Arrays für die Riccati-Bessel-Funktion Psi ihre Ableitung

    Parameters
    ----------
    x : 2d-array
        Argument der Funktion
    n : 1d-array
        Array der Ordnungen

    Returns
    -------
    psi_n : 3d-array
        Riccati-Bessel-Funktion
    d_psi_n : 3d-array
        Ableitung der Riccati-Bessel-Funktion


    """
    N = len(n)
    psi_n = np.zeros((x.shape[0], x.shape[1], N), dtype=complex)
    d_psi_n = np.zeros((x.shape[0], x.shape[1], N), dtype=complex)

    for i in range(N):
        psi_n[:, :, i], d_psi_n[:, :, i] = psi_1d(i + 1, x)

    return psi_n, d_psi_n


def xi_3d(x, n):
    """
    Berechnung 3-dimensionaler Arrays für die Riccati-Bessel-Funktion Xi und ihre Ableitung

    Parameters
    ----------
    x : 2d-array
        Argument der Funktion
    n : 1d-array
        Array der Ordnungen

    Returns
    -------
    xi_n : 3d-array
        Riccati-Bessel-Funktion
    d_xi_n : 3d-array
        Ableitung der Riccati-Bessel-Funktion


    """
    N = len(n)
    xi_n = np.zeros((x.shape[0], x.shape[1], N), dtype=complex)
    d_xi_n = np.zeros((x.shape[0], x.shape[1], N), dtype=complex)

    for i in range(N):
        xi_n[:, :, i], d_xi_n[:, :, i] = xi_1d(i + 1, x)

    return xi_n, d_xi_n


def an_bn_3d(x, m, n):
    """
    Berechnung der Mie-Koeffizienten a_n und b_n

    Parameters
    ----------
    x : 2d-array
        Argument der Funktion
    m : 2d-array
        Verhältnis der komplexen Brechungsindizes
    n : 1d-array
        Array der Ordnungen

    Returns
    -------
    a_n : 3d-array
        Streukoeffizienz a_n
    b_n : 3d-array
        Streukoeffizienz b_n

    """
    psi_n_x, d_psi_n_x = psi_3d(x, n)
    psi_n_mx, d_psi_n_mx = psi_3d(m * x, n)
    xi_n_x, d_xi_n_x = xi_3d(x, n)

    m = m.reshape((m.shape[0], m.shape[1], 1))
    a_n = (m * psi_n_mx * d_psi_n_x - psi_n_x * d_psi_n_mx) / (m * psi_n_mx * d_xi_n_x - xi_n_x * d_psi_n_mx)
    b_n = (psi_n_mx * d_psi_n_x - m * psi_n_x * d_psi_n_mx) / (psi_n_mx * d_xi_n_x - m * xi_n_x * d_psi_n_mx)
    return a_n, b_n


def n_ptcl(Ri, delta_min):
    """
    Funktion zur Bestimmung der maximalen Anzahl an Partikeln je Ringradius

    Parameters
    ----------
    Ri : 1d-array
        Ringradien der Linse
    delta_min : float
        Minimaler Abstand zwischen zwei Partikeln

    Returns
    -------
    n_particle : 1d-array
        Anzahl an Partikeln je Ringradius

    """
    n_max = 2 * np.pi / np.arccos(1 - delta_min ** 2 / (2 * Ri ** 2))
    n_exponent = np.array(np.log2(n_max), dtype=int)
    n_particle = 2 ** n_exponent
    n_particle = np.asarray(n_particle).flatten()
    return n_particle


def prtcls(Ri, delta_min):
    """
    Funktion zur Bestimmung der (x,y,z)-Koordinaten aller Partikel auf Basis des Ringradius

    Parameters
    ----------
    Ri : 1d-array
        Ringradien der Linse
    delta_min : float
        Minimaler Abstand zwischen zwei Partikeln

    Returns
    -------
    particles : 2d-array
        Koordinaten aller Partikel auf der Linse

    """
    n = n_ptcl(Ri, delta_min)
    particles = np.zeros((0, 3))

    for i in range(np.size(Ri)):
        alpha = 2 * np.pi / n[i]
        p = np.zeros((n[i], 3))
        p[:, 0] = Ri[i] * np.sin(alpha * np.arange(0, n[i]))
        p[:, 1] = Ri[i] * np.cos(alpha * np.arange(0, n[i]))
        particles = np.concatenate((particles, p))
    return particles


def pi_tau_3d(theta, n):
    """
    Funktion zur Bestimmung der winkelabhängigen Funktionen pi_n und tau_n

    Parameters
    ----------
    theta : 2d-array
        Polarwinkel Theta
    n : 1d-array
        Array der Ordnungen

    Returns
    -------
    pi_n : 3d-array
        Streuwinkel pi_n
    tau_n : 3d-array
        Streuwinkel tau_n

    """
    N = len(n)
    pi_n = np.zeros((theta.shape[0], theta.shape[1], N))
    tau_n = np.zeros((theta.shape[0], theta.shape[1], N))

    pi_n[:, :, 0] = 1
    pi_n[:, :, 1] = 3 * np.cos(theta)

    for i in range(2, N):
        pi_n[:, :, i] = (2 * i + 1) / i * np.cos(theta) * pi_n[:, :, i - 1] - (i + 1) / i * pi_n[:, :, i - 2]

    tau_n[:, :, 0] = np.cos(theta)
    tau_n[:, :, 1:] = n[1:] * np.cos(theta.reshape(theta.shape[0], theta.shape[1], 1)) * pi_n[:, :, 1:] - (n[1:] + 1) * pi_n[:, :, :-1]

    return pi_n, tau_n


def Mo1n_Ne1n(phi, theta, rho, n):
    """
    Berechnung der Vektorsphärischen Harmonischen M_o1n und N_e1n

    Parameters
    ----------
    phi : 2d-array
        Azimutwinkel Phi
    theta : 2d-array
        Polarwinkel Theta
    rho : 2d-array
        Argument der Funktion
    n : 1d-array
        Array der Ordnungen

    Returns
    -------
    M_o1n : 4d-array
        Vektorspärische Harmonische M_o1n
    N_e1n : 4d-array
        Vektorspärische Harmonische N_e1n

    """
    N = len(n)

    pi_n, tau_n = pi_tau_3d(theta, n)

    xi_n, d_xi_n = xi_3d(rho, n)

    Mr = np.zeros((phi.shape[0], phi.shape[1], N))
    Mtheta = np.cos(phi).reshape((phi.shape[0], phi.shape[1], 1)) * pi_n * xi_n / rho.reshape(rho.shape[0], rho.shape[1], 1)
    Mphi = -np.sin(phi).reshape((phi.shape[0], phi.shape[1], 1)) * tau_n * xi_n / rho.reshape(rho.shape[0], rho.shape[1], 1)

    Nr = np.cos(phi).reshape((phi.shape[0], phi.shape[1], 1)) * n * (n + 1) * np.sin(theta).reshape((phi.shape[0], phi.shape[1], 1)) * pi_n * xi_n / rho.reshape(rho.shape[0], rho.shape[1], 1) ** 2
    Ntheta = np.cos(phi).reshape((phi.shape[0], phi.shape[1], 1)) * tau_n * d_xi_n / rho.reshape(rho.shape[0], rho.shape[1], 1)
    Nphi = -np.sin(phi).reshape((phi.shape[0], phi.shape[1], 1)) * pi_n * d_xi_n / rho.reshape(rho.shape[0], rho.shape[1], 1)

    M_o1n = np.array([Mr, Mtheta, Mphi])
    N_e1n = np.array([Nr, Ntheta, Nphi])

    return M_o1n, N_e1n
