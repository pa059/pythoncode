# Note that lmfit should be Version > 0.9.0. We will check this
import lmfit
from distutils.version import StrictVersion
assert StrictVersion(lmfit.__version__) > StrictVersion("0.9.0")
import comprefimp as com
import numpy as np
import matplotlib.pyplot as plt1
import residualcalc as rs
import LDcalc as LD
import callplotfit
from lmfit import minimize, Parameters, printfuncs

# Choose the model we want to fit it with 
model = 'Drude-Lorentz'

params = Parameters()


params.add('Omega_p', value = 30, min = 29.5, max = 30.69193516) # Omega_p has a value in eV, we will add a starting guess here

# Drude Term
params.add('f0', value=0 , min = 0, max = 0.1) # f  has no units,  we will add a starting guess here 
params.add('Gamma0', value= -9.99956E+306 , min = -9.99956E+307, max = -9.99956E+306 ) # Gamma is damping term and has units in eV, we will add a starting guess here

# Lorentz first Oscillator term
params.add('f1', value=0.104489936, min = 0.1, max = 0.2) # f  has no units,  we will add a starting guess here 
params.add('Gamma1', value= 1.803573504, min = 1.6 ,  max = 2) # Gamma is damping term and has units in eV, we will add a starting guess here
params.add('Omega1', value = -5.593217748 , min = -5.5 , max = -5.4) # Omega_o is oscillator central frequency term  and has units in eV. we will add a starting guess here

# Lorentz Second Oscillator term
params.add('f2', value=-1.886843639, min = -1.9, max = -1.7) # f  has no units,  we will add a starting guess here 
params.add('Gamma2', value= 0.9846178949, min = 0.9,  max = 1 ) # Gamma is damping term and has units in eV, we will add a starting guess here
params.add('Omega2', value = 3.668060488 , min = 3.5 , max = 3.7) # Omega_o is oscillator central frequency term  and has units in eV. we will add a starting guess here

# Lorentz Third Oscillator term
params.add('f3', value=0.01878305559, min = 0.01, max = 0.1) # f  has no units,  we will add a starting guess here 
params.add('Gamma3', value= 0.2145093945, min = 0.2 , max = 0.22) # Gamma is damping term and has units in eV, we will add a starting guess here
params.add('Omega3', value = 3.406791041 , min = 3.3, max = 3.5 ) # Omega_o is oscillator central frequency term  and has units in eV. we will add a starting guess here

# Lorentz Fourth Oscillator term
params.add('f4', value=0.1009430615, min = 0.1, max = 0.12) # f  has no units,  we will add a starting guess here 
params.add('Gamma4', value= 0.5624839897, min = 0.53 , max = 0.57) # Gamma is damping term and has units in eV, we will add a starting guess here
params.add('Omega4', value = -4.224179908 , min = -4.3, max = -4.1 ) # Omega_o is oscillator central frequency term  and has units in eV. we will add a starting guess here

# Lorentz Fifth Oscillator term
params.add('f5', value=1.900783515, min = 1.8, max = 2) # f  has no units,  we will add a starting guess here 
params.add('Gamma5', value= 0.9562973769, min = 0.93 , max = 1) # Gamma is damping term and has units in eV, we will add a starting guess here
params.add('Omega5', value = 3.671618555 , min = 3.5, max = 3.8 ) # Omega_o is oscillator central frequency term  and has units in eV. we will add a starting guess here


# Lets fit to all the data.
w_for_fit = com.w_exp
eps_for_fit = com.eps_exp

# if we were to fit only a part of the data such as fit only below 2 ev
#w_for_fit = w_exp[w_exp<2]
#eps_for_fit = eps_exp[w_exp<2]


# Call the minimize function with required parameters and use differential evolution as a method 
minimizer_results = minimize(rs.complex_residuals, params, args=(model, w_for_fit, eps_for_fit), method = 'differential_evolution', strategy='best1bin',
                             popsize=50, tol=0.01, mutation=(0, 1), recombination=0.9, seed=None, callback=None, disp=True, polish=True, init='latinhypercube')

# If we were to fit ith with least squares methods
# minimizer_results = minimize(complex_residuals, params, args=(model, w_for_fit, eps_for_fit), method = 'leastsq')
                             
#lets see whether the fit exited successfully?
print ("Print exited successfully? :  ", minimizer_results.success)

#lets see the termination status
print ("Termination Status: ", minimizer_results.message)

# lets print the fit report. We dont need lengthy Correlation table
printfuncs.report_fit(minimizer_results, show_correl=False)

print(minimizer_results.params)
# Caluclate the epsilon based on the fit results
eps_fit_result = np.array([LD.drude_lorentz_model(minimizer_results.params, i) for i in w_for_fit])
# eps=(eps_fit_result)**0.5;
# plt1.plot(com.wave_exp,eps.real);
# plt1.plot(com.wave_exp,eps.imag);
# plt1.show();

#print(eps)

#print(eps_fit_result)
# Lets plot the fit data
callplotfit.plot_fit(com.w_exp,com.eps_exp, w_for_fit,eps_fit_result)
