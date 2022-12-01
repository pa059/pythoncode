# Note that lmfit should be Version > 0.9.0. We will check this
import lmfit
from distutils.version import StrictVersion
assert StrictVersion(lmfit.__version__) > StrictVersion("0.9.0")
import comprefimp as com
import numpy as np
import residualcalc as rs
import LDcalc as LD
import callplotfit
from lmfit import minimize, Parameters, printfuncs

# Choose the model we want to fit it with 
model = 'Drude-Lorentz'

params = Parameters()
params.add('Omega_p', value = 30.7, min = 30 , max = 31) # Omega_p has a value in eV, we will add a starting guess here

# Drude Term
params.add('f0', value=0, min = 0, max = 0.1) # f  has no units,  we will add a starting guess here 
params.add('Gamma0', value= float('-inf') , min = -1E+306, max = 0.1 ) # Gamma is damping term and has units in eV, we will add a starting guess here

# Lorentz first Oscillator term
params.add('f1', value=0.1, min = 0.05, max = 0.2) # f  has no units,  we will add a starting guess here 
params.add('Gamma1', value= 1.8, min = 1.5 ,  max = 2) # Gamma is damping term and has units in eV, we will add a starting guess here
params.add('Omega1', value = -5.6 , min = -5.5 , max = -2) # Omega_o is oscillator central frequency term  and has units in eV. we will add a starting guess here

# Lorentz Second Oscillator term
params.add('f2', value=-1.9, min = -2, max = -1) # f  has no units,  we will add a starting guess here 
params.add('Gamma2', value= 0.98, min = 0.96,  max = 1 ) # Gamma is damping term and has units in eV, we will add a starting guess here
params.add('Omega2', value = 3.6 , min = 3.5 , max = 4) # Omega_o is oscillator central frequency term  and has units in eV. we will add a starting guess here

# Lorentz Third Oscillator term
params.add('f3', value=0.0187, min = 0.009, max = 0.02) # f  has no units,  we will add a starting guess here 
params.add('Gamma3', value= 0.21, min = 0.1 , max = 0.3) # Gamma is damping term and has units in eV, we will add a starting guess here
params.add('Omega3', value = 3.4 , min = 3.3, max = 3.6 ) # Omega_o is oscillator central frequency term  and has units in eV. we will add a starting guess here

# Lorentz Fourth Oscillator term
params.add('f4', value=0.1, min = 0.09, max = 0.11) # f  has no units,  we will add a starting guess here 
params.add('Gamma4', value= .56, min = 0.54 , max = 0.58) # Gamma is damping term and has units in eV, we will add a starting guess here
params.add('Omega4', value = -4.2 , min = -4.5, max = -4.1 ) # Omega_o is oscillator central frequency term  and has units in eV. we will add a starting guess here

# Lorentz Fifth Oscillator term
params.add('f5', value=1.9, min = 1.8, max = 2) # f  has no units,  we will add a starting guess here 
params.add('Gamma5', value= 0.95, min = 0.97 , max = 1) # Gamma is damping term and has units in eV, we will add a starting guess here
params.add('Omega5', value = 3.6 , min = 3.5, max = 3.7 ) # Omega_o is oscillator central frequency term  and has units in eV. we will add a starting guess here



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

# Caluclate the epsilon based on the fit results
eps_fit_result = np.array([LD.drude_lorentz_model(minimizer_results.params, i) for i in w_for_fit])

# Lets plot the fit data
callplotfit.plot_fit(com.w_exp,com.eps_exp, w_for_fit,eps_fit_result)
