import LDcalc
def complex_residuals(params, model, w, exp_data):
    '''
    This is the residual function that we will try to minimize.
    It takes the params dict that has the parameters that need to be found. 
   '''
    if model == 'Drude':
        # Our Drude model
        epsilon= drude_model(params,w)
    elif model == 'Drude-Lorentz':
        # Our Drude model
        epsilon= LDcalc.drude_lorentz_model(params,w)
        
    # Lets calculate our complex residual as the way it is done in page 5264
    # Rakic, a D.,et al (1998). Optical properties of metallic films for vertical-cavity optoelectronic devices. Applied Optics, 37(22), 5271–83.
    residual = (abs((epsilon.real - exp_data.real)/exp_data.real) + abs((epsilon.imag - exp_data.imag)/exp_data.imag))**2
    
    # if the residual is being used for least square optimizaiton we should have used
    # residual = (abs((epsilon.real - exp_data.real)/exp_data.real) + abs((epsilon.imag - exp_data.imag)/exp_data.imag))
    # least square method does the square of the residual in its algorithm. see http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.leastsq.html

    return residual