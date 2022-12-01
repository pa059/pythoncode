def drude_lorentz_model(params, w):
    '''
    This is functional representation of drude lorentz model with the first term being Drude
    and the rest terms being lorentz oscillators
    '''
    
    # Get the plasma frequency
    omega_p = params['Omega_p'].value
   
    
    # Get the Drude parms
    f0 = params['f0'].value
    gamma0 = params['Gamma0'].value
    print(gamma0)

    # Get the first Lorentz term parameters,
    f1 = params['f1'].value
    gamma1 = params['Gamma1'].value
    omega1 = params['Omega1'].value
    
    # Get the second Lorentz term parameters
    f2 = params['f2'].value
    gamma2 = params['Gamma2'].value
    omega2 = params['Omega2'].value
    
    # Get the third Lorentz term parameters
    f3 = params['f3'].value
    gamma3 = params['Gamma3'].value
    omega3 = params['Omega3'].value

    # Get the fourth Lorentz term parameters
    f4 = params['f4'].value
    gamma4 = params['Gamma4'].value
    omega4 = params['Omega4'].value
    
    # Get the fifth Lorentz term parameters
    f5 = params['f5'].value
    gamma5 = params['Gamma5'].value
    omega5 = params['Omega5'].value

    # Drude component
    #epsilon_D = 1 - omega_p ** 2 / (w * (w + 1j * gamma0))
    epsilon_D = 1 - (f0 * omega_p ** 2 / (w ** 2 + 1j * (gamma0) * w))

    # Lorentz first oscillator
    epsilon_L1 = (f1 * omega_p ** 2) / (omega1 ** 2 - w ** 2 - 1j * gamma1 * w)

    # Lorentz SEcond oscillator 
    epsilon_L2 = (f2 * omega_p ** 2) / (omega2 ** 2 - w ** 2 - 1j * gamma2 * w)

    # Lorentz Third oscillator 
    epsilon_L3 = (f3 * omega_p ** 2) / (omega3 ** 2 - w ** 2 - 1j * gamma3 * w)
    
    # Lorentz Fourth oscillator 
    epsilon_L4 = (f4 * omega_p ** 2) / (omega4 ** 2 - w ** 2 - 1j * gamma4 * w)
    
    # Lorentz Fifth oscillator 
    epsilon_L5 = (f5 * omega_p ** 2) / (omega5 ** 2 - w ** 2 - 1j * gamma5 * w)

    # Sum all the terms
    epsilon = epsilon_D + epsilon_L1 +epsilon_L2 +  epsilon_L3 +  epsilon_L4 +  epsilon_L5

    return epsilon