#%matplotlib inline
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib import rcParams
def plot_fit(w_exp,eps_exp, w_for_fit = None,eps_fit_result = None):
    '''
    This functions plots the experimental dielectric function as function of wavelength.
    If fitted model data is available it will also plot it
    '''
    #print("in")
    rcParams.update({'font.size': 18}) # Increase the font size
    rcParams['mathtext.default']='regular' # make the mathtext the same size of normal text for better readbility
    
    # Plot the real part of the dielectric function
    fig,ax = plt.subplots(1,2,figsize = (12,5))
    ax[0].scatter(w_exp, abs(eps_exp.real), marker = 'o',facecolor = 'none', edgecolor = 'g')
    
    ax[0].set_xlabel('Energy (eV)')
    ax[0].set_ylabel(r'|$\epsilon^\prime$|')
        

    # Plot the imaginary part of the dielectric function
    ax[1].scatter(w_exp, eps_exp.imag, marker = 'o',facecolor = 'none', edgecolor = 'g')
    ax[1].set_ylabel(r'$\epsilon ^ {\prime \prime}$')
    ax[1].set_xlabel('Energy (eV)')
    
    # If the fit data is available we will plot it with the actual plots
    if (w_for_fit is not None and eps_fit_result is not None):
        ax[0].plot(w_for_fit, abs(eps_fit_result.real),'-r')
        ax[1].plot(w_for_fit,eps_fit_result.imag,'-r')
        
    
    # Set the grid, log and axis formatter
    for axis in ax:
        axis.grid('on')
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.set_xlim([min(w_exp),max(w_exp)])
        for x_yaxis in [axis.xaxis, axis.yaxis]:
            x_yaxis.set_major_formatter(ScalarFormatter())
            
    fig.tight_layout()
    fig.savefig("evplot_si.png")
