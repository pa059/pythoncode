# To make a request and gather data,
import urllib.request as urllib2
import numpy as np
from io import StringIO
import callplotfit
# Query the website

response = urllib2.urlopen('https://refractiveindex.info/database/data/main/Si/Schinke.yml')

# Read the response. Response is in YAML format. Will use PyYAML parser,
yaml_data = response.read()
import yaml
data =  yaml.load(yaml_data,Loader=yaml.Loader)

# Lets see what is inside the data,
print ("Keys in the data:" , data.keys()) 

# # Lets print the reference for this data,
# print ("Reference : ", data['REFERENCES'])

# if 'COMMENTS' in data.keys():
#     print ("Comments : ", data['COMMENTS'])

# Data seems to be stored in 'DATA', if needed uncomment the following line
#print(data['DATA'])
# Some clean up needed. The required data string holding wavelength, n and k, is in the form of dictionary item 
# with key 'data', which is inside a list\n",

cleaned_nk = data['DATA'][0]['data'] 

# Use StringIO to convert String buffer as a file like class\n",
# Read the data into as numpy vectors\n",
wave_micron, n, k  = np.genfromtxt(StringIO(cleaned_nk), unpack=True)

wave_exp = wave_micron * 1000 # Convert wavelength from microns to nm
eps_exp = (n + 1j* k)**2  # Convert n,k into dielectric function
print(eps_exp)
# Lets convert wavelengh in nm to eV and assign it to 'w' 
h = 4.135667516E-15 # plancks's constant in eV-sec
c = 299792458E9 # speed of light in vacuum in nm/sec
w_exp = h*c/wave_exp;  # Convert wavleength in nm to eV
callplotfit.plot_fit(w_exp,eps_exp)