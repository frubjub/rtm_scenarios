#!/usr/bin/python

import math
import numpy as np
from scipy.integrate import fixed_quad
from scipy.integrate import quad
from scipy.integrate import romberg
from scipy.optimize import fsolve
import sys

aot = float(sys.argv[1])
ssa = float(sys.argv[2])
scaleheight = float(sys.argv[3])

g = 9.81    # m.s^-2		gravitational acceleration 
t_0 = 15.0  # celcius            sea level temp in celcius
R_a = 287	# J.kg^-1.K^-1   gas constant for the atmosphere
R = 8.314       # J.mol^-1.K^-1  gas constant
p_0 = 101325	# Pa  		sea level pressure in Pa
N_a = 6.02214e23  # molecules.mol^-1     Avogadro's num

def temperature(z):
	"""returns the temperature in Kelvin at height z (metres) in the standard atmosphere"""
	return 273.15 + (15.0 - 0.0065 * z)

def inv_scale_height(z):
	"""returns the inverse of scale height at height z in the standard atmosphere"""
	return g / (R_a * temperature(z))

def scale_height(z):
	"""returns the scale height at height z in the standard atmosphere"""
	return R_a * temperature(z) / g

def pressure(z):
	"""returns pressure in the standard atmosphere, units will depend on those of p_0"""
	return p_0 * math.exp(-1.0 * quad(inv_scale_height,0,z)[0])

def c(z,c0):
	"""returns an exponential profile of extinction with altitude"""
	return c0 * math.exp(-z/scaleheight)

def integrand(z,c0):
	"""expression aerosol extinction (integrand)"""
	return c(z,c0) 

def func(c0,z0,z1):
	"""gives VCD by integrating between z0 and z1"""
	vcd,err = quad(integrand,z0,z1,(c0,))
	return aot - vcd 

vfunc = np.vectorize(func)
c_guess = 0.2

c_solution = fsolve(vfunc,c_guess,(1400,20000))
#print c_solution

sys.stdout.write("Absorption coefficient\n")
sys.stdout.write("47\n")

step = 100
top = 6000

heights =  range(top,1400-step,-1*step)

# to get the absorption coefficient with single-scattering albedo
# of 0.9 we multiply c(h,c_solution) by 100.0 rather than 1000
# like we would do to get the extinction coefficient

for h in heights:
	sys.stdout.write("%.1f %.3e\n" % (h/1000.0, c(h,c_solution) * 1000.0 * (1.0 - ssa) ))
