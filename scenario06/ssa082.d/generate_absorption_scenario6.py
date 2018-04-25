#!/usr/bin/python

# Copyright (C) 2018 Stephen Broccardo

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import math
import numpy as np
from scipy.integrate import fixed_quad
from scipy.integrate import quad
from scipy.integrate import romberg
from scipy.optimize import fsolve
import sys

aot = float(sys.argv[1])
ssa = float(sys.argv[2])

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
	"""returns concentration as a function of altitude, as a constant"""
	if (z <= 2800.0):
		return c0
	elif (z > 2800.0 and z < 3600.0 ):
		return 0.0
	elif (z >= 3600.0 and z <= 3800.0 ):
		return c0 / 5.0
	else:
		return 0.0

#def c(z,c0):
#	"""returns an exponential profile of altitude"""
#	return c0 * math.exp(-z/1000)

#def integrand(z,c0):
#	"""expression for vertical column density integrand"""
#	return (c(z,c0) * N_a)

def integrand(z,c0):
	"""expression for aerosol extincion  (integrand)"""
	return c(z,c0)

def func(c0,z0,z1):
	"""gives AOT by integrating between z0 and z1"""
	vcd,err = quad(integrand,z0,z1,(c0,))
	return aot - vcd 

vfunc = np.vectorize(func)
c_guess = 0.2

c_solution = fsolve(vfunc,c_guess,(1400,5000))
#print c_solution


sys.stdout.write("Absorption coefficient\n")
sys.stdout.write("16\n")

step = 200
top = 4000

heights =  range(top,1200,-1*step)

for h in heights:
	sys.stdout.write("%.2f %.2e\n" % (h/1000.0, c(h,c_solution)*1000 * (1.0-ssa)))


