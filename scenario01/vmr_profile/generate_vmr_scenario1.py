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
	"""returns vmr as a function of altitude, as a constant block up to 1500m"""
	if (z <= 1400.0):
		return c0
	elif (z > 1400.0):
		return 1.0e-11
#def c(z,c0):
#	"""returns an exponential profile of decreasing vmr with altitude, with scale height 1000"""
#	return c0 * math.exp(-z/1000)

def integrand(z,c0):
	"""expression for VCD integrand"""
	return ( c(z,c0) * N_a * pressure(z) ) / ( R * temperature(z) )

def icd(c0,z0,z1,targetvcd):
	"""gives integrated column density by integrating between z0 and z1"""
	# the factor 10000.0 is to convert molec.m^-2 to molec.cm^-2
	# target is the integrated column density we are aiming for
	# so the equation we solve is 0 = target - vcd
	icd,err = quad(integrand,z0,z1,(c0,))
	return targetvcd - (icd / 10000.0)

vicd = np.vectorize(icd)
c_guess = 5e-9
target = 2e16

c_solution = fsolve(vicd,c_guess,(0,40000,target))
#print c_solution

# print out the vertical profile:
# vertical profile above 10km is taken from climatology

climateprofile = {'60.60': 9.0000e-06,\
		'58.60': 9.0000e-06,\
		'56.50': 1.5000e-05,\
		'54.40': 2.8000e-05,\
		'52.20': 5.2000e-05,\
		'50.00': 1.0340e-05,\
		'47.70': 2.1150e-04,\
		'45.50': 4.4240e-04,\
		'43.30': 9.1800e-04,\
		'41.10': 1.7930e-03,\
		'39.00': 3.1200e-03,\
		'36.90': 4.5620e-03,\
		'34.80': 5.5540e-03,\
		'32.80': 5.7670e-03,\
		'30.90': 5.1500e-03,\
		'29.00': 4.0140e-03,\
		'27.10': 2.7200e-03,\
		'25.20': 1.7090e-03,\
		'23.40': 1.0170e-03,\
		'21.50': 5.5400e-04,\
		'19.70': 1.9260e-04,\
		'17.90': 5.2600e-05,\
		'16.10': 2.7100e-05,\
		'14.40': 1.8000e-05,\
		'12.50': 1.2190e-05,\
		'10.70': 1.0000e-10}

for h, v in reversed(sorted(climateprofile.iteritems())):
	sys.stdout.write("%.2f %.4e\n" % (float(h), float(v)))

#spacing = float(sys.argv[1])
#top = float(sys.argv[2])
step = 200
top = 10000

heights = range(top,0 - step,-1 * step)

for h in heights:
	sys.stdout.write("%.2f %.3e\n" % (float(h)/1000, c(h,c_solution)))


