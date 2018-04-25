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
		return 0.0
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

def intcolumn(c0,z0,z1):
	icd,err = quad(integrand,z0,z1,(c0,))
	return icd/10000.0  # 10000.0 to convert the units!

vicd = np.vectorize(icd)
c_guess = 2.5e-08
target = 100.0000000000000e15

c_solution = fsolve(vicd,c_guess,(0,40000,target))
#print c_solution

step = 200
top = 10000

#heights = range(top,0 - step,-1 * step)

# the lower and upper limits of the sub-column that we integrate over
lower = 0
upper = lower + step


while (upper <= top):

	columnpercent = 100.0 * intcolumn(c_solution,lower,upper)/target
# 	the following gives single-line output
	sys.stdout.write("%.2f %.2f\n" % (lower, columnpercent))

# the following to give multiline output for psxy:
#	sys.stdout.write(">\n")
	#sys.stdout.write(">-G100/100/100\n")
#	sys.stdout.write("0 %.2f\n" % (lower,))
#	sys.stdout.write("%.2f %.2f\n" % (columnpercent,lower) )
#	sys.stdout.write("%.2f %.2f\n" % (columnpercent,upper) )
#	sys.stdout.write("0 %.2f\n" % (upper,) )

	lower = upper
	upper = lower + step





