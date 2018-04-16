#!/usr/bin/python

import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
import re
import sys

# some physical constants:

g = 9.81    # m.s^-2		gravitational acceleration 
t_0 = 15.0  # celcius            sea level temp in celcius
R_a = 287	# J.kg^-1.K^-1   gas constant for the atmosphere
R = 8.314       # J.mol^-1.K^-1  gas constant
p_0 = 101325	# Pa  		sea level pressure in Pa
N_a = 6.02214e23  # molecules.mol^-1     Avogadro's num

# our altitude grid up to 10km

altitudes = [10.000,9.800,9.600,9.400,9.200,9.000,8.800,8.600,8.400,8.200,8.000,7.800,7.600,7.400,7.200,7.000,6.800,6.600,6.400,6.200,6.000,5.800,5.600,5.400,5.200,5.000,4.800,4.600,4.400,4.200,4.000,3.800,3.600,3.400,3.200,3.000,2.800,2.600,2.400,2.200,2.000,1.800,1.600,1.400,1.200,1.000,0.800,0.600,0.400,0.200,0.000]

# a few function defs:

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


# open some files:

if len (sys.argv) < 5:
	raise IOError("The script takes one argument for the chosen SZA and four filename arguments: the output_map.inf and the block_amf.dat file for block amfs, and the output_map.inf and amf.dat file for the amfs. Make sure that the amfs and block amfs are from the same scenario.")

mysza = float(sys.argv[1])
blockinffile = sys.argv[2]
blockamffile = sys.argv[3]
amfinffile = sys.argv[4]
amffile = sys.argv[5]


blockinf = open(blockinffile, 'r')
blockamf = open(blockamffile, 'r')
inf = open(amfinffile, 'r')
amf = open(amffile, 'r')

if not blockinf:
	raise IOError("couldn't open the .inf file %s" % (inffile,))
if not blockamf:
	raise IOError("couldn't open the block_amf file %s" % (blockamffile,))
if not inf:
	raise IOError("couldn't open the .inf file %s" % (amfinffile,))
if not amf:
	raise IOError("couldn't open the amf file %s" % (amffile,))

# make some regular expressions:

# a regular expression to match lines that look like this (in both of the the .inf files):
#  Num.; SZA, LOS angle, azimuth angle (@ Observer position); output altitude
#   1   25.00    0.00   39.00    6.00

infpattern = re.compile(r"""\s+\d\s+
			(?P<sza>\d{2}[.]\d{2})\s+
			(?P<los>\d{1,2}[.]\d{2})\s+
			(?P<azi>\d{2}[.]\d{2})\s+
			(?P<alt>\d{1,2}[.]\d{2})
			.*""", re.VERBOSE)

szalist = []

for line in inf:
	m = infpattern.match(line)
	if m:
		thissza = float(m.group('sza'))
		szalist.append(thissza)

sza_index = szalist.index(mysza)

# a regular expression to match lines from the block_amf file
# that look like this:
#440.0000000                    NaN                    NaN                    NaN 
#               NaN                    NaN                    NaN            
#     NaN                    NaN
# 440.0000000   0.69395535559222D+00   0.70240597064792D+00   0.71128798130466D+00   0.72008562782143D+00   0.72819472951875D+00   0.73494712644128D+00   0.73962767685925D+00   0.73757052011323D+00

blockamfpattern = re.compile(r"""^\s+440.0000000\s
			(?P<blockamfstring>.*)""",re.VERBOSE)

# a regular expression to match lines from the amf file
# that look like this:
# ( we should only have one line matching from the amf file)

amfpattern = re.compile(r"""^\s+440.0000000\s
			(?P<amfstring>.*)""",re.VERBOSE)

# OK, start the logic here:

# get the amf for mysza:

matchcount = 0

for line in amf:
	m = amfpattern.match(line)
	if m:
		matchcount += 1
		amfstring = m.group('amfstring')
		amfs = amfstring.split()

myamf = amfs[sza_index]

# we should never get more than one line matching,  so raise
# an exception if we do:

if matchcount > 1:
	raise ValueError('Too many lines matched in the AMF file')

# get the block amfs for mysza:

blockamfs = []

for line in blockamf:
	m = blockamfpattern.match(line)
	if m:
		blockamfstring = m.group('blockamfstring')
		bamfs = blockamfstring.split()
		blockamfs.append(bamfs[sza_index])


# solve for the concentration:
vicd = np.vectorize(icd)
c_guess = 2.5e-08
target = 100.0000000000000e15

c_solution = fsolve(vicd,c_guess,(0,40000,target))

# we don't start at 0 and 1, so that when we don't go out of range when indexing altitudes[] and
# blockamfs[]
lower = 1
upper = 2


####### 
#we need to calculate the amf, so that we can normalise the block amfs,
#then multiply by the column fraction for each altitude interval


while (upper <= 51):

	altlower = float(altitudes[51-lower])*1000.0
	altupper = float(altitudes[51-upper])*1000.0
	altmean = altlower + (altupper - altlower)/2.0
	blockamflower = float(blockamfs[51-lower].replace('D','e'))
	#blockamfupper = float(blockamfs[51-upper].replace('D','e'))
	#blockamfmean = blockamflower + (blockamfupper - blockamflower)/2.0

	columnpercent = 100.0 * intcolumn(c_solution,altlower,altupper)/target

	contribution = blockamflower * columnpercent / float(myamf.replace('D','e'))

# 	the following gives single-line output
	sys.stdout.write("%.2f %.2f %.2f\n" % (altmean, blockamflower, contribution)) 
#	sys.stdout.write("%.3f %.2f %.2f %.2f %.2f\n" % (columnpercent, altlower, altupper, blockamflower, blockamfupper))
			

# the following to give multiline output for psxy:
#	sys.stdout.write(">\n")
	#sys.stdout.write(">-G100/100/100\n")
#	sys.stdout.write("0 %.2f\n" % (lower,))
#	sys.stdout.write("%.2f %.2f\n" % (columnpercent,lower) )
#	sys.stdout.write("%.2f %.2f\n" % (columnpercent,upper) )
#	sys.stdout.write("0 %.2f\n" % (upper,) )

	lower += 1
	upper += 1 





