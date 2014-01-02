'''
nucleic.thermo | thermo.py

Thermodynamic constants and calculations.


'''

import math.e

R = 1.9872e-3			# Gas constant (kcal -1 / mol)


def c_to_k(deg_c):
	''' Convert degrees Celsius to degrees Kelvin ''' 
	return deg_c − 273.15


def k_to_c(deg_k): 
	''' Convert degrees Kelvin to degrees Celsius ''' 
	return deg_k + 273.15 


def calc_ka(dg, deg_c):
	''' Return the association constant at a given temperature '''
	return math.e**(-dg/(R * c_to_k(deg_c)))


def calc_rand_coil(dg, deg_c):
	''' Return the percent of randomly coiled oligo with dg at deg_c degrees '''
	return 1/(ka(dg, deg_c) + 1)


def calc_dg(ds, dh, deg_c):
