'''
nucleic.thermo | thermo.py

Thermodynamic constants and calculations.

'''

import math

from nucleic.seq import reverse_complement
from nucleic.util import memoize


###############################################################################
# General constants and calculations
###############################################################################

R = 1.9872e-3  # Gas constant (1 / kcal * mol)
KELVIN = 273.15


def c_to_k(deg_c):
    ''' Convert degrees Celsius to degrees Kelvin ''' 
    return deg_c - KELVIN


def k_to_c(deg_k): 
    ''' Convert degrees Kelvin to degrees Celsius ''' 
    return deg_k + KELVIN 


def calc_ka(dg, deg_c):
    ''' Return the association constant at a given temperature '''
    return math.e**(-dg/(R * c_to_k(deg_c)))


def calc_dg(ds, dh, deg_c):
    ''' Return the dg at a given temp using the provided ds and dh '''
    return dh - ds * c_to_k(deg_c)


def calc_rand_coil(dg, deg_c):
    ''' Return the percent of randomly coiled oligo with dg at deg_c degrees '''
    return 1/(calc_ka(dg, deg_c) + 1)


###############################################################################
# DNA Thermodynamics
###############################################################################

# Pythia implementation based on SantaLucia2 2004
# @memoize
# def calc_dna_params(na_conc=50, mg_conc=0.01, dntp_conc=0.0):
#     ''' Return the thermo parameters for DNA under specified salt cond. 

#     '''
#     # SantaLucia 2004 NN data
#     enthalpies = {
#         'AA': -7.6, 'AT': -7.2, 'AG': -7.8, 'AC': -7.2,
#         'TA': -7.2, 'TT': -7.6, 'TG': -8.5, 'TC': -8.2,
#         'GA': -8.2, 'GT': -8.4, 'GG': -8.0, 'GC': -9.8,
#         'CA': -8.5, 'CT': -7.8, 'CG': -10.6, 'CC': -8.0
#     }
#     entropies = { 
#         'AA': -21.3, 'AT': -20.4, 'AG': -21.0, 'AC': -22.4,
#         'TA': -21.3, 'TT': -21.3, 'TG': -22.7, 'TC': -22.2,
#         'GA': -22.2, 'GT': -22.4, 'GG': -19.9, 'GC': -24.4,
#         'CA': -22.7, 'CT': -21.0, 'CG': -27.2, 'CC': -19.9
#     }
#     initiation_enthalpy =           0.2
#     initiation_entropy =           -5.7
#     terminal_at_enthalpy =          2.2
#     terminal_at_entropy =           6.9
#     symmetry_correction_enthalpy =  0.0 
#     symmetry_correction_entropy =  -1.4

#     # From Ahsen, Wittwer, and Schutz 2001
#     na_equivalents = (na_conc + 120. * math.sqrt(mg_conc - dntp_conc))/1000.
#     # From SantaLucia and Hicks, 2004
#     ds_correction = 0.368 * math.log(na_equivalents) / 2

#     thermo_params = {}
#     for k, v in enthalpies.items():
#         thermo_params[k] = {
#             'enthalpy': 1000. * v,
#             'entropy':  entropies[k] + ds_correction
#         }

#     thermo_params.update({
#         'initiation': {
#             'enthalpy':     1000. * initiation_enthalpy,
#             'entropy':      initiation_entropy
#         },
#         'terminal_at': {
#             'enthalpy':     1000. * terminal_at_enthalpy,
#             'entropy':      terminal_at_entropy
#         },
#         'symmetry_correction': {
#             'enthalpy':     1000. * symmetry_correction_enthalpy,
#             'entropy':      symmetry_correction_entropy
#         }
#     })

#     return thermo_params

# def calc_entropy_enthalpy(seq, thermo_params):
#     ''' Return the dS and dH of seq given thermo_params from calc_dna_params
#     '''
#     tp = thermo_params
#     entropy = tp['initiation']['entropy']
#     enthalpy = tp['initiation']['enthalpy']
    
#     for b in range(len(seq)-1):
#         entropy += tp[seq[b:b+2]]['entropy']
#         enthalpy += tp[seq[b:b+2]]['enthalpy']
#     if seq[0] == 'A' or seq[0] == 'T':
#         entropy += tp['terminal_at']['entropy']
#         enthalpy += tp['terminal_at']['enthalpy']
#     if seq[-1] == 'A' or seq[-1] == 'T':
#         entropy += tp['terminal_at']['entropy']
#         enthalpy += tp['terminal_at']['enthalpy']
#     if seq == reverse_complement(seq):
#         entropy += tp['symmetry_correction']['entropy']
#         enthalpy += tp['symmetry_correction']['enthalpy']
#     return entropy, enthalpy

# def calc_tm(seq, thermo_params, oligo_conc):
#     ''' Return the tm of seq.

#     thermo_params should be calculated by calc_dna_params and oligo_conc is in
#     Mol/L.

#     '''
#     entropy, enthalpy = calc_entropy_enthalpy(seq, thermo_params)
#     # Compute tm as per eq. 3 in Santalucia and Hicks, 2004 with additional
#     # compensation for primer concentration as per Cantor and Smith 1999.
#     return enthalpy/(entropy + 1.987 * math.log(oligo_conc)) - 273.15

# Python implementation of primer3 oligo_tm.c method

def divalent_to_monovalent(divalent, dntp):
    if divalent == 0:
        dntp = 0
    if divalent < dntp:
        divalent = dntp
    return 120 * math.sqrt(divalent-dntp)

def calc_thermo(seq, conc_nm=50, monovalent=50, divalent=0.01, dntp=0.0):
    ''' Return the thermo parameters for DNA under specified salt cond. 

    '''

    enthalpies = {
        'AA': 79, 'AT': 72, 'AG': 78, 'AC': 84,
        'TA': 72, 'TT': 79, 'TG': 85, 'TC': 82,
        'GA': 82, 'GT': 84, 'GG': 80, 'GC': 98,
        'CA': 85, 'CT': 78, 'CG': 106, 'CC': 80
    }
    entropies = { 
        'AA': 222, 'AT': 204, 'AG': 210, 'AC': 224,
        'TA': 213, 'TT': 222, 'TG': 227, 'TC': 222,
        'GA': 222, 'GT': 224, 'GG': 199, 'GC': 244,
        'CA': 227, 'CT': 210, 'CG': 272, 'CC': 199
    }
    dH = dS = 0
    # Calculate oligo symmetry
    sym = seq == reverse_complement(seq)
    # Calculate NN uncorrected dS and dH for oligo
    for idx in range(len(seq)-1):
        dH += enthalpies[seq[idx:idx+2]] 
        dS += entropies[seq[idx:idx+2]]
    # Terminal AT penalty and initiation parameters (combined)
    if seq[0] in 'AT':
        dH += -23
        dS += -41
    else:
        dH += -1
        dS += 28
    if seq[-1] in 'AT':
        dH += -23
        dS += -41
    else:
        dH += -1
        dS += 28
    if sym:
        dS += 14
    dH *= -100.0
    dS *= -0.1
    # Convert divalent salt and dntp conc. to monovalent equivalencies
    monovalent += divalent_to_monovalent(divalent, dntp)
    dS = dS + 0.368 * (len(seq) - 1) * math.log(monovalent / 1000.0)
    # Account for oligo symmetry and calculate tm
    if sym:
        tm = dH / (dS + 1.987 * math.log(conc_nm/1.0e9)) - KELVIN
    else:
        tm = dH / (dS + 1.987 * math.log(conc_nm/4.0e9)) - KELVIN
    return dH, dS, tm

def calc_tm(seq, conc_nm=50, monovalent=50, divalent=0.01, dntp=0.0):
    _, _, tm = calc_thermo(seq, conc_nm=conc_nm, monovalent=monovalent,
            divalent=divalent, dntp=dntp)
    return tm
