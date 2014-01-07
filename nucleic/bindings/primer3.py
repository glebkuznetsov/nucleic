'''
nucleic.bindings.primer3 | primer3.py

Low-level subprocess wrappers for the primer3 library. 

Requires that the primer3 executables ntthal, oligotm, primer3_core, ndtpal,
and long_seq_tm_test are present in your $PATH. Additionally, the 
primer3_config folder must be in one of the default locations (see primer3
documentation). The SETUP_UNIX and SETUP_WINDOWS shell scripts provided with
this package ensure that these files are in the proper locations.

'''

import os
import re
import subprocess

from collections import OrderedDict

from nucleic.util import TemporaryDirectory, unique_fn


TMP_DIR = TemporaryDirectory(persist=True).name
DEV_NULL = open(os.devnull, 'wb')


###############################################################################
# Defaults (do not modify)
###############################################################################

# Default low-level calculation parameters
_CALCULATION_PARAMETERS = {

    'DNA_CONC': {
        'FLAG':         'd',
        'UNITS':        'nM',
        'VAL':          50
    },
    'DNTP_CONC': {
        'FLAG':         'n',
        'UNITS':        'mM',
        'VAL':          0
    },
    'MONOVALENT_CONC': {
        'FLAG':         'mv',
        'UNITS':        'mM',
        'VAL':          50
    },
    'DIVALENT_CONC': {
        'FLAG':         'dv',
        'UNITS':        'mM',
        'VAL':          0.01
    },
    'THERMO_MODEL': {
        'FLAG':         'tp',
        'UNITS':        ('[0: Breslauer et al 1986,'
                         ' 1: SantaLucia 1998]'),
        'VAL':          1
    },
    'SALT_MODEL': {
        'FLAG':         'sc',
        'UNITS':        '[0: Schildkraut and Lifson 1965,'
                        ' 1: SantaLucia 1998,'
                        ' 2: Owczarzy et al. 2003]',
        'VAL':          1
    },
    'FOLDING_TEMP': {
        'FLAG':         't',
        'UNITS':        'C',
        'VAL':          37
    }

}


# Calculation parameters grouped by executable
_OLIGOTM_PARAMS = ('DNA_CONC', 'DNTP_CONC', 'MONOVALENT_CONC', 'DIVALENT_CONC',
                   'THERMO_MODEL', 'SALT_MODEL')
_NTTHAL_PARAMS = ('DNA_CONC', 'DNTP_CONC', 'MONOVALENT_CONC', 'DIVALENT_CONC',
                  'FOLDING_TEMP')


###############################################################################
# Config functions (use these to change the current parameters)
###############################################################################


def update_param(param, value):
    _CALCULATION_PARAMETERS[param]['VAL'] = value


def update_params(update_dict):
    for param, value in update_dict.items():
        update_param(param, value)


def _assemble_params(call_type):
    _flag = lambda k: '-' + cp[k]['FLAG']
    _val = lambda k: str(cp[k]['VAL'])
    cp = _CALCULATION_PARAMETERS
    pd = _OLIGOTM_PARAMS if call_type == 'oligotm' else _NTTHAL_PARAMS
    return [f(k) for k in pd for f in (_flag, _val)]


###############################################################################
# Low level bindings
###############################################################################


def calc_tm(seq):
    ''' Return the tm of `seq` as a float. '''
    args = ['oligotm'] + _assemble_params('oligotm') + [seq]
    out = subprocess.check_output(args, stderr=DEV_NULL)
    return float(out.strip())


_ntthal_re = re.compile(r'dS\s+=\s+(\S+)\s+dH\s+=\s+(\S+)\s+' +
                        r'dG\s+=\s+(\S+)\s+t\s+=\s+(\S+)')

def _parse_nthhal(ntthal_output):
    ''' Helper method that uses regex to parse ntthal output. '''
    parsed_vals = re.search(_ntthal_re, ntthal_output)
    return (
        float(parsed_vals.group(1)),    # dS
        float(parsed_vals.group(2)),    # dH
        float(parsed_vals.group(3)),    # dG
        float(parsed_vals.group(4))     # tm
    ) if parsed_vals else None


def calc_hairpin(seq):
    ''' Return a tuple of the dS, dH, dG, and Tm of any hairpin struct present.

    Returns None if the sequences does not form a hairpin returns None.

    '''
    args = ['ntthal', '-a', 'HAIRPIN'] + \
           _assemble_params('ntthal') + ['-s1', seq]
    out = subprocess.check_output(args, cwd=TMP_DIR, stderr=DEV_NULL)
    return _parse_nthhal(out)


def calc_heterodimer(seq, seq2):
    ''' Return a tuple of the dS, dH, dG, and Tm of any predicted heterodimer. 

    Returns None if the sequences do not form a heterodimer.

    '''
    args = ['ntthal'] + _assemble_params('ntthal') + ['-s1', seq, '-s2', seq2]
    out = subprocess.check_output(args, stderr=DEV_NULL, cwd=TMP_DIR)
    return _parse_nthhal(out)


def calc_homodimer(seq):
    ''' Return a tuple of the dS, dH, dG, and Tm of any predicted homodimer. 

    Returns None if the sequence does not form a homodimer.

    '''
    return calc_heterodimer(seq, seq)


def assess_oligo(seq):
    '''
    Return the thermodynamic characteristics of hairpin/homodimer structures.

    Returns a tuple tuples (hairpin stats, homodimer stats) in which each
    individual tuple is structured (dS, dH, dG, Tm).

    '''
    hairpin_out = calc_hairpin(seq)
    homodimer_out = calc_homodimer(seq)
    return (hairpin_out, homodimer_out)


def _write_boulderIO(fp, p3_args):
    with open(fp, 'wb') as fd:
        for k, v in p3_args.items():
            fd.write(k + '=' + v + '\n')
        fd.write('=')


def _parse_boulderIO(fp):
    data_dict = OrderedDict()
    with open(fp, 'rb') as fd:
        for line in fd:
            k,v = line.strip().split('=')
            data_dict[k] = v
    return data_dict


def run_p3_main(p3_args):
    ''' Return the raw primer3_core output for the provided primer3 args.

    Returns an ordered dict of the boulderIO-format primer3 output file
    '''
    fp_in = unique_fn(TMP_DIR, file_ext='.in')
    _write_boulderIO(fp_in, p3_args)
    fp_out = fp_in.rstrip('.in') + '.out'
    subprocess.check_call(['primer3_core', '--output', fp_out, fp_in],
                          stderr=DEV_NULL)
    return _parse_boulderIO(fp_out)

