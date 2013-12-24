'''
nucleic.bindings.unafold | unafold.py

Python bindings for simple UNAfold analysis of oligo sequences

Requires that the UNAfold executables are present in your $PATH. Installing
UNAfold via the Makefile included in its package should be sufficient.

********* This will soon be replaced by C API Nupack bindings *********

'''

import os
import re

from math import log10
from subprocess import check_call

from nucleic.util import TemporaryDirectory, unique_fn

TMP_DIR = TemporaryDirectory(persist=True).name
DEV_NULL = open(os.devnull, 'wb')


###############################################################################
# Defaults (do not modify)
###############################################################################

# Default calculation parameters
_CALCULATION_PARAMETERS = {

    'NA_CONC': {
        'FLAG':         'N',
        'UNITS':        'M',
        'VAL':          0.075
    },
    'MG_CONC': {
        'FLAG':         'M',
        'UNITS':        'M',
        'VAL':          0.003
    },
    'NA_TYPE': {
        'FLAG':         'n',
        'UNITS':        '[DNA, RNA]',
        'VAL':          'RNA'
    }

}


###############################################################################
# Config functions (use these to change the current parameters)
###############################################################################

def update_param(param, value):
    _CALCULATION_PARAMETERS[param]['VAL'] = value


def update_params(update_dict):
    for param, value in update_dict.items():
        updateParam(param, value)


def _assemble_params():
    _flag = lambda k: '-' + cp[k]['FLAG']
    _val = lambda k: str(cp[k]['VAL'])
    cp = _CALCULATION_PARAMETERS
    return [f(k) for k in cp.keys() for f in (_flag, _val)]


###############################################################################
# Main bindings
###############################################################################

_ext_re = re.compile(r'\S+?\t(\S+?)\t*\S+\n')

def _parse_ext(ext_fn):
    '''
    Parse an ext file and return a list of the probabilities that each base
    will be involved in secondary structure (as floats).

    '''
    def _to_ss_int(sprob):
        return (1-float(sprob))
    with open(ext_fn, 'rb') as ext_fd:
        ext_fd.readline() # Discard header
        ssprobs = [_to_ss_int(m[1]) for m in (re.split(r'\t', l)
                  for l in ext_fd)]
    return ssprobs


_det_re = re.compile(r'\s+(\S+)\s+=\s+(\S+)')

def _parse_det(det_fn):
    ''' 
    Parse a det file and return the thermodynamic characteristics of the 
    strongest predicted secondary structure of the input sequence.

    '''
    with open(det_fn, 'rb') as det_fd:
        hd = det_fd.readline()
        parse = re.findall(_det_re, hd)
        out = {k: v for k, v in parse}
    return out


def _write_input(seq, seq2=None):
    fp_in = unique_fn(directory=TMP_DIR, file_ext='.seq')
    with open(fp_in, 'w') as fd:
        fd.write(seq)
    return fp_in


def calc_fold(seq):
    ''' Calculate the predicted secondary structure using full partition 
    functions and return a list of base-pairing probabilities.

    '''

    fp_in = _write_input(seq)
    fp_out = fp_in.rstrip('.seq') + '.37.ext'
    args = _assemble_params() + ['--suffix DAT', '--tracebacks', '100']
    check_call(['hybrid-ss'] + args + [fp_in], cwd=TMP_DIR,
               stdout=DEV_NULL, stderr=DEV_NULL)
    return _parse_ext(fp_out)


def calc_fold_mfe(seq):
    ''' Calculate the predicted mfe secondary structure and return a list of
    base-pairing probabilities.

    This is much faster than the above version, but will also be less accurate.
    See the UNAfold documenation for more information.

    '''
    fp_in = _write_input(seq)
    fp_out = fp_in.rstrip('.seq') + '.37.ext'
    args = _assemble_params()
    check_call(['hybrid-ss-min'] + args + [fp_in], cwd=TMP_DIR,
               stdout=DEV_NULL, stderr=DEV_NULL)
    return _parse_ext(fp_out)


def calc_thermo(seq):
    ''' Calculate the predicted secondary structure using full partition 
    functions and return the thermodynamic characteristics of the strongest
    structure.

    '''
    fp_in = _write_input(seq)
    fp_out = fp_in.rstrip('.seq') + '.det'
    args = _assemble_params()
    check_call(['UNAfold.pl'] + args + [fp_in], cwd=TMP_DIR, stdout=DEV_NULL, 
               stderr=DEV_NULL)
    return _parse_det(fp_out)


def calc_and_write(record):
    out_file, seq_name, seq = record[0], record[1], record[2]
    print 'calculating ss data for:', seq_name
    probs = calc_fold(seq)
    print 'DONE calculating ss data for:', seq_name
    with open(out_file, 'a') as fd:
        fd.write('{}\t{}\n'.format(seq_name, ','.join(probs)))


###############################################################################
# Utility functions
###############################################################################

_bars = (u'\u2581', u'\u2582', u'\u2583', u'\u2584', u'\u2585', u'\u2586',
         u'\u2587', u'\u2588')

def unicode_sparkline(probs, log_scale=True):
    '''
    Takes a list of floats between 0-1 and returns a list of respective
    unicode sparkline characters
    '''
    bars_len = len(_bars) - 1
    if log_scale:
        probs = [max(log10(max(prb, 0.0001) * 10), 0.0) for prb in probs]
    sparkline = [_bars[int(prb*bars_len)].encode('utf8') for prb in probs]
    return sparkline
