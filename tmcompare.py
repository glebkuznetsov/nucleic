from __future__ import print_function

import random

from nucleic import thermo
from nucleic.bindings import primer3


for x in xrange(100):
    seq = ''.join([random.choice('ATGC') for x in range(25)])
    p3tm = primer3.calc_tm(seq)
    ttm = thermo.calc_tm(seq, thermo.calc_dna_params(), 2.5e-7)
    print(p3tm, ttm, p3tm-ttm)
