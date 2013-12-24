'''
fasta.py

'''

import re

from collections import namedtuple


FASTA_RECORD = namedtuple('fasta_record', ['title', 'seq'])


_fasta_re = re.compile(r'>(.*)[\n\r]([^>]+)')

def parse_fasta(fasta_fn):
    ''' Parse a standard fasta file and return a list of FASTA_RECORDs '''
    records = []
    with open(fasta_fn, 'rb') as fd:
        for match in re.finditer(_fasta_re, fd.read()):
            records.append(FASTA_RECORD(match.group(1).strip(),
                           ''.join(match.group(2).strip().split())))
    return records

