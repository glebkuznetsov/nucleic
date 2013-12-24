'''
binaryLUT.py

Generate lookup tables for binary compression of DNA sequences
'''

from binary import Basic
import seq
from itertools import product


# Generate lookup tables

_seq_byte = {''.join(p): bytes(chr(Basic.encode(p)[0])) for p in
             product('ATGC', repeat=4)}
print '{' + ''.join(["'{}':{},".format(k,v) for (k,v) in
                    sorted([(k,hex(ord(v))) for k,v in _seq_byte.items()], key=lambda e: e[0])]) + '}'
_byte_seq = {v: k for k, v in _seq_byte.items()}
_int_seq = {ord(k): v for k,v in _byte_seq.items()}
_LUTis = [0] * 256
for k,v in _int_seq.items():
    _LUTis[k] = v
print "['" + "','".join(_LUTis) + "']"
_compl = {k: _seq_byte[seq.complement(v)] for k, v in _byte_seq.items()}
_int_compl = {ord(k): v for k,v in _compl.items()}
_LUTic = [0] * 256
for k,v in _int_compl.items():
    _LUTic[k] = v
print bytearray(_LUTic).__repr__()
_rev = {v: _seq_byte[k[::-1]] for k, v in _seq_byte.items()}
_int_rev = {ord(k): v for k,v in _rev.items()}
_LUTir = [0] * 256
for k,v in _int_rev.items():
    _LUTir[k] = v
print bytearray(_LUTir).__repr__()
_rev_compl = {k: _seq_byte[seq.reverse_complement(v)] for k, v in _byte_seq.items()}
_int_rev_compl = {ord(k): v for k, v in _rev_compl.items()}
_LUTirc = [0] * 256
for k,v in _int_rev_compl.items():
    _LUTirc[k] = v
print bytearray(_LUTirc).__repr__()