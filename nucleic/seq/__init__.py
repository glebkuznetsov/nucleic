'''
nucleic.seq

Modules related to the creation and manupulation of DNA and RNA sequence
strings, including compact encoding schemes.

'''

from barcode import BarcodeGenerator

import binary

from manip import reverse, complement, reverse_complement

import randseq

# Space efficient abbreviations
r = reverse
c = complement
rc = reverse_complement
