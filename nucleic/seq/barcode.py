'''
barcode.py

Generates unique 5mer barcode sequences randomly or psuedo-rationally

'''
import itertools
import random

import thermo

class BarcodeGenerator(object):

    def __init__(self, length=5, num_needed=None):
        if length > 12:
            raise ValueError('Maximum barcode length is 12 bp')
        self.bc_list = [''.join(b) for b in itertools.product('ATGC', repeat=length)]
        self._filterComplexity(length, (num_needed or 4**5))

    def returnBarcode(self, bc):
        ''' Returns an unused barcode to the pool '''
        self.bc_list.append(bc)

    def randomBarcode(self):
        ''' Grab a random barcode from the list of remaining barcodes '''
        return self.bc_list.pop(random.randint(0, len(self.bc_list)-1))

    def findBarcode(self, pre_seq='', post_seq=''):
        '''
        Find an optimal barcode with the given flanking sequences
        (pre_seq and post_seq)
        '''
        def _reduceTms(attempted_barcode):
            return (attempted_barcode[1]['hrp_out']['t'] +
                    attempted_barcode[1]['homo_out']['t'])

        pot_bcs = self.bc_list[:]
        at_bc_data = []
        i = len(pot_bcs)
        while True:
            i -= 1
            barcode = pot_bcs.pop(random.randint(0, i))
            assess = thermo.assessOligo(pre_seq + barcode + post_seq)
            if assess['pass']:
                bc_assess = (barcode, assess)
                break
            else:
                at_bc_data.append((barcode, assess))
            if i == 0:
                bc_assess = sorted(at_bc_data, key=_reduceTms)[0]
                break
        self.bc_list.remove(bc_assess[0])
        return bc_assess

    def _filterComplexity(self, length, num_needed):
        bases = ['A', 'T', 'G', 'C']
        diff = num_needed - 4**length
        if diff > 4:
            [self.bc_list.remove(length * b) for b in bases]
        if diff > 68:
            [[(self.bc_list.remove((length-1) * b1 + b2),
               self.bc_list.remove(b2 + ((length-1) * b1))) for b1 in bases]
              for b2 in bases]
