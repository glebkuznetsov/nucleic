'''
barcode.py

'''
import itertools
import random

from nucleic.bindings.primer3 import assess_oligo


class BarcodeGenerator(object):
    ''' Generates unique 4-12bp barcode sequences randomly or psuedo-rationally.
    '''

    def __init__(self, length=5, num_needed=None):
        if not 3 < length < 13:
            raise ValueError('Barcode length must be between 4-12bp')
        self.bc_list = [''.join(b) for b in itertools.product('ATGC', 
                                                              repeat=length)]
        self._filter_complexity(length, (num_needed or 4**5))

    def return_bc(self, bc):
        ''' Returns an unused barcode to the pool '''
        self.bc_list.append(bc)

    def random(self):
        ''' Return a random barcode from the list of remaining barcodes '''
        return self.bc_list.pop(random.randint(0, len(self.bc_list)-1))

    def find(self, pre_seq='', post_seq='', max_ss_tm=37):
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
            assess = assess_oligo(pre_seq + barcode + post_seq)
            if max(assess[0][3], assess[1][3]) < max_ss_tm:
                bc_assess = (barcode, assess)
                break
            else:
                at_bc_data.append((barcode, assess))
            if i == 0:
                bc_assess = sorted(at_bc_data, key=_reduceTms)[0]
                break
        self.bc_list.remove(bc_assess[0])
        return bc_assess

    def _filter_complexity(self, length, num_needed):
        bases = ['A', 'T', 'G', 'C']
        diff = num_needed - 4**length
        if diff > 4:
            [self.bc_list.remove(length * b) for b in bases]
        if diff > 68:
            [[(self.bc_list.remove((length-1) * b1 + b2),
               self.bc_list.remove(b2 + ((length-1) * b1))) for b1 in bases]
              for b2 in bases]
