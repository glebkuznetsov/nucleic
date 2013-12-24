'''
seq.py

'''

from string import maketrans


# Translation tables
_dna_in = 'ATGC'
_dna_out = 'TACG'
_DNA_TRANS = maketrans(_dna_in, _dna_out)

_rna_in = 'AUGC'
_rna_out = 'UACG'
_RNA_TRANS = maketrans(_rna_in, _rna_out)


def reverse(seq):
    return seq[::-1]

def complement(seq, na_type='DNA'):
    if na_type == 'DNA':
        return seq.translate(_DNA_TRANS)
    else:
        return seq.translate(_RNA_TRANS)

def reverse_complement(seq, na_type='DNA'):
    return reverse(complement(seq, na_type))


if __name__ == '__main__':
    import random, cProfile
    seq = ''.join([random.choice('ATGC') for x in range(20)])
    print 'Original:\t\t' + seq
    print 'Reverse:\t\t' + reverse(seq)
    print 'Complement:\t\t' + complement(seq)
    print 'Rev Comp:\t\t' + reverse_complement(seq)

    def test():
        seq = ''.join([random.choice('ATGC') for x in range(20)])
        for x in xrange(10000):
            reverse_complement(seq)

    cProfile.run('test()')
