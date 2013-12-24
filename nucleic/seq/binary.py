'''
binary.py

Methods to translate DNA sequences to integers via binary
intermediates (for more efficient memory usage)

'''

# Lookup tables

_LUTsb = {'AAAA':0x0,'AAAC':0x1,'AAAG':0x2,'AAAT':0x3,'AACA':0x4,'AACC':0x5,
'AACG':0x6,'AACT':0x7,'AAGA':0x8,'AAGC':0x9,'AAGG':0xa,'AAGT':0xb,'AATA':0xc,
'AATC':0xd,'AATG':0xe,'AATT':0xf,'ACAA':0x10,'ACAC':0x11,'ACAG':0x12,
'ACAT':0x13,'ACCA':0x14,'ACCC':0x15,'ACCG':0x16,'ACCT':0x17,'ACGA':0x18,
'ACGC':0x19,'ACGG':0x1a,'ACGT':0x1b,'ACTA':0x1c,'ACTC':0x1d,'ACTG':0x1e,
'ACTT':0x1f,'AGAA':0x20,'AGAC':0x21,'AGAG':0x22,'AGAT':0x23,'AGCA':0x24,
'AGCC':0x25,'AGCG':0x26,'AGCT':0x27,'AGGA':0x28,'AGGC':0x29,'AGGG':0x2a,
'AGGT':0x2b,'AGTA':0x2c,'AGTC':0x2d,'AGTG':0x2e,'AGTT':0x2f,'ATAA':0x30,
'ATAC':0x31,'ATAG':0x32,'ATAT':0x33,'ATCA':0x34,'ATCC':0x35,'ATCG':0x36,
'ATCT':0x37,'ATGA':0x38,'ATGC':0x39,'ATGG':0x3a,'ATGT':0x3b,'ATTA':0x3c,
'ATTC':0x3d,'ATTG':0x3e,'ATTT':0x3f,'CAAA':0x40,'CAAC':0x41,'CAAG':0x42,
'CAAT':0x43,'CACA':0x44,'CACC':0x45,'CACG':0x46,'CACT':0x47,'CAGA':0x48,
'CAGC':0x49,'CAGG':0x4a,'CAGT':0x4b,'CATA':0x4c,'CATC':0x4d,'CATG':0x4e,
'CATT':0x4f,'CCAA':0x50,'CCAC':0x51,'CCAG':0x52,'CCAT':0x53,'CCCA':0x54,
'CCCC':0x55,'CCCG':0x56,'CCCT':0x57,'CCGA':0x58,'CCGC':0x59,'CCGG':0x5a,
'CCGT':0x5b,'CCTA':0x5c,'CCTC':0x5d,'CCTG':0x5e,'CCTT':0x5f,'CGAA':0x60,
'CGAC':0x61,'CGAG':0x62,'CGAT':0x63,'CGCA':0x64,'CGCC':0x65,'CGCG':0x66,
'CGCT':0x67,'CGGA':0x68,'CGGC':0x69,'CGGG':0x6a,'CGGT':0x6b,'CGTA':0x6c,
'CGTC':0x6d,'CGTG':0x6e,'CGTT':0x6f,'CTAA':0x70,'CTAC':0x71,'CTAG':0x72,
'CTAT':0x73,'CTCA':0x74,'CTCC':0x75,'CTCG':0x76,'CTCT':0x77,'CTGA':0x78,
'CTGC':0x79,'CTGG':0x7a,'CTGT':0x7b,'CTTA':0x7c,'CTTC':0x7d,'CTTG':0x7e,
'CTTT':0x7f,'GAAA':0x80,'GAAC':0x81,'GAAG':0x82,'GAAT':0x83,'GACA':0x84,
'GACC':0x85,'GACG':0x86,'GACT':0x87,'GAGA':0x88,'GAGC':0x89,'GAGG':0x8a,
'GAGT':0x8b,'GATA':0x8c,'GATC':0x8d,'GATG':0x8e,'GATT':0x8f,'GCAA':0x90,
'GCAC':0x91,'GCAG':0x92,'GCAT':0x93,'GCCA':0x94,'GCCC':0x95,'GCCG':0x96,
'GCCT':0x97,'GCGA':0x98,'GCGC':0x99,'GCGG':0x9a,'GCGT':0x9b,'GCTA':0x9c,
'GCTC':0x9d,'GCTG':0x9e,'GCTT':0x9f,'GGAA':0xa0,'GGAC':0xa1,'GGAG':0xa2,
'GGAT':0xa3,'GGCA':0xa4,'GGCC':0xa5,'GGCG':0xa6,'GGCT':0xa7,'GGGA':0xa8,
'GGGC':0xa9,'GGGG':0xaa,'GGGT':0xab,'GGTA':0xac,'GGTC':0xad,'GGTG':0xae,
'GGTT':0xaf,'GTAA':0xb0,'GTAC':0xb1,'GTAG':0xb2,'GTAT':0xb3,'GTCA':0xb4,
'GTCC':0xb5,'GTCG':0xb6,'GTCT':0xb7,'GTGA':0xb8,'GTGC':0xb9,'GTGG':0xba,
'GTGT':0xbb,'GTTA':0xbc,'GTTC':0xbd,'GTTG':0xbe,'GTTT':0xbf,'TAAA':0xc0,
'TAAC':0xc1,'TAAG':0xc2,'TAAT':0xc3,'TACA':0xc4,'TACC':0xc5,'TACG':0xc6,
'TACT':0xc7,'TAGA':0xc8,'TAGC':0xc9,'TAGG':0xca,'TAGT':0xcb,'TATA':0xcc,
'TATC':0xcd,'TATG':0xce,'TATT':0xcf,'TCAA':0xd0,'TCAC':0xd1,'TCAG':0xd2,
'TCAT':0xd3,'TCCA':0xd4,'TCCC':0xd5,'TCCG':0xd6,'TCCT':0xd7,'TCGA':0xd8,
'TCGC':0xd9,'TCGG':0xda,'TCGT':0xdb,'TCTA':0xdc,'TCTC':0xdd,'TCTG':0xde,
'TCTT':0xdf,'TGAA':0xe0,'TGAC':0xe1,'TGAG':0xe2,'TGAT':0xe3,'TGCA':0xe4,
'TGCC':0xe5,'TGCG':0xe6,'TGCT':0xe7,'TGGA':0xe8,'TGGC':0xe9,'TGGG':0xea,
'TGGT':0xeb,'TGTA':0xec,'TGTC':0xed,'TGTG':0xee,'TGTT':0xef,'TTAA':0xf0,
'TTAC':0xf1,'TTAG':0xf2,'TTAT':0xf3,'TTCA':0xf4,'TTCC':0xf5,'TTCG':0xf6,
'TTCT':0xf7,'TTGA':0xf8,'TTGC':0xf9,'TTGG':0xfa,'TTGT':0xfb,'TTTA':0xfc,
'TTTC':0xfd,'TTTG':0xfe,'TTTT':0xff}

_LUTis = ['AAAA','AAAC','AAAG','AAAT','AACA','AACC','AACG','AACT','AAGA','AAGC',
'AAGG','AAGT','AATA','AATC','AATG','AATT','ACAA','ACAC','ACAG','ACAT','ACCA',
'ACCC','ACCG','ACCT','ACGA','ACGC','ACGG','ACGT','ACTA','ACTC','ACTG','ACTT',
'AGAA','AGAC','AGAG','AGAT','AGCA','AGCC','AGCG','AGCT','AGGA','AGGC','AGGG',
'AGGT','AGTA','AGTC','AGTG','AGTT','ATAA','ATAC','ATAG','ATAT','ATCA','ATCC',
'ATCG','ATCT','ATGA','ATGC','ATGG','ATGT','ATTA','ATTC','ATTG','ATTT','CAAA',
'CAAC','CAAG','CAAT','CACA','CACC','CACG','CACT','CAGA','CAGC','CAGG','CAGT',
'CATA','CATC','CATG','CATT','CCAA','CCAC','CCAG','CCAT','CCCA','CCCC','CCCG',
'CCCT','CCGA','CCGC','CCGG','CCGT','CCTA','CCTC','CCTG','CCTT','CGAA','CGAC',
'CGAG','CGAT','CGCA','CGCC','CGCG','CGCT','CGGA','CGGC','CGGG','CGGT','CGTA',
'CGTC','CGTG','CGTT','CTAA','CTAC','CTAG','CTAT','CTCA','CTCC','CTCG','CTCT',
'CTGA','CTGC','CTGG','CTGT','CTTA','CTTC','CTTG','CTTT','GAAA','GAAC','GAAG',
'GAAT','GACA','GACC','GACG','GACT','GAGA','GAGC','GAGG','GAGT','GATA','GATC',
'GATG','GATT','GCAA','GCAC','GCAG','GCAT','GCCA','GCCC','GCCG','GCCT','GCGA',
'GCGC','GCGG','GCGT','GCTA','GCTC','GCTG','GCTT','GGAA','GGAC','GGAG','GGAT',
'GGCA','GGCC','GGCG','GGCT','GGGA','GGGC','GGGG','GGGT','GGTA','GGTC','GGTG',
'GGTT','GTAA','GTAC','GTAG','GTAT','GTCA','GTCC','GTCG','GTCT','GTGA','GTGC',
'GTGG','GTGT','GTTA','GTTC','GTTG','GTTT','TAAA','TAAC','TAAG','TAAT','TACA',
'TACC','TACG','TACT','TAGA','TAGC','TAGG','TAGT','TATA','TATC','TATG','TATT',
'TCAA','TCAC','TCAG','TCAT','TCCA','TCCC','TCCG','TCCT','TCGA','TCGC','TCGG',
'TCGT','TCTA','TCTC','TCTG','TCTT','TGAA','TGAC','TGAG','TGAT','TGCA','TGCC',
'TGCG','TGCT','TGGA','TGGC','TGGG','TGGT','TGTA','TGTC','TGTG','TGTT','TTAA',
'TTAC','TTAG','TTAT','TTCA','TTCC','TTCG','TTCT','TTGA','TTGC','TTGG','TTGT',
'TTTA','TTTC','TTTG','TTTT']

_LUTic = bytearray(b'\xff\xfe\xfd\xfc\xfb\xfa\xf9\xf8\xf7\xf6\xf5\xf4\xf3\xf2'
'\xf1\xf0\xef\xee\xed\xec\xeb\xea\xe9\xe8\xe7\xe6\xe5\xe4\xe3\xe2\xe1\xe0\xdf'
'\xde\xdd\xdc\xdb\xda\xd9\xd8\xd7\xd6\xd5\xd4\xd3\xd2\xd1\xd0\xcf\xce\xcd\xcc'
'\xcb\xca\xc9\xc8\xc7\xc6\xc5\xc4\xc3\xc2\xc1\xc0\xbf\xbe\xbd\xbc\xbb\xba\xb9'
'\xb8\xb7\xb6\xb5\xb4\xb3\xb2\xb1\xb0\xaf\xae\xad\xac\xab\xaa\xa9\xa8\xa7\xa6'
'\xa5\xa4\xa3\xa2\xa1\xa0\x9f\x9e\x9d\x9c\x9b\x9a\x99\x98\x97\x96\x95\x94\x93'
'\x92\x91\x90\x8f\x8e\x8d\x8c\x8b\x8a\x89\x88\x87\x86\x85\x84\x83\x82\x81\x80'
'\x7f~}|{zyxwvutsrqponmlkjihgfedcba`_^]\\[ZYXWVUTSRQPONMLKJIHGFEDCBA@?>=<;:987'
'6543210/.-,+*)(\'&%$#"! \x1f\x1e\x1d\x1c\x1b\x1a\x19\x18\x17\x16\x15\x14\x13'
'\x12\x11\x10\x0f\x0e\r\x0c\x0b\n\t\x08\x07\x06\x05\x04\x03\x02\x01\x00')

_LUTir = bytearray(b'\x00@\x80\xc0\x10P\x90\xd0 `\xa0\xe00p\xb0\xf0\x04D\x84'
'\xc4\x14T\x94\xd4$d\xa4\xe44t\xb4\xf4\x08H\x88\xc8\x18X\x98\xd8(h\xa8\xe88x'
'\xb8\xf8\x0cL\x8c\xcc\x1c\\\x9c\xdc,l\xac\xec<|\xbc\xfc\x01A\x81\xc1\x11Q\x91'
'\xd1!a\xa1\xe11q\xb1\xf1\x05E\x85\xc5\x15U\x95\xd5%e\xa5\xe55u\xb5\xf5\tI\x89'
'\xc9\x19Y\x99\xd9)i\xa9\xe99y\xb9\xf9\rM\x8d\xcd\x1d]\x9d\xdd-m\xad\xed=}\xbd'
'\xfd\x02B\x82\xc2\x12R\x92\xd2"b\xa2\xe22r\xb2\xf2\x06F\x86\xc6\x16V\x96\xd6&f'
'\xa6\xe66v\xb6\xf6\nJ\x8a\xca\x1aZ\x9a\xda*j\xaa\xea:z\xba\xfa\x0eN\x8e\xce'
'\x1e^\x9e\xde.n\xae\xee>~\xbe\xfe\x03C\x83\xc3\x13S\x93\xd3#c\xa3\xe33s\xb3'
'\xf3\x07G\x87\xc7\x17W\x97\xd7\'g\xa7\xe77w\xb7\xf7\x0bK\x8b\xcb\x1b[\x9b'
'\xdb+k\xab\xeb;{\xbb\xfb\x0fO\x8f\xcf\x1f_\x9f\xdf/o\xaf\xef?\x7f\xbf\xff')


_LUTirc = bytearray(b'\xff\xbf\x7f?\xef\xafo/\xdf\x9f_\x1f\xcf\x8fO\x0f\xfb\xbb'
'{;\xeb\xabk+\xdb\x9b[\x1b\xcb\x8bK\x0b\xf7\xb7w7\xe7\xa7g\'\xd7\x97W\x17\xc7'
'\x87G\x07\xf3\xb3s3\xe3\xa3c#\xd3\x93S\x13\xc3\x83C\x03\xfe\xbe~>\xee\xaen.'
'\xde\x9e^\x1e\xce\x8eN\x0e\xfa\xbaz:\xea\xaaj*\xda\x9aZ\x1a\xca\x8aJ\n\xf6'
'\xb6v6\xe6\xa6f&\xd6\x96V\x16\xc6\x86F\x06\xf2\xb2r2\xe2\xa2b"\xd2\x92R\x12'
'\xc2\x82B\x02\xfd\xbd}=\xed\xadm-\xdd\x9d]\x1d\xcd\x8dM\r\xf9\xb9y9\xe9\xa9i)'
'\xd9\x99Y\x19\xc9\x89I\t\xf5\xb5u5\xe5\xa5e%\xd5\x95U\x15\xc5\x85E\x05\xf1'
'\xb1q1\xe1\xa1a!\xd1\x91Q\x11\xc1\x81A\x01\xfc\xbc|<\xec\xacl,\xdc\x9c\\\x1c'
'\xcc\x8cL\x0c\xf8\xb8x8\xe8\xa8h(\xd8\x98X\x18\xc8\x88H\x08\xf4\xb4t4\xe4'
'\xa4d$\xd4\x94T\x14\xc4\x84D\x04\xf0\xb0p0\xe0\xa0` \xd0\x90P\x10\xc0\x80@'
'\x00')


class Basic:
    ''' Encodes DNA sequences without a lookup table '''
    _seq_bin = {'A':'00', 'T':'11', 'U':'11', 'G':'10', 'C':'01'}
    _bin_seq_dna = {'00':'A', '11':'T', '10':'G', '01':'C'}
    _bin_seq_rna = {'00':'A', '11':'U', '10':'G', '01':'C'}

    @staticmethod
    def _seq_to_bin(seq):
        return ''.join([Basic._seq_bin[b] for b in seq])

    @staticmethod
    def _bin_to_seq(binary, na_type='DNA'):
        to_seq = Basic._bin_seq_dna if na_type == 'DNA' else Basic._bin_seq_rna
        return ''.join([to_seq[binary[x:x+2]] for x in
                        xrange(0, len(binary)-1, 2)])

    @staticmethod
    def decode(integer, length, na_type='DNA'):
        b = bin(integer)[2:]
        b_padded = '0' * (length * 2 - len(b)) + b
        return Basic._bin_to_seq(b_padded, na_type=na_type)

    @staticmethod
    def encode(seq):
        ''' Returns an integer representation of seq as well as its length '''
        return int(Basic._seq_to_bin(seq), 2), len(seq)


class LUT:
    '''
    Encodes DNA sequences with a lookup table and performs basic operations
    (complement, reverse, reverse-complement)

    '''

    @staticmethod
    def encode(seq):
        sb = _LUTsb#_seq_byte
        L = len(seq)
        pad = 4 - ((L % 4) or 4)
        seq += pad * 'A'
        return bytearray([pad] + [sb[seq[b:b+4]] for b in xrange(0, L, 4)])

    @staticmethod
    def decode(seq_ba):
        ts = _LUTis
        offset = seq_ba[0]
        front = offset & 0x04
        idx = offset & 0x03
        # if rear == 4 then offset is from front, else the rear
        if front:
            start = idx
            stop = None
        else:
            start = 0
            stop = -idx or None
        return ''.join([ts[val] for val in seq_ba[1:]])[start:stop]

    @staticmethod
    def rev(seq_ba):
        r = _LUTir
        L = len(seq_ba)
        for i in xrange(L/2-1):
            seq_ba[1+i], seq_ba[-i-1] = r[seq_ba[-i-1]], r[seq_ba[1+i]]
        if L % 2:
            seq_ba[L/2+1] = r[seq_ba[L/2+1]]
        seq_ba[0] = (seq_ba[0] ^ 0x04)
        return seq_ba

    @staticmethod
    def rev_compl(seq_ba):
        rc = _LUTirc
        L = len(seq_ba)
        for i in xrange(L/2-1):
            seq_ba[1+i], seq_ba[-i-1] = rc[seq_ba[-i-1]], rc[seq_ba[1+i]]
        if L % 2:
            seq_ba[L/2+1] = rc[seq_ba[L/2+1]]
        seq_ba[0] = (seq_ba[0] ^ 0x04)
        return seq_ba

    @staticmethod
    def compl(seq_ba):
        c = _LUTic
        for i in xrange(1, len(seq_ba)):
            seq_ba[i] = c[seq_ba[i]]
        return seq_ba


##### Test code below #####

if __name__ == '__main__':

    import random
    import cProfile
    import timeit

    def testBasic():
        global rseq_out
        for x in xrange(10000):
            a = Basic.encode(rseq)
            rseq_out = Basic.decode(*a)

    def testLUT():
        global rseq_out
        for x in xrange(10000):
            enc = LUT.encode(rseq)
            rseq_out = LUT.decode(enc)

    print '-->Timed 10000X for 25mer:'
    rseq = ''.join([random.choice('ATGC') for x in xrange(25)])
    a = timeit.Timer('testBasic()', 'from __main__ import testBasic')
    print 'Basic:', a.timeit(1)
    c = timeit.Timer('testLUT()', 'from __main__ import testLUT')
    print 'LUT:', c.timeit(1)
    # cProfile.run('testBasic()')
    # cProfile.run('testLUT()')
    # cProfile.run('testLUT()')

    print '\n--->Timed 10000X for 500mer:'
    rseq = ''.join([random.choice('ATGC') for x in xrange(500)])
    a = timeit.Timer('testBasic()', 'from __main__ import testBasic')
    print 'Basic:', a.timeit(1)
    c = timeit.Timer('testLUT()', 'from __main__ import testLUT')
    print 'LUT:', c.timeit(1)
    # cProfile.run('testBasic()')
    # cProfile.run('testLUT()')
    # cProfile.run('testLUT()')

    print '\nFull test suite with LUT:'
    rseq = ''.join([random.choice('ATGC') for x in xrange(random.randint(20,30))])
    print 'random_seq_in:\t\t\t', rseq
    enc = LUT.encode(rseq)
    print 'random_seq_out:\t\t\t', LUT.decode(enc)
    enc = LUT.encode(rseq)
    print 'random_seq_compl:\t\t', LUT.decode(LUT.compl(enc))
    enc = LUT.encode(rseq)
    print 'random_seq_rev:\t\t\t', LUT.decode(LUT.rev(enc))
    enc = LUT.encode(rseq)
    print 'random_seq_rev_compl:\t', LUT.decode(LUT.rev_compl(enc))
