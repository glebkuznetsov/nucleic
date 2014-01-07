'''
nucleic.filters | filters.py

Filter DNA sequences by attributes like Tm, GC content, homopolymer runs, etc.
All filter function return True if a seq passes the filter (or False if not)

'''


def gc_percent(seq, gc_min, gc_max):
    ''' Return True if seq has a gc percentage betweeen gc_min and gc_max
    (inclusive)

    gc_min and gc_max should be decimal percentages (i.e., .4 for 40%)
    
    '''
    gc = seq.count('G') + seq.count('C')
    return gc_min <= gc/len(seq) <= gc_max


def gc_run(seq, run_length):
    ''' Return True of seq has a maximum GC run length <= run_length 
    '''
    lrun = 0
    for b in seq:
        if b in 'GC':
            lrun += 1
            if lrun > run_length:
                return False
        else:
            lrun = 0
    return True


def at_run(seq, run_length):
    ''' Return True of seq has a maximum AT run length <= run_length 
    '''
    lrun = 0
    for b in seq:
        if b in 'AT':
            lrun += 1
            if lrun > run_length:
                return False
        else:
            lrun = 0
    return True


def homopol_run(seq, run_length):
    ''' Return True of seq has a maximum homopolymer run length <= run_length 
    '''
    prev = ''
    lrun = 0
    for b in seq:
        if b == prev:
            lrun += 1
            if lrun > run_length:
                return False
        else:
            lrun = 0
    return True


def max_run(seq, max_a=None, max_t=None, max_g=None, max_c=None, 
                 max_at=None, max_gc=None):
    ''' Return True of seq has maximum A, T, G, C, AT, or GC run <= a provided
    maximum.
    '''
    prev = ''
    lrun = 0
    gcrun = 0
    atrun = 0
    for b in seq:
        if b == prev:
            lrun += 1
            if max_a and b == 'A' and lrun > max_a:
                return False        
            if max_t and b == 'T' and lrun > max_t:
                return False  
            if max_g and b == 'G' and lrun > max_g:
                return False  
            if max_c and b == 'C' and lrun > max_c:
                return False          
        else:
            lrun = 0
        if b in 'GC':
            gcrun += 1
            if max_gc and gcrun > max_gc:
                return False
            atrun = 0
        elif b in 'AT':
            atrun += 1
            if max_at and atrun > max_at:
                return False
            gcrun = 0
    return True


def three_prime(seq, max_gcs, max_gc_run):
    tp = seq[-5:]
    gcs = tp.count('G') + tp.count('C')
    if gcs > max_gcs:
        return False
    gc_run = 0
    for b in tp:
        if b in 'GC':
            gc_run += 1
            if gc_run > max_gc_run:
                return False
        else:
            gc_run = 0
    return True 