import random

def weighted_choice(weights):
    choice = random.random() * sum(weights)
    for i, w in enumerate(weights):
        choice -= w
        if choice < 0:
            return i
# end def

def randSeq(length, probs=None):
    """ Probs format {'A':25, 'T':25, 'G':25, 'C':25} (25% liklihood of each)
    """
    if probs is None:
        probs = {'A':25, 'T':25, 'G':25, 'C':25}
    weights = list(probs.values())
    bases = list(probs.keys())
    return ''.join([bases[weighted_choice(weights)] for x in range(length)])
# end def
