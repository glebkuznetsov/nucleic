from nucleic.seq import randseq as rs
from nucleic.bindings import primer3
from nucleic.filters import max_run
from nucleic.seq.manip import reverse_complement as rc
import matplotlib.pyplot as plt
import numpy as np

NUM_SEQS = 10
SEQ_LEN = 30
TM_CUTOFF = 25

_found_seqs = []
_tm_array = np.zeros((10,10))

def _checkAgainstFoundSeqs(s1, found_seqs, tm_cutoff=TM_CUTOFF):
    for s2 in found_seqs:
        out1 = primer3.calc_heterodimer(s1, s2)
        out2 = primer3.calc_heterodimer(rc(s1), s2)
        out3 = primer3.calc_heterodimer(rc(s2), s1)
        out4 = primer3.calc_heterodimer(rc(s2), rc(s1))
        _, _, _, Tm1 = out1 if out1 else (0, 0, 0, 0)
        _, _, _, Tm2 = out2 if out2 else (0, 0, 0, 0)
        _, _, _, Tm3 = out3 if out3 else (0, 0, 0, 0)
        _, _, _, Tm4 = out4 if out4 else (0, 0, 0, 0)
        if max(Tm1, Tm2, Tm3, Tm4) > tm_cutoff:
            return False
    return True

for x in range(NUM_SEQS):
    while True:
        randseq = rs.randSeq(SEQ_LEN, {'A':15, 'T':40, 'G':15, 'C': 30})
        if max_run(randseq, max_a=2, max_t=2, max_g=2, max_c=2, max_at=4, max_gc=2):
            hrp_out, homo_out = primer3.assess_oligo(randseq)
            hrp_tm = hrp_out[3] if hrp_out else 0
            homo_tm = homo_out[3] if homo_out else 0
            if max(hrp_tm, homo_tm) < TM_CUTOFF:
                if _checkAgainstFoundSeqs(randseq, _found_seqs):
                    _found_seqs.append(randseq)
                    break

def _genTmArray():
    for idx1, s1 in enumerate(_found_seqs):
        for idx2, s2 in enumerate(_found_seqs):
            out1 = primer3.calc_heterodimer(s1, s2)
            out2 = primer3.calc_heterodimer(rc(s1), s2)
            out3 = primer3.calc_heterodimer(rc(s2), s1)
            out4 = primer3.calc_heterodimer(rc(s2), rc(s1))
            _, _, _, Tm1 = out1 if out1 else (0, 0, 0, 0)
            _, _, _, Tm2 = out2 if out2 else (0, 0, 0, 0)
            _, _, _, Tm3 = out3 if out3 else (0, 0, 0, 0)
            _, _, _, Tm4 = out4 if out4 else (0, 0, 0, 0)
            _tm_array[idx1][idx2] = max(Tm1, Tm2, Tm3, Tm4)


def displayHeatmap(seq_list, tm_array, title=None, show=True):

    # tm_array = ((tm_array - tm_array.mean()) /
    #               (tm_array.max() - tm_array.min()))

    fig, ax = plt.subplots()
    heatmap = ax.pcolor(tm_array, cmap=plt.cm.Blues, alpha=0.8)
    cbar = plt.colorbar(heatmap)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    ax.set_yticks(range(tm_array.shape[0]), minor=False)
    ax.set_xticks(range(tm_array.shape[0]), minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # set labels + rotation
    ax.set_xticklabels(range(tm_array.shape[0]), minor=False, rotation=0)
    ax.set_yticklabels(range(tm_array.shape[0]), minor=False, rotation=-90)

    # turn off plot grid
    ax.grid(False)

    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    if title:
        fig.canvas.set_window_title(title)
    if show:
        plt.show()
    else:
        plt.savefig('30mer_padding_ortho.pdf')

_genTmArray()

print _found_seqs
with open('30mer_padding_ortho.txt', 'w') as fd:
    for seq in _found_seqs:
        fd.write(seq + '\n')

displayHeatmap(_found_seqs, _tm_array, show=False)
