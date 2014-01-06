'''
findrt.py

Finds gene specific reverse transcription primers for FISSEQ samples

MMuL-V Conditions (for function parameters)
75 mM KCl
3 mM MgCl2
8 uM dNTP

'''

import os
import random
import re
import sys

from math import log, ceil

from nucleic.bindings import unafold, primer3
from nucleic.seq import barcode, manip
from nucleic.io import fasta


CWD = os.path.dirname(os.path.realpath(__file__))


class Primer(object):

    T2S = 'ACTTCAGCTGCCCCGGGTGAAGA'     # 5' to 3'


    def __init__(self, gene_ref, start, seq, init_p3_data, designed=False):
        self.gene = gene_ref           # Respective gene reference string
        self.designed = designed       # True if rationally designed
        self.seq = seq                 # Gene sp. primer seq (no T2s/bc)
        self.barcode = None            # Barcode not yet assigned
        self.barcode_data = None       # Thermodynamic data from barcode gen
        self.start = start             # Primer start (3' end)
        self.stop = self.start + len(self.seq)  # Primer stop (5' end)
        self.init_p3 = init_p3_data    # Initial p3 data
        self.final_p3 = None           # Final p3 data (includes T2S + bc)
        self.full_seq = None           # Full sequence (T2S + barcode + seq)
        self.score = None              # Final weighted score


    def _calcFullP3(self):
        self.final_p3 = primer3.assess_oligo(Primer.T2S + self.seq)
        p3_assess_score = ((self.final_p3[0][3] or 0) +
                          (self.final_p3[1][3] or 0))
        return p3_assess_score


    def _assessLocalSS(self, ss_probs):
        local_start = self.start - 3 # note primer3 has a zero-based idx
        local_end = self.start + len(self.seq) + 3
        self.localSS = (sum(ss_probs[local_start:local_end]) /
                       (local_end - local_start))
        return self.localSS


    def calcWeightedScore(self, ss_probs):
        p3_score = self._calcFullP3()              # Tm, 0-100+
        SS_score = self._assessLocalSS(ss_probs)    # 0 to 1
        p3_rank = self.init_p3.get('RANK', 1)      # 1 to 100
        end_score = max(((self.start+len(self.seq))-200), 0)
        self.score = (p3_score / 100.0 * 2 +
                      SS_score * 5 +
                      p3_rank / 100.0 +
                      end_score * 0.008)
        return self.score


    def genBarcode(self, barcode_gen, random=True):
        if random:
            self.barcode = barcode_gen.random()
        else:
            bc_assess = barcode_gen.find(Primer.T2S, self.seq)
            self.barcode = bc_assess[0]
        self.full_seq = Primer.T2S + self.barcode + self.seq


class Transcript(object):

    def __init__(self, gene_ref, gene_seq, barcode_gen, ssprobs=None,
                 exclude_seqs=None):
        self.gene = gene_ref
        self.seq = gene_seq
        self.bc_gen = barcode_gen
        self.fiveprimespace = 25
        self.tile_spacing = 3
        self.exclude_seqs = exclude_seqs + \
                            [manip.reverse_complement(s) for s in exclude_seqs]
        self._calcSS(ssprobs)


    def genTop10(self):
        pl = self._genPrimerList()
        self._processPrimerList(pl)
        self.primers = self._analyzePrimerList(self.primers, overlap=False)[:10]


    def genTiledPrimers(self):
        st = self.fiveprimespace   # bp at which to start tiling
        length = len(self.seq)
        self.tiled_primers = []
        while st < (length-15):
            self.tiled_primers.append(self._singlePrimer(st, manip.complement(self.seq[st:st+15][::-1])))
            st += 15
        self.tiled_primers = self._analyzePrimerList(self.tiled_primers,
                                                     sort_list=False)


    def _calcSS(self, ssprobs=None):
        '''
        Calculate the probability of secondary structure at each base along the
        gene sequence
        '''
        self.ss_probs = ssprobs or unafold.calc_fold(self.seq)


    def _singlePrimer(self, start, seq):
        p3_args = {
            'SEQUENCE_ID': 'testEX',
            'SEQUENCE_TEMPLATE': self.seq[self.fiveprimespace:],
            'SEQUENCE_PRIMER_REVCOMP': seq,
            'PRIMER_DNTP_CONC': '0.08', # mM in 1X MMuL-V buffer
            'PRIMER_SALT_MONOVALENT': '75', # mM in 1X MMuL-V buffer
            'PRIMER_SALT_DIVALENT': '3', # mM in 1X MMuL-V buffer
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': '0',
            'PRIMER_PICK_RIGHT_PRIMER': '1',
            'PRIMER_PICK_INTERNAL_PRIMER': '0',
            'PRIMER_MAX_SIZE': '16',
            'PRIMER_MIN_SIZE': '12',
            'PRIMER_OPT_SIZE': '14',
            'PRIMER_EXPLAIN_FLAG': '1',
            'PRIMER_MIN_TM': '0',
            'PRIMER_OPT_TM': '47',
            'PRIMER_MAX_TM': '100',
            'PRIMER_MIN_GC': '0',
            'PRIMER_MAX_GC': '100',
            'PRIMER_MAX_POLY_X': '100',
            'PRIMER_MAX_HAIRPIN_TH': '200',
            'PRIMER_MAX_COML_ANY': '100',
            'PRIMER_NUM_RETURN': '1'
        }
        primer_list = primer3.run_p3_main(p3_args)
        pf = 'PRIMER_RIGHT_0'
        try:
            p3_data = {
                'RANK': 1,
                'START_BP': int(primer_list[pf].split(',')[0]) + self.fiveprimespace,
                'SEQ': primer_list[pf + '_SEQUENCE'],
                'TM': float(primer_list[pf + '_TM']),
                'GC_PERCENT': float(primer_list[pf + '_GC_PERCENT']),
                'SELF_ANY_TH': float(primer_list[pf + '_SELF_ANY_TH']),
                'SELF_END_TH': float(primer_list[pf + '_SELF_END_TH']),
                'HARIPIN_TH': float(primer_list[pf + '_HAIRPIN_TH']),
                'END_STABILITY': float(primer_list[pf + '_END_STABILITY']),
                'PENALTY': float(primer_list[pf + '_PENALTY'])
            }
            return Primer(
                gene_ref = self.gene,
                start = p3_data['START_BP'],
                seq = p3_data['SEQ'],
                init_p3_data = p3_data)
        except:
            return Primer(
                gene_ref = self.gene,
                start = start,
                seq = seq,
                init_p3_data = {})


    def _genPrimerList(self):
        p3_args = {
            'SEQUENCE_ID': 'testEX',
            'SEQUENCE_TEMPLATE': self.seq[self.fiveprimespace:],
            'PRIMER_DNTP_CONC': '0.08', # mM in 1X MMuL-V buffer
            'PRIMER_SALT_MONOVALENT': '75', # mM in 1X MMuL-V buffer
            'PRIMER_SALT_DIVALENT': '3', # mM in 1X MMuL-V buffer
            'PRIMER_TASK': 'pick_primer_list',
            'PRIMER_PICK_LEFT_PRIMER': '0',
            'PRIMER_PICK_RIGHT_PRIMER': '1',
            'PRIMER_PICK_INTERNAL_PRIMER': '0',
            'PRIMER_MAX_SIZE': '16',
            'PRIMER_MIN_SIZE': '12',
            'PRIMER_OPT_SIZE': '14',
            'PRIMER_EXPLAIN_FLAG': '1',
            'PRIMER_MIN_TM': '43',
            'PRIMER_OPT_TM': '47',
            'PRIMER_MAX_TM': '51',
            'PRIMER_NUM_RETURN': '100'
        }
        return primer3.run_p3_main(p3_args)


    def _processPrimerList(self, primer_list):
        self.primers = []
        for x in xrange(100):
            if 'PRIMER_RIGHT_' + str(x) in primer_list.keys():
                pf = 'PRIMER_RIGHT_' + str(x)
                p3_data = {
                    'RANK': x,
                    'START_BP': int(primer_list[pf].split(',')[0]) + self.fiveprimespace,
                    'SEQ': primer_list[pf + '_SEQUENCE'],
                    'TM': float(primer_list[pf + '_TM']),
                    'GC_PERCENT': float(primer_list[pf + '_GC_PERCENT']),
                    'SELF_ANY_TH': float(primer_list[pf + '_SELF_ANY_TH']),
                    'SELF_END_TH': float(primer_list[pf + '_SELF_END_TH']),
                    'HARIPIN_TH': float(primer_list[pf + '_HAIRPIN_TH']),
                    'END_STABILITY': float(primer_list[pf + '_END_STABILITY']),
                    'PENALTY': float(primer_list[pf + '_PENALTY'])
                }
                self.primers.append(Primer(
                    gene_ref = self.gene,
                    start = p3_data['START_BP'],
                    seq = p3_data['SEQ'],
                    init_p3_data = p3_data,
                    designed = True)
                )


    def _analyzePrimerList(self, primer_list, sort_list=True, return_bcs=True,
                           random_bcs=True, overlap=True):
        def _exclusionFilter(primer):
            for s in self.exclude_seqs:
                if s in primer.full_seq:
                    return False
            return True
        for primer in primer_list:
            primer.genBarcode(self.bc_gen, random=random_bcs)
            primer.calcWeightedScore(self.ss_probs)
        if self.exclude_seqs:
            primer_list = filter(_exclusionFilter, primer_list)
        if sort_list:
            primer_list.sort(key=lambda p: p.score)
        if return_bcs:
            for primer in primer_list[10:]:
                self.bc_gen.return_bc(primer.barcode)
        if not overlap:
            pl = [primer_list[0]]
            i = 1
            while len(pl) < 10:
                start = primer_list[i-1].start
                stop = len(primer_list[i-1].seq) + start
                if primer_list[i].start not in range(start-3, stop+3):
                    pl.append(primer_list[i])
                i += 1
            return pl
        return primer_list


def estNumBarcodes(to_process):
    bc_num = 0
    for record in to_process:
        bc_num += (len(record.seq) - 25) / 18
        bc_num += 10
    bc_len = int(ceil(log(bc_num, 4)))
    return bc_num, bc_len


def getSSprobs():
    ssprobs = {}
    for fn in os.listdir(os.path.join(CWD, 'ssprobs')):
        fp = os.path.join(CWD, 'ssprobs', fn)
        ssprobs.update({r.title: [int(float(v)) for v in r.seq.split(',')]
                        for r in fasta.parseSSprobs(fp)})
    return ssprobs


def checkSSprobs(to_process, ssprobs):
    for (rec_title, rec_seq) in to_process:
        if not rec_title in ssprobs.keys():
            raise IOError('No ssprobs available for: {}.\n Run calcss.py prior'
                          ' to running findrt.py' .format(rec_title))


def process_file(fn, outfolder=None, exclude_seqs=None):
    to_process = []
    outfolder = outfolder or os.path.join(CWD, 'output')
    tiled_out_fn = os.path.join(outfolder, 'tiled' + os.path.basename(fn))
    designed_out_fn = os.path.join(outfolder, 'designed' + os.path.basename(fn))
    to_process = fasta.parse_fasta(fn)
    print '\nInput file: {}, length: {}, seq lengths: {}'.format(fn, len(to_process),
              [len(rec.seq) for rec in to_process])
    ssprobs = getSSprobs()
    checkSSprobs(to_process, ssprobs)
    bc_num, bc_len = estNumBarcodes(to_process)
    print 'Required num of barcodes:', bc_num
    print 'Required barcode length:', bc_len, '\n'
    bf = barcode.BarcodeGenerator(length=bc_len, num_needed=bc_num)
    tiled_primers = []
    designed_primers = []
    for rec in to_process:
        print 'Finding primers for:', rec.title
        t = Transcript(rec.title, rec.seq, bf, ssprobs=ssprobs.get(rec.title), exclude_seqs=exclude_seqs)
        t.genTop10()
        t.genTiledPrimers()
        designed_primers.extend(t.primers[:10])
        tiled_primers.extend(t.tiled_primers)
    with open(tiled_out_fn, 'wb') as tiled_out_fd, \
         open(designed_out_fn, 'wb') as designed_out_fd:
        print 'Writing tiled output to:', tiled_out_fn
        writeOutput(tiled_out_fd, tiled_primers)
        print 'Writing designed output to:', designed_out_fn
        writeOutput(designed_out_fd, designed_primers)
        print 'Complete.'
        print '*' * 40 + '\n'


def writeOutput(out_fd, primer_list):
    headers = '\t'.join(['gene_ref', 'full_primer', 't2s', 'barcode', 'gs_primer',
                       'start_bp', 'designed', 'weighted_score', 'gs_p3_tm',
                       'gs_p3_gc','gs_p3_self_any_th', 'gs_p3_self_end_th',
                       'gs_p3_harpin_th', 'gs_p3_end_stability',
                       'gs_p3_penalty', 'full_p3_homodimer_dg',
                       'full_p3_homodimer_ds', 'full_p3_homodimer_dh',
                       'full_ps_homodimer_tm', 'full_p3_hairpin_dg',
                       'full_p3_hairpin_ds', 'full_p3_hairpin_dh',
                       'full_ps_hairpin_tm'])
    out_fd.write(headers + '\n')
    for p in primer_list:
        ip3d = p.init_p3
        fp3d = p.final_p3
        primer_data = '\t'.join(str(x) for x in [p.gene, p.full_seq,
                               Primer.T2S, p.barcode, p.seq, p.start,
                               p.designed, p.score, ip3d.get('TM'),
                               ip3d.get('GC_PERCENT'), ip3d.get('SELF_ANY_TH'),
                               ip3d.get('SELF_END_TH'), ip3d.get('HAIRPIN_TH'),
                               ip3d.get('END_STABILITY'), ip3d.get('PENALTY'),
                               fp3d['homo_out']['dG'], fp3d['homo_out']['dS'],
                               fp3d['homo_out']['dH'], fp3d['homo_out']['t'],
                               fp3d['hrp_out']['dG'], fp3d['hrp_out']['dS'],
                               fp3d['hrp_out']['dH'], fp3d['hrp_out']['t']])
        out_fd.write(primer_data + '\n')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Design primers for gene-specific RT')
    parser.add_argument('-o', '--outfolder', help='folder in which output files will be written')
    parser.add_argument('-e', '--exclude', nargs='*', help='subsequences to exclude (e.g., restriction sites)')
    parser.add_argument('infiles', nargs='*', help='input files')
    args = parser.parse_args()
    print args
    for fn in args.infiles:
        print fn
        process_file(fn, outfolder=args.outfolder, exclude_seqs=args.exclude)
