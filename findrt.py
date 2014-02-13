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
from nucleic import filters
from nucleic import thermo

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
        self.linear = None
        self.pcr_f = None
        self.pcr_r = None
        self.full_seq_linear = None
        self.full_seq_pcr = None
        self.score = None              # Final weighted score


    def _calcFullP3(self):
        self.final_p3 = primer3.assess_oligo(Primer.T2S + self.seq)
        p3_assess_score = max((self.final_p3[0][3] or 0),
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
            'PRIMER_MAX_SIZE': '15',
            'PRIMER_MIN_SIZE': '15',
            'PRIMER_OPT_SIZE': '15',
            'PRIMER_EXPLAIN_FLAG': '1',
            'PRIMER_MIN_TM': '40',
            'PRIMER_OPT_TM': '43',
            'PRIMER_MAX_TM': '46',
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
    return {'fn': fn, 'designed': designed_primers,
            'tiled': tiled_primers}


def writeOutput(out_fd, primer_list):
    headers = '\t'.join(['gene_ref', 'full_primer', 't2s', 'barcode', 'gs_primer',
                        'linear_amp_primer', 'to_syntheize_linear', 'pcr_f',
                        'pcr_r', 'to_synthesize_pcr',
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
                               Primer.T2S, p.barcode, p.seq, p.linear,
                               manip.reverse_complement(p.full_seq_linear),
                               p.pcr_f, p.pcr_r, p.full_seq_pcr, p.start,
                               p.designed, p.score, ip3d.get('TM'),
                               ip3d.get('GC_PERCENT'), ip3d.get('SELF_ANY_TH'),
                               ip3d.get('SELF_END_TH'), ip3d.get('HAIRPIN_TH'),
                               ip3d.get('END_STABILITY'), ip3d.get('PENALTY'),
                               fp3d[1][2] if fp3d[1] else 'None',
                               fp3d[1][0] if fp3d[1] else 'None',
                               fp3d[1][1] if fp3d[1] else 'None',
                               fp3d[1][3] if fp3d[1] else 'None',
                               fp3d[0][2] if fp3d[0] else 'None',
                               fp3d[0][0] if fp3d[0] else 'None',
                               fp3d[0][1] if fp3d[0] else 'None',
                               fp3d[0][3] if fp3d[0] else 'None'])
        out_fd.write(primer_data + '\n')


def find_linear_primer(primer_object_list, full_list, full_list_rc):
    def _gen_primer():
        while True:
            p = ''.join([random.choice('ATGC') for x in range(13)]) + 'GAGTC' + \
                ''.join([random.choice('ATGC') for x in range(4)])
            if 'GAGTC' in p[0:15]:
                continue
            if filters.max_run(p, max_a=4, max_t=4, max_g=3, max_c=3,
                               max_at=5, max_gc=4):
                if filters.three_prime(p, max_gcs=3, max_gc_run=2, max_homopolymer=2):
                    if 48 < thermo.calc_tm(p, conc_nm=1000.0, monovalent=20.0,
                                           divalent=2.0, dntp=0.4) < 52:
                        scores = primer3.assess_oligo(p)
                        if not scores[0] or scores[0][3] < 40:
                            if not scores[1] or scores[1][3] < 40:
                                return p
    while True:
        failed = False
        lp = _gen_primer()
        for o in full_list:
            hetero = primer3.calc_heterodimer(lp, o)
            if hetero and hetero[3] > 40:
                failed = True
                break
        if failed:
            continue
        for o in full_list_rc:
            hetro = primer3.calc_heterodimer(lp, o)
            if hetero and hetero[3] > 40:
                failed = True
                break
        if not failed:
            break
    print lp
    for p in primer_object_list:
        p.linear = lp
        p.full_seq_linear = lp + p.full_seq


def find_pcr_primers(primer_object_list, full_list, full_list_rc):
    def _gen_f_primer():
        while True:
            p = ''.join([random.choice('ATGC') for x in range(6)]) + 'GAGGAG' + \
                ''.join([random.choice('ATGC') for x in range(10)])
            if 'GAGGAG' in p[0:11] or 'GAGGAG' in p[9:]:
                continue
            if filters.max_run(p, max_a=5, max_t=5, max_g=3, max_c=3,
                               max_at=5, max_gc=4):
                if filters.three_prime(p, max_gcs=3, max_gc_run=2, max_homopolymer=2):
                    if 58 < thermo.calc_tm(p, conc_nm=500.0, monovalent=50.0,
                                           divalent=1.2, dntp=0.8) < 60:
                        scores = primer3.assess_oligo(p)
                        if not scores[0] or scores[0][3] < 40:
                            if not scores[1] or scores[1][3] < 40:
                                return p
    def _gen_r_primer():
        while True:
            p = ''.join([random.choice('ATGC') for x in range(12)]) + 'CTCTTC' + \
                ''.join([random.choice('ATGC') for x in range(4)])
            if 'CTCTTC' in p[0:17] or 'CTCTTC' in p[14:]:
                continue
            if filters.max_run(p, max_a=5, max_t=5, max_g=3, max_c=3,
                               max_at=5, max_gc=4):
                if filters.three_prime(p, max_gcs=3, max_gc_run=2, max_homopolymer=2):
                    if 58 < thermo.calc_tm(p, conc_nm=500.0, monovalent=50.0,
                                           divalent=1.2, dntp=0.8) < 60:
                        scores = primer3.assess_oligo(p)
                        if not scores[0] or scores[0][3] < 40:
                            if not scores[1] or scores[1][3] < 40:
                                return p
    while True:
        failed = False
        lp = _gen_f_primer()
        for o in full_list:
            hetero = primer3.calc_heterodimer(lp, o)
            if hetero and hetero[3] > 40:
                failed = True
                break
        if failed:
            continue
        for o in full_list_rc:
            hetero = primer3.calc_heterodimer(lp, o)
            if hetero and hetero[3] > 40:
                failed = True
                break
        if not failed:
            break
    for p in primer_object_list:
        p.pcr_f = lp
        p.full_seq_pcr = lp + p.full_seq
    while True:
        failed = False
        rp = _gen_r_primer()
        hetero = primer3.calc_heterodimer(lp, rp)
        if hetero and hetero[3] > 40:
            continue
        for o in full_list:
            hetero = primer3.calc_heterodimer(rp, o)
            if hetero and hetero[3] > 40:
                failed = True
                break
        if failed:
            continue
        for o in full_list_rc:
            hetero = primer3.calc_heterodimer(rp, o)
            if hetero and hetero[3] > 40:
                failed = True
                break
        if not failed:
            break
    for p in primer_object_list:
        p.pcr_r = rp
        p.full_seq_pcr = p.full_seq_pcr + manip.reverse_complement(rp)
    print 'lp', lp, 'rp', rp


def main(infiles, outfolder, exclude_seqs):
    processed = []
    for fn in infiles:
        processed.append(process_file(fn, outfolder=outfolder,
                                      exclude_seqs=exclude_seqs))
    # Find linear amplification oligos
    # SDA takes place in 1X thermopol buffer + 0.4 mM dNTP
    # 10 mM NH4+, 10 mM K+, 2 mM Mg2+ (optimized)
    primer3.update_params({
        'DNA_CONC': 1000, # As per He and Jiang, 2013
        'DNTP_CONC': 0.4,
        'MONOVALENT_CONC': 20,
        'DIVALENT_CONC': 2
        })
    for g in processed:
        for t in ['designed', 'tiled']:
            full_list = []
            full_list += [p.full_seq_linear or p.full_seq for p in g['designed']]
            full_list += [p.full_seq_linear or p.full_seq for p in g['tiled']]
            full_list_rc = [manip.reverse_complement(p) for p in full_list]
            find_linear_primer(g[t], full_list, full_list_rc)
    # Find pcr amplification oligos
    # General PCR conditions: 50 mM K+, 1.5mM Mg2+, 0.5 uM primer, 800 uM dNTP
    primer3.update_params({
        'DNA_CONC': 500, # As per He and Jiang, 2013
        'DNTP_CONC': 0.8,
        'MONOVALENT_CONC': 50,
        'DIVALENT_CONC': 1.2
        })
    for g in processed:
        for t in ['designed', 'tiled']:
            full_list = []
            full_list += [p.full_seq_pcr or p.full_seq for p in g['designed']]
            full_list += [p.full_seq_pcr or p.full_seq for p in g['tiled']]
            full_list_rc = [manip.reverse_complement(p) for p in full_list]
            find_pcr_primers(g[t], full_list, full_list_rc)
    # Write output
    outfolder = outfolder or os.path.join(CWD, 'output')
    for g in processed:
        fn = g['fn']
        tiled_out_fn = os.path.join(outfolder, 'tiled' + os.path.basename(fn))
        designed_out_fn = os.path.join(outfolder, 'designed' + os.path.basename(fn))
        with open(tiled_out_fn, 'wb') as tiled_out_fd, \
              open(designed_out_fn, 'wb') as designed_out_fd:
            writeOutput(tiled_out_fd, g['tiled'])
            writeOutput(designed_out_fd, g['designed'])

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Design primers for gene-specific RT')
    parser.add_argument('-o', '--outfolder', help='folder in which output files will be written')
    parser.add_argument('-e', '--exclude', nargs='*', help='subsequences to exclude (e.g., restriction sites)')
    parser.add_argument('infiles', nargs='*', help='input files')
    args = parser.parse_args()
    main(args.infiles, args.outfolder, args.exclude)
