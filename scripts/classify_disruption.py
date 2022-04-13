import argparse
import pysam
import sys
import re
import gffutils
import pandas as pd
import itertools
import numpy as np
import sys
import binascii

from tqdm import tqdm
from collections import defaultdict
from functools import partial, reduce
from funcy import ilen, take, merge, drop, repeat, walk_keys

def get_length(genome):
    return sum(p['LN'] for p in genome)

def get_nonexist_int_coord(genome):
    return -(get_length(genome) + 1)

def ilen(iterable):
    return reduce(lambda sum, element: sum + 1, iterable, 0)

def hex_to_array(hexstr):
    return np.fromstring(binascii.unhexlify(hexstr), dtype=int)

def hex_to_offtargetinfo(hexstr, delim):
    mainarr = hex_to_array(hexstr)
    index = np.where(mainarr == delim)[0]
    slices = list(zip(np.insert(index, 0, -1), index)) # np.append(index[1:], [len(mainarr)]))
    out = []
    for start, end in slices:
        out += zip(repeat(mainarr[end - 1]), mainarr[start + 1:end - 1])
    return out

def get_offtargets(guide, delim):
    if not guide.has_tag('of'):
        return {
            "0 Off-Targets": 0,
            "1 Off-Targets": 0,
            "2 Off-Targets": 0,
            "3 Off-Targets": 0
        }

    ots = guide.get_tag('of')
    off_targets = hex_to_offtargetinfo(ots, delim=delim)

    return {
        "0 Off-Targets": ilen(filter(lambda x: x[0] == 0, off_targets)),
        "1 Off-Targets": ilen(filter(lambda x: x[0] == 1, off_targets)),
        "2 Off-Targets": ilen(filter(lambda x: x[0] == 2, off_targets)),
        "3 Off-Targets": ilen(filter(lambda x: x[0] == 3, off_targets))
    }
    
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parse_args():
    p = argparse.ArgumentParser(
        description="Classifies gRNA disruption events by their type"
    )

    p.add_argument(
        '--non-specific-db',
        help='gRNA database in SAM/BAM format in original genome'
    )

    p.add_argument(
        '--annotation-db',
        help='Annotation DB used to annotate sgRNAs'
    )

    p.add_argument(
        '--chromosome-mapping',
        help='Mapping file used to map between chromosome names'
    )

    p.add_argument(
        '-o', '--output',
        help='Output CSV file',
        default=sys.stdout
    )

    p.add_argument(
        'snps',
        help='VCF file containing SNPs'
    )

    p.add_argument(
        'indels',
        help='VCF file containing indels'
    )

    p.add_argument(
        'guide_db',
        help='gRNA database in SAM/BAM format'
    )

    return p.parse_args()

def get_indel_size(indels):
    size = 0
    for indel in indels:
        size += len(indel.alleles[1]) - len(indel.alleles[0])
    return size
    
def classify_disruption(indel_db, snp_db, start, end, chrm, strand):
    disruptions = {}
    keys = ['PAM SNP', 'SNP', 'PAM Insertion', 'PAM Deletion', 'Insertion', 'Deletion']
    for key in keys:
        disruptions[key] = 0

    pam_coords = (end, end + 3) if strand == '+' else (start, start + 3)

    het = False
    snps = snp_db.fetch(chrm, start, end)
    for snp in snps:
        if snp.pos >= pam_coords[0] and snp.pos <= pam_coords[1]:
            disruptions['PAM SNP'] += 1
        disruptions['SNP'] += 1

        het = next(snp.filter.iteritems())[0] == 'Het'

    indels = indel_db.fetch(chrm, start, end)
    for indel in indels:
        insertion = get_indel_size([indel]) > 0
        if indel.pos >= pam_coords[0] and indel.pos <= pam_coords[1]:
            disruptions['PAM Insertion'] += 1 if insertion else 0
            disruptions['PAM Deletion']  += 1 if not insertion else 0
        disruptions['Insertion'] += 1 if insertion else 0
        disruptions['Deletion']  += 1 if not insertion else 0

    disruptions['Het'] = het
    return disruptions

def load_non_specific_db(guides, guide_db):
    af = pysam.AlignmentFile(guide_db)
    genome = af.header['SQ']
    delim = get_nonexist_int_coord(genome)
    
    seq_to_guide = {}
    for row, idx in zip(af, itertools.count(1)):
        # debugging
        if idx % 10000 == 0:
            eprint(idx)
        if row.query_sequence in guides:
            seq_to_guide[row.query_sequence] = delim, row

    return seq_to_guide

def load_chromosome_mapping(mapping_file):
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f.readlines():
            chrm, accession = line.split()
            mapping[chrm] = accession
    return mapping

if __name__ == "__main__":
    args = parse_args()

    with pysam.AlignmentFile(args.guide_db) as guide_db:
        guides = set(rec.query_sequence for rec in guide_db)

    guide_db = pysam.AlignmentFile(args.guide_db)
    if args.non_specific_db is not None:
        seq_to_guide = load_non_specific_db(guides, args.non_specific_db)

    snps = pysam.VariantFile(args.snps)
    indels = pysam.VariantFile(args.indels)

    if args.annotation_db:
        annotation_db = gffutils.FeatureDB(args.annotation_db)

    if args.chromosome_mapping:
        chromosome_mapping = load_chromosome_mapping(args.chromosome_mapping)

    classify_disruption = partial(classify_disruption, indels, snps)

    specific_genome = guide_db.header['SQ']
    specific_delim  = get_nonexist_int_coord(specific_genome)

    rows = []
    for rec in tqdm(guide_db):
        start  = rec.reference_start 
        end    = start + 24
        chrm   = rec.reference_name
        strand = '-' if rec.is_reverse else '+'
        sgrna  = rec.query_sequence

        row = {
            'sgRNA'      : sgrna,
            'Chromosome' : chrm,
            'Start'      : start,
            'End'        : end,
            'Strand'     : strand,
        }

        opposite_ots = get_offtargets(rec, specific_delim)
        opposite_ots = walk_keys(lambda k: "Opposite Allele " + k, opposite_ots)
        row = merge(row, opposite_ots)
    
        if args.non_specific_db:
            delim, original_rec = seq_to_guide[sgrna]
            original_ots = get_offtargets(original_rec, delim)
            original_ots = walk_keys(lambda k: "Original Allele " + k, original_ots)
            row = merge(row, original_ots)

        if chrm is None:
            continue

        chrm = re.search(r'\d+|X|Y|MT', chrm)[0]

        event = classify_disruption(start, end, chrm, strand)
        if args.non_specific_db:
            _, guide = seq_to_guide[sgrna]
            #event['Cutting Efficiency'] = guide.get_tag('ds')

        if args.annotation_db:
            if args.chromosome_mapping:
                accession = chromosome_mapping[chrm]
            else:
                accession = chrm

            features = annotation_db.region(f'{accession}:{start}-{end}', featuretype='CDS')
            genes = [f.attributes['gene_name'][0] for f in features]
            row['CDS'] = genes[0] if genes else 'N/A'

        row = merge(row, event)
        rows.append(row)

    pd.DataFrame(rows).to_csv(args.output)
