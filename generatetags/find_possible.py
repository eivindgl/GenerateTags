'''
4C specific. Extracts the 5'->3' DNA sequence between a primer and its
restriction enzyme cut site.
'''
# TODO add check that tags starts with primers.
# output renamed primer.fa (fix _fwd _rev postfix)

import Bio.SeqIO
from pybedtools import BedTool
import collections
import bisect
from bx.intervals import intersection
import argparse
import logbook
from logbook import info, notice, warn
import os
import gzip

def load_genome_dir(genome_dir):
    '''
    Arguments:
    `genome_dir`: directory containing zipped per chrom fasta files

    Returns: A dict of sequence per chromosome
    '''
    notice('called')
    names = os.listdir(genome_dir)
    dna = {}
    for x in names:
        assert x.endswith('.gz')
        name = x.split('.')[0]
        info('loading %s' % name)
        with gzip.open(os.path.join(genome_dir, x)) as f:
            stream = Bio.SeqIO.parse(f, 'fasta')
            dna[name] = stream.next().upper()
    return dna

def load_genome(path):
    '''
    Arguments:
    `genome_dir`: directory containing zipped per chrom fasta files

    Returns: A dict of sequence per chromosome
    '''
    notice('called')
    if path.endswith('.gz'):
        f = gzip.open(path, 'rb')
    else:
        f = open(path, 'r')
    stream = Bio.SeqIO.parse(f, 'fasta')
    dna = {}
    for x in stream:
        dna[x.id] = x.seq.upper()
    f.close()
    return dna

def all_matches_in_genome(needle, genome):
    '''
    Arguments:
    `needle`: sequence we are looking for
    `genome`: dict of sequence per chrom (haystack)

    Returns: An iterable of (chrom, start, end) ntuples for each match.
    '''
    BedEntry = collections.namedtuple('BedEntry', 'chrom start end')
    notice('called with needle "%s"' % needle)
    for chrom, v in genome.items():
        idx = 0
        while True:
            idx = v.seq.find(needle, idx)
            if idx == -1:
                break
            yield BedEntry(chrom, idx, idx + len(needle))
            idx += 1

def save_restriction_sites(outpath, genome, cut_site):
    notice('called')
    rsites = all_matches_in_genome(cut_site, genome)
    with open(outpath, 'w') as f:
        for rsite in rsites:
            f.write('%s\t%s\t%s\n' % (rsite.chrom, rsite.start,
                                      rsite.end))

def save_primer_sites(outpath, primers_fasta_path, genome):
    notice('called')
    primers = []

    for x in Bio.SeqIO.parse(primers_fasta_path, 'fasta'):
        ident = all_matches_in_genome(x.seq, genome)
        rcomp = all_matches_in_genome(x.seq.reverse_complement(), genome)
        primers.append((x.id, '+', ident))
        primers.append((x.id, '-', rcomp))
    with open(outpath, 'w') as f:
        for name, strand, hits in primers:
            for hit in  hits:
                f.write('%s\t%s\t%s\t%s\t.\t%s\n' % (
                    hit.chrom, hit.start, hit.end, name, strand))

class Lookup(object):
    '''
    wrapper around bxpython. Builds tree and
    filters output on chromosome.
    '''

    def __init__(self, items=None):
        self.tree = intersection.IntervalTree()
        if items is not None:
            for x in items:
                self.tree.add(x.start, x.end, x)

    def get_in_interval(self, chrom, start, end):
        return [x for x in self.tree.find(start, end)
                if x.chrom == chrom]

def find_primer_pairs_bordering_fragment(rsites_path, primer_path,
                                         max_dist):
    '''
    Finds pairs of primers flanking the same fragment.
    A fragment is the interval between to adjacent restriction enzyme
    cut sites.
    The leftmost primer site must be on the reverse strand
    and the rightmost primer site must be on the forward strand.

    Arguments:
    `rsites_path`: bed file path to restriction cut sites
    `primer_path`: bed file path to primer matches in the genome
    `max_dist`: the maximum distance between a primer and a cut site.

    Returns:
    An iterable of named tuples containing the flanking rsites and
    the two primers.
    (left_rsite left_primer right_primer right_rsite)
    '''
    notice('called')
    struct = collections.namedtuple(
        'PrimerFragment',
        'left_rsite left_primer right_primer right_rsite')
    mspi_sites = Lookup(BedTool(rsites_path))
    primers = collections.defaultdict(list)

    for p in BedTool(primer_path):
        basename = p.name.split('_')[0]
        primers[basename].append(p)

    primers = {k:sorted(v, key=lambda x: (x.chrom, x.start))
               for k, v in primers.items()}

    tot_frag = 0
    for basename, zs in primers.items():
        info('compting for %s' % basename)
        nfrag = 0
        for l, r in zip(zs[:-1], zs[1:]):
            # RNA POL II moves 3' -> 5' along the template strand
            if (not l.chrom == r.chrom or
                l.name == r.name or
                l.strand == '+' or
                r.strand =='-' or
                l.start - r.end > 10**4):
                continue
            sites = mspi_sites.get_in_interval(l.chrom,
                                               l.start - max_dist,
                                               r.end + max_dist)
            if len(sites) < 2: continue
            starts = [x.start for x in sites]
            lidx = bisect.bisect(starts, l.start)
            ridx = bisect.bisect(starts, r.start)
            if not lidx == ridx: continue
            lsite = sites[lidx-1]
            rsite = sites[ridx]
            if (lsite.start <= l.start and
                r.end <= rsite.end):
                nfrag += 1
                yield struct(lsite, l, r, rsite)
        notice('Stored %d fragments for %s' % (nfrag, basename))
        tot_frag += nfrag
    notice('Stored %d fragments in total' % tot_frag)

tagint = collections.namedtuple('taginterval',
                                'chrom start end name strand')

def bedentry_as_string(x, name=None, extra=None):
    common = '%s\t%s\t%s\t' % (x.chrom, x.start, x.end)
    if name:
        rest = name
    else:
        rest = '%s\t.\t%s' % (x.name, x.strand)
    if extra is not None:
        rest += '\t%s' % extra
    return common + rest + '\n'


def get_tag_interval(primer, re_cut, name=None):
        if primer.start < re_cut.start:
            assert primer.strand == '+'
            return tagint(primer.chrom, primer.start,
                            re_cut.end, name, '+')
        else:
            assert primer.strand == '-'
            return tagint(primer.chrom, re_cut.start,
                            primer.end, name, '-')

def save_tags(filepath, fragments, genome=None):
    def get_tag_intervals():
        for frag in fragments:
            # rev primer
            rname = frag.left_primer.name
            if rname.endswith('_fwd') or rname.endswith('_rev'):
                rname = rname[:-4]
            rname += '_rev'
            rev_tag = get_tag_interval(frag.left_primer,
                                       frag.left_rsite, name=rname)
            yield ('left_primer', rev_tag)
            # fwd primer
            fname = frag.right_primer.name
            if fname.endswith('_fwd') or fname.endswith('_rev'):
                fname = fname[:-4]
            fname += '_fwd'
            fwd_tag = get_tag_interval(frag.right_primer,
                                       frag.right_rsite, name=fname)
            yield ('right_primer', fwd_tag)

    notice('called')
    if genome is None:
        with open(filepath, 'w') as f:
            for prim_loc, x in get_tag_intervals():
                f.write(bedentry_as_string(x, extra=prim_loc))
        return
    z = collections.defaultdict(set)
    for prim_loc, x in get_tag_intervals():
        seq = genome[x.chrom][x.start:x.end]
        if prim_loc == 'left_primer':
            assert x.strand == '-'
            seq = seq.reverse_complement()
        else:
            assert x.strand == '+'
        seq = seq.seq.tostring()
        if seq in z[x.name]:
            warn('%s has multiple identical tag sequences.' % x.name)
        else:
            z[x.name].add(seq)
    with open(filepath, 'w') as f:
        for name in sorted(z):
            v = z[name]
            while len(v):
                f.write('>%s\n' % name)
                f.write('%s\n' % v.pop())

def save_fragments(filepath, fragments):
    notice('called')
    with open(filepath, 'w') as f:
        for x in fragments:
            f.write(bedentry_as_string(x.left_rsite, name='rsite'))
            f.write(bedentry_as_string(x.left_primer))
            f.write(bedentry_as_string(x.right_rsite, name='rsite'))
            f.write(bedentry_as_string(x.right_primer))


def main(out_dir, genome, primers_filepath, re_site,
         max_dist_rsite_primer):
    primer_sites_path = os.path.join(out_dir, 'primers.bed')
    re_cut_sites_path = os.path.join(out_dir, 're_cut_sites.bed')
    tags_bed_path = os.path.join(out_dir, 'tags.bed')
    tags_raw_bed_path = os.path.join(out_dir, 'tags_raw.bed')
    tags_fa_path = os.path.join(out_dir, 'tags.fa')
    if not os.path.isfile(primer_sites_path):
        save_primer_sites(primer_sites_path, primers_filepath, genome)
    else:
        notice('%s exists. using cached version' % primer_sites_path)
    if not os.path.isfile(re_cut_sites_path):
        save_restriction_sites(re_cut_sites_path, genome, re_site)
    else:
        notice('%s exists. using cached version' % re_cut_sites_path)

    fragments = list(find_primer_pairs_bordering_fragment(
        re_cut_sites_path, primer_sites_path, max_dist_rsite_primer))
    save_fragments(tags_raw_bed_path, fragments)
    save_tags(tags_bed_path, fragments)
    save_tags(tags_fa_path, fragments, genome=genome)
