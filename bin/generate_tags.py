import sys 
import argparse
import shutil
import logbook
from logbook import info, notice, warn
from generatetags import find_tags, load_genome


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('genome_dir')
    parser.add_argument('primers_fastq')
    parser.add_argument('output_dir')
    parser.add_argument('-c', '--re_cut_site', default='CCGG')
    parser.add_argument('-m', '--max-distance-rsite-primer',
                        type=int, default=100)
    args = parser.parse_args(sys.argv[1:])
    assert os.path.isdir(args.genome_dir)
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    logfile_name = os.path.join(args.output_dir, 'log.txt')
    fmt = '{record.level_name} {record.func_name}: {record.message}'
    filehandler = logbook.FileHandler(logfile_name, mode='w',
                                      format_string=fmt)
    stderrhandler = logbook.StderrHandler(bubble=True, format_string=fmt)

    with filehandler.applicationbound():
        with stderrhandler.applicationbound():
            notice('Script arguments were:')
            notice(' '.join(sys.argv))
            primer_path = os.path.join(args.output_dir, 'primers.fa')
            notice('restriction enzyme is %s. ' % args.re_cut_site)
            shutil.copy(args.primers_fastq, primer_path)
            notice('Copied primers to %s' % primer_path)
            genome = load_genome(args.genome_dir)
            find_tags(args.output_dir, genome,
                      primer_path, args.re_cut_site,
                      args.max_distance_rsite_primer)
