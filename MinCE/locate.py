"""
Locate

Usage:
    locate <input_file> [-M <M>] [-d <d>] [-c <c>] [-t <t>]

Options:
    -h --help   Show this screen.
    -M=<M>      Megasketch for big fastq files (1 to activate) [default: 0]
    -t=<t>      Threshold t/1000 to show result [default: 998]
    -c=<c>      Candidate set limit [default: 1]
    -d=<d>      Destination directory of sketch(s) [default: 0]
"""

from docopt import docopt
from triangulate import find_genome

        
def main(args):
    new_genome = args['<input_file>']
    mega = args['-M']
    outpath = args['-d']
    cand_limit = args['-c']
    assert cand_limit > 0, "Candidate limit must be positive"
    threshold = args['-t']
    assert threshold <= 1000, "Threshold must be within [0, 1000]"
    finder = find_genome(new_genome, mega, cand_limit, outpath, threshold)
    finder.do_sketch()
    finder.dist_all()
    finder.discern_atom()
    finder.give_results()

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)