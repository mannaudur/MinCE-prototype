"""
Solve features

Usage:
    solve_features <query-file> <atom-file> [-f <f>] [-t <t>] [-r <r>] [-m <m>] [-v <v>]

Options:
    -h --help   Show this screen.
    -t=<t>      Threshold for fastq files [default: 5]
    -f=<f>      Filetype of query file (fasta/fastq) [default: check]
    -r=<r>      Max number of top results to list [default: all]
    -m=<m>      Threshold of missed unique matches for member to be dropped [default: 3]
    -v=<v>      Percentage threshold for found/total for member to be dropped [default: 0.5]
"""

from docopt import docopt
from deBruijn_solver import solve_for_features


def main(args):
    query = args['<query-file>']
    atom = args['<atom-file>']
    k = 31
    t = int(args['-t'])
    ftype = args['-f']
    ratio = float(args['-v'])
    r = args['-r']
    if r != 'all':
        r = int(args['-r'])
    m = int(args['-m'])
    if ftype not in set(['fastq', 'fasta']):
        if query[-1] == 'q':
            ftype = 'fastq'
            print('...assuming filetype is fastq')
        elif query[-1] == 'a':
            ftype = 'fasta'
            print('...assuming filetype is fasta')
        else:
            print('Please specify filetype -f')
            return(-1)

    solver = solve_for_features(query, atom, ftype, k, t, m, ratio)
    solver.process_sequences()
    solver.process_hashmers()
    solver.find_features()
    solver.prune_features()
    solver.get_results(r)
    

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)