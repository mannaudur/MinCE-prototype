"""
Extract features

Usage:
    extract_features <bitmatrix-file> <fasta-file> [-k <k>] [-t <t>] [-S <S>] [-d <d>]

Options:
    -h --help   Show this screen.
    -k=<k>      Size of kmers [default: 31]
    -t=<t>      Target for number of features [default: 10]
    -S=<S>      Seed for hash [default: 42]
    -d=<d>      Destination directory of sketch [default: 0]
"""
from docopt import docopt
from deBruijn_extractor import choose_features

def main(args):
    bitmatrix = args['<bitmatrix-file>']
    fasta_file = args['<fasta-file>']
    k = int(args['-k'])
    t = int(args['-t'])
    seed = int(args['-S'])
    outdir = args['-d']
    extractor = choose_features(bitmatrix, fasta_file, k, seed, t)
    extractor.create_possibility_matrix()
    extractor.create_feature_map()
    extractor.fix_identicals()
    extractor.write_results(outdir)


if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)