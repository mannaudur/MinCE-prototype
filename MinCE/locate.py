"""
Locate

Usage:
    locate <input_file> [-M <M>] [-d <d>] [-c <c>] [-t <t>] [-f <f>]

Options:
    -h --help   Show this screen.
    -M=<M>      Megasketch for big fastq files (1 to activate) [default: 0]
    -t=<t>      Threshold t/1000 to show result [default: 990]
    -c=<c>      Candidate set limit [default: 1]
    -f=<f>      Filetype (fasta/fastq) [default: fasta]
"""

from docopt import docopt
import os
import subprocess
from triangulate import find_genome

def do_sketch(path, cand_limit, mega):
        terminal_command = path+' -c '+cand_limit
        if mega:
            print('Creating Megasketch...')
            terminal_command = 'smash/bin/sketch -Mt '+terminal_command
            # ef beðið er um Megasketch, keyra nýja reikniritið:
            #subprocess.run(terminal_command, shell=True, capture_output=True)
            os.system(terminal_command)
            print('Megasketch successful')
        else:
            print('Creating sketch...')
            terminal_command = 'smash/bin/sketch -t '+terminal_command
            # annars það gamla...
            #subprocess.run(terminal_command, shell=True, capture_output=True)
            os.system(terminal_command)
            print('Sketch successful')
    
        
def main(args):
    new_genome = args['<input_file>']
    mega = int(args['-M'])
    cand_limit = args['-c']
    assert int(cand_limit) > 0, "Candidate limit must be positive"
    threshold = args['-t']
    assert int(threshold) <= 1000, "Threshold must be within [0, 1000]"
    filetype = args['-f']
    if filetype not in set(['fastq', 'fasta']):
        if new_genome[-1] == 'q':
            filetype = 'fastq'
        elif new_genome[-1] == 'a':
            filetype = 'fasta'
        else:
            print('Please specify filetype -f')
            return(-1)
    print('Starting process...')
    do_sketch(new_genome, cand_limit, mega)
    finder = find_genome(new_genome, threshold, filetype)
    finder.read_sketch()
    print('Triangulating sketch...')
    finder.dist_all()
    print('Checking for atoms...')
    finder.discern_atom()
    finder.give_results()

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)