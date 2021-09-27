"""
Get Reverse Dictionaries

Usage:
    get_revDicts <input_file> [-d <d>] [-s <s>]

Options:
    -h --help   Show this screen.
    -d=<d>      Number of splits [default: 50]
    -s=<s>      Size of each dictionary (overrides -d) [default: 0]
"""
from docopt import docopt
import numpy as np
import pandas as pd
import os


        
def main(args):
    hashlocator = args['<input_file>']
    data = pd.read_csv(hashlocator, header = None)
    division = int(args['-d'])
    burst = len(data[0])/division
    if int(args['-s']):
        burst = int(args['-s'])
        division = int(np.ceil(len(data[0])/burst))
    
    os.system('mkdir revDicts')

    rev_dict = []
    for i in range(len(data[0])):
        kash = list(data[0][i].split(' '))[:-1]
        rev_dict.append([int(kash[0]),";".join(kash[1::])])
    rev_dict.sort()

    for i in range(division):
        dictName = 'revDicts/revDict'+str(i)+'.txt'
        with open(dictName, 'w') as f:
            for item in rev_dict[int(np.ceil(i*burst)):min(len(rev_dict), int(np.ceil((i+1)*burst)))]:
                item[0] = str(item[0])
                f.write("%s\n" % ','.join(item))
            f.close()
                
    with open('revDicts/revDict_companion.txt', 'w') as f:
        for j in range(division):
            revDict_j = pd.read_csv('revDicts/revDict'+str(j)+'.txt', header = None, sep=',')
            f.write(str(j) + " %s\n" % list(revDict_j[0])[-1])
        f.close()

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)