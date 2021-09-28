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
    
    os.system('mkdir revDicts')

    rev_dict = []
    for i in range(len(data[0])):
        kash = list(data[0][i].split(' '))[:-1]
        rev_dict.append(",".join(kash[0],";".join(kash[1::])))
    rev_dict.sort()
    with open('MAX_HASH.txt') as max_hash:
        max_hash.write(rev_dict[-1].split(",")[0])
        max_hash.close()

    for i in range(10,100):
        dictName = 'revDicts/revDict'+str(i)+'.txt'
        with open(dictName, 'w') as f:
            for item in rev_dict:
                if int(item[0:2]) < i:
                    print('Something fishy going on...')
                    print(int(item[0:2]))
                    print(i)
                if int(item[0:2]) > i:
                    f.close()
                    break
                else:
                    f.write("%s\n" % item)

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)