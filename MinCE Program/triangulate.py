import os
import collections
import json
import pandas as pd

class find_genome:
    def __init__(self, new_genome, threshold, filetype):
        self.ftype = filetype
        self.unknown = new_genome
        self.sketchname = self.unknown+'.sketch'
        self.threshold = int(threshold)
        self.path_2_revDicts = 'reverseDicts/'
        self.path_2_atom_reference = 'index_map.txt'
        self.current = '10'

    def read_sketch(self):
        self.sketchname = self.unknown+'.sketch'
        self.sketch = list()
        with open(self.sketchname, 'r') as f:
            for i, line in enumerate(f):
                if i > 3:
                    self.sketch.append(line)
        self.sketch.sort()
        print('Sketch size '+str(len(self.sketch)))

    # Print iterations progress
    # Fengið að láni frá https://stackoverflow.com/a/34325723
    def printProgressBar (self, iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        # Print New Line on Complete
        if iteration == total: 
            print()
        
    def dist_all(self):
        # niðurstöðurnar á formi map
        self.result_dict = collections.defaultdict(int)
        self.sub_revDict = pd.read_csv(self.path_2_revDicts+'revDict'+self.current+'.txt', header = None, sep=',')
        self.printProgressBar(0, 1000, prefix = 'Progress:', suffix = 'Complete', length = 50)
        for i, hash_value in enumerate(self.sketch):
            # leitar í öfugri orðabók og færir inn niðurstöður
            if hash_value[0:2] != self.current:
                self.current = hash_value[0:2]
                self.sub_revDict = pd.read_csv(self.path_2_revDicts+'revDict'+self.current+'.txt', header = None, sep=',')
            self.log_sharedKmers(hash_value)
            self.printProgressBar(i + 1, 1000, prefix = 'Progress:', suffix = 'Complete', length = 50)
        return(0)

    def log_sharedKmers(self, val):
        # finnur staðinn með gildinu
        try:
            found = list(self.sub_revDict[0]).index(int(val))
        except ValueError:
            return
        found = list(self.sub_revDict[0]).index(int(val))
        # skilar öllu sem gildið bendir á (öllum rissum á formi index-map talna)
        for result in self.sub_revDict[1][found].split(';'):
            self.result_dict[result] += 1

    def read_atom_reference(self):
        # les inn index-map (vísar á nafn erfðamengis og atóm)
        # formið er [ <tilsvarandi index>, <nafn erfðamengis>, <tilsvarandi atóm>]
        self.atom_reference = list()
        with open(self.path_2_atom_reference, 'r') as f:
            for line in f:
                info = line.split(' ')
                info[-1] = info[-1].rstrip()
                # skrifar út á formið [<nafn erfðamengis>, <tilsvarandi atóm / NULL>]
                self.atom_reference.append([info[1], info[2]])

    def discern_atom(self):
        # heldur utan um hvort atóm sé komið inn
        finalist_atoms = set()
        # niðurstöður á forminu [ <nafn erfðamengis>, <stig úr samanburði>, <atóm>]
        self.finalists = list()
        # skoðar hvort rissur í niðurstöðum séu í atómi
        # formið er [ <tilsvarandi index>, <nafn erfðamengis>, <tilsvarandi atóm>]
        self.read_atom_reference()
        for result in self.result_dict:
            if self.result_dict[result] > self.threshold - 1:
                index_of_genome = int(result)
                if self.atom_reference[index_of_genome][1] == 'NULL':
                    # ef NULL, þá ekki í atómi
                    # formið er [ <nafn erfðamengis>, <stig úr samanburði>, <-1>]
                    self.finalists.append([self.atom_reference[index_of_genome][0], self.result_dict[result], -1])
                else:
                    if self.atom_reference[index_of_genome][1] not in finalist_atoms:
                        # annars í atómi
                        # formið er [ <nafn erfðamengis>, <stig úr samanburði>, <atóm>]
                        self.finalists.append([self.atom_reference[index_of_genome][0], self.result_dict[result], self.atom_reference[index_of_genome][1]])
                        # bætum atóminu í safnið, svo við tvítökum það ekki
                        finalist_atoms.add(self.atom_reference[index_of_genome][1])
                    
        self.results = sorted(self.finalists, key = lambda x: x[1], reverse = True)[0:min(5, len(self.finalists))]

        return(0)

    def give_results(self):
        def read_json(filename):
            with open(filename, "r") as rfile:
                data = json.load(rfile)
            return data

        print('\nResults for input genome '+self.unknown+'...\n\n')
        print('---------------------------------------------------------------------')
        for i, finalist in enumerate(self.results):
            # ef niðurstaðan kemur úr atómi...
            if finalist[2] != -1:
                size = str(len(read_json('features/'+finalist[2]+'.atom.json')['members']))
                print(str(i+1)+':\n Atom ' + finalist[2] + ' of size '+size+' with match at most ' + str(finalist[1]) + '/1000')
                if finalist[1] > 995:
                    cmd = 'python3 deBruijn/solve_features.py ' + self.unknown + ' features/' + finalist[2] + '.atom.json -f '+self.ftype
                    os.system(cmd)
                else:
                    print('\tAtom score too low to justify screening...')
            else:
                print(str(i+1)+':\n Genome ' + finalist[0] + ' with match ' + str(finalist[1]) + '/1000')
        print('---------------------------------------------------------------------')
        return(0)
