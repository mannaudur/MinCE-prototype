import os
import collections
import pandas as pd

class find_genome:
    def __init__(self, new_genome, megasketch, candidate_limit, output_path, threshold):
        self.unknown = new_genome
        self.mega = megasketch
        self.cand_limit = candidate_limit
        self.outpath = output_path
        self.threshold = threshold
        self.path_2_revDicts = 'revDicts/'
        self.path_2_dict_reference = 'revDicts/revDict_companion.txt'
        self.path_2_atom_reference = 'index_map.txt'

    def read_sketch(self):
        self.sketch = list()
        with open(self.sketchname, 'r') as f:
            for i, line in enumerate(f):
                if i > 3:
                    self.sketch.append(int(line))
        

    def do_sketch(self):
        terminal_command = self.unknown+' -c '+self.cand_limit
        # ef beðið er um sérstakan outpath:
        if self.outpath:
            terminal_command = terminal_command + ' -d '+self.outpath
        if self.mega:
            # ef beðið er um Megasketch, keyra nýja reikniritið:
            os.system('sketch -Mt '+terminal_command)
        else:
            # annars það gamla...
            os.system('sketch -t '+terminal_command)
        # getum ýmist bara geymt sketchinn í minni eða skrifað hann niður í skrá
        self.sketchname = self.unknown+'.sketch'
        self.read_sketch()

    def read_dict_reference(self):
        self.dict_reference = list()
        with open(self.path_2_dict_reference, 'r') as f:
            for line in f:
                # tekur inn hæsta gildið fyrir hverja rissu og geymir
                self.dict_reference.append(int(line.split(' ')[-1].rstrip()))
        
    def dist_all(self):
        # les inn hæstu gildi fyrir hverja öfuga orðabók
        self.read_dict_reference()
        # niðurstöðurnar á formi map
        self.result_dict = collections.defaultdict(int)
        # search index segir hvaða skrá á að skoða
        search_index = 0
        self.sub_revDict = pd.read_csv(self.path_2_revDicts+'revDict0.txt', header = None, sep=',')
        dict_changed = False
        for hash_value in self.sketch:
            # finnur rétta öfuga orðabók
            while hash_value > self.dict_reference[search_index]:
                dict_changed = True
                search_index += 1
            if dict_changed:
                self.sub_revDict = pd.read_csv(self.path_2_revDicts+'revDict'+str(search_index)+'.txt', header = None, sep=',')
                dict_changed = False
            # leitar í réttri orðabók og færir inn niðurstöður
            self.log_sharedKmers(hash_value)

    def log_sharedKmers(self, val):
        # finnur staðinn með gildinu
        found = list(self.sub_revDict[0]).index(val)
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
                if self.atom_reference[result][2] == 'NULL':
                    # ef NULL, þá ekki í atómi
                    # formið er [ <nafn erfðamengis>, <stig úr samanburði>, <-1>]
                    self.finalists.append([self.atom_reference[result][1], self.result_dict[result]], -1)
                else:
                    if self.atom_reference[result][2] not in finalist_atoms:
                        # annars í atómi
                        # formið er [ <nafn erfðamengis>, <stig úr samanburði>, <atóm>]
                        self.finalists.append([self.atom_reference[result][1], self.result_dict[result]], self.atom_reference[result][2])
                        # bætum atóminu í safnið, svo við tvítökum það ekki
                        finalist_atoms.add(self.atom_reference[result][2])
                    
        self.finalists.sort(key = lambda x: x[1], reverse = True)

    def give_results(self):
        print('Finding best results for input genome '+self.unknown+'...')
        for i, finalist in enumerate(self.finalists):
            # ef niðurstaðan kemur úr atómi...
            if finalist[2] != -1:
                print(i+': Atom ' + finalist[2] + ' with match at most ' + finalist[1] + '/1000')
                print('Individual results:')
                os.system('python3 deBruijn/solve_features.py ' + self.unknown + ' ' + finalist[2] + '.sketch.json ')
            else:
                print(i+': Genome ' + finalist[0] + ' with match ' + finalist[1] + '/1000')
