import mmh3
import json
import collections
import numpy as np
from Bio.Seq import Seq

class solve_for_features:
    def __init__(self, input_file, atom_file, input_type, k, fasta_threshold, miss_threshold, miss_ratio):
        assert input_type in set(['fastq', 'fasta']), 'Filetype must be \'fasta\' or \'fastq\''
        self.k = k
        self.strikes = miss_threshold - 1
        self.miss_ratio = miss_ratio
        self.threshold = fasta_threshold
        self.input_type = input_type
        self.input_file = input_file
        self.quarks = self.read_json(atom_file)
        self.members = set(self.quarks['members'])
        self.hashmers = set()
        self.pruned = set() # Þeir einstaklingar sem detta út
        self.results = collections.defaultdict(list)
        
    def read_json(self,filename):
        with open(filename, "r") as rfile:
            data = json.load(rfile)
        return data

    def parse_fasta(self,fn):
        desc = []
        seqs = []

        s = []
        for line in open(fn, 'r'):
            if line[0] == '>':
                desc.append(line.strip()[1:])
                if s:
                    seqs.append(''.join(s))
                s = []
            else:
                s.append(line.strip())

        seqs.append(''.join(s))
        
        return {'file': fn, 'desc': desc, 'seqs': seqs}

    def parse_fastq(self,fq):
        seqs = list()
        read = False

        for line in open(fq, 'r'):
            if read:
                seqs.append(line.strip())
                read = False
            if line[0] == '@':
                read = True
        
        return {'seqs': seqs}


    def process_sequences(self):
        self.feature_hashes = collections.defaultdict(list)
        self.features = self.quarks['feature map']
        for feature in self.features:
            self.feature_hashes[int(feature)] = self.features[feature]


    def kmer_hash(self,kmer, seed):
        rmer = str(Seq(kmer).reverse_complement())
        kmer = str(kmer)
        if kmer < rmer:
            Hash = mmh3.hash128(kmer, seed)
        else:
            Hash = mmh3.hash128(rmer, seed)
        return Hash


    def process_hashmers(self):
        if self.input_type == 'fasta':
            self.process_hashmers_fasta()
        elif self.input_type == 'fastq':
            self.process_hashmers_fastq()

    def process_hashmers_fasta(self):
        reference_set = set(self.feature_hashes.keys())
        genome = self.parse_fasta(self.input_file)
        
        for seq in genome['seqs']:
            n = len(seq) - self.k + 1
            for i in range(n):
                kash = self.kmer_hash(seq[i:i + self.k], 42)
                if kash in reference_set:
                    self.hashmers.add(kash)

    def process_hashmers_fastq(self):
        hash_counter = collections.defaultdict(int)
        reference_set = set(self.feature_hashes.keys())
        genome = self.parse_fastq(self.input_file)
        
        for seq in genome['seqs']:
            n = len(seq) - self.k + 1
            for i in range(n):
                kash = self.kmer_hash(seq[i:i + self.k], 42)
                if kash in reference_set:
                    hash_counter[kash] += 1
                    if hash_counter[kash] > self.threshold:
                        self.hashmers.add(kash)


    def find_features(self):
        self.feats = list(self.feature_hashes.keys())
        self.found = np.zeros(len(self.feats))
        for found_feat in self.hashmers.intersection(set(self.feats)):
            self.found[self.feats.index(found_feat)] = 1
        

    def prune_features(self):
        self.eligible = self.members
        self.prune_coarse()
        self.prune_fine()
        

    def prune_coarse(self):
        if len(self.features) == 0:
            return
        self.load_results()
        for contender in self.results.copy():
            hits = self.results[contender][0] + self.results[contender][1]
            missed = self.results[contender][2]
            if hits/(hits - missed) < self.miss_ratio:
                self.pruned.add(contender)
                del self.results[contender]

        for feature in self.feature_hashes:
                self.feature_hashes[feature] = list(set(self.feature_hashes[feature]) - self.pruned)
            
        for feature in self.feature_hashes.copy():
            if not self.feature_hashes[feature]:
                self.feature_hashes.pop(feature)
        
        self.eligible = self.eligible - self.pruned


    def prune_fine(self):
        continue_condition = True
        while continue_condition:
            continue_condition = False
            strike_log = collections.defaultdict(int) # Two strikes and you're out
            not_found = np.where(self.found == 0)[0] # Finnum hvaða features fundust ekki
            for index in not_found: # Fyrir öll einkenni sem fundust ekki...
                if len(self.feature_hashes[self.feats[index]]) == 1: # Fannst ekki feature sem á að vera unique?
                    strike_log[self.feature_hashes[self.feats[index]][0]] += 1 # Merkjum það í orðabók undir einstaklingnum fyrir það feature
            for bad_contender in strike_log: # Skoðum alla einstaklinga sem mældust með fjarverandi uniques
                if strike_log[bad_contender] > self.strikes: # Fengu þeir tvö strike eða fleiri?
                    self.pruned.add(bad_contender) # Bætum þeim í lista til að vinna úr síðar
                    continue_condition = True

            for feature in self.feature_hashes:
                self.feature_hashes[feature] = list(set(self.feature_hashes[feature]) - self.pruned)
            
            for feature in self.feature_hashes.copy():
                if not self.feature_hashes[feature]:
                    self.feature_hashes.pop(feature)
            
            self.eligible = self.eligible - self.pruned
                
    def load_results(self):
        self.results.clear()
        if len(self.eligible) == 0:
            print('\n')
            print('-------------------------------------------')
            print('  Error, no good matches found in cluster  ')
            print('-------------------------------------------')
            print('\n')
            self.eligible = self.members
            print('Printing entire result set: \n')
            return(None)

        if len(self.features) == 0:
            print('\n')
            print('-----------------------------------------------')
            print('  Error, atom is broken (k-mer set identical)  ')
            print('-----------------------------------------------')
            print('\n')
            self.eligible = self.members
            print('Printing members of atom: \n')
            for member in self.members:
                print(member + '\n')
            return(None)

        for member in self.eligible:
            self.results[member] = [0,0,0]
            
        for feat_found in np.where(self.found == 1)[0]:
            for individual in self.feature_hashes[self.feats[feat_found]]:
                self.results[individual][min(len(self.feature_hashes[self.feats[feat_found]])-1, 1)] += 1
        
        for feat_not_found in np.where(self.found == 0)[0]:
            for individual in self.feature_hashes[self.feats[feat_not_found]]:
                if individual in self.eligible:
                    self.results[individual][2] -= 1
                
    def print_results(self, ordered_results, r):
        if len(ordered_results) == 0:
            return
        if r  == 'all':
            r = len(ordered_results)
        print('\n')
        print('NCBI ID: \t \t \t   Unique feat\t Non-unique feat\t Missing feat')
        print('--------------------------------------------------------------------------------------')
        for i,entry in enumerate(ordered_results):
            if i > r-1:
                break
            unq = ordered_results[entry][0]
            mix = ordered_results[entry][1]
            miss = ordered_results[entry][2]
            print(entry.rsplit('_',1)[0] + ':\t\t' + str(unq)+'\t\t'+str(mix)+'\t\t\t'+str(miss))
        print('\n')

    def get_results(self, r):
        self.load_results()
        ordered_results = dict(reversed(sorted(self.results.items(), key=lambda item: item[1])))
        self.print_results(ordered_results, r)
        
        return(ordered_results)
