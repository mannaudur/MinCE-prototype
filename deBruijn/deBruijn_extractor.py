import pandas as pd
import numpy as np
import collections
import json
import copy
from Bio.Seq import Seq
import mmh3

class choose_features:
    def __init__(self, bitmatrix_file, fasta_file, k, seed, target_number):
        self.name = fasta_file.split('/')[-1].split('.')[0]
        self.df = pd.read_csv(bitmatrix_file, sep='\t')
        self.bitmatrix = self.df.values[:,1:]
        self.fasta = self.parse_fasta(fasta_file)
        self.k = k
        self.seed = seed
        self.target = target_number
        self.members = int(np.shape(self.bitmatrix)[1])
        self.M = np.zeros((self.members,self.members))
        self.valid_features = np.where(np.sum(self.bitmatrix, axis=1)<self.members)[0]

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

    def findVector(self,vector):
        found_index = None
        for k in self.valid_features:
            if np.array_equal(self.bitmatrix[k,:], vector):
                found_index = k
                self.bitmatrix[k,:] = np.ones((self.bitmatrix[k,:].shape))
                return(found_index)
        print('Error, no valid index found for ' + np.array2string(vector))
        return(-1)

    def moveGoalPost(self, vector):
        #print('-----')
        #print('Withdrawing vector:')
        #print(vector)
        #print('from GOAL vector:')
        #print(self.GOAL)
        #print('with results:')
        self.GOAL = self.GOAL - vector
        #print(self.GOAL)
        #print('-----')

    def logConnections(self, vector):
        for individual in np.where(vector > 0)[0]:
            for connected in np.where(vector > 0)[0]:
                self.M[individual,connected] += 1

    def kmer_hash(self, kmer, seed):
        rmer = str(Seq(kmer).reverse_complement())
        kmer = str(kmer)
        if kmer < rmer:
            Hash = mmh3.hash128(kmer, seed)
        else:
            Hash = mmh3.hash128(rmer, seed)
        return(Hash)

    def appendIndiv(self, vector, index):
        for indiv in np.where(vector > 0)[0]:
            sequence = self.fasta['seqs'][index]
            kash = self.kmer_hash(sequence[0:31],42) # Ath, alltaf öruggt að taka fyrstu 31?
            self.featureMap[kash].append(list(self.df.columns)[1::][indiv].split('/')[-1])


    def taxCutForUniques(self):
        for member in range(self.M.shape[0]):
            uniques = self.M[member,member]
            self.M[:,member] = np.ones(self.M.shape[0])*-1*uniques//2
            self.M[member,member] = uniques
    
    def create_possibility_matrix(self):
        matrix = self.bitmatrix[self.valid_features,:]
        filtered = copy.deepcopy(matrix)
        filtered_list = list()
        self.count_list = list()
        count = 1

        for i in range(np.shape(filtered)[1]):
            filtered = filtered[filtered[:,i].argsort(kind='stable'),:]

        for j in range(1,np.shape(filtered)[0]):
            if not (np.array_equal(filtered[j,:], filtered[j-1,:])):
                filtered_list.append(j-1)
                self.count_list.append(count)
                count = 1
            else:
                count +=1
            if j == np.shape(filtered)[0] - 1:
                filtered_list.append(j)
                self.count_list.append(count)
        self.buffet = filtered[filtered_list,:]
    
    def create_feature_map(self):

        self.featureMap = collections.defaultdict(list)
        self.get_unique_features()
        self.get_mixed_features()

    def get_unique_features(self):
        
        self.GOAL = np.zeros(self.members)
        for g in range(len(self.GOAL)):
            self.GOAL[g] = min(10,self.buffet[:,g].dot(np.array(self.count_list)))
        
        # Náum í alla unique
        for unique in np.where(np.sum(self.buffet, axis=1) == 1)[0]:
            for i in range(self.count_list[unique]):
                vector = self.buffet[unique,:]
                if self.GOAL[np.where(vector == 1)[0]] > 0:
                    index = self.findVector(vector) # Finnur slíkan feature og skiptir út fyrir np.ones()
                    # Sameina þessi þrjú í eitt fall seinna
                    self.logConnections(vector) #
                    self.appendIndiv(vector, index) # Tengir fasta-skrár við einstaklinga
                    self.moveGoalPost(vector) # Uppfærir hversu marga þarf núna
                
            self.count_list[unique] = 0
            self.buffet[unique,:] = np.zeros(self.members)
        
        self.taxCutForUniques()
    
    def get_mixed_features(self):
        lowest_feature = np.ones((self.members))*2
        loops = np.zeros((self.members))
        

        while(np.any(self.GOAL > 0)):
            
            target = np.argmax(self.GOAL)
            
            workarray = np.where(np.sum(self.buffet, axis = 1) == lowest_feature[target])[0] # Skoðar gráðutölu, t.d. 2 (óþarfi?)
            closer_look = workarray[np.where(self.buffet[workarray,:][:,target] == 1)[0]] # Undirfylki fyrir target

            while len(closer_look) == 0:
                lowest_feature[target] += 1
                workarray = np.where(np.sum(self.buffet, axis = 1) == lowest_feature[target])[0]
                closer_look = workarray[np.where(self.buffet[workarray,:][:,target] == 1)[0]]
                
            workbench = self.buffet[workarray,:]
            
            best = float("inf")
            chosen = None
            
            for entry in range(np.shape(workbench)[0]):
                if workbench[entry,:].dot(self.M[target,:]) < best:
                    chosen = entry
                    best = workbench[entry,:].dot(self.M[:,target])
                
            vector = workbench[chosen,:]
            
            index = self.findVector(vector) # Finnur slíkan feature og skiptir út fyrir np.ones()
            self.logConnections(vector) #Merkir hvaða einstaklingar tengjast núna hverjum (gerist ekkert í 1's)
            self.appendIndiv(vector, index) # Tengir fasta-skrár við einstaklinga
            self.moveGoalPost(vector) # Uppfærir hversu marga þarf núna
        
            self.count_list[workarray[chosen]] -= 1
            if self.count_list[workarray[chosen]] == 0:
                self.buffet[workarray[chosen],:] = np.zeros(self.members)
            

    def fix_identicals(self):
        problems = list()
        
        for c in range(self.M.shape[0]): # Finnur duplicates
            found = False
            group = set()
            top = max(self.M[:,c])
            for r in range(self.M.shape[1]):
                if self.M[r,c] == top and r != c:
                    if np.array_equal(self.M[:,c],self.M[:,r]):
                        group.add(r)
                        group.add(c)
            for old_group in problems:
                if group == old_group:
                    old_group.update(group)
                    found = True
            if found == False and len(group) > 1:
                problems.append(group)
                
        for group in problems:
            for index in group:
                indices1 = np.where(self.buffet[:,index] == 1)[0]
                temp = self.buffet[indices1,:]
                for index2 in group:
                    indices2 = np.where(temp[:,index2] == 0)[0]
                    temp2 = temp[indices2,:]
                    for i in range(temp2.shape[0]):
                        cnt = self.count_list[indices1[indices2[i]]]
                        for j in range(min(2,cnt)):
                            vector = self.buffet[indices1[indices2],:][i,:]
                            index_a = self.findVector(vector)
                            # Sameina þessi þrjú í eitt fall seinna
                            self.logConnections(vector)
                            self.appendIndiv(vector, index_a)



    def write_results(self, outpath):
        cdbg_log = {'members': [x.split('/')[-1] for x in list(self.df.columns)[1::]],
                    'feature map': self.featureMap}
        
        if outpath == '0':
            outpath = self.name
        elif outpath[-1] == '/':
            outpath = outpath + self.name
        
        with open(outpath+'.atom.json', 'w') as wfile:
                json.dump(cdbg_log, wfile, indent = 4)
