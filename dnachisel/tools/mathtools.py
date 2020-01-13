from random import random

#select randomly based on probability, prob_list = [('1',0.1),('2',0.5),('3',0.4)]
def pick_random_tuple(prob_list):
    r, s = random(), 0
    for num in prob_list:
        s += num[1]
        if s >= r:
            return num[0]
        
#select randomly based on probability, prob_list = [0.1,0.5,0.4]
def pick_random(prob_list):
    r, s, i = random(), 0, 0
    for num in prob_list:
        s += num
        if s >= r:
            return i
        i += 1        

def hammingDistance(seq1,seq2):
    score_nt = 0

    for i in range(0,len(seq1)):
        if seq1[i] != seq2[i]:
            score_nt += 1
            
    return score_nt