
def patternToNumber(pattern):
    if len(pattern)== 0:
        return 0 #wenn Länge Null entspricht dann gebe Null aus
    return 4* patternToNumber(pattern[0:-1])+ symbolToNumber(pattern[-1:])

def symbolToNumber(symbol): #Symbol zu Zahl Umwandlung
    if symbol=="A":
        return 0
    if symbol == "C":
        return 1
    if symbol == "G":
        return 2
    if symbol == "T":
        return 3

def numberToPattern(number,k): #"Umrechnung" der Zahl in ein Pattern
    if k==1:
        return numberToSymbol(number)
    return numberToPattern(number//4, k-1) + numberToSymbol(number%4)

def numberToSymbol2(number): #Alternative für Switch
    if number == 0:
        return "A"
    if number == 1:
        return "C"
    if number == 2:
        return "G"
    if number == 3:
        return "T"

def numberToSymbol(number):
    switch={
        O:"A",
        1:"C",
        2:"G",
        3:"T"}
    return switch.get(number)

def profileProbable (text, k, profile):
    maxprob=0
    kmer=text[0:k]
    for i in range(0, len(text)-k+1): #list comprehension
        prob = 1
        pattern= text[i:i+k]
        for j in range(k):
            l=symbolToNumber(pattern[j])
            prob*=profile[l][j] #ich finde den Fehler einfach nicht!!!!!!!!!
        if prob>maxprob:
            maxprob=prob #prob wird ab hier zu neuem Maximum da es größer als das bisherige max ist
            kmer=pattern
        return  kmer

def hemmingDistance(p, q):
    ham=0 #Startwert festlegen für Distanz, damit auf jeden FAll überschritten
    for number, y in zip(p,q):
        if number != y:
            ham +=1
    return ham

def PatternStringDistance(pattern, DNA):
    k= len(pattern)
    distance=0 #Startwert
    for number in DNA:
        hamming=k+1
        for i in range(len(number)-k+1):
            z=hammingDistance(pattern, number[i:i+k])
            if hamming>z:
                hamming=z
        distance += hamming
    return distance

def profileMotifs(motifs):
    k=len(motifs[0])
    profile=[[1 for i in range(k)] for j in range(4)]
    for number in motifs:
        for i in range(len(number)):
            j=symbolToNumber(number[i])
            profile[j][i] +=1
        for number in profile:
            for i in range(len(number)):
                number[i]=number[i]/len(motifs)

def consensus(profile):
    str=""
    for i in range(len(profile[0])):
        max=0
        loc=0
        for j in range(4):
            if profile[j][i]> max:
                loc= j
                max= profile [j][i]
        str += numberToSymbol(loc)
    return str

def score(motifs):
    profile=profileMotifs(motifs)
    cons=consensus(profile)
    score=0
    for number in motifs:
        for i in range(len(number)):
            if cons[i] != number[i]:
                score +=1
    return score

import random

def randomMotifSearch(DNA, k, t):
    winMotifs = []
    motifs= []
    for number in range(t):
        random.seed()
        i= random.randint(0, len(DNA[number])-k)
        motifs.append(DNA[number][i:i+k])
    winMotifs= motifs.copy()
    count=0
    while True: #Schleife
        profile= profileMotifs(motifs)
        for number in range(t):
            motifs[number]=profileProbable(DNA[number], k, profile)
        if score(motifs) < score(winMotifs):
            winMotifs= motifs.copy()
            count +=1
        else:
            print(count)
            return winMotifs

k=8
t=5
DNA= ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
    #input Rosalind
best=randomMotifSearch(DNA, k, t)
min= score(best)
for number in range(1000):
    print(number)
    x= randomMotifSearch(DNA, k, t)
    if score(x) < score(best):
        best= x
        min= score(x)
print(min)
for number in best:
    print(number)

k=15
t=20
DNA=["ACTTATATCTAGAGTAAAGCCCTGATTCCATTGACGCGATCCCTACCTCCATCATACTCCACAGGTTCTTCAATAGAACATGGGGAAAACTGAGGTACACCAGGTCTAACGGAGATTTCTGGCACTAACTACCCAAAATCGAGTGATTGAACTGACTTATATCTAGAGT",
     "AAAGCCCTGATTCCATTGACGCGATCCCTACCTCCATCATACTCCACAGGTTCTTCAATAGAACATGGGGAAAACTGAGGTACACCAGGTCTAACGGAGATTTCTGGCACTAACTACCCAAAATCCTCTCGATCACCGACGAGTGATTGAACTGACTTATATCTAGAGT",
     "CACTCCCGTCCGTCTGACGCCAGGTGCTCTACCCCGCTGATTGTCTGGTACATAGCAGCCTATAGATCACCGATGCAGAAACACTTCGAGGCAGCCGATTTCGCTTATCACAACGTGACGGAATTTGATAAACCACGTACTCTAATACCGTCACGGGCCCATCAACGAA",
     "ACAAGAACTGGTGGGGAGACTATGACACTCTAGCGGTCGCATAAGGGCCGGAAACCAGGACAAATCGATAAGATGAAGCGGGGATATAAGCCTTATACTGCGACTGGTTCCTTATATTATTTAGCCCCGATTGATCACCGATTAAAATATTCTGCGGTTTTCGAGACGG",
     "TAACCACACCTAAAATTTTTCTTGGTGAGATGGACCCCCGCCGTAAATATCAGGATTAAATGTACGGATACCCATGACCCTCCAGTCATCTACCTTCCCGTGGTGGTCGCTCAGCCTTGTGCAGACCGAACTAGCACCTGTCACATACAATGTTGCCCGCATAGATCGT",
     "ATCCGACAGAGGCAGTGAATAAGGTTTCGTTTCCTCAGAGAGTAGAACTGCGTGTGACCTTGCCTTCACCGACATCCGTTTCCAATTGAGCTTTTCAGGACGTTTAGGTAACTGATTGTCATTGCAATTGTCCGGGGGATTTAGATGGCCGGGTACCTCTCGGACTATA",
     "CCTTGTTGCCACCGATTCGCGAGCAACATCGGAGTGCTCTGATTCACGGCGATGCTCCACGAAGAGGACCGCGGCACGACACGCCCTGTACCTACGTTTCTGGATATCCTCCGGCGAGTTAATAGAGCAATACGACCTGGTCGTCGAGATCGTGTATCTAGCCCTACCT",
     "ATAGGTTAACGAATCAGGAGAGTTAATTTTACCTAGCTAGAGCGGACGGTGCCTGGCTGTATTCGCGTTTGACTTTCGGGCTCGCTGATAACTTGTGATCACCTTTTACGCTTACTGGATCCAACGATGGATCAAAGTTGAGAATTTCTGTGCCTTGGGTGTGAGCTGT",
     "CTGACGAAAGGACGGGCGGTGTACTTAGTTTGGGGTAAAATAGTTGGTATAATTCTGTGCGACAGACATTTGGTCAGGCCATACTGCCATATCGTGATGTAACTATCCACACTACGTCATAGGCCCTTGTGATCAATTAAACGTTCCTCATGCCAGGCTATCTGTTTAA"
     "GGCTTCGCGTTTAAGGCTGGATTAAGTACTCCGCCTTGTGATCTGTGATCCTCCGACCTGTGATCAGCAAGATTGGAACCTAGGTAGGCGGCGGGTCTACGCTGGCCCACAATCGTGAGTCCCCCACTCCGTAGGTTGTGGAATTTATAGACCCGCAAGGGGCACCACT"
     "AGGATGACACCCAGGATGAATCTGGATTAGGAACACCAACCCGACATATTTGTTACCGCTGCAGCATTTCGCTCTTGGACGCGTAACCCGAGATCCGTCTCGCGATCGTCACGGATCGGGATTATGCAGGCAATACCTTGTGATCACTCCGCGCTTGGTTTTGCTAGCG"
     "ACATCTCTAGTCACTTTTATTGAGCAGGTGGGCGGATTCATGATCCGGCTCTGTCGTACGTCCAACCACGGTGACATGTTCGGAGCTGTCGCCGTGGAGCAGAGATACATCGGATCTATCAATTTTACTAAGAGCAACTAGCCACGACAAACTGTGATCACCGATTGGA"
     "AATTTGCGTATCTCTAGGACTCCCTCATACAAATCAAAGCTTGGATGGGTAAGATGCCGCAGCAGCAGGTATCTCATATTGGCTATTAAGAGCCAGGCCCTATGGCCTTAGTATCACCGATCAGACGTCGCATGAGCGGGCCCGTTGTCCTATCTCTTTAGCTGCCGCA"
     "GAAGTAAAGGGGTTCCACTGCGTAGAGCGTGCCCCTCTGGTGTGCCGTACTGTTATGGTGATACAGCTTCCTTATACCCCTCGTAAAGCGGCTAATGGTCCTAATGAATGCCCTTGTGAAATCCGAATCGCTTTACAATTGCGTTCGGCGGAATGCAGTCACCAGTGTT"
     "TACACTACGCGTTATTTACTTTTACTGAGTCCTTGTCGCCACCGAACGAGGATTGTTCATTGTATCCGGAGATTAGGAGTTCGCATCGCTGACACAGCCAGTTCGTAGCAAATACCGCTGGCCCTGGGCACTCCAGATCAGAACTACTAGCCCTAAACTCTATGACACA"
     "TTGGGTCTCGATCCCTCTATGTTAAGCTGTTCCGTGGAGAATCTCCTGGGTTTTATGATTTGAATGACGAGAATTGGGAAGTCGGGATGTTGTGATCACCGCCGTTCGCTTTCATAAATGAACCCCTTTTTTTCAGCAGACGGTGGCCTTTCCCTTTCATCATTATACA"
     "TTTCAAGTTACTACCGCCCTCTAGCGATAGAACTGAGGCAAATCATACACCGTGATCACCGACCCATGGAGTTTGACTCAGATTTACACTTTTAGGGGAACATGTTTGTCGGTCAGAGGTGTCAATTATTAGCAGATATCCCCCAACGCAGCGAGAGAGCACGGAGTGA"
     "GATCCATTACCCTACGATATGTATATAGCGCCCTAGTACGGCTTCTCCCTTGCAGACACGCAGGCGCTGTGCGCTATCGGCTTCCTCGGACATTCCTGGATATAAGTAACGGCGAACTGGCTATCACTACCGCCGCTCCTTAAGCCTTGGTTTCACCGACGATTGTCGT"
     "TAGTAGATTATTACCTGTGGACCGTTAGCTTCAAGACCGAAACGTTGGTGATGCTACTTAAATGTCAAGAGTTGCGAAGTTGGGCGAAGCACATCCGTACTCCCAAGTGGACGATCGATAGATCCATGGAGTTTCCATCCATCTTAATCCGCCCTTTGCATCACCGACG"
     "TACAAGGCACAAACGAGACCTGATCGAACGGTGCACGGTCGAGGCAGCGAGATAAATGTACATTGAGAGCACCTTGTGATTTACGACCTGCATCGAAGGTTTCTTGGCACCCACCTGTCGTCCGCCAGGGCAGAGCCGACATTATATGACGCTGATGTACGAAGCCCCT"]
    #extra input rosalind

best=randomMotifSearch(DNA, k, t)
min= score(best)
for number in range(1000):
    print(number)
    x= randomMotifSearch(DNA, k, t)
    if score(x) < score(best):
        best= x
        min= score(x)
print(min)
for number in best:
    print(number)