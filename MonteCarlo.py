# Rawan Abu Alkhayr


# Randomized Search using Monte Carlo 
# source: https://rosalind.info/problems/ba2f/
 

import random

# this func to create profile matrix from given motifs
def Create_Profile(motifs):
    k = len(motifs[0])  
    t = len(motifs)  
    profile = [[0] * k for _ in range(4)]  # profile matrix

    # count nucleotide occurrences in each column
    for motif in motifs:
        for i in range(k):
            nucleotide = motif[i]
            if nucleotide == 'A':
                profile[0][i] += 1
            elif nucleotide == 'C':
                profile[1][i] += 1
            elif nucleotide == 'G':
                profile[2][i] += 1
            elif nucleotide == 'T':
                profile[3][i] += 1

    return profile


# this func to calculate the score of motifs based on mismatches with the consensus string
def Score(motifs):
    mismatch_count = 0  
    consensus_string = ""  # store the consensus string
    profile = Create_Profile(motifs)
    
    # iterate over each position in motifs to find consensus string
    for i in range(len(profile[0])):   
        max_p = 0 
        most_frequent_nucleotide = 0 
          
        # check which nucleotide is most frequent in the current column
        for j in range(4):  
            if profile[j][i] > max_p:
                max_p = profile[j][i]
                most_frequent_nucleotide = j # store the index
        
        # append the nucleotide to the consensus string
        if most_frequent_nucleotide == 0:
            consensus_string += "A"
        elif most_frequent_nucleotide == 1:
            consensus_string += "C"
        elif most_frequent_nucleotide == 2:
            consensus_string += "G"
        elif most_frequent_nucleotide == 3:
            consensus_string += "T"
    
    # count mismatches  
    for motif in motifs:
        for i in range(len(motif)):
            if consensus_string[i] != motif[i]:
                mismatch_count += 1  
    
    return mismatch_count


# this func to find the most probable k-mer 
def profile_most_probable_kmer(text, k, profile):
    max_probability = -1  
    most_probable_kmer = text[:k]  

    # loop through all possible k-mers
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]  # current k-mer
        prob = 1  # current probability  

        # calculate the probability of the k-mer based on the profile matrix
        for j in range(k):
            nucleotide = kmer[j]  # get the nucleotide at position j
            if nucleotide == 'A':
                prob *= profile[0][j]
            elif nucleotide == 'C':
                prob *= profile[1][j]
            elif nucleotide == 'G':
                prob *= profile[2][j]
            elif nucleotide == 'T':
                prob *= profile[3][j]

        # update the maximum probability
        if prob > max_probability:
            max_probability = prob
            most_probable_kmer = kmer

    return most_probable_kmer


# monte carlo search algorthim
def MonteCarlo(Dna, k, t):
    iterations = 10000
    # randomly select initial motifs and calculate their score
    motifs = []
    for i in range(t):
        k_mers = []   
        for j in range(len(Dna[i]) - k + 1):   
            k_mer = Dna[i][j:j+k]   
            k_mers.append(k_mer)  
        
        motif = random.choice(k_mers)   
        motifs.append(motif) 
    
    best_motifs = motifs
    best_score = Score(best_motifs)
    
    no_improvement_count = 0 
    
    # run multiple itereations
    for iteration in range(iterations):
        profile = Create_Profile(motifs)
        # new motifs based on the profile
        new_motifs = []
        for i in range(t):
            new_motif = profile_most_probable_kmer(Dna[i], k, profile)
            new_motifs.append(new_motif)
      
        # update the best motifs 
        current_score = Score(new_motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = new_motifs 
            no_improvement_count = 0  
        else:
            no_improvement_count += 1  # increment
        
        # no improvement threshold
        if no_improvement_count >= 100000: 
            break

    return best_motifs

# input for t, k, and Dna strings
t = int(input("Enter number of Dna strings: "))
k = int(input("Enter k-mer length: "))  
print("Enter the Dna strings:")
Dna = [input().strip().upper() for _ in range(t)]  


# print the result
best_motifs = MonteCarlo(Dna, k, t)
print("Best Motifs:", best_motifs)
