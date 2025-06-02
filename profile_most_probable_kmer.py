# Rawan Abu Alkhayr

# Profile-most Probable k-mer Problem
# This code finds the the most probable k-mer in text

# Source: 
# https://github.com/SinM9/Bioinformatic/blob/master/5.1%20Profile-most%20Probable%20k-mer%20Problem.py


# this func finds the most probable k-mer in DNA sequence using the profile matrix
def profile_most_probable_kmer(text, k, profile):
     
    max_probability = -1  # maximum probability
    most_probable_kmer = text[:k] # most probable k-mer

    # loop through all possible k-mers
    for i in range(len(text) - k + 1):
        
        kmer = text[i:i + k]  # current k-mer
        prob = 1  # current probability  

        # calculate the probability of the k-mer based on the profile matrix
        for j in range(k):

            nucleotide = kmer[j]  # get the nucleotide at position j

            # multiply the probability based on the nucleotide and its position
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


# get the input
text = input("Enter the DNA sequence: ").strip().upper()  

k = int(input("Enter k-mer length: ")) 

print("Enter the profile matrix:") 
matrix_input = [input() for _ in range(4)]  # read the 4 lines of string input

# convert input into a 2D matrix
profile = []
for row in matrix_input:
    row_values = list(map(float, row.split()))  # convert to float values
    profile.append(row_values)

 
# print the most probable k-mer
most_probable_kmer = profile_most_probable_kmer(text, k, profile)
print("\nMost Probable k-mer:", most_probable_kmer)