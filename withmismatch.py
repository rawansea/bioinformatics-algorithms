# Rawan Abu Alkhayr

# Frequent Words with Mismatches Problem
# This code gonna output all most frequent k-mers with allowing d mismatches in Text.
# source; https://github.com/kswang2400/bioinformatics-code-challenges/blob/master/text%20files/Code%20Challenge%201.7:%20Frequent%20Words%20with%20Mismatch

# check if two k-mers match with up to d mismatches
def is_match(s1, s2, k, d):
    mismatch = 0 #mismatch counter
    for i in range(k):
        if s1[i] != s2[i]:
            mismatch += 1
        if mismatch > d: 
            return False
    return True

# func to find the most frequent k-mers with mismatches
def frequent_words(text, k, d):

    # count exact occurrences of each k-mer in the text
    c= {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        if kmer in c:
            c[kmer] += 1
        else:
            c[kmer] = 1

    # count include mismatches
    new_c = {}
    for str1 in c:
        total = 0
        for str2 in c:
            if is_match(str1, str2, k, d): # = True
                total += c[str2] 
        new_c[str1] = total

    # find the maximum frequency 
    max_freq = max(new_c.values())
    # all the maximum frequency
    frequent_kmers = []
    for kmer, count in new_c.items():
        if count == max_freq:
            frequent_kmers.append(kmer)

    return frequent_kmers
    
    
# input
text = input("Enter the string text:")
k = int(input("Enter the k value:"))
d = int(input("Enter the d value:"))

# call the func to find the most frequent k-mers with mismatches
result = frequent_words(text, k, d)

print(result)