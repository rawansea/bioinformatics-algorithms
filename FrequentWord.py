# Rawan Abu Alkhayr

# The Frequent Words Problem
# This program finds the most frequent k-mers in the given data,
# with O(|Text|.k) complexity.

#source:https://stackoverflow.com/questions/65641723/how-to-get-multiple-most-frequent-k-mers-of-a-string-using-python

#to read the dataset files
def read_file(filename):
    file = open(filename, 'r', encoding='utf-8')
    text = file.readline().strip()  #read the first line as text
    k = int(file.readline().strip())  #read the seconde line as K value
    return text, k


Text1, k1 = read_file("dataset_1.txt")
Text2, k2 = read_file("dataset_2.txt")


#find the most frequent k-mers in a given text
def FrequentWords(Text, k):
    count = {} #dictionary to store
    #loop through the text 
    for i in range(len(Text) - k + 1):
        Pattern = Text[i: i + k]
        if Pattern in count:
            count[Pattern] += 1     #already exist +1
        else:
            count[Pattern] = 1     #first time to see 1
            
    #find k-mers with the highest frequency 
    max_count = max(count.values())
    FrequentPatterns = []
    for pattern, freq in count.items():
        if freq == max_count:
            FrequentPatterns.append(pattern)
    
    return FrequentPatterns

print("Frequent patterns in first text:", FrequentWords(Text1, k1))
print("\n")
print("Frequent patterns in second text:", FrequentWords(Text2, k2))

