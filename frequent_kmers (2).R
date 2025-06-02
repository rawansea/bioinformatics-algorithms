# Rawan Abu Alkhayr

# The Frequent Words Problem
# This program finds the most frequent k-mers in the given data,
# with O(|Text|.k) complexity.
 

#find the most frequent k-mers in a given text
find_frequent_kmers <- function(text, k) {
  kmer_counts <- list()
  max_count <- 0
  
  #sliding window to extract and count all k-mers
  for (i in 1:(nchar(text) - k + 1)) {
    kmer <- substr(text, i, i + k - 1)
    #first time to see 1
    if (is.null(kmer_counts[[kmer]])) {
      kmer_counts[[kmer]] <- 1
    } else {
      kmer_counts[[kmer]] <- kmer_counts[[kmer]] + 1 #already exist +1
    }
    #to track the max
    if (kmer_counts[[kmer]] > max_count) {
      max_count <- kmer_counts[[kmer]]
    }
  }
  
  #find k-mers with the highest frequency 
  frequent_kmers <- c()  
  
  for (kmer in names(kmer_counts)) {
    if (kmer_counts[[kmer]] == max_count) {
      frequent_kmers <- c(frequent_kmers, kmer)  
    }
  }
  
  return(frequent_kmers)
}


#to read text from the files 
process_file <- function(file_path, k) {
  if (!file.exists(file_path)) {
    stop(paste("Error: File", file_path, "not found!"))
  }
  
  text <- readLines(file_path, warn = FALSE)
  text <- gsub("\\s+", "", paste(text, collapse = ""))  #remove whitespace and merge lines
  
  return(find_frequent_kmers(text, k))
}

#define k value  
k <- 12

#process dataset files and display results
dataset_files <- c("dataset_1.txt", "dataset_2.txt")
for (file in dataset_files) {
  cat("Processing:", file, "\n")
  result <- process_file(file, k)
  cat("Most Frequent k-mers:", paste(result, collapse = ", "), "\n\n")
}