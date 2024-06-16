# Algorithms and stuff
I often find myself rewriting / optimizing algorithms and data structures implemented in other languages in rust. I am hoping that this repo can help me organize them and perhaps make them easier to import and implement.

## SIMD & GPU optimized functions
Some of the implementations may utilize SIMD / cuda / OpenCL. If you find this repo, it may take a little work to get them to work for you. This repo just has basic implementations for these types of optimizations. When I use these for real world applications, they take a lot of case-by-case tuning, tweaking, and testing.

# Distance algorithms

During some work with bk trees I initially used the [bktree](https://crates.io/crates/bktree) crate as a starting point for another project. The layout was taken from that crate repo, but the only the levenshtein and hamming are from that repo. The different distance metrics are all under the Distance trait.

## A sequence modified Myer's algorithm for fast fixed length distance calculations

There are two distance methods: SequenceLevenshteinDistance & SequenceLevenshteinDistanceSimd
These distance functions are sequenced modified versions of the Myer's algorithm.
They have been adjusted to work in the context of next-gen sequencing (NGS) reads where deletions and insertions do not change the length of the string. This adjustment also allows windowing across strings to find sub-strings without breaking the metric properties of our distance metric. Which just means that our distance measurements will still obey euclidean geometry when windowing, which means we can use it with certain algos and data structures. This opens up some very high performance possibilities.

This is the original [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3853030/) for sequence modified levenshtein distance.

This [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10614987/) is really what helped me see how we can modify Myer's distance for NGS reads. It has nice LD matrix examples for different cases of compared same length sequences.

### How Myer's algo works & how we can modify it for sequence levenshtein distance
With sequencing data, we typically have millions of sequences to process, so we want fast methods for calculating distances (we also want to pair it w/ data structures like bk-trees). Myer's algorithm uses bitwise operations to calculate levenshtein edit distance, so its fast & easy to speed up with SIMD and GPU. If you're familiar with a levenshtein distance (LD) matrix, what Myer's does is calculate each value on the last column. The value we are looking for, LD, is the bottom right value in LD matrix. The values of the last column represent the # of insertions and deletions it would take to transform string B into string A. Since we just need to get to the corner value of the matrix, we can focus on the last column with just insertions and deletions.

Myer's algorithm uses an array of size 256 to index each ASCII character. The character ascii encoding will match its position in vector. It processes each character in first string one by one, inserting the ASCII code into a bit vector (the 256 array) and comparing it to the bit vector of the other string. Since we only care about insertions and deletions, all it has to do is track shifts in the bits to determine if an insertion or deletion occurred. If we insert a value into the array that already exist, no shift. If we insert a character that wasn't in the vector, then it shifts the bits, so we must have an insertion / deletion.

Now to modify for sequence levenshtein-distance, all we have to do is track the lowest score observed. If you look at this [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10614987/) you'll quickly notice that the sequence modified levenshtein distance is always the minimum value between the last row and last column, This means all we have to do is modify myers algorithm to track the minimum value. We also have to perform the calculation twice by swapping order of strings, since we need the last column and last row.

This algorithm can be further improved by using a index array of size 4 (ATGC). This actually helps a lot. The algorithm is CPU memory bound for my use-cases. The biggest time crunch is on moving memory between RAM > Cache > CPU. So only reserving memory for a sized 4 array would be a big improvement of a size 256.

There is also a SequenceLevenshteinDistanceWagner, which is a Wagner-Fischer algorithm modified for sequence levenshtein distance. This is for cases where your embedded sub-string is >32-64 characters.
