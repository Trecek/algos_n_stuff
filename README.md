# Algorithms and Data Structures in Rust

I often find myself rewriting and optimizing algorithms and data structures implemented in other languages into Rust. This repo helps organize these implementations and makes them easier to import and use.

## Algorithms
- **Myer's 1999 algorithm modified for sequence modified Levenshtein distance (seq-lev)**
- **A SIMD variant of sequence modified Levenshtein distance with windowing (a modified Myer's algorithm)**
- **Bit-packing and SIMD-accelerated Hamming distance:** very fast. Bases are encoded into 3 bits and packed continuously into `u64` values, allowing bit-by-bit comparison and summing scores based on index locations within `u64`s.

**TODO**:
- Precompute neighborhood methods
- Mutation methods
- Streaming/channel methods for computing seq-lev distance
- DNA set generations with minimum edit distances using greedy evolutionary algorithms
- BK-tree variant utilizing cosine law (reducing distance calculations) and GPU

## Distance Algorithms

### Sequence Modified Myer's Algorithm for Fast Fixed-Length Distance Calculations

The distance functions `SequenceLevenshteinDistance` and `SequenceLevenshteinDistanceSimd` are sequence-modified versions of Myer's algorithm. 

These versions are adjusted for next-gen sequencing (NGS) reads, where deletions and insertions do not change the string length. This allows windowing across strings to find substrings without breaking the metric properties, enabling use with certain algorithms and data structures for high-performance possibilities.

For more details, refer to the original [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3853030/) on sequence modified Levenshtein distance and another [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10614987/) that explains how to modify Myer's distance for NGS reads.

### How Myer's Algorithm Works and Its Modification for Sequence Levenshtein Distance

With sequencing data, we process millions of sequences, so fast distance calculation methods are crucial, especially when paired with data structures like BK-trees. Myer's algorithm uses bitwise operations to calculate Levenshtein edit distance, which can be sped up with SIMD and GPU.

Myer's algorithm uses an array of size 256 to index each ASCII character, processing characters in the first string one by one. Insertions and deletions are tracked by shifts in the bits, determining if an insertion or deletion occurred.

To modify for sequence Levenshtein distance, we track the lowest score observed. The sequence modified Levenshtein distance is always the minimum value between the last row and column, requiring us to track the minimum value and perform the calculation twice by swapping the order of strings.

This algorithm can be further improved by using an index array of size 4 (ATGC), significantly reducing memory usage and improving performance. For cases where the embedded substring is >32-64 characters, `SequenceLevenshteinDistanceWagner` (a Wagner-Fischer algorithm modified for sequence Levenshtein distance) is used.
