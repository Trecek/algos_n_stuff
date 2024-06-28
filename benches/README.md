# Benchmark Results

## Overview

1. Bit-packed DNA processing (single core)

## Latest Results


### Bit Ham Processor Performance

| Metric | Value |
|--------|-------|
| Sequences | 1536 |
| Sequence Length | 10 bases |
| Total Pairs | 1178880 |
| Total Bases Processed | 11788800 |
| Processing Speed (Bases) | 4597.68 bases/µs |
| Processing Speed (Sequences) | 459.77 sequences/µs |

### Performance Plot

![Benchmark Performance Plot](./benchmark_plot.png)

## Running the Benchmarks

To run the benchmarks and update this README, use the following command:

```
make update_benchmarks
```

This will compile the project with native optimizations, run the benchmarks, generate plots, and update this README with the latest results.

## Methodology

Single core, +avx2, compiled to native cpu
