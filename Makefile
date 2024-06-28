.PHONY: tests precommit setup bench update_benchmarks

setup:
	@command -v rustc >/dev/null 2>&1 || { echo >&2 "Rust is not installed. Installing Rust..."; curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh; }
	cargo install cargo-udeps
	cargo install flamegraph
	pip install matplotlib

# This will check for dead dependencies
check_udeps:
	@echo "Checking dead dependencies without native optimizations:"
	cargo +nightly udeps --all-targets
	@echo "Checking dead dependencies with native optimizations:"
	RUSTFLAGS="-C target-cpu=native" cargo +nightly udeps --all-targets

# This will run benchmarks with native optimizations
bench_native:
	RUSTFLAGS="-C target-cpu=native" cargo bench

# This will run benchmarks with native optimizations
bench:
	cargo bench

build_native:
	RUSTFLAGS="-C target-cpu=native" cargo build --release

bench_native_flame:
	RUSTFLAGS="--emit=llvm-ir -C target-cpu=native" cargo flamegraph --root --bench benchmark

bench_flame:
	RUSTFLAGS="--emit=llvm-ir" cargo flamegraph --root --bench benchmark

# New targets for updating benchmarks and README
update_benchmarks: bench_native generate_plots update_readme

generate_plots:
	python benches/plot_benchmark.py

update_readme:
	python benches/update_readme.py

# This will check for various SIMD features including AVX2, SSE versions, and others
check_simd_features:
	@echo "Checking for SIMD features support..."
	@echo -n "AVX2: "; cat /proc/cpuinfo | grep -q avx2 && echo "Supported" || echo "Not supported"
	@echo -n "SSE: "; cat /proc/cpuinfo | grep -q sse && echo "Supported" || echo "Not supported"
	@echo -n "SSE2: "; cat /proc/cpuinfo | grep -q sse2 && echo "Supported" || echo "Not supported"
	@echo -n "SSE3: "; cat /proc/cpuinfo | grep -q sse3 && echo "Supported" || echo "Not supported"
	@echo -n "SSSE3: "; cat /proc/cpuinfo | grep -q ssse3 && echo "Supported" || echo "Not supported"
	@echo -n "SSE4.1: "; cat /proc/cpuinfo | grep -q sse4_1 && echo "Supported" || echo "Not supported"
	@echo -n "SSE4.2: "; cat /proc/cpuinfo | grep -q sse4_2 && echo "Supported" || echo "Not supported"
	@echo -n "AVX: "; cat /proc/cpuinfo | grep -q avx && echo "Supported" || echo "Not supported"
	@echo -n "FMA: "; cat /proc/cpuinfo | grep -q fma && echo "Supported" || echo "Not supported"
