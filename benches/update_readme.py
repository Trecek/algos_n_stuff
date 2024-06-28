import json
import re


def load_benchmark_data():
    with open("benches/benchmark_results.json", "r") as f:
        return json.load(f)


def update_readme(data):
    with open("benches/README.md", "r") as f:
        content = f.read()

    # Update benchmark results
    content = re.sub(r"\[NUM_SEQUENCES\]", str(data["num_sequences"]), content)
    content = re.sub(r"\[SEQUENCE_LENGTH\]", str(data["sequence_length"]), content)
    content = re.sub(r"\[NUM_PAIRS\]", str(data["num_pairs"]), content)
    content = re.sub(r"\[TOTAL_BASES\]", str(data["total_bases"]), content)
    content = re.sub(
        r"\[BASES_PER_MICROSECOND\]", f"{data['bases_per_microsecond']:.2f}", content
    )
    content = re.sub(
        r"\[SEQUENCES_PER_MICROSECOND\]",
        f"{data['sequences_per_microsecond']:.2f}",
        content,
    )

    with open("benches/README.md", "w") as f:
        f.write(content)


def main():
    data = load_benchmark_data()
    update_readme(data)


if __name__ == "__main__":
    main()
