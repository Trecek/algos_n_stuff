import matplotlib.pyplot as plt
import json
import os


def load_benchmark_data():
    with open("benches/benchmark_results.json", "r") as f:
        return json.load(f)


def generate_performance_plot(data):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

    # Plot bases per microsecond
    ax1.bar(data["algorithm_names"], data["processing_speeds"])
    ax1.set_title("DNA Processing Algorithm Performance (Bases)")
    ax1.set_xlabel("Algorithm")
    ax1.set_ylabel("Processing Speed (bases/µs)")
    ax1.tick_params(axis="x", rotation=45)

    # Calculate sequences per microsecond
    seq_speeds = [
        speed / data["sequence_length"] for speed in data["processing_speeds"]
    ]

    # Plot sequences per microsecond
    ax2.bar(data["algorithm_names"], seq_speeds)
    ax2.set_title("DNA Processing Algorithm Performance (Sequences)")
    ax2.set_xlabel("Algorithm")
    ax2.set_ylabel("Processing Speed (sequences/µs)")
    ax2.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    plt.savefig("benches/benchmark_plot.png")


def main():
    data = load_benchmark_data()
    generate_performance_plot(data)


if __name__ == "__main__":
    main()
