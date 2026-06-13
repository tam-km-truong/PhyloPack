import argparse
import math
from pathlib import Path

def add_batch_parser(subparsers):
    batch_parser = subparsers.add_parser("batch", help="Run batch")
    add_batch_args(batch_parser)
    batch_parser.set_defaults(func=run_batching)

def add_batch_args(parser):
    parser.add_argument(
        "input",
        help="Ordered genome list (one genome per line)"
    )

    parser.add_argument(
        '-o',"--output-dir",
        required=True,
        help="Directory containing batch files"
    )

    parser.add_argument(
        '-s',"--target-size",
        type=int,
        default=10000,
        help="Target number of genomes per batch"
    )

    parser.add_argument(
        "-n", "--batch-name",
        default="batch",
        help="Prefix for batch filenames"
    )

    parser.add_argument(
        "-d","--digits",
        type=int,
        default=3,
        help="Number of digits in batch suffix (like split -d)"
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print batch statistics"
    )

def balanced_partition(n_items, max_batch_size):
    """
    Split n_items into the minimum number of balanced batches
    whose sizes differ by at most 1 and never exceed max_batch_size.
    """
    if n_items == 0:
        return []

    n_batches = math.ceil(n_items / max_batch_size)

    base_size = n_items // n_batches
    remainder = n_items % n_batches

    sizes = []

    for i in range(n_batches):
        if i < remainder:
            sizes.append(base_size + 1)
        else:
            sizes.append(base_size)

    return sizes

def run_batching(args):
    with open(args.input) as f:
        genomes = [line.strip() for line in f if line.strip()]

    n_genomes = len(genomes)

    if args.target_size <= 0:
        raise ValueError("--target-size must be > 0")

    batch_sizes = balanced_partition(
        n_genomes,
        args.target_size
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    start = 0

    for batch_idx, batch_size in enumerate(batch_sizes, start=1):
        end = start + batch_size

        batch = genomes[start:end]

        suffix = f"{batch_idx:0{args.digits}d}"
        batch_path = output_dir / f"{args.batch_name}_{suffix}.txt"

        with open(batch_path, "w") as out:
            out.write("\n".join(batch))
            out.write("\n")

        if args.verbose:
            print(
                f"{args.batch_name}_{suffix}: "
                f"{batch_size} genomes"
            )

        start = end

def main():
    parser = argparse.ArgumentParser(
        description='Split an ordered genome list into balanced contiguous batches.'
    )
    add_batch_args(parser)
    args = parser.parse_args()

    run_batching(args)

if __name__ == "__main__":
    main()