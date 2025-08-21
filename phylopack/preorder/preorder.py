import argparse
import tempfile
import os
import shutil
import json
import csv
from pathlib import Path

from phylopack.preorder.split_cluster import run_split
from phylopack.preorder.py_attotree import run_attotree
from phylopack.preorder.placement import run_placement

def add_preorder_parser(subparsers):
    preorder_parser = subparsers.add_parser("preorder", help="Run full pipeline")
    _add_common_args(preorder_parser)
    preorder_parser.set_defaults(func=run_preorder_pipeline)

def _add_common_args(parser):
    parser.add_argument("input_genomes", help="Path to input genome list")
    parser.add_argument("-o", "--output", required=True, help="Output file for final genome preorder list")
    parser.add_argument(
        "--cut-point", type=float, default=0.01,
        help="Cut size: float <1 for percentage, int ≥1 for number (default: 0.01)"
    )
    parser.add_argument("--seed", type=int, help="Random seed for shuffling (default: current timestamp)")
    parser.add_argument("-k", type=int, default=21, help="K-mer size (default: 21)")
    parser.add_argument("-s-reference", type=int, default=10000, help="Sketch size for reference genomes (default: 10000)")
    parser.add_argument("-s-placement", type=int, default=1000, help="Sketch size placement(default: 1000)")
    parser.add_argument("-t", type=int, default=10, help="Threads (default: 10)")
    parser.add_argument("-m", choices=["nj", "upgma"], default="nj", help="Tree method (default: nj)")
    parser.add_argument("--statistic", action="store_true", help="Enable statistics")
    parser.add_argument(
        "--statistic-file-type", choices=["json", "csv"], default="csv",
        help="Statistics file format (default: json)"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    parser.add_argument("--debug", action="store_true", help="Keep temp files for debugging")

    parser.set_defaults(func=run_preorder_pipeline)

def concat_stat_files(paths, output_path, file_type="json"):
    if file_type == 'json':
        merged = {}
        for path in paths:
            name = Path(path).stem
            with open(path) as f:
                merged[name] = (json.load(f))
        with open(output_path, "w") as out:
            json.dump(merged, out, indent=4)
    elif file_type == 'csv':
        with open(output_path, "w", newline="") as out_csv:
            writer = None
            for path in paths:
                name = Path(path).stem
                with open(path, newline="") as in_csv:
                    reader = csv.reader(in_csv)
                    headers = next(reader)
                    if writer is None:
                        writer = csv.writer(out_csv)
                        writer.writerow(['stat_file'] + headers )
                    for row in reader:
                        writer.writerow([name] + row)

def run_preorder_pipeline(args):
    if args.debug:
        tmpdir = os.path.join(os.path.dirname(args.output), "phylopack_tmp")
        os.makedirs(tmpdir, exist_ok=True)
    else:
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        tmpdir = tempfile.mkdtemp()

    if args.verbose:
        print(f'[INFO] Temp directory at: {tmpdir}')

    ref_file = os.path.join(tmpdir, "references.txt")
    rem_file = os.path.join(tmpdir, "remains.txt")
    leaf_order_file = os.path.join(tmpdir, "leaf_order.txt")
    node_order_file = os.path.join(tmpdir, "node_order.txt")
    output_tree = os.path.join(tmpdir, "tree.nw")
    output_std_tree = os.path.join(tmpdir, "tree_std.nw")
    final_output_tmp = os.path.join(tmpdir, "placement_order.txt")

    split_args = argparse.Namespace(
        input_genomes = args.input_genomes,
        cut_point=args.cut_point,
        output=tmpdir,
        seed=args.seed,
        verbose=args.verbose,
        statistic=args.statistic,
        statistic_file_type=args.statistic_file_type,
        ref_output=ref_file,
        rem_output=rem_file
    )

    run_split(split_args)

    attotree_args = argparse.Namespace(
        input_genomes = ref_file,
        output=tmpdir,
        k=args.k,
        s=args.s_reference,
        t=args.t,
        m=args.m,
        verbose=args.verbose,
        statistic=args.statistic,
        statistic_file_type=args.statistic_file_type,     
        output_tree=output_tree,
        output_std_tree=output_std_tree,
        leaf_order=leaf_order_file,
        node_order=node_order_file
    )

    run_attotree(attotree_args)

    placement_args = argparse.Namespace(
        genomes_list_1=rem_file,
        genomes_list_2=leaf_order_file,
        output=tmpdir,
        k=args.k,
        s=args.s_placement,
        t=args.t,
        verbose=args.verbose,
        statistic=args.statistic,
        statistic_file_type=args.statistic_file_type,  
    )

    run_placement(placement_args)

    shutil.copyfile(final_output_tmp, args.output)

    if args.verbose:
        print(f"[INFO] Preorder written to {args.output}")

    if args.statistic:
        stat_paths = [
            os.path.join(tmpdir, "split_stats." + args.statistic_file_type),
            os.path.join(tmpdir, "tree_stats." + args.statistic_file_type),
            os.path.join(tmpdir, "placement_stats." + args.statistic_file_type),
        ]

        basename = os.path.splitext(os.path.basename(args.input_genomes))[0]

        cutpoint = args.cut_point

        if cutpoint >= 1:
            cutpoint = int(cutpoint)

        merged_stat = os.path.join(os.path.dirname(args.output), f"preorder_stat_{basename}_{cutpoint}.{args.statistic_file_type}")
        concat_stat_files(stat_paths, merged_stat, args.statistic_file_type)
        if args.verbose:
            print(f"[INFO] Statistic written to {merged_stat}")

    if not args.debug:
        shutil.rmtree(tmpdir)

def main():
    parser = argparse.ArgumentParser(description="Run full preorder pipeline")
    _add_common_args(parser)
    args = parser.parse_args()
    run_preorder_pipeline(args)

if __name__ == "__main__":
    main()