import argparse
import os
from pathlib import Path
import shutil
import tempfile
from phylopack.batch.batch import run_batching as batching
from phylopack.preorder.py_attotree import run_attotree as attotree


def add_refine_parser(subparsers):
    parser = subparsers.add_parser("refine", help="Refine genome ordering in windows using attotree")
    add_refine_args(parser)
    parser.set_defaults(func=run_refine)

def add_refine_args(parser):
    parser.add_argument('input_genomes', help='Path to the input list of genomes')
    parser.add_argument('output', help='Output path')
    parser.add_argument('--window', type=int, default=5000, help = 'Window size (number of genomes per batch) for refinement (default: 5000)')
    parser.add_argument('--min-size', type=int, default=10000, help = 'The minimum number of genomes that is in the ordering file to refine by window')

    #attotree setting
    parser.add_argument('-k', type=int, default=21, help='K-mer size (default: 21)')
    parser.add_argument('-s', type=int, default=10000, help='Sketch size (default: 10000)')
    parser.add_argument('-t', type=int, default=10, help='Number of threads (default: 10)')
    parser.add_argument('-m', choices=['nj', 'upgma'], default='nj', help='Tree method: nj or upgma (default: nj)')
    parser.add_argument('-v','--verbose', action='store_true', help='Print logs')
    parser.add_argument("--debug", action="store_true", help="Keep temp files for debugging")


def run_refine(args):

    if args.window < 4:
        raise ValueError("--window must be at least 4")

    input_path = args.input_genomes
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    basename = os.path.splitext(os.path.basename(input_path))[0]
    out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)
    
    with open(input_path, 'r') as f:
        total_genomes = sum(1 for line in f if line.strip())

    if total_genomes < args.min_size:
        if args.verbose:
            print(f"[INFO] Total genomes ({total_genomes}) < min_size ({args.min_size}). Copying directly to {args.output}")
        shutil.copyfile(input_path, args.output)
        return

    if args.debug:
        tmpdir = os.path.join(out_dir, "phylopack_tmp")
        os.makedirs(tmpdir, exist_ok=True)
    else:
        tmpdir = tempfile.mkdtemp()

    try:

        batch_args = argparse.Namespace(
            input = input_path,
            batch_name=basename,
            output_dir=tmpdir,
            target_size = args.window,
            verbose=args.verbose,  
            digits=4,
        )

        batching(batch_args)

        window_lists = sorted([
            os.path.join(tmpdir, f)
            for f in os.listdir(tmpdir)
            if os.path.isfile(os.path.join(tmpdir, f)) and f.startswith(f"{basename}_") and f.endswith(".txt")
        ])

        # Run attotree on each generated batch file

        treetmpdir = os.path.join(tmpdir, "trees")
        os.makedirs(treetmpdir, exist_ok=True)

        leaves_ordered_lists = []
        for window_file in window_lists:
            w_basename = os.path.splitext(os.path.basename(window_file))[0]
            w_tree_dir = os.path.join(treetmpdir, w_basename)
            os.makedirs(w_tree_dir, exist_ok=True)
            w_leaf_order = os.path.join(w_tree_dir, f"{w_basename}_leaf_order.txt")
            attotree_args = argparse.Namespace(
                input_genomes = window_file,
                output=w_tree_dir,
                k=args.k,
                s=args.s,
                t=args.t,
                m=args.m,
                verbose=args.verbose,
                statistic=False,
                statistic_file_type=getattr(args, 'statistic_file_type', 'json'),
                output_tree=None,
                output_std_tree=None,
                leaf_order=w_leaf_order,
                node_order=None,
            )

            attotree(attotree_args)
            leaves_ordered_lists.append(w_leaf_order)

        with open(args.output, "w") as out:
            for window_file in leaves_ordered_lists:
                with open(window_file, "r") as wf:
                    for line in wf:
                        line = line.strip()
                        if line:
                            out.write(line + "\n")

    finally:
        if not args.debug and os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)


def main():
    parser = argparse.ArgumentParser(
        description='Refine genome ordering in windows using attotree'
    )
    add_refine_args(parser)
    args = parser.parse_args()

    run_refine(args)

if __name__ == "__main__":
    main()