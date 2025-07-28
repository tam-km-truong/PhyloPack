import argparse
import os
import subprocess
import time
from datetime import datetime
import json

def get_args():
    parser = argparse.ArgumentParser(
        description='Run attotree with optional parameters and time statistics'
    )
    parser.add_argument('input_genomes', help='Path to the input list of genomes')
    parser.add_argument('-o', '--output', help='Output path (default: current folder)', default='.')
    parser.add_argument('-k', type=int, default=21, help='K-mer size (default: 21)')
    parser.add_argument('-s', type=int, default=10000, help='Sketch size (default: 10000)')
    parser.add_argument('-t', type=int, default=10, help='Number of threads (default: 12)')
    parser.add_argument('-m', choices=['nj', 'upgma'], default='nj', help='Tree method: nj or upgma (default: nj)')
    parser.add_argument('-v','--verbose', action='store_true', help='Print logs')
    return parser.parse_args()

def extract_timestamp(line):
    ts_part = ' '.join(line.split(' ')[1:3])
    return datetime.strptime(ts_part, "%Y-%m-%d %H:%M:%S")

def get_duration(lines, start_keyword, end_keyword):
    start_line = next((l for l in lines if start_keyword in l), None)
    end_line = next((l for l in lines if end_keyword in l), None)
    if start_line and end_line:
        start_time = extract_timestamp(start_line)
        end_time = extract_timestamp(end_line)
        return round((end_time - start_time).total_seconds(), 4)
    return None

def attotree():
    args = get_args()
    input_path = args.input_genomes
    output_basename = os.path.splitext(os.path.basename(input_path))[0]
    output_tree = os.path.join(args.output, f"{output_basename}.nw")
    output_std_tree = os.path.join(args.output, f"{output_basename}_std.nw")
    leaf_order = os.path.join(args.output, f"{output_basename}_leaf_order.txt")
    node_info = os.path.join(args.output, f"{output_basename}_node.txt")


    if args.verbose:
        print(f"Running attotree on {input_path}...")
        wall_start = time.time()
        cpu_start = os.times()

    cmd = [
        "attotree",
        "-L", input_path,
        "-o", output_tree,
        "-k", str(args.k),
        "-s", str(args.s),
        "-t", str(args.t),
        "-m", args.m
    ]


    # Run attotree and capture output
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    attotree_log = result.stdout + result.stderr

    # Run postprocess_tree.py
    cmd_postprocess = [
        "python", "postprocess_tree.py",
        "--standardize", "--midpoint-outgroup", "--ladderize", "--name-internals",
        "-l", leaf_order,
        "-n", node_info,
        output_tree,
        output_std_tree
    ]

    subprocess.run(cmd_postprocess,capture_output=True, check=True)

    # Re-add full paths to leaf_order.txt

    # Build filename full path mapping from original genome list
    path_map = {}
    with open(input_path) as f:
        for line in f:
            basename = os.path.basename(line.strip())
            genome_acession = basename.split('.')[0]
            path_map[genome_acession] = line.strip()

    # Patch leaf_order.txt in-place
    with open(leaf_order) as f:
        leaves = [line.strip() for line in f]

    missing = [leaf for leaf in leaves if leaf not in path_map]
    if missing:
        print("Warning: some leaves not found in input list:")
        for m in missing:
            print("  ", m)

    with open(leaf_order, 'w') as f:
        for leaf in leaves:
            f.write(path_map.get(leaf, leaf) + '\n')

    if args.verbose:
        wall_end = time.time()
        cpu_end = os.times()
        user_time = cpu_end.user - cpu_start.user
        system_time = cpu_end.system - cpu_start.system
        elapsed_time = wall_end - wall_start
        mash_time = get_duration(attotree_log.splitlines(), 'Running Mash', "Finished: 'mash triangle")
        quicktree_time = get_duration(attotree_log.splitlines(), 'Running Quicktree', "Finished: 'quicktree")

        stats = {
            "input": input_path,
            "k": args.k,
            "sketch_size": args.s,
            "threads": args.t,
            "method": args.m,
            "user_time": round(user_time, 4),
            "system_time": round(system_time, 4),
            "elapsed_time": round(elapsed_time, 4),
            "mash_triangle_time": mash_time,
            "quicktree_time": quicktree_time
        }

        print(json.dumps(stats, indent=2))

if __name__ == '__main__':
    attotree()