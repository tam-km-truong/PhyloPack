import argparse
import os
import subprocess
import time
from datetime import datetime
import json
import csv
import resource
from pathlib import Path
import sys
from phylopack.preorder.postprocess_tree import run as postprocesstree

def add_tree_args(parser):
    parser.add_argument('input_genomes', help='Path to the input list of genomes')
    parser.add_argument('-o', '--output', help='Output path (default: current folder)', default='.')
    parser.add_argument('-k', type=int, default=21, help='K-mer size (default: 21)')
    parser.add_argument('-s', type=int, default=10000, help='Sketch size (default: 10000)')
    parser.add_argument('-t', type=int, default=10, help='Number of threads (default: 10)')
    parser.add_argument('-m', choices=['nj', 'upgma'], default='nj', help='Tree method: nj or upgma (default: nj)')
    parser.add_argument('-v','--verbose', action='store_true', help='Print logs')
    parser.add_argument('--statistic', action='store_true', help='Output statistics file')
    parser.add_argument(
        '--statistic-file-type',
        choices=['json', 'csv'],
        default='json',
        help='Output statistics format: json or csv (default: json)'
    )
    parser.add_argument('--output-tree', help='Custom path for the output Newick tree file')
    parser.add_argument('--output-std-tree', help='Custom path for the standardized tree file')
    parser.add_argument('--leaf-order', help='Custom path for the leaf order file')
    parser.add_argument('--node-order', help='Custom path for the internal node info file')

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

def run_attotree(args):

    input_path = args.input_genomes
    output_basename = os.path.splitext(os.path.basename(input_path))[0]
    output_tree = args.output_tree or os.path.join(args.output, f"{output_basename}.nw")
    output_std_tree = args.output_std_tree or os.path.join(args.output, f"{output_basename}_std.nw")
    leaf_order = args.leaf_order or os.path.join(args.output, f"{output_basename}_leaf_order.txt")
    node_order = args.node_order or os.path.join(args.output, f"{output_basename}_node.txt")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    if args.verbose:
        print(f"[INFO] Running attotree on {input_path}...")
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
    # cmd_postprocess = [
    #     "python", "postprocess_tree.py",
    #     "--standardize", "--midpoint-outgroup", "--ladderize", "--name-internals",
    #     "-l", leaf_order,
    #     "-n", node_order,
    #     output_tree,
    #     output_std_tree
    # ]


    # subprocess.run(cmd_postprocess,capture_output=False, check=True)

    postprocesstree(output_tree, output_std_tree, True, True, True, True, leaf_order, node_order)

    # Re-add full paths to leaf_order.txt

    # Build filename full path mapping from original genome list
    path_map = {}
    with open(input_path) as f:
        for line in f:
            basename = os.path.basename(line.strip())
            genome_acession = Path(basename).name.split('.')[0]
            path_map[genome_acession] = line.strip()
    
    # Patch leaf_order.txt in-place
    with open(leaf_order) as f:
        leaves = [line.strip() for line in f]

    missing = [leaf for leaf in leaves if leaf not in path_map]
    if missing:
        print("Error: some leaves not found in input list:")
        for m in missing:
            print("  ", m)
        sys.exit(1)  # Exit with error code 1

    with open(leaf_order, 'w') as f:
        for leaf in leaves:
            f.write(path_map.get(leaf, leaf) + '\n')
    

    wall_end = time.time()
    cpu_end = os.times()
    usage = resource.getrusage(resource.RUSAGE_SELF)

    mash_time = get_duration(attotree_log.splitlines(), 'Running Mash', "Finished: 'mash triangle")
    quicktree_time = get_duration(attotree_log.splitlines(), 'Running Quicktree', "Finished: 'quicktree")

    if args.verbose:
        print(f'[INFO] Tree inference elapsed time: {round(wall_end - wall_start, 4)}s')

    stats = {
        "parameters": {
            "input": input_path,
            "output": args.output,
            "k": args.k,
            "sketch_size": args.s,
            "threads": args.t,
            "method": args.m
        },
        "timings": {
            "total": {
                "wall_time": round(wall_end - wall_start, 4),
                "user_time": round(cpu_end.user - cpu_start.user, 4),
                "system_time": round(cpu_end.system - cpu_start.system, 4)
            },
            "mash_triangle_time": mash_time,
            "quicktree_time": quicktree_time
        },
        "resources": {
            "max_rss_MB": round(usage.ru_maxrss / 1000, 2)
        }
    }

    if args.statistic:
        stats_path = os.path.join(args.output, f"tree_stats.{args.statistic_file_type}")
        if args.statistic_file_type == 'json':
            with open(stats_path, 'w') as f:
                json.dump(stats, f, indent=2)
        elif args.statistic_file_type == 'csv':
            with open(stats_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Category', 'Key', 'Value'])
                for k, v in stats['parameters'].items():
                    writer.writerow(['parameter', k, v])
                for k, v in stats['timings'].items():
                    if isinstance(v, dict):
                        for subk, subv in v.items():
                            writer.writerow(['timing', f"{k}.{subk}", subv])
                    else:
                        writer.writerow(['timing', k, v])
                for k, v in stats['resources'].items():
                    writer.writerow(['resource', k, v])
        if args.verbose:
            print(f"[INFO] Statistics saved to: {stats_path}")

def main():
    parser = argparse.ArgumentParser(
        description='Run attotree with optional parameters and time statistics'
    )
    add_tree_args(parser)
    args = parser.parse_args()
    run_attotree(args)

if __name__ == "__main__":
    main()