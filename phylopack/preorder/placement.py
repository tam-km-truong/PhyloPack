import argparse
import os
import time
import json
import subprocess
from datetime import datetime
from collections import defaultdict
import resource
import csv

SEED_DEFAULT = int(datetime.now().timestamp())

def add_placement_args(parser):
    parser.add_argument('genomes_list_1', help='Path to query genomes (to be placed)')
    parser.add_argument('genomes_list_2', help='Path to reference genomes (anchors)')
    parser.add_argument('-o', '--output', help='Output folder (default: current folder)', default='.')
    parser.add_argument('-v','--verbose', action='store_true', help='Print logs')
    parser.add_argument('-k', type=int, default=21, help='K-mer size (default: 21)')
    parser.add_argument('-s', type=int, default=1000, help='Sketch size (default: 1000)')
    parser.add_argument('-t', type=int, default=10, help='Number of threads (default: 10)')
    parser.add_argument('--statistic', action='store_true', help='Output json statistics file')
    parser.add_argument(
        '--statistic-file-type',
        choices=['json', 'csv'],
        default='json',
        help='Output statistics format: json or csv (default: json)'
    )


def mash_sketch(genomes_list, output, k, s, t, verbose = False):

    if verbose:
        print(f"Sketching genomes from: {genomes_list}")

    basename = os.path.splitext(os.path.basename(genomes_list))[0]
    sketch = os.path.join(output, f"{basename}")


    start = time.time()
    cpu_start = os.times()

    cmd = [
        'mash',
        'sketch',
        '-l', genomes_list,
        '-o', sketch,
        '-k', str(k),
        '-s', str(s),
        '-p', str(t),
    ]    

    print("Running mash sketch with:")
    print(" ".join(cmd))

    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, text=True, check=True)

    end = time.time()
    cpu_end = os.times()

    return sketch, {
        'wall_time': end - start,
        'user_time': cpu_end.user - cpu_start.user,
        'system_time': cpu_end.system - cpu_start.system
    }

def mash_distance(sketch_1, sketch_2, t, output_path, verbose=False):
    
    if verbose:
        print(f"Calculating Mash distance between {sketch_1} and {sketch_2}")

    start = time.time()
    cpu_start = os.times()

    cmd = [
        'mash', 'dist',
        f'{sketch_2}.msh', # reference first
        f'{sketch_1}.msh', # remaining
        '-p', str(t),
        '-t'
    ]    

    distance_file = os.path.join(output_path, f"{os.path.basename(sketch_1)}_{os.path.basename(sketch_2)}_distance.tsv")

    with open(distance_file, 'w') as f:
        subprocess.run(cmd, stdout=f, text=True, check=True)

    end = time.time()
    cpu_end = os.times()

    return distance_file, {
        'wall_time': end - start,
        'user_time': cpu_end.user - cpu_start.user,
        'system_time': cpu_end.system - cpu_start.system
    }

def argmin(distance_file, rows, cols, verbose = False):

    if verbose:
        print(f"Finding argmin from distance file: {distance_file}")

    start = time.time()
    cpu_start = os.times()

    awk_script = (
        'BEGIN { FS="\\t" } '
        'NR > 1 { '
        'min = $2; idx = 0; '
        'for (i = 3; i <= NF; i++) { '
        'if ($i < min) { min = $i; idx = i - 2; } '
        '} print idx; }'
    )

    result = subprocess.run(
        ["awk", awk_script],
        stdin=open(distance_file),
        capture_output=True,
        text=True,
        check=True
    )

    argmin_indices = [int(i) for i in result.stdout.strip().split("\n")]

    groups = defaultdict(list)

    for row_name, idx in zip(rows, argmin_indices):
        col_name = cols[idx]
        groups[col_name].append(row_name)

    end = time.time()
    cpu_end = os.times()

    return groups, {
        'wall_time': end - start,
        'user_time': cpu_end.user - cpu_start.user,
        'system_time': cpu_end.system - cpu_start.system
    }

def run_placement(args):

    full_start = time.time()
    full_cpu_start = os.times()

    ### Sketching

    sketch_1, sketch_time_1 = mash_sketch(args.genomes_list_1, args.output, args.k, args.s, args.t, args.verbose)
    sketch_2, sketch_time_2 = mash_sketch(args.genomes_list_2, args.output, args.k, args.s, args.t, args.verbose)

    ### Calculating distances

    distance_file, dis_time = mash_distance(sketch_1, sketch_2, args.t, args.output, args.verbose)

    ### Grouping and writing to output

    with open(args.genomes_list_1) as f:
        row_names = [os.path.basename(line.strip()).split('.')[0] for line in f]    
    
    with open(args.genomes_list_2) as f:
        col_names = [os.path.basename(line.strip()).split('.')[0] for line in f]    

    argmin_result, grouping_time = argmin(distance_file, row_names, col_names, args.verbose)

    preorder_file = os.path.join(args.output, f"placement_order.txt")
    # Preorder writing
    if args.verbose:
        print(f"Writing preorder result to file")
    with open(preorder_file, 'w') as preorder_file:
        for col in col_names:
            preorder_file.write(col + '\n')
            for val in argmin_result[col]:
                preorder_file.write(val + '\n')

    full_end = time.time()
    full_cpu_end = os.times()

    usage = resource.getrusage(resource.RUSAGE_SELF)

    stats = {
        "parameters": {
            "genomes_list_1": args.genomes_list_1,
            "genomes_list_2": args.genomes_list_2,
            "output": args.output,
            "k": args.k,
            "sketch_size": args.s,
            "threads": args.t,
        },
        "timings": {},
        "resources": {"max_rss_MB": round(usage.ru_maxrss / 1000, 2)}
    }

    stats["timings"]['sketch_list_1'] = sketch_time_1
    stats["timings"]['sketch_list_2'] = sketch_time_2
    stats["timings"]['mash_distance'] = dis_time
    stats["timings"]['grouping'] = grouping_time
    stats["timings"]['total'] = {
        'wall_time': full_end - full_start,
        'user_time': full_cpu_end.user - full_cpu_start.user,
        'system_time': full_cpu_end.system - full_cpu_start.system
    }

    if args.verbose:
        print(f'Placement elapsed time: {round(full_end - full_start, 4)}s')

    if args.statistic:
        stats_path = os.path.join(args.output, f"placement_stats.{args.statistic_file_type}")
        if args.statistic_file_type == 'json':

            with open(stats_path, 'w') as f:
                json.dump(stats, f, indent=2)
        elif args.statistic_file_type == 'csv':
            with open(stats_path, 'w', newline='') as f:
                writer = csv.writer(f)

                writer.writerow(["Category", "Key", "Value"])

                for k, v in stats["parameters"].items():
                    writer.writerow(["parameter", k, v])

                for k, v in stats["timings"].items():
                    for sub_k, sub_v in v.items():
                        writer.writerow(["timing", f"{k}.{sub_k}", sub_v])

                for k, v in stats["resources"].items():
                    writer.writerow(["resource", k, v])      
        if args.verbose:
            print(f"Statistics saved to: {stats_path}")      

def main():
    parser = argparse.ArgumentParser(
        description='Calculate distance from 2 lists of genomes using mash, list 1 is query, list 2 is reference'
    )
    add_placement_args(parser)
    args = parser.parse_args()
    run_placement(args)

if __name__ == "__main__":
    main()