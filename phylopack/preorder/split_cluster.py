import argparse
import os
import sys
import random
from datetime import datetime
import time
import json 
import csv
import resource

SEED_DEFAULT = int(datetime.now().timestamp())

def add_split_args(parser):
    parser.add_argument('input_genomes', help='Path to the input list of genomes')
    parser.add_argument('cut-point', help='Cut size: float < 1 for percentage, int ≥ 1 for number of genomes', type=float)
    parser.add_argument('-o', '--output', help='Output path (default: current folder)', default='.')
    parser.add_argument('--seed', type=int, help='Seed for random function (default: current timestamp)')
    parser.add_argument('-v','--verbose', action='store_true', help='Print logs')
    parser.add_argument('--statistic', action='store_true', help='Output statistics file')
    parser.add_argument(
        '--statistic-file-type',
        choices=['json', 'csv'],
        default='json',
        help='Output statistics format: json or csv (default: json)'
    )
    parser.add_argument('--ref-output', help='Custom output filename for references (overrides default)')
    parser.add_argument('--rem-output', help='Custom output filename for remains (overrides default)')

def run_split(args):
    
    seed = args.seed if args.seed is not None else SEED_DEFAULT
    input_basename = os.path.splitext(os.path.basename(args.input_genomes))[0]

    with open(args.input_genomes, 'r') as infile:
        genome_list = infile.readlines()

    n_genomes = len(genome_list)

    if args.verbose:
        print('Splitting genomes')

    wall_start = time.time()
    cpu_start = os.times()

    ### split operation
    random.seed(seed)
    random.shuffle(genome_list)

    if 0 < args.cut_point < 1:
        cut_point = int(args.cut_point * n_genomes)
        cut_mode = "fraction"
    elif args.cut_point >= 1:
        cut_point = int(args.cut_point)
        cut_mode = "count"
        if cut_point > n_genomes:
            print(f"Warning: requested {cut_point} genomes but only {n_genomes} available. Taking all.")
            cut_point = n_genomes
    else:
        print("Error: cut value must be > 0")
        sys.exit(1)

    reference_genomes = genome_list[:cut_point]
    remaining_genomes = genome_list[cut_point:]

    ref_path = args.ref_output if args.ref_output else os.path.join(args.output, f'references_{input_basename}.txt')
    rem_path = args.rem_output if args.rem_output else os.path.join(args.output, f'remains_{input_basename}.txt')

    print(ref_path)
    with open(ref_path, 'w') as out_reference:
        for gen in reference_genomes:
            out_reference.write(gen)

    with open(rem_path, 'w') as out_remaining:
        for gen in remaining_genomes:
            out_remaining.write(gen)

    wall_end = time.time()
    cpu_end = os.times()
    usage = resource.getrusage(resource.RUSAGE_SELF)    

    if args.verbose:
        print(f'Splitted into {len(reference_genomes)} and {len(remaining_genomes)} lists')
        print(f'Splitting elapsed time: {round(wall_end - wall_start, 4)}s')

    stats = {
        "parameters": {
            "input": os.path.abspath(args.input_genomes),
            "output": os.path.abspath(args.output),
            "cut_point": args.cut_point,
            "cut_mode": cut_mode,
            "seed": seed,
            "reference_count": len(reference_genomes),
            "remaining_count": len(remaining_genomes)
        },
        "timings": {
            "total": {
                "wall_time": round(wall_end - wall_start, 4),
                "user_time": round(cpu_end.user - cpu_start.user, 4),
                "system_time": round(cpu_end.system - cpu_start.system, 4)
            }
        },
        "resources": {
            "max_rss_MB": round(usage.ru_maxrss / 1000, 2)
        }
    }

    if args.statistic:
        stats_path = os.path.join(args.output, f"split_stats.{args.statistic_file_type}")
        if args.statistic_file_type == 'json':
            with open(stats_path, 'w') as f:
                json.dump(stats, f, indent=2)
        else:  # CSV
            with open(stats_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Category', 'Key', 'Value'])
                for k, v in stats['parameters'].items():
                    writer.writerow(['parameter', k, v])
                for k, v in stats['timings']['total'].items():
                    writer.writerow(['timing', f"total.{k}", v])
                for k, v in stats['resources'].items():
                    writer.writerow(['resource', k, v])
        if args.verbose:
            print(f"Statistics saved to: {stats_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Splits a shuffled genome list into reference and remaining genome lists.'
    )
    add_split_args(parser)
    args = parser.parse_args()
    run_split(args)

if __name__ == "__main__":
    main()