import argparse
import os
import sys
import random
from datetime import datetime
import time
import json 

SEED_DEFAULT = int(datetime.now().timestamp())

def get_args():
    parser = argparse.ArgumentParser(
        description='Splits a shuffled genome list into reference and remaining genome lists.'
    )
    parser.add_argument('input_genomes', help='Path to the input list of genomes')
    parser.add_argument('cut_point', help='Cut size: float < 1 for percentage, int ≥ 1 for number of genomes', type=float)
    parser.add_argument('-o', '--output', help='Output path (default: current folder)', default='.')
    parser.add_argument('--seed', type=int, help='Seed for random function (default: current timestamp)')
    parser.add_argument('-v','--verbose', action='store_true', help='Print logs')
    return parser.parse_args()

def split():
    
    args = get_args()
    seed = args.seed if args.seed is not None else SEED_DEFAULT
    input_basename = os.path.splitext(os.path.basename(args.input_genomes))[0]
    cut_arg = args.cut_point

    with open(args.input_genomes, 'r') as infile:
        genome_list = infile.readlines()

    n_genomes = len(genome_list)

    if args.verbose:
        wall_start = time.time()
        cpu_start = os.times()

    ### split operation
    random.seed(seed)
    random.shuffle(genome_list)

    if 0 < cut_arg < 1:
        cut_point = int(cut_arg * n_genomes)
        cut_mode = "fraction"
    elif cut_arg >= 1:
        cut_point = int(cut_arg)
        cut_mode = "count"
        if cut_point > n_genomes:
            print(f"Warning: requested {cut_point} genomes but only {n_genomes} available. Taking all.")
            cut_point = n_genomes
    else:
        print("Error: cut value must be > 0")
        sys.exit(1)
    reference_genomes = genome_list[:cut_point]
    remaining_genomes = genome_list[cut_point:]

    ref_path = os.path.join(args.output, f'references_{input_basename}.txt')
    rem_path = os.path.join(args.output, f'remains_{input_basename}.txt')

    with open(ref_path, 'w') as out_reference:
        for gen in reference_genomes:
            out_reference.write(gen)

    with open(rem_path, 'w') as out_remaining:
        for gen in remaining_genomes:
            out_remaining.write(gen)


    if args.verbose:
        wall_end = time.time()
        cpu_end = os.times()
        user_time = cpu_end.user - cpu_start.user
        system_time = cpu_end.system - cpu_start.system
        elapsed_time = wall_end - wall_start
        timestamp = datetime.now().isoformat()

        log_data = {
            "timestamp": timestamp,
            "seed": seed,
            "cut_mode": cut_mode,
            "input_file": os.path.basename(args.input_genomes),
            "reference_count": len(reference_genomes),
            "remaining_count": len(remaining_genomes),
            "user_time": round(user_time, 4),
            "system_time": round(system_time, 4),
            "elapsed_time": round(elapsed_time, 4)
        }

        print(json.dumps(log_data, indent=2))


if __name__ == "__main__":
    split()