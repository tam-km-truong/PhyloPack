#!/usr/bin/env python3
"""
Generate synthetic multi-species mock genomes for testing, debugging, and benchmarking PhyloPack.

Generates reproducible multi-species FASTA files with controllable phylogenetic divergence
without committing large files to Git. Supports unequal genome distribution across species
while preserving the exact total genome count. Accession IDs remain neutral (without species names),
and an accession-to-species mapping is recorded in a separate TSV metadata file.

Run it: python tests/generate_mock_data.py -u -o tests/data/synthetic -s 10 -n 50 -t accession_species.tsv
"""

import argparse
import csv
import gzip
import os
import random
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional


class MockDatasetResult(NamedTuple):
    manifest_path: Path
    tsv_path: Path
    species_counts: Dict[str, int]


def generate_random_sequence(length: int, gc_content: float = 0.5, rng: Optional[random.Random] = None) -> str:
    """Generate a random DNA sequence with a target GC content."""
    if rng is None:
        rng = random.Random()
    gc_prob = gc_content / 2.0
    at_prob = (1.0 - gc_content) / 2.0
    bases = ['A', 'C', 'G', 'T']
    weights = [at_prob, gc_prob, gc_prob, at_prob]
    return ''.join(rng.choices(bases, weights=weights, k=length))


def mutate_sequence(seq: str, mutation_rate: float, rng: Optional[random.Random] = None) -> str:
    """Introduce random point substitutions into a DNA sequence at a given rate."""
    if rng is None:
        rng = random.Random()
    if mutation_rate <= 0.0:
        return seq

    seq_list = list(seq)
    n = len(seq_list)
    num_muts = int(n * mutation_rate)
    if num_muts == 0 and mutation_rate > 0:
        num_muts = 1

    all_bases = ['A', 'C', 'G', 'T']
    sampled_indices = rng.sample(range(n), min(num_muts, n))
    for idx in sampled_indices:
        current_base = seq_list[idx]
        alternatives = [b for b in all_bases if b != current_base]
        seq_list[idx] = rng.choice(alternatives)

    return ''.join(seq_list)


def split_into_contigs(seq: str, num_contigs: int, rng: Optional[random.Random] = None) -> List[str]:
    """Split a sequence into approximately equal contigs."""
    if num_contigs <= 1 or len(seq) < num_contigs:
        return [seq]
    if rng is None:
        rng = random.Random()

    n = len(seq)
    cut_points = sorted(rng.sample(range(100, n - 100), num_contigs - 1))
    cut_points = [0] + cut_points + [n]
    return [seq[cut_points[i]:cut_points[i + 1]] for i in range(num_contigs)]


def partition_species_counts(
    num_species: int,
    total_genomes: int,
    unequal: bool = False,
    min_count: int = 1,
    rng: Optional[random.Random] = None,
) -> List[int]:
    """
    Distribute total_genomes across num_species.
    Guarantees sum(counts) == total_genomes and every count >= min_count.
    """
    if rng is None:
        rng = random.Random()

    min_required = num_species * min_count
    if total_genomes < min_required:
        raise ValueError(
            f"Total genomes ({total_genomes}) must be >= num_species * min_count "
            f"({num_species} * {min_count} = {min_required})"
        )

    if not unequal:
        base = total_genomes // num_species
        remainder = total_genomes % num_species
        return [base + (1 if i < remainder else 0) for i in range(num_species)]

    # Unequal distribution using Dirichlet-like random weights
    weights = [rng.gammavariate(1.5, 1.0) for _ in range(num_species)]
    total_weight = sum(weights)
    remainder = total_genomes - min_required

    floors = []
    fractions = []
    for i, w in enumerate(weights):
        alloc = (w / total_weight) * remainder
        fl = int(alloc)
        floors.append(fl)
        fractions.append((alloc - fl, i))

    leftover = remainder - sum(floors)
    fractions.sort(reverse=True)
    counts = [min_count + fl for fl in floors]
    for k in range(leftover):
        idx = fractions[k][1]
        counts[idx] += 1

    return counts


def generate_mock_dataset(
    output_dir: Path,
    num_species: int = 4,
    strains_per_species: int = 10,
    total_genomes: Optional[int] = None,
    species_counts: Optional[List[int]] = None,
    unequal_species: bool = False,
    min_strains_per_species: int = 1,
    genome_length: int = 30000,
    species_divergence: float = 0.04,
    strain_divergence: float = 0.002,
    contigs_per_genome: int = 1,
    gc_content: float = 0.5,
    seed: int = 42,
    accession_prefix: str = "SYNTH",
    manifest_name: str = "genomes.txt",
    tsv_name: str = "accession_species.tsv",
    relative_manifest: bool = True,
) -> MockDatasetResult:
    """
    Generate synthetic multi-species mock genomes and write:
      1. FASTA files (.fa.gz) with neutral accession IDs.
      2. Manifest file listing all genome paths.
      3. TSV metadata file mapping accession to species ID and genome path.

    Parameters:
        output_dir: Target directory to save genomes, manifest, and TSV.
        num_species: Number of distinct species clades.
        strains_per_species: Number of strains generated per species if equal.
        total_genomes: Total genomes across all species (defaults to num_species * strains_per_species).
        species_counts: Explicit list of strain counts per species.
        unequal_species: If True, distributes total_genomes unequally across species.
        min_strains_per_species: Minimum strains per species in unequal mode.
        genome_length: Total sequence length per genome (bp).
        species_divergence: Mutation rate separating each species from common ancestor.
        strain_divergence: Mutation rate separating strains within the same species.
        contigs_per_genome: Number of contigs per genome FASTA.
        gc_content: Target GC ratio (0.0 to 1.0).
        seed: Random seed for deterministic reproducibility.
        accession_prefix: Neutral prefix for accession IDs (e.g. 'ACC' -> 'ACC00001').
        manifest_name: Name of the manifest file listing paths.
        tsv_name: Name of TSV file tracking accession and species.
        relative_manifest: If True, manifest contains paths relative to CWD.

    Returns:
        MockDatasetResult(manifest_path, tsv_path, species_counts)
    """
    rng = random.Random(seed)
    output_dir = Path(output_dir).resolve()
    fasta_dir = output_dir / "fasta_files"
    fasta_dir.mkdir(parents=True, exist_ok=True)

    # Determine counts per species
    if species_counts is not None:
        counts = [int(c) for c in species_counts]
        num_species = len(counts)
        total_genomes = sum(counts)
    else:
        if total_genomes is None:
            total_genomes = num_species * strains_per_species
        counts = partition_species_counts(
            num_species=num_species,
            total_genomes=total_genomes,
            unequal=unequal_species,
            min_count=min_strains_per_species,
            rng=rng,
        )

    # 1. Root ancestral genome
    ancestral_root = generate_random_sequence(genome_length, gc_content=gc_content, rng=rng)

    manifest_lines = []
    tsv_records = []
    species_summary: Dict[str, int] = {}
    global_idx = 1

    # 2. Generate species roots and strains
    for sp_idx, count in enumerate(counts, 1):
        species_id = f"species_{sp_idx:02d}"
        species_summary[species_id] = count
        species_root = mutate_sequence(ancestral_root, species_divergence, rng=rng)

        for _ in range(count):
            strain_seq = mutate_sequence(species_root, strain_divergence, rng=rng)
            contigs = split_into_contigs(strain_seq, contigs_per_genome, rng=rng)

            # Neutral accession without species name
            accession = f"{accession_prefix}{global_idx:05d}"
            global_idx += 1

            filename = f"{accession}.fa.gz"
            filepath = fasta_dir / filename

            with gzip.open(filepath, "wt") as f:
                for c_idx, c_seq in enumerate(contigs, 1):
                    f.write(f">{accession}_contig{c_idx:03d} len={len(c_seq)}\n")
                    for i in range(0, len(c_seq), 80):
                        f.write(c_seq[i:i + 80] + "\n")

            if relative_manifest:
                try:
                    rel_path = filepath.relative_to(Path.cwd())
                    formatted_path = f"./{rel_path}"
                except ValueError:
                    formatted_path = str(filepath)
            else:
                formatted_path = str(filepath)

            manifest_lines.append(formatted_path)
            tsv_records.append({
                "accession": accession,
                "species_id": species_id,
                "genome_path": formatted_path
            })

    # 3. Write manifest file
    manifest_path = output_dir / manifest_name
    with open(manifest_path, "w") as f:
        for line in manifest_lines:
            f.write(line + "\n")

    # 4. Write TSV tracking accession and species ID
    tsv_path = output_dir / tsv_name
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["accession", "species_id", "genome_path"], delimiter="\t")
        writer.writeheader()
        writer.writerows(tsv_records)

    return MockDatasetResult(
        manifest_path=manifest_path,
        tsv_path=tsv_path,
        species_counts=species_summary
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate synthetic multi-species mock genomes for PhyloPack testing."
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("tests/data/synthetic"),
        help="Directory to store generated genomes, manifest, and TSV (default: tests/data/synthetic)"
    )
    parser.add_argument(
        "-s", "--num-species",
        type=int,
        default=4,
        help="Number of species clades (default: 4)"
    )
    parser.add_argument(
        "-n", "--strains-per-species",
        type=int,
        default=10,
        help="Number of strains per species if equal distribution (default: 10)"
    )
    parser.add_argument(
        "-T", "--total-genomes",
        type=int,
        default=None,
        help="Total genomes across all species (defaults to num_species * strains_per_species)"
    )
    parser.add_argument(
        "-u", "--unequal",
        action="store_true",
        help="Distribute total genomes unequally across species while preserving the exact total"
    )
    parser.add_argument(
        "--species-counts",
        type=int,
        nargs="+",
        default=None,
        help="Explicit strain counts for each species (e.g. --species-counts 20 12 5 3)"
    )
    parser.add_argument(
        "--min-strains",
        type=int,
        default=1,
        help="Minimum strains per species when using --unequal (default: 1)"
    )
    parser.add_argument(
        "-l", "--genome-length",
        type=int,
        default=30000,
        help="Genome sequence length in bp (default: 30000)"
    )
    parser.add_argument(
        "--accession-prefix",
        type=str,
        default="SYNTH",
        help="Neutral prefix for accession IDs (default: ACC)"
    )
    parser.add_argument(
        "--species-divergence",
        type=float,
        default=0.04,
        help="Inter-species divergence mutation rate (default: 0.04)"
    )
    parser.add_argument(
        "--strain-divergence",
        type=float,
        default=0.002,
        help="Intra-species strain divergence mutation rate (default: 0.002)"
    )
    parser.add_argument(
        "-c", "--contigs",
        type=int,
        default=1,
        help="Number of contigs per genome (default: 1)"
    )
    parser.add_argument(
        "--gc",
        type=float,
        default=0.5,
        help="Target GC content ratio (default: 0.5)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for deterministic generation (default: 42)"
    )
    parser.add_argument(
        "-m", "--manifest-name",
        type=str,
        default="genomes.txt",
        help="Manifest filename (default: genomes.txt)"
    )
    parser.add_argument(
        "-t", "--tsv-name",
        type=str,
        default="accession_species.tsv",
        help="TSV filename tracking accession and species ID (default: accession_species.tsv)"
    )
    parser.add_argument(
        "--absolute-paths",
        action="store_true",
        help="Store absolute paths in manifest instead of relative paths"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    result = generate_mock_dataset(
        output_dir=args.output_dir,
        num_species=args.num_species,
        strains_per_species=args.strains_per_species,
        total_genomes=args.total_genomes,
        species_counts=args.species_counts,
        unequal_species=args.unequal,
        min_strains_per_species=args.min_strains,
        genome_length=args.genome_length,
        species_divergence=args.species_divergence,
        strain_divergence=args.strain_divergence,
        contigs_per_genome=args.contigs,
        gc_content=args.gc,
        seed=args.seed,
        accession_prefix=args.accession_prefix,
        manifest_name=args.manifest_name,
        tsv_name=args.tsv_name,
        relative_manifest=not args.absolute_paths,
    )
    total = sum(result.species_counts.values())
    print(f"[INFO] Successfully generated {total} mock genomes across {len(result.species_counts)} species:")
    for sp_name, count in result.species_counts.items():
        print(f"  - {sp_name}: {count} genomes")
    print(f"[INFO] FASTA directory: {args.output_dir / 'fasta_files'}")
    print(f"[INFO] Manifest written to: {result.manifest_path}")
    print(f"[INFO] Metadata TSV written to: {result.tsv_path}")


if __name__ == "__main__":
    main()
