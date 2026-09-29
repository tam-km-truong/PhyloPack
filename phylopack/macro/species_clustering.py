import argparse
import os
import sys
import heapq
import tempfile
from collections import defaultdict
from collections.abc import Iterator
from contextlib import ExitStack
from typing import Dict, List, Tuple, Optional

def add_species_clustering_parser(parser):
    # Input / Output
    parser.add_argument(
        "metadata",
        help="Path to TSV/CSV metadata file containing accessions, species IDs, and file paths."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output path for the ordered genome file paths (one path per line)."
    )
    parser.add_argument(
        "-d", "--delimiter",
        default="\t",
        help="Metadata field delimiter."
    )
    parser.add_argument(
        "--separators-output",
        dest="separators_output",
        help="Optional TSV output for cluster boundaries: species_id, start_idx, count, representative_path."
    )
    parser.add_argument(
        "--species-order-output",
        help="Optional text file recording the resulting species order (one species_id per line)."
    )

    # Column Mapping
    parser.add_argument(
        "--accession-col",
        default="accession",
        help="Header name or 0-based column index for accession IDs (default: 'accession')."
    )
    parser.add_argument(
        "--species-col",
        default="species_id",
        help="Header name or 0-based column index for species identifiers (default: 'species_id')."
    )
    parser.add_argument(
        "--path-col",
        default=None,
        help="Header name or 0-based column index for genome file paths (default: 'None')."
    )
    parser.add_argument(
        "--path-prefix",
        default="",
        help="Directory prefix prepended to relative genome paths (default: '')."
    )

    # Ordering Strategy
    parser.add_argument(
        "--sort-species-by",
        choices=["desc", "asc"],
        default="desc",
        help="Inter-species ordering: 'desc' (largest species first), 'asc' (smallest first)."
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose progress logging."
    )

# Helper functions
def resolve_columns(header_line: str, delim: str, acc_col: str, sp_col: str, path_col: Optional[str] = None) -> Tuple[int, int, Optional[int]]:    
    fields = [f.strip() for f in header_line.rstrip('\r\n').split(delim)]

    acc_idx = int(acc_col) if acc_col.isdigit() else fields.index(acc_col)
    spc_idx = int(sp_col) if sp_col.isdigit() else fields.index(sp_col)
    if not path_col:
        path_idx = None
    else:    path_idx = int(path_col) if path_col.isdigit() else fields.index(path_col)

    return (acc_idx, spc_idx, path_idx)

def count_and_rank_species(metadata_path: str, delim: str, sp_idx: int, sort_strategy: str) -> Tuple[Dict[str, int], Dict[str, int], List[str]]:
    counts = defaultdict(int)
    max_split = sp_idx + 1

    with open(metadata_path, "rt", buffering=1024 * 1024) as f:
        _ = f.readline()  # skip header
        for line in f:
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            parts = line.split(delim, max_split)
            counts[parts[sp_idx]] += 1

    if sort_strategy == "desc":
        species_order = sorted(counts.keys(), key=lambda sp: (-counts[sp], sp))
    elif sort_strategy == "asc":
        species_order = sorted(counts.keys(), key=lambda sp: (counts[sp], sp))
    else:
        print(f"Error: Unsupported sort strategy '{sort_strategy}'", file=sys.stderr)
        sys.exit(1)

    species_to_rank = {sp: rank for rank, sp in enumerate(species_order)}
    return counts, species_to_rank, species_order

def create_sorted_chunks(
    metadata_path: str,
    delim: str,
    col_indices: Tuple[int, int, Optional[int]],
    species_to_rank: Dict[str, int],
    temp_dir: str,
    chunk_size: int = 1_000_000,
    path_prefix: str = ""
) -> List[str]:

    max_split = max([idx for idx in col_indices if idx is not None]) + 1
    chunk_files = []
    chunk_buffer = []

    def flush_chunk(buffer, chunk_id):
        # Sort in memory by (species_rank, accession)
        buffer.sort(key=lambda row: (row[0], row[1]))
        chunk_path = os.path.join(temp_dir, f"chunk_{chunk_id:04d}.tsv")
        with open(chunk_path, "wt", buffering=256 * 1024) as out:
            for rank, acc, full_path, sp_id in buffer:
                out.write(f"{rank}\t{sp_id}\t{acc}\t{full_path}\n")
        return chunk_path
        
    acc_idx, sp_idx, path_idx = col_indices
    with open(metadata_path, "rt", buffering=1024 * 1024) as f:
        _=f.readline()
        for line in f:
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            parts = line.split(delim, max_split)
            accession = parts[acc_idx]
            species = parts[sp_idx]
            raw_path = "" if path_idx is None else parts[path_idx]
            if path_prefix and raw_path:
                full_path = os.path.join(path_prefix, raw_path)
            else:
                full_path = raw_path if raw_path else accession

            rank = species_to_rank[species]

            chunk_buffer.append((rank, accession, full_path, species))

            if len(chunk_buffer) >= chunk_size:
                chunk_files.append(flush_chunk(chunk_buffer, len(chunk_files)))
                chunk_buffer.clear()

    if chunk_buffer:
        chunk_files.append(flush_chunk(chunk_buffer, len(chunk_files)))
        chunk_buffer.clear()

    return chunk_files

def chunk_line_generator(file_obj) -> Iterator[Tuple[int, str, str, str]]:
    for line in file_obj:
        rank_str, sp_id, acc, path = line.rstrip("\r\n").split("\t", 3)
        yield (int(rank_str), acc, path, sp_id)


def merge_and_emit(
    chunk_paths: List[str],
    output_path: str,
    species_counts: Dict[str, int],
    boundaries_path: Optional[str] = None
) -> List[Tuple[str, int, int, str]]:
    separators = []
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output_path, "wt", buffering=1024 * 1024) as out_f, ExitStack() as stack:
        chunk_files = [
            stack.enter_context(open(cp, "rt", buffering=64 * 1024))
            for cp in chunk_paths
        ]
        streams = [chunk_line_generator(f) for f in chunk_files]
        merged_stream = heapq.merge(*streams, key=lambda row: (row[0], row[1]))
        current_sp = None
        current_start_idx = 0
        current_sp_count = 0
        global_idx = 0
        median_target_idx = -1
        current_rep_path = ""
        for rank, acc, path, sp_id in merged_stream:
            if sp_id != current_sp:
                if current_sp is not None:
                    separators.append((current_sp, current_start_idx, current_sp_count, current_rep_path))
                current_sp = sp_id
                current_start_idx = global_idx
                current_sp_count = 0
                expected_total = species_counts[sp_id]
                median_target_idx = current_start_idx + (expected_total // 2)
            if global_idx == median_target_idx:
                current_rep_path = path
            out_f.write(path + "\n")
            current_sp_count += 1
            global_idx += 1
        if current_sp is not None:
            separators.append((current_sp, current_start_idx, current_sp_count, current_rep_path))
    if boundaries_path:
        b_dir = os.path.dirname(boundaries_path)
        if b_dir:
            os.makedirs(b_dir, exist_ok=True)
        with open(boundaries_path, "wt") as bf:
            bf.write("species_id\tstart_idx\tcount\trepresentative\n")
            for sp_id, s_idx, cnt, rep in separators:
                bf.write(f"{sp_id}\t{s_idx}\t{cnt}\t{rep}\n")
    return separators

def normalize_delimiter(delim: str) -> str:
    """Handles literal escape strings like '\\t' passed from shell."""
    if delim in ("\\t", "tab"):
        return "\t"
    if delim in ("\\s", "space"):
        return " "
    return delim

def run_species_clustering(args):
    delim = normalize_delimiter(args.delimiter)
    with open(args.metadata, "rt") as f:
        header = f.readline()
        if not header:
            raise ValueError(f"Metadata file '{args.metadata}' is empty.")
    col_indices = resolve_columns(
        header, delim, args.accession_col, args.species_col, args.path_col
    )
    species_counts, species_to_rank, species_order = count_and_rank_species(
        args.metadata, delim, col_indices[1], args.sort_species_by
    )
    with tempfile.TemporaryDirectory(prefix="phylopack_sp_") as temp_dir:
        chunk_files = create_sorted_chunks(
            args.metadata, delim, col_indices, species_to_rank,
            temp_dir=temp_dir, path_prefix=args.path_prefix
        )
        separators = merge_and_emit(
            chunk_files, args.output, species_counts, args.separators_output
        )
    if args.species_order_output:
        so_dir = os.path.dirname(args.species_order_output)
        if so_dir:
            os.makedirs(so_dir, exist_ok=True)
        with open(args.species_order_output, "wt") as sof:
            for sp in species_order:
                sof.write(f"{sp}\n")
    return separators, species_order

def main():
    parser = argparse.ArgumentParser(
        description='The description'
    )
    add_species_clustering_parser(parser)
    args = parser.parse_args()

    run_species_clustering(args)

if __name__ == "__main__":
    main()