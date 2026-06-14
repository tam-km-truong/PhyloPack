#!/usr/bin/env bash
set -euo pipefail

# Check for required arguments
if [ "$#" -lt 4 ] || [ "$#" -gt 6 ]; then
    echo "Usage: $0 <input_tsv> <id_column> <species_column> <output_dir> [min_species_size_threshold] [path_column]" >&2
    echo "Default min_species_size_threshold is 100." >&2
    exit 1
fi

INPUT_TSV="$1"
ID_COL="$2"
SPECIES_COL="$3"
OUTPUT_DIR="${4%/}"
MIN_SPECIES_SIZE_THRESHOLD="${5:-100}"
PATH_COL="${6:-}"

# Verify input file exists
if [ ! -f "$INPUT_TSV" ]; then
    echo "Error: Input file '$INPUT_TSV' not found." >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

awk -F'\t' \
    -v id_name="$ID_COL" \
    -v sp_name="$SPECIES_COL" \
    -v path_name="$PATH_COL" \
    -v out_dir="$OUTPUT_DIR/" \
    -v threshold="$MIN_SPECIES_SIZE_THRESHOLD" \
'
NR == 1 {
    # Dynamically locate columns
    for (i = 1; i <= NF; i++) {
        if ($i == id_name) id_col = i
        if ($i == sp_name) sp_col = i
        if (path_name != "" && $i == path_name) path_col = i
    }
    
    if (!id_col || !sp_col) {
        print "Error: Required ID or Species columns not found in header." > "/dev/stderr"
        exit 1
    }
    
    if (path_name != "" && !path_col) {
        print "Error: Specified path column \047" path_name "\047 not found in header." > "/dev/stderr"
        exit 1
    }
    next
}
{
    id = $id_col; sp = $sp_col
    if (id == "" || sp == "") next
    
    prefix = (path_col ? $path_col : "")
    
    gsub(/[[:space:]]+/, "_", sp)
    
    full_path = (prefix != "" ? prefix "/" id : id)
    
    idx = ++count[sp]
    data[sp, idx] = full_path
}
END {
    remainder_file = out_dir "subthreshold_remainder.txt"
    
    for (sp in count) {
        if (count[sp] < threshold) {
            for (j = 1; j <= count[sp]; j++) {
                print data[sp, j] > remainder_file
            }
        } else {
            sp_file = out_dir sp ".txt"
            for (j = 1; j <= count[sp]; j++) {
                print data[sp, j] > sp_file
            }
        }
    }
}' "$INPUT_TSV"