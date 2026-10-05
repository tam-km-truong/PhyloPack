# Snakemake Refinement Workflow

This workflow parallelizes the windowed refinement of large genome orderings (e.g. 700k+ genomes) across local CPU cores or HPC clusters (such as Slurm).

## Architecture

1. **`split_windows` (Scatter checkpoint)**: Divides the global genome ordering into contiguous, balanced windows (default: 5,000 genomes per window) using `phylopack batch`.
2. **`process_window` (Parallel map)**: For each window, runs `attotree` to infer a local phylogenetic tree, standardizes/ladderizes the tree, and extracts leaf ordering with original full paths restored. All windows run in parallel.
3. **`merge_refined_order` (Gather)**: Concatenates all window leaf orders in exact numerical batch sequence into the final refined genome ordering.

## Configuration

Edit `workflow/config.yaml` or override options on the command line:

```yaml
# Input and output paths
input: "path/to/global_ordered_genomes.txt"
output: "results/refined_order.txt"

# Refinement window settings
window_size: 5000     # Number of genomes per local tree
min_size: 10000       # Minimum genome threshold to trigger refinement
work_dir: "results"   # Directory for intermediate batch files and trees

# Attotree inference parameters
attotree:
  k: 21
  sketch_size: 10000
  threads: 8          # CPU cores allocated per window
  method: "nj"        # "nj" or "upgma"
```

## Running the Workflow

### Local Execution (Multi-core Server)

Execute across 16 local CPU cores:

```bash
snakemake -s workflow/Snakefile --cores 16
```

Override configuration directly from the CLI:

```bash
snakemake -s workflow/Snakefile --cores 32 \
    --config input=S_enterica_global.txt output=S_enterica_refined.txt window_size=5000
```

### HPC Execution (Slurm Cluster)

#### Using Snakemake Slurm Executor (Snakemake $\ge$ 8):

```bash
snakemake -s workflow/Snakefile \
    --executor slurm \
    --jobs 50 \
    --default-resources slurm_partition=standard runtime=60 mem_mb=8000
```

#### Using Classic Cluster Submission:

```bash
snakemake -s workflow/Snakefile \
    --jobs 50 \
    --cluster "sbatch --cpus-per-task={threads} --mem=8G --time=02:00:00 -p standard"
```
