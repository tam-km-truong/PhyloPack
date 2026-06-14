# PhyloPack

**PhyloPack** is a protocol to facilitate the handling of million-genomes scale bacterial collections.

It first computes a global genome ordering that places similar genomes close together, then partitions this ordering into batches suitable for downstream applications such as compression, exact k-mer indexing, or distributed processing. The approach uses a skeleton phylogeny built from a subset of reference genomes and subsequently places remaining genomes onto this backbone.

### Key results:

Using the ATB 2024-08 release (approx. 2.4 million genomes):

- **< 15 GB** using MBGC [cite here] for major pathogenic species (~93% of the collection).
- **~33 GB** using AGC [cite here] while retaining fast random-access capabilities.

## Installation

```bash
git clone https://github.com/tam-km-truong/phylopack.git
cd phylopack
pip install .
```

### Requirements

* `Python >= 3.8`
* `attotree`

attotree can be installed through Bioconda:

```bash
conda install -c bioconda -c conda-forge attotree
```

### Input format

Input files should contain one genome path per line:

```
~/data/genome1.fa
~/data/genome2.fa
~/data/genome3.fa
```

### Workflow

1. Compute a global genome ordering using `phylopack preorder`.
2. Split the ordered list into batches using `phylopack batch`.
3. Use the generated ordering or batches in downstream applications.

### Quick Example - Computing order

```bash
phylopack preorder tests/data/genomes.txt --cut-point 0.2 -o ./debug/out.txt -v
```
Output: a text file containing the ordered genome list, one genome path per line.

### Quick Example - Batching

```bash
phylopack batch tests/data/genomes.txt -o debug/ -v
```

Output: multiple batch files containing contiguous segments of the ordered genome list.

## Main command:

### Genome ordering

PhyloPack first splits the input genomes into a skeleton set and a placement set. The skeleton genomes are used to infer a skeleton phylogeny, and the remaining genomes are then placed onto this backbone.

```bash
phylopack preorder [OPTIONS] input_genomes
```

Commonly adjusted options:

```text
-o OUTPUT, --output OUTPUT
    Output file for the final genome ordering.

-c CUT_POINT, --cut-point CUT_POINT
    Controls the size of the skeleton set. Values between 0 and 1 are interpreted as a fraction of the input genomes, while values greater than 1 are interpreted as the number of skeleton genomes.

-m {nj,upgma} [nj]

--splitting-scheme {random,nth-accession,custom}
    Strategy used to select the skeleton genomes (default: random).

--custom-ref CUSTOM_REF
    Path to a custom skeleton genome list, using the same format as the input file.

--exclude-skeleton
    Exclude the skeleton genomes from the final global ordering.
```

Advanced options:

```text
-k K
-s-reference S_REFERENCE
-s-placement S_PLACEMENT
--max-skeleton-genomes MAX_SKELETON_GENOMES
--min-skeleton-genomes MIN_SKELETON_GENOMES
--seed SEED
--statistic
--statistic-file-type {json,csv}
--debug
[...]
```

For the complete list of options:

```bash
phylopack preorder --help
```

### Batching

Batching partitions the ordered genomes into contiguous, balanced batches according to a user-specified maximum batch size.

```bash
phylopack batch [OPTIONS] input
```

Commonly adjusted options:

```text
-o OUTPUT_DIR, --output-dir OUTPUT_DIR
    Directory where batch files are written.

-s TARGET_SIZE, --target-size TARGET_SIZE
    Maximum number of genomes per batch. PhyloPack creates balanced batches while keeping each batch below this limit.

-n BATCH_NAME, --batch-name BATCH_NAME
    Prefix used for generated batch files.
```

Advanced options:

```text
-d DIGITS
[...]
```

For the complete list of options:

```bash
phylopack batch --help
```

## Real-world application: Compressing the AllTheBacteria collection

On the AllTheBacteria (ATB) collection [cite here], containing millions of genomes, PhyloPack was used to generate phylogeny-aware genome orderings and batches prior to compression.

Dataset scale:

- ~2.7 million genomes
- [number] species 
- Largest cluster (*S. enterica*): ~700k genomes

### 1. Species clustering

First, download the ATB metadata table [link here].

Genomes are grouped by species using the helper script in `suppl_scripts/species_clustering`:

```bash
./suppl_scripts/species_clustering.sh \
    <input_tsv> \
    <id_column> \
    <species_column> \
    <output_dir> \
    [min_species_size_threshold] \
    [path_column]
```

The default value of `min_species_size_threshold` is 100.

Example:

```bash
./suppl_scripts/species_clustering.sh \
    file_list_202505.tsv \
    sample \
    sylp_species \
    ~/ATB/clusters/ \
    100 \
    tar_xz
```

This step produces one genome list per species cluster.

### 2. Compute a global ordering for each species

For each species cluster, compute a PhyloPack global ordering.

Example for *S. enterica*:

```bash
phylopack preorder \
    S_enterica_cluster.txt \
    --cut-point 0.05 \
    -t 96 \
    --max-skeleton-genomes 20000 \
    --min-skeleton-genomes 10000 \
    -o S_enterica_cluster_global_ordered.txt \
    -v \
    --statistic \
    --debug
```

In this example, the skeleton set size is controlled by three parameters:

* `--cut-point 0.05` selects 5% of genomes as the skeleton set.
* `--min-skeleton-genomes 10000` ensures that clusters smaller than 10k genomes use all genomes as skeleton genomes, resulting in a full phylogenetic inference.
* `--max-skeleton-genomes 20000` caps the skeleton set size for very large clusters to reduce runtime and memory requirements.

For example, the *S. enterica* cluster contains approximately 700k genomes. With a cut point of 5%, the skeleton set would contain roughly 35k genomes. Because this exceeds the maximum threshold, the skeleton set is capped at 20k genomes.

### 3. Refine the ordering within large clusters

For clusters requiring genome placement, an additional refinement step is performed after the initial global ordering. This improves local ordering accuracy within placement-heavy clusters without requiring phylogenetic inference on the entire dataset.


The ordered genome list is divided into non-overlapping windows of 5k genomes. 

This keeps each local phylogenetic inference trivial while preserving the global-scale ordering produced by PhyloPack.

For each window, a local phylogenetic tree is inferred using Attotree:

```bash
attotree \
    -L ~/ATB/clusters/cluster_window_1.txt \
    -o ~/ATB/tree/cluster_window_1.nw
```

The tree is then converted into a genome ordering:

```bash
python phylopack/preorder/postprocess_tree.py \
    --standardize \
    --midpoint-outgroup \
    --ladderize \
    --name-internals \
    -l ~/leaves/leaves_order.txt \
    ~/tree/tree.nw -
```

The refined window orderings are concatenated to produce the final ordering:

```bash
for cluster in cluster_window*; do
    cat "$cluster" >> refined_order/cluster_refined_order.txt
done
```

### 4. Create balanced batches

The refined ordering is partitioned into balanced batches:

```bash
phylopack batch \
    -o ~/ATB/final_batches/ \
    --target-size 5000 \
    --batch-name cluster_ \
    -v \
    cluster_refined_order.txt
```

This step produces contiguous batches while preserving the locality established by the ordering process.

### 5. Compress each batch with AGC

Each batch is compressed independently using AGC. The first genome in the batch is used as the reference genome.

```bash
for f in ~/ATB/final_batches/*; do
    ref=$(head -n 1 "$f")

    agc create \
        -a \
        -b 500 \
        -s 1500 \
        -t 25 \
        -v 2 \
        -o ~/ATB/results/agc/$(basename "$f" .txt).agc \
        -i "$f" \
        "$ref"
done
```

The resulting AGC archives can then be used for downstream storage and retrieval of bacterial genome collections.
