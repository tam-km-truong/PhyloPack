# PhyloPack

**PhyloPack** is a protocol to facilitate the handling of million-genomes scale bacterial collections by taking advantage of phylogenetic relationships between genomes [1].

Its primary output is a global genome ordering that places similar genomes close together, then this ordering could be partitioned into batches suitable for downstream applications such as compression, exact k-mer indexing, or distributed processing. The approach uses a skeleton phylogeny built from a subset of skeleton genomes and subsequently places remaining genomes onto this backbone.

### Key results

Using PhyloPack, we computed evolutionary-based genome orderings for the ATB 2024-08 release [2] (approximately 2.4 million genomes).

When these orderings were used for genome compression, for major species (~93% of the collection):

- **< 15 GB** using MBGC [4].
- **~33 GB** using AGC [3] which supports random-access capabilities.

These results were obtained using the workflow described below.

## Contents
- Installation
- Quick Example
- Main commands
- Real application: Global ordering of ATB for compression using AGC

## Installation

```bash
git clone https://github.com/tam-km-truong/phylopack.git
cd phylopack
pip install .
```

### Requirements

* `Python >= 3.8`
* `attotree`

attotree [5] can be installed through Bioconda:

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

### High-level protocol

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
phylopack batch tests/data/genomes.txt -o debug/ -s 2 -v
```

Output: multiple batch files containing contiguous segments of the ordered genome list.

## Main commands

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
-k K for Kmers
-s-reference Sketch size for the skeleton genomes set
-s-placement Sketch size for the placement genomes set
--max-skeleton-genomes MAX_SKELETON_GENOMES
--min-skeleton-genomes MIN_SKELETON_GENOMES
[...]
```

For all options:

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

### Real-world application: Computing genome orderings and compression-ready batches for AllTheBacteria v2025-05

On the AllTheBacteria (ATB) collection [2], containing millions of genomes, PhyloPack was used to generate phylogeny-aware genome orderings and batches prior to compression.

Dataset scale - ATB version 2025-05:

- ~2.7 million genomes
- ~12k species 
- Largest cluster (*S. enterica*): ~700k genomes

### Species clustering

First, download the ATB metadata table [[OSF File-list](https://osf.io/zxfmy/files/7dpbh)].

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

The default value of `min_species_size_threshold` is 100. Species clusters containing fewer than 100 genomes are grouped into a pseudo-cluster named `remainder` for easier handling.

Example:

```bash
./suppl_scripts/species_clustering.sh \
    file_list_202505.tsv \
    sample \
    sylph_species \
    ~/ATB/clusters/ \
    100 \
    tar_xz
```

This step produces one genome list per species cluster.

### Compute a global ordering for each species

For each species cluster, compute a PhyloPack global ordering.

Example for *S. enterica*:

```bash
phylopack preorder \
    S_enterica_cluster.txt \
    --cut-point 0.05 \
    -t [Thread] \
    --max-skeleton-genomes 20000 \
    --min-skeleton-genomes 10000 \
    -o S_enterica_cluster_global_ordered.txt \
```

In this example, the skeleton set size is controlled by three parameters:

* `--cut-point 0.05` selects 5% of genomes as the skeleton set.
* `--min-skeleton-genomes 10000` ensures that clusters smaller than 10k genomes use all genomes as skeleton genomes, resulting in a full phylogenetic inference.
* `--max-skeleton-genomes 20000` caps the skeleton set size for very large clusters to reduce runtime and memory requirements.

For example, the *S. enterica* cluster contains approximately 700k genomes. With a cut point of 5%, the skeleton set would contain roughly 35k genomes. Because this exceeds the maximum threshold, the skeleton set is capped at 20k genomes.

### Refine the ordering within large clusters

For clusters requiring genome placement, an additional refinement step is performed after the initial global ordering. This improves local ordering accuracy within placement-heavy clusters without requiring phylogenetic inference on the entire dataset.

Refinement is relevant when the cluster is large enough that genome placement was used during the initial ordering step. For smaller clusters where all genomes were included in the skeleton phylogeny, this step can be skipped.

The ordered genome list is divided into non-overlapping windows of 5k genomes. 

This keeps each local phylogenetic inference computationally manageable while preserving the global-scale ordering produced by PhyloPack.

For each window, a local phylogenetic tree is inferred using Attotree [5]:

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

### Create balanced batches

The refined ordering is partitioned into balanced batches:

```bash
phylopack batch \
    -o ~/ATB/final_batches/ \
    --target-size 5000 \
    --batch-name cluster_ \
    cluster_refined_order.txt
```

The target size is 5000 for most clusters, except for clusters where higher diversity is expected (e.g. unknown-species clusters and the `remainder` pseudo-cluster), where the target size is reduced to 500.
This step aims to create contiguous batches while preserving the locality established by the ordering process.

The target size is application-specific. When random access is prioritized, smaller batches are preferred. When storage efficiency is prioritized, larger batches typically provide better compression ratios.

### Compress each batch with AGC

Each batch is compressed independently using AGC. The first genome in the batch is used as the reference genome.

```bash
for f in ~/ATB/final_batches/*; do
    ref=$(head -n 1 "$f")

    agc create \
        -a \
        -b 500 \
        -s 1500 \
        -t [Thread] \
        -o ~/ATB/results/agc/$(basename "$f" .txt).agc \
        -i "$f" \
        "$ref"
done
```

The resulting AGC archives can then be used for downstream storage and retrieval of bacterial genome collections.

## Bibliography:

[1] Břinda, Karel, Leandro Lima, Simone Pignotti, Natalia Quinones-Olvera, Kamil Salikhov, Rayan Chikhi, Gregory Kucherov, Zamin Iqbal, and Michael Baym. 2025. “Efficient and Robust Search of Microbial Genomes via Phylogenetic Compression.” Nature Methods 22 (4): 692–97.

[2] Hunt, Martin, Leandro Lima, Daniel Anderson, Jane Hawkey, Wei Shen, John Lees, and Zamin Iqbal. 2024. “AllTheBacteria - All Bacterial Genomes Assembled, Available and Searchable.” bioRxiv. https://doi.org/10.1101/2024.03.08.584059.

[3] Deorowicz, Sebastian, Agnieszka Danek, and Heng Li. 2023. “AGC: Compact Representation of Assembled Genomes with Fast Queries and Updates.” Bioinformatics (Oxford, England) 39 (3): btad097.

[4] Kowalski, Tomasz M. 2026. “MBGC2: Boosting Compression via Efficient Encoding of Approximate Matches in Genome Collections.” GigaScience 15 (giag008): giag008.

[5] Břinda, Karel. 2024. attotree 0.1.6. Zenodo.
https://doi.org/10.5281/zenodo.10950480.