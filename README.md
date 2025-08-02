# PhyloPack

**PhyloPack** is a tool for computing phylogenetic-based genome orderings and batching.

## Installation

```bash
git clone https://github.com/tam-km-truong/phylopack.git
cd phylopack
pip install .
````

### Requirements

* `Python >= 3.8`
* `mash`
* `quicktree`
* `attotree`

Make sure the above tools are installed and available in your `PATH`.

## Usage

```bash
phylopack preorder input-genomes.txt -o output.txt
```

For available options:

```bash
phylopack preorder -h
```
