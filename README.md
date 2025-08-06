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
* `attotree`

Make sure the above tools are installed and available in your `PATH`.

## Usage

```bash
phylopack preorder tests/data/genomes.txt --cut-point 0.2 -o ./debug/out.txt -v
```

For available options:

```bash
phylopack preorder -h
```
