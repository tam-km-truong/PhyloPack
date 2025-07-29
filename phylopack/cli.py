import tempfile
import os
import shutil
import argparse

def add_preorder_args(parser):
    parser.add_argument('--input-genomes', required=True, help='Path to the input list of genomes')
    parser.add_argument('--output', default='.', help='Output path (default: current folder)')
    parser.add_argument('-k', type=int, default=21, help='K-mer size (default: 21)')
    parser.add_argument('-s', type=int, default=10000, help='Sketch size (default: 10000)')
    parser.add_argument('-t', type=int, default=10, help='Number of threads (default: 10)')
    parser.add_argument('-m', choices=['nj', 'upgma'], default='nj', help='Tree method: nj or upgma (default: nj)')

    parser.add_argument('--output-tree', help='Custom path to output tree .nw file')
    parser.add_argument('--output-std-tree', help='Custom path to output standardized tree')
    parser.add_argument('--leaf-order', help='Custom path to leaf order output file')
    parser.add_argument('--node-info', help='Custom path to node info file')

    parser.add_argument('--statistic', action='store_true', help='Output statistics file')
    parser.add_argument('--statistic-file-type', choices=['json', 'csv'], default='json', help='Statistics file type (default: json)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print logs')

def run_preorder(args):
    return

def main():
    return

if __name__ == "__main__":
    main()