# cli.py

import argparse
from phylopack.preorder.preorder import add_preorder_parser

def main():
    parser = argparse.ArgumentParser(prog="phylopack", description="Phylopack CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Add the 'preorder' command from the preorder module
    add_preorder_parser(subparsers)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()