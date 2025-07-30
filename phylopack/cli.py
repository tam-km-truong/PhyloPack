# cli.py
import shutil
import sys
import argparse
from phylopack.preorder.preorder import add_preorder_parser


def check_dependencies(tools=["mash", "quicktree", "attotree"]):
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        print(f"Required external tools: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(prog="phylopack", description="Phylopack CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Add the 'preorder' command from the preorder module
    add_preorder_parser(subparsers)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    check_dependencies()
    main()