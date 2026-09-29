import argparse

def add_split_args(parser):
    parser.add_argument('input', help='the input')

def run_macro_ordering(args):

    #task 1: run species clustering

    #task 2: if phylogenetically is on, override the size
    #get the middle genomes of a species
    #run py_attotree
    #rearrange the species order, and the genomes list to the new representative tree
    return

def main():
    parser = argparse.ArgumentParser(
        description='The description'
    )
    add_split_args(parser)
    args = parser.parse_args()

    run_macro_ordering(args)

if __name__ == "__main__":
    main()