import argparse

def add_split_args(parser):
    parser.add_argument('input', help='the input')

def run_species_clustering(args):
    return

def main():
    parser = argparse.ArgumentParser(
        description='The description'
    )
    add_split_args(parser)
    args = parser.parse_args()

    run_species_clustering(args)

if __name__ == "__main__":
    main()