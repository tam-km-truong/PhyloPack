import subprocess
import argparse

def parser():
    parser = argparse.ArgumentParser(description = 'First version for the script of phylopack, currently just do the preordering')
    parser.add_argument('input_genomes', help = 'the input list of genomes')
    parser.add_argument('-o','--output', help = 'Output file destination', default= '-') 

    return parser.parse_args()

def main():
    parser()


if __name__ == "__main__":
    main()