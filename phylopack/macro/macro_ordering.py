import argparse

def add_split_args(parser):
    parser.add_argument('input', help='the input')

def run_macro_ordering(args):
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