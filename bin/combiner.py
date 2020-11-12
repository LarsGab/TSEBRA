#!/usr/bin/env python3
import argparse

def main():
    args = parseCmd()

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--arg', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
