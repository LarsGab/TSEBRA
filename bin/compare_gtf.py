import argparse

def main():
    args = parseCmd()
    gtf1_id = []
    gtf2_id = []

    with open(args.gtf1, 'r') as file:
        for line in file.readlines():
            line=line.split('\t')
            if line[2] == 'transcript':
                id = line[8]
            #else:
                #id = line[8].split('transcript_id "')[1].split('";')[0]
                gtf1_id.append(id)
    with open(args.gtf2, 'r') as file:
        for line in file.readlines():
            line=line.split('\t')
            if line[2] == 'transcript':
                id = line[8]
            #else:
                #id = line[8].split('transcript_id "')[1].split('";"')[0]
                gtf2_id.append(id)

    unique1 = [u for u in gtf1_id if u not in gtf2_id]
    unique2 = [u for u in gtf2_id if u not in gtf1_id]
    with open(args.out + '_1', 'w+') as file:
        file.write('\n'.join(unique1))
    with open(args.out + '_2', 'w+') as file:
        file.write('\n'.join(unique2))

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--gtf1', type=str,
        help='')
    parser.add_argument('--gtf2', type=str,
        help='')
    parser.add_argument('--out', type=str,
        help='')

    return parser.parse_args()

if __name__ == '__main__':
    main()
