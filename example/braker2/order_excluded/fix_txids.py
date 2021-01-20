import argparse


def main():
    args = parseCmd()
    with open(args.gtf, 'r') as file:
        gtf = file.readlines()

    result = ''
    for line in gtf:
        line = line.split('\t')
        if line[2] == 'gene'or line[2] == 'transcript':
                result += '\t'.join(line)
        else:
            temp = line[8].split('transcript_id "')
            transcript_id = temp[1].split('";')[0]
            transcript_id = transcript_id.split('_')
            if len(transcript_id[-1]) == 1:
                transcript_id = '_'.join(transcript_id[-2:])
            else:
                transcript_id = transcript_id[-1]
            line[8] = '{}transcript_id "{}";{}'.format(temp[0], transcript_id, '";'.join(temp[1].split('";')[1:]))
            temp = line[8].split('gene_id "')
            gene_id = temp[1].split('";')[0]
            gene_id = gene_id.split('_')
            if len(gene_id[-1]) == 1:
                gene_id = '_'.join(gene_id[-2:])
            else:
                gene_id = gene_id[-1]
            line[8] = '{}gene_id "{}";{}'.format(temp[0], gene_id, '";'.join(temp[1].split('";')[1:]))
            result += '\t'.join(line)

    with open(args.out, 'w+') as file:
        file.write(result)
def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--gtf', type=str,
        help='')
    parser.add_argument('--out', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
