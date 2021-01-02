import argparse
import os
import subprocess as sp
combiner_bin = "/home/lars/work/combiner/bin/"
def main():
    args = parseCmd()
    braker2_level = ['species_excluded', 'family_excluded', 'order_excluded']
    with open(args.data + '/species.tab', 'r') as file:
        species_list = file.read().split('\n')
    species_list = [s for s in species_list if s]

    full_eval = []
    summary_eval = []
    header = ''
    for species in species_list:
        test_id = "{}_braker1".format(species)
        species_path = "{}/{}".format(args.data, species)
        braker = "{}/braker1/braker_fixed.gtf".format(species_path)
        accuracies = eval("{}/anno/".format(species_path), braker)
        tx_gene = tx_per_gene(braker)
        txt = [a[0] for a in accuracies]
        txt.append(tx_gene[2][0])
        if not header:
            header = '\t'.join(txt)
        elif not header == '\t'.join(txt):
            print('Accuracy assessment output for {}'.format(test_id))
        txt = [a[1] for a in accuracies]
        txt.append(tx_gene[2][1])
        full_eval.append([test_id] + txt)
        txt = list(map(float, txt))
        summary_eval.append([test_id, sum(txt[2:-1]), txt[-1]])
        for level in braker2_level:
            test_id = "{}_braker2_{}".format(species, level)
            braker2 = "{}/braker2/{}/braker_fixed.gtf".format(species_path, level)
            if os.path.exists(braker2):
                accuracies = eval("{}/anno/".format(species_path), braker2)
                tx_gene = tx_per_gene(braker2)
                txt = [a[0] for a in accuracies]
                txt.append(tx_gene[2][0])
                if not header:
                    header = '\t'.join(txt)
                elif not header == '\t'.join(txt):
                    print('Accuracy assessment output for {}'.format(test_id))
                txt = [a[1] for a in accuracies]
                txt.append(tx_gene[2][1])
                full_eval.append([test_id] + txt)
                txt = list(map(float, txt))
                summary_eval.append([test_id, sum(txt[2:-1]), txt[-1]])

    full_eval_out = '# Mode\t{}\n'.format(header)
    for line in full_eval:
        full_eval_out += '# {}\n'.format('\t'.join(line))
    full_eval_out += '\\begin{table}[h]\n\\centering\n\\begin{tabular}{c||c|c||c|c||c|c||c}\n'
    full_eval_out += 'Mode&{}\\\\\n\\hline\\hline\n'.format('&'.join(header.split('\t')))
    for line in full_eval:
        full_eval_out += '{}\\\\\n\\hline\n'.format('&'.join(line))
    full_eval_out += '\\end{tabular}\n\\end{table}'
    full_eval_out = full_eval_out.replace('_', ' ')
    with open(args.out + 'full_evaluation.txt', 'w+') as file:
        file.write(full_eval_out)

    summary_eval = []
    for line in full_eval:
        summary_eval.append([line[0], sum(list(map(float, line[2:-1])))/4, line[-1]])
    summary_eval = sorted(summary_eval, key= lambda s:s[0])
    summary_out = '# {}\n'.format('\t'.join(map(str,[s[0] for s in summary_eval])))
    summary_out += '# {}\n'.format('\t'.join(map(str,[s[1] for s in summary_eval])))
    summary_out += '&'.join(map(str,[s[1] for s in summary_eval]))
    summary_out = summary_out.replace('_', ' ')
    with open(args.out + 'summary_eval.txt', 'w+') as file:
        file.write(summary_out)

def eval(anno, prediction):
    cmd = "{}/compute_accuracies.sh {}/annot.gtf {}/pseudo.gff3 {} gene trans cds".\
        format(combiner_bin, anno, anno, prediction)
    print(cmd)
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if stderr.decode():
        print(stderr.decode())
    stdout = stdout.decode()
    return [s.split('\t') for s in stdout.split('\n') if s]
def tx_per_gene(gtf):
    cmd = "{}/tx_per_gene.py --gtf {}".format(combiner_bin, gtf)
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if stderr.decode():
        print(stderr.decode())
    stdout = stdout.decode()
    stdout = [s for s in stdout.split('\n') if s]
    return [s.split('\t') for s in stdout[-3:]]
def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data', type=str,
        help='')
    parser.add_argument('--out', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
