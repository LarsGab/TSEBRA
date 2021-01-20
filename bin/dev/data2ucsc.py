import argparse
import numpy as np
import subprocess as sp
import os
combiner_bin = "/home/lars/work/combiner/bin/"

def main():
    args = parseCmd()
    braker2_level = ['species_excluded', 'family_excluded', 'order_excluded']
    with open(args.data + '/species.tab', 'r') as file:
        species_list = file.read().split('\n')
    species_list = [s for s in species_list if s]
    for species in species_list:
        species_path = "{}/{}".format(args.data, species)
        anno =  '{}/anno/annot.gtf'.format(species_path)
        anno_out =  '{}/anno/annot_ucsc.gtf'.format(species_path)
        gtf2ucsc(anno, anno_out, 'anno')
        braker = "{}/braker1/braker_fixed.gtf".format(species_path)
        braker_out = "{}/braker1/braker_fixed_ucsc.gtf".format(species_path)
        gtf2ucsc(braker, braker_out, 'braker1')
        hints = "{}/braker1/hintsfile.gff".format(species_path)
        hints_out = "{}/braker1/hintsfile_ucsc.gff".format(species_path)
        gtf2ucsc(hints, hints_out, 'rnaseq', 'prothint')
        for level in braker2_level:
            braker = "{}/braker2/{}/braker_fixed.gtf".format(species_path, level)
            braker_out = "{}/braker2/{}/braker_fixed_ucsc.gtf".format(species_path, level)
            if os.path.exists(braker):
                gtf2ucsc(braker, braker_out, 'braker2_' + level)
                hints = "{}/braker2/{}/hintsfile.gff".format(species_path, level)
                hints_out = "{}/braker2/{}/hintsfile_ucsc.gff".format(species_path, level)
                gtf2ucsc(hints, hints_out, 'prothint', 'prothint')




def gtf2ucsc(gtf, out, name, mode='augustus'):
    color = list(map(str, list(np.random.choice(range(256), size=3))))

    cmd = "{}/gtf2ucsc.py --gtf {} --out {} --name {} --mode {} --color {}".\
        format(combiner_bin, gtf, out, name , mode, ','.join(color))
    print(cmd)
    sp.call(cmd, shell=True)
    if mode == 'augustus':
        cmd = "sed -i -e 's/chrChr/chr/g' " + out
        sp.call(cmd, shell=True)
    elif mode == 'prothint':
        out = '.'.join(out.split('.')[:-1]) + '_'
        for i in ['intron', 'start', 'stop', 'cds']:
            if os.path.exists("{}{}.gff".format(out, i)):
                cmd = "sed -i -e 's/chrChr/chr/g' {}{}.gff".format(out, i)
                sp.call(cmd, shell=True)

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
