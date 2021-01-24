#!/usr/bin/env python3
import argparse
import subprocess as sp
import os


augustus_script = os.path.abspath(os.path.dirname(__file__)) + "/augustus2browser.pl"
def main():
    args = parseCmd()
    mode = args.mode
    color = args.color
    if mode == 'augustus':
        gtf = ""
        temp = args.gtf + '_temp.gtf'
        tool = "augustus2browser.pl"
        cmd = "{} < {} > {}".format(augustus_script, args.gtf, temp)
        sp.call(cmd, shell=True)
        with open(temp, 'r') as file:
            temp_file = file.readlines()

        if args.name:
            gtf += "{}name={} color={} decription=\"{}\"{}".format(temp_file[0].split('name=')[0], \
                args.name, color, args.name, temp_file[0].split('\"')[-1])
        for line in temp_file:
            line = line.split('\t')
            if len(line) == 9:
                '''
                line[0] = line[0].split('_')
                if len(line[0]) == 3:
                    line[0] = line[0][1]
                else:
                    line[0] = line[0][0]
                line[0] = 'c' + line[0][1:]
                '''
                line[0] = 'chr' + line[0]
                gtf += '\t'.join(line)

        with open(args.out, 'w+') as file:
            file.write(gtf)
        os.remove(temp)
    elif mode == 'prothint':
        name = args.name + "_intron"
        out_intron = ["track name={} decription=\"{}\" visibility=3 color={}\n".format(name, name, color)]
        name = args.name + "_start"
        out_start  = ["track name={} decription=\"{}\" visibility=3 color={}\n".format(name, name, color)]
        name = args.name + "_stop"
        out_stop  = ["track name={} decription=\"{}\" visibility=3 color={}\n".format(name, name, color)]
        name = args.name + "_cds"
        out_cds  = ["track name={} decription=\"{}\" visibility=3 color={}\n".format(name, name, color)]
        with open(args.gtf, 'r') as file:
            gff = file.readlines()
        trans_numb = 1
        for line in gff:
            line = line.split('\t')
            line[0] = 'chr' + line[0]
            '''
            line[0] = line[0].split('_')
            if len(line[0]) == 3:
                line[0] = line[0][1]
            else:
                line[0] = line[0][0]
            line[0] = 'c' + line[0][1:]
            '''
            type = line[2]
            line[2] = 'CDS'
            line[8] = '{}; transcript_id=\"t.{}\"\n'.format(line[8].strip('\n').replace(';', '_'), trans_numb)
            trans_numb += 1
            if type == 'intron':
                out_intron.append('\t'.join(line))
            elif type == 'start':
                out_start.append('\t'.join(line))
            elif type == 'stop':
                out_stop.append('\t'.join(line))
            elif type == 'CDSpart':
                out_cds.append('\t'.join(line))
            else:
                print(line)
        out = '.'.join(args.out.split('.')[:-1])
        if len(out_intron) > 1:
            with open (out + "_intron.gff", 'w+') as file:
                file.write(''.join(out_intron))
        if len(out_start) > 1:
            with open (out + "_start.gff", 'w+') as file:
                file.write(''.join(out_start))
        if len(out_stop) > 1:
            with open (out + "_stop.gff", 'w+') as file:
                file.write(''.join(out_stop))
        if len(out_cds) > 1:
            with open (out + "_cds.gff", 'w+') as file:
                file.write(''.join(out_cds))
    elif mode == 'seed_region':
        out = ["track name={} decription=\"{}\" visibility=3 color={}\n".format(args.name, args.name, color)]
        with open(args.gtf, 'r') as file:
            gff = file.readlines()
        trans_numb = 1
        for line in gff:
            line = line.split('\t')
            line[0] = 'chr' + line[0]
            line[2] = 'CDS'
            line[8] = 'transcript_id=\"{}\";\n'.format(line[8].replace('\n', ''))
            trans_numb += 1
            out.append('\t'.join(line))

        with open (args.out, 'w+') as file:
            file.write(''.join(out))
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
    parser.add_argument('--name', type=str,
        help='')
    parser.add_argument('--color', type=str,
        help='')
    parser.add_argument('--mode', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
