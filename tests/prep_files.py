import os
testDir = os.path.abspath(os.path.dirname(__file__))

def genome_anno():
    anno1 = testDir + '/genome_anno/anno1.gtf'
    orig = []
    with open(anno1, 'r') as file:
        for line in file.readlines():
            line = line.strip('\n')
            orig.append(line)
    orig = [f.split('\t') for f in orig]

    anno = orig
    anno[1][8] = 'gene_id "7789_g";'
    anno = ['\t'.join(map(str, line)) for line in anno]
    with open(testDir + '/genome_anno/format_error.gtf', 'w+') as file:
        file.write('\n'.join(anno))

    anno = orig
    anno[1][8] = 'transcript_id "7789_t";'
    anno[6][8] = 'transcript_id "g5980.t1";'
    for line in anno:
        if 'transcript_id "g11080.t1";' in line[8]:
            line[8] = 'transcript_id "g11080.t1";'
    anno = ['\t'.join(map(str, line)) for line in anno]
    with open(testDir + '/genome_anno/missing_gid.gtf', 'w+') as file:
        file.write('\n'.join(anno))




def get_anno(tx_dict, phase):
        template = ['3R', 'AUGUSTUS', '', '', '', phase, '+', '0', '']
        anno = []
        for key in tx_dict:
            coord = tx_dict[key]
            template[8] = 'transcript_id "{}"; gene_id "{}";'.format(key, key + '_g')
            type = 'exon'
            pos = coord[0]
            for c in coord[1:]:
                line = template.copy()
                line[2] = type
                line[3] = pos
                pos += c
                line[4] = pos
                if type == 'intron':
                    line[3] += 1
                    line[4] -= 1
                anno.append(line)
                if type == 'exon':
                    line = line.copy()
                    line[2] = 'CDS'
                    anno.append(line)
                    type = 'intron'
                else:
                    type = 'exon'
            line = template.copy()
            line[2] = 'transcript'
            line[3] = str(coord[0])
            line[4] = str(pos)
            line[8] = key
            anno.append(line)
        return anno

def list2string(gtf):
    gtf = ['\t'.join(map(str, g)) for g in gtf]
    return '\n'.join(gtf)

def graph():
    dir = testDir + '/graph/'
    #example 1
    anno1_txs = { 't1' : [100, 100, 100, 100], \
            't2' : [700, 100, 100, 100, 100, 100], \
            't3' : [1500, 100]}
    anno1 = get_anno(anno1_txs, '0')
    with open(dir + 'ex1_anno1.gtf', 'w+') as file:
        file.write(list2string(anno1))

    anno2_txs = {   't1' : [250, 250, 100, 150],
                    't2' : [1050, 200],
                    't3' : [1700, 100]}
    anno2 = get_anno(anno2_txs, '0')
    with open(dir + 'ex1_anno2.gtf', 'w+') as file:
        file.write(list2string(anno2))

    #example 2
    anno1_txs = { 't1' : [200, 100]}
    anno1 = get_anno(anno1_txs, '0')
    with open(dir + 'ex2_anno1.gtf', 'w+') as file:
        file.write(list2string(anno1))

    anno2_txs = {   't1' : [100, 100], \
                    't2' : [301, 99]}
    anno2 = get_anno(anno2_txs, '0')
    with open(dir + 'ex2_anno2.gtf', 'w+') as file:
        file.write(list2string(anno2))

    #example 3
    anno1_txs = { 't1' : [100, 200, 200, 200, 200, 200]}
    anno1 = get_anno(anno1_txs, '0')
    with open(dir + 'ex3_anno1.gtf', 'w+') as file:
        file.write(list2string(anno1))

    anno2_txs = {   't1' : [110, 90, 600, 200], \
                    't2' : [350, 100]}
    anno2 = get_anno(anno2_txs, '0')
    with open(dir + 'ex3_anno2.gtf', 'w+') as file:
        file.write(list2string(anno2))

    #example 4
    anno1_txs = { 't1' : [100, 100, 100, 100]}
    anno1 = get_anno(anno1_txs, '0')
    with open(dir + 'ex4_anno1.gtf', 'w+') as file:
        file.write(list2string(anno1))

    anno2_txs = { 't1' : [101, 100, 100, 100]}
    anno2 = get_anno(anno2_txs, '1')
    with open(dir + 'ex4_anno2.gtf', 'w+') as file:
        file.write(list2string(anno2))

def evidence():
    dir = testDir + '/evidence/'
    hint_test_file1 = ['3L\tProtHint\tintron\t5812862\t5812941\t24\t-\t.\tsrc=M;mult=24;pri=4\n', \
        '3L\tProtHint\tintron\t12291242\t12291299\t8\t-\t.\ttranscript_id="t1"\n', \
        '3L\tProtHint\tintron\t12291242\t12291299\t8\t-\t.\tsrc=M;pri=4\n',
        '3L\tProtHint\tintron\t12291242\t']
    with open(dir + 'hint1.gff', 'w+') as file:
        file.write(''.join(hint_test_file1))

    hint_test_file2 = ['3L\tProtHint\tintron\t5812862\t5812941\t24\t-\t.\tsrc=M;mult=24;pri=4\n', \
        '3L\tProtHint\tintron\t12291242\t12291299\t8\t-\t.\tsrc=M;mult=8;pri=4\n', \
        '3R\tProtHint\tintron\t17440148\t17440207\t25\t-\t.\tsrc=M;mult=25;pri=4\n', \
        '2R\tProtHint\tintron\t5760114\t5760177\t23\t-\t.\tsrc=M;mult=23;pri=4\n', \
        '2R\tProtHint\tintron\t6210484\t6210546\t21\t-\t.\tsrc=M;mult=21;pri=4\n', \
        '3L\tProtHint\tintron\t20527281\t20527592\t25\t+\t.\tsrc=M;mult=25;pri=4\n', \
        '2L\tProtHint\tintron\t12400752\t12400814\t24\t+\t.\tsrc=M;mult=24;pri=4\n', \
        '2R\tProtHint\tintron\t14988084\t14988142\t25\t-\t.\tsrc=M;mult=25;pri=4\n', \
        '2L\tProtHint\tintron\t6667531\t6667670\t5\t-\t.\tsrc=M;mult=5;pri=4\n', \
        '3R\tProtHint\tintron\t5537551\t5537605\t22\t+\t.\tsrc=M;mult=22;pri=4\n', \
        '3R\tProtHint\tintron\t20813612\t20813665\t12\t-\t.\tsrc=M;mult=12;pri=4\n', \
        'X\tProtHint\tintron\t2145714\t2147174\t25\t+\t.\tsrc=M;mult=25;pri=4\n', \
        '3L\tProtHint\tintron\t8114197\t8114256\t25\t-\t.\tsrc=M;mult=25;pri=4\n', \
        'X\tProtHint\tintron\t11048602\t11048941\t25\t+\t.\tsrc=M;mult=25;pri=4\n', \
        '2L\tProtHint\tintron\t3807462\t3807524\t18\t+\t.\tsrc=M;mult=18;pri=4\n', \
        '3R\tProtHint\tintron\t27059120\t27059364\t19\t-\t.\tsrc=M;mult=19;pri=4\n', \
        '2R\tProtHint\tintron\t13821370\t13821432\t24\t-\t.\tsrc=M;mult=24;pri=4\n', \
        'X\tProtHint\tintron\t8173462\t8173860\t6\t-\t.\tsrc=M;mult=6;pri=4\n', \
        'X\tProtHint\tintron\t13270643\t13271481\t16\t-\t.\tsrc=M;mult=16;pri=4\n', \
        'X\tProtHint\tintron\t2079645\t2079714\t25\t-\t.\tsrc=M;mult=25;pri=4\n']
    with open(dir + 'hint2.gff', 'w+') as file:
        file.write(''.join(hint_test_file2))

    hint_test_file3 = []
    hint_test_file3.append(get_hint(100, 102, 'start_codon'))
    hint_test_file3.append(get_hint(501, 599, 'intron'))
    hint_test_file3.append(get_hint(501, 599, 'intron', src='P', mult=14))
    hint_test_file3.append(get_hint(698, 700, 'stop_codon'))
    hint_test_file3.append(get_hint(801, 899, 'intron'))
    hint_test_file3.append(get_hint(801, 899, 'intron', chr='2L'))
    hint_test_file3.append(get_hint(801, 899, 'intron', src='P', mult=24))
    hint_test_file3.append(get_hint(801, 949, 'intron'))
    hint_test_file3.append(get_hint(801, 899, 'intron', strand='-'))
    hint_test_file3.append(get_hint(1001, 1099, 'intron'))
    hint_test_file3.append(get_hint(1198, 1200, 'stop_codon'))
    hint_test_file3.append(get_hint(1601, 1699, 'intron'))
    with open(dir + 'hint3.gff', 'w+') as file:
        file.write('\n'.join(hint_test_file3))


def get_hint(start, end, type, strand='+', chr='3R', score=10, mult=2, pri=4, src='E'):
    att = 'src={};mult={};pri={}'.format(src,mult,pri)
    template = [chr, 'AUGUSTUS', type, start, end, score, '+', '.', att]
    return '\t'.join(map(str, template))

def get_feature():
    dir = testDir + '/graph/'
    result = []
    with open('/home/lars/work/combiner/example/braker1/braker_fixed.gtf', 'r') as file:
        for line in file.readlines():
            if 'g7604.t1' in line or 'g7603.t1' in line or 'g7605.t1' in line:
                result.append(line)
    with open(dir + 'ex_feature_anno1.gtf', 'w+') as file:
        file.write(''.join(result))

    result = []
    with open('/home/lars/work/combiner/example/braker2/braker.gtf', 'r') as file:
        for line in file.readlines():
            if 'g7700.t1' in line or 'g7701.t1' in line:
                result.append(line)
    with open(dir + 'ex_feature_anno2.gtf', 'w+') as file:
        file.write(''.join(result))

    result = []
    with open('/home/lars/work/combiner/example/braker1/hintsfile.gff', 'r') as file:
        for line in file.readlines():
            line = line.split('\t')
            if len(line) > 8:
                if int(line[3]) >= 21737000 and int(line[4]) <= 21750000 \
                    and line[0] == '3R' and not line[2] == 'CDSpart':
                    result.append(line)
    result = ['\t'.join(r) for r in result]
    with open(dir + 'ex_feature_hint1.gff', 'w+') as file:
        file.write(''.join(result))

    result = []
    with open('/home/lars/work/combiner/example/braker2/hintsfile.gff', 'r') as file:
        for line in file.readlines():
            line = line.split('\t')
            if len(line) > 8:
                if int(line[3]) >= 21737000 and int(line[4]) <= 21750000 \
                    and line[0] == '3R' and not line[2] == 'CDSpart':
                    result.append(line)
    result = ['\t'.join(r) for r in result]
    with open(dir + 'ex_feature_hint2.gff', 'w+') as file:
        file.write(''.join(result))

if __name__ == '__main__':
    #genome_anno()
    #graph()
    #evidence()
    get_feature()
