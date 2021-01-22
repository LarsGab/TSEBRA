#!/usr/bin/env python3
import os
import sys
import pytest

sys.path.append('/home/lars/work/prevco/bin')

from genome_anno import Transcript, Anno, NotGtfFormat


testDir = os.path.abspath(os.path.dirname(__file__))
anno1 = testDir + '/genome_anno/anno1.gtf'
anno_format_error = testDir + '/genome_anno/format_error.gtf'
anno_missing_gid = testDir + '/genome_anno/missing_gid.gtf'
tx1 = testDir + '/genome_anno/tx1.gtf'
tx1_args = ('tx1', 'g.tx1', '3L', 'GeneMark.hmm')

@pytest.fixture
def transcript():
    return Transcript(*tx1_args)

@pytest.fixture
def file_tx1():
    result = []
    with open(tx1, 'r') as file:
        for line in file.readlines():
            line = line.strip('\n').split('\t')
            result.append(line)
    return result

@pytest.fixture
def file_anno1():
    result = []
    with open(anno1, 'r') as file:
        for line in file.readlines():
            line = line.strip('\n')
            result.append(line)
    return result

@pytest.fixture
def transcript_tx1(file_tx1):
    t = Transcript(*tx1_args)
    for line in file_tx1:
        t.add_line(line)
    return t

@pytest.fixture
def anno_anno1():
    anno = Anno(anno1, 'anno1')
    anno.addGtf()
    return anno

def test_transcript_defaults(transcript):
    assert transcript.id == tx1_args[0]
    assert transcript.gene_id == tx1_args[1]
    assert transcript.chr == tx1_args[2]
    assert transcript.source_anno == tx1_args[3]

def test_transcript_add_lines(transcript_tx1, file_tx1):
    list = []
    for key in transcript_tx1.transcript_lines.keys():
        list += transcript_tx1.transcript_lines[key]
    assert len(list) == len(file_tx1)
    for line in list:
        assert line in file_tx1

def test_transcript_find_lines(transcript_tx1):
    missing = {"intron" : [['3L', 'GeneMark.hmm', 'intron', 18462541, 18462718, \
            '0', '-', '.', \
            'gene_id "g.tx1"; transcript_id "tx1";']], \
        "start_codon" : [['3L', 'GeneMark.hmm', 'start_codon', 18463066, 18463068, \
            '.', '-', '.', 'gene_id "g.tx1"; transcript_id "tx1";']], \
        "transcript" : [['3L', 'GeneMark.hmm', 'transcript', 18462228, 18463068, \
            '0', '-', '.', 'tx1']]}
    transcript_tx1.add_missing_lines()
    for key in missing.keys():
        for line in missing[key]:
            assert line in transcript_tx1.transcript_lines[key]

def test_anno_read_file(anno_anno1, file_anno1):
    gtf_anno = anno_anno1.get_gtf().split('\n')
    gtf_anno = [g.split('\t')[:8] for g in gtf_anno]
    file_anno1 = [f.split('\t')[:8] for f in file_anno1]
    assert len(gtf_anno) == len(file_anno1)
    for line in file_anno1:
        assert line in gtf_anno

def test_format_error():
    anno = Anno(anno_format_error, 'error_anno')
    with pytest.raises(NotGtfFormat):
        anno.addGtf()

def test_missing_gid(file_anno1):
    anno = Anno(anno_missing_gid, 'anno1')
    anno.addGtf()
    gtf_anno = anno.get_gtf().split('\n')
    gtf_anno = [g.split('\t')[:8] for g in gtf_anno]
    file_anno1 = [f.split('\t')[:8] for f in file_anno1]
    assert len(gtf_anno) == len(file_anno1)
    for line in gtf_anno:
        assert line in file_anno1




if __name__ == '__main__':
    os.mkdir(tempDir)
    #sys.path.append(testDir + "/../bin")
