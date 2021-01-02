import sys
import os
import pytest

sys.path.append('/home/lars/work/')

from combiner.bin.evidence import NotGtfFormat, AttributeMissing, Hint, Hintfile

testDir = os.path.abspath(os.path.dirname(__file__))

@pytest.fixture
def hints1():
    with open(testDir + '/evidence/hint1.gff') as file:
        hints = file.readlines()
    return hints


def test_hint(hints1):
    hint = Hint(hints1[0])
    assert list(map(str, hint.hint2list())) == hints1[0].split('\t')
    hint = Hint(hints1[2])
    assert  list(map(str, hint.hint2list())) == hints1[2].split('\t')

def test_hint_error(hints1):
    with pytest.raises(AttributeMissing):
        Hint(hints1[1])
    with pytest.raises(NotGtfFormat):
        Hint(hints1[3])

def test_hint_range():
    hint3 = Hintfile(testDir + '/evidence/hint3.gff')
    h_list = hint3.hints_in_range(550, 899, '3R')
    with open(testDir + '/evidence/hint3.gff') as file:
        hints = file.read().split('\n')
    correct = hints[3:5]
    correct.append(hints[6])
    correct.append(hints[8])
    correct = [c.split('\t') for c in correct]
    h_list = [list(map(str, l)) for l in h_list]
    print(correct)
    print(h_list)
    assert h_list == correct
