import sys
import os
import pytest

sys.path.append('/home/lars/work/prevco/bin')

from evidence import NotGtfFormat, AttributeMissing, Hint, Evidence

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

def test_get_hint():
    evi = Evidence()
    evi.add_hintfile(testDir + '/evidence/hint3.gff')
    mult = evi.get_hint('3R','801','899','intron','+')
    assert mult == 28
