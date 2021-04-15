import sys
import os
import pytest
import csv
sys.path.append('/home/lars/work/prevco/bin')

from evidence import NotGtfFormat, AttributeMissing, Hint, Evidence

testDir = os.path.abspath(os.path.dirname(__file__))

@pytest.fixture
def hints1():
    hints = []
    with open(testDir + '/evidence/hint1.gff') as file:
        hints_tab = csv.reader(file, delimiter='\t')
        for line in hints_tab:
            hints.append(line)
    return hints

def test_hint(hints1):
    hint = Hint(hints1[0])
    assert list(map(str, hint.hint2list())) == hints1[0]
    hint = Hint(hints1[2])
    assert  list(map(str, hint.hint2list())) == hints1[2]

def test_hint_error(hints1):
    with pytest.raises(AttributeMissing):
        Hint(hints1[1])
    with pytest.raises(NotGtfFormat):
        Hint(hints1[3])

def test_get_hint():
    evi = Evidence()
    evi.add_hintfile(testDir + '/evidence/hint3.gff')
    mult = evi.get_hint('3R','801','899','intron','+')
    assert sum(mult.values()) == 28
