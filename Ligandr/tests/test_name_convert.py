import pytest
import sys
from Ligandr import NameToSMILES

def test_xml_grab():
    assert NameToSMILES.getSMILE("SAM") == "C[S@@+](CC[C@@H](C(=O)[O-])N)C[C@@H]1[C@H]([C@H]([C@@H](O1)n2cnc3c2ncnc3N)O)O"


def test_pH():
    assert ("O-" in NameToSMILES.addH("CC(=O)O")) == True
