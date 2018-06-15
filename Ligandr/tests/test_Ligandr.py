"""
Unit and regression test for the Ligandr package.
"""

# Import package, test suite, and other packages as needed
import Ligandr
import pytest
import sys
from Ligandr import NameToSMILES


def test_Ligandr_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "Ligandr" in sys.modules
