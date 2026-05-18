"""Make the Sentinel package importable when running pytest from anywhere."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
if SRC_ROOT not in sys.path:
    sys.path.insert(0, SRC_ROOT)
