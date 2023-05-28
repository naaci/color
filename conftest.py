import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent))

from numpy import asarray, isclose
from numpy.random import random

# sys.path.insert(0, "./colorspaces")


def _is_inverse(g, f):
    c = random(3)
    assert isclose(c, g(f(c))).all()
    assert isclose(c, f(g(c))).all()


def _test_delta(f, value):
    c = asarray([0.1, 0.2, 0.3])
    assert f(c, c * 2) == value
