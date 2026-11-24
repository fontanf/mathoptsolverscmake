"""Python bindings for MathOptSolversCMake's solver-independent part: the
MathOptModel structure and its MPS / LP / JSON / AMPL .nl readers and
writers.

    >>> import mathoptmodel as mom
    >>> model = mom.read("p0033.mps.gz")
    >>> model.number_of_variables, model.number_of_constraints
    (33, 16)
"""
import bz2
import gzip
import os

from mathoptmodel._version import __version__
from mathoptmodel._core import (
    ConstraintClass,
    ConstraintSense,
    Model,
    ObjectiveDirection,
    VariableType,
    dual,
    expand_defined_variables,
    read_json_bytes,
    read_lp_bytes,
    read_mps_bytes,
    read_nl_bytes,
    write_json_string,
    write_lp_string,
    write_mps_string,
    write_nl_string,
)

__all__ = [
    "__version__",
    "ConstraintClass", "ConstraintSense", "Model", "ObjectiveDirection", "VariableType",
    "FORMATS", "detect_format", "dual", "expand_defined_variables", "read", "write",
    "read_json_bytes", "read_lp_bytes", "read_mps_bytes", "read_nl_bytes",
    "write_json_string", "write_lp_string", "write_mps_string", "write_nl_string",
]

_READERS = {
    "mps": read_mps_bytes,
    "lp": read_lp_bytes,
    "json": read_json_bytes,
    "nl": read_nl_bytes,
}
_WRITERS = {
    "mps": write_mps_string,
    "lp": write_lp_string,
    "json": write_json_string,
    "nl": write_nl_string,
}
_OPENERS = {".gz": gzip.open, ".bz2": bz2.open}

FORMATS = tuple(_READERS)


def detect_format(path):
    """The format of `path` ("mps", "lp", "json" or "nl") from its extension,
    ignoring a trailing .gz/.bz2, or None if it isn't recognized."""
    base, ext = os.path.splitext(path.lower())
    if ext in _OPENERS:
        base, ext = os.path.splitext(base)
    ext = ext.lstrip(".")
    return ext if ext in _READERS else None


def read(path, format=None):
    """Read a model file, optionally .gz/.bz2 compressed. The format is
    guessed from the extension unless given ("mps", "lp", "json" or "nl")."""
    path = os.fspath(path)
    format = format or detect_format(path)
    if format not in _READERS:
        raise ValueError(f"Unrecognized model file format: {path}")
    opener = _OPENERS.get(os.path.splitext(path.lower())[1], open)
    with opener(path, "rb") as f:
        data = f.read()
    return _READERS[format](data)


def write(model, path, format=None):
    """Write a model file, .gz/.bz2 compressed if its name says so. The format
    is guessed from the extension unless given ("mps", "lp", "json" or "nl")."""
    path = os.fspath(path)
    format = format or detect_format(path)
    if format not in _WRITERS:
        raise ValueError(f"Unrecognized model file format: {path}")
    opener = _OPENERS.get(os.path.splitext(path.lower())[1], open)
    with opener(path, "wt") as f:
        f.write(_WRITERS[format](model))
