"""Python bindings for MathOptSolversCMake: the MathOptModel structure, its
MPS / LP / JSON / AMPL .nl readers and writers, and solving with the solvers
available in the build.

    >>> import mathoptmodel as mom
    >>> model = mom.read("p0033.mps.gz")
    >>> model.number_of_variables, model.number_of_constraints
    (33, 16)
    >>> mom.solve(model, "highs", verbose=False).objective_value
    3089.0

Models can also be built in batches of NumPy arrays: see Model.add_variables,
Model.add_constraints, Model.add_defined_variables and Model.set_objective.
"""
import bz2
import functools
import gzip
import os

import numpy as np

from mathoptmodel._version import __version__
from mathoptmodel._core import (
    ConstraintClass,
    ConstraintSense,
    Model,
    ObjectiveDirection,
    SolveResult,
    Solver,
    VariableType,
    available_solvers,
    dual,
    expand_defined_variables,
    read_json_bytes,
    read_lp_bytes,
    read_mps_bytes,
    read_nl_bytes,
    write_json_file,
    write_json_string,
    write_lp_file,
    write_lp_string,
    write_mps_file,
    write_mps_string,
    write_nl_file,
    write_nl_string,
    solve_with as _solve_with,
)

__all__ = [
    "__version__",
    "ConstraintClass", "ConstraintSense", "Model", "ObjectiveDirection", "VariableType",
    "SolveResult", "Solver", "available_solvers", "solve",
    "FORMATS", "detect_format", "dual", "expand_defined_variables", "read", "write",
    "read_json_bytes", "read_lp_bytes", "read_mps_bytes", "read_nl_bytes",
    "write_json_string", "write_lp_string", "write_mps_string", "write_nl_string",
    "write_json_file", "write_lp_file", "write_mps_file", "write_nl_file",
]


def _accept_sequences(method):
    """Let a Model building method take lists/tuples where it takes arrays
    (nanobind's array arguments only accept array objects). Lists of str
    (names) are left as they are."""

    def convert(value):
        if isinstance(value, (list, tuple)) and not (value and isinstance(value[0], str)):
            return np.asarray(value)
        return value

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        args = [convert(value) for value in args]
        kwargs = {key: convert(value) for key, value in kwargs.items()}
        return method(self, *args, **kwargs)

    return wrapper


for _name in ("add_variables", "add_constraints", "add_defined_variables", "set_objective"):
    setattr(Model, _name, _accept_sequences(getattr(Model, _name)))
del _name

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
_FILE_WRITERS = {
    "mps": write_mps_file,
    "lp": write_lp_file,
    "json": write_json_file,
    "nl": write_nl_file,
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
    opener = _OPENERS.get(os.path.splitext(path.lower())[1])
    # An uncompressed file is written by the C++ writer directly, rather than
    # through a str (but on Windows, only for ASCII paths, which the C++
    # writers open with the narrow-string API).
    if opener is None and (os.name != "nt" or path.isascii()):
        _FILE_WRITERS[format](model, path)
        return
    with (opener or open)(path, "wt") as f:
        f.write(_WRITERS[format](model))


_SOLVERS = {solver.name.lower(): solver for solver in Solver}


def solve(model, solver="highs", time_limit=None, verbose=True):
    """Solve a model and return a SolveResult.

    solver is a Solver or its name ("cbc", "highs", "xpress", "cplex",
    "gurobi", "knitro", "ipopt" or "scip", in any case); it must be available
    in this build (see available_solvers()). time_limit is in seconds (None
    for no limit). With verbose=False, the solver's log is not written to the
    standard output.

    The GIL is released while the solver runs.
    """
    if isinstance(solver, str):
        try:
            solver = _SOLVERS[solver.lower()]
        except KeyError:
            raise ValueError(
                f"Unknown solver {solver!r}; expected one of {', '.join(_SOLVERS)}") from None
    if time_limit is None:
        time_limit = float("inf")
    return _solve_with(model, solver, float(time_limit), bool(verbose))
