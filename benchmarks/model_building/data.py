"""Benchmark model data, generated identically in Python and Java (Bench.java).

The model is a random sparse MILP with n variables, n constraints and K
nonzeros per constraint:

    max  sum_j obj[j] x[j]
    s.t. sum_k coef[i, k] x[cols[i, k]] <= rhs[i]     for each constraint i
         0 <= x[j] <= 10, x[j] integer for j < n / 2, continuous otherwise

The values come from splitmix64 hashes of (stream, index), so that both
languages build exactly the same model without sharing data files. The K
columns of a row are distinct: cols[i, k] = (base[i] + k * (n // K)) % n.
"""
import numpy as np

K = 5
UPPER_BOUND = 10.0

_GOLDEN = np.uint64(0x9E3779B97F4A7C15)
_MIX_1 = np.uint64(0xBF58476D1CE4E5B9)
_MIX_2 = np.uint64(0x94D049BB133111EB)


def _hash(stream, indices):
    """splitmix64((stream << 40) + index), for an array of indices."""
    z = (np.uint64(stream) << np.uint64(40)) + indices.astype(np.uint64) + _GOLDEN
    z = (z ^ (z >> np.uint64(30))) * _MIX_1
    z = (z ^ (z >> np.uint64(27))) * _MIX_2
    return z ^ (z >> np.uint64(31))


def generate(n):
    """Return (obj, cols, coefs, rhs): obj (n,), cols (n, K) int32,
    coefs (n, K) and rhs (n,), all as NumPy arrays."""
    m = n
    rows = np.arange(m, dtype=np.uint64)
    base = (_hash(1, rows) % np.uint64(n)).astype(np.int64)
    offsets = np.arange(K, dtype=np.int64) * (n // K)
    cols = ((base[:, None] + offsets[None, :]) % n).astype(np.int32)
    coefs = (1 + _hash(2, np.arange(m * K, dtype=np.uint64)) % np.uint64(9)).astype(np.float64).reshape(m, K)
    rhs = (10 + _hash(3, rows) % np.uint64(90)).astype(np.float64)
    obj = (1 + _hash(4, np.arange(n, dtype=np.uint64)) % np.uint64(9)).astype(np.float64)
    return obj, cols, coefs, rhs


def number_of_integers(n):
    return n // 2
