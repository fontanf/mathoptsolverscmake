# Model-building benchmark: results

Results of the benchmark in this directory (see [README.md](README.md) for the
model, the tools and how to run it), run on 2026-11-24.

## Setup

| | |
|---|---|
| Machine | Intel Core Ultra 5 135U (14 logical CPUs), 15 GB RAM, laptop with other applications running |
| OS | Ubuntu 24.04.5 LTS, Linux 7.0.0 |
| Python | 3.12.3: mathoptmodel 0.2.0, Pyomo 6.11.0, pyoptinterface 0.6.1 (HiGHS 1.15.1 from highsbox), NumPy 2.5.5, SciPy 1.18.1 |
| Java | OpenJDK 25.0.5 (Temurin): mathoptmodel 0.2.0, OR-Tools 9.15.6755 |

The model is a random sparse MILP with `n` variables (half integer), `n`
constraints and 5 nonzeros per constraint; `n` = 10,000, 100,000 and
1,000,000. Times are the best of 5, 3 and 2 repetitions respectively, and each
size runs in its own process, capped at 4 GB. The LP files written for
`n` = 10,000 were read back with mathoptmodel: all the tools wrote the correct
model.

With other applications running on the machine, differences under about 20%
are not significant.

## Build time

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| Python | mathoptmodel, one batch | 0.2 ms | 1.6 ms | **37 ms** |
| Python | mathoptmodel, one call per constraint | 17 ms | 166 ms | 2.12 s |
| Python | Pyomo | 95 ms | 1.30 s | 14.71 s |
| Python | pyoptinterface, one call per constraint | 127 ms | 1.27 s | 13.26 s |
| Python | pyoptinterface, matrix API | 146 ms | 1.49 s | 15.98 s |
| Java | mathoptmodel, one batch | 0.6 ms | 2.1 ms | **38 ms** |
| Java | mathoptmodel, one call per constraint | 3.8 ms | 18 ms | 247 ms |
| Java | OR-Tools ModelBuilder | 9.1 ms | 107 ms | 1.01 s |
| Java | OR-Tools MPSolver | 7.6 ms | 118 ms | 1.90 s |

- Building in batches takes 37–38 ms at `n` = 1,000,000 in both languages,
  essentially copying the arrays: about 400 times faster than Pyomo,
  360 times faster than pyoptinterface and 27 times faster than OR-Tools
  ModelBuilder.
- One call per constraint stays faster than the other tools:
  - in Java, 4 times faster than OR-Tools ModelBuilder;
  - in Python, 7 times faster than Pyomo or pyoptinterface, although each
    call costs about 2 µs, mostly to convert its arguments.
- pyoptinterface's matrix API (`add_m_linear_constraints`) loops over the rows
  in Python, so it's no faster than its per-constraint API.

## LP write time

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| Python | mathoptmodel 0.2.0 | 18 ms | 169 ms | 2.47 s |
| Python | Pyomo | 91 ms | 1.39 s | 18.87 s |
| Python | pyoptinterface (HiGHS's writer) | 18 ms | 217 ms | 3.73 s |
| Java | mathoptmodel 0.2.0 | 13 ms | 151 ms | 2.06 s |
| Java | OR-Tools ModelBuilder | 9.0 ms | 105 ms | 1.76 s |
| Java | OR-Tools MPSolver | 16 ms | 281 ms | 4.40 s |
| Python | mathoptmodel, new LP writer (unreleased) | 2.6 ms | 16 ms | **0.18 s** |
| Java | mathoptmodel, new LP writer (unreleased) | 1.8 ms | 19 ms | **0.17 s** |

With 0.2.0's LP writer, writing the file took 60 times longer than building
the model. The new writers are 5 to 12 times faster for this model and write
the same bytes; see "The new writers" below.

## Peak memory

Peak resident memory of the whole process (the benchmark's input arrays,
Python or the JVM included):

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| Python | mathoptmodel, one batch | 30 MB | 75 MB | **422 MB** |
| Python | mathoptmodel, one call per constraint | 45 MB | 112 MB | 835 MB |
| Python | Pyomo | 70 MB | 283 MB | 2,115 MB |
| Python | pyoptinterface, one call per constraint | 61 MB | 153 MB | 1,185 MB |
| Python | pyoptinterface, matrix API | 63 MB | 126 MB | 722 MB |
| Java | mathoptmodel, one batch | 80 MB | 147 MB | **767 MB** |
| Java | mathoptmodel, one call per constraint | 92 MB | 184 MB | 944 MB |
| Java | OR-Tools ModelBuilder | 110 MB | 250 MB | 1,356 MB |
| Java | OR-Tools MPSolver | 115 MB | 333 MB | 2,236 MB |

## Caveats

- pyoptinterface builds the model directly into HiGHS, so its numbers include
  HiGHS's own storage of the model.
- OR-Tools' MPSolver is created with the SCIP backend, but the model is only
  passed to SCIP when solving, which isn't timed.
- Only linear models are compared, since OR-Tools has no nonlinear modeling.

## The new writers

Timed in C++ on the 1,000,000-variable model, with integer coefficients (as
in the benchmark) and with the coefficients divided by 7 (fractional):

| Format | Size | 0.2.0, integer | New, integer | 0.2.0, fractional | New, fractional |
|---|---|---|---|---|---|
| LP | 98 MB | 2.16 s | **0.18 s** | 2.25 s | 0.78 s |
| MPS | 164 MB | 3.29 s | **0.60 s** | 3.19 s | 1.05 s |
| `.nl` | 94 MB | 2.17 s | **0.45 s** | 2.23 s | 1.08 s |
| JSON | 101 MB | 1.82 s | **0.19 s** | 1.76 s | 0.75 s |

(The fractional files are 70–75 MB larger.)

The four writers now share a buffered output (`src/text_output.hpp`):

- **Buffered output:** the text is built in a 1 MB buffer and written in large
  chunks, instead of one `std::ostream` insertion per token.
- **Number formatting:**
  - integers, and numbers equal to an integer (most coefficients and bounds in
    practice), are formatted directly;
  - other numbers use `snprintf("%.17g")`, which gives the same text as
    `std::setprecision(17)`, and which is now most of the time with
    fractional data.
- **Names:** default names (`x<i>`, `c<i>`) are written directly, without
  building a string per variable and constraint first.
- **Per-format changes:**
  - **LP:** the name check that runs the LP tokenizer now only runs for names
    that aren't plain identifiers.
  - **MPS:**
    - the matrix is transposed into flat arrays, rather than one vector per
      column holding copies of the row names;
    - the objective row name's clash check scans the constraint names, rather
      than putting them all in an `std::set`.
  - **JSON:** variable type names are computed once per type, rather than
    with an `std::ostringstream` per variable.
  - **`.nl`:** the J segments reuse one vector across constraints.
- **Direct file writing:** the Python and Java bindings write uncompressed
  files with the C++ writers directly, instead of building the whole file as
  a string, then copying it to Python or Java and writing it. In the
  benchmark, writing the LP file at `n` = 1,000,000 went from 2.47 s to
  0.18 s in Python, and from 2.06 s to 0.17 s in Java.

The output of all four writers is byte-identical to 0.2.0's:
- on the benchmark model, with integer and fractional coefficients;
- on 1,883 instances (Netlib and set-covering MPS files, and MINLPLib 2024
  `.nl` files), each written in the four formats: all 7,512 files are
  identical. The 3,165 writes that fail, mostly nonlinear models written as
  LP or MPS, fail with the same errors.
