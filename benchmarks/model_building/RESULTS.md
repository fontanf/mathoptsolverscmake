# Model-building benchmark: results

Results of the benchmark in this directory (see [README.md](README.md) for the
model, the tools and how to run it), run on 2026-11-25 with mathoptmodel
0.2.1.

## Setup

| | |
|---|---|
| Machine | Intel Core Ultra 5 135U (14 logical CPUs), 15 GB RAM, laptop with other applications running |
| OS | Ubuntu 24.04.5 LTS, Linux 7.0.0 |
| Python | 3.12.3: mathoptmodel 0.2.1, Pyomo 6.11.0, pyoptinterface 0.6.1 (HiGHS 1.15.1 from highsbox), NumPy 2.5.5, SciPy 1.18.1 |
| Java | OpenJDK 25.0.5 (Temurin): mathoptmodel 0.2.1, OR-Tools 9.15.6755 |

The model is a random sparse MILP with `n` variables (half integer), `n`
constraints and 5 nonzeros per constraint; `n` = 10,000, 100,000 and
1,000,000. Times are the best of 5, 3 and 2 repetitions respectively, and each
size runs in its own process, capped at 4 GB. The LP files written for
`n` = 10,000 were read back with mathoptmodel: all the tools wrote the correct
model.

With other applications running on the machine, differences under about 20%
are not significant. For example, building in one batch at
`n` = 1,000,000 took between 31 and 50 ms over several runs.

For Pyomo and OR-Tools, the benchmark has two kinds of rows:

- the tool's usual API, as its documentation shows it: Pyomo's `sum(...)`
  expressions, OR-Tools' `ModelBuilder` and `MPSolver`;
- the fastest way we found to build the same model with the tool:
  - Pyomo's `LinearExpression`, which takes the coefficients and the
    variables directly;
  - OR-Tools' `ModelBuilderHelper`, the low-level layer under `ModelBuilder`,
    with one native call per variable attribute and per term;
  - an OR-Tools `MPModelProto`, built in Java and then loaded into `MPSolver`.

pyoptinterface's per-constraint API is already its fastest: most of its time
goes into adding the integer variables to HiGHS (about 20 µs each, against
1 µs for a continuous variable), with its vector API `add_m_variables` too.

## Build time

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| Python | mathoptmodel, one batch | 0.2 ms | 3.2 ms | **41 ms** |
| Python | mathoptmodel, one call per constraint | 17 ms | 172 ms | 1.73 s |
| Python | Pyomo | 95 ms | 1.31 s | 14.72 s |
| Python | Pyomo, `LinearExpression` | 59 ms | 861 ms | 9.94 s |
| Python | pyoptinterface, one call per constraint | 128 ms | 1.29 s | 13.50 s |
| Python | pyoptinterface, matrix API | 146 ms | 1.50 s | 16.22 s |
| Java | mathoptmodel, one batch | 0.5 ms | 4.2 ms | **51 ms** |
| Java | mathoptmodel, one call per constraint | 6.6 ms | 19 ms | 287 ms |
| Java | OR-Tools ModelBuilder | 10 ms | 87 ms | 1.00 s |
| Java | OR-Tools ModelBuilder, helper | 4.8 ms | 37 ms | 262 ms |
| Java | OR-Tools MPSolver | 8.6 ms | 123 ms | 1.94 s |
| Java | OR-Tools MPSolver, `MPModelProto` | 10 ms | 114 ms | 1.52 s |

- Building in batches takes 31–51 ms at `n` = 1,000,000 in both languages,
  essentially copying the arrays:
  - in Python, 240 times faster than Pyomo with `LinearExpression`, and
    330 times faster than pyoptinterface;
  - in Java, 5 to 8 times faster than OR-Tools' `ModelBuilderHelper`, and
    20 to 30 times faster than its `ModelBuilder` API.
- One call per constraint:
  - in Python, 6 to 8 times faster than Pyomo with `LinearExpression` or
    pyoptinterface, although each call costs about 2 µs, mostly to convert
    its arguments;
  - in Java, about as fast as OR-Tools' `ModelBuilderHelper`.
- Pyomo's `LinearExpression` builds 1.5 times faster than `sum(...)`, and
  OR-Tools' `ModelBuilderHelper` 4 times faster than the `ModelBuilder` API.
  Loading an `MPModelProto` into `MPSolver` is only about 25% faster than
  adding the constraints one by one.
- pyoptinterface's matrix API (`add_m_linear_constraints`) loops over the rows
  in Python, so it's no faster than its per-constraint API.

## LP write time

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| Python | mathoptmodel | 1.9 ms | 17 ms | **0.18 s** |
| Python | Pyomo | 93 ms | 1.43 s | 19.07 s |
| Python | Pyomo, `LinearExpression` | 91 ms | 1.39 s | 19.32 s |
| Python | pyoptinterface (HiGHS's writer) | 19 ms | 223 ms | 3.91 s |
| Java | mathoptmodel | 1.8 ms | 19 ms | **0.19 s** |
| Java | OR-Tools ModelBuilder | 10 ms | 107 ms | 1.83 s |
| Java | OR-Tools ModelBuilder, helper | 8.3 ms | 104 ms | 1.74 s |
| Java | OR-Tools MPSolver | 16 ms | 279 ms | 4.24 s |
| Java | OR-Tools MPSolver, `MPModelProto` | 16 ms | 254 ms | 4.24 s |
| Python | mathoptmodel 0.2.0 | 18 ms | 169 ms | 2.47 s |
| Java | mathoptmodel 0.2.0 | 13 ms | 151 ms | 2.06 s |

(The mathoptmodel rows are the batch-built model; the per-constraint one
writes the same file in the same time. pyoptinterface's row is the
matrix-built model, its per-constraint one is within 5%.)

The last two rows are 0.2.0's writer, measured on 2026-11-24 for the same
model: 0.2.1's is 11 to 14 times faster here, and now writes the file about 10
times faster than OR-Tools and 20 times faster than pyoptinterface (HiGHS),
4 times longer than building the model in one batch; see "The new writers"
below.

## Peak memory

Peak resident memory of the whole process (the benchmark's input arrays,
Python or the JVM included), writing the LP file included:

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| Python | mathoptmodel, one batch | 29 MB | 46 MB | **213 MB** |
| Python | mathoptmodel, one call per constraint | 45 MB | 85 MB | 600 MB |
| Python | Pyomo | 71 MB | 283 MB | 2,115 MB |
| Python | Pyomo, `LinearExpression` | 71 MB | 280 MB | 2,153 MB |
| Python | pyoptinterface, one call per constraint | 71 MB | 153 MB | 1,184 MB |
| Python | pyoptinterface, matrix API | 71 MB | 127 MB | 722 MB |
| Java | mathoptmodel, one batch | 77 MB | 103 MB | **268 MB** |
| Java | mathoptmodel, one call per constraint | 91 MB | 149 MB | 420 MB |
| Java | OR-Tools ModelBuilder | 110 MB | 239 MB | 1,280 MB |
| Java | OR-Tools ModelBuilder, helper | 95 MB | 172 MB | 876 MB |
| Java | OR-Tools MPSolver | 115 MB | 332 MB | 2,193 MB |
| Java | OR-Tools MPSolver, `MPModelProto` | 135 MB | 380 MB | 2,339 MB |

With 0.2.0, the bindings built the whole LP file as a string before writing
it: the peak at `n` = 1,000,000 was 422 MB in Python and 767 MB in Java for
the batch-built model. 0.2.1 writes the file directly, which halves it in
Python and divides it by 3 in Java.

## Caveats

- pyoptinterface builds the model directly into HiGHS, so its numbers include
  HiGHS's own storage of the model.
- OR-Tools' MPSolver is created with the SCIP backend, but the model is only
  passed to SCIP when solving, which isn't timed.
- Only linear models are compared, since OR-Tools has no nonlinear modeling.

## The new writers

0.2.1 makes the four writers faster. Timed in C++ on the 1,000,000-variable
model, with integer coefficients (as in the benchmark) and with the
coefficients divided by 7 (fractional):

| Format | Size | 0.2.0, integer | 0.2.1, integer | 0.2.0, fractional | 0.2.1, fractional |
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
  a string, then copying it to Python or Java and writing it.

The output of all four writers is byte-identical to 0.2.0's:
- on the benchmark model, with integer and fractional coefficients;
- on 1,883 instances (Netlib and set-covering MPS files, and MINLPLib 2024
  `.nl` files), each written in the four formats: all 7,512 files are
  identical. The 3,165 writes that fail, mostly nonlinear models written as
  LP or MPS, fail with the same errors.
