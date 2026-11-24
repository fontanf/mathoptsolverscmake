# mathoptmodel

`mathoptmodel` is the Python package of
[MathOptSolversCMake](https://github.com/fontanf/mathoptsolverscmake): Python
bindings for the solver-independent part of the library, its
`MathOptModel` optimization-model structure and its MPS, LP, JSON and AMPL
`.nl` readers and writers. No solver is bundled. Requires Python 3.12 or
newer (a single stable-ABI wheel per platform).

```shell
pip install mathoptmodel
```

```python
import mathoptmodel as mom

model = mom.read("p0033.mps.gz")  # format from the extension; .gz/.bz2 handled
print(model.number_of_variables, model.number_of_constraints)

# Arrays are read-only NumPy views onto the C++ model (no copy).
starts = model.constraints_starts
for i in range(model.number_of_constraints):
    begin, end = starts[i], model.constraint_end(i)
    variables = model.elements_variables[begin:end]
    coefficients = model.elements_coefficients[begin:end]
    print(model.constraint_name(i), model.constraint_class(i), dict(zip(variables, coefficients)))

mom.write(model, "p0033.lp")
```

The fields mirror `MathOptModel` in
[`mathopt.hpp`](https://github.com/fontanf/mathoptsolverscmake/blob/main/include/mathoptsolverscmake/mathopt.hpp),
which documents each of them, including the expression-tree layout of
nonlinear objectives and constraints. `variables_types` and
`constraints_classes` are int8 arrays of `VariableType` / `ConstraintClass`
values.

## Building a model

A model can also be built from Python, in batches of NumPy arrays laid out
like `MathOptModel`'s own fields, so that building it costs a few array
copies rather than one Python call per variable, constraint or term:

```python
import numpy as np
import mathoptmodel as mom

model = mom.Model()

# max x + 2 y + 3 z, x in [0, 10] continuous, y binary, z in [0, 5] integer.
model.add_variables(
    lower_bounds=0.0,                      # a scalar applies to the whole batch
    upper_bounds=np.array([10.0, 1.0, 5.0]),
    types=np.array([0, 1, 2], dtype=np.int8),  # VariableType values
    names=["x", "y", "z"])
model.set_objective(mom.ObjectiveDirection.Maximize, coefficients=np.array([1.0, 2.0, 3.0]))

# x + y + z <= 10 and x - y >= -3, in CSR form: the terms of all the rows
# concatenated, and the offset of each row's first term.
model.add_constraints(
    lower_bounds=np.array([-np.inf, -3.0]),
    upper_bounds=np.array([10.0, np.inf]),
    starts=np.array([0, 3]),
    variables=np.array([0, 1, 2, 0, 1]),
    coefficients=np.array([1.0, 1.0, 1.0, 1.0, -1.0]))

# exp(x) + y^2 <= 4, as a pre-order expression tree (see mathopt.hpp).
model.add_constraints(
    upper_bounds=4.0,
    nonlinear_starts=[0],
    nonlinear_operators=b"+evpvk",
    nonlinear_number_of_children=[2, 1, 0, 2, 0, 0],
    nonlinear_variables=[-1, -1, 0, -1, 1, -1],
    nonlinear_values=[0, 0, 0, 0, 0, 2])

assert model.check()
mom.write(model, "model.nl")
```

- `add_variables()` and `add_constraints()` return the index of the first
  variable / constraint of the batch.
- Quadratic constraint terms go through `quadratic_starts`,
  `quadratic_variables_1`, `quadratic_variables_2` and
  `quadratic_coefficients`, laid out the same way.
- Common subexpressions go through `add_defined_variables()`.
- `set_objective()` replaces only the parts it is given: the direction, the
  linear coefficients, the quadratic terms and the expression tree.
- Arrays of another dtype, non-contiguous arrays and lists are converted.
- Each batch is validated before the model is modified (`ValueError`
  otherwise), and the copies run without the GIL.

Growing a model may reallocate its arrays: views such as
`model.variables_lower_bounds` obtained before a modification must not be used
after it (get them again).

## Building from source

```shell
pip install .
pytest python/tests
```
