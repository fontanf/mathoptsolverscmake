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

## Building from source

```shell
pip install .
pytest python/tests
```
