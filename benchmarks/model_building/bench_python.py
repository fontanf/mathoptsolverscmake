"""Time building (and writing as LP) the benchmark model with a Python tool.

    python bench_python.py TOOL N REPETITIONS OUTPUT_DIRECTORY

TOOL is one of TOOLS below. Prints one JSON line with the best build and
write times over the repetitions and the process's peak memory. The model
data is generated before the timings (see data.py); for the tools that add
constraints one at a time, converting the data to the tool's input types is
done beforehand too, so that only the tool's own work is timed.
"""
import gc
import json
import os
import resource
import sys
import time

import numpy as np

import data


def mathoptmodel_batch(n, obj, cols, coefs, rhs):
    import mathoptmodel as mom

    def build():
        model = mom.Model()
        types = np.zeros(n, dtype=np.int8)
        types[:data.number_of_integers(n)] = int(mom.VariableType.Integer)
        model.add_variables(lower_bounds=0.0, upper_bounds=data.UPPER_BOUND, types=types,
                            objective_coefficients=obj)
        model.objective_direction = mom.ObjectiveDirection.Maximize
        model.add_constraints(upper_bounds=rhs,
                              starts=np.arange(0, n * data.K, data.K, dtype=np.int32),
                              variables=cols.reshape(-1), coefficients=coefs.reshape(-1))
        return model

    def write(model, path):
        mom.write(model, path)

    return build, write


def mathoptmodel_row(n, obj, cols, coefs, rhs):
    import mathoptmodel as mom

    starts = np.zeros(1, dtype=np.int32)
    rows = [(cols[i], coefs[i]) for i in range(n)]
    rhs_list = rhs.tolist()

    def build():
        model = mom.Model()
        integers = data.number_of_integers(n)
        model.add_variables(integers, lower_bounds=0.0, upper_bounds=data.UPPER_BOUND,
                            types=mom.VariableType.Integer)
        model.add_variables(n - integers, lower_bounds=0.0, upper_bounds=data.UPPER_BOUND)
        model.set_objective(mom.ObjectiveDirection.Maximize, coefficients=obj)
        for (row_cols, row_coefs), row_rhs in zip(rows, rhs_list):
            model.add_constraints(upper_bounds=row_rhs, starts=starts,
                                  variables=row_cols, coefficients=row_coefs)
        return model

    def write(model, path):
        mom.write(model, path)

    return build, write


def pyomo(n, obj, cols, coefs, rhs):
    import pyomo.environ as pyo

    integers = data.number_of_integers(n)
    cols_list = cols.tolist()
    coefs_list = coefs.tolist()
    rhs_list = rhs.tolist()
    obj_list = obj.tolist()

    def build():
        model = pyo.ConcreteModel()
        model.x = pyo.Var(
            range(n), bounds=(0, data.UPPER_BOUND),
            within=lambda model, j: pyo.Integers if j < integers else pyo.Reals)
        x = model.x
        model.c = pyo.Constraint(
            range(n),
            rule=lambda model, i: sum(c * x[j] for c, j in zip(coefs_list[i], cols_list[i])) <= rhs_list[i])
        model.obj = pyo.Objective(expr=sum(c * x[j] for j, c in enumerate(obj_list)), sense=pyo.maximize)
        return model

    def write(model, path):
        model.write(path, format="lp", io_options={"symbolic_solver_labels": False})

    return build, write


def _pyoptinterface_model():
    from pyoptinterface import highs

    if not highs.is_library_loaded() and not highs.autoload_library():
        raise RuntimeError("HiGHS library not found (pip install highsbox)")
    model = highs.Model()
    model.set_model_attribute(__import__("pyoptinterface").ModelAttribute.Silent, True)
    return model


def pyoptinterface_row(n, obj, cols, coefs, rhs):
    import pyoptinterface as poi

    integers = data.number_of_integers(n)
    cols_list = cols.tolist()
    coefs_list = coefs.tolist()
    rhs_list = rhs.tolist()
    obj_list = obj.tolist()

    def build():
        model = _pyoptinterface_model()
        x = [model.add_variable(lb=0, ub=data.UPPER_BOUND,
                                domain=poi.VariableDomain.Integer if j < integers else poi.VariableDomain.Continuous)
             for j in range(n)]
        # ScalarAffineFunction's fast constructor takes the variables' indices.
        index = [v.index for v in x]
        for row_cols, row_coefs, row_rhs in zip(cols_list, coefs_list, rhs_list):
            model.add_linear_constraint(
                poi.ScalarAffineFunction(row_coefs, [index[j] for j in row_cols]), poi.Leq, row_rhs)
        model.set_objective(poi.ScalarAffineFunction(obj_list, index), poi.ObjectiveSense.Maximize)
        return model

    def write(model, path):
        model.write(path)

    return build, write


def pyoptinterface_matrix(n, obj, cols, coefs, rhs):
    import pyoptinterface as poi
    from scipy.sparse import csr_array

    integers = data.number_of_integers(n)
    matrix = csr_array((coefs.reshape(-1), cols.reshape(-1),
                        np.arange(0, n * data.K + 1, data.K)), shape=(n, n))
    obj_list = obj.tolist()

    def build():
        model = _pyoptinterface_model()
        x_integer = model.add_m_variables(integers, domain=poi.VariableDomain.Integer,
                                          lb=0, ub=data.UPPER_BOUND)
        x_continuous = model.add_m_variables(n - integers, lb=0, ub=data.UPPER_BOUND)
        x = np.concatenate([x_integer, x_continuous])
        model.add_m_linear_constraints(matrix, x, poi.Leq, rhs)
        model.set_objective(poi.ScalarAffineFunction(obj_list, [v.index for v in x]),
                            poi.ObjectiveSense.Maximize)
        return model

    def write(model, path):
        model.write(path)

    return build, write


TOOLS = {
    "mathoptmodel-batch": mathoptmodel_batch,
    "mathoptmodel-row": mathoptmodel_row,
    "pyomo": pyomo,
    "pyoptinterface-row": pyoptinterface_row,
    "pyoptinterface-matrix": pyoptinterface_matrix,
}


def main():
    tool, n, repetitions, output_directory = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    obj, cols, coefs, rhs = data.generate(n)
    build, write = TOOLS[tool](n, obj, cols, coefs, rhs)
    path = os.path.join(output_directory, f"{tool}-{n}.lp")
    build_times, write_times = [], []
    for _ in range(repetitions):
        gc.collect()
        start = time.perf_counter()
        model = build()
        build_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        write(model, path)
        write_times.append(time.perf_counter() - start)
        del model
    print(json.dumps({
        "language": "python",
        "tool": tool,
        "n": n,
        "build": min(build_times),
        "write": min(write_times),
        "peak_memory_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
        "lp": path,
    }))


if __name__ == "__main__":
    main()
