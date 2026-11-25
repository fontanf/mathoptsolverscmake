"""Run the model-building benchmark: every tool, every size, each in its own
process under a memory cap, then check the models written and print a table.

    python run.py [--sizes 10000 100000 1000000] [--tools ...] [--memory 4G]

Run it with a Python environment that has the Python tools installed (see
README.md); the Java tools need java/target/classpath.txt, made by Maven.
Results are written to results/results.json and results/results.md.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np

import data

HERE = os.path.dirname(os.path.abspath(__file__))

PYTHON_TOOLS = ["mathoptmodel-batch", "mathoptmodel-row", "pyomo", "pyomo-linear-expression", "pyoptinterface-row", "pyoptinterface-matrix"]
JAVA_TOOLS = ["mathoptmodel-batch", "mathoptmodel-row", "ortools-modelbuilder", "ortools-modelbuilder-helper",
              "ortools-mpsolver", "ortools-mpsolver-proto"]


def repetitions(n):
    return 5 if n <= 10_000 else 3 if n <= 100_000 else 2


def command(language, tool, n, output, memory):
    if language == "python":
        run = [sys.executable, os.path.join(HERE, "bench_python.py"), tool, str(n), str(repetitions(n)), output]
    else:
        classpath_file = os.path.join(HERE, "java", "target", "classpath.txt")
        with open(classpath_file) as f:
            classpath = os.path.join(HERE, "java", "target", "classes") + os.pathsep + f.read().strip()
        java = os.path.join(os.environ["JAVA_HOME"], "bin", "java") if "JAVA_HOME" in os.environ else "java"
        run = [java, "-Xmx3g", "--enable-native-access=ALL-UNNAMED", "-cp", classpath, "bench.Bench",
               tool, str(n), str(repetitions(n)), output]
    if shutil.which("systemd-run"):
        # A cgroup memory cap, which also covers the JVM's native memory.
        run = ["systemd-run", "--user", "--scope", "--quiet", "-p", f"MemoryMax={memory}",
               "-p", "MemorySwapMax=0"] + run
    return run


def expected_signature(n):
    obj, cols, coefs, rhs = data.generate(n)
    return {
        "variables": n,
        "constraints": n,
        "elements": n * data.K,
        "coefficients": float(coefs.sum()),
        "objective": float(obj.sum()),
        "rhs": float(rhs.sum()),
        "integers": data.number_of_integers(n),
    }


def signature(path):
    """Read a written LP file back with mathoptmodel and summarize it."""
    import mathoptmodel as mom

    model = mom.read(path)
    types = model.variables_types
    # Some writers add a constant variable or an objective-constant row:
    # only count the constraints that have terms.
    starts = model.constraints_starts
    ends = np.append(starts[1:], model.number_of_elements)
    with_terms = ends > starts
    upper_bounds = model.constraints_upper_bounds
    return {
        "variables": model.number_of_variables,
        "constraints": int(with_terms.sum()),
        "elements": model.number_of_elements,
        "coefficients": float(np.abs(model.elements_coefficients).sum()),
        "objective": float(np.abs(model.objective_coefficients).sum()),
        "rhs": float(np.abs(upper_bounds[with_terms & np.isfinite(upper_bounds)]).sum()),
        "integers": int((types != int(mom.VariableType.Continuous)).sum()),
    }


def check(path, n):
    try:
        actual = signature(path)
    except Exception as e:
        return f"unreadable ({str(e).splitlines()[0][:80]})"
    expected = expected_signature(n)
    differences = [key for key in expected if not np.isclose(actual[key], expected[key])]
    return "ok" if not differences else "mismatch: " + ", ".join(
        f"{key} {actual[key]} != {expected[key]}" for key in differences)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", type=int, nargs="+", default=[10_000, 100_000, 1_000_000])
    parser.add_argument("--languages", nargs="+", default=["python", "java"])
    parser.add_argument("--tools", nargs="+", help="only these tools (default: all)")
    parser.add_argument("--memory", default="4G", help="memory cap of each run (systemd-run)")
    parser.add_argument("--timeout", type=int, default=1200, help="time limit of each run, in seconds")
    parser.add_argument("--check-size", type=int, default=10_000,
                        help="read back the LP files written for this size")
    args = parser.parse_args()

    results_directory = os.path.join(HERE, "results")
    lp_directory = os.path.join(results_directory, "lp")
    os.makedirs(lp_directory, exist_ok=True)
    results = []
    for language in args.languages:
        tools = PYTHON_TOOLS if language == "python" else JAVA_TOOLS
        for tool in tools:
            if args.tools and tool not in args.tools and f"{language}:{tool}" not in args.tools:
                continue
            for n in args.sizes:
                label = f"{language} {tool} n={n}"
                print(f"{label}...", file=sys.stderr, flush=True)
                try:
                    process = subprocess.run(command(language, tool, n, lp_directory, args.memory),
                                             capture_output=True, text=True, timeout=args.timeout, cwd=HERE)
                except subprocess.TimeoutExpired:
                    results.append({"language": language, "tool": tool, "n": n, "error": "timeout"})
                    continue
                lines = [line for line in process.stdout.splitlines() if line.startswith("{")]
                if process.returncode != 0 or not lines:
                    error = "out of memory" if process.returncode in (-9, 137) else \
                        (process.stderr.strip().splitlines() or ["failed"])[-1][:120]
                    results.append({"language": language, "tool": tool, "n": n, "error": error})
                    print(f"  {error}", file=sys.stderr)
                    continue
                result = json.loads(lines[-1])
                if n == args.check_size:
                    result["check"] = check(result["lp"], n)
                lp = result.pop("lp")
                if n != args.check_size and os.path.exists(lp):
                    os.remove(lp)
                results.append(result)
                print(f"  build {result['build']:.3f} s, write {result['write']:.3f} s, "
                      f"{result['peak_memory_mb']:.0f} MB {result.get('check', '')}", file=sys.stderr)

    with open(os.path.join(results_directory, "results.json"), "w") as f:
        json.dump(results, f, indent=1)
    table = markdown(results, args.sizes)
    with open(os.path.join(results_directory, "results.md"), "w") as f:
        f.write(table)
    print(table)


def markdown(results, sizes):
    def cell(result, key):
        if result is None:
            return ""
        if "error" in result:
            return result["error"]
        value = result[key]
        return f"{value:.0f} MB" if key == "peak_memory_mb" else \
            f"{value * 1000:.1f} ms" if value < 1 else f"{value:.2f} s"

    lines = []
    for key, title in [("build", "Build time"), ("write", "LP write time"), ("peak_memory_mb", "Peak memory")]:
        lines.append(f"### {title}\n")
        lines.append("| Language | Tool | " + " | ".join(f"n = {n:,}" for n in sizes) + " |")
        lines.append("|---|---|" + "---|" * len(sizes))
        seen = []
        for result in results:
            if (result["language"], result["tool"]) not in seen:
                seen.append((result["language"], result["tool"]))
        for language, tool in seen:
            by_size = {r["n"]: r for r in results if r["language"] == language and r["tool"] == tool}
            lines.append(f"| {language} | {tool} | "
                         + " | ".join(cell(by_size.get(n), key) for n in sizes) + " |")
        lines.append("")
    checks = [r for r in results if "check" in r]
    if checks:
        lines.append("### Models written (read back with mathoptmodel)\n")
        lines.append("| Language | Tool | Check |")
        lines.append("|---|---|---|")
        for r in checks:
            lines.append(f"| {r['language']} | {r['tool']} | {r['check']} |")
        lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":
    main()
