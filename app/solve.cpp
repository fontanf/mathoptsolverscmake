#include "mathoptsolverscmake/mathopt.hpp"
#include "mathoptsolverscmake/dual.hpp"
#include "mathoptsolverscmake/mathopt_json.hpp"
#include "mathoptsolverscmake/mathopt_lp.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"
#include "mathoptsolverscmake/mathopt_nl.hpp"

#ifdef CBC_FOUND
#include "mathoptsolverscmake/mathopt_cbc.hpp"
#endif
#ifdef HIGHS_FOUND
#include "mathoptsolverscmake/mathopt_highs.hpp"
#endif
#ifdef XPRESS_FOUND
#include "mathoptsolverscmake/mathopt_xpress.hpp"
#endif
#ifdef KNITRO_FOUND
#include "mathoptsolverscmake/mathopt_knitro.hpp"
#endif

#include <boost/program_options.hpp>

#include <cctype>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>

using namespace mathoptsolverscmake;
namespace po = boost::program_options;

namespace
{

enum class FileFormat
{
    Mps,
    Lp,
    Json,
    Nl,
};

std::istream& operator>>(std::istream& in, FileFormat& format)
{
    std::string token;
    in >> token;
    for (char& c: token)
        c = (char)std::tolower((unsigned char)c);
    if (token == "mps") {
        format = FileFormat::Mps;
    } else if (token == "lp") {
        format = FileFormat::Lp;
    } else if (token == "json") {
        format = FileFormat::Json;
    } else if (token == "nl") {
        format = FileFormat::Nl;
    } else {
        in.setstate(std::ios_base::failbit);
    }
    return in;
}

FileFormat format_from_extension(const std::string& file_path)
{
    std::string extension = std::filesystem::path(file_path).extension().string();
    for (char& c: extension)
        c = (char)std::tolower((unsigned char)c);
    if (extension == ".mps") return FileFormat::Mps;
    if (extension == ".lp") return FileFormat::Lp;
    if (extension == ".json") return FileFormat::Json;
    if (extension == ".nl") return FileFormat::Nl;
    throw std::invalid_argument(
            "unable to guess the file format from the extension of \"" + file_path + "\"; "
            "pass --format explicitly (mps, lp, json or nl).");
}

MathOptModel read_model(const std::string& file_path, FileFormat format)
{
    switch (format) {
    case FileFormat::Mps: return read_mps(file_path);
    case FileFormat::Lp: return read_lp(file_path);
    case FileFormat::Json: return read_json(file_path);
    case FileFormat::Nl: return read_nl(file_path);
    }
    throw std::logic_error("unreachable");
}

/** The solver-specific dispatch; only Cbc, HiGHS, XPRESS and Knitro have a
 * generic load()/solve()/get_solution() interface suited to "load an
 * arbitrary file and solve it" (dlib/ConicBundle/SMS++ target narrower use
 * cases -- box-constrained, convex Lagrangian relaxation, black-box-only --
 * meant for embedding in a larger algorithm, not standalone solving). */
struct SolveResult
{
    bool found_solution = false;
    std::vector<double> solution;
    double objective_value = std::numeric_limits<double>::quiet_NaN();
    bool has_bound = false;
    double bound = std::numeric_limits<double>::quiet_NaN();
};

SolveResult solve_model(
        const MathOptModel& model,
        SolverName solver,
        bool has_time_limit,
        double time_limit)
{
    SolveResult result;

    switch (solver) {
    case SolverName::Cbc: {
#ifdef CBC_FOUND
        OsiCbcSolverInterface osi_solver;
        CbcModel cbc_model(osi_solver);
        if (has_time_limit)
            set_time_limit(cbc_model, time_limit);
        load(cbc_model, model);
        solve(cbc_model);
        result.solution = get_solution(cbc_model);
        result.found_solution = !result.solution.empty();
        result.bound = get_bound(cbc_model);
        result.has_bound = true;
#else
        throw std::invalid_argument("solver \"Cbc\" is not available in this build.");
#endif
        break;
    } case SolverName::Highs: {
#ifdef HIGHS_FOUND
        Highs highs;
        if (has_time_limit)
            set_time_limit(highs, time_limit);
        load(highs, model);
        solve(highs);
        result.solution = get_solution(highs);
        result.found_solution = !result.solution.empty();
        result.bound = get_bound(highs);
        result.has_bound = true;
#else
        throw std::invalid_argument("solver \"HiGHS\" is not available in this build.");
#endif
        break;
    } case SolverName::Xpress: {
#ifdef XPRESS_FOUND
        XPRSprob xpress_model;
        XPRScreateprob(&xpress_model);
        if (has_time_limit)
            set_time_limit(xpress_model, time_limit);
        load(xpress_model, model);
        solve(xpress_model);
        // Unlike Cbc/HiGHS, XPRESS's get_solution() throws rather than
        // returning empty when no solution is available.
        try {
            result.solution = get_solution(xpress_model);
            result.found_solution = true;
            result.bound = get_bound(xpress_model);
            result.has_bound = true;
        } catch (const std::exception&) {
            result.found_solution = false;
        }
        XPRSdestroyprob(xpress_model);
#else
        throw std::invalid_argument("solver \"XPRESS\" is not available in this build.");
#endif
        break;
    } case SolverName::Knitro: {
#ifdef KNITRO_FOUND
        KnitroContext knitro;
        if (has_time_limit)
            set_time_limit(knitro, time_limit);
        load(knitro, model);
        solve(knitro);
        if (!is_infeasible(knitro)) {
            result.solution = get_solution(knitro);
            result.objective_value = get_solution_value(knitro);
            result.found_solution = true;
        }
#else
        throw std::invalid_argument("solver \"Knitro\" is not available in this build.");
#endif
        break;
    } default: {
        std::ostringstream solver_name;
        solver_name << solver;
        throw std::invalid_argument(
                "solver \"" + solver_name.str() + "\" is not supported by this tool "
                "(only Cbc, HiGHS, XPRESS and Knitro are).");
    }
    }

    if (result.found_solution && solver != SolverName::Knitro)
        result.objective_value = model.evaluate_objective(result.solution);

    return result;
}

}

int main(int argc, char* argv[])
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i", po::value<std::string>(), "input model file (required)")
        ("format,f", po::value<FileFormat>(), "input file format: mps, lp, json or nl (guessed from the file extension if not given)")
        ("solver,s", po::value<SolverName>(), "solver: cbc, highs, xpress or knitro (required)")
        ("dual,d", po::bool_switch(), "solve the dual of the model instead of the model itself (LPs only, see mathoptsolverscmake::dual()); --output then writes the dual solution")
        ("time-limit,t", po::value<double>(), "time limit in seconds (no limit if not given)")
        ("output,o", po::value<std::string>(), "solution output file (see MathOptModel::write_solution())")
        ("verbosity-level,v", po::value<int>()->default_value(1), "model summary verbosity level (0: none)")
        ("tolerance", po::value<double>()->default_value(1e-6), "feasibility/integrality tolerance used when checking the returned solution (real solvers never return exactly-integral/exactly-on-bound values)");

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
        po::notify(vm);
    } catch (const po::error& e) {
        std::cerr << "Error: " << e.what() << "\n\n" << desc << std::endl;
        return 1;
    }

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }
    if (!vm.count("input")) {
        std::cerr << "Error: --input is required.\n\n" << desc << std::endl;
        return 1;
    }
    if (!vm.count("solver")) {
        std::cerr << "Error: --solver is required.\n\n" << desc << std::endl;
        return 1;
    }

    std::string input_path = vm["input"].as<std::string>();
    SolverName solver = vm["solver"].as<SolverName>();
    bool has_time_limit = vm.count("time-limit") > 0;
    double time_limit = has_time_limit? vm["time-limit"].as<double>(): 0.0;

    MathOptModel model;
    try {
        FileFormat format = vm.count("format")?
            vm["format"].as<FileFormat>():
            format_from_extension(input_path);
        model = read_model(input_path, format);
    } catch (const std::exception& e) {
        std::cerr << "Error reading \"" << input_path << "\": " << e.what() << std::endl;
        return 1;
    }
    if (vm["dual"].as<bool>()) {
        try {
            model = dual(model);
        } catch (const std::exception& e) {
            std::cerr << "Error building the dual of \"" << input_path << "\": " << e.what() << std::endl;
            return 1;
        }
    }
    model.feasibility_tolerance = vm["tolerance"].as<double>();
    model.integrality_tolerance = vm["tolerance"].as<double>();

    int verbosity_level = vm["verbosity-level"].as<int>();
    if (verbosity_level > 0) {
        std::cout
            << "===========================" << std::endl
            << "          MathOpt          " << std::endl
            << "===========================" << std::endl
            << std::endl
            << "Model" << std::endl
            << "-----" << std::endl;
        model.format(std::cout, verbosity_level);
        std::cout << std::endl;
    }

    SolveResult result;
    try {
        result = solve_model(model, solver, has_time_limit, time_limit);
    } catch (const std::exception& e) {
        std::cerr << "Error solving the model: " << e.what() << std::endl;
        return 1;
    }

    if (!result.found_solution) {
        std::cout << "No feasible solution found "
            "(the model may be infeasible, or the solver may have stopped "
            "before finding one)." << std::endl;
        if (result.has_bound)
            std::cout << "Bound: " << result.bound << std::endl;
        return 1;
    }

    if (verbosity_level > 0) {
        std::cout
            << std::endl
            << "Final statistics" << std::endl
            << "----------------" << std::endl;
        std::cout << "Objective value:  " << result.objective_value << std::endl;
        if (result.has_bound)
            std::cout << "Bound:            " << result.bound << std::endl;
        std::cout << "Feasible:         " << model.check_solution(result.solution) << std::endl;
    }

    if (vm.count("output"))
        model.write_solution(result.solution, vm["output"].as<std::string>());

    return 0;
}
