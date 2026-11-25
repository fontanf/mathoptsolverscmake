#include "mathoptsolverscmake/mathopt.hpp"
#include "mathoptsolverscmake/dual.hpp"
#include "mathoptsolverscmake/mathopt_json.hpp"
#include "mathoptsolverscmake/mathopt_lp.hpp"
#include "mathoptsolverscmake/mathopt_mps.hpp"
#include "mathoptsolverscmake/mathopt_nl.hpp"
#include "mathoptsolverscmake/mathopt_solve.hpp"

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

}

int main(int argc, char* argv[])
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i", po::value<std::string>(), "input model file (required)")
        ("format,f", po::value<FileFormat>(), "input file format: mps, lp, json or nl (guessed from the file extension if not given)")
        ("solver,s", po::value<SolverName>(), "solver: cbc, highs, xpress, cplex, gurobi, knitro, ipopt or scip (required)")
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

    SolveParameters parameters;
    parameters.solver = solver;
    if (has_time_limit)
        parameters.time_limit = time_limit;
    SolveOutput result;
    try {
        result = solve(model, parameters);
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
