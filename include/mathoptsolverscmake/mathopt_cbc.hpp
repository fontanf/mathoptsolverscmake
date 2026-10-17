#pragma once

#include "mathoptsolverscmake/mathopt.hpp"

#include "coin/OsiCbcSolverInterface.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mathoptsolverscmake
{

struct CbcModel
{
    /**
     * Message handler which writes the log to a file (see set_log_file).
     *
     * Cbc doesn't take ownership of the message handlers passed to it, so it
     * is owned here. It is declared before the model so that it outlives it.
     */
    std::unique_ptr<CoinMessageHandler> log_handler;

    ::CbcModel model;

    CbcModel();
    CbcModel(const CbcModel&) = delete;
    CbcModel& operator=(const CbcModel&) = delete;
};

void load(
        CbcModel& cbc_model,
        const MathOptModel& model);

void reduce_printout(
        CbcModel& cbc_model);

void set_time_limit(
        CbcModel& cbc_model,
        double time_limit);

/**
 * Write the log of Cbc to a file.
 *
 * It must be called after reduce_printout, which would otherwise disable the
 * log.
 */
void set_log_file(
        CbcModel& cbc_model,
        const std::string& log_file);

void solve(
        CbcModel& cbc_model);

double get_solution_value(
        const CbcModel& cbc_model);

std::vector<double> get_solution(
        const CbcModel& cbc_model);

double get_bound(
        const CbcModel& cbc_model);

int get_number_of_nodes(
        const CbcModel& cbc_model);

/*
 * The following functions take the model of Cbc (::CbcModel) directly, so
 * that they can be used in a CbcEventHandler, which only has access to it.
 */

double get_solution_value(
        const ::CbcModel& cbc_model);

std::vector<double> get_solution(
        const ::CbcModel& cbc_model);

double get_bound(
        const ::CbcModel& cbc_model);

int get_number_of_nodes(
        const ::CbcModel& cbc_model);

}
