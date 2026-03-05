#pragma once
#include <vector>

namespace LibHUPF {

// Returns a reference to the thread-local cache for the last resultant polynomial
// computed by SolveForAngles::calculateRoots().
// Coefficients are in ascending order: coefficient[k] is the coefficient of t^k.
// This is overwritten on every call to calculateRoots(); use primary_poly_ref()
// for the snapshot taken after the primary (non-perturbed) solve.
inline std::vector<double>& poly_capture_ref() {
    thread_local std::vector<double> v;
    return v;
}

// Snapshot of poly_capture_ref() taken immediately after the primary HuPf solve
// (before any perturbation fallback).  Set by IKDH::Solver::solve().
inline std::vector<double>& primary_poly_ref() {
    thread_local std::vector<double> v;
    return v;
}

} // namespace LibHUPF
