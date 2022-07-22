
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "CommonOptions.hpp"
#include "Options.hpp"

#include <iostream>

#include <boost/program_options.hpp>

//#include "PropagationAlgorithm.hpp"

namespace ActsExamples {

  namespace Options {

    /// @brief PropagationAlgorithm options
    ///
    /// @tparam aopt_t Type of the options class from boost
    inline void addPropagationOptions(boost::program_options::options_description& opt) {
      namespace po = boost::program_options;
      using namespace Acts::UnitLiterals;
      opt.add_options()("prop-debug", po::value<bool>()->default_value(false),
                        "Run in debug mode, will create propagation screen output.")(
          "prop-step-collection", po::value<std::string>()->default_value("propagation-steps"),
          "Propgation step collection.")("prop-stepper", po::value<int>()->default_value(1),
                                         "Propgation type: 0 (StraightLine), 1 (Eigen), 2 (Atlas).")(
          "prop-mode", po::value<int>()->default_value(0), "Propgation modes: 0 (inside-out), 1 (surface to surface).")(
          "prop-cov", po::value<bool>()->default_value(false), "Propagate (random) test covariances.")(
          "prop-energyloss", po::value<bool>()->default_value(true),
          "Apply energy loss correction - in extrapolation mode only.")(
          "prop-scattering", po::value<bool>()->default_value(true),
          "Apply scattering correction - in extrapolation mode only.")(
          "prop-record-material", po::value<bool>()->default_value(true),
          "Record the material interaction and - in extrapolation mode only.")(
          "prop-material-collection", po::value<std::string>()->default_value("propagation-material"),
          "Propagation material collection.")("prop-ntests", po::value<size_t>()->default_value(1000),
                                              "Number of tests performed.")(
          "prop-resolve-material", po::value<bool>()->default_value(true), "Resolve all smaterial surfaces.")(
          "prop-resolve-passive", po::value<bool>()->default_value(false), "Resolve all passive surfaces.")(
          "prop-resolve-sensitive", po::value<bool>()->default_value(true), "Resolve all sensitive surfaces.")(
          "prop-d0-sigma", po::value<double>()->default_value(15_um),
          "Sigma of the transverse impact parameter [in mm].")("prop-z0-sigma",
                                                               po::value<double>()->default_value(55_mm),
                                                               "Sigma of the longitudinal impact parameter [in mm].")(
          "prop-phi-sigma", po::value<double>()->default_value(0.001), "Sigma of the azimuthal angle [in rad].")(
          "prop-theta-sigma", po::value<double>()->default_value(0.001), "Sigma of the polar angle [in rad].")(
          "prop-qp-sigma", po::value<double>()->default_value(0.0001 / 1_GeV),
          "Sigma of the signed inverse momentum [in GeV^{-1}].")(
          "prop-t-sigma", po::value<double>()->default_value(1_ns), "Sigma of the time parameter [in ns].")(
          "prop-corr-offd", po::value<Reals<15>>(),
          "The 15 off-diagonal correlation rho(d0,z0), rho(d0,phi), [...], "
          "rho(z0,phi), rho(z0, theta), [...], rho(qop,t). Row-wise.")(
          "prop-phi-range", po::value<Reals<2>>()->default_value({{-M_PI, M_PI}}),
          "Azimutal angle phi range for proprapolated tracks.")("prop-eta-range",
                                                                po::value<Reals<2>>()->default_value({{-4., 4.}}),
                                                                "Pseudorapidity range for proprapolated tracks.")(
          "prop-pt-range", po::value<Reals<2>>()->default_value({{100_MeV, 100_GeV}}),
          "Transverse momentum range for proprapolated tracks [in GeV].")(
          "prop-max-stepsize", po::value<double>()->default_value(3_m),
          "Maximum step size for the propagation [in mm].")(
          "prop-pt-loopers", po::value<double>()->default_value(500_MeV),
          "Transverse momentum below which loops are being detected [in GeV].");
    }

  }  // namespace Options
}  // namespace ActsExamples
