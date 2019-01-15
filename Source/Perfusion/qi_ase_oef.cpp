/*
 *  qi_ase_oef.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Dense>
#include <iostream>

// #define QI_DEBUG_BUILD

#include "Args.h"
#include "FitFunction.h"
#include "ImageIO.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "MultiEchoSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct ASEModel {
    using SequenceType  = QI::MultiEchoFlexSequence;
    using DataType      = double;
    using ParameterType = double;

    static const int NV = 4;
    static const int NF = 3;
    using VaryingArray  = QI_ARRAYN(ParameterType, NV);
    using FixedArray    = QI_ARRAYN(ParameterType, NF);
    const std::array<const std::string, NV> varying_names{{"S0"s, "R2p"s, "DBV"s, "dT"s}};
    const std::array<const std::string, NF> fixed_names{{}};
    const FixedArray                        fixed_defaults{};

    const SequenceType &sequence;
    const double        B0, Tc;
    VaryingArray        start, bounds_lo, bounds_hi;

    static constexpr double kappa    = 0.03;         // Conversion factor
    static constexpr double gamma    = 42.577e6;     // Gyromagnetic Ratio
    static constexpr double delta_X0 = 0.264e-6;     // Susc diff oxy and fully de-oxy blood
    static constexpr double Hb       = 0.34 / kappa; // Hct = 0.34;
    ASEModel(const SequenceType &s, const double B0in)
        : sequence{s}, B0{B0in}, Tc{0.015 / (B0in / 3)}
    // Nic Blockley uses Tc = 15 ms for 3T, scale for other field-strengths
    {
        start << 1., 1., 0.1, 0.;
        bounds_lo << 0.1, -20., 1.e-6, -5.e-3;
        bounds_hi << 10., 20., 0.99, 5.e-3;
    }

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &varying, const FixedArray & /* Unused */) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T      = typename Derived::Scalar;
        const T &S0  = varying[0];
        const T &R2p = varying[1];
        const T &DBV = varying[2];
        // const T &OEF = varying[3];
        const T &dT = varying[3];

        // Calculate characteristic frequency
        // const auto dHb = 3. * R2p / (DBV * 4. * gamma * M_PI * delta_X0 * kappa * B0);
        const auto dw = R2p / DBV;
        const auto Tc = 1. / dw;
        // Adjust TE to find peak
        const auto aTE = (sequence.TE + dT).abs();
        QI_ARRAY(T) S(sequence.size());
        for (int i = 0; i < sequence.size(); i++) {
            const auto tau    = aTE(i);
            const auto S_quad = S0 * exp(-(2. / 9.) * (R2p * tau) * (R2p * tau) / DBV);
            const auto S_lin  = S0 * exp(-tau * R2p + DBV);
            std::cout << "< 1.5Tc " << (tau < (1.5 * Tc)) << "\ttau " << tau << "\tSl " << S_lin
                      << "\tSq " << S_quad << std::endl;
            if (tau < (1.5 * Tc)) {
                S(i) = S_quad;
            } else {
                S(i) = S_lin;
            }
        }
        QI_DB(Tc)

        QI_DB(S0)
        QI_DB(R2p)
        QI_DB(DBV)
        QI_DB(dT)
        QI_DBVEC(sequence.TE)
        QI_DBVEC(aTE)
        QI_DBVEC(S)
        return S;
    }
};

using ASEFit = QI::ScaledNLLSFitFunction<ASEModel>;

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "Calculates the OEF from ASE data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "ASE_FILE", "Input ASE file");
    args::HelpFlag                help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename",
                                        {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::ValueFlag<double> B0(parser, "B0", "Field-strength (Tesla), default 3", {'B', "B0"}, 3.0);
    args::ValueFlag<std::string> gradx(
        parser, "GRADX", "Gradient of field-map in x-direction for MFG correction", {'x', "gradx"});
    args::ValueFlag<std::string> grady(
        parser, "GRADY", "Gradient of field-map in y-direction for MFG correction", {'y', "grady"});
    args::ValueFlag<std::string> gradz(
        parser, "GRADZ", "Gradient of field-map in z-direction for MFG correction", {'z', "gradz"});
    args::ValueFlag<double> slice_arg(
        parser, "SLICE THICKNESS",
        "Slice-thickness for MFG calculation (useful if there was a slice gap)", {'s', "slice"});
    args::ValueFlag<std::string> subregion(
        parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::ValueFlag<std::string> json_file(parser, "FILE",
                                           "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)",
        {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    rapidjson::Document json = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    QI::MultiEchoFlexSequence sequence(json["MultiEchoFlex"]);
    ASEModel                  model{sequence, B0.Get()};
    if (simulate) {
        QI::SimulateModel<ASEModel, false>(json, model, {gradx.Get(), grady.Get(), gradz.Get()},
                                           {input_path.Get()}, verbose, simulate.Get());
    } else {
        ASEFit fit(model);
        auto   fit_filter = itk::ModelFitFilter<ASEFit>::New(&fit);
        fit_filter->SetVerbose(verbose);
        fit_filter->SetOutputAllResiduals(false);
        QI_LOG(verbose, "Reading ASE data from: " << QI::CheckPos(input_path));
        auto input = QI::ReadVectorImage(QI::CheckPos(input_path));
        // QI::VolumeF::SpacingType  vox_size = input->GetSpacing();
        // if (slice_arg) {
        //     vox_size[2] = slice_arg.Get();
        // }
        fit_filter->SetInput(0, input);
        if (mask)
            fit_filter->SetMask(QI::ReadImage(mask.Get()));
        if (gradx)
            fit_filter->SetFixed(0, QI::ReadImage(gradx.Get()));
        if (grady)
            fit_filter->SetFixed(1, QI::ReadImage(grady.Get()));
        if (gradz)
            fit_filter->SetFixed(2, QI::ReadImage(gradz.Get()));
        if (subregion) {
            fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
        }
        QI_LOG(verbose, "Processing");
        if (verbose) {
            auto monitor = QI::GenericMonitor::New();
            fit_filter->AddObserver(itk::ProgressEvent(), monitor);
        }
        fit_filter->Update();
        QI_LOG(verbose, "Elapsed time was " << fit_filter->GetTotalTime() << "s");

        const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
        for (size_t i = 0; i < model.NV; i++) {
            const std::string fname = outPrefix + "_" + model.varying_names[i] + QI::OutExt();
            QI_LOG(verbose, "Writing file: " << fname);
            QI::WriteImage(fit_filter->GetOutput(i), fname);
        }
    }
    return EXIT_SUCCESS;
}
