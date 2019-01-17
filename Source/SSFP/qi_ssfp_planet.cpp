/*
 *  qdespot1.cpp - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <array>
#include <iostream>

#include <Eigen/Core>

#include "Args.h"
#include "ImageIO.h"
#include "ModelFitFilter.h"
#include "SSFPSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct PLANETModel {
    using DataType      = double;
    using ParameterType = double;
    using SequenceType  = QI::SSFPSequence;

    static constexpr int NV = 3;
    static constexpr int ND = 0;
    static constexpr int NF = 1;

    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const QI_ARRAYN(double, NF) fixed_defaults;
    const SequenceType &sequence;

    size_t num_outputs() const { return 3; }
    int    output_size(int /* Unused */) { return sequence.size(); }

    auto signals(const Eigen::ArrayBase<QI_ARRAYN(double, NV)> &v,
                 const QI_ARRAYN(double, NF) & f) const -> std::vector<QI_ARRAY(double)> {
        const double &                PD      = v[0];
        const double &                T1      = v[1];
        const double &                T2      = v[2];
        const double &                B1      = f[0];
        const double                  E1      = exp(-sequence.TR / T1);
        const double                  E2      = exp(-sequence.TR / T2);
        const Eigen::ArrayXd          alpha   = B1 * sequence.FA;
        const Eigen::ArrayXd          d       = (1. - E1 * E2 * E2 - (E1 - E2 * E2) * cos(alpha));
        const Eigen::ArrayXd          G       = PD * sqrt(E2) * (1 - E1) * sin(alpha) / d;
        const Eigen::ArrayXd          a       = Eigen::ArrayXd::Constant(sequence.size(), E2);
        const Eigen::ArrayXd          b       = E2 * (1. - E1) * (1. + cos(alpha)) / d;
        std::vector<QI_ARRAY(double)> outputs = {G, a, b};
        return outputs;
    }
};
std::array<const std::string, 3> PLANETModel::varying_names{{"PD"s, "T1"s, "T2"s}};
std::array<const std::string, 1> PLANETModel::fixed_names{{"B1"s}};
const QI_ARRAYN(double, 1) PLANETModel::fixed_defaults{1.0};

struct PLANETFit {
    static const bool Blocked = true;
    static const bool Indexed = false;
    using InputType           = double;
    using OutputType          = double;
    using ResidualType        = double;
    using FlagType            = int;
    using ModelType           = PLANETModel;
    ModelType model;

    int n_inputs() const { return 3; }
    int input_size(const int /* Unused */) const { return 1; }
    int n_fixed() const { return 1; }
    int n_outputs() const { return 3; }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs, const Eigen::ArrayXd &fixed,
                          QI_ARRAYN(OutputType, PLANETModel::NV) & out, ResidualType & /* Unused */,
                          std::vector<Eigen::ArrayXd> & /* Unused */, FlagType & /* Unused */,
                          const int block) const {
        const double &G    = inputs[0][0];
        const double &a    = inputs[1][0];
        const double &b    = inputs[2][0];
        const double  b1   = fixed[0];
        const double  cosa = cos(b1 * model.sequence.FA(block));
        const double  sina = sin(b1 * model.sequence.FA(block));
        const double  T1   = -model.sequence.TR / log((a * (1. + cosa - a * b * cosa) - b) /
                                                   (a * (1. + cosa - a * b) - b * cosa));
        const double  T2   = -model.sequence.TR / log(a);
        const double  E1   = exp(-model.sequence.TR / T1);
        const double  E2   = a; // For simplicity copying formulas
        const double  PD =
            G * (1. - E1 * cosa - E2 * E2 * (E1 - cosa)) / (sqrt(E2) * (1. - E1) * sina);
        out[0] = PD;
        out[1] = T1;
        out[2] = T2;
        return std::make_tuple(true, "");
    }
};

int main(int argc, char **argv) {
    args::ArgumentParser parser(
        "Calculates T1&T2 from SSFP Ellipse Parameters.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> G_filename(parser, "G", "Ellipse parameter G");
    args::Positional<std::string> a_filename(parser, "a", "Ellipse parameter a");
    args::Positional<std::string> b_filename(parser, "b", "Ellipse parameter b");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> out_prefix(parser, "OUTPREFIX", "Add a prefix to output filenames",
                                            {'o', "out"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::ValueFlag<std::string> subregion(
        parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::ValueFlag<std::string> seq_arg(parser, "FILE",
                                         "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)",
        {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(G_filename);
    QI::CheckPos(a_filename);
    QI::CheckPos(b_filename);

    QI_LOG(verbose, "Reading sequence information");
    rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::SSFPSequence    ssfp(QI::GetMember(input, "SSFP"));
    PLANETModel         model{ssfp};
    if (simulate) {
        QI::SimulateModel<PLANETModel, true>(input, model, {B1.Get()},
                                             {G_filename.Get(), a_filename.Get(), b_filename.Get()},
                                             verbose, simulate.Get());
    } else {
        PLANETFit fit{model};
        auto      fit_filter = itk::ModelFitFilter<PLANETFit>::New(&fit);
        fit_filter->SetVerbose(verbose);
        fit_filter->SetInput(0, QI::ReadVectorImage(G_filename.Get(), verbose));
        fit_filter->SetInput(1, QI::ReadVectorImage(a_filename.Get(), verbose));
        fit_filter->SetInput(2, QI::ReadVectorImage(b_filename.Get(), verbose));
        fit_filter->SetBlocks(ssfp.size());
        if (B1)
            fit_filter->SetFixed(0, QI::ReadImage(B1.Get(), verbose));
        if (mask)
            fit_filter->SetMask(QI::ReadImage(mask.Get(), verbose));
        if (subregion)
            fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
        QI_LOG(verbose, "Processing");
        if (verbose) {
            auto monitor = QI::GenericMonitor::New();
            fit_filter->AddObserver(itk::ProgressEvent(), monitor);
        }
        fit_filter->Update();
        QI_LOG(verbose, "Elapsed time was " << fit_filter->GetTotalTime() << "s\n"
                                            << "Writing results files.");
        std::string outPrefix = out_prefix.Get() + "PLANET_";
        for (int i = 0; i < PLANETModel::NV; i++) {
            QI::WriteVectorImage(fit_filter->GetOutput(i),
                                 outPrefix + PLANETModel::varying_names.at(i) + QI::OutExt());
        }
        QI_LOG(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
