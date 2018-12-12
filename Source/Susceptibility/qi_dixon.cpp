/*
 *  qi_cestasym.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Core>

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "Util.h"
#include "Args.h"
#include "ImageIO.h"
#include "JSON.h"

int main(int argc, char **argv)
{
    Eigen::initParallel();
    args::ArgumentParser parser("Classic 3-point Dixon\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> input_path(parser, "INPUT", "Input Z-spectrum file");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI_LOG(verbose, "Opening file: " << QI::CheckPos(input_path));
    auto input = QI::ReadVectorImage<std::complex<float>>(QI::CheckPos(input_path), verbose);

    auto phiVol = QI::VolumeF::New();
    phiVol->CopyInformation(input);
    phiVol->SetRegions(input->GetBufferedRegion());
    phiVol->Allocate(true);

    // auto T2sVol = QI::VolumeF::New();
    // T2sVol->CopyInformation(input);
    // T2sVol->SetRegions(input->GetBufferedRegion());
    // T2sVol->Allocate(true);

    auto WfVol = QI::VolumeF::New();
    WfVol->CopyInformation(input);
    WfVol->SetRegions(input->GetBufferedRegion());
    WfVol->Allocate(true);

    auto FfVol = QI::VolumeF::New();
    FfVol->CopyInformation(input);
    FfVol->SetRegions(input->GetBufferedRegion());
    FfVol->Allocate(true);

    auto mt = itk::MultiThreaderBase::New();
    QI_LOG(verbose, "Processing");
    mt->ParallelizeImageRegion<3>(input->GetBufferedRegion(),
        [&](const QI::VectorVolumeF::RegionType &region) {
            itk::ImageRegionConstIterator<QI::VectorVolumeXF> in_it(input, region);
            itk::ImageRegionIterator<QI::VolumeF> phi_it(phiVol, region);
            itk::ImageRegionIterator<QI::VolumeF> Wf_it(WfVol, region);
            itk::ImageRegionIterator<QI::VolumeF> Ff_it(FfVol, region);
            // itk::ImageRegionIterator<QI::VolumeF>  T2s(T2sVol, region);
            for (in_it.GoToBegin(); !in_it.IsAtEnd(); ++in_it, ++phi_it, ++Wf_it, ++Ff_it)
            {
                const Eigen::Map<const Eigen::ArrayXcf> data(in_it.Get().GetDataPointer(), 3);
                const auto phi = std::arg(data[2] / data[0]) / 2.;
                const auto psi = cos(std::arg(data[1] / data[0]) - phi);
                Wf_it.Set(0.5 * (0.5 * (std::abs(data[0]) + std::abs(data[2])) + psi * std::abs(data[1])));
                Ff_it.Set(0.5 * (0.5 * (std::abs(data[0]) + std::abs(data[2])) - psi * std::abs(data[1])));
                phi_it.Set(phi);
            }
        },
        nullptr);
    QI::WriteImage(phiVol, outarg.Get() + "DIX_Phi" + QI::OutExt(), verbose);
    QI::WriteImage(WfVol, outarg.Get() + "DIX_Fw" + QI::OutExt(), verbose);
    QI::WriteImage(FfVol, outarg.Get() + "DIX_Ff" + QI::OutExt(), verbose);
    QI_LOG(verbose, "Finished.");
    return EXIT_SUCCESS;
}
