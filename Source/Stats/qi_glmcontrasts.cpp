/*
 *  qi_glmdiffs.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <Eigen/Dense>

#include "itkTileImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include "MeanImageFilter.h"

#include "Util.h"
#include "Args.h"
#include "ImageIO.h"

class ContrastsFilter : public itk::ImageToImageFilter<QI::VectorVolumeF, QI::VolumeF> {
public:
    /** Standard class typedefs. */
    typedef ContrastsFilter                    Self;
    typedef ImageToImageFilter<QI::VectorVolumeF, QI::VolumeF> Superclass;
    typedef itk::SmartPointer<Self>            Pointer;
    typedef typename QI::VolumeF::RegionType   RegionType;
    
    itkNewMacro(Self);
    itkTypeMacro(Self, Superclass);

    void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        auto input     = this->GetInput(0);
        auto region    = input->GetLargestPossibleRegion();
        auto spacing   = input->GetSpacing();
        auto origin    = input->GetOrigin();
        auto direction = input->GetDirection();
        for (int i = 0; i < m_mat.rows(); i++) {
            auto op = this->GetOutput(i);
            op->SetRegions(region);
            op->SetSpacing(spacing);
            op->SetOrigin(origin);
            op->SetDirection(direction);
            op->Allocate(true);
        }
    }

    void SetMatrix(const Eigen::MatrixXd &d, const Eigen::MatrixXd &c, const bool s = false)
    {
        m_scale = s;
        m_mat = c * ((d.transpose() * d).inverse()) * d.transpose();
        this->SetNumberOfRequiredOutputs(m_mat.rows());
    }

protected:
    bool m_scale = false;
    Eigen::MatrixXd m_mat;

    ContrastsFilter() {
        this->SetNumberOfRequiredInputs(1);
    }
    ~ContrastsFilter() {}

    void DynamicThreadedGenerateData(const RegionType &region) ITK_OVERRIDE {
        const auto input_image = this->GetInput(0);
        itk::ImageRegionConstIterator<QI::VectorVolumeF> input_iter(input_image, region);
        input_iter.GoToBegin();
        std::vector<itk::ImageRegionIterator<QI::VolumeF>> out_iters(m_mat.rows());
        for (int i = 0; i < m_mat.rows(); i++) {
            out_iters.at(i) = itk::ImageRegionIterator<QI::VolumeF>(this->GetOutput(i), region);
            out_iters.at(i).GoToBegin();
        }
        
        while(!input_iter.IsAtEnd()) {
            const auto input_vec = input_iter.Get();

            Eigen::Map<const Eigen::VectorXf> indata(input_vec.GetDataPointer(), input_vec.Size());
            Eigen::VectorXd c = m_mat * indata.cast<double>();
            if (m_scale) {
                c /= indata.mean();
            }
            for (int i = 0; i < m_mat.rows(); i++) {
                out_iters[i].Set(c[i]);
                ++out_iters[i];
            }
            ++input_iter;
        }
    }

private:
    ContrastsFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("A utility for calculating group means, differences etc.\n"
                                "One output file will be generated for each contrast.\n"
                                "\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "IMAGE", "The combined image file from qi_glmsetup");
    args::Positional<std::string> design_path(parser, "DESIGN", "GLM Design matrix from qi_glmsetup");
    args::Positional<std::string> contrasts_path(parser, "CONTRASTS", "Contrasts matrix from qi_glmsetup");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag fraction(parser, "FRACTION", "Output contrasts as fraction of grand mean", {'F',"frac"});
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    QI::ParseArgs(parser, argc, argv, verbose);

    QI_LOG(verbose, "Reading input file " << QI::CheckPos(input_path));
    QI::VectorVolumeF::Pointer merged = QI::ReadVectorImage<float>(QI::CheckPos(input_path));
    QI_LOG(verbose, "Reading design matrix" << QI::CheckPos(design_path));
    Eigen::ArrayXXd design_matrix = QI::ReadArrayFile(QI::CheckPos(design_path));
    QI_LOG(verbose, "Reading contrasts file" << QI::CheckPos(contrasts_path));
    Eigen::ArrayXXd contrasts     = QI::ReadArrayFile(QI::CheckPos(contrasts_path));
    if (design_matrix.rows() != merged->GetNumberOfComponentsPerPixel()) {
        QI_FAIL("Number of rows in design matrix (" << design_matrix.rows() <<
                ") does not match number of volumes in image (" << merged->GetNumberOfComponentsPerPixel() << ")");
    }
    if (design_matrix.cols() != contrasts.cols()) {
        QI_FAIL("Number of columns in design matrix (" << design_matrix.cols() <<
                ") does not match contrasts (" << contrasts.cols() << ")");
    }

    auto con_filter = ContrastsFilter::New();
    con_filter->SetMatrix(design_matrix, contrasts, fraction);
    con_filter->SetInput(0, merged);
    QI_LOG(verbose, "Calculating contrasts" );
    con_filter->Update();
    for (int c = 0; c < contrasts.rows(); c++) {
        QI_LOG(verbose, "Writing contrast " << (c + 1));
        QI::WriteImage(con_filter->GetOutput(c), outarg.Get() + "con" + std::to_string(c + 1) + QI::OutExt());
    }
}

