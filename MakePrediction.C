
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <utility>

#include <algorithm>
#include <functional>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TChain.h"
#include "TBranch.h"
#include "TSpline.h"
#include "TVirtualFFT.h"
#include "Math/GSLMinimizer.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

ClassImp(TSplinePoly);
ClassImp(TSplinePoly3);
ClassImp(TSplinePoly5);
ClassImp(TSpline3);
ClassImp(TSpline5);
ClassImp(TSpline);

#include <math.h>
#include <complex>
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
#include <complex>
#include <cmath>
using namespace Eigen;

#include "PrepareDarkData.C"

MatrixXcd mtx_data;
MatrixXcd mtx_dark;
MatrixXcd mtx_data_bkgd;
MatrixXcd mtx_gamma_raw;
MatrixXcd mtx_gamma;
MatrixXcd mtx_eigenvalue_gamma;
MatrixXcd mtx_eigenvalue_data;
MatrixXcd mtx_eigenvalue_dark;
MatrixXcd mtx_eigenvalue_init;
MatrixXcd mtx_eigenvalue_vari;
MatrixXcd mtx_eigenvector_gamma;
MatrixXcd mtx_eigenvector_data;
MatrixXcd mtx_eigenvector_dark;
MatrixXcd mtx_eigenvector_init;
MatrixXcd mtx_eigenvector_vari;
MatrixXcd mtx_eigenvector_inv_gamma;
MatrixXcd mtx_eigenvector_inv_data;
MatrixXcd mtx_eigenvector_inv_dark;
MatrixXcd mtx_eigenvector_inv_init;
MatrixXcd mtx_eigenvector_inv_vari;
MatrixXcd mtx_eigenvector;
MatrixXcd mtx_eigenvector_inv;
MatrixXcd mtx_eigenvalue;
MatrixXcd mtx_eigenval_data_redu;
MatrixXcd mtx_data_redu;
ComplexEigenSolver<MatrixXcd> eigensolver_gamma;
ComplexEigenSolver<MatrixXcd> eigensolver_gamma_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_bkgd;
ComplexEigenSolver<MatrixXcd> eigensolver_bkgd_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_dark;
ComplexEigenSolver<MatrixXcd> eigensolver_dark_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_data;
ComplexEigenSolver<MatrixXcd> eigensolver_data_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_init;
ComplexEigenSolver<MatrixXcd> eigensolver_init_transpose;
int binx_blind_global;
int biny_blind_global;
int VaryLeftOrRightVector;
int VaryNthVector;

double GammaScale = 0.;
double MinChi2Unblind = 1e10;

void ResetMatrixDimension()
{
    mtx_data = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_dark = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_data_bkgd = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_gamma_raw = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_gamma = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvalue_gamma = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvalue_data = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvalue_dark = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvalue_init = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvalue_vari = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_gamma = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_data = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_dark = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_init = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_vari = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_inv_gamma = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_inv_data = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_inv_dark = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_inv_init = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_inv_vari = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvector_inv = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvalue = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenval_data_redu = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_data_redu = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
}
MatrixXcd fillMatrix(TH2D* hist)
{
    MatrixXcd matrix(hist->GetNbinsX(),hist->GetNbinsY());
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            matrix(binx,biny) = hist->GetBinContent(binx+1,biny+1);
        }
    }
    return matrix;
}
void fill2DHistogramAbs(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            hist->SetBinContent(binx+1,biny+1,max(0.,mtx(binx,biny).real()));
        }
    }
}
void fill2DHistogram(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            hist->SetBinContent(binx+1,biny+1,mtx(binx,biny).real());
        }
    }
}
void fill1DHistogram(TH1D* hist,MatrixXcd mtx, int rank, double weight)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        double old_content = hist->GetBinContent(binx+1);
        double new_content = mtx(binx,mtx.cols()-1-rank).real()*weight;
        hist->SetBinContent(binx+1,old_content+new_content);
    }
    double integral = hist->Integral();
    hist->Scale(abs(integral)/integral);
}
MatrixXcd MakeRealEigenvectors(MatrixXcd mtx_input)
{
    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    for (int col=0;col<mtx_input.cols();col++)
    {
        if (col<mtx_input.cols()-NumberOfRealEigenvectors) continue;
        for (int row=0;row<mtx_input.rows();row++)
        {
            double real = mtx_input(row,col).real();
            mtx_output(row,col) = real/abs(real)*std::abs(mtx_input(row,col));
        }
    }
    return mtx_output;
}
MatrixXcd GetLambdaNuMatrix_v2(MatrixXcd mtx_input_lambda, MatrixXcd mtx_input_eta)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_input_lambda;
    eigensolver_input_lambda = ComplexEigenSolver<MatrixXcd>(mtx_input_lambda);
    ComplexEigenSolver<MatrixXcd> eigensolver_input_eta;
    eigensolver_input_eta = ComplexEigenSolver<MatrixXcd>(mtx_input_eta);
    ComplexEigenSolver<MatrixXcd> eigensolver_input_eta_transpose;
    eigensolver_input_eta_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input_eta.transpose());

    MatrixXcd mtx_U_r = eigensolver_input_eta.eigenvectors();
    MatrixXcd mtx_U_l = eigensolver_input_eta_transpose.eigenvectors();
    mtx_U_r = MakeRealEigenvectors(mtx_U_r);
    mtx_U_l = MakeRealEigenvectors(mtx_U_l);
    MatrixXcd mtx_H = mtx_U_l.transpose()*mtx_U_r;
    MatrixXcd mtx_lambdanu(mtx_input_lambda.rows(),mtx_input_lambda.cols());
    for (int row=0;row<mtx_input_lambda.rows();row++)
    {
        for (int col=0;col<mtx_input_lambda.cols();col++)
        {
            mtx_lambdanu(row,col) = 0.;
            if (row==col)
            {
                if (col>=mtx_input_lambda.cols()-NumberOfRealEigenvectors)
                {
                    mtx_lambdanu(row,col) = eigensolver_input_lambda.eigenvalues()(col).real()/mtx_H(row,col).real();
                }
                else
                {
                    mtx_lambdanu(row,col) = eigensolver_input_lambda.eigenvalues()(col)/mtx_H(row,col);
                }
            }
        }
    }
    return mtx_lambdanu;
}
MatrixXcd CutoffEigenvalueMatrix(MatrixXcd mtx_input, int cutoff_entry)
{
    MatrixXcd mtx_output = mtx_input;
    for (int entry=1;entry<=mtx_input.rows();entry++)
    {
        if (entry>cutoff_entry)
        {
            mtx_output(mtx_input.rows()-entry,mtx_input.rows()-entry) = 0.;
        }
    }
    return mtx_output;
}

MatrixXd SmoothingRealVectorsFFT(MatrixXd mtx_input)
{
    MatrixXd mtx_output = mtx_input;
    double *in = new double[2*((mtx_input.rows()+1)/2+1)];
    double *re_full = new double[mtx_input.rows()];
    double *im_full = new double[mtx_input.rows()];
    double *out = new double[mtx_input.rows()];
    for (int col=0;col<mtx_input.cols();col++)
    {
        if (col<mtx_input.cols()-NumberOfRealEigenvectors) continue;
        for (int row=0;row<mtx_input.rows();row++)
        {
            in[row] = mtx_input(row,col);
        }
        int n_size = mtx_input.rows()+1;
        TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &n_size, "R2C ES K");
        fft_own->SetPoints(in);
        fft_own->Transform();
        fft_own->GetPointsComplex(re_full,im_full); //Copy all the output points
        // filtering high frequency
        for (int row=0;row<mtx_input.rows();row++)
        {
            if (row>=6)
            {
                re_full[row] = 0.;
                im_full[row] = 0.;
            }
        }
        int n = mtx_input.rows();
        TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
        fft_back->SetPointsComplex(re_full,im_full);
        fft_back->Transform();
        fft_back->GetPoints(out);
        for (int row=0;row<mtx_input.rows();row++)
        {
            // the y-axes has to be rescaled (factor 1/bins)
            mtx_output(row,col) = out[row]/double(mtx_input.rows());
        }
    }
    return mtx_output;
}

MatrixXd SmoothingRealVectors(MatrixXd mtx_input)
{
    MatrixXd mtx_output = mtx_input;
    //for (int col=0;col<mtx_input.cols()-1;col++)
    for (int col=0;col<mtx_input.cols();col++)
    {
        for (int row=1;row<mtx_input.rows()-1;row++)
        {
            mtx_output(row,col) = (0.5*mtx_input(row+1,col)+mtx_input(row,col)+0.5*mtx_input(row-1,col))/2.;
        }
    }
    return mtx_output;
}

MatrixXd SmoothingRealVectorsSpline(MatrixXd mtx_input)
{
    MatrixXd mtx_output = mtx_input;
    for (int col=0;col<mtx_input.cols();col++)
    {
        TH1D Hist_temp = TH1D("Hist_temp","",mtx_input.rows(),0,mtx_input.rows());
        for (int row=0;row<mtx_input.rows();row++)
        {
            Hist_temp.SetBinContent(row+1,mtx_output(row,col));
        }
        int n_rebin = 2;
        Hist_temp.Rebin(n_rebin);
        TSpline3 spline_eigenvec_real(&Hist_temp);
        for (int row=0;row<mtx_input.rows();row++)
        {
            double xx = Hist_temp.GetBinCenter(row+1);
            mtx_output(row,col) = spline_eigenvec_real.Eval(xx)/double(n_rebin); 
        }
    }
    return mtx_output;
}

void SetInitialSpectralvectors(int binx_blind, int biny_blind, MatrixXcd mtx_input)
{

    std::cout << "SetInitialSpectralvectors..." << std::endl;

    eigensolver_init = ComplexEigenSolver<MatrixXcd>(mtx_input);
    eigensolver_init_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_U_r_init = eigensolver_init.eigenvectors();
    MatrixXcd mtx_U_l_init = eigensolver_init_transpose.eigenvectors();
    mtx_U_r_init = MakeRealEigenvectors(mtx_U_r_init);
    mtx_U_l_init = MakeRealEigenvectors(mtx_U_l_init);

    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    eigensolver_dark_transpose = ComplexEigenSolver<MatrixXcd>(mtx_dark.transpose());
    MatrixXcd mtx_U_r_dark = eigensolver_dark.eigenvectors();
    MatrixXcd mtx_U_l_dark = eigensolver_dark_transpose.eigenvectors();
    mtx_U_r_dark = MakeRealEigenvectors(mtx_U_r_dark);
    mtx_U_l_dark = MakeRealEigenvectors(mtx_U_l_dark);

    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
    eigensolver_data_transpose = ComplexEigenSolver<MatrixXcd>(mtx_data.transpose());
    MatrixXcd mtx_U_r_data = eigensolver_data.eigenvectors();
    MatrixXcd mtx_U_l_data = eigensolver_data_transpose.eigenvectors();
    mtx_U_r_data = MakeRealEigenvectors(mtx_U_r_data);
    mtx_U_l_data = MakeRealEigenvectors(mtx_U_l_data);

    eigensolver_gamma = ComplexEigenSolver<MatrixXcd>(mtx_gamma_raw);
    eigensolver_gamma_transpose = ComplexEigenSolver<MatrixXcd>(mtx_gamma_raw.transpose());
    MatrixXcd mtx_U_r_gamma = eigensolver_gamma.eigenvectors();
    MatrixXcd mtx_U_l_gamma = eigensolver_gamma_transpose.eigenvectors();
    mtx_U_r_gamma = MakeRealEigenvectors(mtx_U_r_gamma);
    mtx_U_l_gamma = MakeRealEigenvectors(mtx_U_l_gamma);

    MatrixXcd mtx_lambdanu = GetLambdaNuMatrix_v2(mtx_input,mtx_input);
    mtx_lambdanu = CutoffEigenvalueMatrix(mtx_lambdanu, NumberOfRealEigenvectors);

    MatrixXcd mtx_lambdanu_dark = GetLambdaNuMatrix_v2(mtx_dark,mtx_dark);
    mtx_lambdanu_dark = CutoffEigenvalueMatrix(mtx_lambdanu_dark, NumberOfRealEigenvectors);

    MatrixXcd mtx_lambdanu_data = GetLambdaNuMatrix_v2(mtx_data,mtx_data);
    mtx_lambdanu_data = CutoffEigenvalueMatrix(mtx_lambdanu_data, NumberOfRealEigenvectors);

    MatrixXcd mtx_lambdanu_gamma = GetLambdaNuMatrix_v2(mtx_gamma_raw,mtx_gamma_raw);
    mtx_lambdanu_gamma = CutoffEigenvalueMatrix(mtx_lambdanu_gamma, NumberOfRealEigenvectors);

    const std::complex<double> If(0.0, 1.0);
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;

    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector_init(row,col) = mtx_U_r_init(row,col);
            mtx_eigenvector_inv_init(row,col) = mtx_U_l_init.transpose()(row,col);
            mtx_eigenvector_vari(row,col) = 0.;
            mtx_eigenvector_inv_vari(row,col) = 0.;
            mtx_eigenvalue_init(row,col) = mtx_lambdanu(row,col);
            mtx_eigenvalue_vari(row,col) = 0.;
            mtx_eigenvector_dark(row,col) = mtx_U_r_dark(row,col);
            mtx_eigenvector_inv_dark(row,col) = mtx_U_l_dark.transpose()(row,col);
            mtx_eigenvalue_dark(row,col) = mtx_lambdanu_dark(row,col);
            mtx_eigenvector_data(row,col) = mtx_U_r_data(row,col);
            mtx_eigenvector_inv_data(row,col) = mtx_U_l_data.transpose()(row,col);
            mtx_eigenvalue_data(row,col) = mtx_lambdanu_data(row,col);
            mtx_eigenvector_gamma(row,col) = mtx_U_r_gamma(row,col);
            mtx_eigenvector_inv_gamma(row,col) = mtx_U_l_gamma.transpose()(row,col);
            mtx_eigenvalue_gamma(row,col) = mtx_lambdanu_gamma(row,col);
        }
    }

    std::cout << "SetInitialSpectralvectors...Done" << std::endl;

}
MatrixXcd GetErrorMap()
{
    MatrixXcd mtx_error(N_bins_for_deconv,N_bins_for_deconv);
    mtx_error = mtx_dark-mtx_eigenvector_dark*mtx_eigenvalue_dark*mtx_eigenvector_inv_dark;
    return mtx_error;
}
void RingBackgroundSkymap(TH2D* hist_input, TH2D* hist_output, TH2D* hist_response)
{
    double radius_inner = 0.2;
    double radius_outer = 0.4;
    for (int bx_target=1;bx_target<=hist_input->GetNbinsX();bx_target++)
    {
        for (int by_target=1;by_target<=hist_input->GetNbinsY();by_target++)
        {
            double target_x = hist_input->GetXaxis()->GetBinCenter(bx_target);
            double target_y = hist_input->GetYaxis()->GetBinCenter(by_target);
            double target_norm = hist_response->GetBinContent(bx_target,by_target);
            double ring_count = 0.;
            double ring_norm = 0.;
            for (int bx_ring=1;bx_ring<=hist_input->GetNbinsX();bx_ring++)
            {
                for (int by_ring=1;by_ring<=hist_input->GetNbinsY();by_ring++)
                {
                    double ring_x = hist_input->GetXaxis()->GetBinCenter(bx_ring);
                    double ring_y = hist_input->GetYaxis()->GetBinCenter(by_ring);
                    if (abs(target_x-ring_x)>radius_outer) continue;
                    if (abs(target_y-ring_y)>radius_outer) continue;
                    double distance = (target_x-ring_x)*(target_x-ring_x) + (target_y-ring_y)*(target_y-ring_y);
                    if (distance<radius_inner*radius_inner) continue;
                    if (distance>radius_outer*radius_outer) continue;
                    ring_norm += hist_response->GetBinContent(bx_ring,by_ring);
                    ring_count += hist_input->GetBinContent(bx_ring,by_ring);
                }
            }
            if (ring_norm==0.)
            {
                hist_output->SetBinContent(bx_target,by_target,0.);
            }
            else 
            {
                hist_output->SetBinContent(bx_target,by_target,ring_count*target_norm/ring_norm);
            }
            double data_count = hist_input->GetBinContent(bx_target,by_target);
            double bkgd_count = hist_output->GetBinContent(bx_target,by_target);
            //if (data_count<100) continue;
            //std::cout << "data_count = " << data_count << ", bkgd_count = " << bkgd_count << std::endl;
        }
    }
}
double BlindedChi2(TH2D* hist_data, TH2D* hist_dark, TH2D* hist_model, TH2D* hist_error, int LeftOrRight)
{
    int binx_blind_upper = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind_upper = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    int binx_blind_lower = hist_data->GetXaxis()->FindBin(-MSCL_cut_blind);
    int biny_blind_lower = hist_data->GetYaxis()->FindBin(-MSCW_cut_blind);
    int binx_upper = hist_data->GetNbinsX();
    int biny_upper = hist_data->GetNbinsY();
    double chi2 = 0.;
    double nbins = 0.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            if (bx>binx_upper || by>biny_upper) continue;
            if (LeftOrRight==2)
            {
                if (by<biny_blind_upper) continue;
            }
            if (LeftOrRight==1)
            {
                if (bx<binx_blind_upper) continue;
            }
            double data = hist_data->GetBinContent(bx,by);
            double dark = hist_dark->GetBinContent(bx,by);
            double dark_err = abs(hist_error->GetBinContent(bx,by));
            double model = hist_model->GetBinContent(bx,by);
            double weight = 1.;
            double dx = hist_data->GetXaxis()->GetBinCenter(bx)-(1.);
            double dy = hist_data->GetYaxis()->GetBinCenter(by)-(1.);
            double width = 0.1;
            double data_err = max(1.0,pow(data,0.5));
            double model_err = max(1.0,pow(abs(model),0.5));
            //weight = 1./(data_err*data_err+model_err*model_err);
            //weight = 1./(data_err*data_err);
            double chi2_this = weight*pow(data-model,2);
            if (isnan(chi2_this))
            {
                if (isnan(weight)) std::cout << "weight==nan!!!" << std::endl;
                if (isnan(pow(data-model,2))) std::cout << "pow(data-model,2)==nan!!!" << std::endl;
                std::cout << "model = " << model << std::endl;
                std::cout << "data = " << data << std::endl;
                return 1e10;
            }
            if (bx<binx_blind_upper && by<biny_blind_upper)
            //if (bx<binx_blind_upper && by<biny_blind_upper && bx>=binx_blind_lower && by>=biny_blind_lower)
            {
                //if ((data-model)<0.) chi2 += chi2_this;
                continue;
            }
            else
            {
                chi2 += chi2_this;
            }
            nbins += 1.;
        }
    }
    
    return chi2;
}
double GetChi2Function(MatrixXcd mtx_model, bool isBlind)
{

    MatrixXcd mtx_error = GetErrorMap();

    TH2D hist_data = TH2D("hist_data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_dark = TH2D("hist_dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_model = TH2D("hist_model","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_error = TH2D("hist_error","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);

    fill2DHistogram(&hist_data,mtx_data);
    fill2DHistogram(&hist_dark,mtx_dark);
    fill2DHistogram(&hist_model,mtx_model);
    fill2DHistogram(&hist_error,mtx_error);

    double chi2 = 0.;
    int binx_blind_upper = hist_data.GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind_upper = hist_data.GetYaxis()->FindBin(MSCW_cut_blind);
    int binx_blind_lower = hist_data.GetXaxis()->FindBin(-MSCL_cut_blind);
    int biny_blind_lower = hist_data.GetYaxis()->FindBin(-MSCW_cut_blind);
    int binx_upper = hist_data.GetNbinsX();
    int biny_upper = hist_data.GetNbinsY();
    for (int bx=1;bx<=hist_data.GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data.GetNbinsY();by++)
        {
            if (bx>binx_upper || by>biny_upper) continue;
            double data = hist_data.GetBinContent(bx,by);
            double model = hist_model.GetBinContent(bx,by);
            double weight = 1.;
            double chi2_this = weight*pow(data-model,2);
            if (isBlind && bx<binx_blind_upper && by<biny_blind_upper)
            {
                continue;
            }
            else
            {
                chi2 += chi2_this;
            }
        }
    }

    return chi2;

}
pair<MatrixXd,VectorXd> RemoveEmptyRows(MatrixXd mtx_input,VectorXd vtr_input)
{
    int n_empty_rows = 0;
    vector<bool> empty_rows;
    for (int row=0;row<mtx_input.rows();row++)
    {
        bool row_is_empty = true;
        for (int col=0;col<mtx_input.cols();col++)
        {
            if (mtx_input(row,col)!=0.)
            {
                row_is_empty = false;
            }
        }
        empty_rows.push_back(row_is_empty);
        if (row_is_empty)
        {
            n_empty_rows += 1;
        }
    }
    std::cout << "n_empty_rows = " << n_empty_rows << std::endl;
    MatrixXd mtx_output = MatrixXd::Zero(mtx_input.rows()-n_empty_rows,mtx_input.cols());
    VectorXd vtr_output = VectorXd::Zero(mtx_input.rows()-n_empty_rows);
    int new_row = -1;
    for (int row=0;row<mtx_input.rows();row++)
    {
        if (empty_rows.at(row)) continue;
        new_row += 1;
        vtr_output(new_row) = vtr_input(row);
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_output(new_row,col) = mtx_input(row,col);
        }
    }
    return std::make_pair(mtx_output,vtr_output);
}
VectorXd SolutionWithConstraints(MatrixXd mtx_big, MatrixXd mtx_constraints_input, VectorXd vtr_delta, VectorXd vtr_constraints_delta_input)
{

    MatrixXd mtx_constraints = RemoveEmptyRows(mtx_constraints_input,vtr_constraints_delta_input).first;
    VectorXd vtr_constraints_delta = RemoveEmptyRows(mtx_constraints_input,vtr_constraints_delta_input).second;
    std::cout << "mtx_constraints.rows() = " << mtx_constraints.rows() << std::endl;
    if (mtx_constraints.rows()==0)
    {
        BDCSVD<MatrixXd> svd(mtx_big, ComputeThinU | ComputeThinV);
        svd.setThreshold(svd_threshold);  // size of singular value to be considered as nonzero.
        VectorXd vtr_vari_big = VectorXd::Zero(mtx_big.cols());
        vtr_vari_big = svd.solve(vtr_delta);
        return vtr_vari_big;
    }

    //MatrixXd mtx_Bigger = MatrixXd::Zero(mtx_big.rows()+mtx_constraints.rows(),mtx_big.cols());
    //mtx_Bigger.block(0,0,mtx_big.rows(),mtx_big.cols()) = mtx_big;
    //mtx_Bigger.block(mtx_big.rows(),0,mtx_constraints.rows(),mtx_big.cols()) = mtx_constraints;
    //VectorXd vtr_bigger_delta = VectorXd::Zero(mtx_big.rows()+mtx_constraints.rows());
    //vtr_bigger_delta.segment(0,mtx_big.cols()) = vtr_delta;
    //vtr_bigger_delta.segment(vtr_delta.cols(),vtr_constraints_delta.size()) = vtr_constraints_delta;
    //BDCSVD<MatrixXd> svd(mtx_Bigger, ComputeThinU | ComputeThinV);
    //svd.setThreshold(svd_threshold);  // size of singular value to be considered as nonzero.
    //VectorXd vtr_vari_big = VectorXd::Zero(mtx_big.cols());
    //vtr_vari_big = svd.solve(vtr_bigger_delta);
    //return vtr_vari_big;

    MatrixXd BTB = mtx_big.transpose()*mtx_big;
    VectorXd BTD = mtx_big.transpose()*vtr_delta;

    MatrixXd mtx_Bigger = MatrixXd::Zero(BTB.rows()+mtx_constraints.rows(),BTB.cols()+mtx_constraints.rows());
    mtx_Bigger.block(0,0,BTB.rows(),BTB.cols()) = 2.*BTB;
    mtx_Bigger.block(BTB.rows(),0,mtx_constraints.rows(),mtx_constraints.cols()) = mtx_constraints;
    mtx_Bigger.block(0,BTB.cols(),mtx_constraints.cols(),mtx_constraints.rows()) = mtx_constraints.transpose();

    VectorXd vtr_bigger_delta = VectorXd::Zero(BTB.rows()+mtx_constraints.rows());
    vtr_bigger_delta.segment(0,BTB.cols()) = 2.*BTD;
    vtr_bigger_delta.segment(BTB.cols(),vtr_constraints_delta.size()) = vtr_constraints_delta;

    ComplexEigenSolver<MatrixXcd> eigensolver_bigger = ComplexEigenSolver<MatrixXcd>(mtx_Bigger);
    std::cout << "eigensolver_bigger.eigenvalues()(mtx_Bigger.cols()-1)  = " << eigensolver_bigger.eigenvalues()(mtx_Bigger.cols()-1) << std::endl;
    //svd_threshold = eigensolver_bigger.eigenvalues()(mtx_Bigger.cols()-1).real()/10000.;

    VectorXd vtr_vari_bigger = VectorXd::Zero(BTB.cols()+mtx_constraints.rows());
    //vtr_vari_bigger = mtx_Bigger.bdcSvd(ComputeThinU | ComputeThinV).solve(vtr_bigger_delta);
    //JacobiSVD<MatrixXd> svd(mtx_Bigger, ComputeThinU | ComputeThinV);
    BDCSVD<MatrixXd> svd(mtx_Bigger, ComputeThinU | ComputeThinV);
    //std::cout << "svd.singularValues()" << std::endl;
    //std::cout << svd.singularValues() << std::endl;
    svd.setThreshold(svd_threshold);  // size of singular value to be considered as nonzero.
    vtr_vari_bigger = svd.solve(vtr_bigger_delta);

    bool bad_result = false;
    for (int i=0;i<vtr_vari_bigger.size();i++)
    {
        if (isnan(vtr_vari_bigger(i))) bad_result = true;
    }
    if (bad_result)
    {
        svd.setThreshold(1e-3);
        vtr_vari_bigger = svd.solve(vtr_bigger_delta);
    }

    return vtr_vari_bigger.segment(0,BTB.cols());

}

MatrixXcd SmoothingComplexVectors(MatrixXcd mtx_input, int entry_start)
{
    MatrixXcd mtx_output = mtx_input;
    int col = mtx_input.cols()-entry_start;
    for (int row=1;row<mtx_input.rows()-1;row++)
    {
        mtx_output(row,col) = (mtx_input(row+1,col)+mtx_input(row,col)+mtx_input(row-1,col))/3.;
    }
    return mtx_output;
}

MatrixXcd GetTruncatedMatrix(MatrixXcd mtx_input, int rank_cutoff)
{
    MatrixXcd mtx_temp = mtx_input;
    ComplexEigenSolver<MatrixXcd> eigensolver_temp = ComplexEigenSolver<MatrixXcd>(mtx_temp);
    MatrixXcd mtx_eigenvector_temp = eigensolver_temp.eigenvectors();
    MatrixXcd mtx_eigenvector_inv_temp = eigensolver_temp.eigenvectors().inverse();
    MatrixXcd mtx_eigenvalue_temp = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    for (int row=0;row<mtx_temp.rows();row++)
    {
        for (int col=0;col<mtx_temp.cols();col++)
        {
            mtx_eigenvalue_temp(row,col) = 0.;
            if (row==col)
            {
                mtx_eigenvalue_temp(row,col) = eigensolver_temp.eigenvalues()(col);
            }
        }
    }
    mtx_eigenvalue_temp = CutoffEigenvalueMatrix(mtx_eigenvalue_temp, rank_cutoff);
    mtx_temp = mtx_eigenvector_temp*mtx_eigenvalue_temp*mtx_eigenvector_inv_temp;
    return mtx_temp;
}
int DetermineStableNumberOfEigenvalues(MatrixXcd mtx_input, MatrixXcd mtx_eigenvalue_input, MatrixXcd mtx_eigenvector_input, MatrixXcd mtx_eigenvector_inv_input)
{
    //return 1;
    ComplexEigenSolver<MatrixXcd> eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    int stable_number_max = 3;
    int stable_number = 0;
    MatrixXcd mtx_H_input = mtx_eigenvector_inv_input*mtx_eigenvector_input;
    for (int cutoff=1;cutoff<=stable_number_max;cutoff++)
    {
        double real_lambda = eigensolver_input.eigenvalues()(mtx_eigenvalue_input.rows()-cutoff).real();
        double imag_lambda = eigensolver_input.eigenvalues()(mtx_eigenvalue_input.rows()-cutoff).imag();
        double real_S = mtx_eigenvalue_input(mtx_eigenvalue_input.rows()-cutoff,mtx_eigenvalue_input.cols()-cutoff).real();
        double imag_S = mtx_eigenvalue_input(mtx_eigenvalue_input.rows()-cutoff,mtx_eigenvalue_input.cols()-cutoff).imag();
        double real_eta = mtx_H_input(mtx_eigenvalue_input.rows()-cutoff,mtx_eigenvalue_input.cols()-cutoff).real();
        if (abs(imag_lambda)/abs(real_lambda)>0.01) break;
        if (real_lambda<0.) break;
        if (real_S*real_eta<0.) break;
        stable_number = cutoff;
    }
    return stable_number;
}
MatrixXcd SpectralDecompositionMethod_v3(MatrixXcd mtx_input, int entry_start, int entry_size, double step_frac, bool isBlind)
{

    //SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_input);

    MatrixXcd mtx_q_init = mtx_eigenvector_init;
    MatrixXcd mtx_S = mtx_eigenvalue_init;
    MatrixXcd mtx_p_init = mtx_eigenvector_inv_init.transpose();
    MatrixXcd mtx_H_init = mtx_eigenvector_inv_init*mtx_eigenvector_init;

    MatrixXcd mtx_q_gamma = mtx_eigenvector_gamma;
    MatrixXcd mtx_p_gamma = mtx_eigenvector_inv_gamma.transpose();

    MatrixXcd mtx_q_dark = mtx_eigenvector_dark;
    MatrixXcd mtx_p_dark = mtx_eigenvector_inv_dark.transpose();

    MatrixXcd mtx_S_data = mtx_eigenvalue_data;
    MatrixXcd mtx_H_data = mtx_eigenvector_inv_data*mtx_eigenvector_data;

    int truncate_cutoff = 3;
    int dark_cutoff = DetermineStableNumberOfEigenvalues(mtx_input,mtx_eigenvalue_init, mtx_eigenvector_init, mtx_eigenvector_inv_init);
    int data_cutoff = DetermineStableNumberOfEigenvalues(mtx_data,mtx_eigenvalue_data, mtx_eigenvector_data, mtx_eigenvector_inv_data);
    truncate_cutoff = min(dark_cutoff,data_cutoff);
    MatrixXcd mtx_data_truncated = mtx_data;
    MatrixXcd mtx_input_truncated = mtx_input;
    if (TruncateNoise)
    {
        mtx_data_truncated = GetTruncatedMatrix(mtx_data, truncate_cutoff);
        mtx_input_truncated = GetTruncatedMatrix(mtx_input, truncate_cutoff);
    }

    TH2D Hist_Temp_Data = TH2D("Hist_Temp_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    fill2DHistogram(&Hist_Temp_Data,mtx_data_truncated);
    TH2D Hist_Temp_Dark = TH2D("Hist_Temp_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    fill2DHistogram(&Hist_Temp_Dark,mtx_input_truncated);
    int binx_lower = Hist_Temp_Data.GetXaxis()->FindBin(MSCL_cut_lower);
    int binx_blind = Hist_Temp_Data.GetXaxis()->FindBin(MSCL_cut_blind)-1;
    int binx_upper = Hist_Temp_Data.GetXaxis()->FindBin(1.)-1;
    int biny_lower = Hist_Temp_Data.GetYaxis()->FindBin(MSCW_cut_lower);
    int biny_blind = Hist_Temp_Data.GetYaxis()->FindBin(MSCW_cut_blind)-1;
    int biny_upper = Hist_Temp_Data.GetYaxis()->FindBin(1.)-1;
    double Data_CR_Integral = Hist_Temp_Data.Integral(binx_blind+1,binx_upper+1,biny_lower,biny_blind);
    double Dark_CR_Integral = Hist_Temp_Dark.Integral(binx_blind+1,binx_upper+1,biny_lower,biny_blind);


    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());

    bool DoRegularization = solution_w_regularizations;
    if (!isBlind) DoRegularization = false;
    int n_regularization = 1;
    if (!DoRegularization) n_regularization = 0;
    int row_size_big = mtx_input.rows()*mtx_input.cols() + n_regularization*entry_size;
    VectorXd vtr_Delta = VectorXd::Zero(row_size_big);
    MatrixXd mtx_Big = MatrixXd::Zero(row_size_big,2*entry_size*mtx_input.cols());
    //MatrixXd mtx_Constraint = MatrixXd::Zero((3*(NumberOfEigenvectors_Stable))*entry_size,2*entry_size*mtx_input.cols());
    //VectorXd vtr_Constraint_Delta = VectorXd::Zero((3*(NumberOfEigenvectors_Stable))*entry_size);
    MatrixXd mtx_Constraint = MatrixXd::Zero(2*3*entry_size,2*entry_size*mtx_input.cols());
    VectorXd vtr_Constraint_Delta = VectorXd::Zero(2*3*entry_size);
    for (int idx_k=0; idx_k<entry_size; idx_k++)
    {
        int nth_entry = idx_k + entry_start;
        for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
        {
            for (int idx_j=0; idx_j<mtx_input.cols(); idx_j++)
            {
                int idx_m = idx_j + mtx_input.cols()*idx_i;
                double weight = 1.;
                //weight = 1./max(1.,pow(mtx_data(idx_i,idx_j).real(),0.5)); // you need to freeze eigenvalues
                //if (!isBlind) weight = 1./max(1.,pow(mtx_data(idx_i,idx_j).real(),0.5));
                if (idx_i<binx_blind_global && idx_j<biny_blind_global)
                {
                    if (isBlind) weight = 0.; // blind gamma-ray region
                }
                vtr_Delta(idx_m) = weight*(mtx_data_truncated-mtx_input_truncated)(idx_i,idx_j).real();
                int idx_n = idx_j + mtx_input.cols()*idx_k;
                int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;

                mtx_Big(idx_m,idx_n) = weight*mtx_S(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry).real()*mtx_q_init(idx_i,mtx_input.rows()-nth_entry).real();
                mtx_Big(idx_m,idx_w) = weight*mtx_S(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry).real()*mtx_p_init(idx_j,mtx_input.rows()-nth_entry).real();
            }
        }
        if (DoRegularization)
        {
            double weight_s = 0.1*mtx_S(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry).real();
            double weight_p = 0.;
            double weight_q = 0.;
            double weight_norm = 0.;
            for (int idx_j=0; idx_j<mtx_input.cols(); idx_j++)
            {
                weight_norm += 1.;
                weight_p += pow(mtx_p_init(idx_j,mtx_input.rows()-nth_entry).real(),2);
                weight_q += pow(mtx_q_init(idx_j,mtx_input.rows()-nth_entry).real(),2);
            }
            weight_p = 1.0*weight_s*pow(weight_p/weight_norm,0.5);
            weight_q = 1.0*weight_s*pow(weight_q/weight_norm,0.5);
            for (int idx_j=0; idx_j<mtx_input.cols(); idx_j++)
            {
                int idx_m = mtx_input.cols()*mtx_input.rows() + n_regularization*idx_k;
                int idx_n = idx_j + mtx_input.cols()*idx_k;
                int idx_w = idx_j + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
                for (int idx_r=0; idx_r<n_regularization; idx_r++)
                {
                    vtr_Delta(idx_m+idx_r) = 0.;
                    mtx_Big(idx_m+idx_r,idx_w) = weight_s*mtx_p_init(idx_j,mtx_input.rows()-nth_entry).real();
                    mtx_Big(idx_m+idx_r,idx_n) = weight_s*mtx_q_init(idx_j,mtx_input.rows()-nth_entry).real();
                }
            }
        }
        for (int idx_l=0; idx_l<3; idx_l++)
        {
            int nth_entry2 = idx_l + 1;
            double weight_s1 = mtx_S(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry).real();
            double weight_s2 = mtx_S(mtx_input.rows()-nth_entry2,mtx_input.rows()-nth_entry2).real();
            double weight_s = 100.*mtx_S(mtx_input.rows()-1,mtx_input.rows()-1).real();
            for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
            {
                int idx_u = 2*idx_l + idx_k*2*3;
                int idx_n = idx_i + mtx_input.cols()*idx_k;
                int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
                vtr_Constraint_Delta(idx_u) = 0.;
                //if (nth_entry2==nth_entry) continue;
                //mtx_Constraint(idx_u+0,idx_w) = weight_s*mtx_q_gamma(idx_i,mtx_input.rows()-nth_entry2).real();
                //mtx_Constraint(idx_u+1,idx_n) = weight_s*mtx_p_gamma(idx_i,mtx_input.rows()-nth_entry2).real();
                if (nth_entry2>NumberOfEigenvectors_Stable) continue;
                if (nth_entry2==nth_entry)
                {
                    mtx_Constraint(idx_u+0,idx_w) = weight_s*mtx_q_init(idx_i,mtx_input.rows()-nth_entry2).real();
                    //mtx_Constraint(idx_u+1,idx_n) = weight_s*mtx_p_init(idx_i,mtx_input.rows()-nth_entry2).real();
                    //mtx_Constraint(idx_u,idx_w) = weight_s*mtx_p_init(idx_i,mtx_input.rows()-nth_entry2).real();
                    //mtx_Constraint(idx_u,idx_n) = weight_s*mtx_q_init(idx_i,mtx_input.rows()-nth_entry2).real();
                }
            }
        }
    }

    double lambda_0 = mtx_S(mtx_input.rows()-1,mtx_input.rows()-1).real();
    svd_threshold = 1./abs(1.*lambda_0);

    VectorXd vtr_vari_big = VectorXd::Zero(2*entry_size*mtx_input.cols());
    if (solution_w_constraints && entry_start>0) 
    {
        vtr_vari_big = SolutionWithConstraints(mtx_Big, mtx_Constraint, vtr_Delta, vtr_Constraint_Delta);
        //VectorXd vtr_should_be_zero = mtx_Constraint*vtr_vari_big;
        //for (int i=0;i<vtr_should_be_zero.size();i++)
        //{
        //    if (isnan(vtr_should_be_zero(i))) return mtx_input;
        //    if (abs(vtr_should_be_zero(i))>0.1) return mtx_input;
        //}
    }
    else
    {
        //vtr_vari_big = (mtx_Big.transpose()*mtx_Big).ldlt().solve(mtx_Big.transpose()*vtr_Delta);
        //vtr_vari_big = mtx_Big.colPivHouseholderQr().solve(vtr_Delta);
        //vtr_vari_big = mtx_Big.bdcSvd(ComputeThinU | ComputeThinV).solve(vtr_Delta);
        BDCSVD<MatrixXd> svd(mtx_Big, ComputeThinU | ComputeThinV);
        svd.setThreshold(svd_threshold);  // size of singular value to be considered as nonzero.
        vtr_vari_big = svd.solve(vtr_Delta);
    }
    //vtr_vari_big = (mtx_Big.transpose()*mtx_Big).inverse()*mtx_Big.transpose()*vtr_Delta;


    MatrixXcd mtx_p_final = mtx_p_init;
    MatrixXcd mtx_q_final = mtx_q_init;
    MatrixXd mtx_p_vari = MatrixXd::Zero(mtx_p_init.rows(),mtx_p_init.cols());
    MatrixXd mtx_q_vari = MatrixXd::Zero(mtx_q_init.rows(),mtx_q_init.cols());
    for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
    {
        for (int idx_j=0; idx_j<mtx_input.cols(); idx_j++)
        {
            for (int idx_k=0; idx_k<entry_size; idx_k++)
            {
                int idx_m = idx_j + mtx_input.cols()*idx_i;
                int idx_n = idx_j + mtx_input.cols()*idx_k;
                int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
                int nth_entry = idx_k + entry_start;
                mtx_p_vari(idx_j,mtx_p_vari.cols()-nth_entry) = vtr_vari_big(idx_n);
                mtx_q_vari(idx_i,mtx_q_vari.cols()-nth_entry) = vtr_vari_big(idx_w);
            }
        }
    }

    //for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
    //{
    //    for (int idx_j=0; idx_j<mtx_input.cols(); idx_j++)
    //    {
    //        for (int idx_k=0; idx_k<entry_size; idx_k++)
    //        {
    //            int idx_m = idx_j + mtx_input.cols()*idx_i;
    //            int idx_n = idx_j + mtx_input.cols()*idx_k;
    //            int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
    //            int nth_entry = idx_k + entry_start;
    //            if (isnan(mtx_p_vari(idx_j,mtx_p_vari.cols()-nth_entry))) return mtx_input;
    //            if (isnan(mtx_q_vari(idx_i,mtx_q_vari.cols()-nth_entry))) return mtx_input;
    //        }
    //    }
    //}

    if (isBlind)
    {
        std::cout << "mtx_S:" << std::endl;
        std::cout << mtx_S.block(mtx_S.rows()-3,mtx_S.cols()-3,3,3) << std::endl;

        std::cout << "mtx_H_init:" << std::endl;
        std::cout << mtx_H_init.block(mtx_H_init.rows()-3,mtx_H_init.cols()-3,3,3) << std::endl;

        MatrixXcd mtx_H_vari = mtx_p_vari.transpose()*mtx_q_init + mtx_p_init.transpose()*mtx_q_vari;
        std::cout << "mtx_H_vari:" << std::endl;
        std::cout << mtx_H_vari.block(mtx_H_vari.rows()-3,mtx_H_vari.cols()-3,3,3) << std::endl;

        MatrixXcd mtx_P_vari = mtx_p_vari.transpose()*mtx_p_init;
        std::cout << "mtx_P_vari:" << std::endl;
        std::cout << mtx_P_vari.block(mtx_P_vari.rows()-3,mtx_P_vari.cols()-3,3,3) << std::endl;

        MatrixXcd mtx_Q_vari = mtx_q_vari.transpose()*mtx_q_init;
        std::cout << "mtx_Q_vari:" << std::endl;
        std::cout << mtx_Q_vari.block(mtx_Q_vari.rows()-3,mtx_Q_vari.cols()-3,3,3) << std::endl;
    }

    //mtx_p_vari = SmoothingRealVectors(mtx_p_vari);
    //mtx_q_vari = SmoothingRealVectors(mtx_q_vari);

    mtx_output = mtx_input;
    mtx_output += mtx_q_init*mtx_S*(mtx_p_vari.transpose()*step_frac);
    mtx_output += (mtx_q_vari*step_frac)*mtx_S*mtx_p_init.transpose();

    mtx_eigenvector_init += mtx_q_vari*step_frac;
    mtx_eigenvector_inv_init += mtx_p_vari.transpose()*step_frac;

    //MatrixXcd mtx_H_final = mtx_eigenvector_inv_init*mtx_eigenvector_init;
    //std::cout << "mtx_H_final:" << std::endl;
    //std::cout << mtx_H_final.block(mtx_H_init.rows()-4,mtx_H_init.cols()-4,4,4) << std::endl;

    return mtx_output;

}
bool CheckIfEigenvalueMakeSense(MatrixXcd mtx_input, double init_chi2, int rank, bool isBlind)
{
    //ComplexEigenSolver<MatrixXcd> eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    //for (int i=1;i<=rank;i++)
    //{
    //    double lambda_real = eigensolver_input.eigenvalues()(mtx_input.cols()-i).real();
    //    double lambda_imag = eigensolver_input.eigenvalues()(mtx_input.cols()-i).imag();
    //    if (lambda_real<0.0)
    //    {
    //        std::cout << "break (lambda_real is negative)" << std::endl;
    //        return false;
    //    }
    //    if (abs(lambda_imag)/abs(lambda_real)>0.01)
    //    {
    //        std::cout << "break (lambda_imag too big.)" << std::endl;
    //        return false;
    //    }
    //}

    //init_chi2 = GetChi2Function(mtx_data_bkgd,isBlind);
    double current_chi2 = GetChi2Function(mtx_input,isBlind);
    //if (current_chi2>1.*init_chi2) 
    //{
    //    std::cout << "break (chi2 increasing.)" << std::endl;
    //    return false;
    //}
    if (isnan(current_chi2))
    {
        std::cout << "break (chi2 is nan.)" << std::endl;
        return false;
    }

    return true;
}
void GetTruncatedHistogram(TH2D* hist_input, TH2D* hist_output, int rank_cutoff)
{
    MatrixXcd mtx_temp = fillMatrix(hist_input);
    ComplexEigenSolver<MatrixXcd> eigensolver_temp = ComplexEigenSolver<MatrixXcd>(mtx_temp);
    MatrixXcd mtx_eigenvector_temp = eigensolver_temp.eigenvectors();
    MatrixXcd mtx_eigenvector_inv_temp = eigensolver_temp.eigenvectors().inverse();
    MatrixXcd mtx_eigenvalue_temp = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    for (int row=0;row<mtx_temp.rows();row++)
    {
        for (int col=0;col<mtx_temp.cols();col++)
        {
            mtx_eigenvalue_temp(row,col) = 0.;
            if (row==col)
            {
                mtx_eigenvalue_temp(row,col) = eigensolver_temp.eigenvalues()(col);
            }
        }
    }
    mtx_eigenvalue_temp = CutoffEigenvalueMatrix(mtx_eigenvalue_temp, rank_cutoff);
    mtx_temp = mtx_eigenvector_temp*mtx_eigenvalue_temp*mtx_eigenvector_inv_temp;
    fill2DHistogram(hist_output,mtx_temp);
}
void LeastSquareSolutionMethod(bool DoSequential, bool isBlind)
{

    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_data_bkgd(row,col) = 0.;
            mtx_eigenvector(row,col) = 0.;
            mtx_eigenvector_inv(row,col) = 0.;
            mtx_eigenvalue(row,col) = 0.;
        }
    }

    SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_dark);
    mtx_data_bkgd = mtx_dark;

    if (eigensolver_data.eigenvalues()(mtx_data.cols()-1)==0.) return;
    if (eigensolver_dark.eigenvalues()(mtx_dark.cols()-1)==0.) return;

    mtx_eigenvector = mtx_eigenvector_init;
    mtx_eigenvalue = mtx_eigenvalue_init;
    mtx_eigenvector_inv = mtx_eigenvector_inv_init;
    //mtx_data_bkgd = mtx_eigenvector_init*mtx_eigenvalue_init*mtx_eigenvector_inv_init;

    std::cout << "initial chi2 in CR = " << GetChi2Function(mtx_dark,isBlind) << std::endl;
    double init_chi2 = GetChi2Function(mtx_dark,isBlind);

    MatrixXcd mtx_temp = mtx_dark;
    double eigenvalue_dark_real = 0.;
    double eigenvalue_dark_imag = 0.;
    double eigenvalue_data_real = 0.;
    double eigenvalue_data_imag = 0.;
    double step_frac= 1.0;
    double imag_real_ratio = 1./100.;

    int dark_cutoff = DetermineStableNumberOfEigenvalues(mtx_dark,mtx_eigenvalue_init, mtx_eigenvector_init, mtx_eigenvector_inv_init);
    int data_cutoff = DetermineStableNumberOfEigenvalues(mtx_data,mtx_eigenvalue_data, mtx_eigenvector_data, mtx_eigenvector_inv_data);
    NumberOfEigenvectors_Stable = min(dark_cutoff,data_cutoff);
    std::cout << "NumberOfEigenvectors_Stable = " << NumberOfEigenvectors_Stable << std::endl;
    if (DoSequential)
    {
        for (int rank=1; rank<=NumberOfEigenvectors_Stable; rank++)
        {
            mtx_temp = SpectralDecompositionMethod_v3(mtx_data_bkgd, rank, 1, 1.0, isBlind);
            if (CheckIfEigenvalueMakeSense(mtx_temp, init_chi2, rank, isBlind))
            {
                mtx_data_bkgd = mtx_temp;
            }
        }
    }
    else
    {
        if (NumberOfEigenvectors_Stable==0)
        {
            return;
        }
        mtx_temp = SpectralDecompositionMethod_v3(mtx_data_bkgd, 1, NumberOfEigenvectors_Stable, 1.0, isBlind);
        if (CheckIfEigenvalueMakeSense(mtx_temp, init_chi2, NumberOfEigenvectors_Stable, isBlind))
        {
            mtx_data_bkgd = mtx_temp;
            return;
        }
        if (NumberOfEigenvectors_Stable-1>0)
        {
            mtx_temp = SpectralDecompositionMethod_v3(mtx_data_bkgd, 1, NumberOfEigenvectors_Stable-1, 1.0, isBlind);
            if (CheckIfEigenvalueMakeSense(mtx_temp, init_chi2, NumberOfEigenvectors_Stable-1, isBlind))
            {
                mtx_data_bkgd = mtx_temp;
                return;
            }
        }
        if (NumberOfEigenvectors_Stable-2>0)
        {
            mtx_temp = SpectralDecompositionMethod_v3(mtx_data_bkgd, 1, NumberOfEigenvectors_Stable-2, 1.0, isBlind);
            if (CheckIfEigenvalueMakeSense(mtx_temp, init_chi2, NumberOfEigenvectors_Stable-2, isBlind))
            {
                mtx_data_bkgd = mtx_temp;
                return;
            }
        }
    }

}
void GetReducedEigenvalueMatrix(int rank_cutoff)
{
    for (int i=0;i<mtx_data.cols();i++)
    {
        for (int j=0;j<mtx_data.rows();j++)
        {
            if (i==j && i==mtx_data.cols()-rank_cutoff-1) 
            {
                mtx_eigenval_data_redu(i,j) = eigensolver_data.eigenvalues()(i);
            }
            else
            {
                mtx_eigenval_data_redu(i,j) = 0;
            }
        }
    }
}
void GetCRReplacedMatrix(TH2D* hist_data, TH2D* hist_dark)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    TH2D hist_temp = TH2D("hist_temp","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            if (bx<binx_blind && by<biny_blind)
            {
                double dark = hist_dark->GetBinContent(bx,by);
                hist_temp.SetBinContent(bx,by,dark);
            }
            else
            {
                double data = hist_data->GetBinContent(bx,by);
                hist_temp.SetBinContent(bx,by,data);
            }
        }
    }
    hist_dark->Reset();
    hist_dark->Add(&hist_temp);
}
void GetEigenvaluesReplacedMatrix(TH2D* hist_input, TH2D* hist_output)
{
    MatrixXcd mtx_input = fillMatrix(hist_input);
    ComplexEigenSolver<MatrixXcd> eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    MatrixXcd mtx_output = fillMatrix(hist_output);
    ComplexEigenSolver<MatrixXcd> eigensolver_output = ComplexEigenSolver<MatrixXcd>(mtx_output);
    MatrixXcd mtx_eigenvector_output = eigensolver_output.eigenvectors();
    MatrixXcd mtx_eigenvector_inv_output = eigensolver_output.eigenvectors().inverse();
    MatrixXcd mtx_eigenvalue_output = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    for (int row=0;row<mtx_output.rows();row++)
    {
        for (int col=0;col<mtx_output.cols();col++)
        {
            mtx_eigenvalue_output(row,col) = 0.;
            if (row==col)
            {
                double ratio = eigensolver_input.eigenvalues()(col).real()/eigensolver_input.eigenvalues()(mtx_output.rows()-1).real();
                double new_eigenvalue = eigensolver_output.eigenvalues()(mtx_output.rows()-1).real()*ratio;
                mtx_eigenvalue_output(row,col) = new_eigenvalue;
            }
        }
    }
    mtx_output = mtx_eigenvector_output*mtx_eigenvalue_output*mtx_eigenvector_inv_output;
    fill2DHistogram(hist_output,mtx_output);
}
void GetNoiseReplacedMatrix(TH2D* hist_data, TH2D* hist_dark, bool isTight)
{
    TH2D hist_data_temp = TH2D("hist_data_temp","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_data_noise = TH2D("hist_data_noise","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_dark_temp = TH2D("hist_dark_temp","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);

    int cutoff = 3;
    //MatrixXcd mtx_S = mtx_eigenvalue_dark;
    //MatrixXcd mtx_H_init = mtx_eigenvector_inv_dark*mtx_eigenvector_dark;
    //double sign_3 = mtx_S(mtx_S.rows()-3,mtx_S.rows()-3).real()*mtx_H_init(mtx_S.rows()-3,mtx_S.rows()-3).real();
    //double sign_2 = mtx_S(mtx_S.rows()-2,mtx_S.rows()-2).real()*mtx_H_init(mtx_S.rows()-2,mtx_S.rows()-2).real();
    //double sign_1 = mtx_S(mtx_S.rows()-1,mtx_S.rows()-1).real()*mtx_H_init(mtx_S.rows()-1,mtx_S.rows()-1).real();
    //if (sign_3<0.) cutoff = 2;
    //if (sign_2<0.) cutoff = 1;
    //if (sign_1<0.) cutoff = 0;

    GetTruncatedHistogram(hist_data,&hist_data_temp,cutoff);
    GetTruncatedHistogram(hist_dark,&hist_dark_temp,cutoff);
    hist_data_noise.Add(hist_data);
    hist_data_noise.Add(&hist_data_temp,-1.);
    hist_dark->Reset();
    hist_dark->Add(&hist_dark_temp);
    hist_dark->Add(&hist_data_noise);
}
void NormalizaDarkMatrix(TH2D* hist_data, TH2D* hist_dark)
{
    int binx_lower = hist_data->GetXaxis()->FindBin(MSCL_cut_lower);
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind)-1;
    int binx_upper = hist_data->GetNbinsX();
    int biny_lower = hist_data->GetYaxis()->FindBin(MSCW_cut_lower);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind)-1;
    int biny_upper = hist_data->GetNbinsY();
    double Data_SR_Integral = hist_data->Integral(binx_lower,binx_blind,biny_lower,biny_blind);
    double Data_Integral = hist_data->Integral();
    double Data_CR_Integral_1 = hist_data->Integral(binx_blind,binx_upper,biny_blind,biny_upper);
    //double Data_CR_Integral_2 = Data_Integral-Data_SR_Integral-Data_CR_Integral_1;
    double Data_CR_Integral_2 = Data_Integral-Data_SR_Integral;
    double Dark_SR_Integral = hist_dark->Integral(binx_lower,binx_blind,biny_lower,biny_blind);
    double Dark_Integral = hist_dark->Integral();
    double Dark_CR_Integral_1 = hist_dark->Integral(binx_blind,binx_upper,biny_blind,biny_upper);
    //double Dark_CR_Integral_2 = Dark_Integral-Dark_SR_Integral-Dark_CR_Integral_1;
    double Dark_CR_Integral_2 = Dark_Integral-Dark_SR_Integral;
    if (Dark_CR_Integral_1!=0)
    {
        double dark_scale = Data_CR_Integral_2/Dark_CR_Integral_2;
        hist_dark->Scale(dark_scale);
    }
    else
    {
        hist_data->Scale(0.);
        hist_dark->Scale(0.);
    }
}
void AlterDarkMatrix(TH2D* hist_data, TH2D* hist_dark, TH2D* hist_dark_alter)
{

    hist_dark_alter->Reset();
    hist_dark_alter->Add(hist_dark);
    if (UseTruncatedONData)
    {
        GetTruncatedHistogram(hist_data,hist_dark_alter,3);
    }
    if (UseReplacedONData)
    {
        GetCRReplacedMatrix(hist_data,hist_dark_alter);
    }

}
MatrixXcd MinChi2Method(TH1D* Hist_Chi2, TH2D* Hist_Data, TH2D* Hist_Dark)
{

    mtx_dark = fillMatrix(Hist_Dark);
    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    mtx_data_bkgd = fillMatrix(Hist_Dark);
    MatrixXcd mtx_best_bkgd = mtx_dark;

    TH2D Hist_Temp_Dark = TH2D("Hist_Temp_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D Hist_Temp_Dark_Alt = TH2D("Hist_Temp_Dark_Alt","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);

    TH2D Hist_Temp2_Gamma = TH2D("Hist_Temp2_Gamma","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    fill2DHistogram(&Hist_Temp2_Gamma,mtx_gamma_raw);
    std::cout << "Hist_Temp2_Gamma.Integral() = " << Hist_Temp2_Gamma.Integral() << std::endl;
    TH2D Hist_Temp_Data = TH2D("Hist_Temp_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    Hist_Temp_Data.Add(Hist_Data);
    Hist_Temp_Dark.Add(Hist_Dark);

    int binx_lower = Hist_Temp_Data.GetXaxis()->FindBin(MSCL_cut_lower);
    int binx_blind = Hist_Temp_Data.GetXaxis()->FindBin(MSCL_cut_blind)-1;
    int binx_upper = Hist_Temp_Data.GetXaxis()->FindBin(1.)-1;
    int biny_lower = Hist_Temp_Data.GetYaxis()->FindBin(MSCW_cut_lower);
    int biny_blind = Hist_Temp_Data.GetYaxis()->FindBin(MSCW_cut_blind)-1;
    int biny_upper = Hist_Temp_Data.GetYaxis()->FindBin(1.)-1;

    TH2D Hist_Temp_Bkgd = TH2D("Hist_Temp_Bkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    double Data_SR_Integral = Hist_Temp_Data.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
    double Data_Integral = Hist_Temp_Data.Integral();
    double Dark_Integral = Hist_Temp_Dark.Integral();
    double Gamma_SR_Integral = Hist_Temp2_Gamma.Integral(binx_lower,binx_blind,biny_lower,biny_blind);

    if (Data_Integral==0.) return mtx_best_bkgd;
    if (Dark_Integral==0.) return mtx_best_bkgd;
    if (Gamma_SR_Integral==0.) return mtx_best_bkgd;

    Hist_Temp_Bkgd.Reset();
    double Dark_SR_Integral = Hist_Temp_Dark.Integral(binx_lower,binx_blind,biny_lower,biny_blind);

    double chi2_best = 1e20;
    for (int bin=0;bin<Hist_Chi2->GetNbinsX();bin++)
    {
        Hist_Temp_Dark_Alt.Reset();
        Hist_Temp_Data.Reset();
        Hist_Temp_Data.Add(Hist_Data);
        Hist_Temp_Dark.Reset();
        Hist_Temp_Dark.Add(Hist_Dark);
        Hist_Temp_Bkgd.Reset();
        Hist_Temp2_Gamma.Reset();
        fill2DHistogram(&Hist_Temp2_Gamma,mtx_gamma_raw);
        double scale = Hist_Chi2->GetBinCenter(bin+1);
        double delta = Dark_SR_Integral*scale;
        Hist_Temp2_Gamma.Scale((delta)/Gamma_SR_Integral);
        Hist_Temp_Data.Add(&Hist_Temp2_Gamma,-1.);
        mtx_gamma = fillMatrix(&Hist_Temp2_Gamma);
        mtx_data = fillMatrix(&Hist_Temp_Data);
        NormalizaDarkMatrix(&Hist_Temp_Data,&Hist_Temp_Dark);
        AlterDarkMatrix(Hist_Data, &Hist_Temp_Dark, &Hist_Temp_Dark_Alt);
        mtx_dark = fillMatrix(&Hist_Temp_Dark_Alt);
        eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
        LeastSquareSolutionMethod(false, true);
        fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
        double chi2 = GetChi2Function(mtx_data_bkgd,true);
        if (chi2<chi2_best)
        {
            mtx_best_bkgd = mtx_data_bkgd;
        }
        Hist_Chi2->SetBinContent(bin+1,chi2);
    }
    int min_error_bin = Hist_Chi2->GetMinimumBin();
    Hist_Temp_Dark_Alt.Reset();
    Hist_Temp_Data.Reset();
    Hist_Temp_Data.Add(Hist_Data);
    Hist_Temp_Dark.Reset();
    Hist_Temp_Dark.Add(Hist_Dark);
    Hist_Temp_Bkgd.Reset();
    Hist_Temp2_Gamma.Reset();
    fill2DHistogram(&Hist_Temp2_Gamma,mtx_gamma_raw);
    double scale = Hist_Chi2->GetBinCenter(min_error_bin);
    double delta = Dark_SR_Integral*scale;
    Hist_Temp2_Gamma.Scale((delta)/Gamma_SR_Integral);
    Hist_Temp_Data.Add(&Hist_Temp2_Gamma,-1.);
    mtx_gamma = fillMatrix(&Hist_Temp2_Gamma);
    mtx_data = fillMatrix(&Hist_Temp_Data);
    NormalizaDarkMatrix(&Hist_Temp_Data,&Hist_Temp_Dark);
    AlterDarkMatrix(Hist_Data, &Hist_Temp_Dark, &Hist_Temp_Dark_Alt);
    mtx_dark = fillMatrix(&Hist_Temp_Dark_Alt);
    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    LeastSquareSolutionMethod(false, true);
    fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);

    return mtx_data_bkgd;

}
void MakePrediction(string target_data, double tel_elev_lower_input, double tel_elev_upper_input, int MJD_start_cut, int MJD_end_cut, double input_theta2_cut_lower, double input_theta2_cut_upper, double SignalNorm, bool isON)
{

    TH1::SetDefaultSumw2();
    sprintf(target, "%s", target_data.c_str());

    TString ONOFF_tag;
    if (isON) 
    {
        source_theta2_cut = 0.;
        ONOFF_tag = "ON";
    }
    else
    {
        ONOFF_tag = "OFF";
    }

    if (MJD_start_cut!=0 || MJD_end_cut!=0)
    {
        sprintf(mjd_cut_tag, "_MJD%dto%d", MJD_start_cut, MJD_end_cut);
    }
    camera_theta2_cut_lower = input_theta2_cut_lower;
    camera_theta2_cut_upper = input_theta2_cut_upper;
    sprintf(theta2_cut_tag, "_Theta2%dto%d", int(camera_theta2_cut_lower), int(camera_theta2_cut_upper));
    sprintf(signal_tag, "_S%d", int(SignalNorm*100.));
    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;
    sprintf(elev_cut_tag, "_TelElev%dto%d", int(TelElev_lower), int(TelElev_upper));
    MSCW_cut_blind = MSCW_cut_moderate;
    MSCL_cut_blind = MSCL_cut_moderate;
    if (TString(target).Contains("Crab"))
    {
        if (source_theta2_cut==0.)
        {
            MSCW_cut_blind = MSCW_cut_loose;
            MSCL_cut_blind = MSCL_cut_loose;
        }
    }
    if (TString(target).Contains("Mrk421"))
    {
        if (source_theta2_cut==0.)
        {
            MSCW_cut_blind = MSCW_cut_loose;
            MSCL_cut_blind = MSCL_cut_loose;
        }
    }

    mean_tele_point_ra = 0.;
    mean_tele_point_dec = 0.;
    pair<double,double> source_ra_dec = GetSourceRaDec(TString(target));
    mean_tele_point_ra = source_ra_dec.first;
    mean_tele_point_dec = source_ra_dec.second;
    pair<double,double> tele_point_l_b = ConvertRaDecToGalactic(mean_tele_point_ra, mean_tele_point_dec);
    mean_tele_point_l = tele_point_l_b.first;
    mean_tele_point_b = tele_point_l_b.second;


    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    TH1D Hist_ErecS_fine = TH1D("Hist_ErecS_fine","",N_energy_fine_bins,energy_fine_bins);
    TH1D Hist_EffArea = TH1D("Hist_EffArea","",N_energy_fine_bins,energy_fine_bins);

    SizeSecondMax_Cut = 600.;
    if (TString(target).Contains("V4")) SizeSecondMax_Cut = 400.;
    if (TString(target).Contains("V5")) SizeSecondMax_Cut = 400.;

    vector<int> group_size;
    for (int e=0;e<N_energy_bins;e++) 
    {
        group_size.push_back(0);
    }

    vector<string>* Data_runlist_name_ptr = new std::vector<string>(10);
    vector<int>* Data_runlist_number_ptr = new std::vector<int>(10);
    vector<int>* Data_runlist_MJD_ptr = new std::vector<int>(10);
    vector<double>* Data_runlist_elev_ptr = new std::vector<double>(10);
    vector<double>* Data_runlist_exposure_ptr = new std::vector<double>(10);
    vector<string>* roi_name_ptr = new std::vector<string>(10);
    TFile InputDataFile("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+"_"+ONOFF_tag+".root");
    TTree* InfoTree_ptr = nullptr;
    InfoTree_ptr = (TTree*) InputDataFile.Get("InfoTree");
    InfoTree_ptr->SetBranchAddress("Data_runlist_name",&Data_runlist_name_ptr);
    InfoTree_ptr->SetBranchAddress("Data_runlist_number",&Data_runlist_number_ptr);
    InfoTree_ptr->SetBranchAddress("Data_runlist_MJD",&Data_runlist_MJD_ptr);
    InfoTree_ptr->SetBranchAddress("Data_runlist_elev",&Data_runlist_elev_ptr);
    InfoTree_ptr->SetBranchAddress("Data_runlist_exposure",&Data_runlist_exposure_ptr);
    InfoTree_ptr->SetBranchAddress("roi_name",&roi_name_ptr);
    InfoTree_ptr->GetEntry(0);

    int FirstRun = 0;
    int LastRun = Data_runlist_name_ptr->size();


    std::cout << "Working on OFF data..." << std::endl;
    vector<vector<TH2D>> Hist_OffData_MSCLW;
    vector<vector<TH2D>> Hist_OffBkgd_MSCLW;
    vector<vector<TH2D>> Hist_OnSyst_MSCLW;
    vector<vector<TH1D>> Hist_OnSyst_Chi2;
    vector<vector<TH1D>> Hist_OneGroup_OnSyst_Chi2;
    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
    {
        std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "sample = " << nth_sample << std::endl;
        char sample_tag[50];
        sprintf(sample_tag, "%i", nth_sample);
        vector<TH2D> Hist_OffData_OneSample_MSCLW;
        vector<TH2D> Hist_OffBkgd_OneSample_MSCLW;
        vector<TH2D> Hist_OnSyst_OneSample_MSCLW;
        vector<TH1D> Hist_OnSyst_OneSample_Chi2;
        vector<TH1D> Hist_OneGroup_OnSyst_OneSample_Chi2;
        for (int e=0;e<N_energy_bins;e++) 
        {
            std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
            std::cout << "energy = " << energy_bins[e] << std::endl;
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            NumberOfEigenvectors = N_max_ranks_func_E[e];

            MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
            MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
            N_bins_for_deconv = N_bins_for_deconv_func_E[e];
            ResetMatrixDimension();

            Hist_OffData_OneSample_MSCLW.push_back(TH2D("Hist_OffData2_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_OffBkgd_OneSample_MSCLW.push_back(TH2D("Hist_OffBkgd_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_OnSyst_OneSample_MSCLW.push_back(TH2D("Hist_OnSyst_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_OnSyst_OneSample_Chi2.push_back(TH1D("Hist_OnSyst_Chi2_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",40,0.,0.2));
            Hist_OneGroup_OnSyst_OneSample_Chi2.push_back(TH1D("Hist_OneGroup_OnSyst_Chi2_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",40,0.,0.2));
        }
        Hist_OffData_MSCLW.push_back(Hist_OffData_OneSample_MSCLW);
        Hist_OffBkgd_MSCLW.push_back(Hist_OffBkgd_OneSample_MSCLW);
        Hist_OnSyst_MSCLW.push_back(Hist_OnSyst_OneSample_MSCLW);
        Hist_OnSyst_Chi2.push_back(Hist_OnSyst_OneSample_Chi2);
        Hist_OneGroup_OnSyst_Chi2.push_back(Hist_OneGroup_OnSyst_OneSample_Chi2);
    }

    vector<vector<TH2D>> Hist_OneGroup_OffData_MSCLW;
    vector<vector<TH2D>> Hist_OneGroup_OffDark_MSCLW;
    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
    {
        char sample_tag[50];
        sprintf(sample_tag, "%i", nth_sample);
        vector<TH2D> Hist_OneSample_OffData_MSCLW;
        vector<TH2D> Hist_OneSample_OffDark_MSCLW;
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            NumberOfEigenvectors = N_max_ranks_func_E[e];

            MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
            MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
            N_bins_for_deconv = N_bins_for_deconv_func_E[e];
            ResetMatrixDimension();

            Hist_OneSample_OffData_MSCLW.push_back(TH2D("Hist_OneSample_OffData_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_OneSample_OffDark_MSCLW.push_back(TH2D("Hist_OneSample_OffDark_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        }
        Hist_OneGroup_OffData_MSCLW.push_back(Hist_OneSample_OffData_MSCLW);
        Hist_OneGroup_OffDark_MSCLW.push_back(Hist_OneSample_OffDark_MSCLW);
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "energy = " << energy_bins[e] << std::endl;
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));

        NumberOfEigenvectors = N_max_ranks_func_E[e];

        MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
        MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
        N_bins_for_deconv = N_bins_for_deconv_func_E[e];
        ResetMatrixDimension();

        for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
        {
            char sample2_tag[50];
            sprintf(sample2_tag, "%i", nth_sample);

            for (int on_run=FirstRun;on_run<LastRun;on_run++)
            {

                bool use_this_run = true;
                if (MJD_start_cut!=0 || MJD_end_cut!=0)
                {
                    if (MJD_start_cut>Data_runlist_MJD_ptr->at(on_run)) use_this_run = false;
                    if (MJD_end_cut<Data_runlist_MJD_ptr->at(on_run)) use_this_run = false;
                }

                if (use_this_run)
                {
                    char sample_tag[50];
                    sprintf(sample_tag, "%i", on_run);

                    int binx_lower = Hist_OffData_MSCLW.at(0).at(e).GetXaxis()->FindBin(MSCL_cut_lower);
                    binx_blind_global = Hist_OffData_MSCLW.at(0).at(e).GetXaxis()->FindBin(MSCL_cut_blind)-1;
                    int binx_upper = Hist_OffData_MSCLW.at(0).at(e).GetXaxis()->FindBin(1.)-1;
                    int biny_lower = Hist_OffData_MSCLW.at(0).at(e).GetYaxis()->FindBin(MSCW_cut_lower);
                    biny_blind_global = Hist_OffData_MSCLW.at(0).at(e).GetYaxis()->FindBin(MSCW_cut_blind)-1;
                    int biny_upper = Hist_OffData_MSCLW.at(0).at(e).GetYaxis()->FindBin(1.)-1;

                    TString hist_name;
                    hist_name  = "Hist_OnDark_MSCLW_R"+TString(sample_tag)+"_V"+TString(sample2_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                    Hist_OneGroup_OffData_MSCLW.at(nth_sample).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                    group_size.at(e) += ((TH2D*)InputDataFile.Get(hist_name))->Integral();
                    for (int nth_other_sample=0;nth_other_sample<n_dark_samples;nth_other_sample++)
                    {
                        if (nth_sample==nth_other_sample) continue;
                        char sample2_tag[50];
                        sprintf(sample2_tag, "%i", nth_other_sample);
                        hist_name  = "Hist_OnDark_MSCLW_R"+TString(sample_tag)+"_V"+TString(sample2_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                        Hist_OneGroup_OffDark_MSCLW.at(nth_sample).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                    }

                }

                if (group_size.at(e)>=group_size_limit[e]*1000 || on_run==Data_runlist_number_ptr->size()-1)
                {
                    mtx_data = fillMatrix(&Hist_OneGroup_OffData_MSCLW.at(nth_sample).at(e));
                    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
                    NormalizaDarkMatrix(&Hist_OneGroup_OffData_MSCLW.at(nth_sample).at(e), &Hist_OneGroup_OffDark_MSCLW.at(nth_sample).at(e));
                    TH2D hist_dark_alter = TH2D("hist_dark_alter","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
                    AlterDarkMatrix(&Hist_OneGroup_OffData_MSCLW.at(nth_sample).at(e), &Hist_OneGroup_OffDark_MSCLW.at(nth_sample).at(e), &hist_dark_alter);
                    mtx_dark = fillMatrix(&hist_dark_alter);
                    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
                    mtx_data_bkgd = mtx_dark;
                    //LeastSquareSolutionMethod(false, true);
                    TH2D Hist_Temp_Bkgd = TH2D("Hist_Temp_Bkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
                    fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
                    Hist_OffData_MSCLW.at(nth_sample).at(e).Add(&Hist_OneGroup_OffData_MSCLW.at(nth_sample).at(e));
                    Hist_OffBkgd_MSCLW.at(nth_sample).at(e).Add(&Hist_Temp_Bkgd);

                    Hist_OneGroup_OffData_MSCLW.at(nth_sample).at(e).Reset();
                    Hist_OneGroup_OffDark_MSCLW.at(nth_sample).at(e).Reset();
                    group_size.at(e) = 0;

                }
            }

        }
    }

    std::cout << "Working on ON data..." << std::endl;


    vector<vector<TH2D>> Hist_OneGroup_Dark_MSCLW;
    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
    {
        char sample_tag[50];
        sprintf(sample_tag, "%i", nth_sample);
        vector<TH2D> Hist_OneSample_Dark_MSCLW;
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            NumberOfEigenvectors = N_max_ranks_func_E[e];

            MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
            MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
            N_bins_for_deconv = N_bins_for_deconv_func_E[e];
            ResetMatrixDimension();

            Hist_OneSample_Dark_MSCLW.push_back(TH2D("Hist_OneSample_Dark_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        }
        Hist_OneGroup_Dark_MSCLW.push_back(Hist_OneSample_Dark_MSCLW);
    }
    vector<TH2D> Hist_OneGroup_Data_MSCLW;
    vector<TH2D> Hist_OneGroup_Bkgd_MSCLW;
    vector<TH2D> Hist_OneGroup_Bkgd_Unblind_wGamma_MSCLW;
    vector<TH2D> Hist_OneGroup_Bkgd_Unblind_woGamma_MSCLW;
    vector<TH2D> Hist_OneGroup_Gamma_MSCLW;
    vector<TH2D> Hist_Gamma_MSCLW;
    vector<TH2D> Hist_OnData_MSCLW;
    vector<TH2D> Hist_OnBkgd_MSCLW;
    vector<TH2D> Hist_OnBkgd_Unblind_wGamma_MSCLW;
    vector<TH2D> Hist_OnBkgd_Unblind_woGamma_MSCLW;
    vector<TH2D> Hist_OnDark_MSCLW;
    vector<TH2D> Hist_Rank0_MSCLW;
    vector<TH2D> Hist_Rank1_MSCLW;
    vector<TH2D> Hist_Rank2_MSCLW;
    vector<TH2D> Hist_Rank3_MSCLW;
    vector<TH1D> Hist_Data_Rank0_LeftVector;
    vector<TH1D> Hist_Data_Rank1_LeftVector;
    vector<TH1D> Hist_Data_Rank2_LeftVector;
    vector<TH1D> Hist_Data_Rank3_LeftVector;
    vector<TH1D> Hist_Data_Rank0_RightVector;
    vector<TH1D> Hist_Data_Rank1_RightVector;
    vector<TH1D> Hist_Data_Rank2_RightVector;
    vector<TH1D> Hist_Data_Rank3_RightVector;
    vector<TH1D> Hist_Bkgd_Rank0_LeftVector;
    vector<TH1D> Hist_Bkgd_Rank1_LeftVector;
    vector<TH1D> Hist_Bkgd_Rank2_LeftVector;
    vector<TH1D> Hist_Bkgd_Rank3_LeftVector;
    vector<TH1D> Hist_Bkgd_Rank0_RightVector;
    vector<TH1D> Hist_Bkgd_Rank1_RightVector;
    vector<TH1D> Hist_Bkgd_Rank2_RightVector;
    vector<TH1D> Hist_Bkgd_Rank3_RightVector;
    vector<TH1D> Hist_Dark_Rank0_LeftVector;
    vector<TH1D> Hist_Dark_Rank1_LeftVector;
    vector<TH1D> Hist_Dark_Rank2_LeftVector;
    vector<TH1D> Hist_Dark_Rank3_LeftVector;
    vector<TH1D> Hist_Dark_Rank0_RightVector;
    vector<TH1D> Hist_Dark_Rank1_RightVector;
    vector<TH1D> Hist_Dark_Rank2_RightVector;
    vector<TH1D> Hist_Dark_Rank3_RightVector;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));

        NumberOfEigenvectors = N_max_ranks_func_E[e];

        MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
        MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
        N_bins_for_deconv = N_bins_for_deconv_func_E[e];
        ResetMatrixDimension();

        Hist_OneGroup_Data_MSCLW.push_back(TH2D("Hist_OneGroup_Data_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OneGroup_Bkgd_MSCLW.push_back(TH2D("Hist_OneGroup_Bkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OneGroup_Bkgd_Unblind_wGamma_MSCLW.push_back(TH2D("Hist_OneGroup_Bkgd_Unblind_wGamma_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OneGroup_Bkgd_Unblind_woGamma_MSCLW.push_back(TH2D("Hist_OneGroup_Bkgd_Unblind_woGamma_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OneGroup_Gamma_MSCLW.push_back(TH2D("Hist_OneGroup_Gamma_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Gamma_MSCLW.push_back(TH2D("Hist_Gamma_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnData_MSCLW.push_back(TH2D("Hist_OnData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnBkgd_MSCLW.push_back(TH2D("Hist_OnBkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnBkgd_Unblind_wGamma_MSCLW.push_back(TH2D("Hist_OnBkgd_Unblind_wGamma_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnBkgd_Unblind_woGamma_MSCLW.push_back(TH2D("Hist_OnBkgd_Unblind_woGamma_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnDark_MSCLW.push_back(TH2D("Hist_OnDark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));

        Hist_Rank0_MSCLW.push_back(TH2D("Hist_Rank0_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank1_MSCLW.push_back(TH2D("Hist_Rank1_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank2_MSCLW.push_back(TH2D("Hist_Rank2_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank3_MSCLW.push_back(TH2D("Hist_Rank3_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));

        Hist_Data_Rank0_LeftVector.push_back(TH1D("Hist_Data_Rank0_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank1_LeftVector.push_back(TH1D("Hist_Data_Rank1_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank2_LeftVector.push_back(TH1D("Hist_Data_Rank2_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank3_LeftVector.push_back(TH1D("Hist_Data_Rank3_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank0_RightVector.push_back(TH1D("Hist_Data_Rank0_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank1_RightVector.push_back(TH1D("Hist_Data_Rank1_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank2_RightVector.push_back(TH1D("Hist_Data_Rank2_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank3_RightVector.push_back(TH1D("Hist_Data_Rank3_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank0_LeftVector.push_back(TH1D("Hist_Bkgd_Rank0_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank1_LeftVector.push_back(TH1D("Hist_Bkgd_Rank1_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank2_LeftVector.push_back(TH1D("Hist_Bkgd_Rank2_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank3_LeftVector.push_back(TH1D("Hist_Bkgd_Rank3_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank0_RightVector.push_back(TH1D("Hist_Bkgd_Rank0_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank1_RightVector.push_back(TH1D("Hist_Bkgd_Rank1_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank2_RightVector.push_back(TH1D("Hist_Bkgd_Rank2_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank3_RightVector.push_back(TH1D("Hist_Bkgd_Rank3_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Rank0_LeftVector.push_back(TH1D("Hist_Dark_Rank0_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Rank1_LeftVector.push_back(TH1D("Hist_Dark_Rank1_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Rank2_LeftVector.push_back(TH1D("Hist_Dark_Rank2_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Rank3_LeftVector.push_back(TH1D("Hist_Dark_Rank3_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Rank0_RightVector.push_back(TH1D("Hist_Dark_Rank0_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Rank1_RightVector.push_back(TH1D("Hist_Dark_Rank1_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Rank2_RightVector.push_back(TH1D("Hist_Dark_Rank2_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Rank3_RightVector.push_back(TH1D("Hist_Dark_Rank3_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
    }

    vector<TH1D> Hist_OnData_SR_Energy;
    vector<TH1D> Hist_OnData_SR_Energy_OneGroup;
    vector<TH1D> Hist_OnData_CR_Energy;
    vector<TH1D> Hist_OnData_CR_Energy_OneGroup;
    vector<TH1D> Hist_OnData_CR_Energy_Raw_OneGroup;
    vector<TH1D> Hist_OnData_SR_Zenith;
    vector<TH1D> Hist_OnData_SR_Zenith_OneGroup;
    vector<TH1D> Hist_OnData_CR_Zenith;
    vector<TH1D> Hist_OnData_CR_Zenith_OneGroup;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));

        NumberOfEigenvectors = N_max_ranks_func_E[e];

        MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
        MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
        N_bins_for_deconv = N_bins_for_deconv_func_E[e];
        ResetMatrixDimension();

        Hist_OnData_SR_Energy.push_back(TH1D("Hist_OnData_SR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_SR_Energy_OneGroup.push_back(TH1D("Hist_OnData_SR_Energy_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_CR_Energy.push_back(TH1D("Hist_OnData_CR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_CR_Energy_OneGroup.push_back(TH1D("Hist_OnData_CR_Energy_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_CR_Energy_Raw_OneGroup.push_back(TH1D("Hist_OnData_CR_Energy_Raw_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_SR_Zenith.push_back(TH1D("Hist_OnData_SR_Zenith_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
        Hist_OnData_SR_Zenith_OneGroup.push_back(TH1D("Hist_OnData_SR_Zenith_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
        Hist_OnData_CR_Zenith.push_back(TH1D("Hist_OnData_CR_Zenith_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
        Hist_OnData_CR_Zenith_OneGroup.push_back(TH1D("Hist_OnData_CR_Zenith_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
    }

    vector<TH1D> Hist_OnData_SR_Skymap_Theta2;
    vector<TH1D> Hist_OnData_SR_Skymap_Theta2_OneGroup;
    vector<TH1D> Hist_OnData_CR_Skymap_Theta2;
    vector<TH1D> Hist_OnData_CR_Skymap_Theta2_OneGroup;
    vector<TH1D> Hist_OnDark_CR_Skymap_Theta2;
    vector<TH1D> Hist_OnDark_CR_Skymap_Theta2_OneGroup;
    vector<TH2D> Hist_OnData_SR_Skymap;
    vector<TH2D> Hist_OnData_SR_Skymap_OneGroup;
    vector<TH2D> Hist_OnData_VR_Skymap;
    vector<TH2D> Hist_OnData_VR_Skymap_OneGroup;
    vector<TH2D> Hist_OnData_CR_Skymap;
    vector<TH2D> Hist_OnData_CR_Skymap_OneGroup;
    vector<TH2D> Hist_OnDark_CR_Skymap;
    vector<TH2D> Hist_OnDark_CR_Skymap_OneGroup;
    vector<TH2D> Hist_OnData_CR_Skymap_Syst;
    vector<TH2D> Hist_OnData_CR_Skymap_Syst_OneGroup;
    vector<TH2D> Hist_OnDark_CR_Skymap_Syst;
    vector<TH2D> Hist_OnDark_CR_Skymap_Syst_OneGroup;
    vector<TH2D> Hist_OnData_SR_Skymap_Galactic;
    vector<TH2D> Hist_OnData_SR_Skymap_Galactic_OneGroup;
    vector<TH2D> Hist_OnData_VR_Skymap_Galactic;
    vector<TH2D> Hist_OnData_VR_Skymap_Galactic_OneGroup;
    vector<TH2D> Hist_OnData_CR_Skymap_Galactic;
    vector<TH2D> Hist_OnData_CR_Skymap_Galactic_OneGroup;
    vector<TH2D> Hist_OnData_CR_Skymap_Galactic_Syst;
    vector<TH2D> Hist_OnData_CR_Skymap_Galactic_Syst_OneGroup;
    vector<TH2D> Hist_OnDark_CR_Skymap_Galactic_Syst;
    vector<TH2D> Hist_OnDark_CR_Skymap_Galactic_Syst_OneGroup;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));

        NumberOfEigenvectors = N_max_ranks_func_E[e];

        MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
        MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
        N_bins_for_deconv = N_bins_for_deconv_func_E[e];
        ResetMatrixDimension();

        Hist_OnData_SR_Skymap_Theta2.push_back(TH1D("Hist_OnData_SR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_SR_Skymap_Theta2_OneGroup.push_back(TH1D("Hist_OnData_SR_Skymap_Theta2_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_CR_Skymap_Theta2.push_back(TH1D("Hist_OnData_CR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_CR_Skymap_Theta2_OneGroup.push_back(TH1D("Hist_OnData_CR_Skymap_Theta2_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnDark_CR_Skymap_Theta2.push_back(TH1D("Hist_OnDark_CR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnDark_CR_Skymap_Theta2_OneGroup.push_back(TH1D("Hist_OnDark_CR_Skymap_Theta2_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_SR_Skymap.push_back(TH2D("Hist_OnData_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_SR_Skymap_OneGroup.push_back(TH2D("Hist_OnData_SR_Skymap_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_VR_Skymap.push_back(TH2D("Hist_OnData_VR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_VR_Skymap_OneGroup.push_back(TH2D("Hist_OnData_VR_Skymap_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_CR_Skymap.push_back(TH2D("Hist_OnData_CR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_CR_Skymap_OneGroup.push_back(TH2D("Hist_OnData_CR_Skymap_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnDark_CR_Skymap.push_back(TH2D("Hist_OnDark_CR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnDark_CR_Skymap_OneGroup.push_back(TH2D("Hist_OnDark_CR_Skymap_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_CR_Skymap_Syst.push_back(TH2D("Hist_OnData_Skymap_Syst_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_CR_Skymap_Syst_OneGroup.push_back(TH2D("Hist_OnData_CR_Skymap_Syst_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnDark_CR_Skymap_Syst.push_back(TH2D("Hist_OnDark_Skymap_Syst_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnDark_CR_Skymap_Syst_OneGroup.push_back(TH2D("Hist_OnDark_CR_Skymap_Syst_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,75,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));

        Hist_OnData_SR_Skymap_Galactic.push_back(TH2D("Hist_OnData_SR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,75,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_SR_Skymap_Galactic_OneGroup.push_back(TH2D("Hist_OnData_SR_Skymap_Galactic_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,75,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_VR_Skymap_Galactic.push_back(TH2D("Hist_OnData_VR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,75,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_VR_Skymap_Galactic_OneGroup.push_back(TH2D("Hist_OnData_VR_Skymap_Galactic_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,75,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_CR_Skymap_Galactic.push_back(TH2D("Hist_OnData_CR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,75,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_CR_Skymap_Galactic_OneGroup.push_back(TH2D("Hist_OnData_CR_Skymap_Galactic_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,75,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_CR_Skymap_Galactic_Syst.push_back(TH2D("Hist_OnData_Skymap_Galactic_Syst_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,75,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_CR_Skymap_Galactic_Syst_OneGroup.push_back(TH2D("Hist_OnData_CR_Skymap_Galactic_Syst_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,75,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnDark_CR_Skymap_Galactic_Syst.push_back(TH2D("Hist_OnDark_Skymap_Galactic_Syst_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,75,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnDark_CR_Skymap_Galactic_Syst_OneGroup.push_back(TH2D("Hist_OnDark_CR_Skymap_Galactic_Syst_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",75,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,75,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));

    }

    vector<vector<TH1D>> Hist_OnData_SR_RoI_Energy;
    vector<vector<TH1D>> Hist_OnData_SR_RoI_Energy_OneGroup;
    vector<vector<TH1D>> Hist_OnData_CR_RoI_Energy;
    vector<vector<TH1D>> Hist_OnData_CR_RoI_Energy_OneGroup;
    vector<vector<TH1D>> Hist_OnData_SR_Skymap_RoI_Theta2;
    vector<vector<TH1D>> Hist_OnData_SR_Skymap_RoI_Theta2_OneGroup;
    vector<vector<TH1D>> Hist_OnData_CR_Skymap_RoI_Theta2;
    vector<vector<TH1D>> Hist_OnData_CR_Skymap_RoI_Theta2_OneGroup;
    vector<vector<TH1D>> Hist_OnData_SR_RoI_MJD;
    vector<vector<TH1D>> Hist_OnData_SR_RoI_MJD_OneGroup;
    vector<vector<TH1D>> Hist_OnData_CR_RoI_MJD;
    vector<vector<TH1D>> Hist_OnData_CR_RoI_MJD_OneGroup;
    for (int nth_roi=0;nth_roi<roi_name_ptr->size();nth_roi++)
    {
        char roi_tag[50];
        sprintf(roi_tag, "%i", nth_roi);
        vector<TH1D> Hist_OnData_OneRoI_SR_RoI_Energy;
        vector<TH1D> Hist_OnData_OneRoI_SR_RoI_Energy_OneGroup;
        vector<TH1D> Hist_OnData_OneRoI_CR_RoI_Energy;
        vector<TH1D> Hist_OnData_OneRoI_CR_RoI_Energy_OneGroup;
        vector<TH1D> Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2;
        vector<TH1D> Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2_OneGroup;
        vector<TH1D> Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2;
        vector<TH1D> Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2_OneGroup;
        vector<TH1D> Hist_OnData_OneRoI_SR_RoI_MJD;
        vector<TH1D> Hist_OnData_OneRoI_SR_RoI_MJD_OneGroup;
        vector<TH1D> Hist_OnData_OneRoI_CR_RoI_MJD;
        vector<TH1D> Hist_OnData_OneRoI_CR_RoI_MJD_OneGroup;
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            NumberOfEigenvectors = N_max_ranks_func_E[e];

            MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
            MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
            N_bins_for_deconv = N_bins_for_deconv_func_E[e];
            ResetMatrixDimension();

            Hist_OnData_OneRoI_SR_RoI_Energy.push_back(TH1D("Hist_OnData_SR_RoI_Energy_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
            Hist_OnData_OneRoI_SR_RoI_Energy_OneGroup.push_back(TH1D("Hist_OnData_SR_RoI_Energy_OneGroup_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
            Hist_OnData_OneRoI_CR_RoI_Energy.push_back(TH1D("Hist_OnData_CR_RoI_Energy_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
            Hist_OnData_OneRoI_CR_RoI_Energy_OneGroup.push_back(TH1D("Hist_OnData_CR_RoI_Energy_OneGroup_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        }
        Hist_OnData_SR_RoI_Energy.push_back(Hist_OnData_OneRoI_SR_RoI_Energy);
        Hist_OnData_SR_RoI_Energy_OneGroup.push_back(Hist_OnData_OneRoI_SR_RoI_Energy_OneGroup);
        Hist_OnData_CR_RoI_Energy.push_back(Hist_OnData_OneRoI_CR_RoI_Energy);
        Hist_OnData_CR_RoI_Energy_OneGroup.push_back(Hist_OnData_OneRoI_CR_RoI_Energy_OneGroup);
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            NumberOfEigenvectors = N_max_ranks_func_E[e];

            MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
            MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
            N_bins_for_deconv = N_bins_for_deconv_func_E[e];
            ResetMatrixDimension();

            Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2.push_back(TH1D("Hist_OnData_SR_Skymap_RoI_Theta2_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,0.5));
            Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2_OneGroup.push_back(TH1D("Hist_OnData_SR_Skymap_RoI_Theta2_OneGroup_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,0.5));
            Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2.push_back(TH1D("Hist_OnData_CR_Skymap_RoI_Theta2_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,0.5));
            Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2_OneGroup.push_back(TH1D("Hist_OnData_CR_Skymap_RoI_Theta2_OneGroup_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,0.5));
            Hist_OnData_OneRoI_SR_RoI_MJD.push_back(TH1D("Hist_OnData_SR_RoI_MJD_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",800,56200-4000,56200+4000));
            Hist_OnData_OneRoI_SR_RoI_MJD_OneGroup.push_back(TH1D("Hist_OnData_SR_RoI_MJD_OneGroup_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",800,56200-4000,56200+4000));
            Hist_OnData_OneRoI_CR_RoI_MJD.push_back(TH1D("Hist_OnData_CR_RoI_MJD_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",800,56200-4000,56200+4000));
            Hist_OnData_OneRoI_CR_RoI_MJD_OneGroup.push_back(TH1D("Hist_OnData_CR_RoI_MJD_OneGroup_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",800,56200-4000,56200+4000));
        }
        Hist_OnData_SR_Skymap_RoI_Theta2.push_back(Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2);
        Hist_OnData_SR_Skymap_RoI_Theta2_OneGroup.push_back(Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2_OneGroup);
        Hist_OnData_CR_Skymap_RoI_Theta2.push_back(Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2);
        Hist_OnData_CR_Skymap_RoI_Theta2_OneGroup.push_back(Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2_OneGroup);
        Hist_OnData_SR_RoI_MJD.push_back(Hist_OnData_OneRoI_SR_RoI_MJD);
        Hist_OnData_SR_RoI_MJD_OneGroup.push_back(Hist_OnData_OneRoI_SR_RoI_MJD_OneGroup);
        Hist_OnData_CR_RoI_MJD.push_back(Hist_OnData_OneRoI_CR_RoI_MJD);
        Hist_OnData_CR_RoI_MJD_OneGroup.push_back(Hist_OnData_OneRoI_CR_RoI_MJD_OneGroup);
    }

    exposure_hours = 0.;
    MJD_Start = 2147483647;
    MJD_End = 0;
    for (int e=0;e<N_energy_bins;e++) 
    {
        std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "energy = " << energy_bins[e] << std::endl;
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));

        NumberOfEigenvectors = N_max_ranks_func_E[e];

        MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
        MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
        N_bins_for_deconv = N_bins_for_deconv_func_E[e];
        ResetMatrixDimension();

        TString hist_name;

        for (int on_run=FirstRun;on_run<LastRun;on_run++)
        {

            bool use_this_run = true;
            if (MJD_start_cut!=0 || MJD_end_cut!=0)
            {
                if (MJD_start_cut>Data_runlist_MJD_ptr->at(on_run)) use_this_run = false;
                if (MJD_end_cut<Data_runlist_MJD_ptr->at(on_run)) use_this_run = false;
                if (use_this_run) std::cout << Data_runlist_MJD_ptr->at(on_run) << std::endl;
            }

            if (use_this_run)
            {
                if (e==0)
                {
                    exposure_hours += Data_runlist_exposure_ptr->at(on_run);
                }
                if (MJD_Start>Data_runlist_MJD_ptr->at(on_run)) MJD_Start = Data_runlist_MJD_ptr->at(on_run);
                if (MJD_End<Data_runlist_MJD_ptr->at(on_run)) MJD_End = Data_runlist_MJD_ptr->at(on_run);
                char sample_tag[50];
                sprintf(sample_tag, "%i", on_run);

                int binx_lower = Hist_OneGroup_Data_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_lower);
                binx_blind_global = Hist_OneGroup_Data_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind)-1;
                int binx_upper = Hist_OneGroup_Data_MSCLW.at(e).GetXaxis()->FindBin(1.)-1;
                int biny_lower = Hist_OneGroup_Data_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_lower);
                biny_blind_global = Hist_OneGroup_Data_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind)-1;
                int biny_upper = Hist_OneGroup_Data_MSCLW.at(e).GetYaxis()->FindBin(1.)-1;
                hist_name  = "Hist_OnData_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OneGroup_Data_MSCLW.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                group_size.at(e) += ((TH2D*)InputDataFile.Get(hist_name))->Integral();

                //hist_name  = "Hist_GammaData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                hist_name  = "Hist_GammaMC_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OneGroup_Gamma_MSCLW.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                mtx_gamma_raw = fillMatrix(&Hist_OneGroup_Gamma_MSCLW.at(e));
                hist_name  = "Hist_OnData_SR_Skymap_Theta2_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_SR_Skymap_Theta2_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_CR_Skymap_Theta2_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_CR_Skymap_Theta2_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                Hist_OnDark_CR_Skymap_Theta2_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_SR_Skymap_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_SR_Skymap_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_VR_Skymap_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_VR_Skymap_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_CR_Skymap_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_CR_Skymap_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                Hist_OnDark_CR_Skymap_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnDark_CR_Skymap_Syst_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnDark_CR_Skymap_Syst_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_SR_Skymap_Galactic_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_SR_Skymap_Galactic_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_VR_Skymap_Galactic_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_VR_Skymap_Galactic_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_CR_Skymap_Galactic_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_CR_Skymap_Galactic_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnDark_CR_Skymap_Galactic_Syst_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnDark_CR_Skymap_Galactic_Syst_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_SR_Energy_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_SR_Energy_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_CR_Energy_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_CR_Energy_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_CR_Energy_Raw_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_CR_Energy_Raw_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_SR_Zenith_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_SR_Zenith_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                hist_name  = "Hist_OnData_CR_Zenith_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OnData_CR_Zenith_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                for (int nth_roi=0;nth_roi<roi_name_ptr->size();nth_roi++)
                {
                    char roi_tag[50];
                    sprintf(roi_tag, "%i", nth_roi);
                    hist_name  = "Hist_OnData_SR_RoI_Energy_R"+TString(sample_tag)+"_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                    Hist_OnData_SR_RoI_Energy_OneGroup.at(nth_roi).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                    hist_name  = "Hist_OnData_SR_Skymap_RoI_Theta2_R"+TString(sample_tag)+"_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                    Hist_OnData_SR_Skymap_RoI_Theta2_OneGroup.at(nth_roi).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                    hist_name  = "Hist_OnData_SR_RoI_MJD_R"+TString(sample_tag)+"_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                    Hist_OnData_SR_RoI_MJD_OneGroup.at(nth_roi).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                    hist_name  = "Hist_OnData_CR_RoI_Energy_R"+TString(sample_tag)+"_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                    Hist_OnData_CR_RoI_Energy_OneGroup.at(nth_roi).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                    hist_name  = "Hist_OnData_CR_Skymap_RoI_Theta2_R"+TString(sample_tag)+"_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                    Hist_OnData_CR_Skymap_RoI_Theta2_OneGroup.at(nth_roi).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                    hist_name  = "Hist_OnData_CR_RoI_MJD_R"+TString(sample_tag)+"_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                    Hist_OnData_CR_RoI_MJD_OneGroup.at(nth_roi).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                }
                for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
                {
                    char sample2_tag[50];
                    sprintf(sample2_tag, "%i", nth_sample);
                    hist_name  = "Hist_OnDark_MSCLW_R"+TString(sample_tag)+"_V"+TString(sample2_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                    Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                }

            }

            if (group_size.at(e)>=group_size_limit[e] || on_run==Data_runlist_number_ptr->size()-1)
            {

                TH2D Hist_Temp_Gamma = TH2D("Hist_Temp_Gamma","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
                int binx_lower = Hist_Temp_Gamma.GetXaxis()->FindBin(MSCL_cut_lower);
                int binx_blind = Hist_Temp_Gamma.GetXaxis()->FindBin(MSCL_cut_blind)-1;
                int binx_upper = Hist_Temp_Gamma.GetXaxis()->FindBin(1.)-1;
                int biny_lower = Hist_Temp_Gamma.GetYaxis()->FindBin(MSCW_cut_lower);
                int biny_blind = Hist_Temp_Gamma.GetYaxis()->FindBin(MSCW_cut_blind)-1;
                int biny_upper = Hist_Temp_Gamma.GetYaxis()->FindBin(1.)-1;
                fill2DHistogram(&Hist_Temp_Gamma,mtx_gamma_raw);
                double Gamma_SR_Integral = Hist_Temp_Gamma.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
                double Data_SR_Integral = Hist_OneGroup_Data_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
                Hist_Temp_Gamma.Scale((Data_SR_Integral*SignalNorm)/Gamma_SR_Integral);
                Hist_OneGroup_Data_MSCLW.at(e).Add(&Hist_Temp_Gamma);

                mtx_data = fillMatrix(&Hist_OneGroup_Data_MSCLW.at(e));
                eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);

                for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
                {

                    NormalizaDarkMatrix(&Hist_OneGroup_Data_MSCLW.at(e), &Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e));
                    TH2D hist_dark_alter = TH2D("hist_dark_alter","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
                    AlterDarkMatrix(&Hist_OneGroup_Data_MSCLW.at(e), &Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e), &hist_dark_alter);

                    fill2DHistogram(&Hist_Temp_Gamma,mtx_gamma_raw);
                    Gamma_SR_Integral = Hist_Temp_Gamma.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
                    Hist_Temp_Gamma.Scale((0.0)/Gamma_SR_Integral);
                    mtx_gamma = fillMatrix(&Hist_Temp_Gamma);

                    mtx_dark = fillMatrix(&hist_dark_alter);
                    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
                    LeastSquareSolutionMethod(false, true);
                    TH2D Hist_Temp_Bkgd = TH2D("Hist_Temp_Bkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
                    fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
                    
                    if (UseMinChi2)
                    {
                        mtx_data_bkgd = MinChi2Method(&Hist_OneGroup_OnSyst_Chi2.at(nth_sample).at(e), &Hist_OneGroup_Data_MSCLW.at(e), &Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e));
                        fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
                    }

                    double Bkgd_SR_Integral = Hist_Temp_Bkgd.Integral(binx_lower,binx_blind,biny_lower,biny_blind);

                    //fill2DHistogram(&Hist_Temp_Gamma,mtx_gamma_raw);
                    //Gamma_SR_Integral = Hist_Temp_Gamma.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
                    //Hist_Temp_Gamma.Scale(max(0.,Data_SR_Integral-Bkgd_SR_Integral)/Gamma_SR_Integral);
                    //mtx_gamma = fillMatrix(&Hist_Temp_Gamma);
                    //LeastSquareSolutionMethod(false, true);
                    //fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
                    
                    double dark_weight = 1./double(n_dark_samples);
                    int error_bin = Hist_OneGroup_OnSyst_Chi2.at(nth_sample).at(e).GetMinimumBin();
                    GetNoiseReplacedMatrix(&Hist_OneGroup_Data_MSCLW.at(e),&Hist_Temp_Bkgd,false);
                    Hist_OneGroup_Bkgd_MSCLW.at(e).Add(&Hist_Temp_Bkgd,dark_weight);
                    Hist_OnSyst_MSCLW.at(nth_sample).at(e).Add(&Hist_Temp_Bkgd);
                    GetNoiseReplacedMatrix(&Hist_OneGroup_Data_MSCLW.at(e),&hist_dark_alter,false);
                    Hist_OnDark_MSCLW.at(e).Add(&hist_dark_alter,dark_weight);
                    Hist_OnSyst_Chi2.at(nth_sample).at(e).Add(&Hist_OneGroup_OnSyst_Chi2.at(nth_sample).at(e));
                    Hist_OneGroup_OnSyst_Chi2.at(nth_sample).at(e).Reset();
                    
                    //fill2DHistogram(&Hist_Temp_Gamma,mtx_gamma_raw);
                    //Gamma_SR_Integral = Hist_Temp_Gamma.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
                    //Hist_Temp_Gamma.Scale((Bkgd_SR_Integral*0.0)/Gamma_SR_Integral);
                    //mtx_gamma = fillMatrix(&Hist_Temp_Gamma);
                    //LeastSquareSolutionMethod(false, false);
                    //fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
                    //Hist_OneGroup_Bkgd_Unblind_woGamma_MSCLW.at(e).Add(&Hist_Temp_Bkgd,dark_weight);
                    //fill2DHistogram(&Hist_Temp_Gamma,mtx_gamma_raw);
                    //Gamma_SR_Integral = Hist_Temp_Gamma.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
                    //Hist_Temp_Gamma.Scale((Bkgd_SR_Integral*0.1)/Gamma_SR_Integral);
                    //mtx_gamma = fillMatrix(&Hist_Temp_Gamma);
                    //LeastSquareSolutionMethod(false, false);
                    //fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
                    //Hist_OneGroup_Bkgd_Unblind_wGamma_MSCLW.at(e).Add(&Hist_Temp_Bkgd,dark_weight);
                    //Hist_OneGroup_Gamma_MSCLW.at(e).Reset();
                    //Hist_OneGroup_Gamma_MSCLW.at(e).Add(&Hist_Temp_Gamma);
                }
                if (!isnan(Hist_OneGroup_Bkgd_MSCLW.at(e).Integral()) && !isnan(Hist_OneGroup_Data_MSCLW.at(e).Integral()))
                {
                    Hist_OnData_MSCLW.at(e).Add(&Hist_OneGroup_Data_MSCLW.at(e));
                    Hist_OnBkgd_MSCLW.at(e).Add(&Hist_OneGroup_Bkgd_MSCLW.at(e));
                    Hist_OnBkgd_Unblind_wGamma_MSCLW.at(e).Add(&Hist_OneGroup_Bkgd_Unblind_wGamma_MSCLW.at(e));
                    Hist_OnBkgd_Unblind_woGamma_MSCLW.at(e).Add(&Hist_OneGroup_Bkgd_Unblind_woGamma_MSCLW.at(e));
                    Hist_Gamma_MSCLW.at(e).Add(&Hist_OneGroup_Gamma_MSCLW.at(e));

                    int binx_lower = Hist_OnData_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_lower);
                    int binx_blind = Hist_OnData_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind)-1;
                    int binx_validate = Hist_OnData_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind+0.2)-1;
                    int binx_upper = Hist_OnData_MSCLW.at(e).GetXaxis()->FindBin(1.)-1;
                    int biny_lower = Hist_OnData_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_lower);
                    int biny_blind = Hist_OnData_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind)-1;
                    int biny_validate = Hist_OnData_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind+0.2)-1;
                    int biny_upper = Hist_OnData_MSCLW.at(e).GetYaxis()->FindBin(1.)-1;
                    double Bkgd_SR_Integral = Hist_OneGroup_Bkgd_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
                    double Dark_SR_Integral = 0.;
                    double dark_weight = 1./double(n_dark_samples);
                    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
                    {
                        Dark_SR_Integral += dark_weight*Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
                    }
                    double Bkgd_VR_Integral = Hist_OnData_VR_Skymap_OneGroup.at(e).Integral();
                    double Old_Integral = Hist_OnData_CR_Energy_OneGroup.at(e).Integral();
                    double Old_Integral_Raw = Hist_OnData_CR_Energy_Raw_OneGroup.at(e).Integral();
                    double Old_Integral_VR = Hist_OnDark_CR_Skymap_Syst_OneGroup.at(e).Integral();
                    double scale = 0.;
                    double scale_dark = 0.;
                    double scale_raw = 0.;
                    double scale_vr = 0.;
                    double scale_vrsr = 0.;
                    if (Old_Integral>0.)
                    {
                        scale = Bkgd_SR_Integral/Old_Integral;
                        scale_dark = Dark_SR_Integral/Old_Integral;
                        scale_raw = Bkgd_SR_Integral/Old_Integral_Raw;
                    }
                    if (Bkgd_VR_Integral>0.)
                    {
                        scale_vrsr = Bkgd_SR_Integral/Bkgd_VR_Integral;
                    }
                    if (Old_Integral_VR>0.)
                    {
                        scale_vr = Bkgd_VR_Integral/Old_Integral_VR;
                    }
                    std::cout << "Old_Integral_VR = " << Old_Integral_VR << std::endl;
                    std::cout << "Bkgd_VR_Integral = " << Bkgd_VR_Integral << std::endl;
                    std::cout << "scale_vr = " << scale_vr << std::endl;
                    std::cout << "scale_vrsr = " << scale_vrsr << std::endl;
                    Hist_OnData_CR_Energy_OneGroup.at(e).Scale(scale);
                    Hist_OnData_CR_Zenith_OneGroup.at(e).Scale(scale);
                    Hist_OnData_CR_Skymap_Theta2_OneGroup.at(e).Scale(scale);
                    Hist_OnDark_CR_Skymap_Theta2_OneGroup.at(e).Scale(scale_dark);
                    Hist_OnData_CR_Skymap_OneGroup.at(e).Scale(scale);
                    Hist_OnDark_CR_Skymap_OneGroup.at(e).Scale(scale_dark);
                    Hist_OnData_CR_Skymap_Galactic_OneGroup.at(e).Scale(scale);
                    //RingBackgroundSkymap(&Hist_OnData_VR_Skymap_OneGroup.at(e),&Hist_OnDark_CR_Skymap_Syst_OneGroup.at(e));
                    //Hist_OnDark_CR_Skymap_Syst_OneGroup.at(e).Scale(scale_vr);
                    //Hist_OnDark_CR_Skymap_Syst_OneGroup.at(e).Add(&Hist_OnData_VR_Skymap_OneGroup.at(e),-1.);
                    //Hist_OnDark_CR_Skymap_Syst_OneGroup.at(e).Scale(scale_vrsr);
                    //RingBackgroundSkymap(&Hist_OnData_VR_Skymap_Galactic_OneGroup.at(e),&Hist_OnDark_CR_Skymap_Galactic_Syst_OneGroup.at(e));
                    //Hist_OnDark_CR_Skymap_Galactic_Syst_OneGroup.at(e).Scale(scale_vr);
                    //Hist_OnDark_CR_Skymap_Galactic_Syst_OneGroup.at(e).Add(&Hist_OnData_VR_Skymap_Galactic_OneGroup.at(e),-1.);
                    //Hist_OnData_CR_Skymap_Galactic_Syst_OneGroup.at(e).Scale(scale_vrsr);
                    for (int nth_roi=0;nth_roi<roi_name_ptr->size();nth_roi++)
                    {
                        Hist_OnData_CR_RoI_Energy_OneGroup.at(nth_roi).at(e).Scale(scale);
                        Hist_OnData_CR_Skymap_RoI_Theta2_OneGroup.at(nth_roi).at(e).Scale(scale);
                        Hist_OnData_CR_RoI_MJD_OneGroup.at(nth_roi).at(e).Scale(scale_raw);
                    }
                    Hist_OnData_SR_Energy.at(e).Add(&Hist_OnData_SR_Energy_OneGroup.at(e));
                    Hist_OnData_CR_Energy.at(e).Add(&Hist_OnData_CR_Energy_OneGroup.at(e));
                    Hist_OnData_SR_Zenith.at(e).Add(&Hist_OnData_SR_Zenith_OneGroup.at(e));
                    Hist_OnData_CR_Zenith.at(e).Add(&Hist_OnData_CR_Zenith_OneGroup.at(e));
                    Hist_OnData_SR_Skymap_Theta2.at(e).Add(&Hist_OnData_SR_Skymap_Theta2_OneGroup.at(e));
                    Hist_OnData_CR_Skymap_Theta2.at(e).Add(&Hist_OnData_CR_Skymap_Theta2_OneGroup.at(e));
                    Hist_OnDark_CR_Skymap_Theta2.at(e).Add(&Hist_OnDark_CR_Skymap_Theta2_OneGroup.at(e));
                    Hist_OnData_SR_Skymap.at(e).Add(&Hist_OnData_SR_Skymap_OneGroup.at(e));
                    Hist_OnData_VR_Skymap.at(e).Add(&Hist_OnData_VR_Skymap_OneGroup.at(e));
                    Hist_OnData_CR_Skymap.at(e).Add(&Hist_OnData_CR_Skymap_OneGroup.at(e));
                    Hist_OnDark_CR_Skymap.at(e).Add(&Hist_OnDark_CR_Skymap_OneGroup.at(e));
                    Hist_OnDark_CR_Skymap_Syst.at(e).Add(&Hist_OnDark_CR_Skymap_Syst_OneGroup.at(e));
                    Hist_OnData_SR_Skymap_Galactic.at(e).Add(&Hist_OnData_SR_Skymap_Galactic_OneGroup.at(e));
                    Hist_OnData_VR_Skymap_Galactic.at(e).Add(&Hist_OnData_VR_Skymap_Galactic_OneGroup.at(e));
                    Hist_OnData_CR_Skymap_Galactic.at(e).Add(&Hist_OnData_CR_Skymap_Galactic_OneGroup.at(e));
                    Hist_OnData_CR_Skymap_Galactic_Syst.at(e).Add(&Hist_OnData_CR_Skymap_Galactic_Syst_OneGroup.at(e));
                    Hist_OnDark_CR_Skymap_Galactic_Syst.at(e).Add(&Hist_OnDark_CR_Skymap_Galactic_Syst_OneGroup.at(e));
                    for (int nth_roi=0;nth_roi<roi_name_ptr->size();nth_roi++)
                    {
                        Hist_OnData_SR_RoI_Energy.at(nth_roi).at(e).Add(&Hist_OnData_SR_RoI_Energy_OneGroup.at(nth_roi).at(e));
                        Hist_OnData_CR_RoI_Energy.at(nth_roi).at(e).Add(&Hist_OnData_CR_RoI_Energy_OneGroup.at(nth_roi).at(e));
                        Hist_OnData_SR_Skymap_RoI_Theta2.at(nth_roi).at(e).Add(&Hist_OnData_SR_Skymap_RoI_Theta2_OneGroup.at(nth_roi).at(e));
                        Hist_OnData_CR_Skymap_RoI_Theta2.at(nth_roi).at(e).Add(&Hist_OnData_CR_Skymap_RoI_Theta2_OneGroup.at(nth_roi).at(e));
                        Hist_OnData_SR_RoI_MJD.at(nth_roi).at(e).Add(&Hist_OnData_SR_RoI_MJD_OneGroup.at(nth_roi).at(e));
                        Hist_OnData_CR_RoI_MJD.at(nth_roi).at(e).Add(&Hist_OnData_CR_RoI_MJD_OneGroup.at(nth_roi).at(e));
                    }
                }
                for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
                {
                    Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e).Reset();
                }
                Hist_OneGroup_Data_MSCLW.at(e).Reset();
                Hist_OneGroup_Bkgd_MSCLW.at(e).Reset();
                Hist_OneGroup_Bkgd_Unblind_wGamma_MSCLW.at(e).Reset();
                Hist_OneGroup_Bkgd_Unblind_woGamma_MSCLW.at(e).Reset();
                Hist_OneGroup_Gamma_MSCLW.at(e).Reset();
                Hist_OnData_SR_Energy_OneGroup.at(e).Reset();
                Hist_OnData_CR_Energy_OneGroup.at(e).Reset();
                Hist_OnData_CR_Energy_Raw_OneGroup.at(e).Reset();
                Hist_OnData_SR_Zenith_OneGroup.at(e).Reset();
                Hist_OnData_CR_Zenith_OneGroup.at(e).Reset();
                Hist_OnData_SR_Skymap_Theta2_OneGroup.at(e).Reset();
                Hist_OnData_CR_Skymap_Theta2_OneGroup.at(e).Reset();
                Hist_OnDark_CR_Skymap_Theta2_OneGroup.at(e).Reset();
                Hist_OnData_SR_Skymap_OneGroup.at(e).Reset();
                Hist_OnData_VR_Skymap_OneGroup.at(e).Reset();
                Hist_OnData_CR_Skymap_OneGroup.at(e).Reset();
                Hist_OnDark_CR_Skymap_OneGroup.at(e).Reset();
                Hist_OnDark_CR_Skymap_Syst_OneGroup.at(e).Reset();
                Hist_OnData_SR_Skymap_Galactic_OneGroup.at(e).Reset();
                Hist_OnData_VR_Skymap_Galactic_OneGroup.at(e).Reset();
                Hist_OnData_CR_Skymap_Galactic_OneGroup.at(e).Reset();
                Hist_OnDark_CR_Skymap_Galactic_Syst_OneGroup.at(e).Reset();
                for (int nth_roi=0;nth_roi<roi_name_ptr->size();nth_roi++)
                {
                    Hist_OnData_SR_RoI_Energy_OneGroup.at(nth_roi).at(e).Reset();
                    Hist_OnData_CR_RoI_Energy_OneGroup.at(nth_roi).at(e).Reset();
                    Hist_OnData_SR_Skymap_RoI_Theta2_OneGroup.at(nth_roi).at(e).Reset();
                    Hist_OnData_CR_Skymap_RoI_Theta2_OneGroup.at(nth_roi).at(e).Reset();
                    Hist_OnData_SR_RoI_MJD_OneGroup.at(nth_roi).at(e).Reset();
                    Hist_OnData_CR_RoI_MJD_OneGroup.at(nth_roi).at(e).Reset();
                }
                mtx_data = fillMatrix(&Hist_OnData_MSCLW.at(e));
                eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
                mtx_dark = fillMatrix(&Hist_OnDark_MSCLW.at(e));
                eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
                mtx_data_bkgd = fillMatrix(&Hist_OnBkgd_MSCLW.at(e));
                eigensolver_bkgd = ComplexEigenSolver<MatrixXcd>(mtx_data_bkgd);

                std::cout << "=====================================================================" << std::endl;
                std::cout << "eigensolver_data.eigenvalues()(mtx_data.cols()-1):" << std::endl;
                std::cout << eigensolver_data.eigenvalues()(mtx_data.cols()-1) << std::endl;
                std::cout << "eigensolver_data.eigenvalues()(mtx_data.cols()-2):" << std::endl;
                std::cout << eigensolver_data.eigenvalues()(mtx_data.cols()-2) << std::endl;
                std::cout << "eigensolver_data.eigenvalues()(mtx_data.cols()-3):" << std::endl;
                std::cout << eigensolver_data.eigenvalues()(mtx_data.cols()-3) << std::endl;
                std::cout << "eigensolver_data.eigenvalues()(mtx_data.cols()-4):" << std::endl;
                std::cout << eigensolver_data.eigenvalues()(mtx_data.cols()-4) << std::endl;
                std::cout << "=====================================================================" << std::endl;
                std::cout << "eigensolver_dark.eigenvalues()(mtx_dark.cols()-1):" << std::endl;
                std::cout << eigensolver_dark.eigenvalues()(mtx_dark.cols()-1) << std::endl;
                std::cout << "eigensolver_dark.eigenvalues()(mtx_dark.cols()-2):" << std::endl;
                std::cout << eigensolver_dark.eigenvalues()(mtx_dark.cols()-2) << std::endl;
                std::cout << "eigensolver_dark.eigenvalues()(mtx_dark.cols()-3):" << std::endl;
                std::cout << eigensolver_dark.eigenvalues()(mtx_dark.cols()-3) << std::endl;
                std::cout << "eigensolver_dark.eigenvalues()(mtx_dark.cols()-4):" << std::endl;
                std::cout << eigensolver_dark.eigenvalues()(mtx_dark.cols()-4) << std::endl;
                std::cout << "=====================================================================" << std::endl;
                std::cout << "eigensolver_bkgd.eigenvalues()(mtx_data.cols()-1):" << std::endl;
                std::cout << eigensolver_bkgd.eigenvalues()(mtx_data.cols()-1) << std::endl;
                std::cout << "eigensolver_bkgd.eigenvalues()(mtx_data.cols()-2):" << std::endl;
                std::cout << eigensolver_bkgd.eigenvalues()(mtx_data.cols()-2) << std::endl;
                std::cout << "eigensolver_bkgd.eigenvalues()(mtx_data.cols()-3):" << std::endl;
                std::cout << eigensolver_bkgd.eigenvalues()(mtx_data.cols()-3) << std::endl;
                std::cout << "eigensolver_bkgd.eigenvalues()(mtx_data.cols()-4):" << std::endl;
                std::cout << eigensolver_bkgd.eigenvalues()(mtx_data.cols()-4) << std::endl;
                std::cout << "=====================================================================" << std::endl;

                GetReducedEigenvalueMatrix(0);
                mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_redu*eigensolver_data.eigenvectors().inverse();
                fill2DHistogram(&Hist_Rank0_MSCLW.at(e),mtx_data_redu);
                GetReducedEigenvalueMatrix(1);
                mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_redu*eigensolver_data.eigenvectors().inverse();
                fill2DHistogram(&Hist_Rank1_MSCLW.at(e),mtx_data_redu);
                GetReducedEigenvalueMatrix(2);
                mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_redu*eigensolver_data.eigenvectors().inverse();
                fill2DHistogram(&Hist_Rank2_MSCLW.at(e),mtx_data_redu);
                GetReducedEigenvalueMatrix(3);
                mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_redu*eigensolver_data.eigenvectors().inverse();
                fill2DHistogram(&Hist_Rank3_MSCLW.at(e),mtx_data_redu);


                SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data);
                fill1DHistogram(&Hist_Data_Rank0_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),0,1.0);
                fill1DHistogram(&Hist_Data_Rank1_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),1,1.0);
                fill1DHistogram(&Hist_Data_Rank2_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),2,1.0);
                fill1DHistogram(&Hist_Data_Rank3_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),2,1.0);
                fill1DHistogram(&Hist_Data_Rank0_RightVector.at(e),mtx_eigenvector_init,0,1.0);
                fill1DHistogram(&Hist_Data_Rank1_RightVector.at(e),mtx_eigenvector_init,1,1.0);
                fill1DHistogram(&Hist_Data_Rank2_RightVector.at(e),mtx_eigenvector_init,2,1.0);
                fill1DHistogram(&Hist_Data_Rank3_RightVector.at(e),mtx_eigenvector_init,3,1.0);

                SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
                fill1DHistogram(&Hist_Bkgd_Rank0_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),0,1.0);
                fill1DHistogram(&Hist_Bkgd_Rank1_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),1,1.0);
                fill1DHistogram(&Hist_Bkgd_Rank2_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),2,1.0);
                fill1DHistogram(&Hist_Bkgd_Rank3_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),3,1.0);
                fill1DHistogram(&Hist_Bkgd_Rank0_RightVector.at(e),mtx_eigenvector_init,0,1.0);
                fill1DHistogram(&Hist_Bkgd_Rank1_RightVector.at(e),mtx_eigenvector_init,1,1.0);
                fill1DHistogram(&Hist_Bkgd_Rank2_RightVector.at(e),mtx_eigenvector_init,2,1.0);
                fill1DHistogram(&Hist_Bkgd_Rank3_RightVector.at(e),mtx_eigenvector_init,3,1.0);

                SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_dark);
                fill1DHistogram(&Hist_Dark_Rank0_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),0,1.0);
                fill1DHistogram(&Hist_Dark_Rank1_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),1,1.0);
                fill1DHistogram(&Hist_Dark_Rank2_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),2,1.0);
                fill1DHistogram(&Hist_Dark_Rank3_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),3,1.0);
                fill1DHistogram(&Hist_Dark_Rank0_RightVector.at(e),mtx_eigenvector_init,0,1.0);
                fill1DHistogram(&Hist_Dark_Rank1_RightVector.at(e),mtx_eigenvector_init,1,1.0);
                fill1DHistogram(&Hist_Dark_Rank2_RightVector.at(e),mtx_eigenvector_init,2,1.0);
                fill1DHistogram(&Hist_Dark_Rank3_RightVector.at(e),mtx_eigenvector_init,3,1.0);
            
                group_size.at(e) = 0;

            }
        }
    }
    InputDataFile.Close();

    for (int e=0;e<N_energy_bins;e++) 
    {
        int binx_lower = Hist_OnData_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_lower);
        int binx_blind = Hist_OnData_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind)-1;
        int biny_lower = Hist_OnData_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_lower);
        int biny_blind = Hist_OnData_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind)-1;
        double Bkgd_SR_Integral = Hist_OnBkgd_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
        double Bkgd_VR_Integral = Hist_OnData_VR_Skymap.at(e).Integral();
        double scale_vrsr = 0.;
        if (Bkgd_VR_Integral>0.)
        {
            scale_vrsr = Bkgd_SR_Integral/Bkgd_VR_Integral;
        }
        RingBackgroundSkymap(&Hist_OnData_VR_Skymap.at(e),&Hist_OnData_CR_Skymap_Syst.at(e),&Hist_OnDark_CR_Skymap_Syst.at(e));
        Hist_OnData_CR_Skymap_Syst.at(e).Add(&Hist_OnData_VR_Skymap.at(e),-1.);
        Hist_OnData_CR_Skymap_Syst.at(e).Scale(scale_vrsr);
        RingBackgroundSkymap(&Hist_OnData_VR_Skymap_Galactic.at(e),&Hist_OnData_CR_Skymap_Galactic_Syst.at(e),&Hist_OnDark_CR_Skymap_Galactic_Syst.at(e));
        Hist_OnData_CR_Skymap_Galactic_Syst.at(e).Add(&Hist_OnData_VR_Skymap_Galactic.at(e),-1.);
        Hist_OnData_CR_Skymap_Galactic_Syst.at(e).Scale(scale_vrsr);
    }


    TFile InputFile("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+"_"+ONOFF_tag+".root");
    TTree* InfoTree = nullptr;
    InfoTree = (TTree*) InputFile.Get("InfoTree");
    TTree* StarTree = nullptr;
    StarTree = (TTree*) InputFile.Get("StarTree");
    TTree* FaintStarTree = nullptr;
    FaintStarTree = (TTree*) InputFile.Get("FaintStarTree");
    TString hist_name;
    hist_name  = "Hist_EffArea";
    Hist_EffArea.Add( (TH2D*)InputFile.Get(hist_name) );

    TFile OutputFile("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+"_"+TString(output_file2_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+TString(signal_tag)+TString(mjd_cut_tag)+"_"+ONOFF_tag+".root","recreate");

    TTree *newtree = InfoTree->CloneTree();
    newtree->Write();
    TTree *newstartree = StarTree->CloneTree();
    newstartree->Write();
    TTree *newfaintstartree = FaintStarTree->CloneTree();
    newfaintstartree->Write();

    TTree NewInfoTree("NewInfoTree","new info tree");
    NewInfoTree.Branch("MJD_Start",&MJD_Start,"MJD_Start/I");
    NewInfoTree.Branch("MJD_End",&MJD_End,"MJD_End/I");
    NewInfoTree.Branch("exposure_hours",&exposure_hours,"exposure_hours/D");
    NewInfoTree.Fill();
    NewInfoTree.Write();

    Hist_EffArea.Write();
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_OnData_MSCLW.at(e).Write();
        Hist_OnBkgd_MSCLW.at(e).Write();
        Hist_OnBkgd_Unblind_wGamma_MSCLW.at(e).Write();
        Hist_OnBkgd_Unblind_woGamma_MSCLW.at(e).Write();
        Hist_Gamma_MSCLW.at(e).Write();
        Hist_OnDark_MSCLW.at(e).Write();
        Hist_Rank0_MSCLW.at(e).Write();
        Hist_Rank1_MSCLW.at(e).Write();
        Hist_Rank2_MSCLW.at(e).Write();
        Hist_Rank3_MSCLW.at(e).Write();
        Hist_Data_Rank0_LeftVector.at(e).Write();
        Hist_Data_Rank1_LeftVector.at(e).Write();
        Hist_Data_Rank2_LeftVector.at(e).Write();
        Hist_Data_Rank3_LeftVector.at(e).Write();
        Hist_Data_Rank0_RightVector.at(e).Write();
        Hist_Data_Rank1_RightVector.at(e).Write();
        Hist_Data_Rank2_RightVector.at(e).Write();
        Hist_Data_Rank3_RightVector.at(e).Write();
        Hist_Bkgd_Rank0_LeftVector.at(e).Write();
        Hist_Bkgd_Rank1_LeftVector.at(e).Write();
        Hist_Bkgd_Rank2_LeftVector.at(e).Write();
        Hist_Bkgd_Rank3_LeftVector.at(e).Write();
        Hist_Bkgd_Rank0_RightVector.at(e).Write();
        Hist_Bkgd_Rank1_RightVector.at(e).Write();
        Hist_Bkgd_Rank2_RightVector.at(e).Write();
        Hist_Bkgd_Rank3_RightVector.at(e).Write();
        Hist_Dark_Rank0_LeftVector.at(e).Write();
        Hist_Dark_Rank1_LeftVector.at(e).Write();
        Hist_Dark_Rank2_LeftVector.at(e).Write();
        Hist_Dark_Rank3_LeftVector.at(e).Write();
        Hist_Dark_Rank0_RightVector.at(e).Write();
        Hist_Dark_Rank1_RightVector.at(e).Write();
        Hist_Dark_Rank2_RightVector.at(e).Write();
        Hist_Dark_Rank3_RightVector.at(e).Write();
    }
    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
    {
        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OffData_MSCLW.at(nth_sample).at(e).Write();
            Hist_OffBkgd_MSCLW.at(nth_sample).at(e).Write();
            Hist_OnSyst_MSCLW.at(nth_sample).at(e).Write();
            Hist_OnSyst_Chi2.at(nth_sample).at(e).Write();
        }
    }
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_OnData_SR_Energy.at(e).Write();
        Hist_OnData_CR_Energy.at(e).Write();
        Hist_OnData_SR_Zenith.at(e).Write();
        Hist_OnData_CR_Zenith.at(e).Write();
    }
    for (int nth_roi=0;nth_roi<Hist_OnData_SR_RoI_Energy.size();nth_roi++)
    {
        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OnData_SR_RoI_Energy.at(nth_roi).at(e).Write();
            Hist_OnData_CR_RoI_Energy.at(nth_roi).at(e).Write();
        }
        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OnData_SR_RoI_MJD.at(nth_roi).at(e).Write();
            Hist_OnData_CR_RoI_MJD.at(nth_roi).at(e).Write();
            Hist_OnData_SR_Skymap_RoI_Theta2.at(nth_roi).at(e).Write();
            Hist_OnData_CR_Skymap_RoI_Theta2.at(nth_roi).at(e).Write();
        }
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        Hist_OnData_SR_Skymap_Theta2.at(e).Write();
        Hist_OnData_CR_Skymap_Theta2.at(e).Write();
        Hist_OnDark_CR_Skymap_Theta2.at(e).Write();
        Hist_OnData_SR_Skymap.at(e).Write();
        Hist_OnData_VR_Skymap.at(e).Write();
        Hist_OnData_CR_Skymap.at(e).Write();
        Hist_OnDark_CR_Skymap.at(e).Write();
        Hist_OnData_CR_Skymap_Syst.at(e).Write();
        Hist_OnData_SR_Skymap_Galactic.at(e).Write();
        Hist_OnData_VR_Skymap_Galactic.at(e).Write();
        Hist_OnData_CR_Skymap_Galactic.at(e).Write();
        Hist_OnData_CR_Skymap_Galactic_Syst.at(e).Write();
    }
    OutputFile.Close();
    InputFile.Close();

    std::cout << "Done." << std::endl;

}
