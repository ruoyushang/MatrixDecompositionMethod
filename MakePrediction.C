
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

MatrixXcd mtx_regularization;
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
MatrixXcd mtx_eigenval_dark_redu;
MatrixXcd mtx_dark_redu;
MatrixXcd mtx_eigenval_diff_redu;
MatrixXcd mtx_diff_redu;
ComplexEigenSolver<MatrixXcd> eigensolver_gamma;
ComplexEigenSolver<MatrixXcd> eigensolver_gamma_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_bkgd;
ComplexEigenSolver<MatrixXcd> eigensolver_bkgd_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_dark;
ComplexEigenSolver<MatrixXcd> eigensolver_dark_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_data;
ComplexEigenSolver<MatrixXcd> eigensolver_data_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_diff;
ComplexEigenSolver<MatrixXcd> eigensolver_diff_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_init;
ComplexEigenSolver<MatrixXcd> eigensolver_init_transpose;
VectorXcd vtr_eigenval_vvv;
int binx_blind_global;
int biny_blind_global;
int VaryLeftOrRightVector;
int VaryNthVector;
int n_iterations = 100;

double GammaScale = 0.;
double MinChi2Unblind = 1e10;

//double svd_threshold = 1e-6; // size of singular value to be considered as nonzero.
double svd_threshold = 1e-20; // size of singular value to be considered as nonzero.
double svd_threshold_scale = 1.0;

vector<double> data_gamma_count;
vector<double> dark_gamma_count;
vector<double> bkgd_gamma_count;

void ResetMatrixDimension()
{
    mtx_regularization = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
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
    mtx_eigenval_dark_redu = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_dark_redu = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenval_diff_redu = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    mtx_diff_redu = MatrixXcd(N_bins_for_deconv,N_bins_for_deconv);
    vtr_eigenval_vvv = VectorXcd(N_bins_for_deconv*N_bins_for_deconv);
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
            //if (binx<hist->GetNbinsX()-2 && biny<hist->GetNbinsY()-2) continue;
            hist->SetBinContent(binx+1,biny+1,abs(mtx(binx,biny)));
            //hist->SetBinContent(binx+1,biny+1,mtx(binx,biny).real());
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
    //mtx_U_r = MakeRealEigenvectors(mtx_U_r);
    //mtx_U_l = MakeRealEigenvectors(mtx_U_l);
    MatrixXcd mtx_H = mtx_U_l.transpose()*mtx_U_r;
    MatrixXcd mtx_lambdanu(mtx_input_lambda.rows(),mtx_input_lambda.cols());
    for (int row=0;row<mtx_input_lambda.rows();row++)
    {
        for (int col=0;col<mtx_input_lambda.cols();col++)
        {
            mtx_lambdanu(row,col) = 0.;
            if (row==col)
            {
                mtx_lambdanu(row,col) = eigensolver_input_lambda.eigenvalues()(col)/mtx_H(row,col);
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

    eigensolver_init = ComplexEigenSolver<MatrixXcd>(mtx_input);
    eigensolver_init_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_U_r_init = eigensolver_init.eigenvectors();
    MatrixXcd mtx_U_l_init = eigensolver_init_transpose.eigenvectors();
    //mtx_U_r_init = MakeRealEigenvectors(mtx_U_r_init);
    //mtx_U_l_init = MakeRealEigenvectors(mtx_U_l_init);

    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    eigensolver_dark_transpose = ComplexEigenSolver<MatrixXcd>(mtx_dark.transpose());
    MatrixXcd mtx_U_r_dark = eigensolver_dark.eigenvectors();
    MatrixXcd mtx_U_l_dark = eigensolver_dark_transpose.eigenvectors();
    //mtx_U_r_dark = MakeRealEigenvectors(mtx_U_r_dark);
    //mtx_U_l_dark = MakeRealEigenvectors(mtx_U_l_dark);

    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
    eigensolver_data_transpose = ComplexEigenSolver<MatrixXcd>(mtx_data.transpose());
    MatrixXcd mtx_U_r_data = eigensolver_data.eigenvectors();
    MatrixXcd mtx_U_l_data = eigensolver_data_transpose.eigenvectors();
    //mtx_U_r_data = MakeRealEigenvectors(mtx_U_r_data);
    //mtx_U_l_data = MakeRealEigenvectors(mtx_U_l_data);

    eigensolver_gamma = ComplexEigenSolver<MatrixXcd>(mtx_gamma_raw);
    eigensolver_gamma_transpose = ComplexEigenSolver<MatrixXcd>(mtx_gamma_raw.transpose());
    MatrixXcd mtx_U_r_gamma = eigensolver_gamma.eigenvectors();
    MatrixXcd mtx_U_l_gamma = eigensolver_gamma_transpose.eigenvectors();
    //mtx_U_r_gamma = MakeRealEigenvectors(mtx_U_r_gamma);
    //mtx_U_l_gamma = MakeRealEigenvectors(mtx_U_l_gamma);

    MatrixXcd mtx_lambdanu = GetLambdaNuMatrix_v2(mtx_input,mtx_input);
    //mtx_lambdanu = CutoffEigenvalueMatrix(mtx_lambdanu, NumberOfRealEigenvectors);

    MatrixXcd mtx_lambdanu_dark = GetLambdaNuMatrix_v2(mtx_dark,mtx_dark);
    //mtx_lambdanu_dark = CutoffEigenvalueMatrix(mtx_lambdanu_dark, NumberOfRealEigenvectors);

    MatrixXcd mtx_lambdanu_data = GetLambdaNuMatrix_v2(mtx_data,mtx_data);
    //mtx_lambdanu_data = CutoffEigenvalueMatrix(mtx_lambdanu_data, NumberOfRealEigenvectors);

    MatrixXcd mtx_lambdanu_gamma = GetLambdaNuMatrix_v2(mtx_gamma_raw,mtx_gamma_raw);
    //mtx_lambdanu_gamma = CutoffEigenvalueMatrix(mtx_lambdanu_gamma, NumberOfRealEigenvectors);

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
            //mtx_eigenvector_vari(row,col) = 0.;
            //mtx_eigenvector_inv_vari(row,col) = 0.;
            mtx_eigenvalue_init(row,col) = mtx_lambdanu(row,col);
            //mtx_eigenvalue_vari(row,col) = 0.;
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
MatrixXcd GetPseudoInverse(MatrixXcd mtx_input, double threshold)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_temp = ComplexEigenSolver<MatrixXcd>(mtx_input);
    MatrixXcd mtx_U = eigensolver_temp.eigenvectors().conjugate().transpose();
    MatrixXcd mtx_V = eigensolver_temp.eigenvectors().inverse().conjugate().transpose();
    MatrixXcd mtx_S = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_S(row,col) = 0.;
            if (row==col && abs(eigensolver_temp.eigenvalues()(col))>threshold)
            {
                mtx_S(row,col) = 1./eigensolver_temp.eigenvalues()(col);
            }
        }
    }
    MatrixXcd mtx_temp = mtx_V*mtx_S*mtx_U;
    return mtx_temp;
}
MatrixXcd GetTruncatedMatrixThreshold(MatrixXcd mtx_input, double threshold)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_temp2 = ComplexEigenSolver<MatrixXcd>(mtx_input);
    MatrixXcd mtx_eigenvector_temp2 = eigensolver_temp2.eigenvectors();
    MatrixXcd mtx_eigenvector_inv_temp2 = eigensolver_temp2.eigenvectors().inverse();
    MatrixXcd mtx_eigenvalue_temp2 = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_eigenvalue_temp2(row,col) = 0.;
            if (row==col && abs(eigensolver_temp2.eigenvalues()(col))>threshold)
            {
                mtx_eigenvalue_temp2(row,col) = eigensolver_temp2.eigenvalues()(col);
            }
        }
    }
    MatrixXcd mtx_temp2 = mtx_eigenvector_temp2*mtx_eigenvalue_temp2*mtx_eigenvector_inv_temp2;
    return mtx_temp2;
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
double GetChi2Function(MatrixXcd mtx_model, MatrixXcd mtx_ref, bool isBlind)
{

    MatrixXcd mtx_data_truncated = mtx_ref;
    MatrixXcd mtx_model_truncated = mtx_model;
    //if (TruncateNoise)
    //{
    //    mtx_data_truncated = GetTruncatedMatrix(mtx_data, NumberOfEigenvectors_Stable);
    //    mtx_model_truncated = GetTruncatedMatrix(mtx_model, NumberOfEigenvectors_Stable);
    //}
    double chi2 = 0.;
    for (int row=0;row<mtx_ref.rows();row++)
    {
        for (int col=0;col<mtx_ref.cols();col++)
        {
            double data = mtx_data_truncated(row,col).real();
            double model = mtx_model_truncated(row,col).real();
            double chi2_this = pow(data-model,2);
            if (isBlind && row<binx_blind_global && col<biny_blind_global)
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
pair<MatrixXcd,VectorXcd> RemoveEmptyRows(MatrixXcd mtx_input,VectorXcd vtr_input)
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
    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_input.rows()-n_empty_rows,mtx_input.cols());
    VectorXcd vtr_output = VectorXcd::Zero(mtx_input.rows()-n_empty_rows);
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
VectorXcd SolutionWithConstraints(MatrixXcd mtx_big, MatrixXcd mtx_constraints_input, VectorXcd vtr_delta, VectorXcd vtr_constraints_delta_input)
{

    MatrixXcd mtx_constraints = RemoveEmptyRows(mtx_constraints_input,vtr_constraints_delta_input).first;
    VectorXcd vtr_constraints_delta = RemoveEmptyRows(mtx_constraints_input,vtr_constraints_delta_input).second;
    std::cout << "mtx_constraints.rows() = " << mtx_constraints.rows() << std::endl;
    if (mtx_constraints.rows()==0)
    {
        BDCSVD<MatrixXcd> svd(mtx_big, ComputeThinU | ComputeThinV);
        svd.setThreshold(svd_threshold);  // size of singular value to be considered as nonzero.
        VectorXcd vtr_vari_big = VectorXcd::Zero(mtx_big.cols());
        vtr_vari_big = svd.solve(vtr_delta);
        return vtr_vari_big;
    }

    MatrixXcd BTB = mtx_big.transpose()*mtx_big;
    VectorXcd BTD = mtx_big.transpose()*vtr_delta;

    MatrixXcd mtx_Bigger = MatrixXcd::Zero(BTB.rows()+mtx_constraints.rows(),BTB.cols()+mtx_constraints.rows());
    mtx_Bigger.block(0,0,BTB.rows(),BTB.cols()) = 2.*BTB;
    mtx_Bigger.block(BTB.rows(),0,mtx_constraints.rows(),mtx_constraints.cols()) = mtx_constraints;
    mtx_Bigger.block(0,BTB.cols(),mtx_constraints.cols(),mtx_constraints.rows()) = mtx_constraints.transpose();

    VectorXcd vtr_bigger_delta = VectorXcd::Zero(BTB.rows()+mtx_constraints.rows());
    vtr_bigger_delta.segment(0,BTB.cols()) = 2.*BTD;
    vtr_bigger_delta.segment(BTB.cols(),vtr_constraints_delta.size()) = vtr_constraints_delta;

    ComplexEigenSolver<MatrixXcd> eigensolver_bigger = ComplexEigenSolver<MatrixXcd>(mtx_Bigger);
    std::cout << "eigensolver_bigger.eigenvalues()(0)  = " << eigensolver_bigger.eigenvalues()(0) << std::endl;
    std::cout << "eigensolver_bigger.eigenvalues()(mtx_Bigger.cols()-1)  = " << eigensolver_bigger.eigenvalues()(mtx_Bigger.cols()-1) << std::endl;
    //svd_threshold = eigensolver_bigger.eigenvalues()(mtx_Bigger.cols()-1).real()/10000.;

    VectorXcd vtr_vari_bigger = VectorXcd::Zero(BTB.cols()+mtx_constraints.rows());
    //vtr_vari_bigger = mtx_Bigger.bdcSvd(ComputeThinU | ComputeThinV).solve(vtr_bigger_delta);
    //JacobiSVD<MatrixXcd> svd(mtx_Bigger, ComputeThinU | ComputeThinV);
    BDCSVD<MatrixXcd> svd(mtx_Bigger, ComputeThinU | ComputeThinV);
    //std::cout << "svd.singularValues()" << std::endl;
    //std::cout << svd.singularValues() << std::endl;
    svd.setThreshold(svd_threshold);  // size of singular value to be considered as nonzero.
    vtr_vari_bigger = svd.solve(vtr_bigger_delta);

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

int DetermineStableNumberOfEigenvalues(MatrixXcd mtx_input)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    int stable_number_max = NumberOfRealEigenvectors;
    int stable_number = 0;
    for (int cutoff=1;cutoff<=stable_number_max;cutoff++)
    {
        double real_lambda = eigensolver_input.eigenvalues()(mtx_input.rows()-cutoff).real();
        double imag_lambda = eigensolver_input.eigenvalues()(mtx_input.rows()-cutoff).imag();
        if (abs(imag_lambda)/abs(real_lambda)>0.01) break;
        if (real_lambda<0.) break;
        stable_number = cutoff;
    }
    return stable_number;
}
int CountEffectiveEquations(MatrixXcd mtx_input, int region)
{
    double threshold = 30000.;
    double evt_count = 0.;
    int eqn_count = 0;
    for (int col=0;col<mtx_input.cols();col++)
    {
        for (int row=0;row<mtx_input.rows();row++)
        {
            if (region==0) // gamma
            {
                if (row>=binx_blind_global) continue;
                if (col>=biny_blind_global) continue;
                if (mtx_input(row,col).real()>threshold)
                {
                    eqn_count += 1;
                }
                else
                {
                    evt_count += mtx_input(row,col).real();
                }
            }
            if (region==1)
            {
                if (row<binx_blind_global) continue;
                if (col>=biny_blind_global) continue;
                if (mtx_input(row,col).real()>threshold)
                {
                    eqn_count += 1;
                }
                else
                {
                    evt_count += mtx_input(row,col).real();
                }
            }
            if (region==2)
            {
                if (row>=binx_blind_global) continue;
                if (col<biny_blind_global) continue;
                if (mtx_input(row,col).real()>threshold)
                {
                    eqn_count += 1;
                }
                else
                {
                    evt_count += mtx_input(row,col).real();
                }
            }
            if (region==3)
            {
                if (row<binx_blind_global) continue;
                if (col<biny_blind_global) continue;
                if (mtx_input(row,col).real()>threshold)
                {
                    eqn_count += 1;
                }
                else
                {
                    evt_count += mtx_input(row,col).real();
                }
            }
            if (region==4)
            {
                if (row>=binx_blind_global && col>=biny_blind_global) continue;
                if (mtx_input(row,col).real()>threshold)
                {
                    eqn_count += 1;
                }
                else
                {
                    evt_count += mtx_input(row,col).real();
                }
            }
        }
    } 
    eqn_count += int(evt_count/threshold);
    return eqn_count;
}
void PrintEigenvalues(MatrixXcd mtx_input)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    for (int entry=1;entry<=mtx_input.rows();entry++)
    {
        std::cout << "eigenvalue " << entry << " = " << abs(eigensolver_input.eigenvalues()(mtx_input.cols()-entry)) << std::endl;
    }
}
void CheckRuoSolution(int entry_size, MatrixXcd mtx_r, MatrixXcd mtx_l, MatrixXcd mtx_dM, VectorXcd vtr_t)
{

    double eta_threshold = 0.0001;
    MatrixXcd mtx_H = mtx_r.transpose()*mtx_l;
    MatrixXcd mtx_FdM = MatrixXcd::Zero(mtx_dM.rows(),mtx_dM.cols());
    MatrixXcd mtx_FFt = MatrixXcd::Zero(mtx_dM.rows(),mtx_dM.cols());
    for (int idx_s=0;idx_s<mtx_dM.cols();idx_s++)
    {
        int sth_entry = idx_s+1;
        std::complex<double> eta_s = mtx_H(mtx_dM.cols()-sth_entry,mtx_dM.cols()-sth_entry);
        if (abs(eta_s)<eta_threshold) continue;
        for (int idx_m=0;idx_m<mtx_dM.cols();idx_m++)
        {
            int mth_entry = idx_m+1;
            std::complex<double> eta_m = mtx_H(mtx_dM.cols()-mth_entry,mtx_dM.cols()-mth_entry);
            if (abs(eta_m)<eta_threshold) continue;
            if (sth_entry>entry_size && mth_entry>entry_size) continue;
            if (sth_entry>5 || mth_entry>5) continue;
            for (int idx_i=0;idx_i<mtx_dM.cols();idx_i++)
            {
                for (int idx_j=0;idx_j<mtx_dM.cols();idx_j++)
                {
                    if (idx_i<binx_blind_global && idx_j<biny_blind_global)
                    {
                        continue;
                    }
                    mtx_FdM(mtx_dM.cols()-sth_entry,mtx_dM.cols()-mth_entry) +=  mtx_dM(idx_i,idx_j)*std::conj(mtx_r(idx_i,mtx_dM.cols()-sth_entry))*std::conj(mtx_l(idx_j,mtx_dM.cols()-mth_entry))/(std::conj(eta_s)*std::conj(eta_m));
                    for (int idx_k=0;idx_k<mtx_dM.cols();idx_k++)
                    {
                        int kth_entry = idx_k+1;
                        std::complex<double> eta_k = mtx_H(mtx_dM.cols()-kth_entry,mtx_dM.cols()-kth_entry);
                        if (abs(eta_k)<eta_threshold) continue;
                        for (int idx_n=0;idx_n<mtx_dM.cols();idx_n++)
                        {
                            int nth_entry = idx_n+1;
                            std::complex<double> eta_n = mtx_H(mtx_dM.cols()-nth_entry,mtx_dM.cols()-nth_entry);
                            if (abs(eta_n)<eta_threshold) continue;
                            if (kth_entry>entry_size && nth_entry>entry_size) continue;
                            if (kth_entry>5 || nth_entry>5) continue;
                            int idx_v = idx_k*mtx_dM.cols() + idx_n;
                            mtx_FFt(mtx_dM.cols()-sth_entry,mtx_dM.cols()-mth_entry) += vtr_t(idx_v)*mtx_r(idx_i,mtx_dM.cols()-kth_entry)*mtx_l(idx_j,mtx_dM.cols()-nth_entry)/(eta_k*eta_n)*std::conj(mtx_r(idx_i,mtx_dM.cols()-sth_entry))*std::conj(mtx_l(idx_j,mtx_dM.cols()-mth_entry))/(std::conj(eta_s)*std::conj(eta_m));
                        }
                    }
                }
            }
        }
    }
    std::cout << "mtx_FdM = " << std::endl;
    std::cout << (mtx_FdM).cwiseAbs() << std::endl;
    std::cout << "mtx_FFt = " << std::endl;
    std::cout << (mtx_FFt).cwiseAbs() << std::endl;

    MatrixXcd mtx_R_full = MatrixXcd::Zero(mtx_r.rows(),mtx_r.cols());
    MatrixXcd mtx_R_gamma = MatrixXcd::Zero(mtx_r.rows(),mtx_r.cols());
    MatrixXcd mtx_L_full = MatrixXcd::Zero(mtx_r.rows(),mtx_r.cols());
    MatrixXcd mtx_L_gamma = MatrixXcd::Zero(mtx_r.rows(),mtx_r.cols());
    for (int idx_k=0;idx_k<mtx_r.rows();idx_k++)
    {
        int kth_entry = idx_k+1;
        for (int idx_n=0;idx_n<mtx_r.cols();idx_n++)
        {
            int nth_entry = idx_n+1;
            //if (kth_entry>entry_size && nth_entry>entry_size) continue;
            //if (kth_entry>5 || nth_entry>5) continue;
            for (int idx_i=0;idx_i<mtx_r.rows();idx_i++)
            {
                mtx_R_full(mtx_r.cols()-kth_entry,mtx_r.cols()-nth_entry) += std::conj(mtx_r(idx_i,mtx_r.cols()-kth_entry))*mtx_r(idx_i,mtx_r.cols()-nth_entry);
                mtx_L_full(mtx_r.cols()-kth_entry,mtx_r.cols()-nth_entry) += std::conj(mtx_l(idx_i,mtx_r.cols()-kth_entry))*mtx_l(idx_i,mtx_r.cols()-nth_entry);
                if (idx_i<binx_blind_global)
                {
                    mtx_R_gamma(mtx_r.cols()-kth_entry,mtx_r.cols()-nth_entry) += std::conj(mtx_r(idx_i,mtx_r.cols()-kth_entry))*mtx_r(idx_i,mtx_r.cols()-nth_entry);
                    mtx_L_gamma(mtx_r.cols()-kth_entry,mtx_r.cols()-nth_entry) += std::conj(mtx_l(idx_i,mtx_r.cols()-kth_entry))*mtx_l(idx_i,mtx_r.cols()-nth_entry);
                }
            }
        }
    }

    MatrixXcd mtx_DM = MatrixXcd::Zero(mtx_r.rows(),mtx_r.cols());
    MatrixXcd mtx_RLDM = MatrixXcd::Zero(mtx_r.rows(),mtx_r.cols());
    MatrixXcd mtx_tRL = MatrixXcd::Zero(mtx_r.rows(),mtx_r.cols());
    for (int idx_s=0;idx_s<mtx_dM.cols();idx_s++)
    {
        int sth_entry = idx_s+1;
        std::complex<double> eta_s = mtx_H(mtx_dM.cols()-sth_entry,mtx_dM.cols()-sth_entry);
        if (abs(eta_s)<eta_threshold) continue;
        for (int idx_m=0;idx_m<mtx_dM.cols();idx_m++)
        {
            int mth_entry = idx_m+1;
            std::complex<double> eta_m = mtx_H(mtx_dM.cols()-mth_entry,mtx_dM.cols()-mth_entry);
            if (abs(eta_m)<eta_threshold) continue;
            if (sth_entry>entry_size && mth_entry>entry_size) continue;
            if (sth_entry>5 || mth_entry>5) continue;
            for (int idx_i=0;idx_i<mtx_dM.cols();idx_i++)
            {
                for (int idx_j=0;idx_j<mtx_dM.cols();idx_j++)
                {
                    if (idx_i<binx_blind_global && idx_j<biny_blind_global)
                    {
                        continue;
                    }
                    mtx_DM(mtx_dM.cols()-sth_entry,mtx_dM.cols()-mth_entry) +=  mtx_dM(idx_i,idx_j)*mtx_l(idx_i,mtx_dM.cols()-sth_entry)*mtx_r(idx_j,mtx_dM.cols()-mth_entry);
                }
            }
        }
    }
    for (int idx_s=0;idx_s<mtx_dM.cols();idx_s++)
    {
        int sth_entry = idx_s+1;
        std::complex<double> eta_s = mtx_H(mtx_dM.cols()-sth_entry,mtx_dM.cols()-sth_entry);
        if (abs(eta_s)<eta_threshold) continue;
        for (int idx_m=0;idx_m<mtx_dM.cols();idx_m++)
        {
            int mth_entry = idx_m+1;
            std::complex<double> eta_m = mtx_H(mtx_dM.cols()-mth_entry,mtx_dM.cols()-mth_entry);
            if (abs(eta_m)<eta_threshold) continue;
            if (sth_entry>entry_size && mth_entry>entry_size) continue;
            if (sth_entry>5 || mth_entry>5) continue;
            for (int idx_k=0;idx_k<mtx_dM.cols();idx_k++)
            {
                int kth_entry = idx_k+1;
                std::complex<double> eta_k = mtx_H(mtx_dM.cols()-kth_entry,mtx_dM.cols()-kth_entry);
                if (abs(eta_k)<eta_threshold) continue;
                for (int idx_n=0;idx_n<mtx_dM.cols();idx_n++)
                {
                    int nth_entry = idx_n+1;
                    std::complex<double> eta_n = mtx_H(mtx_dM.cols()-nth_entry,mtx_dM.cols()-nth_entry);
                    if (abs(eta_n)<eta_threshold) continue;
                    if (kth_entry>entry_size && nth_entry>entry_size) continue;
                    if (kth_entry>5 || nth_entry>5) continue;
                    int idx_v = idx_k*mtx_dM.cols() + idx_n;
                    mtx_RLDM(mtx_dM.cols()-sth_entry,mtx_dM.cols()-mth_entry) +=  mtx_DM(mtx_dM.cols()-kth_entry,mtx_dM.cols()-nth_entry)*mtx_R_full(mtx_dM.cols()-sth_entry,mtx_dM.cols()-kth_entry)*mtx_L_full(mtx_dM.cols()-mth_entry,mtx_dM.cols()-nth_entry)/(eta_k*std::conj(eta_s)*eta_n*std::conj(eta_m));
                    mtx_tRL(mtx_dM.cols()-sth_entry,mtx_dM.cols()-mth_entry) += vtr_t(idx_v)*(mtx_R_full(mtx_dM.cols()-sth_entry,mtx_dM.cols()-kth_entry)*mtx_L_full(mtx_dM.cols()-mth_entry,mtx_dM.cols()-nth_entry)-mtx_R_gamma(mtx_dM.cols()-sth_entry,mtx_dM.cols()-kth_entry)*mtx_L_gamma(mtx_dM.cols()-mth_entry,mtx_dM.cols()-nth_entry))/(eta_k*std::conj(eta_s)*eta_n*std::conj(eta_m));
                }
            }
        }
    }
    std::cout << "mtx_RLDM = " << std::endl;
    std::cout << (mtx_RLDM).cwiseAbs() << std::endl;
    std::cout << "mtx_tRL = " << std::endl;
    std::cout << (mtx_tRL).cwiseAbs() << std::endl;

}
MatrixXcd SortEigenvectors(VectorXcd eigenval, MatrixXcd mtx_r, MatrixXcd mtx_l)
{
    double eta_threshold = 0.0001;
    MatrixXcd mtx_l_output = mtx_l;
    MatrixXcd mtx_H = mtx_r.transpose()*mtx_l;
    for (int entry=1;entry<=mtx_r.cols();entry++)
    {
        std::complex<double> eta = mtx_H(mtx_r.cols()-entry,mtx_r.cols()-entry);
        if (abs(eta)<eta_threshold && entry<mtx_r.cols())
        {
            if (abs(abs(eigenval(mtx_r.cols()-entry))-abs(eigenval(mtx_r.cols()-entry-1)))<0.1)
            {
                mtx_l_output.col(mtx_r.cols()-entry) = mtx_l.col(mtx_r.cols()-entry-1);
                mtx_l_output.col(mtx_r.cols()-entry-1) = mtx_l.col(mtx_r.cols()-entry);
            }
        }
    }
    return mtx_l_output;
}

MatrixXcd NuclearNormMinimization(MatrixXcd mtx_init_input, MatrixXcd mtx_data_input, MatrixXcd mtx_dark_input, int entry_start, int entry_size, double step_frac, bool isBlind, int isLeft, double alpha)
{
    JacobiSVD<MatrixXd> svd_data(mtx_data_input.real(), ComputeFullU | ComputeFullV);
    MatrixXd mtx_U_data = svd_data.matrixU();
    MatrixXd mtx_V_data = svd_data.matrixV();
    MatrixXd mtx_S_data = MatrixXd::Zero(mtx_data_input.rows(),mtx_data_input.cols());
    for (int entry=0;entry<svd_data.singularValues().size();entry++)
    {
        mtx_S_data(entry,entry) = svd_data.singularValues()(entry);
    }

    JacobiSVD<MatrixXd> svd_init(mtx_init_input.real(), ComputeFullU | ComputeFullV);
    MatrixXd mtx_U_init = svd_init.matrixU();
    MatrixXd mtx_V_init = svd_init.matrixV();
    MatrixXd mtx_S_init = MatrixXd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int entry=0;entry<svd_init.singularValues().size();entry++)
    {
        mtx_S_init(entry,entry) = svd_init.singularValues()(entry);
    }

    JacobiSVD<MatrixXd> svd_dark(mtx_dark_input.real(), ComputeFullU | ComputeFullV);
    MatrixXd mtx_U_dark = svd_dark.matrixU();
    MatrixXd mtx_V_dark = svd_dark.matrixV();
    MatrixXd mtx_S_dark = MatrixXd::Zero(mtx_dark_input.rows(),mtx_dark_input.cols());
    for (int entry=0;entry<svd_dark.singularValues().size();entry++)
    {
        mtx_S_dark(entry,entry) = svd_dark.singularValues()(entry);
    }
    //MatrixXd mtx_UTV_dark = mtx_U_dark.transpose()*mtx_V_dark;
    //std::cout << "mtx_UTV_dark:" << std::endl;
    //std::cout << mtx_UTV_dark << std::endl;
    //MatrixXd mtx_UTU_dark = mtx_U_dark.transpose()*mtx_U_dark;
    //std::cout << "mtx_UTU_dark:" << std::endl;
    //std::cout << mtx_UTU_dark << std::endl;
    //MatrixXd mtx_VTV_dark = mtx_V_dark.transpose()*mtx_V_dark;
    //std::cout << "mtx_VTV_dark:" << std::endl;
    //std::cout << mtx_VTV_dark << std::endl;

    int size_k = mtx_init_input.cols();
    int size_n = mtx_init_input.cols();
    int length_tkn = size_k*size_n;
    int regularization_size = length_tkn;
    if (RegularizationType==1 || RegularizationType==2 || RegularizationType==3)
    {
        regularization_size = 0;
    }
    if (!isBlind)
    {
        regularization_size = 0;
    }
    VectorXcd vtr_Delta = VectorXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols()+regularization_size);
    MatrixXcd mtx_A = MatrixXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols()+regularization_size,length_tkn);
    for (int idx_i=0;idx_i<mtx_init_input.rows();idx_i++)
    {
        for (int idx_j=0;idx_j<mtx_init_input.rows();idx_j++)
        {
            int idx_u = idx_j*mtx_init_input.rows() + idx_i;
            if (isBlind && idx_i<binx_blind_global && idx_j<biny_blind_global)
            {
                continue;
            }
            vtr_Delta(idx_u) = (mtx_data_input-mtx_init_input)(idx_i,idx_j);
            for (int idx_k=0;idx_k<size_k;idx_k++)
            {
                int kth_entry = idx_k+1;
                for (int idx_n=0;idx_n<size_n;idx_n++)
                {
                    int nth_entry = idx_n+1;
                    if (kth_entry>entry_size && nth_entry>entry_size) continue;
                    if (RegularizationType==1)
                    {
                        if (kth_entry>entry_size && nth_entry>entry_size) continue;
                        if (kth_entry==3 && nth_entry==3) continue;
                        if (kth_entry>3 || nth_entry>3) continue;
                    }
                    if (RegularizationType==2)
                    {
                        if (kth_entry>entry_size && nth_entry>entry_size) continue;
                        if (kth_entry==3 && nth_entry==3) continue;
                        if (kth_entry>5 || nth_entry>5) continue;
                    }
                    if (isLeft==1 && kth_entry>entry_size) continue;
                    if (isLeft==2 && nth_entry>entry_size) continue;
                    int idx_v = idx_k*size_n + idx_n;
                    mtx_A(idx_u,idx_v) = mtx_U_dark(idx_i,idx_k)*mtx_V_dark(idx_j,idx_n);
                }
            }
        }
    }
    if (regularization_size==length_tkn)
    {
        for (int idx_k=0;idx_k<size_k;idx_k++)
        {
            int kth_entry = idx_k+1;
            for (int idx_n=0;idx_n<size_n;idx_n++)
            {
                int nth_entry = idx_n+1;
                int idx_v = idx_k*size_n + idx_n;
                int idx_u = idx_v + mtx_init_input.rows()*mtx_init_input.cols();
                if (RegularizationType==0)
                {
                    double regularization_weight = abs(mtx_regularization(idx_k,idx_n));
                    if (regularization_weight==0.) continue;
                    mtx_A(idx_u,idx_v) = pow(alpha*1./regularization_weight,1);
                }
                else if (RegularizationType==4)
                {
                    mtx_A(idx_u,idx_v) = alpha*1e5;
                }
            }
        }
    }
    VectorXcd vtr_t = VectorXcd::Zero(length_tkn);
    BDCSVD<MatrixXcd> bdc_svd(mtx_A, ComputeThinU | ComputeThinV);
    //bdc_svd.setThreshold(1e-10);  // size of singular value to be considered as nonzero.
    //std::cout << "bdc_svd.singularValues():" << std::endl;
    //std::cout << bdc_svd.singularValues() << std::endl;
    if (RegularizationType==3)
    {
        MatrixXcd mtx_U = bdc_svd.matrixU();
        MatrixXcd mtx_V = bdc_svd.matrixV();
        MatrixXcd mtx_S = MatrixXcd::Zero(mtx_A.rows(),mtx_A.cols());
        MatrixXcd mtx_S_inv = MatrixXcd::Zero(mtx_A.cols(),mtx_A.rows());
        for (int entry=0;entry<bdc_svd.singularValues().size();entry++)
        {
            if (abs(bdc_svd.singularValues()(entry))==0.) continue;
            if (entry>int(alpha-1.)) continue;
            mtx_S(entry,entry) = bdc_svd.singularValues()(entry);
            mtx_S_inv(entry,entry) = 1./(bdc_svd.singularValues()(entry));
        }
        vtr_t = mtx_V*mtx_S_inv*mtx_U.conjugate().transpose()*vtr_Delta;
    }
    else
    {
        vtr_t = bdc_svd.solve(vtr_Delta);
    }
    vtr_eigenval_vvv = bdc_svd.singularValues();

    MatrixXcd mtx_t = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_C = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_D = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_E = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int idx_k=0;idx_k<size_k;idx_k++)
    {
        int kth_entry = idx_k+1;
        std::complex<double> sigma_k = mtx_S_dark(idx_k,idx_k);
        for (int idx_n=0;idx_n<size_n;idx_n++)
        {
            int nth_entry = idx_n+1;
            std::complex<double> sigma_n = mtx_S_dark(idx_n,idx_n);
            int idx_kn = idx_k*size_n + idx_n;
            int idx_nk = idx_n*size_k + idx_k;
            mtx_t(idx_k,idx_n) = vtr_t(idx_kn); 
            if (kth_entry!=nth_entry)
            {
                mtx_C(idx_k,idx_n) = (vtr_t(idx_kn)/sigma_k+vtr_t(idx_nk)/sigma_n)/(sigma_n/sigma_k-sigma_k/sigma_n); 
                mtx_D(idx_k,idx_n) = (vtr_t(idx_kn)/sigma_n+vtr_t(idx_nk)/sigma_k)/(sigma_n/sigma_k-sigma_k/sigma_n); 
            }
            else
            {
                mtx_E(idx_k,idx_n) = vtr_t(idx_kn); 
            }
        }
    }
    if (!isBlind)
    {
        return mtx_t;
    }

    MatrixXcd mtx_U_vari = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_V_vari = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_S_vari = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int idx_i=0;idx_i<mtx_init_input.cols();idx_i++)
    {
        for (int idx_k=0;idx_k<size_k;idx_k++)
        {
            int kth_entry = idx_k+1;
            mtx_S_vari(idx_k,idx_k) = mtx_E(idx_k,idx_k);
            for (int idx_n=0;idx_n<size_n;idx_n++)
            {
                int nth_entry = idx_n+1;
                std::complex<double> coeff_u = mtx_C(idx_k,idx_n);
                mtx_U_vari(idx_i,idx_n) += coeff_u*mtx_U_dark(idx_i,idx_k);
                std::complex<double> coeff_v = mtx_D(idx_k,idx_n);
                mtx_V_vari(idx_i,idx_n) += coeff_v*mtx_V_dark(idx_i,idx_k);
            }
        }
    }


    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    mtx_output = mtx_U_init*mtx_S_init*mtx_V_init.transpose();
    MatrixXcd mtx_vari = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    mtx_vari += step_frac*mtx_U_init*mtx_S_init*mtx_V_vari.transpose();
    mtx_vari += step_frac*mtx_U_vari*mtx_S_init*mtx_V_init.transpose();
    mtx_vari += step_frac*mtx_U_init*mtx_S_vari*mtx_V_init.transpose();
    mtx_output += mtx_vari;

    return mtx_output;

}

MatrixXcd LowRankOptimizationMethod(MatrixXcd mtx_init_input, MatrixXcd mtx_data_input, MatrixXcd mtx_dark_input, int entry_start, int entry_size, double step_frac, bool isBlind, int isLeft, double alpha)
{

    bool doPrint = false;
    double eta_threshold = 0.0001;

    eigensolver_init = ComplexEigenSolver<MatrixXcd>(mtx_init_input);
    eigensolver_init_transpose = ComplexEigenSolver<MatrixXcd>(mtx_init_input.transpose());
    MatrixXcd mtx_r_init = eigensolver_init.eigenvectors();
    MatrixXcd mtx_l_init = eigensolver_init_transpose.eigenvectors();
    mtx_l_init = SortEigenvectors(eigensolver_init.eigenvalues(), mtx_r_init, mtx_l_init);
    MatrixXcd mtx_H_init = mtx_r_init.transpose()*mtx_l_init;
    MatrixXcd mtx_lambda_init(mtx_init_input.rows(),mtx_init_input.cols());
    for (int row=0;row<mtx_init_input.rows();row++)
    {
        for (int col=0;col<mtx_init_input.cols();col++)
        {
            mtx_lambda_init(row,col) = 0.;
            if (row==col)
            {
                mtx_lambda_init(row,col) = eigensolver_init.eigenvalues()(col);
            }
        }
    }

    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark_input);
    eigensolver_dark_transpose = ComplexEigenSolver<MatrixXcd>(mtx_dark_input.transpose());
    MatrixXcd mtx_r_dark = eigensolver_dark.eigenvectors();
    MatrixXcd mtx_l_dark = eigensolver_dark_transpose.eigenvectors();
    mtx_l_dark = SortEigenvectors(eigensolver_dark.eigenvalues(), mtx_r_dark, mtx_l_dark);
    MatrixXcd mtx_H_dark = mtx_r_dark.transpose()*mtx_l_dark;
    MatrixXcd mtx_lambda_dark(mtx_init_input.rows(),mtx_init_input.cols());
    for (int row=0;row<mtx_init_input.rows();row++)
    {
        for (int col=0;col<mtx_init_input.cols();col++)
        {
            mtx_lambda_dark(row,col) = 0.;
            if (row==col)
            {
                mtx_lambda_dark(row,col) = eigensolver_dark.eigenvalues()(col);
            }
        }
    }
    MatrixXcd mtx_S_init = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int col=0;col<mtx_init_input.cols();col++)
    {
        mtx_S_init(col,col) = mtx_lambda_init(col,col)/mtx_H_dark(col,col);
    }

    MatrixXcd mtx_Delta = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int idx_k=0;idx_k<mtx_init_input.rows();idx_k++)
    {
        int kth_entry = idx_k+1;
        for (int idx_n=0;idx_n<mtx_init_input.cols();idx_n++)
        {
            int nth_entry = idx_n+1;
            //if (kth_entry>entry_size && nth_entry>entry_size) continue;
            //if (kth_entry>5 || nth_entry>5) continue;
            for (int idx_i=0;idx_i<mtx_init_input.rows();idx_i++)
            {
                for (int idx_j=0;idx_j<mtx_init_input.rows();idx_j++)
                {
                    if (isBlind && idx_i<binx_blind_global && idx_j<biny_blind_global)
                    {
                        continue;
                    }
                    mtx_Delta(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-nth_entry) += mtx_l_dark(idx_i,mtx_init_input.cols()-kth_entry)*(mtx_data_input-mtx_init_input)(idx_i,idx_j)*mtx_r_dark(idx_j,mtx_init_input.cols()-nth_entry);
                }
            }
        }
    }

    MatrixXcd mtx_W = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int idx_i=0;idx_i<mtx_init_input.rows();idx_i++)
    {
        if (idx_i<binx_blind_global)
        {
            mtx_W(idx_i,idx_i) = 1.;
        }
    }
    MatrixXcd mtx_R_full = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_Ur2 = mtx_r_dark.conjugate().transpose()*mtx_r_dark;
    MatrixXcd mtx_WUr2 = (mtx_W*mtx_r_dark).conjugate().transpose()*(mtx_W*mtx_r_dark);
    MatrixXcd mtx_UrWUr = mtx_r_dark.conjugate().transpose()*mtx_W*mtx_r_dark;
    MatrixXcd mtx_R_gamma = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_L_full = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_Ul2 = mtx_l_dark.conjugate().transpose()*mtx_l_dark;
    MatrixXcd mtx_WUl2 = (mtx_W*mtx_l_dark).conjugate().transpose()*(mtx_W*mtx_l_dark);
    MatrixXcd mtx_UlWUl = mtx_l_dark.conjugate().transpose()*mtx_W*mtx_l_dark;
    MatrixXcd mtx_L_gamma = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int idx_k=0;idx_k<mtx_init_input.rows();idx_k++)
    {
        int kth_entry = idx_k+1;
        for (int idx_n=0;idx_n<mtx_init_input.cols();idx_n++)
        {
            int nth_entry = idx_n+1;
            //if (kth_entry>entry_size && nth_entry>entry_size) continue;
            //if (kth_entry>5 || nth_entry>5) continue;
            for (int idx_i=0;idx_i<mtx_init_input.rows();idx_i++)
            {
                mtx_R_full(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-nth_entry) += std::conj(mtx_r_dark(idx_i,mtx_init_input.cols()-kth_entry))*mtx_r_dark(idx_i,mtx_init_input.cols()-nth_entry);
                mtx_L_full(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-nth_entry) += std::conj(mtx_l_dark(idx_i,mtx_init_input.cols()-kth_entry))*mtx_l_dark(idx_i,mtx_init_input.cols()-nth_entry);
                if (isBlind) 
                {
                    if (idx_i<binx_blind_global)
                    {
                        mtx_R_gamma(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-nth_entry) += std::conj(mtx_r_dark(idx_i,mtx_init_input.cols()-kth_entry))*mtx_r_dark(idx_i,mtx_init_input.cols()-nth_entry);
                        mtx_L_gamma(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-nth_entry) += std::conj(mtx_l_dark(idx_i,mtx_init_input.cols()-kth_entry))*mtx_l_dark(idx_i,mtx_init_input.cols()-nth_entry);
                    }
                }
            }
        }
    }
    if (doPrint)
    {
        std::cout << "mtx_R_full:" << std::endl;
        std::cout << mtx_R_full.cwiseAbs() << std::endl;
        std::cout << "mtx_Ur2:" << std::endl;
        std::cout << mtx_Ur2.cwiseAbs() << std::endl;
        std::cout << "mtx_L_full:" << std::endl;
        std::cout << mtx_L_full.cwiseAbs() << std::endl;
        std::cout << "mtx_Ul2:" << std::endl;
        std::cout << mtx_Ul2.cwiseAbs() << std::endl;
        std::cout << "mtx_R_gamma:" << std::endl;
        std::cout << mtx_R_gamma.cwiseAbs() << std::endl;
        std::cout << "mtx_WUr2:" << std::endl;
        std::cout << mtx_WUr2.cwiseAbs() << std::endl;
        std::cout << "mtx_UrWUr:" << std::endl;
        std::cout << mtx_UrWUr.cwiseAbs() << std::endl;
        std::cout << "mtx_L_gamma:" << std::endl;
        std::cout << mtx_L_gamma.cwiseAbs() << std::endl;
        std::cout << "mtx_WUl2:" << std::endl;
        std::cout << mtx_WUl2.cwiseAbs() << std::endl;
        std::cout << "mtx_UlWUl:" << std::endl;
        std::cout << mtx_UlWUl.cwiseAbs() << std::endl;
        MatrixXcd mtx_O_full = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
        MatrixXcd mtx_O_gamma_P = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
        MatrixXcd mtx_O_gamma_Q = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
        MatrixXcd mtx_LgR = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
        MatrixXcd mtx_LgRg = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
        MatrixXcd mtx_LRg = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
        MatrixXcd mtx_UlWUr = mtx_l_dark.conjugate().transpose()*mtx_W*mtx_r_dark;
        for (int idx_k=0;idx_k<mtx_init_input.rows();idx_k++)
        {
            for (int idx_n=0;idx_n<mtx_init_input.cols();idx_n++)
            {
                for (int idx_m=0;idx_m<mtx_init_input.cols();idx_m++)
                {
                    mtx_O_full(idx_k,idx_n) += mtx_L_full(idx_m,idx_k)*mtx_R_full(idx_m,idx_n)/(mtx_H_dark(idx_k,idx_k)*std::conj(mtx_H_dark(idx_m,idx_m)));
                    mtx_O_gamma_P(idx_k,idx_n) += mtx_L_gamma(idx_m,idx_k)*mtx_R_gamma(idx_m,idx_n)/(mtx_H_dark(idx_k,idx_k)*std::conj(mtx_H_dark(idx_m,idx_m)));
                    mtx_O_gamma_Q(idx_k,idx_n) += mtx_L_gamma(idx_m,idx_k)*mtx_R_gamma(idx_m,idx_n)/(mtx_H_dark(idx_n,idx_n)*std::conj(mtx_H_dark(idx_m,idx_m)));
                    mtx_LgRg(idx_k,idx_n) += mtx_L_gamma(idx_m,idx_k)*mtx_R_gamma(idx_m,idx_n)/(std::conj(mtx_H_dark(idx_m,idx_m)));
                    mtx_LgR(idx_k,idx_n) += mtx_L_gamma(idx_m,idx_k)*mtx_R_full(idx_m,idx_n)/(std::conj(mtx_H_dark(idx_m,idx_m)));
                    mtx_LRg(idx_k,idx_n) += mtx_L_full(idx_m,idx_k)*mtx_R_gamma(idx_m,idx_n)/(std::conj(mtx_H_dark(idx_m,idx_m)));
                }
            }
        }
        std::cout << "mtx_O_full:" << std::endl;
        std::cout << mtx_O_full.cwiseAbs() << std::endl;
        std::cout << "mtx_O_gamma_P:" << std::endl;
        std::cout << mtx_O_gamma_P.cwiseAbs() << std::endl;
        std::cout << "mtx_O_gamma_Q:" << std::endl;
        std::cout << mtx_O_gamma_Q.cwiseAbs() << std::endl;
        std::cout << "mtx_LgRg:" << std::endl;
        std::cout << mtx_LgRg.cwiseAbs() << std::endl;
        std::cout << "mtx_LgR:" << std::endl;
        std::cout << mtx_LgR.cwiseAbs() << std::endl;
        std::cout << "mtx_LRg:" << std::endl;
        std::cout << mtx_LRg.cwiseAbs() << std::endl;
        std::cout << "mtx_UlWUr:" << std::endl;
        std::cout << mtx_UlWUr.cwiseAbs() << std::endl;
    }

    MatrixXcd mtx_P = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_Q = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_P_swap = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_Q_swap = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int idx_k=0;idx_k<mtx_init_input.rows();idx_k++)
    {
        std::complex<double> eta_k = mtx_H_dark(idx_k,idx_k);
        if (abs(eta_k)<eta_threshold) continue;
        for (int idx_n=0;idx_n<mtx_init_input.rows();idx_n++)
        {
            std::complex<double> eta_n = mtx_H_dark(idx_n,idx_n);
            if (abs(eta_n)<eta_threshold) continue;
            for (int idx_s=0;idx_s<mtx_init_input.rows();idx_s++)
            {
                std::complex<double> eta_s = mtx_H_dark(idx_s,idx_s);
                if (abs(eta_s)<eta_threshold) continue;
                mtx_P(idx_n,idx_k) += mtx_L_gamma(idx_s,idx_n)*mtx_R_full(idx_s,idx_k)/(eta_n*std::conj(eta_s)); 
                mtx_Q(idx_n,idx_k) += mtx_L_full(idx_s,idx_n)*mtx_R_gamma(idx_s,idx_k)/(eta_k*std::conj(eta_s)); 
                mtx_P_swap(idx_n,idx_k) += mtx_L_full(idx_s,idx_n)*mtx_R_gamma(idx_s,idx_k)/(eta_n*std::conj(eta_s)); 
                mtx_Q_swap(idx_n,idx_k) += mtx_L_gamma(idx_s,idx_n)*mtx_R_full(idx_s,idx_k)/(eta_k*std::conj(eta_s)); 
            }
        }
    }
    if (doPrint)
    {
        std::cout << "mtx_P:" << std::endl;
        std::cout << mtx_P.cwiseAbs() << std::endl;
        std::cout << "mtx_Q:" << std::endl;
        std::cout << mtx_Q.cwiseAbs() << std::endl;
        std::cout << "mtx_P_swap:" << std::endl;
        std::cout << mtx_P_swap.cwiseAbs() << std::endl;
        std::cout << "mtx_Q_Swap:" << std::endl;
        std::cout << mtx_Q_swap.cwiseAbs() << std::endl;
        MatrixXcd mtx_PP = mtx_P*mtx_P-mtx_P;
        MatrixXcd mtx_QQ = mtx_Q*mtx_Q-mtx_Q;
        std::cout << "mtx_PP:" << std::endl;
        std::cout << mtx_PP.cwiseAbs() << std::endl;
        std::cout << "mtx_QQ:" << std::endl;
        std::cout << mtx_QQ.cwiseAbs() << std::endl;
    }

    int size_k = mtx_init_input.cols();
    int size_n = mtx_init_input.cols();
    int length_tkn = size_k*size_n;
    VectorXcd vtr_Delta = VectorXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols());
    MatrixXcd mtx_F = MatrixXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols(),length_tkn);
    for (int idx_i=0;idx_i<mtx_init_input.rows();idx_i++)
    {
        int ith_entry = idx_i+1;
        for (int idx_j=0;idx_j<mtx_init_input.rows();idx_j++)
        {
            int jth_entry = idx_j+1;
            int idx_u = idx_j*mtx_init_input.rows() + idx_i;
            vtr_Delta(idx_u) = mtx_Delta(mtx_init_input.cols()-jth_entry,mtx_init_input.cols()-ith_entry);
            for (int idx_k=0;idx_k<size_k;idx_k++)
            {
                int kth_entry = idx_k+1;
                std::complex<double> eta_k = mtx_H_dark(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-kth_entry);
                for (int idx_n=0;idx_n<size_n;idx_n++)
                {
                    int nth_entry = idx_n+1;
                    std::complex<double> eta_n = mtx_H_dark(mtx_init_input.cols()-nth_entry,mtx_init_input.cols()-nth_entry);
                    //if (kth_entry>entry_size && nth_entry>entry_size) continue;
                    //if (kth_entry>5 || nth_entry>5) continue;
                    int idx_v = idx_k*size_n + idx_n;
                    mtx_F(idx_u,idx_v) += -1.*mtx_Q(mtx_init_input.cols()-jth_entry,mtx_init_input.cols()-kth_entry)*mtx_P(mtx_init_input.cols()-ith_entry,mtx_init_input.cols()-nth_entry);
                    if (jth_entry==kth_entry && ith_entry==nth_entry)
                    {
                        mtx_F(idx_u,idx_v) += 1.;
                    }
                }
            }
        }
    }

    int regularization_size = length_tkn;
    if (RegularizationType==1 || RegularizationType==2 || RegularizationType==3 || RegularizationType==5)
    {
        regularization_size = 0;
    }
    VectorXcd vtr_Delta_Ruo = VectorXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols()+regularization_size);
    MatrixXcd mtx_F_Ruo = MatrixXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols()+regularization_size,length_tkn);
    for (int idx_i=0;idx_i<mtx_init_input.rows();idx_i++)
    {
        int ith_entry = idx_i+1;
        for (int idx_j=0;idx_j<mtx_init_input.rows();idx_j++)
        {
            int jth_entry = idx_j+1;
            int idx_u = idx_j*mtx_init_input.rows() + idx_i;
            if (isBlind && idx_i<binx_blind_global && idx_j<biny_blind_global)
            {
                continue;
            }
            vtr_Delta_Ruo(idx_u) = (mtx_data_input-mtx_init_input)(idx_i,idx_j);
            for (int idx_k=0;idx_k<size_k;idx_k++)
            {
                int kth_entry = idx_k+1;
                std::complex<double> eta_k = mtx_H_dark(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-kth_entry);
                if (abs(eta_k)<eta_threshold) continue;
                for (int idx_n=0;idx_n<size_n;idx_n++)
                {
                    int nth_entry = idx_n+1;
                    std::complex<double> eta_n = mtx_H_dark(mtx_init_input.cols()-nth_entry,mtx_init_input.cols()-nth_entry);
                    if (abs(eta_n)<eta_threshold) continue;
                    if (RegularizationType!=5) 
                    {
                        if (kth_entry>entry_size && nth_entry>entry_size) continue;
                    }
                    if (RegularizationType==1)
                    {
                        if (kth_entry>entry_size && nth_entry>entry_size) continue;
                        if (kth_entry==3 && nth_entry==3) continue;
                        if (kth_entry>3 || nth_entry>3) continue;
                    }
                    if (RegularizationType==2)
                    {
                        if (kth_entry>entry_size && nth_entry>entry_size) continue;
                        if (kth_entry==3 && nth_entry==3) continue;
                        if (kth_entry>5 || nth_entry>5) continue;
                    }
                    //if (isLeft==1 && kth_entry>entry_size) continue;
                    //if (isLeft==2 && nth_entry>entry_size) continue;
                    int idx_v = idx_k*size_n + idx_n;
                    mtx_F_Ruo(idx_u,idx_v) = mtx_r_dark(idx_i,mtx_init_input.cols()-kth_entry)*mtx_l_dark(idx_j,mtx_init_input.cols()-nth_entry)/(eta_k*eta_n);
                }
            }
        }
    }
    if (regularization_size==length_tkn)
    {
        for (int idx_k=0;idx_k<size_k;idx_k++)
        {
            int kth_entry = idx_k+1;
            for (int idx_n=0;idx_n<size_n;idx_n++)
            {
                int nth_entry = idx_n+1;
                std::complex<double> eta_n = mtx_H_dark(mtx_init_input.cols()-nth_entry,mtx_init_input.cols()-nth_entry);
                int idx_v = idx_k*size_n + idx_n;
                int idx_u = idx_v + mtx_init_input.rows()*mtx_init_input.cols();
                if (RegularizationType==0)
                {
                    double regularization_weight = abs(mtx_regularization(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-nth_entry));
                    //regularization_weight = regularization_weight*pow(exposure_hours/100.,2);
                    if (regularization_weight==0.) continue;
                    mtx_F_Ruo(idx_u,idx_v) = alpha*1./regularization_weight;
                    //mtx_F_Ruo(idx_u,idx_v) = pow(alpha*1./regularization_weight,2);
                }
                else if (RegularizationType==4)
                {
                    mtx_F_Ruo(idx_u,idx_v) = alpha*1e5;
                }
            }
        }
    }
    VectorXcd vtr_t_Ruo = VectorXcd::Zero(length_tkn);
    BDCSVD<MatrixXcd> bdc_svd(mtx_F_Ruo, ComputeThinU | ComputeThinV);
    //std::cout << "bdc_svd.singularValues():" << std::endl;
    //std::cout << bdc_svd.singularValues() << std::endl;
    if (RegularizationType==3)
    {
        MatrixXcd mtx_U = bdc_svd.matrixU();
        MatrixXcd mtx_V = bdc_svd.matrixV();
        MatrixXcd mtx_S = MatrixXcd::Zero(mtx_F_Ruo.rows(),mtx_F_Ruo.cols());
        MatrixXcd mtx_S_inv = MatrixXcd::Zero(mtx_F_Ruo.cols(),mtx_F_Ruo.rows());
        for (int entry=0;entry<bdc_svd.singularValues().size();entry++)
        {
            if (abs(bdc_svd.singularValues()(entry))==0.) continue;
            if (entry>int(alpha-1.)) continue;
            mtx_S(entry,entry) = bdc_svd.singularValues()(entry);
            mtx_S_inv(entry,entry) = 1./(bdc_svd.singularValues()(entry));
        }
        vtr_t_Ruo = mtx_V*mtx_S_inv*mtx_U.conjugate().transpose()*vtr_Delta_Ruo;
        vtr_eigenval_vvv = bdc_svd.singularValues();
    }
    else if (RegularizationType==5)
    {
        MatrixXcd mtx_FTF_Ruo = mtx_F_Ruo.conjugate().transpose()*mtx_F_Ruo;
        ComplexEigenSolver<MatrixXcd> eigensolver_vvv = ComplexEigenSolver<MatrixXcd>(mtx_FTF_Ruo);
        MatrixXcd mtx_U = eigensolver_vvv.eigenvectors();
        MatrixXcd mtx_U_inv = eigensolver_vvv.eigenvectors().inverse();
        MatrixXcd mtx_S = MatrixXcd::Zero(mtx_F_Ruo.rows(),mtx_F_Ruo.cols());
        MatrixXcd mtx_S_inv = MatrixXcd::Zero(mtx_F_Ruo.cols(),mtx_F_Ruo.rows());
        //std::cout << "eigensolver_vvv.eigenvalues().size() = " << eigensolver_vvv.eigenvalues().size() << std::endl;
        //std::cout << "eigensolver_vvv.eigenvalues() = " << eigensolver_vvv.eigenvalues() << std::endl;
        for (int entry=0;entry<eigensolver_vvv.eigenvalues().size();entry++)
        {
            vtr_eigenval_vvv(entry) = abs(eigensolver_vvv.eigenvalues()(entry));
            if (abs(eigensolver_vvv.eigenvalues()(entry))==0.) continue;
            if (entry<eigensolver_vvv.eigenvalues().size()-int(alpha)) continue;
            mtx_S(entry,entry) = eigensolver_vvv.eigenvalues()(entry);
            mtx_S_inv(entry,entry) = 1./eigensolver_vvv.eigenvalues()(entry);
        }
        vtr_t_Ruo = mtx_U*mtx_S_inv*mtx_U_inv*mtx_F_Ruo.conjugate().transpose()*vtr_Delta_Ruo;
    }
    else
    {
        vtr_t_Ruo = bdc_svd.solve(vtr_Delta_Ruo);
        vtr_eigenval_vvv = bdc_svd.singularValues();
    }
    //CheckRuoSolution(entry_size, mtx_r_dark,mtx_l_dark, mtx_data_input-mtx_init_input, vtr_t_Ruo);

    double sum_F_abs = 0.;
    MatrixXcd mtx_FFmF = mtx_F*mtx_F-mtx_F;
    for (int idx_i=0;idx_i<mtx_F.cols();idx_i++)
    {
        for (int idx_j=0;idx_j<mtx_F.cols();idx_j++)
        {
            sum_F_abs += abs(mtx_FFmF(idx_i,idx_j));
        }
    }
    if (doPrint)
    {
        std::cout << "sum abs(mtx_F*mtx_F - mtx_F):" << std::endl;
        std::cout << sum_F_abs << std::endl;
    }

    VectorXcd vtr_t = VectorXcd::Zero(length_tkn);
    vtr_t = vtr_t_Ruo;

    MatrixXcd mtx_C = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int idx_k=0;idx_k<size_k;idx_k++)
    {
        int kth_entry = idx_k+1;
        std::complex<double> lambda_k = mtx_lambda_dark(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-kth_entry);
        std::complex<double> eta_k = mtx_H_dark(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-kth_entry);
        for (int idx_n=0;idx_n<size_n;idx_n++)
        {
            int nth_entry = idx_n+1;
            //if (kth_entry>entry_size && nth_entry>entry_size) continue;
            //if (kth_entry>5 || nth_entry>5) continue;
            std::complex<double> lambda_n = mtx_lambda_dark(mtx_init_input.cols()-nth_entry,mtx_init_input.cols()-nth_entry);
            std::complex<double> eta_n = mtx_H_dark(mtx_init_input.cols()-nth_entry,mtx_init_input.cols()-nth_entry);
            int idx_v = idx_k*size_n + idx_n;
            if (kth_entry!=nth_entry)
            {
                mtx_C(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-nth_entry) = vtr_t(idx_v)/(lambda_k-lambda_n); 
            }
            else
            {
                if (abs(eta_k)<eta_threshold) continue;
                mtx_C(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-nth_entry) = vtr_t(idx_v)/(eta_k); 
            }
        }
    }
    if (doPrint)
    {
        int print_size = 6;
        std::cout << "mtx_C:" << std::endl;
        std::cout << mtx_C.block(mtx_init_input.rows()-print_size,mtx_init_input.cols()-print_size,print_size,print_size).cwiseAbs() << std::endl;
        std::cout << "mtx_regularization:" << std::endl;
        std::cout << mtx_regularization.block(mtx_init_input.rows()-print_size,mtx_init_input.cols()-print_size,print_size,print_size).cwiseAbs() << std::endl;
    }

    MatrixXcd mtx_l_vari = MatrixXcd::Zero(mtx_l_init.rows(),mtx_l_init.cols());
    MatrixXcd mtx_r_vari = MatrixXcd::Zero(mtx_r_init.rows(),mtx_r_init.cols());
    MatrixXcd mtx_S_vari = MatrixXcd::Zero(mtx_r_init.rows(),mtx_r_init.cols());
    for (int idx_i=0;idx_i<mtx_init_input.cols();idx_i++)
    {
        for (int idx_k=0;idx_k<size_k;idx_k++)
        {
            int kth_entry = idx_k+1;
            //if (kth_entry>entry_size) continue;
            std::complex<double> eta_k = mtx_H_dark(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-kth_entry);
            if (abs(eta_k)<eta_threshold) continue;
            mtx_S_vari(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-kth_entry) = mtx_C(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-kth_entry)/eta_k;
            for (int idx_n=0;idx_n<size_n;idx_n++)
            {
                int nth_entry = idx_n+1;
                //if (kth_entry>entry_size && nth_entry>entry_size) continue;
                //if (kth_entry>5 || nth_entry>5) continue;
                //if (kth_entry==nth_entry) continue;
                std::complex<double> eta_n = mtx_H_dark(mtx_init_input.cols()-nth_entry,mtx_init_input.cols()-nth_entry);
                if (abs(eta_n)<eta_threshold) continue;
                std::complex<double> coeff_l = mtx_C(mtx_init_input.cols()-kth_entry,mtx_init_input.cols()-nth_entry)/eta_n;
                mtx_l_vari(idx_i,mtx_init_input.cols()-kth_entry) += coeff_l*mtx_l_dark(idx_i,mtx_init_input.cols()-nth_entry);
                std::complex<double> coeff_r = -1.*mtx_C(mtx_init_input.cols()-nth_entry,mtx_init_input.cols()-kth_entry)/eta_n;
                mtx_r_vari(idx_i,mtx_init_input.cols()-kth_entry) += coeff_r*mtx_r_dark(idx_i,mtx_init_input.cols()-nth_entry);
            }
        }
    }

    double step_scale = step_frac;
    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    mtx_output = mtx_init_input;
    MatrixXcd mtx_vari = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    mtx_vari += step_scale*mtx_r_init*mtx_S_init*mtx_l_vari.transpose();
    mtx_vari += step_scale*mtx_r_vari*mtx_S_init*mtx_l_init.transpose();
    mtx_vari += step_scale*mtx_r_init*mtx_S_vari*mtx_l_init.transpose();
    //mtx_vari = GetTruncatedMatrix(mtx_vari, truncate_cutoff);
    mtx_output += mtx_vari;

    return mtx_output;

}
MatrixXcd SpectralDecompositionMethod_v3(MatrixXcd mtx_input, int entry_start, int entry_size, double step_frac, bool isBlind, int isLeft)
{

    //std::cout << "entry_size = " << entry_size << std::endl;

    //int effective_eqns = CountEffectiveEquations(mtx_input, 4);
    //std::cout << "effective_eqns = " << effective_eqns << std::endl;
    //int effective_eqns_r = CountEffectiveEquations(mtx_input, 1);
    //int effective_eqns_l = CountEffectiveEquations(mtx_input, 2);
    //std::cout << "effective_eqns_r = " << effective_eqns_r << std::endl;
    //std::cout << "effective_eqns_l = " << effective_eqns_l << std::endl;

    eigensolver_init = ComplexEigenSolver<MatrixXcd>(mtx_input);
    eigensolver_init_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_r_init = eigensolver_init.eigenvectors();
    MatrixXcd mtx_l_init = eigensolver_init_transpose.eigenvectors();
    MatrixXcd mtx_H_init = mtx_r_init.transpose()*mtx_l_init;
    MatrixXcd mtx_lambda_init(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_lambda_init(row,col) = 0.;
            if (row==col)
            {
                mtx_lambda_init(row,col) = eigensolver_init.eigenvalues()(col);
            }
        }
    }

    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    eigensolver_dark_transpose = ComplexEigenSolver<MatrixXcd>(mtx_dark.transpose());
    MatrixXcd mtx_r_dark = eigensolver_dark.eigenvectors();
    MatrixXcd mtx_l_dark = eigensolver_dark_transpose.eigenvectors();
    MatrixXcd mtx_H_dark = mtx_r_dark.transpose()*mtx_l_dark;
    MatrixXcd mtx_lambda_dark(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_lambda_dark(row,col) = 0.;
            if (row==col)
            {
                mtx_lambda_dark(row,col) = eigensolver_dark.eigenvalues()(col);
            }
        }
    }
    MatrixXcd mtx_S_dark = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    for (int col=0;col<mtx_input.cols();col++)
    {
        if (col<mtx_input.cols()-entry_size) continue;
        mtx_S_dark(col,col) = mtx_lambda_dark(col,col)/mtx_H_dark(col,col);
    }
    MatrixXcd mtx_S_init = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    for (int col=0;col<mtx_input.cols();col++)
    {
        if (col<mtx_input.cols()-entry_size) continue;
        mtx_S_init(col,col) = mtx_lambda_init(col,col)/mtx_H_dark(col,col);
    }

    MatrixXcd mtx_r_gamma = mtx_eigenvector_gamma;
    MatrixXcd mtx_l_gamma = mtx_eigenvector_inv_gamma.transpose();

    MatrixXcd mtx_r_data = mtx_eigenvector_data;
    MatrixXcd mtx_l_data = mtx_eigenvector_inv_data.transpose();
    MatrixXcd mtx_S_data = mtx_eigenvalue_data;
    MatrixXcd mtx_H_data = mtx_eigenvector_inv_data*mtx_eigenvector_data;

    int truncate_cutoff = entry_size;
    MatrixXcd mtx_data_truncated = mtx_data;
    MatrixXcd mtx_input_truncated = mtx_input;
    if (TruncateNoise)
    {
        //mtx_data_truncated = GetTruncatedMatrix(mtx_data, truncate_cutoff);
        mtx_input_truncated = GetTruncatedMatrix(mtx_input, truncate_cutoff);
    }
    MatrixXcd mtx_data_noise = mtx_data - GetTruncatedMatrix(mtx_data, truncate_cutoff);
    MatrixXcd mtx_input_noise = mtx_input - GetTruncatedMatrix(mtx_input, truncate_cutoff);

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

    int row_size_big = mtx_input.rows()*mtx_input.cols();
    int regularization_size = 0;
    //if (UseRegularization)
    //{
    //    //regularization_size += entry_size*mtx_input.cols();
    //    regularization_size += mtx_input.rows()*mtx_input.cols();
    //}
    VectorXcd vtr_vari_big = VectorXcd::Zero(2*entry_size*mtx_input.cols()+entry_size);
    VectorXcd vtr_Delta = VectorXcd::Zero(row_size_big+regularization_size);
    MatrixXcd mtx_Big = MatrixXcd::Zero(row_size_big+regularization_size,2*entry_size*mtx_input.cols()+entry_size);
    MatrixXcd mtx_Constraint = MatrixXcd::Zero(2*entry_size*mtx_input.cols(),2*entry_size*mtx_input.cols());
    VectorXcd vtr_Constraint_Delta = VectorXcd::Zero(2*entry_size*mtx_input.cols());
    for (int idx_k=0; idx_k<entry_size; idx_k++)
    {
        int nth_entry = idx_k + entry_start;
        for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
        {
            for (int idx_j=0; idx_j<mtx_input.cols(); idx_j++)
            {
                int idx_m = idx_j + mtx_input.cols()*idx_i;
                double weight = 1.;
                if (idx_i<binx_blind_global && idx_j<biny_blind_global)
                {
                    if (isBlind) weight = 0.; // blind gamma-ray region
                }
                vtr_Delta(idx_m) = weight*(mtx_data_truncated-mtx_input_truncated)(idx_i,idx_j);
                int idx_n = idx_j + mtx_input.cols()*idx_k;
                int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
                int idx_t = idx_k + mtx_input.cols()*entry_size + mtx_input.cols()*entry_size;

                mtx_Big(idx_m,idx_n) = weight*mtx_S_init(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry)*mtx_r_init(idx_i,mtx_input.rows()-nth_entry);
                mtx_Big(idx_m,idx_w) = weight*mtx_S_init(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry)*mtx_l_init(idx_j,mtx_input.rows()-nth_entry);
                mtx_Big(idx_m,idx_t) = weight*mtx_r_init(idx_i,mtx_input.rows()-nth_entry)*mtx_l_init(idx_j,mtx_input.rows()-nth_entry);
            }
        }
    }

    //if (UseRegularization)
    //{
    //    for (int idx_l=0; idx_l<mtx_input.rows(); idx_l++)
    //    {
    //        int nth_entry2 = idx_l + 1;
    //        for (int idx_m=0; idx_m<mtx_input.rows(); idx_m++)
    //        {
    //            int nth_entry3 = idx_m + 1;
    //            double weight_s1 = 1.*eigensolver_dark.eigenvalues()(mtx_dark.cols()-1).real();
    //            double weight_s2 = 1.*eigensolver_dark.eigenvalues()(mtx_dark.cols()-nth_entry2).real();
    //            double weight_s3 = 1.*eigensolver_dark.eigenvalues()(mtx_dark.cols()-nth_entry3).real();
    //            for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
    //            {
    //                int idx_u = idx_l + mtx_input.rows()*idx_m + row_size_big;
    //                vtr_Delta(idx_u) = 0.;
    //                if (nth_entry2!=nth_entry3 || (nth_entry2==entry_size && nth_entry3==entry_size))
    //                {
    //                    if (nth_entry2==entry_size && nth_entry3==entry_size)
    //                    {
    //                        int idx_t = idx_l + mtx_input.cols()*entry_size + mtx_input.cols()*entry_size;
    //                        mtx_Big(idx_u,idx_t) = weight_s1;
    //                    }
    //                    else if (nth_entry2<=entry_size)
    //                    {
    //                        int idx_n = idx_i + mtx_input.rows()*idx_l;
    //                        mtx_Big(idx_u,idx_n) = weight_s1*mtx_r_init(idx_i,mtx_input.rows()-nth_entry3);
    //                    }
    //                    else if (nth_entry3<=entry_size)
    //                    {
    //                        int idx_w = idx_i + mtx_input.rows()*idx_m + mtx_input.cols()*entry_size;
    //                        mtx_Big(idx_u,idx_w) = weight_s1*mtx_l_init(idx_i,mtx_input.rows()-nth_entry2);
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}

    //int coeff_size = mtx_input.rows();
    //int coeff_size = entry_size;
    int coeff_size = 5;
    //if (UseRegularization)
    //{
    //    coeff_size = mtx_input.rows();
    //}
    VectorXcd vtr_Coeff;
    MatrixXcd mtx_Basis;

    // orthogonal and normalization
    vtr_Coeff = VectorXcd::Zero(2*entry_size*coeff_size-entry_size*entry_size);
    mtx_Basis = MatrixXcd::Zero(2*entry_size*mtx_input.cols()+entry_size,2*entry_size*coeff_size-entry_size*entry_size);
    for (int idx_k=0; idx_k<entry_size; idx_k++)
    {
        int nth_entry = idx_k + 1;
        for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
        {
            int idx_n = idx_i + mtx_input.cols()*idx_k;
            int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
            int idx_t = idx_k + mtx_input.cols()*entry_size + mtx_input.cols()*entry_size;
            for (int idx_l=0; idx_l<entry_size; idx_l++)
            {
                int nth_entry2 = idx_l + 1;
                int idx_u = idx_l + entry_size*idx_k;
                int idx_v = idx_l*entry_size + idx_k;
                std::complex<double> eta_1 = mtx_H_dark(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry);
                std::complex<double> eta_2 = mtx_H_dark(mtx_input.rows()-nth_entry2,mtx_input.rows()-nth_entry2);
                //std::complex<double> eta_1 = mtx_H_init(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry);
                //std::complex<double> eta_2 = mtx_H_init(mtx_input.rows()-nth_entry2,mtx_input.rows()-nth_entry2);
                std::complex<double> lambda_1 = mtx_lambda_dark(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry);
                std::complex<double> lambda_2 = mtx_lambda_dark(mtx_input.rows()-nth_entry2,mtx_input.rows()-nth_entry2);
                if (nth_entry!=nth_entry2)
                {
                    //if (nth_entry==entry_size) continue;
                    //mtx_Basis(idx_n,idx_u) = 1./eta_2*1./(lambda_1-lambda_2)*mtx_l_dark(idx_i,mtx_input.rows()-nth_entry2);
                    //mtx_Basis(idx_w,idx_v) = -1./eta_2*1./(lambda_1-lambda_2)*mtx_r_dark(idx_i,mtx_input.rows()-nth_entry2);
                    mtx_Basis(idx_n,idx_u) = 1./eta_2*mtx_l_dark(idx_i,mtx_input.rows()-nth_entry2);
                    mtx_Basis(idx_w,idx_v) = -1./eta_2*mtx_r_dark(idx_i,mtx_input.rows()-nth_entry2);
                }
                else
                {
                    //mtx_Basis(idx_t,idx_u) = 1.;
                    if (nth_entry<3)
                    {
                        mtx_Basis(idx_t,idx_u) = 1.;
                    }
                }
            }
        }
    }
    for (int idx_k=0; idx_k<entry_size; idx_k++)
    {
        int nth_entry = idx_k + 1;
        for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
        {
            int idx_n = idx_i + mtx_input.cols()*idx_k;
            int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
            for (int idx_l=0; idx_l<coeff_size-entry_size; idx_l++)
            {
                int nth_entry2 = idx_l + entry_size + 1;
                int idx_u = idx_l + (coeff_size-entry_size)*idx_k + entry_size*entry_size;
                int idx_v = idx_l + (coeff_size-entry_size)*idx_k + (coeff_size-entry_size)*entry_size + entry_size*entry_size;
                std::complex<double> eta_1 = mtx_H_dark(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry);
                std::complex<double> eta_2 = mtx_H_dark(mtx_input.rows()-nth_entry2,mtx_input.rows()-nth_entry2);
                //std::complex<double> eta_1 = mtx_H_init(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry);
                //std::complex<double> eta_2 = mtx_H_init(mtx_input.rows()-nth_entry2,mtx_input.rows()-nth_entry2);
                std::complex<double> lambda_1 = mtx_lambda_dark(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry);
                std::complex<double> lambda_2 = mtx_lambda_dark(mtx_input.rows()-nth_entry2,mtx_input.rows()-nth_entry2);
                double lambda_threshold = 1.0*mtx_lambda_dark(mtx_input.rows()-entry_size,mtx_input.rows()-entry_size).real();
                if (abs(eta_2)<0.01) continue;
                if (isLeft==0)
                {
                    //mtx_Basis(idx_n,idx_u) = 1./eta_2*1./(lambda_1-lambda_2)*mtx_l_dark(idx_i,mtx_input.rows()-nth_entry2); 
                    //mtx_Basis(idx_w,idx_v) = 1./eta_2*1./(lambda_1-lambda_2)*mtx_r_dark(idx_i,mtx_input.rows()-nth_entry2); 
                    mtx_Basis(idx_n,idx_u) = mtx_l_dark(idx_i,mtx_input.rows()-nth_entry2); 
                    mtx_Basis(idx_w,idx_v) = mtx_r_dark(idx_i,mtx_input.rows()-nth_entry2); 
                }
                else if (isLeft==1)
                {
                    //mtx_Basis(idx_n,idx_u) = 1./eta_2*1./(lambda_1-lambda_2)*mtx_l_dark(idx_i,mtx_input.rows()-nth_entry2); 
                    mtx_Basis(idx_n,idx_u) = mtx_l_dark(idx_i,mtx_input.rows()-nth_entry2); 
                }
                else if (isLeft==2)
                {
                    //mtx_Basis(idx_w,idx_v) = 1./eta_2*1./(lambda_1-lambda_2)*mtx_r_dark(idx_i,mtx_input.rows()-nth_entry2); 
                    mtx_Basis(idx_w,idx_v) = mtx_r_dark(idx_i,mtx_input.rows()-nth_entry2); 
                }
            }
        }
    }

    //for (int idx_k=0; idx_k<entry_size; idx_k++)
    //{
    //    int nth_entry = idx_k + 1;
    //    for (int idx_l=0; idx_l<mtx_input.cols(); idx_l++)
    //    {
    //        int nth_entry2 = idx_l + 1;
    //        double weight_s = 1000.*eigensolver_dark.eigenvalues()(mtx_dark.cols()-1).real();
    //        for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
    //        {
    //            int idx_u = 2*idx_l + idx_k*2*mtx_input.cols();
    //            int idx_n = idx_i + mtx_input.cols()*idx_k;
    //            int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
    //            vtr_Constraint_Delta(idx_u+0) = 0.;
    //            vtr_Constraint_Delta(idx_u+1) = 0.;
    //            if (nth_entry2!=nth_entry)
    //            {
    //                mtx_Constraint(idx_u+0,idx_w) = weight_s*mtx_l_init(idx_i,mtx_input.rows()-nth_entry2);
    //                mtx_Constraint(idx_u+1,idx_n) = weight_s*mtx_r_init(idx_i,mtx_input.rows()-nth_entry2);
    //            }
    //        }
    //    }
    //}

    //int entry_threshold = NumberOfEigenvectors_Stable;
    //if (entry_threshold==0) entry_threshold = 1;
    //entry_threshold = 1;
    //double lambda_noise = eigensolver_dark.eigenvalues()(mtx_dark.cols()-entry_threshold).real();
    //svd_threshold = 1./abs(svd_threshold_scale*lambda_noise);
    //std::cout << "svd_threshold = " << svd_threshold << std::endl;

    //BDCSVD<MatrixXcd> svd(mtx_Big*mtx_Basis, ComputeThinU | ComputeThinV);
    //svd.setThreshold(svd_threshold);  // size of singular value to be considered as nonzero.
    //vtr_Coeff = svd.solve(vtr_Delta);
    //vtr_vari_big = mtx_Basis*vtr_Coeff;

    //vtr_Coeff = SolutionWithConstraints(mtx_Big*mtx_Basis, mtx_Constraint*mtx_Basis, vtr_Delta, vtr_Constraint_Delta);
    //vtr_vari_big = mtx_Basis*vtr_Coeff;


    MatrixXcd mtx_l_final = mtx_l_init;
    MatrixXcd mtx_r_final = mtx_r_init;
    MatrixXcd mtx_l_vari = MatrixXcd::Zero(mtx_l_init.rows(),mtx_l_init.cols());
    MatrixXcd mtx_r_vari = MatrixXcd::Zero(mtx_r_init.rows(),mtx_r_init.cols());
    MatrixXcd mtx_S_vari = MatrixXcd::Zero(mtx_r_init.rows(),mtx_r_init.cols());

    MatrixXcd mtx_H_vari_true = (mtx_l_data-mtx_l_init).transpose()*mtx_r_init + mtx_l_init.transpose()*(mtx_r_data-mtx_r_init);

    MatrixXcd mtx_H_vari = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());

    bool isGood = false;
    while (!isGood)
    {
        isGood = true;
        BDCSVD<MatrixXcd> svd(mtx_Big*mtx_Basis, ComputeThinU | ComputeThinV);
        svd.setThreshold(svd_threshold);  // size of singular value to be considered as nonzero.
        vtr_Coeff = svd.solve(vtr_Delta);
        //MatrixXcd mtx_A = mtx_Big*mtx_Basis;
        //JacobiSVD<MatrixXcd> svd(mtx_A, ComputeFullU | ComputeFullV);
        //MatrixXcd mtx_U = svd.matrixU();
        //MatrixXcd mtx_V = svd.matrixV();
        //MatrixXcd mtx_S = MatrixXcd::Zero(mtx_A.rows(),mtx_A.cols());
        //MatrixXcd mtx_S_inv = MatrixXcd::Zero(mtx_A.cols(),mtx_A.rows());
        //std::cout << "svd.singularValues():" << std::endl;
        //std::cout << svd.singularValues() << std::endl;
        //for (int entry=0;entry<svd.singularValues().size();entry++)
        //{
        //    if (abs(svd.singularValues()(entry))<1e-10) continue;
        //    mtx_S(entry,entry) = svd.singularValues()(entry);
        //    //mtx_S_inv(entry,entry) = 1./(svd.singularValues()(entry)+0.1*svd.singularValues()(0));
        //    mtx_S_inv(entry,entry) = 1./(svd.singularValues()(entry));
        //}
        //vtr_Coeff = mtx_V*mtx_S_inv*mtx_U.conjugate().transpose()*vtr_Delta;
        vtr_vari_big = mtx_Basis*vtr_Coeff;
        for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
        {
            for (int idx_j=0; idx_j<mtx_input.cols(); idx_j++)
            {
                for (int idx_k=0; idx_k<entry_size; idx_k++)
                {
                    int idx_m = idx_j + mtx_input.cols()*idx_i;
                    int idx_n = idx_j + mtx_input.cols()*idx_k;
                    int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
                    int idx_t = idx_k + mtx_input.cols()*entry_size + mtx_input.cols()*entry_size;
                    int nth_entry = idx_k + entry_start;
                    mtx_l_vari(idx_j,mtx_l_vari.cols()-nth_entry) = vtr_vari_big(idx_n);
                    mtx_r_vari(idx_i,mtx_r_vari.cols()-nth_entry) = vtr_vari_big(idx_w);
                    mtx_S_vari(mtx_S_vari.cols()-nth_entry,mtx_S_vari.cols()-nth_entry) = vtr_vari_big(idx_t);
                }
            }
        }
        mtx_H_vari = mtx_l_vari.transpose()*mtx_r_init + mtx_l_init.transpose()*mtx_r_vari;
        //if (UseRegularization)
        //{
        //    for (int idx_k=0; idx_k<entry_size; idx_k++)
        //    {
        //        int nth_entry = idx_k + 1;
        //        for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
        //        {
        //            int idx_n = idx_i + mtx_input.cols()*idx_k;
        //            int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
        //            for (int idx_l=0; idx_l<coeff_size-entry_size; idx_l++)
        //            {
        //                int nth_entry2 = idx_l + entry_size + 1;
        //                int idx_u = idx_l + (coeff_size-entry_size)*idx_k + entry_size*entry_size;
        //                int idx_v = idx_l + (coeff_size-entry_size)*idx_k + (coeff_size-entry_size)*entry_size + entry_size*entry_size;
        //                if (abs(mtx_H_vari(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry2))>0.005)
        //                {
        //                    if (mtx_Basis(idx_n,idx_u)!=0.)
        //                    {
        //                        isGood = false;
        //                        mtx_Basis(idx_n,idx_u) = 0.;
        //                    }
        //                }
        //                if (abs(mtx_H_vari(mtx_input.rows()-nth_entry2,mtx_input.rows()-nth_entry))>0.005)
        //                {
        //                    if (mtx_Basis(idx_w,idx_v)!=0.)
        //                    {
        //                        isGood = false;
        //                        mtx_Basis(idx_w,idx_v) = 0.;
        //                    }
        //                }
        //            }
        //        }
        //    }
        //}
        //else
        //{
        //    isGood = true;
        //}
    }

    //int which_entry =  max(1,entry_size-1);
    //double largest_S_init = mtx_S_init(mtx_S_init.rows()-which_entry,mtx_S_init.rows()-which_entry).real();
    //double largest_S_vari = mtx_S_vari(mtx_S_init.rows()-which_entry,mtx_S_init.rows()-which_entry).real();
    //double step_size = abs(largest_S_init*step_frac);
    //double step_scale = abs(step_size/largest_S_vari);
    //step_scale = min(1.,step_scale);
    double step_scale = step_frac;

    if (isBlind)
    {
        //std::cout << "mtx_S:" << std::endl;
        //std::cout << mtx_S.block(mtx_S.rows()-4,mtx_S.cols()-4,4,4) << std::endl;

        //std::cout << "mtx_S_vari:" << std::endl;
        //std::cout << mtx_S_vari.block(mtx_S_vari.rows()-5,mtx_S_vari.cols()-5,5,5).cwiseAbs() << std::endl;
        for (int i=1;i<=5;i++)
        {
            std::cout << step_scale*mtx_S_vari(mtx_S_vari.rows()-i,mtx_S_vari.rows()-i).real() << ", " ;
        }
        std::cout << std::endl;

        //std::cout << "mtx_H_data:" << std::endl;
        //std::cout << mtx_H_data.block(mtx_H_init.rows()-5,mtx_H_init.cols()-5,5,5).cwiseAbs() << std::endl;
        //std::cout << mtx_H_data.cwiseAbs() << std::endl;

        //std::cout << "mtx_H_vari:" << std::endl;
        //std::cout << mtx_H_vari.block(mtx_H_vari.rows()-5,mtx_H_vari.cols()-5,5,5).cwiseAbs() << std::endl;

        //std::cout << "mtx_H_vari_true:" << std::endl;
        //std::cout << mtx_H_vari_true.block(mtx_H_vari_true.rows()-5,mtx_H_vari_true.cols()-5,5,5).cwiseAbs() << std::endl;

        //MatrixXcd mtx_R = mtx_r_init.transpose()*mtx_r_init;
        //std::cout << "mtx_R:" << std::endl;
        //std::cout << mtx_R.cwiseAbs() << std::endl;
    }

    mtx_output = mtx_input;
    MatrixXcd mtx_delta = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    mtx_delta += step_scale*mtx_r_init*mtx_S_init*mtx_l_vari.transpose();
    mtx_delta += step_scale*mtx_r_vari*mtx_S_init*mtx_l_init.transpose();
    mtx_delta += step_scale*mtx_r_init*mtx_S_vari*mtx_l_init.transpose();
    //mtx_delta = GetTruncatedMatrix(mtx_delta, truncate_cutoff);
    mtx_output += mtx_delta;

    mtx_eigenvector_init += mtx_r_vari;
    mtx_eigenvector_inv_init += mtx_l_vari.transpose();
    mtx_eigenvector_vari = mtx_r_vari;
    mtx_eigenvector_inv_vari = mtx_l_vari.transpose();

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
    double current_chi2 = GetChi2Function(mtx_input,mtx_data,isBlind);
    //if (current_chi2>2.*init_chi2) 
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
double CountGammaRegion(MatrixXcd mtx_input)
{
    double count = 0.;
    for (int col=0;col<mtx_input.cols();col++)
    {
        for (int row=0;row<mtx_input.rows();row++)
        {
            if (row>=binx_blind_global) continue;
            if (col>=biny_blind_global) continue;
            count += mtx_input(row,col).real();
        }
    }
    return count;
}
double CountCosmicRayRegion(MatrixXcd mtx_input)
{
    double count = 0.;
    for (int col=0;col<mtx_input.cols();col++)
    {
        for (int row=0;row<mtx_input.rows();row++)
        {
            if (row<binx_blind_global && col<biny_blind_global) continue;
            count += mtx_input(row,col).real();
        }
    }
    return count;
}
double Chi2CosmicRayRegion(MatrixXcd mtx_input, MatrixXcd mtx_model)
{
    double count = 0.;
    for (int col=0;col<mtx_input.cols();col++)
    {
        for (int row=0;row<mtx_input.rows();row++)
        {
            count += mtx_input(row,col).real();
        }
    }
    double chi2 = 0.;
    for (int col=0;col<mtx_input.cols();col++)
    {
        for (int row=0;row<mtx_input.rows();row++)
        {
            if (row<binx_blind_global && col<biny_blind_global) continue;
            chi2 += pow((mtx_input(row,col)-mtx_model(row,col)).real()/count,2);
        }
    }
    return chi2;
}
double SwapCRRegion(MatrixXcd mtx_input, MatrixXcd& mtx_output)
{
    for (int col=0;col<mtx_input.cols();col++)
    {
        for (int row=0;row<mtx_input.rows();row++)
        {
            if (col>=binx_blind_global || row>=biny_blind_global)
            {
                mtx_output(row,col) = mtx_input(row,col);
            }
        }
    }
}
void FindPerturbationSolution(MatrixXcd mtx_ref, MatrixXcd mtx_input_fix, MatrixXcd& mtx_input, MatrixXcd mtx_H_init, MatrixXcd mtx_r_init, MatrixXcd mtx_lambda_init, MatrixXcd mtx_l_init, int entry_size, bool doPrint)
{

    double cutoff = 0.001;

    MatrixXcd mtx_bkgd = mtx_input;

    MatrixXcd mtx_delta = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());

    mtx_delta = mtx_ref-mtx_bkgd;
    //mtx_delta = mtx_bkgd-mtx_input_fix;
    //mtx_delta.block(0,0,binx_blind_global,biny_blind_global) = (mtx_ref-mtx_bkgd).block(0,0,binx_blind_global,biny_blind_global);
    //mtx_delta.block(0,0,binx_blind_global,biny_blind_global) = (mtx_bkgd-mtx_input_fix).block(0,0,binx_blind_global,biny_blind_global);
    //mtx_delta.block(0,0,binx_blind_global,biny_blind_global) = (mtx_bkgd-1.1*mtx_input_fix).block(0,0,binx_blind_global,biny_blind_global);

    MatrixXcd mtx_coeff = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    MatrixXcd mtx_proj = mtx_l_init.transpose()*mtx_delta*mtx_r_init;
    for (int col=0;col<mtx_input.cols();col++)
    {
        for (int row=0;row<mtx_input.rows();row++)
        {
            if (row<mtx_input.cols()-entry_size && col<mtx_input.cols()-entry_size) continue;
            if (row==col) continue;
            //if (abs(mtx_lambda_init(row,row)-mtx_lambda_init(col,col))<cutoff) continue;
            mtx_coeff(row,col) = 1.*mtx_proj(row,col)/(mtx_lambda_init(row,row)-mtx_lambda_init(col,col));
        }
    }

    MatrixXcd mtx_r_next = mtx_r_init;
    MatrixXcd mtx_l_next = mtx_l_init;
    MatrixXcd mtx_lambda_next = mtx_lambda_init;
    for (int idx_j=0;idx_j<mtx_input.cols();idx_j++)
    {
        if (idx_j<mtx_input.cols()-entry_size) continue;
        mtx_lambda_next(idx_j,idx_j) += mtx_proj(idx_j,idx_j)/mtx_H_init(idx_j,idx_j);
        for (int idx_i=0;idx_i<mtx_input.rows();idx_i++)
        {
            for (int idx_k=0;idx_k<mtx_input.rows();idx_k++)
            {
                if (idx_k<mtx_input.cols()-entry_size) continue;
                mtx_r_next(idx_i,idx_j) += -1.*mtx_coeff(idx_k,idx_j)*mtx_r_init(idx_i,idx_k)/mtx_H_init(idx_k,idx_k);
                mtx_l_next(idx_i,idx_j) += mtx_coeff(idx_j,idx_k)*mtx_l_init(idx_i,idx_k)/mtx_H_init(idx_k,idx_k);
            }
        }
    }

    MatrixXcd mtx_H_next = mtx_r_next.transpose()*mtx_l_next;
    if (doPrint)
    {
        std::cout << "mtx_coeff:" << std::endl;
        //std::cout << mtx_coeff.block(mtx_input.rows()-6,mtx_input.cols()-6,6,6).cwiseAbs() << std::endl;
        std::cout << mtx_coeff.cwiseAbs() << std::endl;
        std::cout << "mtx_lambda_init:" << std::endl;
        std::cout << mtx_lambda_init.block(mtx_input.rows()-6,mtx_input.cols()-6,6,6).cwiseAbs() << std::endl;
        std::cout << "mtx_lambda_next:" << std::endl;
        std::cout << mtx_lambda_next.block(mtx_input.rows()-6,mtx_input.cols()-6,6,6).cwiseAbs() << std::endl;
        std::cout << "mtx_H_init:" << std::endl;
        std::cout << mtx_H_init.block(mtx_input.rows()-6,mtx_input.cols()-6,6,6).cwiseAbs() << std::endl;
        std::cout << "mtx_H_next:" << std::endl;
        std::cout << mtx_H_next.block(mtx_input.rows()-6,mtx_input.cols()-6,6,6).cwiseAbs() << std::endl;
    }

    MatrixXcd mtx_S_next = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    for (int col=0;col<mtx_input.cols();col++)
    {
        if (col<mtx_input.cols()-entry_size) continue;
        //if (abs(mtx_H_init(col,col))<cutoff) continue;
        mtx_S_next(col,col) = mtx_lambda_next(col,col)/mtx_H_init(col,col);
    }

    mtx_input = mtx_r_next*mtx_S_next*mtx_l_next.transpose();

}
MatrixXcd PerturbationMethod(MatrixXcd mtx_ref, MatrixXcd mtx_input, bool doSwap)
{

    int entry_size = DetermineStableNumberOfEigenvalues(mtx_input);
    entry_size = min(entry_size,3);
    //entry_size = min(entry_size,2);

    MatrixXcd mtx_input_swap = 1.0*mtx_input;
    if (doSwap)
    {
        SwapCRRegion(mtx_ref, mtx_input_swap);
    }

    MatrixXcd mtx_bkgd = mtx_input_swap;
    eigensolver_init = ComplexEigenSolver<MatrixXcd>(mtx_bkgd);
    eigensolver_init_transpose = ComplexEigenSolver<MatrixXcd>(mtx_bkgd.transpose());
    MatrixXcd mtx_r_init = eigensolver_init.eigenvectors();
    MatrixXcd mtx_l_init = eigensolver_init_transpose.eigenvectors();
    MatrixXcd mtx_H_init = mtx_r_init.transpose()*mtx_l_init;
    MatrixXcd mtx_lambda_init(mtx_bkgd.rows(),mtx_bkgd.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_lambda_init(row,col) = 0.;
            if (row==col)
            {
                mtx_lambda_init(row,col) = eigensolver_init.eigenvalues()(col);
            }
        }
    }
    MatrixXcd mtx_S_init = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    for (int col=0;col<mtx_input.cols();col++)
    {
        if (col<mtx_input.cols()-entry_size) continue;
        mtx_S_init(col,col) = mtx_lambda_init(col,col)/mtx_H_init(col,col);
    }
    mtx_bkgd = mtx_r_init*mtx_S_init*mtx_l_init.transpose();

    FindPerturbationSolution(mtx_ref, mtx_input_swap, mtx_bkgd, mtx_H_init, mtx_r_init, mtx_lambda_init, mtx_l_init, entry_size, false);
    //for (int iter=0;iter<20;iter++)
    //{
    //    bool doPrint = false;
    //    //if (iter<=3) doPrint = true;
    //    //std::cout << "=========== " << iter+2 << "-th iteration ===========" << std::endl;
    //    FindPerturbationSolution(mtx_ref, mtx_input_swap, mtx_bkgd, mtx_H_init, mtx_r_init, mtx_lambda_init, mtx_l_init, entry_size, doPrint);
    //}
    //std::cout << "=========== final iteration ===========" << std::endl;
    //FindPerturbationSolution(mtx_ref, mtx_input_swap, mtx_bkgd, mtx_H_init, mtx_r_init, mtx_lambda_init, mtx_l_init, entry_size, true);

    return mtx_bkgd;

}
void LeastSquareSolutionMethod(bool DoSequential, bool isBlind, TH1D* Hist_Converge, TH1D* Hist_Optimization, TH1D* Hist_CosmicRayChi2, double alpha)
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

    std::cout << "initial chi2 in CR = " << GetChi2Function(mtx_dark,mtx_data,isBlind) << std::endl;
    double init_chi2 = GetChi2Function(mtx_dark,mtx_data,isBlind);

    MatrixXcd mtx_temp = mtx_dark;
    double eigenvalue_dark_real = 0.;
    double eigenvalue_dark_imag = 0.;
    double eigenvalue_data_real = 0.;
    double eigenvalue_data_imag = 0.;
    double step_frac= 1.0;
    double imag_real_ratio = 1./100.;

    std::cout << "NumberOfEigenvectors_Stable = " << NumberOfEigenvectors_Stable << std::endl;
    if (DoSequential)
    {
        for (int rank=1; rank<=NumberOfEigenvectors_Stable; rank++)
        {
            mtx_temp = SpectralDecompositionMethod_v3(mtx_data_bkgd, rank, 1, 1.0, isBlind, 0);
            if (CheckIfEigenvalueMakeSense(mtx_temp, init_chi2, rank, isBlind))
            {
                mtx_data_bkgd = mtx_temp;
            }
        }
    }
    else
    {
        if (NumberOfEigenvectors_Stable<=0)
        {
            return;
        }
        double step_frac = 0.1;
        int entry_size = NumberOfEigenvectors_Stable;
        double data_count = CountGammaRegion(mtx_data);
        double data_cr_count = CountCosmicRayRegion(mtx_data);
        double bkgd_count = CountGammaRegion(mtx_data_bkgd);
        Hist_Converge->SetBinContent(1,1.-bkgd_count/data_count);
        data_cr_count = 1.;
        MatrixXcd mtx_data_rescale = mtx_data/data_cr_count;
        MatrixXcd mtx_dark_rescale = mtx_dark/data_cr_count;
        MatrixXcd mtx_init_rescale = mtx_data_bkgd/data_cr_count;
        for (int binx=1; binx<=Hist_Optimization->GetNbinsX(); binx++)
        {
            double alpha_temp = pow(10.,Hist_Optimization->GetBinCenter(binx));
            if (RegularizationType==3 || RegularizationType==5)
            {
                alpha_temp = double(binx);
            }
            if (EigenDecomposition)
            {
                mtx_temp = LowRankOptimizationMethod(mtx_init_rescale, mtx_data_rescale, mtx_dark_rescale, 1, entry_size, 1.0, true, 0, alpha_temp);
            }
            else
            {
                mtx_temp = NuclearNormMinimization(mtx_init_rescale, mtx_data_rescale, mtx_dark_rescale, 1, entry_size, 1.0, true, 0, alpha_temp);
            }
            mtx_temp = mtx_temp*data_cr_count;
            if (!CheckIfEigenvalueMakeSense(mtx_temp, init_chi2, entry_size, isBlind))
            {
                mtx_temp = mtx_dark;
            }
            bkgd_count = CountGammaRegion(mtx_temp);
            Hist_Optimization->SetBinContent(binx,abs(1.-bkgd_count/data_count));
            double model_chi2 = Chi2CosmicRayRegion(mtx_data, mtx_temp);
            Hist_CosmicRayChi2->SetBinContent(binx,model_chi2);
        }
        if (RegularizationType==3 || RegularizationType==4 || RegularizationType==5)
        {
            int min_bin = 1;
            double min_val = 1e6;
            for (int binx=1; binx<=Hist_Optimization->GetNbinsX(); binx++)
            {
                if (min_val>Hist_Optimization->GetBinContent(binx))
                {
                    min_val = Hist_Optimization->GetBinContent(binx);
                    min_bin = binx;
                }
            }
            std::cout << "min_val = " << min_val << std::endl;
            std::cout << "min_bin = " << min_bin << std::endl;
            std::cout << "Hist_Optimization->GetBinCenter(min_bin) = " << Hist_Optimization->GetBinCenter(min_bin) << std::endl;
            alpha = double(min_bin);
            if (RegularizationType==4)
            {
                alpha = pow(10.,Hist_Optimization->GetBinCenter(min_bin));
            }
        }
        //mtx_temp = SpectralDecompositionMethod_v3(mtx_data_bkgd, 1, entry_size, 1.0, true, 0);
        if (EigenDecomposition)
        {
            mtx_temp = LowRankOptimizationMethod(mtx_init_rescale, mtx_data_rescale, mtx_dark_rescale, 1, entry_size, 1.0, true, 0, alpha);
        }
        else
        {
            mtx_temp = NuclearNormMinimization(mtx_init_rescale, mtx_data_rescale, mtx_dark_rescale, 1, entry_size, 1.0, true, 0, alpha);
        }
        mtx_temp = mtx_temp*data_cr_count;
        if (CheckIfEigenvalueMakeSense(mtx_temp, init_chi2, entry_size, isBlind))
        {
            mtx_data_bkgd = mtx_temp;
            bkgd_count = CountGammaRegion(mtx_data_bkgd);
        }

        return;
    }

}
void GetReducedEigenvalueMatrix(int rank_cutoff)
{
    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
    JacobiSVD<MatrixXd> svd_data(mtx_data.real(), ComputeFullU | ComputeFullV);
    JacobiSVD<MatrixXd> svd_dark(mtx_dark.real(), ComputeFullU | ComputeFullV);
    for (int i=0;i<mtx_data.cols();i++)
    {
        for (int j=0;j<mtx_data.rows();j++)
        {
            if (EigenDecomposition)
            {
                if (i==j && i==mtx_data.cols()-rank_cutoff-1) 
                {
                    mtx_eigenval_data_redu(i,j) = eigensolver_data.eigenvalues()(i);
                    mtx_eigenval_dark_redu(i,j) = eigensolver_dark.eigenvalues()(i);
                }
                else
                {
                    mtx_eigenval_data_redu(i,j) = 0;
                    mtx_eigenval_dark_redu(i,j) = 0;
                }
            }
            else
            {
                if (i==j && i==rank_cutoff) 
                {
                    mtx_eigenval_data_redu(i,j) = svd_data.singularValues()(i);
                    mtx_eigenval_dark_redu(i,j) = svd_dark.singularValues()(i);
                }
                else
                {
                    mtx_eigenval_data_redu(i,j) = 0;
                    mtx_eigenval_dark_redu(i,j) = 0;
                }
            }
        }
    }
}
void GetCRReplacedHistogram(TH2D* hist_data, TH2D* hist_dark)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    TH2D hist_temp = TH2D("hist_temp","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            if (bx<binx_blind && by<biny_blind)
            //if (bx<binx_blind || by<biny_blind)
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
MatrixXcd GetCRReplacedMatrix(MatrixXcd mtx_ref, MatrixXcd mtx_input)
{
    MatrixXcd mtx_output = mtx_input;
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            if (col<binx_blind_global && row<biny_blind_global)
            {
                mtx_output(row,col) = mtx_input(row,col);
            }
            else
            {
                mtx_output(row,col) = mtx_ref(row,col);
            }
        }
    }
    return mtx_output;
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

    GetTruncatedHistogram(hist_data,&hist_data_temp,NumberOfEigenvectors_Stable);
    GetTruncatedHistogram(hist_dark,&hist_dark_temp,NumberOfEigenvectors_Stable);
    hist_data_noise.Add(hist_data);
    hist_data_noise.Add(&hist_data_temp,-1.);
    hist_dark->Reset();
    hist_dark->Add(&hist_dark_temp);
    hist_dark->Add(&hist_data_noise);
}
void NormalizeDarkMatrix(TH2D* hist_data, TH2D* hist_dark)
{
    int binx_lower = hist_data->GetXaxis()->FindBin(MSCL_cut_lower);
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int binx_upper = hist_data->GetNbinsX();
    int biny_lower = hist_data->GetYaxis()->FindBin(MSCW_cut_lower);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    int biny_upper = hist_data->GetNbinsY();
    double Data_CR_Integral = 0.;
    double Dark_CR_Integral = 0.;
    for (int binx=1;binx<=hist_data->GetNbinsX();binx++)
    {
        for (int biny=1;biny<=hist_data->GetNbinsY();biny++)
        {
            if (binx<binx_blind && biny<biny_blind) continue;
            //if (binx>binx_blind+hist_data->GetNbinsX()/4) continue;
            //if (biny>biny_blind+hist_data->GetNbinsY()/4) continue;
            //if (binx>=binx_blind) continue;
            Data_CR_Integral += hist_data->GetBinContent(binx,biny);
            Dark_CR_Integral += hist_dark->GetBinContent(binx,biny);
        }
    }
    if (Dark_CR_Integral!=0)
    {
        double dark_scale = Data_CR_Integral/Dark_CR_Integral;
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
        GetCRReplacedHistogram(hist_data,hist_dark_alter);
    }

}
//MatrixXcd MinChi2Method(TH1D* Hist_Chi2, TH2D* Hist_Data, TH2D* Hist_Dark)
//{
//
//    mtx_dark = fillMatrix(Hist_Dark);
//    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
//    mtx_data_bkgd = fillMatrix(Hist_Dark);
//    MatrixXcd mtx_best_bkgd = mtx_dark;
//
//    TH2D Hist_Temp_Dark = TH2D("Hist_Temp_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
//    TH2D Hist_Temp_Dark_Alt = TH2D("Hist_Temp_Dark_Alt","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
//
//    TH2D Hist_Temp2_Gamma = TH2D("Hist_Temp2_Gamma","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
//    fill2DHistogram(&Hist_Temp2_Gamma,mtx_gamma_raw);
//    std::cout << "Hist_Temp2_Gamma.Integral() = " << Hist_Temp2_Gamma.Integral() << std::endl;
//    TH2D Hist_Temp_Data = TH2D("Hist_Temp_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
//    Hist_Temp_Data.Add(Hist_Data);
//    Hist_Temp_Dark.Add(Hist_Dark);
//
//    int binx_lower = Hist_Temp_Data.GetXaxis()->FindBin(MSCL_cut_lower);
//    int binx_blind = Hist_Temp_Data.GetXaxis()->FindBin(MSCL_cut_blind)-1;
//    int binx_upper = Hist_Temp_Data.GetXaxis()->FindBin(1.)-1;
//    int biny_lower = Hist_Temp_Data.GetYaxis()->FindBin(MSCW_cut_lower);
//    int biny_blind = Hist_Temp_Data.GetYaxis()->FindBin(MSCW_cut_blind)-1;
//    int biny_upper = Hist_Temp_Data.GetYaxis()->FindBin(1.)-1;
//
//    TH2D Hist_Temp_Bkgd = TH2D("Hist_Temp_Bkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
//    double Data_SR_Integral = Hist_Temp_Data.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
//    double Data_Integral = Hist_Temp_Data.Integral();
//    double Dark_Integral = Hist_Temp_Dark.Integral();
//    double Gamma_SR_Integral = Hist_Temp2_Gamma.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
//
//    if (Data_Integral==0.) return mtx_best_bkgd;
//    if (Dark_Integral==0.) return mtx_best_bkgd;
//    if (Gamma_SR_Integral==0.) return mtx_best_bkgd;
//
//    Hist_Temp_Bkgd.Reset();
//    double Dark_SR_Integral = Hist_Temp_Dark.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
//
//    double chi2_best = 1e20;
//    for (int bin=0;bin<Hist_Chi2->GetNbinsX();bin++)
//    {
//        Hist_Temp_Dark_Alt.Reset();
//        Hist_Temp_Data.Reset();
//        Hist_Temp_Data.Add(Hist_Data);
//        Hist_Temp_Dark.Reset();
//        Hist_Temp_Dark.Add(Hist_Dark);
//        Hist_Temp_Bkgd.Reset();
//        Hist_Temp2_Gamma.Reset();
//        fill2DHistogram(&Hist_Temp2_Gamma,mtx_gamma_raw);
//        double scale = Hist_Chi2->GetBinCenter(bin+1);
//        double delta = Dark_SR_Integral*scale;
//        Hist_Temp2_Gamma.Scale((delta)/Gamma_SR_Integral);
//        Hist_Temp_Data.Add(&Hist_Temp2_Gamma,-1.);
//        mtx_gamma = fillMatrix(&Hist_Temp2_Gamma);
//        mtx_data = fillMatrix(&Hist_Temp_Data);
//        NormalizeDarkMatrix(&Hist_Temp_Data,&Hist_Temp_Dark);
//        AlterDarkMatrix(Hist_Data, &Hist_Temp_Dark, &Hist_Temp_Dark_Alt);
//        mtx_dark = fillMatrix(&Hist_Temp_Dark_Alt);
//        eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
//        LeastSquareSolutionMethod(false, true);
//        //fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
//        //GetNoiseReplacedMatrix(&Hist_Temp_Data,&Hist_Temp_Bkgd,false);
//        //mtx_data_bkgd = fillMatrix(&Hist_Temp_Bkgd);
//        double chi2 = GetChi2Function(mtx_data_bkgd,mtx_data,true);
//        if (chi2<chi2_best)
//        {
//            mtx_best_bkgd = mtx_data_bkgd;
//        }
//        Hist_Chi2->SetBinContent(bin+1,chi2);
//    }
//    int min_error_bin = Hist_Chi2->GetMinimumBin();
//    Hist_Temp_Dark_Alt.Reset();
//    Hist_Temp_Data.Reset();
//    Hist_Temp_Data.Add(Hist_Data);
//    Hist_Temp_Dark.Reset();
//    Hist_Temp_Dark.Add(Hist_Dark);
//    Hist_Temp_Bkgd.Reset();
//    Hist_Temp2_Gamma.Reset();
//    fill2DHistogram(&Hist_Temp2_Gamma,mtx_gamma_raw);
//    double scale = Hist_Chi2->GetBinCenter(min_error_bin);
//    double delta = Dark_SR_Integral*scale;
//    Hist_Temp2_Gamma.Scale((delta)/Gamma_SR_Integral);
//    Hist_Temp_Data.Add(&Hist_Temp2_Gamma,-1.);
//    mtx_gamma = fillMatrix(&Hist_Temp2_Gamma);
//    mtx_data = fillMatrix(&Hist_Temp_Data);
//    NormalizeDarkMatrix(&Hist_Temp_Data,&Hist_Temp_Dark);
//    AlterDarkMatrix(Hist_Data, &Hist_Temp_Dark, &Hist_Temp_Dark_Alt);
//    mtx_dark = fillMatrix(&Hist_Temp_Dark_Alt);
//    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
//    LeastSquareSolutionMethod(false, true);
//    fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
//
//    return mtx_data_bkgd;
//
//}
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
    //if (TString(target).Contains("Crab"))
    //{
    //    if (source_theta2_cut==0.)
    //    {
    //        MSCW_cut_blind = MSCW_cut_loose;
    //        MSCL_cut_blind = MSCL_cut_loose;
    //    }
    //}
    //if (TString(target).Contains("Mrk421"))
    //{
    //    if (source_theta2_cut==0.)
    //    {
    //        MSCW_cut_blind = MSCW_cut_loose;
    //        MSCL_cut_blind = MSCL_cut_loose;
    //    }
    //}

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
    TH2D Hist_Data_ElevNSB = TH2D("Hist_Data_ElevNSB","",20,0,10,90,0,90);
    TH2D Hist_Dark_ElevNSB = TH2D("Hist_Dark_ElevNSB","",20,0,10,90,0,90);
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
    TString regularization_name;
    regularization_name  = "../Regularization_svd.root";
    if (EigenDecomposition)
    {
        regularization_name  = "../Regularization_eigen.root";
    }
    TFile RegularizationFile(regularization_name);
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
                    NormalizeDarkMatrix(&Hist_OneGroup_OffData_MSCLW.at(nth_sample).at(e), &Hist_OneGroup_OffDark_MSCLW.at(nth_sample).at(e));
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


    TH2D Hist_Regularization = TH2D("Hist_Regularization","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv);
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
    vector<TH2D> Hist_H_Vari;
    vector<TH2D> Hist_H_Vari_Bkgd;
    vector<TH2D> Hist_Coeff_Data;
    vector<TH2D> Hist_Coeff_Bkgd;
    vector<TH2D> Hist_Rank0_MSCLW_Data;
    vector<TH2D> Hist_Rank1_MSCLW_Data;
    vector<TH2D> Hist_Rank2_MSCLW_Data;
    vector<TH2D> Hist_Rank3_MSCLW_Data;
    vector<TH2D> Hist_Rank4_MSCLW_Data;
    vector<TH2D> Hist_Rank0_MSCLW_Dark;
    vector<TH2D> Hist_Rank1_MSCLW_Dark;
    vector<TH2D> Hist_Rank2_MSCLW_Dark;
    vector<TH2D> Hist_Rank3_MSCLW_Dark;
    vector<TH2D> Hist_Rank4_MSCLW_Dark;
    vector<TH1D> Hist_Bkgd_Optimization;
    vector<TH1D> Hist_Bkgd_Chi2;
    vector<TH1D> Hist_Bkgd_Converge_Blind;
    vector<TH1D> Hist_Bkgd_Converge_Unblind;
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
    vector<TH1D> Hist_Data_Eigenvalues;
    vector<TH1D> Hist_VVV_Eigenvalues;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));


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

        Hist_H_Vari.push_back(TH2D("Hist_H_Vari_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_H_Vari_Bkgd.push_back(TH2D("Hist_H_Vari_Bkgd_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Coeff_Data.push_back(TH2D("Hist_Coeff_Data_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Coeff_Bkgd.push_back(TH2D("Hist_Coeff_Bkgd_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv));

        Hist_Rank0_MSCLW_Data.push_back(TH2D("Hist_Rank0_MSCLW_Data_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank1_MSCLW_Data.push_back(TH2D("Hist_Rank1_MSCLW_Data_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank2_MSCLW_Data.push_back(TH2D("Hist_Rank2_MSCLW_Data_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank3_MSCLW_Data.push_back(TH2D("Hist_Rank3_MSCLW_Data_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank4_MSCLW_Data.push_back(TH2D("Hist_Rank4_MSCLW_Data_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank0_MSCLW_Dark.push_back(TH2D("Hist_Rank0_MSCLW_Dark_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank1_MSCLW_Dark.push_back(TH2D("Hist_Rank1_MSCLW_Dark_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank2_MSCLW_Dark.push_back(TH2D("Hist_Rank2_MSCLW_Dark_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank3_MSCLW_Dark.push_back(TH2D("Hist_Rank3_MSCLW_Dark_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank4_MSCLW_Dark.push_back(TH2D("Hist_Rank4_MSCLW_Dark_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));

        //double optimiz_lower = -6.;
        //double optimiz_upper = 6.;
        double optimiz_lower = -10.;
        double optimiz_upper = 0.;
        Hist_Bkgd_Optimization.push_back(TH1D("Hist_Bkgd_Optimization_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv*N_bins_for_deconv,optimiz_lower,optimiz_upper));
        Hist_Bkgd_Chi2.push_back(TH1D("Hist_Bkgd_Chi2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv*N_bins_for_deconv,optimiz_lower,optimiz_upper));
        Hist_Bkgd_Converge_Blind.push_back(TH1D("Hist_Bkgd_Converge_Blind_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",n_iterations,0,n_iterations));
        Hist_Bkgd_Converge_Unblind.push_back(TH1D("Hist_Bkgd_Converge_Unblind_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",n_iterations,0,n_iterations));

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
        Hist_Data_Eigenvalues.push_back(TH1D("Hist_Data_Eigenvalues_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_VVV_Eigenvalues.push_back(TH1D("Hist_VVV_Eigenvalues_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv*N_bins_for_deconv,0,N_bins_for_deconv*N_bins_for_deconv));
    }

    vector<TH1D> Hist_OnData_SR_Energy;
    vector<TH1D> Hist_OnData_SR_Energy_OneGroup;
    vector<TH1D> Hist_OnData_CR_Energy;
    vector<TH1D> Hist_OnData_CR_Energy_OneGroup;
    vector<TH1D> Hist_OnDark_CR_Energy;
    vector<TH1D> Hist_OnDark_CR_Energy_OneGroup;
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


        MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
        MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
        N_bins_for_deconv = N_bins_for_deconv_func_E[e];
        ResetMatrixDimension();

        Hist_OnData_SR_Energy.push_back(TH1D("Hist_OnData_SR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_SR_Energy_OneGroup.push_back(TH1D("Hist_OnData_SR_Energy_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_CR_Energy.push_back(TH1D("Hist_OnData_CR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_CR_Energy_OneGroup.push_back(TH1D("Hist_OnData_CR_Energy_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnDark_CR_Energy.push_back(TH1D("Hist_OnDark_CR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnDark_CR_Energy_OneGroup.push_back(TH1D("Hist_OnDark_CR_Energy_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
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
                biny_blind_global = Hist_OneGroup_Data_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind)-1;
                //binx_blind_global = Hist_OneGroup_Data_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind)-5;
                //biny_blind_global = Hist_OneGroup_Data_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind)-5;
                int binx_upper = Hist_OneGroup_Data_MSCLW.at(e).GetXaxis()->FindBin(1.)-1;
                int biny_lower = Hist_OneGroup_Data_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_lower);
                int biny_upper = Hist_OneGroup_Data_MSCLW.at(e).GetYaxis()->FindBin(1.)-1;

                hist_name  = "Hist2D_Regularization_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_Regularization.Reset();
                Hist_Regularization.Add( (TH2D*)RegularizationFile.Get(hist_name) );

                hist_name  = "Hist_OnData_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OneGroup_Data_MSCLW.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                group_size.at(e) += ((TH2D*)InputDataFile.Get(hist_name))->Integral();

                hist_name  = "Hist_GammaData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                //hist_name  = "Hist_GammaMC_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
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
                Hist_OnDark_CR_Energy_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
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

                mtx_regularization = fillMatrix(&Hist_Regularization);

                for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
                {

                    NormalizeDarkMatrix(&Hist_OneGroup_Data_MSCLW.at(e), &Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e));
                    TH2D hist_dark_alter = TH2D("hist_dark_alter","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
                    AlterDarkMatrix(&Hist_OneGroup_Data_MSCLW.at(e), &Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e), &hist_dark_alter);

                    fill2DHistogram(&Hist_Temp_Gamma,mtx_gamma_raw);
                    Gamma_SR_Integral = Hist_Temp_Gamma.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
                    Hist_Temp_Gamma.Scale((0.0)/Gamma_SR_Integral);
                    mtx_gamma = fillMatrix(&Hist_Temp_Gamma);

                    mtx_dark = fillMatrix(&hist_dark_alter);
                    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);

                    SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_dark);

                    if (UsePerturbation)
                    {
                        mtx_data_bkgd = PerturbationMethod(mtx_data,mtx_dark,true);
                    }
                    else
                    {
                        //double chi2_best = 1e20;
                        //int cutoff_max = NumberOfEigenvectors_Stable;
                        //int cutoff_best = 0;
                        //for (int cutoff=1;cutoff<=cutoff_max;cutoff++)
                        //{
                        //    NumberOfEigenvectors_Stable = cutoff;
                        //    LeastSquareSolutionMethod(false, true);
                        //    double chi2 = GetChi2Function(mtx_data_bkgd,mtx_data,true);
                        //    if (chi2<chi2_best)
                        //    {
                        //        chi2_best = chi2;
                        //        cutoff_best = cutoff;
                        //    }
                        //}
                        //NumberOfEigenvectors_Stable = cutoff_best;
                        
                        NumberOfEigenvectors_Stable = RankTruncation[e];
                        std::cout << "NumberOfEigenvectors_Stable = " << NumberOfEigenvectors_Stable << std::endl;
                        double optimized_alpha = pow(10.,Log10_alpha[e]);
                        if (RegularizationType==3)
                        {
                            optimized_alpha = double(30);
                        }
                        else if (RegularizationType==5)
                        {
                            optimized_alpha = double(256-30);
                        }
                        else if (RegularizationType==4)
                        {
                            optimized_alpha = 1e-1;
                        }
                        LeastSquareSolutionMethod(false, true, &Hist_Bkgd_Converge_Blind.at(e), &Hist_Bkgd_Optimization.at(e), &Hist_Bkgd_Chi2.at(e), optimized_alpha);
                        //LeastSquareSolutionMethod(false, false, &Hist_Bkgd_Converge_Blind.at(e));
                        for (int entry=0;entry<vtr_eigenval_vvv.size();entry++)
                        {
                            //int hist_entry = vtr_eigenval_vvv.size()-entry;
                            int hist_entry = entry+1;
                            Hist_VVV_Eigenvalues.at(e).SetBinContent(hist_entry,abs(vtr_eigenval_vvv(entry)));
                        }
                    }
                    
                    TH2D Hist_Temp_Bkgd = TH2D("Hist_Temp_Bkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
                    fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);

                    //if (UseMinChi2)
                    //{
                    //    mtx_data_bkgd = MinChi2Method(&Hist_OneGroup_OnSyst_Chi2.at(nth_sample).at(e), &Hist_OneGroup_Data_MSCLW.at(e), &Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e));
                    //    fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
                    //}

                    double Bkgd_SR_Integral = Hist_Temp_Bkgd.Integral(binx_lower,binx_blind,biny_lower,biny_blind);

                    //LeastSquareSolutionMethod(false, false, &Hist_Bkgd_Converge_Unblind.at(e));
                    //MatrixXcd mtx_dark_truncate = mtx_dark;
                    //if (TruncateNoise)
                    //{
                    //    mtx_dark_truncate = GetTruncatedMatrix(mtx_dark, NumberOfEigenvectors_Stable);
                    //}
                    //MatrixXcd mtx_diff = mtx_data_bkgd - mtx_dark_truncate;
                    //TH2D hist_diff = TH2D("hist_diff","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
                    //fill2DHistogram(&hist_diff,mtx_diff);
                    //Hist_Trunc_MSCLW_Diff.at(e).Add(&hist_diff);

                    double dark_weight = 1./double(n_dark_samples);
                    int error_bin = Hist_OneGroup_OnSyst_Chi2.at(nth_sample).at(e).GetMinimumBin();
                    Hist_OneGroup_Bkgd_MSCLW.at(e).Add(&Hist_Temp_Bkgd,dark_weight);
                    Hist_OnSyst_MSCLW.at(nth_sample).at(e).Add(&Hist_Temp_Bkgd);
                    if (TruncateNoise)
                    {
                        mtx_dark = GetTruncatedMatrix(mtx_dark, NumberOfEigenvectors_Stable);
                    }
                    fill2DHistogram(&hist_dark_alter,mtx_dark);
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
                    Hist_OnDark_CR_Energy_OneGroup.at(e).Scale(scale_dark);
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
                    Hist_OnDark_CR_Energy.at(e).Add(&Hist_OnDark_CR_Energy_OneGroup.at(e));
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
                Hist_OnDark_CR_Energy_OneGroup.at(e).Reset();
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

                std::cout << "Total data counts: " << Hist_OnData_MSCLW.at(e).Integral() << std::endl;
                double lambda_data_1 =  eigensolver_data.eigenvalues()(mtx_data.cols()-1).real();
                double lambda_data_2 =  eigensolver_data.eigenvalues()(mtx_data.cols()-2).real();
                double lambda_data_3 =  eigensolver_data.eigenvalues()(mtx_data.cols()-3).real();
                double lambda_dark_1 =  eigensolver_dark.eigenvalues()(mtx_data.cols()-1).real();
                double lambda_dark_2 =  eigensolver_dark.eigenvalues()(mtx_data.cols()-2).real();
                double lambda_dark_3 =  eigensolver_dark.eigenvalues()(mtx_data.cols()-3).real();
                double lambda_bkgd_1 =  eigensolver_bkgd.eigenvalues()(mtx_data.cols()-1).real();
                double lambda_bkgd_2 =  eigensolver_bkgd.eigenvalues()(mtx_data.cols()-2).real();
                double lambda_bkgd_3 =  eigensolver_bkgd.eigenvalues()(mtx_data.cols()-3).real();
                std::cout << "=====================================================================" << std::endl;
                std::cout << "lambda_data_1:" << std::endl;
                std::cout << lambda_data_1 << std::endl;
                std::cout << "lambda_data_2:" << std::endl;
                std::cout << lambda_data_2 << std::endl;
                std::cout << "lambda_data_3:" << std::endl;
                std::cout << lambda_data_3 << std::endl;
                std::cout << "=====================================================================" << std::endl;
                std::cout << "lambda_dark_1:" << std::endl;
                std::cout << lambda_dark_1 << " (" << (lambda_dark_1-lambda_data_1)/lambda_data_1 << ")" << std::endl;
                std::cout << "lambda_dark_2:" << std::endl;
                std::cout << lambda_dark_2 << " (" << (lambda_dark_2-lambda_data_2)/lambda_data_2 << ")" << std::endl;
                std::cout << "lambda_dark_3:" << std::endl;
                std::cout << lambda_dark_3 << " (" << (lambda_dark_3-lambda_data_3)/lambda_data_3 << ")" << std::endl;
                std::cout << "=====================================================================" << std::endl;
                std::cout << "lambda_bkgd_1:" << std::endl;
                std::cout << lambda_bkgd_1 << " (" << (lambda_bkgd_1-lambda_data_1)/lambda_data_1 << ")" << std::endl;
                std::cout << "lambda_bkgd_2:" << std::endl;
                std::cout << lambda_bkgd_2 << " (" << (lambda_bkgd_2-lambda_data_2)/lambda_data_2 << ")" << std::endl;
                std::cout << "lambda_bkgd_3:" << std::endl;
                std::cout << lambda_bkgd_3 << " (" << (lambda_bkgd_3-lambda_data_3)/lambda_data_3 << ")" << std::endl;
                std::cout << "=====================================================================" << std::endl;

                MatrixXcd mtx_data_truncate = GetTruncatedMatrix(mtx_data, NumberOfEigenvectors_Stable);
                double dark_full_count = CountGammaRegion(mtx_dark);
                double data_full_count = CountGammaRegion(mtx_data);
                double data_redu_count = CountGammaRegion(mtx_data_truncate);
                double epsilon_redu  = data_redu_count/data_full_count -1.;
                std::cout << "N^{ON} = " << data_full_count << ", N^{ON}_{3} = " << data_redu_count << ", epsilon^{ON}_{3} = " << epsilon_redu << std::endl;
                MatrixXcd mtx_data_replace = GetCRReplacedMatrix(mtx_data, mtx_dark);
                MatrixXcd mtx_data_replace_truncate = GetTruncatedMatrix(mtx_data_replace, NumberOfEigenvectors_Stable);
                double data_replace_redu_count = CountGammaRegion(mtx_data_replace_truncate);
                double epsilon_replace_redu  = data_replace_redu_count/data_full_count -1.;
                std::cout << "N^{ON} = " << data_full_count << ", N^{Swap}_{3} = " << data_replace_redu_count << ", epsilon^{Swap}_{3} = " << epsilon_replace_redu << std::endl;
                MatrixXcd mtx_bkgd_truncate = GetTruncatedMatrix(mtx_data_bkgd, NumberOfEigenvectors_Stable);
                double bkgd_redu_count = CountGammaRegion(mtx_bkgd_truncate);
                double epsilon_bkgd_redu  = bkgd_redu_count/data_full_count -1.;
                std::cout << "N^{ON} = " << data_full_count << ", N^{Chi2}_{3} = " << bkgd_redu_count << ", epsilon^{Chi2}_{3} = " << epsilon_bkgd_redu << std::endl;
                data_gamma_count.push_back(data_full_count);
                dark_gamma_count.push_back(dark_full_count);
                bkgd_gamma_count.push_back(bkgd_redu_count);
                MatrixXcd mtx_vvv_bkgd = PerturbationMethod(mtx_data,mtx_dark,true);
                MatrixXcd mtx_vvv_truncate = GetTruncatedMatrix(mtx_vvv_bkgd, NumberOfEigenvectors_Stable);
                double vvv_redu_count = CountGammaRegion(mtx_vvv_truncate);
                double epsilon_vvv_redu  = vvv_redu_count/data_full_count -1.;
                std::cout << "N^{ON} = " << data_full_count << ", N^{Pert,swap}_{3} = " << vvv_redu_count << ", epsilon^{Pert,swap}_{3} = " << epsilon_vvv_redu << std::endl;
                MatrixXcd mtx_vvv_noswap_bkgd = PerturbationMethod(mtx_data,mtx_dark,false);
                MatrixXcd mtx_vvv_noswap_truncate = GetTruncatedMatrix(mtx_vvv_noswap_bkgd, NumberOfEigenvectors_Stable);
                double vvv_noswap_redu_count = CountGammaRegion(mtx_vvv_noswap_truncate);
                double epsilon_vvv_noswap_redu  = vvv_noswap_redu_count/data_full_count -1.;
                std::cout << "N^{ON} = " << data_full_count << ", N^{Pert}_{3} = " << vvv_noswap_redu_count << ", epsilon^{Pert}_{3} = " << epsilon_vvv_noswap_redu << std::endl;

                double data_cr_count = CountCosmicRayRegion(mtx_data);
                MatrixXcd mtx_data_rescale = mtx_data/data_cr_count;
                MatrixXcd mtx_dark_rescale = mtx_dark/data_cr_count;
                MatrixXcd mtx_bkgd_rescale = mtx_data_bkgd/data_cr_count;

                eigensolver_bkgd = ComplexEigenSolver<MatrixXcd>(mtx_bkgd_rescale);
                eigensolver_bkgd_transpose = ComplexEigenSolver<MatrixXcd>(mtx_bkgd_rescale.transpose());
                MatrixXcd mtx_r_bkgd = eigensolver_bkgd.eigenvectors();
                MatrixXcd mtx_l_bkgd = eigensolver_bkgd_transpose.eigenvectors();
                MatrixXcd mtx_H_bkgd = mtx_r_bkgd.transpose()*mtx_l_bkgd;
                eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark_rescale);
                eigensolver_dark_transpose = ComplexEigenSolver<MatrixXcd>(mtx_dark_rescale.transpose());
                MatrixXcd mtx_r_dark = eigensolver_dark.eigenvectors();
                MatrixXcd mtx_l_dark = eigensolver_dark_transpose.eigenvectors();
                mtx_l_dark = SortEigenvectors(eigensolver_dark.eigenvalues(), mtx_r_dark, mtx_l_dark);
                MatrixXcd mtx_H_dark = mtx_r_dark.transpose()*mtx_l_dark;
                eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data_rescale);
                eigensolver_data_transpose = ComplexEigenSolver<MatrixXcd>(mtx_data_rescale.transpose());
                MatrixXcd mtx_r_data = eigensolver_data.eigenvectors();
                MatrixXcd mtx_l_data = eigensolver_data_transpose.eigenvectors();
                mtx_l_data = SortEigenvectors(eigensolver_data.eigenvalues(), mtx_r_data, mtx_l_data);
                MatrixXcd mtx_H_data = mtx_r_data.transpose()*mtx_l_data;
                MatrixXcd mtx_H_vari_true = (mtx_r_data-mtx_r_dark).transpose()*mtx_l_dark+mtx_r_dark.transpose()*(mtx_l_data-mtx_l_dark);
                MatrixXcd mtx_H_vari_bkgd = (mtx_r_bkgd-mtx_r_dark).transpose()*mtx_l_dark+mtx_r_dark.transpose()*(mtx_l_bkgd-mtx_l_dark);
                fill2DHistogramAbs(&Hist_H_Vari.at(e),mtx_H_vari_true);
                fill2DHistogramAbs(&Hist_H_Vari_Bkgd.at(e),mtx_H_vari_bkgd);

                MatrixXcd mtx_lambda_dark(mtx_dark.rows(),mtx_dark.cols());
                for (int row=0;row<mtx_dark.rows();row++)
                {
                    for (int col=0;col<mtx_dark.cols();col++)
                    {
                        mtx_lambda_dark(row,col) = 0.;
                        if (row==col)
                        {
                            mtx_lambda_dark(row,col) = eigensolver_dark.eigenvalues()(col);
                        }
                    }
                }

                MatrixXcd mtx_coeff_data = MatrixXcd::Zero(mtx_dark.rows(),mtx_dark.cols());
                MatrixXcd mtx_proj_data = mtx_l_dark.transpose()*(mtx_data_rescale-mtx_dark_rescale)*mtx_r_dark;
                MatrixXcd mtx_coeff_bkgd = MatrixXcd::Zero(mtx_dark.rows(),mtx_dark.cols());
                MatrixXcd mtx_proj_bkgd = mtx_l_dark.transpose()*(mtx_bkgd_rescale-mtx_dark_rescale)*mtx_r_dark;
                for (int col=0;col<mtx_data.cols();col++)
                {
                    for (int row=0;row<mtx_data.rows();row++)
                    {
                        if (row==col)
                        {
                            mtx_coeff_data(row,col) = mtx_proj_data(row,col);
                            mtx_coeff_bkgd(row,col) = mtx_proj_bkgd(row,col);
                            //mtx_coeff_data(row,col) = mtx_proj_data(row,col)/mtx_H_dark(row,col);
                            //mtx_coeff_bkgd(row,col) = mtx_proj_bkgd(row,col)/mtx_H_dark(row,col);
                        }
                        else
                        {
                            mtx_coeff_data(row,col) = mtx_proj_data(row,col);
                            mtx_coeff_bkgd(row,col) = mtx_proj_bkgd(row,col);
                            //mtx_coeff_data(row,col) = mtx_proj_data(row,col)/(mtx_lambda_dark(row,row)-mtx_lambda_dark(col,col));
                            //mtx_coeff_bkgd(row,col) = mtx_proj_bkgd(row,col)/(mtx_lambda_dark(row,row)-mtx_lambda_dark(col,col));
                        }
                    }
                }

                if (!EigenDecomposition)
                {
                    mtx_coeff_data = NuclearNormMinimization(mtx_dark_rescale, mtx_data_rescale, mtx_dark_rescale, 1, mtx_dark.rows(), 1.0, false, 0, 0.);
                    mtx_coeff_bkgd = NuclearNormMinimization(mtx_dark_rescale, mtx_bkgd_rescale, mtx_dark_rescale, 1, mtx_dark.rows(), 1.0, false, 0, 0.);
                }

                fill2DHistogramAbs(&Hist_Coeff_Data.at(e),mtx_coeff_data);
                fill2DHistogramAbs(&Hist_Coeff_Bkgd.at(e),mtx_coeff_bkgd);
                //int print_size = mtx_data.rows();
                int print_size = 10;
                std::cout << "mtx_H_data:" << std::endl;
                std::cout << mtx_H_data.block(mtx_data.rows()-print_size,mtx_data.cols()-print_size,print_size,print_size).cwiseAbs() << std::endl;
                std::cout << "mtx_H_dark:" << std::endl;
                std::cout << mtx_H_dark.block(mtx_data.rows()-print_size,mtx_data.cols()-print_size,print_size,print_size).cwiseAbs() << std::endl;
                std::cout << "mtx_H_vari_true:" << std::endl;
                std::cout << mtx_H_vari_true.block(mtx_data.rows()-print_size,mtx_data.cols()-print_size,print_size,print_size).cwiseAbs() << std::endl;
                std::cout << "mtx_coeff_data:" << std::endl;
                std::cout << mtx_coeff_data.block(mtx_data.rows()-print_size,mtx_data.cols()-print_size,print_size,print_size).cwiseAbs() << std::endl;
                std::cout << "mtx_coeff_bkgd:" << std::endl;
                std::cout << mtx_coeff_bkgd.block(mtx_data.rows()-print_size,mtx_data.cols()-print_size,print_size,print_size).cwiseAbs() << std::endl;

                if (EigenDecomposition)
                {
                    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
                    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
                    mtx_eigenvector_dark = eigensolver_dark.eigenvectors();
                    mtx_eigenvector_inv_dark = eigensolver_dark.eigenvectors().inverse();
                    mtx_eigenvector_data = eigensolver_data.eigenvectors();
                    mtx_eigenvector_inv_data = eigensolver_data.eigenvectors().inverse();
                    for (int i=0;i<mtx_data.cols();i++)
                    {
                        for (int j=0;j<mtx_data.rows();j++)
                        {
                            if (i==j) 
                            {
                                mtx_eigenval_data_redu(i,j) = eigensolver_data.eigenvalues()(i);
                                mtx_eigenval_dark_redu(i,j) = eigensolver_dark.eigenvalues()(i);
                            }
                        }
                    }
                }
                else
                {
                    JacobiSVD<MatrixXd> svd_dark(mtx_dark.real(), ComputeFullU | ComputeFullV);
                    mtx_eigenvector_dark = svd_dark.matrixU();
                    mtx_eigenvector_inv_dark = svd_dark.matrixV().transpose();
                    JacobiSVD<MatrixXd> svd_data(mtx_data.real(), ComputeFullU | ComputeFullV);
                    mtx_eigenvector_data = svd_data.matrixU();
                    mtx_eigenvector_inv_data = svd_data.matrixV().transpose();
                    for (int i=0;i<mtx_data.cols();i++)
                    {
                        for (int j=0;j<mtx_data.rows();j++)
                        {
                            if (i==j) 
                            {
                                mtx_eigenval_data_redu(i,j) = svd_data.singularValues()(i);
                                mtx_eigenval_dark_redu(i,j) = svd_dark.singularValues()(i);
                            }
                        }
                    }
                }
                for (int entry=0;entry<mtx_data.cols();entry++)
                {
                    int hist_entry = entry+1;
                    Hist_Data_Eigenvalues.at(e).SetBinContent(hist_entry,abs(mtx_eigenval_data_redu(entry,entry)));
                }
                GetReducedEigenvalueMatrix(0);
                mtx_data_redu = mtx_eigenvector_data*mtx_eigenval_data_redu*mtx_eigenvector_inv_data;
                fill2DHistogram(&Hist_Rank0_MSCLW_Data.at(e),mtx_data_redu);
                mtx_dark_redu = mtx_eigenvector_dark*mtx_eigenval_dark_redu*mtx_eigenvector_inv_dark;
                fill2DHistogram(&Hist_Rank0_MSCLW_Dark.at(e),mtx_dark_redu);
                GetReducedEigenvalueMatrix(1);
                mtx_data_redu = mtx_eigenvector_data*mtx_eigenval_data_redu*mtx_eigenvector_inv_data;
                fill2DHistogram(&Hist_Rank1_MSCLW_Data.at(e),mtx_data_redu);
                mtx_dark_redu = mtx_eigenvector_dark*mtx_eigenval_dark_redu*mtx_eigenvector_inv_dark;
                fill2DHistogram(&Hist_Rank1_MSCLW_Dark.at(e),mtx_dark_redu);
                GetReducedEigenvalueMatrix(2);
                mtx_data_redu = mtx_eigenvector_data*mtx_eigenval_data_redu*mtx_eigenvector_inv_data;
                fill2DHistogram(&Hist_Rank2_MSCLW_Data.at(e),mtx_data_redu);
                mtx_dark_redu = mtx_eigenvector_dark*mtx_eigenval_dark_redu*mtx_eigenvector_inv_dark;
                fill2DHistogram(&Hist_Rank2_MSCLW_Dark.at(e),mtx_dark_redu);
                GetReducedEigenvalueMatrix(3);
                mtx_data_redu = mtx_eigenvector_data*mtx_eigenval_data_redu*mtx_eigenvector_inv_data;
                fill2DHistogram(&Hist_Rank3_MSCLW_Data.at(e),mtx_data_redu);
                mtx_dark_redu = mtx_eigenvector_dark*mtx_eigenval_dark_redu*mtx_eigenvector_inv_dark;
                fill2DHistogram(&Hist_Rank3_MSCLW_Dark.at(e),mtx_dark_redu);
                GetReducedEigenvalueMatrix(4);
                mtx_data_redu = mtx_eigenvector_data*mtx_eigenval_data_redu*mtx_eigenvector_inv_data;
                fill2DHistogram(&Hist_Rank4_MSCLW_Data.at(e),mtx_data_redu);
                mtx_dark_redu = mtx_eigenvector_dark*mtx_eigenval_dark_redu*mtx_eigenvector_inv_dark;
                fill2DHistogram(&Hist_Rank4_MSCLW_Dark.at(e),mtx_dark_redu);


                SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data);
                fill1DHistogram(&Hist_Data_Rank0_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),0,1.0);
                fill1DHistogram(&Hist_Data_Rank1_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),1,1.0);
                fill1DHistogram(&Hist_Data_Rank2_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),2,1.0);
                fill1DHistogram(&Hist_Data_Rank3_LeftVector.at(e),mtx_eigenvector_inv_init.transpose(),3,1.0);
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
    hist_name  = "Hist_Data_ElevNSB";
    Hist_Data_ElevNSB.Add( (TH2D*)InputFile.Get(hist_name) );
    hist_name  = "Hist_Dark_ElevNSB";
    Hist_Dark_ElevNSB.Add( (TH2D*)InputFile.Get(hist_name) );
    hist_name  = "Hist_EffArea";
    Hist_EffArea.Add( (TH1D*)InputFile.Get(hist_name) );
    double NSB_mean_data = Hist_Data_ElevNSB.GetMean(1);
    double NSB_RMS_data = Hist_Data_ElevNSB.GetRMS(1);
    double NSB_mean_dark = Hist_Dark_ElevNSB.GetMean(1);
    double NSB_RMS_dark = Hist_Dark_ElevNSB.GetRMS(1);
    double Zenith_mean_data = Hist_Data_ElevNSB.GetMean(2);
    double Zenith_RMS_data = Hist_Data_ElevNSB.GetRMS(2);
    double Zenith_mean_dark = Hist_Dark_ElevNSB.GetMean(2);
    double Zenith_RMS_dark = Hist_Dark_ElevNSB.GetRMS(2);

    char lowrank_tag[50] = "";
    if (EigenDecomposition)
    {
        sprintf(lowrank_tag, "eigen");
    }
    else
    {
        sprintf(lowrank_tag, "svd");
    }
    TFile OutputFile("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+"_"+TString(output_file2_tag)+"_"+TString(lowrank_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+TString(signal_tag)+TString(mjd_cut_tag)+"_"+ONOFF_tag+".root","recreate");

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
    NewInfoTree.Branch("NSB_mean_data",&NSB_mean_data,"NSB_mean_data/D");
    NewInfoTree.Branch("NSB_RMS_data",&NSB_RMS_data,"NSB_RMS_data/D");
    NewInfoTree.Branch("NSB_mean_dark",&NSB_mean_dark,"NSB_mean_dark/D");
    NewInfoTree.Branch("NSB_RMS_dark",&NSB_RMS_dark,"NSB_RMS_dark/D");
    NewInfoTree.Branch("Zenith_mean_data",&Zenith_mean_data,"Zenith_mean_data/D");
    NewInfoTree.Branch("Zenith_RMS_data",&Zenith_RMS_data,"Zenith_RMS_data/D");
    NewInfoTree.Branch("Zenith_mean_dark",&Zenith_mean_dark,"Zenith_mean_dark/D");
    NewInfoTree.Branch("Zenith_RMS_dark",&Zenith_RMS_dark,"Zenith_RMS_dark/D");
    NewInfoTree.Branch("data_gamma_count","std::vector<double>",&data_gamma_count);
    NewInfoTree.Branch("dark_gamma_count","std::vector<double>",&dark_gamma_count);
    NewInfoTree.Branch("bkgd_gamma_count","std::vector<double>",&bkgd_gamma_count);
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
        Hist_H_Vari.at(e).Write();
        Hist_H_Vari_Bkgd.at(e).Write();
        Hist_Coeff_Data.at(e).Write();
        Hist_Coeff_Bkgd.at(e).Write();
        Hist_Rank0_MSCLW_Data.at(e).Write();
        Hist_Rank1_MSCLW_Data.at(e).Write();
        Hist_Rank2_MSCLW_Data.at(e).Write();
        Hist_Rank3_MSCLW_Data.at(e).Write();
        Hist_Rank4_MSCLW_Data.at(e).Write();
        Hist_Rank0_MSCLW_Dark.at(e).Write();
        Hist_Rank1_MSCLW_Dark.at(e).Write();
        Hist_Rank2_MSCLW_Dark.at(e).Write();
        Hist_Rank3_MSCLW_Dark.at(e).Write();
        Hist_Rank4_MSCLW_Dark.at(e).Write();
        Hist_Bkgd_Optimization.at(e).Write();
        Hist_Bkgd_Chi2.at(e).Write();
        Hist_Bkgd_Converge_Blind.at(e).Write();
        Hist_Bkgd_Converge_Unblind.at(e).Write();
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
        Hist_Data_Eigenvalues.at(e).Write();
        Hist_VVV_Eigenvalues.at(e).Write();
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
        Hist_OnDark_CR_Energy.at(e).Write();
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
