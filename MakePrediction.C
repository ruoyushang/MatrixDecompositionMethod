
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

MatrixXcd mtx_data(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_dark(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_data_bkgd(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvalue_dark(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvalue_init(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvalue_vari(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_dark(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_init(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_vari(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv_dark(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv_init(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv_vari(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvalue(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenval_data_redu(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_data_redu(N_bins_for_deconv,N_bins_for_deconv);
ComplexEigenSolver<MatrixXcd> eigensolver_bkgd;
ComplexEigenSolver<MatrixXcd> eigensolver_bkgd_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_data;
ComplexEigenSolver<MatrixXcd> eigensolver_data_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_dark;
ComplexEigenSolver<MatrixXcd> eigensolver_dark_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_init;
ComplexEigenSolver<MatrixXcd> eigensolver_init_transpose;
int binx_blind_global;
int biny_blind_global;
int VaryLeftOrRightVector;
int VaryNthVector;

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
void fill1DHistogram(TH1D* hist,MatrixXcd mtx, int rank)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        hist->SetBinContent(binx+1,mtx(binx,mtx.cols()-1-rank).real());
    }
    double integral = hist->Integral();
    hist->Scale(abs(integral)/integral);
}
MatrixXcd MakeRealEigenvectors(MatrixXcd mtx_input)
{
    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    for (int col=0;col<mtx_input.cols();col++)
    {
        if (col<mtx_input.cols()-NumberOfEigenvectors) continue;
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
                if (col>=mtx_input_lambda.cols()-NumberOfEigenvectors)
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
void SetInitialSpectralvectors(int binx_blind, int biny_blind, MatrixXcd mtx_input)
{

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

    MatrixXcd mtx_lambdanu = GetLambdaNuMatrix_v2(mtx_input,mtx_input);
    mtx_lambdanu = CutoffEigenvalueMatrix(mtx_lambdanu, NumberOfEigenvectors);

    MatrixXcd mtx_lambdanu_dark = GetLambdaNuMatrix_v2(mtx_dark,mtx_dark);
    mtx_lambdanu_dark = CutoffEigenvalueMatrix(mtx_lambdanu_dark, NumberOfEigenvectors);

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
            //if (solution_w_constraints && row<N_bins_for_deconv-1 && col<N_bins_for_deconv-1) 
            //{
            //    double scale = mtx_lambdanu(N_bins_for_deconv-1,N_bins_for_deconv-1).real()/mtx_lambdanu_dark(N_bins_for_deconv-1,N_bins_for_deconv-1).real();
            //    mtx_eigenvalue_init(row,col) = mtx_lambdanu_dark(row,col)*scale;
            //}
            mtx_eigenvalue_vari(row,col) = 0.;
            mtx_eigenvector_dark(row,col) = mtx_U_r_dark(row,col);
            mtx_eigenvector_inv_dark(row,col) = mtx_U_l_dark.transpose()(row,col);
            mtx_eigenvalue_dark(row,col) = mtx_lambdanu_dark(row,col);
        }
    }

}
MatrixXcd GetErrorMap()
{
    MatrixXcd mtx_error(N_bins_for_deconv,N_bins_for_deconv);
    mtx_error = mtx_dark-mtx_eigenvector_dark*mtx_eigenvalue_dark*mtx_eigenvector_inv_dark;
    return mtx_error;
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
double GetChi2Function(MatrixXcd mtx_model, int LeftOrRight)
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
    chi2 += BlindedChi2(&hist_data,&hist_dark,&hist_model,&hist_error,LeftOrRight);

    return chi2;

}
VectorXd SolutionWithConstraints(MatrixXd mtx_big, MatrixXd mtx_constraints, VectorXd vtr_delta)
{

    MatrixXd BTB = mtx_big.transpose()*mtx_big;
    VectorXd BTD = mtx_big.transpose()*vtr_delta;

    MatrixXd mtx_Bigger = MatrixXd::Zero(BTB.rows()+mtx_constraints.rows(),BTB.cols()+mtx_constraints.rows());
    mtx_Bigger.block(0,0,BTB.rows(),BTB.cols()) = BTB;
    mtx_Bigger.block(BTB.rows(),0,mtx_constraints.rows(),mtx_constraints.cols()) = mtx_constraints;
    mtx_Bigger.block(0,BTB.cols(),mtx_constraints.cols(),mtx_constraints.rows()) = mtx_constraints.transpose();

    VectorXd vtr_bigger_delta = VectorXd::Zero(BTB.rows()+mtx_constraints.rows());
    vtr_bigger_delta.segment(0,BTB.rows()) = BTD;

    VectorXd vtr_vari_bigger = VectorXd::Zero(BTB.cols()+mtx_constraints.rows());
    vtr_vari_bigger = mtx_Bigger.bdcSvd(ComputeThinU | ComputeThinV).solve(vtr_bigger_delta);

    return vtr_vari_bigger.segment(0,BTB.cols());

}
MatrixXcd SpectralDecompositionMethod_v3(int entry_start, int entry_size)
{

    double step_frac = 1.0;

    MatrixXcd mtx_q_init = mtx_eigenvector_init;
    MatrixXcd mtx_S = mtx_eigenvalue_init;
    MatrixXcd mtx_p_init = mtx_eigenvector_inv_init.transpose();
    MatrixXcd mtx_input = mtx_eigenvector_init*mtx_eigenvalue_init*mtx_eigenvector_inv_init;
    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());

    int row_size_big = mtx_input.rows()*mtx_input.cols();
    VectorXd vtr_Delta = VectorXd::Zero(row_size_big);
    MatrixXd mtx_Big = MatrixXd::Zero(row_size_big,2*entry_size*mtx_input.cols());
    MatrixXd mtx_Constraint = MatrixXd::Zero(NumberOfEigenvectors*entry_size,2*entry_size*mtx_input.cols());
    //MatrixXd mtx_Constraint = MatrixXd::Zero(1*entry_size,2*entry_size*mtx_input.cols());
    for (int idx_i=0; idx_i<mtx_input.rows(); idx_i++)
    {
        for (int idx_j=0; idx_j<mtx_input.cols(); idx_j++)
        {
            int idx_m = idx_j + idx_i*mtx_input.cols();
            double weight = 1.;
            //weight = 1./max(1.,pow(mtx_data(idx_i,idx_j).real(),0.5));
            if (idx_i<binx_blind_global && idx_j<biny_blind_global)
            {
                weight = 0.; // blind gamma-ray region
            }
            vtr_Delta(idx_m) = weight*(mtx_data-mtx_input)(idx_i,idx_j).real();
            for (int idx_k=0; idx_k<entry_size; idx_k++)
            {
                int idx_m = idx_j + mtx_input.cols()*idx_i;
                int idx_n = idx_j + mtx_input.cols()*idx_k;
                int idx_w = idx_i + mtx_input.cols()*idx_k + mtx_input.cols()*entry_size;
                int nth_entry = idx_k + entry_start;

                mtx_Big(idx_m,idx_n) = weight*mtx_S(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry).real()*mtx_q_init(idx_i,mtx_input.rows()-nth_entry).real();
                mtx_Big(idx_m,idx_w) = weight*mtx_S(mtx_input.rows()-nth_entry,mtx_input.rows()-nth_entry).real()*mtx_p_init(idx_j,mtx_input.rows()-nth_entry).real();
                for (int idx_l=0; idx_l<NumberOfEigenvectors; idx_l++)
                //for (int idx_l=0; idx_l<1; idx_l++)
                {
                    int nth_entry2 = idx_l + 1;
                    //int nth_entry2 = idx_k + entry_start;
                    int idx_u = idx_l + idx_k*NumberOfEigenvectors;
                    //int idx_u = idx_l + idx_k*1;
                    int idx_x = idx_l + idx_i*entry_size;
                    int idx_y = idx_k + idx_i*entry_size + mtx_input.cols()*entry_size;
                    if (nth_entry>=nth_entry2)
                    {
                        mtx_Constraint(idx_u,idx_x) = mtx_q_init(idx_i,mtx_input.rows()-nth_entry).real();
                        mtx_Constraint(idx_u,idx_y) = mtx_p_init(idx_i,mtx_input.rows()-nth_entry2).real();
                    }
                }
            }
        }
    }
    VectorXd vtr_vari_big = VectorXd::Zero(2*entry_size*mtx_input.cols());
    if (solution_w_constraints && entry_start>1) 
    {
        vtr_vari_big = SolutionWithConstraints(mtx_Big, mtx_Constraint, vtr_Delta);
        VectorXd vtr_should_be_zero = mtx_Constraint*vtr_vari_big;
        //std::cout << "entry_start = " << entry_start << std::endl;
        //std::cout << "mtx_Constraint = " << std::endl;
        //std::cout << mtx_Constraint << std::endl;
        //std::cout << "vtr_should_be_zero = " << std::endl;
        //std::cout << vtr_should_be_zero << std::endl;
        for (int i=0;i<vtr_should_be_zero.size();i++)
        {
            if (isnan(vtr_should_be_zero(i))) return mtx_input;
            if (abs(vtr_should_be_zero(i))>0.1) return mtx_input;
        }
    }
    else
    {
        vtr_vari_big = mtx_Big.bdcSvd(ComputeThinU | ComputeThinV).solve(vtr_Delta);
    }
    //VectorXd vtr_vari_big = VectorXd::Zero(2*entry_size*mtx_input.cols());
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
    mtx_p_final += mtx_p_vari*step_frac;
    mtx_q_final += mtx_q_vari*step_frac;

    mtx_output = mtx_eigenvector_init*mtx_eigenvalue_init*mtx_eigenvector_inv_init;
    mtx_output += mtx_q_vari*mtx_eigenvalue_init*mtx_eigenvector_inv_init;
    mtx_output += mtx_eigenvector_init*mtx_eigenvalue_init*mtx_p_vari.transpose();

    return mtx_output;

}
bool CheckIfEigenvalueMakeSense(MatrixXcd mtx_input, int rank)
{
    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
    double lambda_data = eigensolver_data.eigenvalues()(mtx_data.cols()-rank).real();
    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    double lambda_dark = eigensolver_dark.eigenvalues()(mtx_dark.cols()-rank).real();
    ComplexEigenSolver<MatrixXcd> eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    double lambda_input = eigensolver_input.eigenvalues()(mtx_input.cols()-rank).real();
    double diff_data = (lambda_data-lambda_dark)/lambda_dark;
    double diff_input = (lambda_input-lambda_dark)/lambda_dark;
    //std::cout << "rank = " << rank << std::endl;
    //std::cout << "lambda_data = " << lambda_data << std::endl;
    //std::cout << "lambda_dark = " << lambda_dark << std::endl;
    //std::cout << "lambda_input = " << lambda_input << std::endl;
    //if (rank<10)
    //{
    //    if (diff_data<0. && diff_input>0.) 
    //    {
    //        std::cout << "break at rank " << rank << " (unphysical result.)" << std::endl;
    //        return false;
    //    }
    //    //if (abs(diff_input)>0.2) 
    //    //{
    //    //    std::cout << "break at rank " << rank << " (deviate from initial by 20%.)" << std::endl;
    //    //    return false;
    //    //}
    //}

    double init_chi2 = GetChi2Function(mtx_data_bkgd,0);
    double current_chi2 = GetChi2Function(mtx_input,0);
    if (current_chi2>init_chi2) 
    {
        std::cout << "break at rank " << rank << " (chi2 increasing.)" << std::endl;
        return false;
    }

    return true;
}
void LeastSquareSolutionMethod(int rank_variation)
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

    std::cout << "initial chi2 = " << GetChi2Function(mtx_dark,0) << std::endl;

    MatrixXcd mtx_temp = mtx_dark;
    int n_iterations = 5;
    //int n_iterations = 2;
    for (int iteration=0;iteration<n_iterations;iteration++)
    {
        if (rank_variation==1)
        {
            mtx_temp = SpectralDecompositionMethod_v3(1, 1);
            if (!CheckIfEigenvalueMakeSense(mtx_temp,1)) break;
            mtx_data_bkgd = mtx_temp;
            SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
            //std::cout << "chi2 (rank 1) = " << GetChi2Function(mtx_data_bkgd,0) << std::endl;
        }
        else if (rank_variation==2)
        {
            mtx_temp = SpectralDecompositionMethod_v3(1, 1);
            if (!CheckIfEigenvalueMakeSense(mtx_temp,1)) break;
            mtx_data_bkgd = mtx_temp;
            SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
            //std::cout << "chi2 (rank 1) = " << GetChi2Function(mtx_data_bkgd,0) << std::endl;
            
            mtx_temp = SpectralDecompositionMethod_v3(2, 1);
            if (!CheckIfEigenvalueMakeSense(mtx_temp,2)) break;
            mtx_data_bkgd = mtx_temp;
            SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
            //std::cout << "chi2 (rank 2) = " << GetChi2Function(mtx_data_bkgd,0) << std::endl;
        }
        else
        {
            mtx_temp = SpectralDecompositionMethod_v3(1, 1);
            if (CheckIfEigenvalueMakeSense(mtx_temp,1))
            {
                mtx_data_bkgd = mtx_temp;
                SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
            }
            //
            mtx_temp = SpectralDecompositionMethod_v3(2, 1);
            if (CheckIfEigenvalueMakeSense(mtx_temp,2))
            {
                mtx_data_bkgd = mtx_temp;
                SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
            }
            //
            mtx_temp = SpectralDecompositionMethod_v3(3, 1);
            if (CheckIfEigenvalueMakeSense(mtx_temp,3))
            {
                mtx_data_bkgd = mtx_temp;
                SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
            }
        }
    }

    if (!CheckIfEigenvalueMakeSense(mtx_data_bkgd,1))
    {
        mtx_data_bkgd = mtx_dark;
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
void NormalizaDarkMatrix(TH2D* hist_data, TH2D* hist_dark)
{
    int binx_lower = hist_data->GetXaxis()->FindBin(MSCL_cut_lower);
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind)-1;
    int binx_upper = hist_data->GetXaxis()->FindBin(1.)-1;
    int biny_lower = hist_data->GetYaxis()->FindBin(MSCW_cut_lower);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind)-1;
    int biny_upper = hist_data->GetYaxis()->FindBin(1.)-1;
    double Data_SR_Integral = hist_data->Integral(binx_lower,binx_blind,biny_lower,biny_blind);
    double Data_Integral = hist_data->Integral();
    double Data_CR_Integral = Data_Integral-Data_SR_Integral;
    double Dark_SR_Integral = hist_dark->Integral(binx_lower,binx_blind,biny_lower,biny_blind);
    double Dark_Integral = hist_dark->Integral();
    double Dark_CR_Integral = Dark_Integral-Dark_SR_Integral;
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
void MakePrediction(string target_data, double tel_elev_lower_input, double tel_elev_upper_input, int MJD_start_cut, int MJD_end_cut, bool isON)
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

    //int rank_variation = 1;
    int rank_variation = NumberOfEigenvectors;

    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;
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
    MSCW_plot_upper = gamma_hadron_dim_ratio*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
    MSCL_plot_upper = gamma_hadron_dim_ratio*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;

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

    vector<int> group_size_limit;
    vector<int> group_size;
    for (int e=0;e<N_energy_bins;e++) 
    {
        group_size_limit.push_back(20);
        if (energy_bins[e]>=pow(10,2.6))
        {
            group_size_limit.at(e) = 20;
        }
        if (energy_bins[e]>=pow(10,3.0))
        {
            group_size_limit.at(e) = 100;
        }
        group_size.push_back(0);
    }

    vector<string>* Data_runlist_name_ptr = new std::vector<string>(10);
    vector<int>* Data_runlist_number_ptr = new std::vector<int>(10);
    vector<string>* roi_name_ptr = new std::vector<string>(10);
    TFile InputDataFile("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+TString(mjd_cut_tag)+"_"+ONOFF_tag+".root");
    TTree* InfoTree_ptr = nullptr;
    InfoTree_ptr = (TTree*) InputDataFile.Get("InfoTree");
    InfoTree_ptr->SetBranchAddress("Data_runlist_number",&Data_runlist_number_ptr);
    InfoTree_ptr->SetBranchAddress("Data_runlist_name",&Data_runlist_name_ptr);
    InfoTree_ptr->SetBranchAddress("roi_name",&roi_name_ptr);
    InfoTree_ptr->GetEntry(0);

    int FirstRun = 0;
    int LastRun = Data_runlist_name_ptr->size();
    //FirstRun = int(double(Data_runlist_name_ptr->size())*0.5);
    //LastRun = int(double(Data_runlist_name_ptr->size())*0.5);


    std::cout << "Working on OFF data..." << std::endl;
    vector<vector<TH2D>> Hist_OffData_MSCLW;
    vector<vector<TH2D>> Hist_OffBkgd_MSCLW;
    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
    {
        std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "sample = " << nth_sample << std::endl;
        char sample_tag[50];
        sprintf(sample_tag, "%i", nth_sample);
        vector<TH2D> Hist_OffData_OneSample_MSCLW;
        vector<TH2D> Hist_OffBkgd_OneSample_MSCLW;
        for (int e=0;e<N_energy_bins;e++) 
        {
            std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
            std::cout << "energy = " << energy_bins[e] << std::endl;
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));
            Hist_OffData_OneSample_MSCLW.push_back(TH2D("Hist_OffData2_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_OffBkgd_OneSample_MSCLW.push_back(TH2D("Hist_OffBkgd_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        }
        Hist_OffData_MSCLW.push_back(Hist_OffData_OneSample_MSCLW);
        Hist_OffBkgd_MSCLW.push_back(Hist_OffBkgd_OneSample_MSCLW);
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
            Hist_OneSample_OffData_MSCLW.push_back(TH2D("Hist_OneSample_OffData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_OneSample_OffDark_MSCLW.push_back(TH2D("Hist_OneSample_OffDark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
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

        for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
        {
            char sample2_tag[50];
            sprintf(sample2_tag, "%i", nth_sample);

            for (int on_run=FirstRun;on_run<LastRun;on_run++)
            {

                char sample_tag[50];
                sprintf(sample_tag, "%i", on_run);

                int binx_lower = Hist_OffData_MSCLW.at(0).at(0).GetXaxis()->FindBin(MSCL_cut_lower);
                binx_blind_global = Hist_OffData_MSCLW.at(0).at(0).GetXaxis()->FindBin(MSCL_cut_blind)-1;
                int binx_upper = Hist_OffData_MSCLW.at(0).at(0).GetXaxis()->FindBin(1.)-1;
                int biny_lower = Hist_OffData_MSCLW.at(0).at(0).GetYaxis()->FindBin(MSCW_cut_lower);
                biny_blind_global = Hist_OffData_MSCLW.at(0).at(0).GetYaxis()->FindBin(MSCW_cut_blind)-1;
                int biny_upper = Hist_OffData_MSCLW.at(0).at(0).GetYaxis()->FindBin(1.)-1;

                TString hist_name;
                hist_name  = "Hist_OnDark_MSCLW_R"+TString(sample_tag)+"_V"+TString(sample2_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                Hist_OneGroup_OffData_MSCLW.at(nth_sample).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                for (int nth_other_sample=0;nth_other_sample<n_dark_samples;nth_other_sample++)
                {
                    if (nth_sample==nth_other_sample) continue;
                    char sample2_tag[50];
                    sprintf(sample2_tag, "%i", nth_other_sample);
                    hist_name  = "Hist_OnDark_MSCLW_R"+TString(sample_tag)+"_V"+TString(sample2_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                    Hist_OneGroup_OffDark_MSCLW.at(nth_sample).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                }

                group_size.at(e) += 1;

                if (group_size.at(e)==group_size_limit.at(e) || on_run==Data_runlist_number_ptr->size()-1)
                {
                    mtx_data = fillMatrix(&Hist_OneGroup_OffData_MSCLW.at(nth_sample).at(e));
                    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
                    NormalizaDarkMatrix(&Hist_OneGroup_OffData_MSCLW.at(nth_sample).at(e), &Hist_OneGroup_OffDark_MSCLW.at(nth_sample).at(e));
                    mtx_dark = fillMatrix(&Hist_OneGroup_OffDark_MSCLW.at(nth_sample).at(e));
                    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
                    LeastSquareSolutionMethod(rank_variation);
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
            Hist_OneSample_Dark_MSCLW.push_back(TH2D("Hist_OneSample_Dark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        }
        Hist_OneGroup_Dark_MSCLW.push_back(Hist_OneSample_Dark_MSCLW);
    }
    vector<TH2D> Hist_OneGroup_Data_MSCLW;
    vector<TH2D> Hist_OneGroup_Bkgd_MSCLW;
    vector<TH2D> Hist_OnData_MSCLW;
    vector<TH2D> Hist_OnBkgd_MSCLW;
    vector<TH2D> Hist_OnDark_MSCLW;
    vector<TH2D> Hist_Rank0_MSCLW;
    vector<TH2D> Hist_Rank1_MSCLW;
    vector<TH2D> Hist_Rank2_MSCLW;
    vector<TH1D> Hist_Data_Rank0_LeftVector;
    vector<TH1D> Hist_Data_Rank1_LeftVector;
    vector<TH1D> Hist_Data_Rank2_LeftVector;
    vector<TH1D> Hist_Data_Rank0_RightVector;
    vector<TH1D> Hist_Data_Rank1_RightVector;
    vector<TH1D> Hist_Data_Rank2_RightVector;
    vector<TH1D> Hist_Bkgd_Rank0_LeftVector;
    vector<TH1D> Hist_Bkgd_Rank1_LeftVector;
    vector<TH1D> Hist_Bkgd_Rank2_LeftVector;
    vector<TH1D> Hist_Bkgd_Rank0_RightVector;
    vector<TH1D> Hist_Bkgd_Rank1_RightVector;
    vector<TH1D> Hist_Bkgd_Rank2_RightVector;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        Hist_OneGroup_Data_MSCLW.push_back(TH2D("Hist_OneGroup_Data_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OneGroup_Bkgd_MSCLW.push_back(TH2D("Hist_OneGroup_Bkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnData_MSCLW.push_back(TH2D("Hist_OnData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnBkgd_MSCLW.push_back(TH2D("Hist_OnBkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnDark_MSCLW.push_back(TH2D("Hist_OnDark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));

        Hist_Rank0_MSCLW.push_back(TH2D("Hist_Rank0_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank1_MSCLW.push_back(TH2D("Hist_Rank1_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank2_MSCLW.push_back(TH2D("Hist_Rank2_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));

        Hist_Data_Rank0_LeftVector.push_back(TH1D("Hist_Data_Rank0_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank1_LeftVector.push_back(TH1D("Hist_Data_Rank1_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank2_LeftVector.push_back(TH1D("Hist_Data_Rank2_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank0_RightVector.push_back(TH1D("Hist_Data_Rank0_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank1_RightVector.push_back(TH1D("Hist_Data_Rank1_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Rank2_RightVector.push_back(TH1D("Hist_Data_Rank2_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank0_LeftVector.push_back(TH1D("Hist_Bkgd_Rank0_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank1_LeftVector.push_back(TH1D("Hist_Bkgd_Rank1_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank2_LeftVector.push_back(TH1D("Hist_Bkgd_Rank2_LeftVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank0_RightVector.push_back(TH1D("Hist_Bkgd_Rank0_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank1_RightVector.push_back(TH1D("Hist_Bkgd_Rank1_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Bkgd_Rank2_RightVector.push_back(TH1D("Hist_Bkgd_Rank2_RightVector_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
    }

    vector<TH1D> Hist_OnData_SR_Energy;
    vector<TH1D> Hist_OnData_SR_Energy_OneGroup;
    vector<TH1D> Hist_OnData_CR_Energy;
    vector<TH1D> Hist_OnData_CR_Energy_OneGroup;
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
        Hist_OnData_SR_Energy.push_back(TH1D("Hist_OnData_SR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_SR_Energy_OneGroup.push_back(TH1D("Hist_OnData_SR_Energy_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_CR_Energy.push_back(TH1D("Hist_OnData_CR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_CR_Energy_OneGroup.push_back(TH1D("Hist_OnData_CR_Energy_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_SR_Zenith.push_back(TH1D("Hist_OnData_SR_Zenith_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
        Hist_OnData_SR_Zenith_OneGroup.push_back(TH1D("Hist_OnData_SR_Zenith_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
        Hist_OnData_CR_Zenith.push_back(TH1D("Hist_OnData_CR_Zenith_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
        Hist_OnData_CR_Zenith_OneGroup.push_back(TH1D("Hist_OnData_CR_Zenith_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
    }

    vector<TH1D> Hist_OnData_SR_Skymap_Theta2;
    vector<TH1D> Hist_OnData_SR_Skymap_Theta2_OneGroup;
    vector<TH1D> Hist_OnData_CR_Skymap_Theta2;
    vector<TH1D> Hist_OnData_CR_Skymap_Theta2_OneGroup;
    vector<TH2D> Hist_OnData_SR_Skymap;
    vector<TH2D> Hist_OnData_SR_Skymap_OneGroup;
    vector<TH2D> Hist_OnData_CR_Skymap;
    vector<TH2D> Hist_OnData_CR_Skymap_OneGroup;
    vector<TH2D> Hist_OnData_SR_Skymap_Galactic;
    vector<TH2D> Hist_OnData_SR_Skymap_Galactic_OneGroup;
    vector<TH2D> Hist_OnData_CR_Skymap_Galactic;
    vector<TH2D> Hist_OnData_CR_Skymap_Galactic_OneGroup;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        Hist_OnData_SR_Skymap_Theta2.push_back(TH1D("Hist_OnData_SR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_SR_Skymap_Theta2_OneGroup.push_back(TH1D("Hist_OnData_SR_Skymap_Theta2_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_CR_Skymap_Theta2.push_back(TH1D("Hist_OnData_CR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_CR_Skymap_Theta2_OneGroup.push_back(TH1D("Hist_OnData_CR_Skymap_Theta2_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_SR_Skymap.push_back(TH2D("Hist_OnData_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,150,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_SR_Skymap_OneGroup.push_back(TH2D("Hist_OnData_SR_Skymap_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,150,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_CR_Skymap.push_back(TH2D("Hist_OnData_CR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,150,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_CR_Skymap_OneGroup.push_back(TH2D("Hist_OnData_CR_Skymap_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,150,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));

        Hist_OnData_SR_Skymap_Galactic.push_back(TH2D("Hist_OnData_SR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,150,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_SR_Skymap_Galactic_OneGroup.push_back(TH2D("Hist_OnData_SR_Skymap_Galactic_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,150,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_CR_Skymap_Galactic.push_back(TH2D("Hist_OnData_CR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,150,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_CR_Skymap_Galactic_OneGroup.push_back(TH2D("Hist_OnData_CR_Skymap_Galactic_OneGroup_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,150,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));

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

        for (int on_run=FirstRun;on_run<LastRun;on_run++)
        {

            char sample_tag[50];
            sprintf(sample_tag, "%i", on_run);

            int binx_lower = Hist_OneGroup_Data_MSCLW.at(0).GetXaxis()->FindBin(MSCL_cut_lower);
            binx_blind_global = Hist_OneGroup_Data_MSCLW.at(0).GetXaxis()->FindBin(MSCL_cut_blind)-1;
            int binx_upper = Hist_OneGroup_Data_MSCLW.at(0).GetXaxis()->FindBin(1.)-1;
            int biny_lower = Hist_OneGroup_Data_MSCLW.at(0).GetYaxis()->FindBin(MSCW_cut_lower);
            biny_blind_global = Hist_OneGroup_Data_MSCLW.at(0).GetYaxis()->FindBin(MSCW_cut_blind)-1;
            int biny_upper = Hist_OneGroup_Data_MSCLW.at(0).GetYaxis()->FindBin(1.)-1;
            TString hist_name;
            hist_name  = "Hist_OnData_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            Hist_OneGroup_Data_MSCLW.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
            hist_name  = "Hist_OnData_SR_Skymap_Theta2_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            Hist_OnData_SR_Skymap_Theta2_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
            hist_name  = "Hist_OnData_CR_Skymap_Theta2_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            Hist_OnData_CR_Skymap_Theta2_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
            hist_name  = "Hist_OnData_SR_Skymap_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            Hist_OnData_SR_Skymap_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
            hist_name  = "Hist_OnData_CR_Skymap_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            Hist_OnData_CR_Skymap_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
            hist_name  = "Hist_OnData_SR_Skymap_Galactic_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            Hist_OnData_SR_Skymap_Galactic_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
            hist_name  = "Hist_OnData_CR_Skymap_Galactic_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            Hist_OnData_CR_Skymap_Galactic_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
            hist_name  = "Hist_OnData_SR_Energy_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            Hist_OnData_SR_Energy_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
            hist_name  = "Hist_OnData_CR_Energy_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            Hist_OnData_CR_Energy_OneGroup.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
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

            group_size.at(e) += 1;

            if (group_size.at(e)==group_size_limit.at(e) || on_run==Data_runlist_number_ptr->size()-1)
            {
                mtx_data = fillMatrix(&Hist_OneGroup_Data_MSCLW.at(e));
                eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
                for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
                {
                    NormalizaDarkMatrix(&Hist_OneGroup_Data_MSCLW.at(e), &Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e));
                    mtx_dark = fillMatrix(&Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e));
                    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
                    LeastSquareSolutionMethod(rank_variation);
                    TH2D Hist_Temp_Bkgd = TH2D("Hist_Temp_Bkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
                    fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
                    Hist_OneGroup_Bkgd_MSCLW.at(e).Add(&Hist_Temp_Bkgd,1./double(n_dark_samples));
                    Hist_OnDark_MSCLW.at(e).Add(&Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e),1./double(n_dark_samples));
                }
                if (!isnan(Hist_OneGroup_Bkgd_MSCLW.at(e).Integral()) && !isnan(Hist_OneGroup_Data_MSCLW.at(e).Integral()))
                {
                    Hist_OnData_MSCLW.at(e).Add(&Hist_OneGroup_Data_MSCLW.at(e));
                    Hist_OnBkgd_MSCLW.at(e).Add(&Hist_OneGroup_Bkgd_MSCLW.at(e));

                    int binx_lower = Hist_OnData_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_lower);
                    int binx_blind = Hist_OnData_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind)-1;
                    int binx_upper = Hist_OnData_MSCLW.at(e).GetXaxis()->FindBin(1.)-1;
                    int biny_lower = Hist_OnData_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_lower);
                    int biny_blind = Hist_OnData_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind)-1;
                    int biny_upper = Hist_OnData_MSCLW.at(e).GetYaxis()->FindBin(1.)-1;
                    double Bkgd_SR_Integral = Hist_OneGroup_Bkgd_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
                    double Old_Integral = Hist_OnData_CR_Energy_OneGroup.at(e).Integral();
                    double scale = 0.;
                    if (Old_Integral>0.)
                    {
                        scale = Bkgd_SR_Integral/Old_Integral;
                    }
                    Hist_OnData_CR_Energy_OneGroup.at(e).Scale(scale);
                    Hist_OnData_CR_Zenith_OneGroup.at(e).Scale(scale);
                    Hist_OnData_CR_Skymap_Theta2_OneGroup.at(e).Scale(scale);
                    Hist_OnData_CR_Skymap_OneGroup.at(e).Scale(scale);
                    Hist_OnData_CR_Skymap_Galactic_OneGroup.at(e).Scale(scale);
                    for (int nth_roi=0;nth_roi<roi_name_ptr->size();nth_roi++)
                    {
                        Hist_OnData_CR_RoI_Energy_OneGroup.at(nth_roi).at(e).Scale(scale);
                        Hist_OnData_CR_Skymap_RoI_Theta2_OneGroup.at(nth_roi).at(e).Scale(scale);
                        Hist_OnData_CR_RoI_MJD_OneGroup.at(nth_roi).at(e).Scale(scale);
                    }
                    Hist_OnData_SR_Energy.at(e).Add(&Hist_OnData_SR_Energy_OneGroup.at(e));
                    Hist_OnData_CR_Energy.at(e).Add(&Hist_OnData_CR_Energy_OneGroup.at(e));
                    Hist_OnData_SR_Zenith.at(e).Add(&Hist_OnData_SR_Zenith_OneGroup.at(e));
                    Hist_OnData_CR_Zenith.at(e).Add(&Hist_OnData_CR_Zenith_OneGroup.at(e));
                    Hist_OnData_SR_Skymap_Theta2.at(e).Add(&Hist_OnData_SR_Skymap_Theta2_OneGroup.at(e));
                    Hist_OnData_CR_Skymap_Theta2.at(e).Add(&Hist_OnData_CR_Skymap_Theta2_OneGroup.at(e));
                    Hist_OnData_SR_Skymap.at(e).Add(&Hist_OnData_SR_Skymap_OneGroup.at(e));
                    Hist_OnData_CR_Skymap.at(e).Add(&Hist_OnData_CR_Skymap_OneGroup.at(e));
                    Hist_OnData_SR_Skymap_Galactic.at(e).Add(&Hist_OnData_SR_Skymap_Galactic_OneGroup.at(e));
                    Hist_OnData_CR_Skymap_Galactic.at(e).Add(&Hist_OnData_CR_Skymap_Galactic_OneGroup.at(e));
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
                Hist_OnData_SR_Energy_OneGroup.at(e).Reset();
                Hist_OnData_CR_Energy_OneGroup.at(e).Reset();
                Hist_OnData_SR_Zenith_OneGroup.at(e).Reset();
                Hist_OnData_CR_Zenith_OneGroup.at(e).Reset();
                Hist_OnData_SR_Skymap_Theta2_OneGroup.at(e).Reset();
                Hist_OnData_CR_Skymap_Theta2_OneGroup.at(e).Reset();
                Hist_OnData_SR_Skymap_OneGroup.at(e).Reset();
                Hist_OnData_CR_Skymap_OneGroup.at(e).Reset();
                Hist_OnData_SR_Skymap_Galactic_OneGroup.at(e).Reset();
                Hist_OnData_CR_Skymap_Galactic_OneGroup.at(e).Reset();
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
                std::cout << "=====================================================================" << std::endl;
                std::cout << "eigensolver_dark.eigenvalues()(mtx_dark.cols()-1):" << std::endl;
                std::cout << eigensolver_dark.eigenvalues()(mtx_dark.cols()-1) << std::endl;
                std::cout << "eigensolver_dark.eigenvalues()(mtx_dark.cols()-2):" << std::endl;
                std::cout << eigensolver_dark.eigenvalues()(mtx_dark.cols()-2) << std::endl;
                std::cout << "eigensolver_dark.eigenvalues()(mtx_dark.cols()-3):" << std::endl;
                std::cout << eigensolver_dark.eigenvalues()(mtx_dark.cols()-3) << std::endl;
                std::cout << "=====================================================================" << std::endl;
                std::cout << "eigensolver_bkgd.eigenvalues()(mtx_data.cols()-1):" << std::endl;
                std::cout << eigensolver_bkgd.eigenvalues()(mtx_data.cols()-1) << std::endl;
                std::cout << "eigensolver_bkgd.eigenvalues()(mtx_data.cols()-2):" << std::endl;
                std::cout << eigensolver_bkgd.eigenvalues()(mtx_data.cols()-2) << std::endl;
                std::cout << "eigensolver_bkgd.eigenvalues()(mtx_data.cols()-3):" << std::endl;
                std::cout << eigensolver_bkgd.eigenvalues()(mtx_data.cols()-3) << std::endl;
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

                fill1DHistogram(&Hist_Data_Rank0_LeftVector.at(e),eigensolver_data.eigenvectors().inverse().transpose(),0);
                fill1DHistogram(&Hist_Data_Rank1_LeftVector.at(e),eigensolver_data.eigenvectors().inverse().transpose(),1);
                fill1DHistogram(&Hist_Data_Rank2_LeftVector.at(e),eigensolver_data.eigenvectors().inverse().transpose(),2);
                fill1DHistogram(&Hist_Data_Rank0_RightVector.at(e),eigensolver_data.eigenvectors(),0);
                fill1DHistogram(&Hist_Data_Rank1_RightVector.at(e),eigensolver_data.eigenvectors(),1);
                fill1DHistogram(&Hist_Data_Rank2_RightVector.at(e),eigensolver_data.eigenvectors(),2);
                fill1DHistogram(&Hist_Bkgd_Rank0_LeftVector.at(e),eigensolver_bkgd.eigenvectors().inverse().transpose(),0);
                fill1DHistogram(&Hist_Bkgd_Rank1_LeftVector.at(e),eigensolver_bkgd.eigenvectors().inverse().transpose(),1);
                fill1DHistogram(&Hist_Bkgd_Rank2_LeftVector.at(e),eigensolver_bkgd.eigenvectors().inverse().transpose(),2);
                fill1DHistogram(&Hist_Bkgd_Rank0_RightVector.at(e),eigensolver_bkgd.eigenvectors(),0);
                fill1DHistogram(&Hist_Bkgd_Rank1_RightVector.at(e),eigensolver_bkgd.eigenvectors(),1);
                fill1DHistogram(&Hist_Bkgd_Rank2_RightVector.at(e),eigensolver_bkgd.eigenvectors(),2);
            
                group_size.at(e) = 0;

            }
        }
    }
    InputDataFile.Close();


    for (int e=0;e<N_energy_fine_bins;e++) 
    {
        Hist_EffArea.SetBinContent(e+1,Hist_EffArea.GetBinContent(e+1)/(3600.*exposure_hours));
    }

    TFile InputFile("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+TString(mjd_cut_tag)+"_"+ONOFF_tag+".root");
    TTree* InfoTree = nullptr;
    InfoTree = (TTree*) InputFile.Get("InfoTree");
    TTree* StarTree = nullptr;
    StarTree = (TTree*) InputFile.Get("StarTree");
    TTree* FaintStarTree = nullptr;
    FaintStarTree = (TTree*) InputFile.Get("FaintStarTree");

    TFile OutputFile("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+"_"+TString(output_file2_tag)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+TString(mjd_cut_tag)+"_"+ONOFF_tag+".root","recreate");

    TTree *newtree = InfoTree->CloneTree();
    newtree->Write();
    TTree *newstartree = StarTree->CloneTree();
    newstartree->Write();
    TTree *newfaintstartree = FaintStarTree->CloneTree();
    newfaintstartree->Write();

    Hist_EffArea.Write();
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_OnData_MSCLW.at(e).Write();
        Hist_OnBkgd_MSCLW.at(e).Write();
        Hist_OnDark_MSCLW.at(e).Write();
        Hist_Rank0_MSCLW.at(e).Write();
        Hist_Rank1_MSCLW.at(e).Write();
        Hist_Rank2_MSCLW.at(e).Write();
        Hist_Data_Rank0_LeftVector.at(e).Write();
        Hist_Data_Rank1_LeftVector.at(e).Write();
        Hist_Data_Rank2_LeftVector.at(e).Write();
        Hist_Data_Rank0_RightVector.at(e).Write();
        Hist_Data_Rank1_RightVector.at(e).Write();
        Hist_Data_Rank2_RightVector.at(e).Write();
        Hist_Bkgd_Rank0_LeftVector.at(e).Write();
        Hist_Bkgd_Rank1_LeftVector.at(e).Write();
        Hist_Bkgd_Rank2_LeftVector.at(e).Write();
        Hist_Bkgd_Rank0_RightVector.at(e).Write();
        Hist_Bkgd_Rank1_RightVector.at(e).Write();
        Hist_Bkgd_Rank2_RightVector.at(e).Write();
    }
    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
    {
        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OffData_MSCLW.at(nth_sample).at(e).Write();
            Hist_OffBkgd_MSCLW.at(nth_sample).at(e).Write();
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
        Hist_OnData_SR_Skymap.at(e).Write();
        Hist_OnData_CR_Skymap.at(e).Write();
        Hist_OnData_SR_Skymap_Galactic.at(e).Write();
        Hist_OnData_CR_Skymap_Galactic.at(e).Write();
    }
    OutputFile.Close();
    InputFile.Close();

    std::cout << "Done." << std::endl;

}
