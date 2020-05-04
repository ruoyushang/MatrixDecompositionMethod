
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

#include "NetflixMethodGetShowerImage.C"

MatrixXcd mtx_gamma(N_bins_for_deconv,N_bins_for_deconv);
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
int NumberOfEigenvectors = 0;
int VaryLeftOrRightVector;
int VaryNthVector;

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
void NetflixParametrizeEigenvectors(const double *par)
{

    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector(row,col) = 0.;
            mtx_eigenvector_inv(row,col) = 0.;
            mtx_eigenvalue(row,col) = 0.;
        }
    }

    const std::complex<double> If(0.0, 1.0);
    // build eigenvector matrix
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {
        col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-2)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector_vari(row,col_fix) = par[first_index+row];
        }

        row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-1)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            mtx_eigenvector_inv_vari(row_fix,col) = par[first_index+col];
        }

        // build eigenvalue matrix
        first_index = 1+(2*NthEigenvector-0)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
        {
            mtx_eigenvalue_vari(N_bins_for_deconv-NthEigenvector,N_bins_for_deconv-NthEigenvalue) = par[first_index+NthEigenvalue-1];
        }
    }
}
MatrixXcd BuildModelMatrix()
{

    // build model
    MatrixXcd mtx_model(N_bins_for_deconv,N_bins_for_deconv);
    mtx_eigenvalue = mtx_eigenvalue_init+mtx_eigenvalue_vari;
    mtx_eigenvector = mtx_eigenvector_init+mtx_eigenvector_vari;
    mtx_eigenvector_inv = mtx_eigenvector_inv_init+mtx_eigenvector_inv_vari;

    mtx_model = mtx_eigenvector*mtx_eigenvalue*mtx_eigenvector_inv;

    return mtx_model;
}
MatrixXcd GetErrorMap()
{
    MatrixXcd mtx_error(N_bins_for_deconv,N_bins_for_deconv);
    mtx_error = mtx_dark-mtx_eigenvector_dark*mtx_eigenvalue_dark*mtx_eigenvector_inv_dark;
    return mtx_error;
}
double BlindedChi2(TH2D* hist_data, TH2D* hist_dark, TH2D* hist_model, TH2D* hist_error)
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
                continue;
            }
            if (bx<binx_blind_upper && by<biny_blind_upper)
            //if (bx<binx_blind_upper && by<biny_blind_upper && bx>=binx_blind_lower && by>=biny_blind_lower)
            {
                if ((data-model)<0.) chi2 += chi2_this;
                //continue;
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
double NetflixChi2Function(const double *par)
{

    NetflixParametrizeEigenvectors(par);
    
    // build model
    MatrixXcd mtx_model(N_bins_for_deconv,N_bins_for_deconv);
    mtx_model = BuildModelMatrix();

    MatrixXcd mtx_error = GetErrorMap();

    TH2D hist_gamma = TH2D("hist_gamma","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_diff = TH2D("hist_diff","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_data = TH2D("hist_data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_dark = TH2D("hist_dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_model = TH2D("hist_model","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_error = TH2D("hist_error","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);

    fill2DHistogramAbs(&hist_gamma,mtx_gamma);
    fill2DHistogram(&hist_data,mtx_data);
    fill2DHistogram(&hist_dark,mtx_dark);
    fill2DHistogram(&hist_model,mtx_model);
    fill2DHistogram(&hist_error,mtx_error);

    int binx_lower = hist_diff.GetXaxis()->FindBin(MSCL_cut_lower);
    int binx_blind = hist_diff.GetXaxis()->FindBin(MSCL_cut_blind)-1;
    int binx_upper = hist_diff.GetXaxis()->FindBin(1.)-1;
    int biny_lower = hist_diff.GetYaxis()->FindBin(MSCW_cut_lower);
    int biny_blind = hist_diff.GetYaxis()->FindBin(MSCW_cut_blind)-1;
    int biny_upper = hist_diff.GetYaxis()->FindBin(1.)-1;
    hist_diff.Add(&hist_data);
    hist_diff.Add(&hist_model,-1.);
    double gamma_total = hist_diff.Integral(binx_lower,binx_upper,biny_lower,biny_upper);
    gamma_total = max(0.,gamma_total);
    double scale = gamma_total/double(hist_gamma.Integral(binx_lower,binx_upper,biny_lower,biny_upper));
    hist_gamma.Scale(scale);

    double chi2 = 0.;
    chi2 += BlindedChi2(&hist_data,&hist_dark,&hist_model,&hist_error);
    //if (signal_model) chi2 += SignalChi2(&hist_data,&hist_gamma,&hist_model);

    return chi2;

}
void NetflixSetInitialVariables(ROOT::Math::GSLMinimizer* Chi2Minimizer, int binx_blind, int biny_blind, int fix_which, int which_to_fit)
{

    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);

    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    double scale = 0.;
    double limit = 0.; 

    ROOT::Math::GSLMinimizer Chi2Minimizer1D( ROOT::Math::kVectorBFGS );
    Chi2Minimizer1D.SetMaxFunctionCalls(200); // for Minuit/Minuit2
    Chi2Minimizer1D.SetMaxIterations(200); // for GSL
    Chi2Minimizer1D.SetTolerance(0.001);
    const double *par_1D;
    

    Chi2Minimizer->SetVariable(0, "par["+std::to_string(int(0))+"]", 1., 0.001);
    Chi2Minimizer->SetVariableLimits(0,1.,1.);
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        limit = 0.0;
        int col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-2)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            Chi2Minimizer->SetVariable(first_index+row,"par["+std::to_string(int(first_index+row))+"]",0.,0.01);
            if (fix_which==0) Chi2Minimizer->SetVariableLimits(first_index+row,0.,0.);
            if (NthEigenvector>which_to_fit) Chi2Minimizer->SetVariableLimits(first_index+row,0.,0.);
            if (which_to_fit!=NthEigenvector && which_to_fit!=-1) Chi2Minimizer->SetVariableLimits(first_index+row,0.,0.);
        }

        limit = 0.0;
        int row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-1)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            Chi2Minimizer->SetVariable(first_index+col,"par["+std::to_string(int(first_index+col))+"]",0.,0.01);
            if (fix_which==1) Chi2Minimizer->SetVariableLimits(first_index+col,0.,0.);
            if (NthEigenvector>which_to_fit) Chi2Minimizer->SetVariableLimits(first_index+col,0.,0.);
            if (which_to_fit!=NthEigenvector && which_to_fit!=-1) Chi2Minimizer->SetVariableLimits(first_index+col,0.,0.);
        }

        // eigenvalues
        double input_value = 0.;
        first_index = 1+(2*NthEigenvector-0)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
        {
            input_value = 0.;
            limit = 0.;
            Chi2Minimizer->SetVariable(first_index+NthEigenvalue-1, "par["+std::to_string(int(first_index+NthEigenvalue-1))+"]", input_value, 0.01*limit);
            Chi2Minimizer->SetVariableLimits(first_index+NthEigenvalue-1,input_value-limit,input_value+limit);
        }
    }

}
void SingleTimeMinimization(int fix_which, int which_to_fit)
{

    VaryLeftOrRightVector = fix_which;
    VaryNthVector = which_to_fit;

    ROOT::Math::Functor Chi2Func;
    //std::cout << "total n paramters = " << 1+2*NumberOfEigenvectors*(N_bins_for_deconv)+NumberOfEigenvectors << std::endl;
    Chi2Func = ROOT::Math::Functor(&NetflixChi2Function,1+2*NumberOfEigenvectors*(N_bins_for_deconv)+NumberOfEigenvectors*NumberOfEigenvectors); 

    // Choose method upon creation between:
    // kConjugateFR, kConjugatePR, kVectorBFGS,
    // kVectorBFGS2, kSteepestDescent

    //ROOT::Math::GSLMinimizer Chi2Minimizer_1st( ROOT::Math::kSteepestDescent );
    ROOT::Math::GSLMinimizer Chi2Minimizer_1st( ROOT::Math::kVectorBFGS2 );
    Chi2Minimizer_1st.SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    Chi2Minimizer_1st.SetMaxIterations(100); // for GSL
    Chi2Minimizer_1st.SetTolerance(0.1);
    Chi2Minimizer_1st.SetPrintLevel(0);
    Chi2Minimizer_1st.SetFunction(Chi2Func);

    NetflixSetInitialVariables(&Chi2Minimizer_1st,binx_blind_global,biny_blind_global,fix_which,which_to_fit);

    const double *par_1st = Chi2Minimizer_1st.X();
    //std::cout << "initial chi2 = " << NetflixChi2Function(par_1st) << std::endl;
    Chi2Minimizer_1st.SetTolerance(0.01*double(N_bins_for_deconv*N_bins_for_deconv));
    Chi2Minimizer_1st.Minimize();
    par_1st = Chi2Minimizer_1st.X();
    //std::cout << "final chi2 = " << NetflixChi2Function(par_1st) << std::endl;

    NetflixParametrizeEigenvectors(par_1st);

    mtx_data_bkgd = BuildModelMatrix();

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
            mtx_eigenvalue_vari(row,col) = 0.;
            mtx_eigenvector_dark(row,col) = mtx_U_r_dark(row,col);
            mtx_eigenvector_inv_dark(row,col) = mtx_U_l_dark.transpose()(row,col);
            mtx_eigenvalue_dark(row,col) = mtx_lambdanu_dark(row,col);
        }
    }

}
void MatrixFactorizationMethod()
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
    mtx_data_bkgd = mtx_eigenvector_init*mtx_eigenvalue_init*mtx_eigenvector_inv_init;

    for (int iteration=0;iteration<20;iteration++)
    {
        //std::cout << "iteration = " << iteration << std::endl;
        SingleTimeMinimization(0,1);
        SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
        SingleTimeMinimization(1,1);
        SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
    }
    for (int iteration=0;iteration<20;iteration++)
    {
        //std::cout << "iteration = " << iteration << std::endl;
        SingleTimeMinimization(0,2);
        SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
        SingleTimeMinimization(1,2);
        SetInitialSpectralvectors(binx_blind_global,biny_blind_global,mtx_data_bkgd);
    }

}
void NetflixMethodPrediction(string target_data, double PercentCrab, double tel_elev_lower_input, double tel_elev_upper_input, double MSCW_cut_input, double MSCL_cut_input, int rank)
{

    TH1::SetDefaultSumw2();
    sprintf(target, "%s", target_data.c_str());
    MSCW_cut_blind = MSCW_cut_input;
    MSCL_cut_blind = MSCL_cut_input;
    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;
    MSCW_plot_upper = gamma_hadron_dim_ratio*(MSCW_cut_input-MSCW_plot_lower)+MSCW_cut_input;
    MSCL_plot_upper = gamma_hadron_dim_ratio*(MSCL_cut_input-MSCL_plot_lower)+MSCL_cut_input;
    NumberOfEigenvectors = rank;

    TFile InputDataFile("../Netflix_"+TString(target_data)+"_Crab"+std::to_string(int(PercentCrab))+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+".root");
    TTree* InfoTree = nullptr;
    InfoTree = (TTree*) InputDataFile.Get("InfoTree");

    std::cout << "Working on OFF data..." << std::endl;
    vector<vector<TH2D>> Hist_OffBkgd_MSCLW;
    vector<vector<TH2D>> Hist_OffGamma_MSCLW;
    for (int nth_sample=0;nth_sample<n_control_samples-1;nth_sample++)
    {
        std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "sample = " << nth_sample << std::endl;
        char sample_tag[50];
        sprintf(sample_tag, "%i", nth_sample);
        vector<TH2D> Hist_OffBkgd_OneSample_MSCLW;
        vector<TH2D> Hist_OffGamma_OneSample_MSCLW;
        for (int e=0;e<N_energy_bins;e++) 
        {
            std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
            std::cout << "energy = " << energy_bins[e] << std::endl;
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            TString filename_gamma  = "Hist_GammaData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            if (PercentCrab>0.) filename_gamma  = "Hist_GammaMC_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            TString filename_data  = "Hist_OffData_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            TString filename_dark  = "Hist_OffDark_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
            TH2D* Hist_GammaMC = (TH2D*)InputDataFile.Get(filename_gamma);
            TH2D* Hist_Data = (TH2D*)InputDataFile.Get(filename_data);
            TH2D* Hist_Dark = (TH2D*)InputDataFile.Get(filename_dark);
            Hist_OffBkgd_OneSample_MSCLW.push_back(TH2D("Hist_OffBkgd_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_OffGamma_OneSample_MSCLW.push_back(TH2D("Hist_OffGamma_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));

            mtx_gamma = fillMatrix(Hist_GammaMC);
            mtx_data = fillMatrix(Hist_Data);
            mtx_dark = fillMatrix(Hist_Dark);
            eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
            eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);

            int binx_lower = Hist_Data->GetXaxis()->FindBin(MSCL_cut_lower);
            binx_blind_global = Hist_Data->GetXaxis()->FindBin(MSCL_cut_blind)-1;
            int binx_upper = Hist_Data->GetXaxis()->FindBin(1.)-1;
            int biny_lower = Hist_Data->GetYaxis()->FindBin(MSCW_cut_lower);
            biny_blind_global = Hist_Data->GetYaxis()->FindBin(MSCW_cut_blind)-1;
            int biny_upper = Hist_Data->GetYaxis()->FindBin(1.)-1;

            MatrixFactorizationMethod();

            fill2DHistogram(&Hist_OffBkgd_OneSample_MSCLW.at(e),mtx_data_bkgd);
            eigensolver_bkgd = ComplexEigenSolver<MatrixXcd>(mtx_data_bkgd);

            fill2DHistogramAbs(&Hist_OffGamma_OneSample_MSCLW.at(e),mtx_gamma);
            double gamma_total = Hist_Data->Integral(binx_lower,binx_upper,biny_lower,biny_upper)-Hist_OffBkgd_OneSample_MSCLW.at(e).Integral(binx_lower,binx_upper,biny_lower,biny_upper);
            gamma_total = max(0.,gamma_total);
            double scale_gamma = double(gamma_total)/double(Hist_OffGamma_OneSample_MSCLW.at(e).Integral(binx_lower,binx_upper,biny_lower,biny_upper));
            Hist_OffGamma_OneSample_MSCLW.at(e).Scale(scale_gamma);

        }
        Hist_OffBkgd_MSCLW.push_back(Hist_OffBkgd_OneSample_MSCLW);
        Hist_OffGamma_MSCLW.push_back(Hist_OffGamma_OneSample_MSCLW);
    }

    std::cout << "Working on ON data..." << std::endl;
    vector<TH2D> Hist_OnBkgd_MSCLW;
    vector<TH2D> Hist_OnGamma_MSCLW;
    for (int e=0;e<N_energy_bins;e++) 
    {
        std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "energy = " << energy_bins[e] << std::endl;
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        TString filename_gamma  = "Hist_GammaData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        if (PercentCrab>0.) filename_gamma  = "Hist_GammaMC_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TString filename_data  = "Hist_OnData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TString filename_dark  = "Hist_OnDark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TH2D* Hist_GammaMC = (TH2D*)InputDataFile.Get(filename_gamma);
        TH2D* Hist_Data = (TH2D*)InputDataFile.Get(filename_data);
        TH2D* Hist_Dark = (TH2D*)InputDataFile.Get(filename_dark);
        Hist_OnBkgd_MSCLW.push_back(TH2D("Hist_OnBkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnGamma_MSCLW.push_back(TH2D("Hist_OnGamma_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));

        mtx_gamma = fillMatrix(Hist_GammaMC);
        mtx_data = fillMatrix(Hist_Data);
        mtx_dark = fillMatrix(Hist_Dark);
        eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
        eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);

        int binx_lower = Hist_Data->GetXaxis()->FindBin(MSCL_cut_lower);
        binx_blind_global = Hist_Data->GetXaxis()->FindBin(MSCL_cut_blind)-1;
        int binx_upper = Hist_Data->GetXaxis()->FindBin(1.)-1;
        int biny_lower = Hist_Data->GetYaxis()->FindBin(MSCW_cut_lower);
        biny_blind_global = Hist_Data->GetYaxis()->FindBin(MSCW_cut_blind)-1;
        int biny_upper = Hist_Data->GetYaxis()->FindBin(1.)-1;

        MatrixFactorizationMethod();

        fill2DHistogram(&Hist_OnBkgd_MSCLW.at(e),mtx_data_bkgd);
        eigensolver_bkgd = ComplexEigenSolver<MatrixXcd>(mtx_data_bkgd);

        fill2DHistogramAbs(&Hist_OnGamma_MSCLW.at(e),mtx_gamma);
        double gamma_total = Hist_Data->Integral(binx_lower,binx_upper,biny_lower,biny_upper)-Hist_OnBkgd_MSCLW.at(e).Integral(binx_lower,binx_upper,biny_lower,biny_upper);
        gamma_total = max(0.,gamma_total);
        double scale_gamma = double(gamma_total)/double(Hist_OnGamma_MSCLW.at(e).Integral(binx_lower,binx_upper,biny_lower,biny_upper));
        Hist_OnGamma_MSCLW.at(e).Scale(scale_gamma);

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

    }

    InputDataFile.Close();

    TFile OutputFile("../Netflix_"+TString(target_data)+"_Crab"+std::to_string(int(PercentCrab))+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+".root","update");
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_OnBkgd_MSCLW.at(e).Write();
        Hist_OnGamma_MSCLW.at(e).Write();
    }
    for (int nth_sample=0;nth_sample<n_control_samples-1;nth_sample++)
    {
        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OffBkgd_MSCLW.at(nth_sample).at(e).Write();
            Hist_OffGamma_MSCLW.at(nth_sample).at(e).Write();
        }
    }

    std::cout << "Done." << std::endl;
}
