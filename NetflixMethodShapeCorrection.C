
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
#include <math.h> 

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

#include <complex>
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
using namespace Eigen;

#include "NetflixMethodGetShowerImage.C"

void NetflixMethodShapeCorrection(string target_data, double tel_elev_lower_input, double tel_elev_upper_input, int MJD_start_cut, int MJD_end_cut, bool isON)
{

    TH1::SetDefaultSumw2();

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
    sprintf(target, "%s", target_data.c_str());
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
    if (TString(target).Contains("SgrA"))
    {
        n_control_samples = 2;
    }
    MSCW_plot_upper = gamma_hadron_dim_ratio*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
    MSCL_plot_upper = gamma_hadron_dim_ratio*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;

    SizeSecondMax_Cut = 600.;
    if (TString(target).Contains("V4")) SizeSecondMax_Cut = 400.;
    if (TString(target).Contains("V5")) SizeSecondMax_Cut = 400.;


    mean_tele_point_ra = 0.;
    mean_tele_point_dec = 0.;
    pair<double,double> source_ra_dec = GetSourceRaDec(TString(target));
    mean_tele_point_ra = source_ra_dec.first;
    mean_tele_point_dec = source_ra_dec.second;

    GetBrightStars();
    GetGammaSources();


    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    TH1D Hist_ErecS_fine = TH1D("Hist_ErecS_fine","",N_energy_fine_bins,energy_fine_bins);

    vector<double>* roi_ra_ptr = new std::vector<double>(10);
    vector<double>* roi_dec_ptr = new std::vector<double>(10);
    vector<double>* roi_radius_ptr = new std::vector<double>(10);
    vector<string>* Dark_runlist_name_ptr = new std::vector<string>(10);
    vector<int>* Dark_runlist_number_ptr = new std::vector<int>(10);
    vector<string>* Data_runlist_name_ptr = new std::vector<string>(10);
    vector<int>* Data_runlist_number_ptr = new std::vector<int>(10);
    TFile InputDataFile("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+"_"+TString(output_file2_tag)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+TString(mjd_cut_tag)+"_"+ONOFF_tag+".root");
    TTree* InfoTree = nullptr;
    InfoTree = (TTree*) InputDataFile.Get("InfoTree");
    InfoTree->SetBranchAddress("roi_ra",&roi_ra_ptr);
    InfoTree->SetBranchAddress("roi_dec",&roi_dec_ptr);
    InfoTree->SetBranchAddress("roi_radius",&roi_radius_ptr);
    InfoTree->SetBranchAddress("Dark_runlist_number",&Dark_runlist_number_ptr);
    InfoTree->SetBranchAddress("Dark_runlist_name",&Dark_runlist_name_ptr);
    InfoTree->SetBranchAddress("Data_runlist_number",&Data_runlist_number_ptr);
    InfoTree->SetBranchAddress("Data_runlist_name",&Data_runlist_name_ptr);
    InfoTree->GetEntry(0);
    std::cout << "Dark_runlist_number_ptr->size() = " << Dark_runlist_number_ptr->size() << std::endl;
    std::cout << "Data_runlist_number_ptr->size() = " << Data_runlist_number_ptr->size() << std::endl;

    roi_ra.clear();
    roi_dec.clear();
    roi_radius.clear();
    for (int nth_roi=0;nth_roi<roi_ra_ptr->size();nth_roi++)
    {
        roi_ra.push_back(roi_ra_ptr->at(nth_roi));
        roi_dec.push_back(roi_dec_ptr->at(nth_roi));
        roi_radius.push_back(roi_radius_ptr->at(nth_roi));
    }


    vector<TH1D> Hist_OnDark_SR_Energy;
    vector<TH1D> Hist_OnDark_CR_Energy;
    vector<TH1D> Hist_OnData_SR_Energy;
    vector<TH1D> Hist_OnData_CR_Energy;

    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        Hist_OnDark_SR_Energy.push_back(TH1D("Hist_OnDark_SR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnDark_CR_Energy.push_back(TH1D("Hist_OnDark_CR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_SR_Energy.push_back(TH1D("Hist_OnData_SR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_CR_Energy.push_back(TH1D("Hist_OnData_CR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
    }

    vector<TH2D> Hist_OnDark_SR_CameraFoV;
    vector<TH2D> Hist_OnDark_CR_CameraFoV;
    vector<TH1D> Hist_OnDark_SR_Theta2;
    vector<TH1D> Hist_OnDark_CR_Theta2;
    vector<TH1D> Hist_OnData_SR_Skymap_Theta2;
    vector<TH1D> Hist_OnData_CR_Skymap_Theta2;
    vector<TH1D> Hist_OnData_CR_Skymap_Theta2_Raw;
    vector<TH1D> Hist_OnData_SR_CameraFoV_Theta2;
    vector<TH1D> Hist_OnData_CR_CameraFoV_Theta2;
    vector<TH1D> Hist_OnData_CR_CameraFoV_Theta2_Raw;
    vector<TH2D> Hist_OnData_SR_Skymap;
    vector<TH2D> Hist_OnData_CR_Skymap;
    vector<TH2D> Hist_OnData_CR_Skymap_Raw;
    vector<TH2D> Hist_OnData_SR_Skymap_Galactic;
    vector<TH2D> Hist_OnData_CR_Skymap_Galactic;
    vector<TH2D> Hist_OnData_SR_CameraFoV;
    vector<TH2D> Hist_OnData_CR_CameraFoV;
    vector<TH2D> Hist_OnData_CR_CameraFoV_Raw;
    for (int e=0;e<N_energy_fine_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_fine_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_fine_bins[e+1]));
        Hist_OnDark_SR_CameraFoV.push_back(TH2D("Hist_OnDark_SR_CameraFoV_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",40,0,10,6,0,2*M_PI));
        Hist_OnDark_CR_CameraFoV.push_back(TH2D("Hist_OnDark_CR_CameraFoV_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",40,0,10,6,0,2*M_PI));
        Hist_OnDark_SR_Theta2.push_back(TH1D("Hist_OnDark_SR_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,10));
        Hist_OnDark_CR_Theta2.push_back(TH1D("Hist_OnDark_CR_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,10));
        Hist_OnData_SR_Skymap_Theta2.push_back(TH1D("Hist_OnData_SR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_CR_Skymap_Theta2.push_back(TH1D("Hist_OnData_CR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_CR_Skymap_Theta2_Raw.push_back(TH1D("Hist_OnData_CR_Skymap_Theta2_Raw_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_SR_CameraFoV_Theta2.push_back(TH1D("Hist_OnData_SR_CameraFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_CR_CameraFoV_Theta2.push_back(TH1D("Hist_OnData_CR_CameraFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_CR_CameraFoV_Theta2_Raw.push_back(TH1D("Hist_OnData_CR_CameraFoV_Theta2_Raw_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_SR_Skymap.push_back(TH2D("Hist_OnData_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,mean_tele_point_ra-3,mean_tele_point_ra+3,150,mean_tele_point_dec-3,mean_tele_point_dec+3));
        Hist_OnData_CR_Skymap.push_back(TH2D("Hist_OnData_CR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,mean_tele_point_ra-3,mean_tele_point_ra+3,150,mean_tele_point_dec-3,mean_tele_point_dec+3));
        Hist_OnData_CR_Skymap_Raw.push_back(TH2D("Hist_OnData_CR_Skymap_Raw_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,mean_tele_point_ra-3,mean_tele_point_ra+3,150,mean_tele_point_dec-3,mean_tele_point_dec+3));

        pair<double,double> tele_point_l_b = ConvertRaDecToGalactic(mean_tele_point_ra, mean_tele_point_dec);
        mean_tele_point_l = tele_point_l_b.first;
        mean_tele_point_b = tele_point_l_b.second;
        Hist_OnData_SR_Skymap_Galactic.push_back(TH2D("Hist_OnData_SR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,tele_point_l_b.first-3,tele_point_l_b.first+3,150,tele_point_l_b.second-3,tele_point_l_b.second+3));
        Hist_OnData_CR_Skymap_Galactic.push_back(TH2D("Hist_OnData_CR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,tele_point_l_b.first-3,tele_point_l_b.first+3,150,tele_point_l_b.second-3,tele_point_l_b.second+3));

        Hist_OnData_SR_CameraFoV.push_back(TH2D("Hist_OnData_SR_CameraFoV_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,-3,3,150,-3,3));
        Hist_OnData_CR_CameraFoV.push_back(TH2D("Hist_OnData_CR_CameraFoV_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,-3,3,150,-3,3));
        Hist_OnData_CR_CameraFoV_Raw.push_back(TH2D("Hist_OnData_CR_CameraFoV_Raw_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,-3,3,150,-3,3));
    }

    vector<vector<TH1D>> Hist_OnData_SR_RoI_Energy;
    vector<vector<TH1D>> Hist_OnData_CR_RoI_Energy;
    vector<vector<TH1D>> Hist_OnData_SR_Skymap_RoI_Theta2;
    vector<vector<TH1D>> Hist_OnData_CR_Skymap_RoI_Theta2;
    vector<vector<TH1D>> Hist_OnData_SR_RoI_MJD;
    vector<vector<TH1D>> Hist_OnData_CR_RoI_MJD;
    for (int nth_roi=0;nth_roi<roi_ra.size();nth_roi++)
    {
        char roi_tag[50];
        sprintf(roi_tag, "%i", nth_roi);
        vector<TH1D> Hist_OnData_OneRoI_SR_RoI_Energy;
        vector<TH1D> Hist_OnData_OneRoI_CR_RoI_Energy;
        vector<TH1D> Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2;
        vector<TH1D> Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2;
        vector<TH1D> Hist_OnData_OneRoI_SR_RoI_MJD;
        vector<TH1D> Hist_OnData_OneRoI_CR_RoI_MJD;
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));
            Hist_OnData_OneRoI_SR_RoI_Energy.push_back(TH1D("Hist_OnData_SR_RoI_Energy_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
            Hist_OnData_OneRoI_CR_RoI_Energy.push_back(TH1D("Hist_OnData_CR_RoI_Energy_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        }
        Hist_OnData_SR_RoI_Energy.push_back(Hist_OnData_OneRoI_SR_RoI_Energy);
        Hist_OnData_CR_RoI_Energy.push_back(Hist_OnData_OneRoI_CR_RoI_Energy);
        for (int e=0;e<N_energy_fine_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_fine_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_fine_bins[e+1]));
            Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2.push_back(TH1D("Hist_OnData_SR_Skymap_RoI_Theta2_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,0.5));
            Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2.push_back(TH1D("Hist_OnData_CR_Skymap_RoI_Theta2_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,0.5));
            Hist_OnData_OneRoI_SR_RoI_MJD.push_back(TH1D("Hist_OnData_SR_RoI_MJD_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",800,56200-4000,56200+4000));
            Hist_OnData_OneRoI_CR_RoI_MJD.push_back(TH1D("Hist_OnData_CR_RoI_MJD_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",800,56200-4000,56200+4000));
        }
        Hist_OnData_SR_Skymap_RoI_Theta2.push_back(Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2);
        Hist_OnData_CR_Skymap_RoI_Theta2.push_back(Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2);
        Hist_OnData_SR_RoI_MJD.push_back(Hist_OnData_OneRoI_SR_RoI_MJD);
        Hist_OnData_CR_RoI_MJD.push_back(Hist_OnData_OneRoI_CR_RoI_MJD);
    }

    // Load MSCL/MSCW histograms
    vector<TH2D> Hist_OnBkgd_MSCLW;
    vector<TH2D> Hist_OnDark_MSCLW;
    for (int e=0;e<N_energy_bins;e++) 
    {
        std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "energy = " << energy_bins[e] << std::endl;
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));

        Hist_OnBkgd_MSCLW.push_back(TH2D("Hist_OnBkgd_MSCLW_new_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnDark_MSCLW.push_back(TH2D("Hist_OnDark_MSCLW_new_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));

        TString filename_bkgd  = "Hist_OnBkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TH2D* Hist_Bkgd = (TH2D*)InputDataFile.Get(filename_bkgd);
        for (int binx=0;binx<Hist_Bkgd->GetNbinsX();binx++)
        {
            for (int biny=0;biny<Hist_Bkgd->GetNbinsY();biny++)
            {
                //std::cout << "Hist_Bkgd->GetBinContent(binx+1,biny+1) = " << Hist_Bkgd->GetBinContent(binx+1,biny+1) << std::endl;
                Hist_OnBkgd_MSCLW.at(e).SetBinContent(binx+1,biny+1,Hist_Bkgd->GetBinContent(binx+1,biny+1));
            }
        }
    }
    InputDataFile.Cp("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+"_"+TString(output_file2_tag)+"_Modified_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+TString(mjd_cut_tag)+"_"+ONOFF_tag+".root");
    InputDataFile.Close();

    std::cout << "Prepare dark run samples..." << std::endl;
    for (int run=0;run<Dark_runlist_name_ptr->size();run++)
    {

        char run_number[50];
        char Dark_observation[50];
        sprintf(run_number, "%i", int(Dark_runlist_number_ptr->at(run)));
        sprintf(Dark_observation, "%s", Dark_runlist_name_ptr->at(run).c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Dark_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        //std::cout << "Get telescope pointing RA and Dec..." << std::endl;
        pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
        if (!TString(target).Contains("Proton")) tele_point_ra_dec = GetRunRaDec(filename,int(Dark_runlist_number_ptr->at(run)));
        run_tele_point_ra = tele_point_ra_dec.first;
        run_tele_point_dec = tele_point_ra_dec.second;

        TFile*  input_file = TFile::Open(filename.c_str());
        TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
        TTree* Dark_tree = (TTree*) input_file->Get(root_file);
        Dark_tree->SetBranchAddress("Xoff",&Xoff);
        Dark_tree->SetBranchAddress("Yoff",&Yoff);
        Dark_tree->SetBranchAddress("Xoff_derot",&Xoff_derot);
        Dark_tree->SetBranchAddress("Yoff_derot",&Yoff_derot);
        Dark_tree->SetBranchAddress("theta2",&theta2);
        Dark_tree->SetBranchAddress("ra",&ra_sky);
        Dark_tree->SetBranchAddress("dec",&dec_sky);
        Dark_tree->SetBranchAddress("ErecS",&ErecS);
        Dark_tree->SetBranchAddress("EChi2S",&EChi2S);
        Dark_tree->SetBranchAddress("MSCW",&MSCW);
        Dark_tree->SetBranchAddress("MSCL",&MSCL);
        Dark_tree->SetBranchAddress("NImages",&NImages);
        Dark_tree->SetBranchAddress("Xcore",&Xcore);
        Dark_tree->SetBranchAddress("Ycore",&Ycore);
        Dark_tree->SetBranchAddress("SizeSecondMax",&SizeSecondMax);
        Dark_tree->SetBranchAddress("EmissionHeight",&EmissionHeight);
        Dark_tree->SetBranchAddress("Time",&Time);
        Dark_tree->SetBranchAddress("Shower_Ze",&Shower_Ze);
        Dark_tree->SetBranchAddress("Shower_Az",&Shower_Az);

        for (int entry=0;entry<Dark_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            Xcore = 0;
            Ycore = 0;
            SizeSecondMax = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Dark_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            Phioff = atan2(Yoff,Xoff)+M_PI;
            ra_sky = tele_point_ra_dec.first+Xoff_derot;
            dec_sky = tele_point_ra_dec.second+Yoff_derot;
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            if (!SelectNImages(3,4)) continue;
            if (SizeSecondMax<SizeSecondMax_Cut) continue;
            if (EmissionHeight<6.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            //if (R2off>4.) continue;
            if (DarkFoV())
            {
                Hist_OnDark_MSCLW.at(energy).Fill(MSCL,MSCW);
            }
        }
        input_file->Close();
    }

    for (int run=0;run<Dark_runlist_number_ptr->size();run++)
    {

        char run_number[50];
        char Dark_observation[50];
        sprintf(run_number, "%i", int(Dark_runlist_number_ptr->at(run)));
        sprintf(Dark_observation, "%s", Dark_runlist_name_ptr->at(run).c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Dark_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        //std::cout << "Get telescope pointing RA and Dec..." << std::endl;
        pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
        if (!TString(target).Contains("Proton")) tele_point_ra_dec = GetRunRaDec(filename,int(Dark_runlist_number_ptr->at(run)));
        run_tele_point_ra = tele_point_ra_dec.first;
        run_tele_point_dec = tele_point_ra_dec.second;

        TFile*  input_file = TFile::Open(filename.c_str());
        TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
        TTree* Dark_tree = (TTree*) input_file->Get(root_file);
        Dark_tree->SetBranchAddress("Xoff",&Xoff);
        Dark_tree->SetBranchAddress("Yoff",&Yoff);
        Dark_tree->SetBranchAddress("Xoff_derot",&Xoff_derot);
        Dark_tree->SetBranchAddress("Yoff_derot",&Yoff_derot);
        Dark_tree->SetBranchAddress("theta2",&theta2);
        Dark_tree->SetBranchAddress("ra",&ra_sky);
        Dark_tree->SetBranchAddress("dec",&dec_sky);
        Dark_tree->SetBranchAddress("ErecS",&ErecS);
        Dark_tree->SetBranchAddress("EChi2S",&EChi2S);
        Dark_tree->SetBranchAddress("MSCW",&MSCW);
        Dark_tree->SetBranchAddress("MSCL",&MSCL);
        Dark_tree->SetBranchAddress("NImages",&NImages);
        Dark_tree->SetBranchAddress("Xcore",&Xcore);
        Dark_tree->SetBranchAddress("Ycore",&Ycore);
        Dark_tree->SetBranchAddress("SizeSecondMax",&SizeSecondMax);
        Dark_tree->SetBranchAddress("EmissionHeight",&EmissionHeight);
        Dark_tree->SetBranchAddress("Time",&Time);
        Dark_tree->SetBranchAddress("Shower_Ze",&Shower_Ze);
        Dark_tree->SetBranchAddress("Shower_Az",&Shower_Az);

        for (int entry=0;entry<Dark_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            Xcore = 0;
            Ycore = 0;
            SizeSecondMax = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Dark_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            Phioff = atan2(Yoff,Xoff)+M_PI;
            ra_sky = tele_point_ra_dec.first+Xoff_derot;
            dec_sky = tele_point_ra_dec.second+Yoff_derot;
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            if (!SelectNImages(3,4)) continue;
            if (SizeSecondMax<SizeSecondMax_Cut) continue;
            if (EmissionHeight<6.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            //if (R2off>4.) continue;
            if (DarkFoV())
            {
                int binx = Hist_OnDark_MSCLW.at(energy).GetXaxis()->FindBin(MSCL);
                int biny = Hist_OnDark_MSCLW.at(energy).GetYaxis()->FindBin(MSCW);
                double dark_content = Hist_OnDark_MSCLW.at(energy).GetBinContent(binx,biny);
                double bkgd_content = Hist_OnBkgd_MSCLW.at(energy).GetBinContent(binx,biny);
                double weight = 1.;
                if (dark_content>0.) weight = bkgd_content/dark_content;
                if (SignalSelectionTheta2())
                {
                    Hist_OnDark_SR_CameraFoV.at(energy_fine).Fill(R2off,Phioff,weight);
                    Hist_OnDark_SR_Theta2.at(energy_fine).Fill(Xoff*Xoff+Yoff*Yoff,weight);
                    Hist_OnDark_SR_Energy.at(energy).Fill(ErecS*1000.,weight);
                }
                if (ControlSelectionTheta2())
                {
                    Hist_OnDark_CR_CameraFoV.at(energy_fine).Fill(R2off,Phioff,weight);
                    Hist_OnDark_CR_Theta2.at(energy_fine).Fill(Xoff*Xoff+Yoff*Yoff,weight);
                    Hist_OnDark_CR_Energy.at(energy).Fill(ErecS*1000.,weight);
                }
            }
        }
        input_file->Close();
    }
    

    std::cout << "Prepare ON run samples..." << std::endl;
    for (int run=0;run<Data_runlist_name_ptr->size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist_number_ptr->at(run)));
        sprintf(Data_observation, "%s", Data_runlist_name_ptr->at(run).c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Data_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        //std::cout << "Get telescope pointing RA and Dec..." << std::endl;
        pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
        if (!TString(target).Contains("Proton")) tele_point_ra_dec = GetRunRaDec(filename,int(Data_runlist_number_ptr->at(run)));
        run_tele_point_ra = tele_point_ra_dec.first;
        run_tele_point_dec = tele_point_ra_dec.second;

        TFile*  input_file = TFile::Open(filename.c_str());
        TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
        TTree* Data_tree = (TTree*) input_file->Get(root_file);
        Data_tree->SetBranchAddress("Xoff",&Xoff);
        Data_tree->SetBranchAddress("Yoff",&Yoff);
        Data_tree->SetBranchAddress("theta2",&theta2);
        Data_tree->SetBranchAddress("Xoff_derot",&Xoff_derot);
        Data_tree->SetBranchAddress("Yoff_derot",&Yoff_derot);
        Data_tree->SetBranchAddress("ErecS",&ErecS);
        Data_tree->SetBranchAddress("EChi2S",&EChi2S);
        Data_tree->SetBranchAddress("MSCW",&MSCW);
        Data_tree->SetBranchAddress("MSCL",&MSCL);
        Data_tree->SetBranchAddress("NImages",&NImages);
        Data_tree->SetBranchAddress("Xcore",&Xcore);
        Data_tree->SetBranchAddress("Ycore",&Ycore);
        Data_tree->SetBranchAddress("SizeSecondMax",&SizeSecondMax);
        Data_tree->SetBranchAddress("EmissionHeight",&EmissionHeight);
        Data_tree->SetBranchAddress("Time",&Time);
        Data_tree->SetBranchAddress("Shower_Ze",&Shower_Ze);
        Data_tree->SetBranchAddress("Shower_Az",&Shower_Az);
        Data_tree->SetBranchAddress("MJD",&MJD);

        for (int entry=0;entry<Data_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            Xcore = 0;
            Ycore = 0;
            SizeSecondMax = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Data_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            Phioff = atan2(Yoff,Xoff)+M_PI;
            ra_sky = tele_point_ra_dec.first+Xoff_derot;
            dec_sky = tele_point_ra_dec.second+Yoff_derot;
            // redefine theta2
            theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
            pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            if (!SelectNImages(3,4)) continue;
            if (SizeSecondMax<SizeSecondMax_Cut) continue;
            if (EmissionHeight<6.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            //if (TString(target).Contains("Crab") && theta2<0.3) continue;
            //if (TString(target).Contains("Mrk421") && theta2<0.3) continue;
            //if (R2off>4.) continue;
            if (SignalSelectionTheta2())
            {
                if (FoV(true))
                {
                    Hist_OnData_SR_Energy.at(energy).Fill(ErecS*1000.);
                    for (int nth_roi=0;nth_roi<roi_ra.size();nth_roi++)
                    {
                        if (RoIFoV(nth_roi)) 
                        {
                            Hist_OnData_SR_RoI_Energy.at(nth_roi).at(energy).Fill(ErecS*1000.);
                            Hist_OnData_SR_RoI_MJD.at(nth_roi).at(energy_fine).Fill(MJD);
                        }
                        theta2_roi = pow(ra_sky-roi_ra.at(nth_roi),2)+pow(dec_sky-roi_dec.at(nth_roi),2);
                        Hist_OnData_SR_Skymap_RoI_Theta2.at(nth_roi).at(energy_fine).Fill(theta2_roi);
                    }
                    Hist_OnData_SR_Skymap_Theta2.at(energy_fine).Fill(theta2);
                    Hist_OnData_SR_Skymap.at(energy_fine).Fill(ra_sky,dec_sky);
                    Hist_OnData_SR_CameraFoV.at(energy_fine).Fill(Xoff,Yoff);
                    Hist_OnData_SR_Skymap_Galactic.at(energy_fine).Fill(evt_l_b.first,evt_l_b.second);
                    Hist_OnData_SR_CameraFoV_Theta2.at(energy_fine).Fill(R2off);
                }
            }
            if (ControlSelectionTheta2())
            {
                if (FoV(true))
                {
                    int binx = Hist_OnDark_SR_CameraFoV.at(energy_fine).GetXaxis()->FindBin(R2off);
                    int biny = Hist_OnDark_SR_CameraFoV.at(energy_fine).GetYaxis()->FindBin(Phioff);
                    double dark_cr_content = Hist_OnDark_CR_CameraFoV.at(energy_fine).GetBinContent(binx,biny);
                    double dark_sr_content = Hist_OnDark_SR_CameraFoV.at(energy_fine).GetBinContent(binx,biny);
                    double weight = 0.;
                    if (dark_cr_content>0.) weight = dark_sr_content/dark_cr_content;
                    int bin_e = Hist_OnDark_SR_Energy.at(energy).FindBin(ErecS*1000.);
                    double dark_cr_content_e = Hist_OnDark_CR_Energy.at(energy).GetBinContent(bin_e);
                    double dark_sr_content_e = Hist_OnDark_SR_Energy.at(energy).GetBinContent(bin_e);
                    double weight_e = 0.;
                    if (dark_cr_content_e>0.) weight_e = dark_sr_content_e/dark_cr_content_e;
                    Hist_OnData_CR_Energy.at(energy).Fill(ErecS*1000.,weight);
                    for (int nth_roi=0;nth_roi<roi_ra.size();nth_roi++)
                    {
                        if (RoIFoV(nth_roi)) 
                        {
                            Hist_OnData_CR_RoI_Energy.at(nth_roi).at(energy).Fill(ErecS*1000.,weight);
                            Hist_OnData_CR_RoI_MJD.at(nth_roi).at(energy_fine).Fill(MJD,weight);
                        }
                        theta2_roi = pow(ra_sky-roi_ra.at(nth_roi),2)+pow(dec_sky-roi_dec.at(nth_roi),2);
                        Hist_OnData_CR_Skymap_RoI_Theta2.at(nth_roi).at(energy_fine).Fill(theta2_roi,weight);
                    }
                    Hist_OnData_CR_Skymap_Theta2.at(energy_fine).Fill(theta2,weight);
                    Hist_OnData_CR_Skymap_Theta2_Raw.at(energy_fine).Fill(theta2,1.);
                    Hist_OnData_CR_CameraFoV_Theta2.at(energy_fine).Fill(R2off,weight);
                    Hist_OnData_CR_CameraFoV_Theta2_Raw.at(energy_fine).Fill(R2off,1.);
                    Hist_OnData_CR_Skymap.at(energy_fine).Fill(ra_sky,dec_sky,weight);
                    Hist_OnData_CR_Skymap_Raw.at(energy_fine).Fill(ra_sky,dec_sky,1.);
                    Hist_OnData_CR_Skymap_Galactic.at(energy_fine).Fill(evt_l_b.first,evt_l_b.second,weight);
                    Hist_OnData_CR_CameraFoV.at(energy_fine).Fill(Xoff,Yoff,weight);
                    Hist_OnData_CR_CameraFoV_Raw.at(energy_fine).Fill(Xoff,Yoff,1.);

                    //vector<vector<double>> new_locations = FillBrightStarHoles(ra_sky,dec_sky);
                    //for (int star=0;star<new_locations.size();star++)
                    //{
                    //    double new_ra_sky = new_locations.at(star).at(0);
                    //    double new_dec_sky = new_locations.at(star).at(1);
                    //    double new_weight_sky = new_locations.at(star).at(2);
                    //    if (new_weight_sky==0.) continue;
                    //    pair<double,double> new_evt_l_b = ConvertRaDecToGalactic(new_ra_sky,new_dec_sky);
                    //    Hist_OnData_CR_Energy.at(energy).Fill(ErecS*1000.,weight*new_weight_sky);
                    //    for (int nth_roi=0;nth_roi<roi_ra.size();nth_roi++)
                    //    {
                    //        theta2_roi = pow(new_ra_sky-roi_ra.at(nth_roi),2)+pow(new_dec_sky-roi_dec.at(nth_roi),2);
                    //        if (theta2_roi<bright_star_radius_cut*bright_star_radius_cut) 
                    //        {
                    //            Hist_OnData_CR_RoI_Energy.at(nth_roi).at(energy).Fill(ErecS*1000.,weight*new_weight_sky);
                    //            Hist_OnData_CR_RoI_MJD.at(nth_roi).at(energy_fine).Fill(MJD,weight*new_weight_sky);
                    //        }
                    //        Hist_OnData_CR_Skymap_RoI_Theta2.at(nth_roi).at(energy_fine).Fill(theta2_roi,weight*new_weight_sky);
                    //    }
                    //    double new_theta2 = pow(new_ra_sky-mean_tele_point_ra,2)+pow(new_dec_sky-mean_tele_point_dec,2);
                    //    Hist_OnData_CR_Skymap_Theta2.at(energy_fine).Fill(new_theta2,weight*new_weight_sky);
                    //    Hist_OnData_CR_Skymap.at(energy_fine).Fill(new_ra_sky,new_dec_sky,weight*new_weight_sky);
                    //    Hist_OnData_CR_Skymap_Galactic.at(energy_fine).Fill(new_evt_l_b.first,new_evt_l_b.second,weight*new_weight_sky);
                    //    Hist_OnData_CR_CameraFoV_Theta2.at(energy_fine).Fill(R2off,weight*new_weight_sky);
                    //    Hist_OnData_CR_CameraFoV.at(energy_fine).Fill(Xoff,Yoff,weight*new_weight_sky);
                    //}

                }
            }
        }
        input_file->Close();
    }

    TFile OutputFile("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+"_"+TString(output_file2_tag)+"_Modified_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+TString(mjd_cut_tag)+"_"+ONOFF_tag+".root","update");
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_OnDark_SR_Energy.at(e).Write();
        Hist_OnDark_CR_Energy.at(e).Write();
        Hist_OnData_SR_Energy.at(e).Write();
        Hist_OnData_CR_Energy.at(e).Write();
    }
    for (int nth_roi=0;nth_roi<roi_ra.size();nth_roi++)
    {
        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OnData_SR_RoI_Energy.at(nth_roi).at(e).Write();
            Hist_OnData_CR_RoI_Energy.at(nth_roi).at(e).Write();
        }
        for (int e=0;e<N_energy_fine_bins;e++) 
        {
            Hist_OnData_SR_RoI_MJD.at(nth_roi).at(e).Write();
            Hist_OnData_CR_RoI_MJD.at(nth_roi).at(e).Write();
            Hist_OnData_SR_Skymap_RoI_Theta2.at(nth_roi).at(e).Write();
            Hist_OnData_CR_Skymap_RoI_Theta2.at(nth_roi).at(e).Write();
        }
    }
    for (int e=0;e<N_energy_fine_bins;e++) 
    {
        Hist_OnDark_SR_CameraFoV.at(e).Write();
        Hist_OnDark_CR_CameraFoV.at(e).Write();
        Hist_OnDark_SR_Theta2.at(e).Write();
        Hist_OnDark_CR_Theta2.at(e).Write();
        Hist_OnData_SR_Skymap_Theta2.at(e).Write();
        Hist_OnData_CR_Skymap_Theta2.at(e).Write();
        Hist_OnData_CR_Skymap_Theta2_Raw.at(e).Write();
        Hist_OnData_SR_CameraFoV_Theta2.at(e).Write();
        Hist_OnData_CR_CameraFoV_Theta2.at(e).Write();
        Hist_OnData_CR_CameraFoV_Theta2_Raw.at(e).Write();
        Hist_OnData_SR_Skymap.at(e).Write();
        Hist_OnData_CR_Skymap.at(e).Write();
        Hist_OnData_CR_Skymap_Raw.at(e).Write();
        Hist_OnData_SR_Skymap_Galactic.at(e).Write();
        Hist_OnData_CR_Skymap_Galactic.at(e).Write();
        Hist_OnData_SR_CameraFoV.at(e).Write();
        Hist_OnData_CR_CameraFoV.at(e).Write();
        Hist_OnData_CR_CameraFoV_Raw.at(e).Write();
    }
    OutputFile.Close();

    std::cout << "Done." << std::endl;
}
