
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

#include "PrepareDarkData.C"

void ObservingEffect()
{

    SMI_INPUT = string(std::getenv("SMI_INPUT"));
    SMI_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    SMI_DIR = string(std::getenv("SMI_DIR"));
    SMI_AUX = string(std::getenv("SMI_AUX"));

    TH1::SetDefaultSumw2();

    camera_theta2_cut_lower = 0.;
    camera_theta2_cut_upper = 9.;
    TelElev_lower = 25.;
    TelElev_upper = 85.;
    MSCW_cut_blind = MSCW_cut_moderate;
    MSCL_cut_blind = MSCL_cut_moderate;
    source_theta2_cut = 0.03;

    SizeSecondMax_Cut = 600.;
    if (TString(target).Contains("V4")) SizeSecondMax_Cut = 400.;
    if (TString(target).Contains("V5")) SizeSecondMax_Cut = 400.;

    vector<string> source_list;
    source_list.push_back("OJ287V6");
    source_list.push_back("M82V6");
    source_list.push_back("M82V5");
    source_list.push_back("Segue1V6");
    source_list.push_back("Segue1V5");
    source_list.push_back("1ES0647V6");
    source_list.push_back("1ES1011V6");
    source_list.push_back("3C264V6");
    source_list.push_back("PKS1424V6");
    source_list.push_back("PKS1424V5");
    source_list.push_back("H1426V6");
    source_list.push_back("1ES0229V6");
    source_list.push_back("1ES0229V5");
    source_list.push_back("SNR_G150p3Plus04p5_V6");
    source_list.push_back("DracoV6");
    source_list.push_back("DracoV5");
    source_list.push_back("1ES0502V6");
    source_list.push_back("1ES0502V5");
    source_list.push_back("BLLacV6");
    source_list.push_back("BLLacV5");
    source_list.push_back("RGBJ0710V5");
    source_list.push_back("1ES0414V5");
    source_list.push_back("PG1553V6");
    source_list.push_back("PG1553V5");

    TH1D Hist_Elev = TH1D("Hist_Elev","",N_elev_bins,elev_bins);
    TH1D Hist_MJD = TH1D("Hist_MJD","",N_MJD_bins,MJD_bins);
    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);

    vector<vector<TH1D>> Hist_OnData_PerElev_MSCW;
    vector<vector<TH1D>> Hist_OnData_PerElev_MSCL;
    vector<vector<TH1D>> Hist_OnData_PerElev_Xoff;
    vector<vector<TH1D>> Hist_OnData_PerElev_Yoff;
    vector<vector<TH1D>> Hist_CRData_PerElev_Yoff;
    for (int elev=0;elev<N_elev_bins;elev++)
    {
        char elev_tag[50];
        sprintf(elev_tag, "%i", elev);
        vector<TH1D> Hist_OnData_ThisElev_MSCW;
        vector<TH1D> Hist_OnData_ThisElev_MSCL;
        vector<TH1D> Hist_OnData_ThisElev_Xoff;
        vector<TH1D> Hist_OnData_ThisElev_Yoff;
        vector<TH1D> Hist_CRData_ThisElev_Yoff;
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
            MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
            N_bins_for_deconv = N_bins_for_deconv_func_E[e];

            Hist_OnData_ThisElev_MSCW.push_back(TH1D("Hist_OnData_ThisElev_MSCW_V"+TString(elev_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_OnData_ThisElev_MSCL.push_back(TH1D("Hist_OnData_ThisElev_MSCL_V"+TString(elev_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
            Hist_OnData_ThisElev_Xoff.push_back(TH1D("Hist_OnData_ThisElev_Xoff_V"+TString(elev_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins/6,-Skymap_size,Skymap_size));
            Hist_OnData_ThisElev_Yoff.push_back(TH1D("Hist_OnData_ThisElev_Yoff_V"+TString(elev_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins/6,-Skymap_size,Skymap_size));
            Hist_CRData_ThisElev_Yoff.push_back(TH1D("Hist_CRData_ThisElev_Yoff_V"+TString(elev_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins/6,-Skymap_size,Skymap_size));
        }
        Hist_OnData_PerElev_MSCW.push_back(Hist_OnData_ThisElev_MSCW);
        Hist_OnData_PerElev_MSCL.push_back(Hist_OnData_ThisElev_MSCL);
        Hist_OnData_PerElev_Xoff.push_back(Hist_OnData_ThisElev_Xoff);
        Hist_OnData_PerElev_Yoff.push_back(Hist_OnData_ThisElev_Yoff);
        Hist_CRData_PerElev_Yoff.push_back(Hist_CRData_ThisElev_Yoff);
    }

    vector<vector<TH1D>> Hist_OnData_PerYear_MSCW;
    vector<vector<TH1D>> Hist_OnData_PerYear_MSCL;
    vector<vector<TH1D>> Hist_OnData_PerYear_Xoff;
    vector<vector<TH1D>> Hist_OnData_PerYear_Yoff;
    for (int year=0;year<N_MJD_bins;year++)
    {
        char year_tag[50];
        sprintf(year_tag, "%i", year);
        vector<TH1D> Hist_OnData_ThisYear_MSCW;
        vector<TH1D> Hist_OnData_ThisYear_MSCL;
        vector<TH1D> Hist_OnData_ThisYear_Xoff;
        vector<TH1D> Hist_OnData_ThisYear_Yoff;
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
            MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
            N_bins_for_deconv = N_bins_for_deconv_func_E[e];

            Hist_OnData_ThisYear_MSCW.push_back(TH1D("Hist_OnData_ThisYear_MSCW_V"+TString(year_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_OnData_ThisYear_MSCL.push_back(TH1D("Hist_OnData_ThisYear_MSCL_V"+TString(year_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
            Hist_OnData_ThisYear_Xoff.push_back(TH1D("Hist_OnData_ThisYear_Xoff_V"+TString(year_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins/6,-Skymap_size,Skymap_size));
            Hist_OnData_ThisYear_Yoff.push_back(TH1D("Hist_OnData_ThisYear_Yoff_V"+TString(year_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins/6,-Skymap_size,Skymap_size));
        }
        Hist_OnData_PerYear_MSCW.push_back(Hist_OnData_ThisYear_MSCW);
        Hist_OnData_PerYear_MSCL.push_back(Hist_OnData_ThisYear_MSCL);
        Hist_OnData_PerYear_Xoff.push_back(Hist_OnData_ThisYear_Xoff);
        Hist_OnData_PerYear_Yoff.push_back(Hist_OnData_ThisYear_Yoff);
    }

    vector<vector<TH1D>> Hist_Source_Theta2;
    vector<vector<TH1D>> Hist_Source_Theta2_Weighted;
    for (int source=0;source<source_list.size();source++)
    {
        vector<TH1D> Hist_ThisSource_Theta2;
        vector<TH1D> Hist_ThisSource_Theta2_Weighted;
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            Hist_ThisSource_Theta2.push_back(TH1D("Hist_Source_Theta2_"+TString(source_list[source])+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",5,0.,5.*source_theta2_cut));
            Hist_ThisSource_Theta2_Weighted.push_back(TH1D("Hist_Source_Theta2_Weighted_"+TString(source_list[source])+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",5,0.,5.*source_theta2_cut));
        }
        Hist_Source_Theta2.push_back(Hist_ThisSource_Theta2);
        Hist_Source_Theta2_Weighted.push_back(Hist_ThisSource_Theta2_Weighted);
    }

    std::cout << "Get a list of target observation runs" << std::endl;
    for (int source=0;source<source_list.size();source++)
    {
        vector<pair<string,int>> Data_runlist_init = GetRunList(source_list[source]);
        vector<pair<string,int>> Data_runlist;
        if (!TString(target).Contains("Proton")) Data_runlist = SelectONRunList(Data_runlist_init,TelElev_lower,TelElev_upper,0,360,0,0);
        else Data_runlist = Data_runlist_init;

        std::cout << "Prepare ON run samples... " << source_list[source] << std::endl;
        std::cout << "Data_runlist size = " << Data_runlist.size() << std::endl;
        for (int run=0;run<Data_runlist.size();run++)
        {
            std::cout << "Get run... " << Data_runlist[run].second << std::endl;
            char run_number[50];
            char Data_observation[50];
            sprintf(run_number, "%i", int(Data_runlist[run].second));
            sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
            string filename;
            filename = TString(TString(SMI_INPUT)+"/"+TString(run_number)+".anasum.root");

            mean_tele_point_ra = 0.;
            mean_tele_point_dec = 0.;
            pair<double,double> source_ra_dec = GetSourceRaDec(TString(Data_runlist[run].first));
            mean_tele_point_ra = source_ra_dec.first;
            mean_tele_point_dec = source_ra_dec.second;

            //std::cout << "Get telescope pointing RA and Dec..." << std::endl;
            pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
            tele_point_ra_dec = GetRunRaDec(filename,int(Data_runlist[run].second));
            run_tele_point_ra = tele_point_ra_dec.first;
            run_tele_point_dec = tele_point_ra_dec.second;

            TFile*  input_file = TFile::Open(filename.c_str());
            TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
            TTree* Data_tree = (TTree*) input_file->Get(root_file);
            SetEventDisplayTreeBranch(Data_tree);
            vector<pair<double,double>> timecut_thisrun = GetRunTimecuts(Data_runlist[run].second);

            double tele_elev = GetRunElevAzim(filename,int(Data_runlist[run].second)).first;

            Data_tree->GetEntry(0);
            double time_0 = Time;
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
                int elevation = Hist_Elev.FindBin(tele_elev)-1;
                int year = Hist_MJD.FindBin(MJD)-1;
                if (energy<0) continue;
                if (energy>=N_energy_bins) continue;
                if (elevation<0) continue;
                if (elevation>=N_elev_bins) continue;
                if (!SelectNImages()) continue;
                if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
                if (SizeSecondMax<SizeSecondMax_Cut) continue;
                if (EmissionHeight<6.) continue;
                if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
                if (MSCL<MSCL_cut_blind && MSCW<MSCW_cut_blind)
                {
                    Hist_Source_Theta2.at(source).at(energy).Fill(theta2);
                }
            }
            input_file->Close();
        }

        vector<double> source_weight;
        for (int e=0;e<N_energy_bins;e++) 
        {
            int nbins = Hist_Source_Theta2.at(source).at(e).GetNbinsX();
            double background_counts = 0.;
            for (int bin=2;bin<=nbins;bin++) // first bin is the signal region
            {
                background_counts += Hist_Source_Theta2.at(source).at(e).GetBinContent(bin);
            }
            double background_avg = background_counts/double(nbins-1);
            if (Hist_Source_Theta2.at(source).at(e).GetBinContent(1)>0.)
            {
                source_weight.push_back(background_avg/Hist_Source_Theta2.at(source).at(e).GetBinContent(1));
            }
            else
            {
                source_weight.push_back(0.);
            }
        }

        for (int run=0;run<Data_runlist.size();run++)
        {
            char run_number[50];
            char Data_observation[50];
            sprintf(run_number, "%i", int(Data_runlist[run].second));
            sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
            string filename;
            filename = TString(TString(SMI_INPUT)+"/"+TString(run_number)+".anasum.root");

            mean_tele_point_ra = 0.;
            mean_tele_point_dec = 0.;
            pair<double,double> source_ra_dec = GetSourceRaDec(TString(Data_runlist[run].first));
            mean_tele_point_ra = source_ra_dec.first;
            mean_tele_point_dec = source_ra_dec.second;

            //std::cout << "Get telescope pointing RA and Dec..." << std::endl;
            pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
            tele_point_ra_dec = GetRunRaDec(filename,int(Data_runlist[run].second));
            run_tele_point_ra = tele_point_ra_dec.first;
            run_tele_point_dec = tele_point_ra_dec.second;

            TFile*  input_file = TFile::Open(filename.c_str());
            TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
            TTree* Data_tree = (TTree*) input_file->Get(root_file);
            SetEventDisplayTreeBranch(Data_tree);
            vector<pair<double,double>> timecut_thisrun = GetRunTimecuts(Data_runlist[run].second);

            double tele_elev = GetRunElevAzim(filename,int(Data_runlist[run].second)).first;

            Data_tree->GetEntry(0);
            double time_0 = Time;
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
                int elevation = Hist_Elev.FindBin(tele_elev)-1;
                int year = Hist_MJD.FindBin(MJD)-1;
                if (energy<0) continue;
                if (energy>=N_energy_bins) continue;
                if (elevation<0) continue;
                if (elevation>=N_elev_bins) continue;
                if (!SelectNImages()) continue;
                if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
                if (SizeSecondMax<SizeSecondMax_Cut) continue;
                if (EmissionHeight<6.) continue;
                if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
                double weight = 1.;
                if (theta2<source_theta2_cut) weight = source_weight.at(energy);
                if (MSCL<MSCL_cut_blind)
                {
                    Hist_OnData_PerElev_MSCW.at(elevation).at(energy).Fill(MSCW,weight);
                    Hist_OnData_PerYear_MSCW.at(year).at(energy).Fill(MSCW,weight);
                }
                if (MSCW<MSCW_cut_blind)
                {
                    Hist_OnData_PerElev_MSCL.at(elevation).at(energy).Fill(MSCL,weight);
                    Hist_OnData_PerYear_MSCL.at(year).at(energy).Fill(MSCL,weight);
                }
                if (MSCL<MSCL_cut_blind && MSCW<MSCW_cut_blind)
                {
                    Hist_Source_Theta2_Weighted.at(source).at(energy).Fill(theta2,weight);
                    Hist_OnData_PerElev_Xoff.at(elevation).at(energy).Fill(Xoff,weight);
                    Hist_OnData_PerElev_Yoff.at(elevation).at(energy).Fill(Yoff,weight);
                    Hist_OnData_PerYear_Xoff.at(year).at(energy).Fill(Xoff,weight);
                    Hist_OnData_PerYear_Yoff.at(year).at(energy).Fill(Yoff,weight);
                }
                if (ControlSelectionTheta2())
                {
                    Hist_CRData_PerElev_Yoff.at(elevation).at(energy).Fill(Yoff);
                }
            }
            input_file->Close();
        }

    }

    TFile OutputFile(TString(SMI_OUTPUT)+"/ObservingEffect.root","recreate");
    for (int source=0;source<source_list.size();source++)
    {
        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_Source_Theta2.at(source).at(e).Write();
            Hist_Source_Theta2_Weighted.at(source).at(e).Write();
        }
    }
    for (int elev=0;elev<N_elev_bins;elev++)
    {
        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OnData_PerElev_MSCL.at(elev).at(e).Write();
            Hist_OnData_PerElev_MSCW.at(elev).at(e).Write();
            Hist_OnData_PerElev_Xoff.at(elev).at(e).Write();
            Hist_OnData_PerElev_Yoff.at(elev).at(e).Write();
            Hist_CRData_PerElev_Yoff.at(elev).at(e).Write();
        }
    }
    for (int year=0;year<N_MJD_bins;year++)
    {
        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OnData_PerYear_MSCL.at(year).at(e).Write();
            Hist_OnData_PerYear_MSCW.at(year).at(e).Write();
            Hist_OnData_PerYear_Xoff.at(year).at(e).Write();
            Hist_OnData_PerYear_Yoff.at(year).at(e).Write();
        }
    }

    OutputFile.Close();

    std::cout << "Done." << std::endl;

}
