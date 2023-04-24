
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
#include "TProfile2D.h"
#include "TGraph.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TChain.h"
#include "TBranch.h"
#include "TRandom.h"
#include "TCut.h"
#include "TMath.h"

//#include "../../EventDisplay/VEvndispRunParameter.h"
#include "/home/rshang/MatrixDecompositionMethod/EventDisplay/VEvndispRunParameter.h"

#include "GetRunList.h"
#include "ResetPublicVariables.C"

#include <complex>
//#include "../../Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
//#include "../../Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
#include "/home/rshang/MatrixDecompositionMethod/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/MatrixDecompositionMethod/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
using namespace Eigen;

double TelElev_lower = 70.;
double TelElev_upper = 80.;
char target[50] = "";
char elev_cut_tag[50] = "";
char theta2_cut_tag[50] = "";
char mjd_cut_tag[50] = "";
char group_tag[50] = "";
char map_x_tag[50] = "";
char map_y_tag[50] = "";
double SizeSecondMax_Cut = 0.;

// EventDisplay variables
double TelElevation = 0;
double TelAzimuth = 0;
double TelRAJ2000 = 0;
double TelDecJ2000 = 0;
double ErecS = 0;
double EChi2S = 0;
double MSCW = 0;
int NImages = 0;
double Xcore = 0.;
double Ycore = 0.;
double SizeSecondMax = 0;
double MSCL = 0;
double Time = 0;
int MJD = 0;
UInt_t MJD_UInt_t = 0;
double Shower_Ze = 0;
double Shower_Az = 0;
double SlantDepth = 0;
float EmissionHeight = 0;
float EmissionHeightChi2 = 0;
double Xoff = 0;
double Yoff = 0;
double Xoff_derot = 0;
double Yoff_derot = 0;
double R2off = 0;
double Phioff = 0;
double theta2 = 0;
double theta2_dark = 0;
double theta2_roi = 0;
double proj_x_roi = 0;
double proj_y_roi = 0;
double ra_sky = 0;
double dec_sky = 0;
double ra_sky_dark = 0;
double dec_sky_dark = 0;
double exposure_hours = 0.;
double exposure_hours_usable = 0.;
double exposure_hours_ref = 0.;
int MJD_Start = 2147483647;
int MJD_End = 0;

int n_expect_matches = 0;
int n_good_matches = 0;
double map_center_x = 0.;
double map_center_y = 0.;
double mean_tele_point_ra = 0.;
double mean_tele_point_dec = 0.;
double mean_tele_point_l = 0.;
double mean_tele_point_b = 0.;
double run_tele_point_ra = 0.;
double run_tele_point_dec = 0.;
vector<string> roi_name;
vector<double> roi_ra;
vector<double> roi_dec;
vector<double> roi_radius_inner;
vector<double> roi_radius_outer;
vector<vector<double>> BrightStars_Data;
vector<vector<double>> FaintStars_Data;
vector<vector<double>> GammaSource_Data;
vector<vector<double>> Dark_weight;
vector<bool> did_i_find_a_match;

string SMI_INPUT;
string SMI_OUTPUT;
string SMI_DIR;
string SMI_AUX;

double GetRunPedestalVar(int run_number)
{
    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_version = "";
    string acc_nsb = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;
    double NSB = 0.;

    ifstream myfile (SMI_AUX+"/NSB_allruns.txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==0)
                {
                    acc_runnumber += line[i];
                }
            }
            if (std::stoi(acc_runnumber,nullptr,10)==run_number)
            {
                //std::cout << "find run " << run_number << std::endl;
                nth_delimiter = 0;
                acc_nsb = "";
                for(int i = 0; i < line.size(); i++)
                {
                    if (line[i] != delimiter)
                    {
                        acc_nsb += line[i];
                    }
                    else
                    {
                        acc_nsb = "";
                        nth_delimiter += 1;
                    }
                }
                NSB = std::stod(acc_nsb,&sz);
                break;
            }
        }
        myfile.close();
    }
    else std::cout << "Unable to open file NSB_allruns.txt" << std::endl; 

    return NSB;
}

int RunTypeCategory(int run_number, bool doPrint)
{
    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_runtype = "";
    int nth_delimiter = 0;
    std::string::size_type sz;
    int runtype = 2;

    // run veritas_db_query.py to get this file
    ifstream myfile (SMI_AUX+"/runtype_allruns.txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==0)
                {
                    acc_runnumber += line[i];
                }
            }
            if (std::stoi(acc_runnumber,nullptr,10)==run_number)
            {
                //std::cout << "find run " << run_number << std::endl;
                nth_delimiter = 0;
                acc_runtype = "";
                for(int i = 0; i < line.size(); i++)
                {
                    if (line[i] != delimiter)
                    {
                        acc_runtype += line[i];
                    }
                    else
                    {
                        acc_runtype = "";
                        nth_delimiter += 1;
                    }
                }
                if (acc_runtype=="observing")
                {
                    runtype = 0;
                }
                else if (acc_runtype=="obsLowHV")
                {
                    runtype = 1;
                    if (RHVData)
                    {
                        runtype = 0;
                    }
                }
                else
                {
                    runtype = 2;
                    //if (doPrint)
                    //{
                    //    std::cout << "Run " << run_number << " rejected. Run Type " << acc_runtype << std::endl;
                    //}
                }
                break;
            }
        }
        myfile.close();
    }
    else std::cout << "Unable to open file category_allruns.txt" << std::endl; 

    return runtype;
}

double GetRunL3Rate(int run_number)
{
    return 0.;

    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_rate = "";
    int nth_delimiter = 0;
    std::string::size_type sz;
    double L3_rate = 0.;

    // run GetRunL3Rates.py to get this file
    ifstream myfile (SMI_AUX+"/L3rate_allruns.txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==0)
                {
                    acc_runnumber += line[i];
                }
            }
            if (std::stoi(acc_runnumber,nullptr,10)==run_number)
            {
                //std::cout << "find run " << run_number << std::endl;
                nth_delimiter = 0;
                acc_rate = "";
                for(int i = 0; i < line.size(); i++)
                {
                    if (line[i] != delimiter)
                    {
                        acc_rate += line[i];
                    }
                    else
                    {
                        acc_rate = "";
                        nth_delimiter += 1;
                    }
                }
                L3_rate = std::stod(acc_rate,&sz);
                break;
            }
        }
        myfile.close();
    }
    else std::cout << "Unable to open file L3rate_allruns.txt" << std::endl; 

    return L3_rate;
}

double GetRunUsableTime(string file_name,int run_number)
{
    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_time = "";
    int nth_delimiter = 0;
    std::string::size_type sz;
    double usable_time = 0.;

    // run veritas_db_query.py to get this file
    ifstream myfile (SMI_AUX+"/usable_time_allruns.txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==0)
                {
                    acc_runnumber += line[i];
                }
            }
            if (std::stoi(acc_runnumber,nullptr,10)==run_number)
            {
                //std::cout << "find run " << run_number << std::endl;
                nth_delimiter = 0;
                acc_time = "";
                for(int i = 0; i < line.size(); i++)
                {
                    if (line[i] != delimiter)
                    {
                        acc_time += line[i];
                    }
                    else
                    {
                        acc_time = "";
                        nth_delimiter += 1;
                    }
                }
                usable_time = std::stod(acc_time,&sz);
                //size_t pos = 0;
                //string time_delimiter = ":";
                //if ((pos = acc_time.find(time_delimiter)) != std::string::npos)
                //{
                //    pos = acc_time.find(time_delimiter);
                //    string time_hour = acc_time.substr(0, pos);
                //    acc_time.erase(0, pos + time_delimiter.length());
                //    pos = acc_time.find(time_delimiter);
                //    string time_minute = acc_time.substr(0, pos);
                //    acc_time.erase(0, pos + time_delimiter.length());
                //    string time_second = acc_time;
                //    usable_time = std::stod(time_hour,&sz)*60.*60.;
                //    usable_time += std::stod(time_minute,&sz)*60.;
                //    usable_time += std::stod(time_second,&sz);
                //}
                break;
            }
        }
        myfile.close();
    }
    else std::cout << "Unable to open file usable_time_allruns.txt" << std::endl; 

    //if (usable_time==0.)
    //{
    //    char run_number_char[50];
    //    sprintf(run_number_char, "%i", run_number);
    //    TFile*  input_file = TFile::Open(file_name.c_str());
    //    TTree* pointing_tree = nullptr;
    //    pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number_char)+"/stereo/pointingDataReduced"));
    //    pointing_tree->SetBranchStatus("*",0);
    //    pointing_tree->SetBranchStatus("Time",1);
    //    pointing_tree->SetBranchAddress("Time",&Time);
    //    double total_entries = (double)pointing_tree->GetEntries();
    //    pointing_tree->GetEntry(0);
    //    int time_start = Time;
    //    pointing_tree->GetEntry(total_entries-1);
    //    int time_end = Time;
    //    input_file->Close();
    //    usable_time = time_end-time_start;
    //}
    if (usable_time==0.)
    {
        usable_time = 30.*60.;
    }

    return usable_time;
}

vector<pair<double,double>> GetRunTimecuts(int run_number)
{
    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_timecuts = "";
    vector<pair<double,double>> timecuts;
    int nth_delimiter = 0;
    std::string::size_type sz;

    // run veritas_db_query.py to get this file
    ifstream myfile (SMI_AUX+"/timecuts_allruns.txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==0)
                {
                    acc_runnumber += line[i];
                }
            }
            if (std::stoi(acc_runnumber,nullptr,10)==run_number)
            {
                //std::cout << "find run " << run_number << std::endl;
                nth_delimiter = 0;
                acc_timecuts = "";
                for(int i = 0; i < line.size(); i++)
                {
                    if (line[i] != delimiter)
                    {
                        acc_timecuts += line[i];
                    }
                    else
                    {
                        size_t pos = 0;
                        string time_delimiter = "/";
                        if ((pos = acc_timecuts.find(time_delimiter)) != std::string::npos)
                        {
                            string timecut_start = acc_timecuts.substr(0, pos);
                            acc_timecuts.erase(0, pos + time_delimiter.length());
                            string timecut_end = acc_timecuts;
                            //std::cout << "timecut_start " << timecut_start << std::endl;
                            //std::cout << "timecut_end " << timecut_end << std::endl;
                            double val_timecut_start = std::stod(timecut_start,&sz);
                            double val_timecut_end = std::stod(timecut_end,&sz);
                            timecuts.push_back(std::make_pair(val_timecut_start,val_timecut_end));
                        }
                        acc_timecuts = "";
                        nth_delimiter += 1;
                    }
                }
                size_t pos = 0;
                string time_delimiter = "/";
                if ((pos = acc_timecuts.find(time_delimiter)) != std::string::npos)
                {
                    string timecut_start = acc_timecuts.substr(0, pos);
                    acc_timecuts.erase(0, pos + time_delimiter.length());
                    string timecut_end = acc_timecuts;
                    //std::cout << "timecut_start " << timecut_start << std::endl;
                    //std::cout << "timecut_end " << timecut_end << std::endl;
                    double val_timecut_start = std::stod(timecut_start,&sz);
                    double val_timecut_end = std::stod(timecut_end,&sz);
                    timecuts.push_back(std::make_pair(val_timecut_start,val_timecut_end));
                }
                break;
            }
        }
        myfile.close();
    }
    else std::cout << "Unable to open file timecuts_allruns.txt" << std::endl; 

    return timecuts;
}

int GetRunMJD(string file_name,int run)
{

    int run_mjd = 0;
    return run_mjd;

    //char run_number[50];
    //sprintf(run_number, "%i", int(run));
    //TFile*  input_file = TFile::Open(file_name.c_str());
    //TTree* pointing_tree = nullptr;
    //pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/pointingDataReduced"));
    //pointing_tree->SetBranchStatus("*",0);
    //pointing_tree->SetBranchStatus("MJD",1);
    //pointing_tree->SetBranchAddress("MJD",&MJD_UInt_t);
    //double total_entries = (double)pointing_tree->GetEntries();
    //pointing_tree->GetEntry(0);
    //run_mjd = MJD_UInt_t;
    //input_file->Close();
    ////std::cout << "root file MJD = " << run_mjd << std::endl;

    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_version = "";
    string acc_nsb = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;

    ifstream myfile (SMI_AUX+"/diagnostics.txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            acc_version = "";
            acc_nsb = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if (nth_line<84) continue;
                if(line[i] == delimiter)
                {
                    if (nth_delimiter==4 && std::stoi(acc_runnumber,nullptr,10)==run) 
                    {
                        run_mjd = std::stoi(acc_nsb,&sz);
                        if (std::stoi(acc_version,nullptr,10)!=2) run_mjd = 0;
                    }
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==0)
                {
                    acc_runnumber += line[i];
                }
                else if (nth_delimiter==1)
                {
                    acc_version += line[i];
                }
                else if (nth_delimiter==4)
                {
                    acc_nsb += line[i];
                }
            }
            nth_line += 1;
        }
        myfile.close();
    }
    else std::cout << "Unable to open file diagnostics.txt" << std::endl; 
    //std::cout << "diag file MJD = " << run_mjd << std::endl;

    return run_mjd;
}

pair<double,double> GetRunElevAzim(string file_name, int run)
{
    double TelElevation_avg = 0.;
    double TelAzimuth_avg = 0.;

    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_version = "";
    string acc_elev = "";
    string acc_az = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;
    
    if (!UseDBOnly)
    {

        if (!gSystem->AccessPathName(file_name.c_str()))
        {
            char run_number[50];
            sprintf(run_number, "%i", int(run));
            TFile*  input_file = TFile::Open(file_name.c_str());
            TTree* pointing_tree = nullptr;
            pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/pointingDataReduced"));
            pointing_tree->SetBranchStatus("*",0);
            pointing_tree->SetBranchStatus("TelElevation",1);
            pointing_tree->SetBranchStatus("TelAzimuth",1);
            pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
            pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
            double total_entries = (double)pointing_tree->GetEntries();
            pointing_tree->GetEntry(int(total_entries/2.));
            TelElevation_avg = TelElevation;
            TelAzimuth_avg = TelAzimuth;
            input_file->Close();
            std::cout << "root file elev = " << TelElevation_avg << " azim = " << TelAzimuth_avg << std::endl;
        }

    }

    return std::make_pair(TelElevation_avg,TelAzimuth_avg);
}

bool PointingSelection(string file_name,int run, double Elev_cut_lower, double Elev_cut_upper)
{
    //if (run>100000) return true;
    char run_number[50];
    sprintf(run_number, "%i", int(run));
    //TFile*  input_file = TFile::Open(file_name.c_str());

    //TTree* pointing_tree = nullptr;
    //pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/pointingDataReduced"));
    //pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
    //pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
    //pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
    //pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
    //double total_entries = (double)pointing_tree->GetEntries();
    //pointing_tree->GetEntry(int(total_entries/2.));

    TelElevation = GetRunElevAzim(file_name,run).first;
    TelAzimuth = GetRunElevAzim(file_name,run).second;
    int MJD = GetRunMJD(file_name,run);
    if (TelElevation<Elev_cut_lower) 
    {
        //input_file->Close();
        std::cout << run << ", elev = " << TelElevation << ", pointing rejected." << std::endl;
        return false;
    }
    if (TelElevation>Elev_cut_upper)
    {
        //input_file->Close();
        std::cout << run << ", elev = " << TelElevation << ", pointing rejected." << std::endl;
        return false;
    }

    //double period = 162.5;
    //double MJD_ref = 59246 + 10.*period; // 2021-02-01
    //double X = ( MJD_ref - double(MJD) )/period;
    //double Phase=X-floor(X);
    //if (TString(target).Contains("SS433Half1"))
    //{
    //    if (Phase>0.5)
    //    {
    //        return false;
    //    }
    //}
    //if (TString(target).Contains("SS433Half2"))
    //{
    //    if (Phase<=0.5)
    //    {
    //        return false;
    //    }
    //}

    //if (Azim_region=="North")
    //{
    //    if (!(TelAzimuth<45. || TelAzimuth>315.))
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (Azim_region=="South")
    //{
    //    if (!(TelAzimuth>180.-45. && TelAzimuth<180.+45.))
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (Azim_region=="East")
    //{
    //    if (!(TelAzimuth>90.-45. && TelAzimuth<90.+45.))
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (Azim_region=="West")
    //{
    //    if (!(TelAzimuth>270.-45. && TelAzimuth<270.+45.))
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (TString(target).Contains("SS433"))
    //{
    //    double TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
    //    double TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
    //    double delta_RA = TelRAJ2000_tmp - 287.9565;
    //    double delta_Dec = TelDecJ2000_tmp - 4.9827;
    //    double distance = pow(delta_RA*delta_RA + delta_Dec*delta_Dec,0.5);
    //    //if (distance>2.0)
    //    //{
    //    //    input_file->Close();
    //    //    return false;
    //    //}
    //}
    //else if (TString(target).Contains("MGRO_J2031"))
    //{
    //    double TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
    //    double TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
    //    double delta_RA = TelRAJ2000_tmp - 308.041666667;
    //    double delta_Dec = TelDecJ2000_tmp - 41.4594444444;
    //    double distance = pow(delta_RA*delta_RA + delta_Dec*delta_Dec,0.5);
    //    if (distance>2.0)
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //else if (TString(target).Contains("MGRO_J1908"))
    //{
    //    double TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
    //    double TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
    //    double delta_RA = TelRAJ2000_tmp - 287.188333333;
    //    double delta_Dec = TelDecJ2000_tmp - 6.16233333333;
    //    double distance = pow(delta_RA*delta_RA + delta_Dec*delta_Dec,0.5);
    //    if (distance>2.0)
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //else if (TString(target).Contains("WComae"))
    //{
    //    double TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
    //    double TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
    //    //double delta_RA = TelRAJ2000_tmp - 185.382;
    //    //double delta_Dec = TelDecJ2000_tmp - 28.233;
    //    double delta_RA = TelRAJ2000_tmp - (185.360+184.616)/2.;
    //    double delta_Dec = TelDecJ2000_tmp - (30.191+30.130)/2.;
    //    double distance = pow(delta_RA*delta_RA + delta_Dec*delta_Dec,0.5);
    //    if (distance>2.0)
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //    //if (MJD>57250)
    //    //{
    //    //    input_file->Close();
    //    //    return false;
    //    //}
    //}
    //int MJD_start = 59246; // 2021-02-01
    //int MJD_end = 59327; // 2021-04-23
    //double X = ( double(MJD_start) - double(MJD) )/162.5;
    //double Phase=X-floor(X);
    //if (TString(target).Contains("SS433Half1"))
    //{
    //    if (Phase>0.5)
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (TString(target).Contains("SS433Half2"))
    //{
    //    if (Phase<=0.5)
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (TString(target).Contains("MGRO_J1908Half1"))
    //{
    //    if (Phase>0.5)
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (TString(target).Contains("MGRO_J1908Half2"))
    //{
    //    if (Phase<=0.5)
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (TString(target).Contains("SS433Quad1"))
    //{
    //    if (!(Phase<0.25))
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (TString(target).Contains("SS433Quad2"))
    //{
    //    if (!(Phase>=0.25 && Phase<0.5))
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (TString(target).Contains("SS433Quad3"))
    //{
    //    if (!(Phase>=0.5 && Phase<0.75))
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}
    //if (TString(target).Contains("SS433Quad4"))
    //{
    //    if (!(Phase>=0.75))
    //    {
    //        input_file->Close();
    //        return false;
    //    }
    //}

    //input_file->Close();
    return true;
}

pair<double,double> ConvertRaDecToGalactic(double RA, double Dec)
{
    double delta = Dec*M_PI/180.;
    double delta_G = 27.12825*M_PI/180.;
    double alpha = RA*M_PI/180.;
    double alpha_G = 192.85948*M_PI/180.;
    double l_NCP = 122.93192*M_PI/180.;
    double sin_b = TMath::Sin(delta)*TMath::Sin(delta_G)+TMath::Cos(delta)*TMath::Cos(delta_G)*TMath::Cos(alpha-alpha_G);
    double cos_b = TMath::Cos(TMath::ASin(sin_b));
    double sin_l_NCP_m_l = TMath::Cos(delta)*TMath::Sin(alpha-alpha_G)/cos_b;
    double cos_l_NCP_m_l = (TMath::Cos(delta_G)*TMath::Sin(delta)-TMath::Sin(delta_G)*TMath::Cos(delta)*TMath::Cos(alpha-alpha_G))/cos_b;
    double b = (TMath::ASin(sin_b))*180./M_PI;
    double l = (l_NCP-TMath::ATan2(sin_l_NCP_m_l,cos_l_NCP_m_l))*180./M_PI;
    return std::make_pair(l,b);
}

pair<double,double> GetRunRaDec(string file_name, int run)
{

    double TelRAJ2000_tmp = 0.;
    double TelDecJ2000_tmp = 0.;

    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_version = "";
    string acc_elev = "";
    string acc_az = "";
    string acc_ra = "";
    string acc_dec = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;

    if (!UseDBOnly)
    {
        if (!gSystem->AccessPathName(file_name.c_str()))
        {
            char run_number[50];
            sprintf(run_number, "%i", int(run));
            TFile*  input_file = TFile::Open(file_name.c_str());
            TTree* pointing_tree = nullptr;
            pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/pointingDataReduced"));
            pointing_tree->SetBranchStatus("*",0);
            pointing_tree->SetBranchStatus("TelRAJ2000",1);
            pointing_tree->SetBranchStatus("TelDecJ2000",1);
            pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
            pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
            double total_entries = (double)pointing_tree->GetEntries();
            pointing_tree->GetEntry(int(total_entries/2.));
            TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
            TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
            input_file->Close();

            std::cout << "root file RA = " << TelRAJ2000_tmp << " Dec = " << TelDecJ2000_tmp << std::endl;
        }
    }

    return std::make_pair(TelRAJ2000_tmp,TelDecJ2000_tmp);
}

bool MJDSelection(string file_name,int run)
{
    return true;

    if (MJD_start_cut==0 && MJD_end_cut==0) return true;
    //if (run>100000) return true;
    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());

    //TTree* pointing_tree = nullptr;
    //pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/pointingDataReduced"));
    //pointing_tree->SetBranchAddress("MJD",&MJD_UInt_t);
    //double total_entries = (double)pointing_tree->GetEntries();
    //pointing_tree->GetEntry(0);

    MJD_UInt_t = GetRunMJD(file_name,run);
    if (MJD_UInt_t<MJD_start_cut || MJD_UInt_t>MJD_end_cut) 
    {
        input_file->Close();
        return false;
    }
    input_file->Close();
    return true;
}
vector<double> GetRunElevationList(vector<pair<string,int>> list)
{
    vector<double> list_elev;
    for (int run=0;run<list.size();run++)
    {
        double run_elev = GetRunElevAzim(list.at(run).first,int(list.at(run).second)).first;
        list_elev.push_back(run_elev);
    }
    return list_elev;
}
vector<double> GetRunRaDecList(vector<pair<string,int>> list)
{
    vector<double> list_radec;
    for (int run=0;run<list.size();run++)
    {
        double run_RA_Dec = GetRunRaDec(list.at(run).first,int(list.at(run).second)).second;
        list_radec.push_back(run_RA_Dec);
    }
    return list_radec;
}
void SortingList(vector<pair<string,int>>* list, vector<double>* list_pointing)
{
    pair<string,int> temp_list;
    double temp_list_pointing;
    for (int run1=0;run1<list->size()-1;run1++)
    {
        for (int run2=run1+1;run2<list->size();run2++)
        {
            if (list_pointing->at(run1)>list_pointing->at(run2))
            {
                temp_list = list->at(run1);
                temp_list_pointing = list_pointing->at(run1);
                list->at(run1) = list->at(run2);
                list_pointing->at(run1) = list_pointing->at(run2);
                list->at(run2) = temp_list;
                list_pointing->at(run2) = temp_list_pointing;
            }
        }
    }
}

int FindAMatchedRunFromPairList(int ON_runnumber, string source, vector<int> exclusion_list, bool isImposter, int true_ON_runnumber, double elev_offset) 
{

    string line;
    char delimiter = ' ';
    string txt_on_runnumber = "";
    string txt_off_runnumber = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;
    vector<int> offrun_list;

    string SMI_DIR;
    SMI_DIR = string(std::getenv("SMI_DIR"));
    ifstream myfile (SMI_DIR+"/runlist_backup/PairList_"+source+".txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            txt_on_runnumber = "";
            txt_off_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==1)
                {
                    txt_off_runnumber += line[i];
                }
                else if (nth_delimiter==0)
                {
                    txt_on_runnumber += line[i];
                }
                if (i==line.size()-1)
                {
                    double this_on_runnumber = std::stod(txt_on_runnumber,&sz);
                    double this_off_runnumber = std::stod(txt_off_runnumber,&sz);
                    if (!isImposter)
                    {
                        if (int(this_on_runnumber)==ON_runnumber)
                        {
                            offrun_list.push_back(int(this_off_runnumber));
                        }
                    }
                    else
                    {
                        if (int(this_on_runnumber)==true_ON_runnumber)
                        {
                            offrun_list.push_back(int(this_off_runnumber));
                        }
                    }
                }
            }
            nth_line += 1;
        }
        myfile.close();
    }

    if (offrun_list.size()==0)
    {
        std::cout << "ON run " << ON_runnumber << " does not have matched off runs." << std::endl;
        return 0;
    }

    char char_onrun_number[50];
    sprintf(char_onrun_number, "%i", ON_runnumber);
    string ON_filename;
    ON_filename = TString(SMI_INPUT+"/"+string(char_onrun_number)+".anasum.root");
    std::pair<double,double> onrun_pointing = GetRunElevAzim(ON_filename,ON_runnumber);

    double chi2 = 1e10;
    int best_off_runnumber = 0;
    for (int off_run=0;off_run<offrun_list.size();off_run++)
    {
        int OFF_runnumber = offrun_list[off_run];

        bool do_not_use = false;
        for (int exclude_run=0;exclude_run<exclusion_list.size();exclude_run++)
        {
            if (OFF_runnumber==exclusion_list[exclude_run])
            {
                do_not_use = true;
            }
        }
        if (do_not_use) continue;

        char char_offrun_number[50];
        sprintf(char_offrun_number, "%i", OFF_runnumber);
        string OFF_filename;
        OFF_filename = TString(SMI_INPUT+"/"+string(char_offrun_number)+".anasum.root");
        std::pair<double,double> offrun_pointing = GetRunElevAzim(OFF_filename,OFF_runnumber);

        double delta_elev = onrun_pointing.first-offrun_pointing.first;
        double delta_airmass = (1./sin(onrun_pointing.first*M_PI/180.)-1./sin(offrun_pointing.first*M_PI/180.));
        double delta_azim = onrun_pointing.second-offrun_pointing.second;
        if (delta_azim>180.) delta_azim = 360.-delta_azim;

        if (!isImposter)
        {
            if (elev_offset>0.)
            {
                if (delta_elev>0.) continue;
            }
            else
            {
                if (delta_elev<0.) continue;
            }
        }
        if (abs(delta_airmass)>MatchRun_dElev) continue;
        if (abs(delta_azim)>MatchRun_dAzim) continue;

        double this_chi2 = pow(delta_elev,2);
        if (this_chi2<chi2)
        {
            chi2 = this_chi2;
            best_off_runnumber = OFF_runnumber;
        }
    }

    return best_off_runnumber;

}

int FindAMatchedRun(int ON_runnumber, pair<double,double> ON_pointing, double ON_NSB, double ON_L3Rate, double ON_MJD, vector<int> OFF_runnumber, vector<pair<double,double>> OFF_pointing, vector<double> OFF_NSB, vector<double> OFF_L3Rate, vector<double> OFF_MJD, vector<int> exclusion_list, bool isImposter, double elev_offset) 
{
    int matched_runnumber = 0;
    double match_chi2 = 1e10;
    double threshold_dElev = 2.0;
    double threshold_dNSB = 0.4; // for J1908 analysis
    double threshold_dAzim = 180.;
    double threshold_dMJD = 2.*365.;
    double threshold_dL3Rate = 50.;

    threshold_dElev = 10.;
    threshold_dAzim = 45.; // default, do not use for J1908
    threshold_dNSB = 1.0; // default, do not use for J1908
    threshold_dMJD = 1000.*365.;
    threshold_dL3Rate = 10000.;

    if (RHVData)
    {
        threshold_dElev = 2.*MatchRun_dElev; // default 
        threshold_dAzim = 45.; // default 
        threshold_dNSB = 10.; // default 
        threshold_dMJD = 100.*365.;
    }
    else
    {
        threshold_dElev = 2.*MatchRun_dElev; // default 
        threshold_dAzim = 2.*MatchRun_dAzim; // default 
        threshold_dNSB = 2.*MatchRun_dNSB; // default 
        threshold_dMJD = 100.*365.;
    }

    if (!isImposter)
    {
        threshold_dElev = MatchRun_dElev; // default 
        threshold_dAzim = MatchRun_dAzim; // default 
        threshold_dNSB = MatchRun_dNSB; // default 
        threshold_dMJD = 100.*365.;

        if (MatchingSelection==1)
        {
            threshold_dElev = 1000.;
        }
        if (MatchingSelection==2)
        {
            threshold_dAzim = 1000.;
        }
        if (MatchingSelection==3)
        {
            threshold_dNSB = 1000.;
        }
        if (MatchingSelection==4)
        {
            threshold_dMJD = 100.*365.;
        }
        if (MatchingSelection==5)
        {
            threshold_dL3Rate = 1000.;
        }
    }

    for (int off_run=0;off_run<OFF_runnumber.size();off_run++)
    {
        bool do_not_use = false;
        for (int exclude_run=0;exclude_run<exclusion_list.size();exclude_run++)
        {
            if (OFF_runnumber[off_run]==exclusion_list[exclude_run])
            {
                do_not_use = true;
            }
        }
        if (do_not_use) continue;
        double delta_l3rate = abs(ON_L3Rate-OFF_L3Rate[off_run]);
        double delta_mjd = abs(ON_MJD-OFF_MJD[off_run]);
        double delta_nsb = abs(ON_NSB-OFF_NSB[off_run]);
        double delta_elev = (ON_pointing.first-OFF_pointing[off_run].first);
        double delta_airmass = (1./sin(ON_pointing.first*M_PI/180.)-1./sin(OFF_pointing[off_run].first*M_PI/180.));
        double signed_delta_airmass = (1./sin(ON_pointing.first*M_PI/180.)-1./sin(OFF_pointing[off_run].first*M_PI/180.));

        double delta_azim = abs(ON_pointing.second-OFF_pointing[off_run].second);
        if (delta_azim>180.) delta_azim = 360.-delta_azim;

        if (!isImposter)
        {
            //if (ON_runnumber % 2 == 1)
            if (elev_offset>0.)
            {
                if (delta_elev>0.) continue;
            }
            else
            {
                if (delta_elev<0.) continue;
            }
        }

        if (abs(delta_airmass)>threshold_dElev) continue;
        //if (abs(delta_elev)>threshold_dElev) continue;
        if (delta_nsb>threshold_dNSB) continue;
        if (delta_azim>threshold_dAzim) continue; //default, do not use for J1908

        if (delta_mjd>threshold_dMJD) continue;
        if (delta_l3rate>threshold_dL3Rate) continue;

        double chi2 = pow(delta_azim,2); // for J1908 analysis

        //chi2 = pow(delta_airmass,2);  // default, do not use for J1908
        chi2 = pow(delta_elev,2);
        //chi2 = pow(2.*delta_elev,2)+pow(delta_azim,2);
        //chi2 = pow(delta_mjd,2);

        if (chi2<match_chi2)
        {
            match_chi2 = chi2;
            matched_runnumber = OFF_runnumber[off_run];
        }
    }
    return matched_runnumber;
}

bool CheckBrightGammaRaySource(double input_ra, double input_dec)
{
    double delta_ra = 0.;
    double delta_dec = 0.;
    double delta_angle = 0.;

    // geminga
    delta_ra = 98.476-input_ra;
    delta_dec = 17.770-input_dec;
    delta_angle = pow(delta_ra*delta_ra+delta_dec*delta_dec,0.5);
    if (delta_angle<5.)
    {
        return true;
    }
    // mrk 421
    delta_ra = 166.079-input_ra;
    delta_dec = 38.195-input_dec;
    delta_angle = pow(delta_ra*delta_ra+delta_dec*delta_dec,0.5);
    if (delta_angle<3.)
    {
        return true;
    }
    // 1es 1218
    delta_ra = 185.360-input_ra;
    delta_dec = 30.191-input_dec;
    delta_angle = pow(delta_ra*delta_ra+delta_dec*delta_dec,0.5);
    if (delta_angle<5.)
    {
        return true;
    }

    return false;
}

void SelectImposterRunList(TTree * RunListTree, TString root_file_name, vector<pair<string,int>> Data_runlist, vector<pair<string,int>> OFF_runlist_input, double Elev_cut_lower, double Elev_cut_upper, int iteration, pair<double,double> source_ra_dec, string source_strip)
{
    std::cout << "initial runs = " << Data_runlist.size() << std::endl;

    vector<int> new_list;
    vector<pair<double,double>> new_list_pointing;
    vector<double> ON_time;
    vector<double> ON_NSB;
    vector<double> ON_L3Rate;
    for (int run=0;run<Data_runlist.size();run++)
    {

        if (RunTypeCategory(Data_runlist[run].second,true)!=0) 
        {
            std::cout << int(Data_runlist[run].second) << " category rejected." << std::endl;
            continue;
        }
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        double NSB_thisrun = GetRunPedestalVar(int(Data_runlist[run].second));
        string filename;
        filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

        if (!PointingSelection(filename,int(Data_runlist[run].second),Elev_cut_lower,Elev_cut_upper))
        {
            //std::cout << int(Data_runlist[run].second) << " pointing rejected." << std::endl;
            continue;
        }
        if (!MJDSelection(filename,int(Data_runlist[run].second))) 
        {
            //std::cout << int(Data_runlist[run].second) << " MJD rejected." << std::endl;
            continue;
        }
        double L3_rate = GetRunL3Rate(Data_runlist[run].second);
        //if (NSB_thisrun>7.0)
        //{
        //    std::cout << int(Data_runlist[run].second) << " NSB rejected." << std::endl;
        //    continue;
        //}
        if (UseGalacticCoord)
        {
            std::pair<double,double> on_run_RA_Dec = GetRunRaDec(filename,int(Data_runlist[run].second));
            pair<double,double> run_galactic_coords =  ConvertRaDecToGalactic(on_run_RA_Dec.first,on_run_RA_Dec.second);
            pair<double,double> mean_galactic_coords =  ConvertRaDecToGalactic(source_ra_dec.first,source_ra_dec.second);
            if (abs(run_galactic_coords.first-mean_galactic_coords.first)>Skymap_size_x+2.)
            {
                std::cout << int(Data_runlist[run].second) << " sky coord. rejected." << std::endl;
                continue;
            }
        }
        else
        {
            std::pair<double,double> on_run_RA_Dec = GetRunRaDec(filename,int(Data_runlist[run].second));
            if (abs(on_run_RA_Dec.first-source_ra_dec.first)>Skymap_size_x+1.)
            {
                std::cout << int(Data_runlist[run].second) << " sky coord. rejected." << std::endl;
                continue;
            }
            if (abs(on_run_RA_Dec.second-source_ra_dec.second)>Skymap_size_y+1.)
            {
                std::cout << int(Data_runlist[run].second) << " sky coord. rejected." << std::endl;
                continue;
            }
        }

        new_list.push_back(Data_runlist[run].second);
        new_list_pointing.push_back(GetRunElevAzim(filename,int(Data_runlist[run].second)));
        ON_NSB.push_back(NSB_thisrun);
        double exposure_thisrun = GetRunUsableTime(filename,Data_runlist[run].second);
        ON_time.push_back(exposure_thisrun);
        ON_L3Rate.push_back(GetRunL3Rate(Data_runlist[run].second));

    }
    std::cout << "ON runs = " << new_list.size() << std::endl;
    if (new_list.size()==0)
    {
        return;
    }

    std::cout << "Prepare OFF run list..." << std::endl;
    std::cout << "Initial OFF run list = " << OFF_runlist_input.size() << std::endl;
    vector<int> OFF_runlist;
    vector<pair<double,double>> OFF_pointing;
    vector<pair<double,double>> OFF_pointing_radec;
    vector<double> OFF_time;
    vector<double> OFF_L3Rate;
    vector<double> OFF_NSB;
    vector<double> OFF_MJD;
    for (int off_run=0;off_run<OFF_runlist_input.size();off_run++)
    {
        //std::cout << "Load OFF run " << OFF_runlist_input[off_run].second << std::endl;
        char OFF_runnumber[50];
        char OFF_observation[50];
        sprintf(OFF_runnumber, "%i", int(OFF_runlist_input[off_run].second));
        sprintf(OFF_observation, "%s", OFF_runlist_input[off_run].first.c_str());
        string OFF_filename;
        OFF_filename = TString(SMI_INPUT+"/"+string(OFF_runnumber)+".anasum.root");
        std::pair<double,double> run_pointing = GetRunElevAzim(OFF_filename,int(OFF_runlist_input[off_run].second));
        double run_elevation = run_pointing.first;
        if (run_elevation<Elev_cut_lower-5.) continue;
        if (run_elevation>Elev_cut_upper+5.) continue;
        if (RunTypeCategory(OFF_runlist_input[off_run].second,true)!=0) 
        {
            //std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " category rejected." << std::endl;
            continue;
        }
        double L3_rate = GetRunL3Rate(OFF_runlist_input[off_run].second);
        std::pair<double,double> off_run_RA_Dec = GetRunRaDec(OFF_filename,int(OFF_runlist_input[off_run].second));
        double delta_ra = source_ra_dec.first-off_run_RA_Dec.first;
        double delta_dec = source_ra_dec.second-off_run_RA_Dec.second;
        double delta_angle = pow(delta_ra*delta_ra+delta_dec*delta_dec,0.5);
        if (delta_angle<5.)
        {
            //std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " coordinate rejected." << std::endl;
            continue;
        }
        if (CheckBrightGammaRaySource(off_run_RA_Dec.first,off_run_RA_Dec.second))
        {
            //std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " coordinate rejected." << std::endl;
            continue;
        }
        std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " selected." << std::endl;
        OFF_runlist.push_back(OFF_runlist_input[off_run].second);
        OFF_pointing.push_back(run_pointing);
        OFF_pointing_radec.push_back(off_run_RA_Dec);
        double NSB_thisrun = GetRunPedestalVar(int(OFF_runlist_input[off_run].second));
        OFF_NSB.push_back(NSB_thisrun);
        double exposure_thisrun = GetRunUsableTime(OFF_filename,OFF_runlist_input[off_run].second);
        OFF_time.push_back(exposure_thisrun);
        OFF_L3Rate.push_back(L3_rate);
        double MJD_thisrun = double(GetRunMJD(OFF_filename,OFF_runlist_input[off_run].second));
        OFF_MJD.push_back(MJD_thisrun);

    }


    int ON_runnumber = 0;
    double ON_pointing_RA = 0.;
    double ON_pointing_Dec = 0.;
    double ON_exposure_hour = 0.;
    vector<pair<double,double>> ON_timecut;
    vector<int> OFF_runnumber;
    vector<double> OFF_exposure_hour;

    vector<int> exclusion_list = new_list;
    for (int off_run=0;off_run<OFF_runlist.size();off_run++)
    {
        double dist_crab = pow(OFF_pointing_radec.at(off_run).first-83.633,2)+pow(OFF_pointing_radec.at(off_run).second-22.014,2);
        double dist_mrk421 = pow(OFF_pointing_radec.at(off_run).first-166.079,2)+pow(OFF_pointing_radec.at(off_run).second-38.195,2);
        double dist_mrk501 = pow(OFF_pointing_radec.at(off_run).first-253.467,2)+pow(OFF_pointing_radec.at(off_run).second-39.76,2);
        double dist_geminga = pow(OFF_pointing_radec.at(off_run).first-98.117,2)+pow(OFF_pointing_radec.at(off_run).second-17.367,2);
        if (dist_crab<2.*2.)
        {
            exclusion_list.push_back(OFF_runlist.at(off_run));
        }
        if (dist_mrk421<2.*2.)
        {
            exclusion_list.push_back(OFF_runlist.at(off_run));
        }
        if (dist_mrk501<2.*2.)
        {
            exclusion_list.push_back(OFF_runlist.at(off_run));
        }
        if (dist_geminga<2.*2.)
        {
            exclusion_list.push_back(OFF_runlist.at(off_run));
        }
    }
    vector<std::pair<int,int>> imposter_list;
    vector<pair<double,double>> imposter_list_radec;
    std::cout << "Find matched imposter run list..." << std::endl;
    std::cout << "ON runs = " << new_list.size() << std::endl;
    for (int trial=0;trial<iteration;trial++)
    {
        if (root_file_name.Contains("MGRO_J1908")) continue;
        for (int on_run=0;on_run<new_list.size();on_run++)
        {
            ON_runnumber = new_list.at(on_run);
            ON_timecut = GetRunTimecuts(ON_runnumber);
            std::cout << "Selected ON run " << ON_runnumber << std::endl;
            std::cout << "time cut: ";
            for (int entry=0;entry<ON_timecut.size();entry++)
            {
                std::cout << ON_timecut.at(entry).first << "-" << ON_timecut.at(entry).second << ", ";
            }
            std::cout << std::endl;

            char char_runnumber[50];
            sprintf(char_runnumber, "%i", ON_runnumber);
            string ON_filename;
            ON_filename = TString(SMI_INPUT+"/"+string(char_runnumber)+".anasum.root");
            std::pair<double,double> on_run_elev_azim = GetRunElevAzim(ON_filename,ON_runnumber);
            std::pair<double,double> on_run_RA_Dec = GetRunRaDec(ON_filename,ON_runnumber);
            ON_pointing_RA = on_run_RA_Dec.first;
            ON_pointing_Dec = on_run_RA_Dec.second;
            double on_run_NSB = GetRunPedestalVar(ON_runnumber);
            ON_exposure_hour = GetRunUsableTime(ON_filename,ON_runnumber)/3600.;
            double ON_L3Rate = GetRunL3Rate(ON_runnumber);
            pair<double,double> ON_pointing = new_list_pointing[on_run];
            double on_run_MJD = double(GetRunMJD(ON_filename,ON_runnumber));

            int matched_runnumber = FindAMatchedRunFromPairList(ON_runnumber, source_strip, exclusion_list, false, 0, 0.);
            if (matched_runnumber==0)
            {
                matched_runnumber = FindAMatchedRun(ON_runnumber, ON_pointing, on_run_NSB, ON_L3Rate, on_run_MJD, OFF_runlist, OFF_pointing, OFF_NSB, OFF_L3Rate, OFF_MJD, exclusion_list, true, 0.); 
            }
            if (matched_runnumber==0)
            {
                std::cout << "ON run " << ON_runnumber << " failed to find an imposter." << std::endl;
            }
            else
            {
                std::cout << "ON run " << ON_runnumber << " found an imposter " << matched_runnumber << std::endl;
                if (trial==iteration-1)
                {
                    imposter_list.push_back(std::make_pair(matched_runnumber,ON_runnumber));
                    imposter_list_radec.push_back(std::make_pair(ON_pointing_RA,ON_pointing_Dec));
                }
                exclusion_list.push_back(matched_runnumber);
            }
        }
    }


    RunListTree->Branch("ON_runnumber",&ON_runnumber,"ON_runnumber/I");
    RunListTree->Branch("ON_exposure_hour",&ON_exposure_hour,"ON_exposure_hour/D");
    RunListTree->Branch("ON_pointing_RA",&ON_pointing_RA,"ON_pointing_RA/D");
    RunListTree->Branch("ON_pointing_Dec",&ON_pointing_Dec,"ON_pointing_Dec/D");
    RunListTree->Branch("ON_timecut","std::vector<std::pair<double,double>>",&ON_timecut);
    RunListTree->Branch("OFF_runnumber","std::vector<int>",&OFF_runnumber);
    RunListTree->Branch("OFF_exposure_hour","std::vector<double>",&OFF_exposure_hour);


    exclusion_list = new_list;
    for (int on_run=0;on_run<imposter_list.size();on_run++)
    {
        exclusion_list.push_back(imposter_list.at(on_run).first);
    }
    for (int on_run=0;on_run<imposter_list.size();on_run++)
    {

        ON_runnumber = imposter_list.at(on_run).first;
        int true_ON_runnumber = imposter_list.at(on_run).second;
        ON_timecut = GetRunTimecuts(ON_runnumber);
        std::cout << "Selected imposter run " << ON_runnumber << std::endl;
        std::cout << "time cut: ";
        for (int entry=0;entry<ON_timecut.size();entry++)
        {
            std::cout << ON_timecut.at(entry).first << "-" << ON_timecut.at(entry).second << ", ";
        }
        std::cout << std::endl;

        char char_runnumber[50];
        sprintf(char_runnumber, "%i", ON_runnumber);
        string ON_filename;
        ON_filename = TString(SMI_INPUT+"/"+string(char_runnumber)+".anasum.root");
        std::pair<double,double> on_run_elev_azim = GetRunElevAzim(ON_filename,ON_runnumber);
        std::pair<double,double> on_run_RA_Dec = imposter_list_radec.at(on_run);
        ON_pointing_RA = on_run_RA_Dec.first;
        ON_pointing_Dec = on_run_RA_Dec.second;
        double on_run_NSB = GetRunPedestalVar(ON_runnumber);
        ON_exposure_hour = GetRunUsableTime(ON_filename,ON_runnumber)/3600.;
        double ON_L3Rate = GetRunL3Rate(ON_runnumber);
        pair<double,double> ON_pointing = GetRunElevAzim(ON_filename,ON_runnumber);
        double on_run_MJD = double(GetRunMJD(ON_filename,ON_runnumber));

        double total_off_time = 0.;
        OFF_runnumber.clear();
        OFF_exposure_hour.clear();
        bool continue_to_find_match = true;
        double sum_on_elev = 0.;
        double sum_off_elev = 0.;

        int matched_runnumber = FindAMatchedRunFromPairList(ON_runnumber, source_strip, exclusion_list, true, true_ON_runnumber, sum_on_elev-sum_off_elev);
        if (matched_runnumber==0)
        {
            matched_runnumber = FindAMatchedRun(ON_runnumber, ON_pointing, on_run_NSB, ON_L3Rate, on_run_MJD, OFF_runlist, OFF_pointing, OFF_NSB, OFF_L3Rate, OFF_MJD, exclusion_list, false, sum_on_elev-sum_off_elev); 
        }
        if (matched_runnumber==0)
        {
            std::cout << "Imposter run " << ON_runnumber << " failed to find a match." << std::endl;
            ON_runnumber = 0;
            OFF_runnumber.push_back(0);
            OFF_exposure_hour.push_back(0);
            continue_to_find_match = false;
        }
        else
        {
            std::cout << "Imposter run " << ON_runnumber << " found a match " << matched_runnumber << std::endl;
            exclusion_list.push_back(matched_runnumber);
            string OFF_filename;
            char char_off_runnumber[50];
            sprintf(char_off_runnumber, "%i", matched_runnumber);
            OFF_filename = TString(SMI_INPUT+"/"+string(char_off_runnumber)+".anasum.root");
            OFF_runnumber.push_back(matched_runnumber);
            double off_time = GetRunUsableTime(OFF_filename,matched_runnumber)/3600.;
            OFF_exposure_hour.push_back(off_time);
            total_off_time += off_time;
            std::pair<double,double> off_run_pointing = GetRunElevAzim(OFF_filename,matched_runnumber);
            sum_on_elev += ON_pointing.first*off_time;
            sum_off_elev += off_run_pointing.first*off_time;
        }

        if (ON_runnumber==0) continue;
        RunListTree->Fill();
    }

}

void SelectONRunList(TTree * RunListTree, vector<pair<string,int>> Data_runlist, vector<pair<string,int>> OFF_runlist_input, double Elev_cut_lower, double Elev_cut_upper, pair<double,double> source_ra_dec, string source_strip)
{
    std::cout << "initial runs = " << Data_runlist.size() << std::endl;

    vector<int> new_list;
    vector<pair<double,double>> new_list_pointing;
    vector<double> ON_time;
    vector<double> ON_NSB;
    vector<double> ON_L3Rate;
    std::cout << "ON runs (initial) = " << Data_runlist.size() << std::endl;
    for (int run=0;run<Data_runlist.size();run++)
    {

        //if (RunTypeCategory(Data_runlist[run].second,true)!=0) 
        //{
        //    std::cout << int(Data_runlist[run].second) << " category rejected." << std::endl;
        //    continue;
        //}
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        double NSB_thisrun = GetRunPedestalVar(int(Data_runlist[run].second));
        string filename;
        filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

        if (!PointingSelection(filename,int(Data_runlist[run].second),Elev_cut_lower,Elev_cut_upper))
        {
            //std::cout << int(Data_runlist[run].second) << " pointing rejected." << std::endl;
            continue;
        }
        if (!MJDSelection(filename,int(Data_runlist[run].second))) 
        {
            std::cout << int(Data_runlist[run].second) << " MJD rejected." << std::endl;
            continue;
        }
        double L3_rate = GetRunL3Rate(Data_runlist[run].second);
        //if (NSB_thisrun>7.0)
        //{
        //    std::cout << int(Data_runlist[run].second) << " NSB rejected." << std::endl;
        //    continue;
        //}
        if (UseGalacticCoord)
        {
            std::pair<double,double> on_run_RA_Dec = GetRunRaDec(filename,int(Data_runlist[run].second));
            pair<double,double> run_galactic_coords =  ConvertRaDecToGalactic(on_run_RA_Dec.first,on_run_RA_Dec.second);
            pair<double,double> mean_galactic_coords =  ConvertRaDecToGalactic(source_ra_dec.first,source_ra_dec.second);
            if (abs(run_galactic_coords.first-mean_galactic_coords.first)>Skymap_size_x+0.5)
            {
                std::cout << int(Data_runlist[run].second) << " sky coord. rejected." << std::endl;
                continue;
            }
        }
        else
        {
            std::pair<double,double> on_run_RA_Dec = GetRunRaDec(filename,int(Data_runlist[run].second));
            if (abs(on_run_RA_Dec.first-source_ra_dec.first)>Skymap_size_x+0.5)
            {
                std::cout << int(Data_runlist[run].second) << " sky coord. rejected." << std::endl;
                continue;
            }
            if (abs(on_run_RA_Dec.second-source_ra_dec.second)>Skymap_size_y+0.5)
            {
                std::cout << int(Data_runlist[run].second) << " sky coord. rejected." << std::endl;
                continue;
            }
        }

        new_list.push_back(Data_runlist[run].second);
        new_list_pointing.push_back(GetRunElevAzim(filename,int(Data_runlist[run].second)));
        ON_NSB.push_back(NSB_thisrun);
        double exposure_thisrun = GetRunUsableTime(filename,Data_runlist[run].second);
        ON_time.push_back(exposure_thisrun);
        ON_L3Rate.push_back(GetRunL3Rate(Data_runlist[run].second));

    }
    std::cout << "ON runs (final) = " << new_list.size() << std::endl;
    if (new_list.size()==0)
    {
        return;
    }

    std::cout << "Prepare OFF run list..." << std::endl;
    std::cout << "Initial OFF run list = " << OFF_runlist_input.size() << std::endl;
    vector<int> OFF_runlist;
    vector<pair<double,double>> OFF_pointing;
    vector<pair<double,double>> OFF_pointing_radec;
    vector<double> OFF_time;
    vector<double> OFF_L3Rate;
    vector<double> OFF_NSB;
    vector<double> OFF_MJD;
    for (int off_run=0;off_run<OFF_runlist_input.size();off_run++)
    {
        //std::cout << "Load OFF run " << OFF_runlist_input[off_run].second << std::endl;
        char OFF_runnumber[50];
        char OFF_observation[50];
        sprintf(OFF_runnumber, "%i", int(OFF_runlist_input[off_run].second));
        sprintf(OFF_observation, "%s", OFF_runlist_input[off_run].first.c_str());
        string OFF_filename;
        OFF_filename = TString(SMI_INPUT+"/"+string(OFF_runnumber)+".anasum.root");
        std::pair<double,double> run_pointing = GetRunElevAzim(OFF_filename,int(OFF_runlist_input[off_run].second));
        double run_elevation = run_pointing.first;
        if (run_elevation<Elev_cut_lower-5.) continue;
        if (run_elevation>Elev_cut_upper+5.) continue;
        if (RunTypeCategory(OFF_runlist_input[off_run].second,true)!=0) 
        {
            //std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " category rejected." << std::endl;
            continue;
        }
        double L3_rate = GetRunL3Rate(OFF_runlist_input[off_run].second);
        std::pair<double,double> off_run_RA_Dec = GetRunRaDec(OFF_filename,int(OFF_runlist_input[off_run].second));
        double delta_ra = source_ra_dec.first-off_run_RA_Dec.first;
        double delta_dec = source_ra_dec.second-off_run_RA_Dec.second;
        double delta_angle = pow(delta_ra*delta_ra+delta_dec*delta_dec,0.5);
        if (delta_angle<5.)
        {
            //std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " coordinate rejected." << std::endl;
            continue;
        }
        if (CheckBrightGammaRaySource(off_run_RA_Dec.first,off_run_RA_Dec.second))
        {
            //std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " coordinate rejected." << std::endl;
            continue;
        }
        std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " selected." << std::endl;
        OFF_runlist.push_back(OFF_runlist_input[off_run].second);
        OFF_pointing.push_back(run_pointing);
        OFF_pointing_radec.push_back(off_run_RA_Dec);
        double NSB_thisrun = GetRunPedestalVar(int(OFF_runlist_input[off_run].second));
        OFF_NSB.push_back(NSB_thisrun);
        double exposure_thisrun = GetRunUsableTime(OFF_filename,OFF_runlist_input[off_run].second);
        OFF_time.push_back(exposure_thisrun);
        OFF_L3Rate.push_back(L3_rate);
        double MJD_thisrun = double(GetRunMJD(OFF_filename,OFF_runlist_input[off_run].second));
        OFF_MJD.push_back(MJD_thisrun);

    }


    int ON_runnumber = 0;
    double ON_pointing_RA = 0.;
    double ON_pointing_Dec = 0.;
    double ON_exposure_hour = 0.;
    vector<pair<double,double>> ON_timecut;
    vector<int> OFF_runnumber;
    vector<double> OFF_exposure_hour;

    RunListTree->Branch("ON_runnumber",&ON_runnumber,"ON_runnumber/I");
    RunListTree->Branch("ON_exposure_hour",&ON_exposure_hour,"ON_exposure_hour/D");
    RunListTree->Branch("ON_pointing_RA",&ON_pointing_RA,"ON_pointing_RA/D");
    RunListTree->Branch("ON_pointing_Dec",&ON_pointing_Dec,"ON_pointing_Dec/D");
    RunListTree->Branch("ON_timecut","std::vector<std::pair<double,double>>",&ON_timecut);
    RunListTree->Branch("OFF_runnumber","std::vector<int>",&OFF_runnumber);
    RunListTree->Branch("OFF_exposure_hour","std::vector<double>",&OFF_exposure_hour);

    vector<int> exclusion_list = new_list;
    for (int off_run=0;off_run<OFF_runlist.size();off_run++)
    {
        double dist_crab = pow(OFF_pointing_radec.at(off_run).first-83.633,2)+pow(OFF_pointing_radec.at(off_run).second-22.014,2);
        double dist_mrk421 = pow(OFF_pointing_radec.at(off_run).first-166.079,2)+pow(OFF_pointing_radec.at(off_run).second-38.195,2);
        double dist_mrk501 = pow(OFF_pointing_radec.at(off_run).first-253.467,2)+pow(OFF_pointing_radec.at(off_run).second-39.76,2);
        double dist_geminga = pow(OFF_pointing_radec.at(off_run).first-98.117,2)+pow(OFF_pointing_radec.at(off_run).second-17.367,2);
        if (dist_crab<2.*2.)
        {
            exclusion_list.push_back(OFF_runlist.at(off_run));
        }
        if (dist_mrk421<2.*2.)
        {
            exclusion_list.push_back(OFF_runlist.at(off_run));
        }
        if (dist_mrk501<2.*2.)
        {
            exclusion_list.push_back(OFF_runlist.at(off_run));
        }
        if (dist_geminga<2.*2.)
        {
            exclusion_list.push_back(OFF_runlist.at(off_run));
        }
    }
    std::cout << "Find matched OFF run list..." << std::endl;
    std::cout << "ON runs = " << new_list.size() << std::endl;
    for (int on_run=0;on_run<new_list.size();on_run++)
    {

        ON_runnumber = new_list.at(on_run);
        ON_timecut = GetRunTimecuts(ON_runnumber);
        std::cout << "Selected ON run " << ON_runnumber << std::endl;
        std::cout << "time cut: ";
        for (int entry=0;entry<ON_timecut.size();entry++)
        {
            std::cout << ON_timecut.at(entry).first << "-" << ON_timecut.at(entry).second << ", ";
        }
        std::cout << std::endl;

        char char_runnumber[50];
        sprintf(char_runnumber, "%i", ON_runnumber);
        string ON_filename;
        ON_filename = TString(SMI_INPUT+"/"+string(char_runnumber)+".anasum.root");
        std::pair<double,double> on_run_elev_azim = GetRunElevAzim(ON_filename,ON_runnumber);
        std::pair<double,double> on_run_RA_Dec = GetRunRaDec(ON_filename,ON_runnumber);
        ON_pointing_RA = on_run_RA_Dec.first;
        ON_pointing_Dec = on_run_RA_Dec.second;
        double on_run_NSB = GetRunPedestalVar(ON_runnumber);
        ON_exposure_hour = GetRunUsableTime(ON_filename,ON_runnumber)/3600.;
        double ON_L3Rate = GetRunL3Rate(ON_runnumber);
        pair<double,double> ON_pointing = new_list_pointing[on_run];
        double on_run_MJD = double(GetRunMJD(ON_filename,ON_runnumber));

        double total_off_time = 0.;
        OFF_runnumber.clear();
        OFF_exposure_hour.clear();
        bool continue_to_find_match = true;
        double ratio_exposure = 1.0;
        if (UseDBOnly)
        {
            ratio_exposure = 2.0;
        }
        double sum_on_elev = 0.;
        double sum_off_elev = 0.;
        int matched_runnumber = FindAMatchedRunFromPairList(ON_runnumber, source_strip, exclusion_list, false, 0, sum_on_elev-sum_off_elev);
        if (matched_runnumber==0)
        {
            matched_runnumber = FindAMatchedRun(ON_runnumber, ON_pointing, on_run_NSB, ON_L3Rate, on_run_MJD, OFF_runlist, OFF_pointing, OFF_NSB, OFF_L3Rate, OFF_MJD, exclusion_list, false, sum_on_elev-sum_off_elev); 
        }
        if (matched_runnumber==0)
        {
            std::cout << "ON run " << ON_runnumber << " failed to find a match." << std::endl;
            std::cout << "ON run elev = " << ON_pointing.first << std::endl;
            std::cout << "ON run azim = " << ON_pointing.second << std::endl;
            std::cout << "ON run NSB = " << on_run_NSB << std::endl;
            std::cout << "ON run L3 rate = " << ON_L3Rate << std::endl;
            ON_runnumber = 0;
            OFF_runnumber.push_back(0);
            OFF_exposure_hour.push_back(0);
            continue_to_find_match = false;
        }
        else
        {
            std::cout << "ON run " << ON_runnumber << " found a match " << matched_runnumber << std::endl;
            exclusion_list.push_back(matched_runnumber);
            string OFF_filename;
            char char_off_runnumber[50];
            sprintf(char_off_runnumber, "%i", matched_runnumber);
            OFF_filename = TString(SMI_INPUT+"/"+string(char_off_runnumber)+".anasum.root");
            OFF_runnumber.push_back(matched_runnumber);
            double off_time = GetRunUsableTime(OFF_filename,matched_runnumber)/3600.;
            OFF_exposure_hour.push_back(off_time);
            total_off_time += off_time;
            std::pair<double,double> off_run_pointing = GetRunElevAzim(OFF_filename,matched_runnumber);
            sum_on_elev += ON_pointing.first*off_time;
            sum_off_elev += off_run_pointing.first*off_time;
        }


        if (ON_runnumber==0) continue;
        RunListTree->Fill();
    }

}

double GetDistance(double x1, double y1, double x2, double y2)
{
    double distance = pow((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2),0.5);
    return distance;
}

vector<std::pair<string,int>> GetDBRunList(string epoch, double elev_low, double elev_up) 
{

    vector<std::pair<string,int>> list;

    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_version = "";
    string acc_ra = "";
    string acc_dec = "";
    string acc_elev = "";
    string acc_az = "";
    string acc_fir = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;
    int runnumber = 0;
    double TelRAJ2000_tmp = 0.;
    double TelDecJ2000_tmp = 0.;
    double TelElevation_avg = 0.;
    double TelAzimuth_avg = 0.;
    double Tel0Fir_avg = 0.;

    ifstream myfile (SMI_AUX+"/diagnostics.txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            acc_version = "";
            acc_ra = "";
            acc_dec = "";
            acc_elev = "";
            acc_az = "";
            acc_fir = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if (nth_line<84) continue;
                if(line[i] == delimiter)
                {
                    if (nth_delimiter==0)
                    {
                        runnumber = std::stoi(acc_runnumber,nullptr,10);
                        //std::cout << "Processing run " << runnumber << std::endl;
                    }
                    if (nth_delimiter==7) 
                    {
                        TelElevation_avg = std::stod(acc_elev,&sz);
                        if (std::stoi(acc_version,nullptr,10)!=2) TelElevation_avg = 0.;
                    }
                    if (nth_delimiter==8) 
                    {
                        TelAzimuth_avg = std::stod(acc_az,&sz);
                        if (std::stoi(acc_version,nullptr,10)!=2) TelAzimuth_avg = 0.;
                    }
                    if (nth_delimiter==106) 
                    {
                        Tel0Fir_avg = std::stod(acc_fir,&sz);
                        if (std::stoi(acc_version,nullptr,10)!=2) Tel0Fir_avg = 0.;
                    }
                    if (nth_delimiter==45) 
                    {
                        TelRAJ2000_tmp = std::stod(acc_ra,&sz);
                        if (std::stoi(acc_version,nullptr,10)!=2) TelRAJ2000_tmp = 0.;
                    }
                    if (nth_delimiter==46) 
                    {
                        TelDecJ2000_tmp = std::stod(acc_dec,&sz);
                        if (std::stoi(acc_version,nullptr,10)!=2) TelDecJ2000_tmp = 0.;
                    }
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==0)
                {
                    acc_runnumber += line[i];
                }
                else if (nth_delimiter==1)
                {
                    acc_version += line[i];
                }
                else if (nth_delimiter==7)
                {
                    acc_elev += line[i];
                }
                else if (nth_delimiter==8)
                {
                    acc_az += line[i];
                }
                else if (nth_delimiter==106)
                {
                    acc_fir += line[i];
                }
                else if (nth_delimiter==45)
                {
                    acc_ra += line[i];
                }
                else if (nth_delimiter==46)
                {
                    acc_dec += line[i];
                }
            }
            if (nth_line>84)
            {
                pair<double,double> galactic_coords =  ConvertRaDecToGalactic(TelRAJ2000_tmp,TelDecJ2000_tmp);
                bool qualified = true;
                if (abs(galactic_coords.second)<5.) qualified = false;
                if (GetDistance(TelRAJ2000_tmp,TelDecJ2000_tmp,83.633,22.014)<3.) qualified = false; //Crab
                if (GetDistance(TelRAJ2000_tmp,TelDecJ2000_tmp,166.079,38.195)<3.) qualified = false; //Mrk 421
                if (GetDistance(TelRAJ2000_tmp,TelDecJ2000_tmp,98.117,17.367)<3.) qualified = false; //Geminga
                if (RunTypeCategory(runnumber,false)!=0) qualified = false;
                double exposure_minute = GetRunUsableTime("",runnumber)/60.;
                if (exposure_minute<10.) qualified = false;
                double L3_rate = GetRunL3Rate(runnumber);
                if (TelElevation_avg<elev_low-5.) qualified = false;
                if (TelElevation_avg>elev_up+5.) qualified = false;
                if (Tel0Fir_avg>0.3) qualified = false;
                if (qualified)
                {
                    if (epoch=="V6")
                    {
                        if (runnumber>=63373)
                        {
                            //std::cout << "Adding run " << runnumber << std::endl;
                            //std::cout << "T0 FIR temp = " << Tel0Fir_avg << std::endl;
                            list.push_back(std::make_pair("",runnumber));
                        }
                    }
                    else if (epoch=="V5")
                    {
                        if (runnumber>=46642 && runnumber<63373)
                        {
                            list.push_back(std::make_pair("",runnumber));
                        }
                    }
                }
            }
            nth_line += 1;
            runnumber = 0;
        }
        myfile.close();
    }
    else std::cout << "Unable to open file diagnostics.txt" << std::endl; 
    //std::cout << "diag file RA = " << TelRAJ2000_tmp << " Dec = " << TelDecJ2000_tmp << std::endl;

    return list;

}

vector<pair<string,int>> GetNewRunListFromTree(TTree * RunListTree)
{
    RunListTree->SetBranchStatus("*", 0);
    RunListTree->SetBranchStatus("ON_runnumber", 1);
    int ON_runnumber_ptr;
    RunListTree->SetBranchAddress("ON_runnumber",&ON_runnumber_ptr);

    vector<pair<string,int>> Data_runlist_new;
    for (int on_run=0;on_run<RunListTree->GetEntries();on_run++)
    {
        RunListTree->GetEntry(on_run);
        if (ON_runnumber_ptr==0) continue;
        std::cout << "ON run " << ON_runnumber_ptr << " is good for finding an imposter." << std::endl;
        Data_runlist_new.push_back(std::make_pair("ON run",ON_runnumber_ptr));
    }
    return Data_runlist_new;
}

void PrepareRunList(string target_data, double tel_elev_lower_input, double tel_elev_upper_input, double input_theta2_cut_lower, double input_theta2_cut_upper, bool isON, bool doImposter, int GammaModel)
{

    sprintf(target, "%s", target_data.c_str());
    ResetPublicVariables(TString(target));

    SMI_INPUT = string(std::getenv("SMI_INPUT"));
    SMI_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    SMI_DIR = string(std::getenv("SMI_DIR"));
    SMI_AUX = string(std::getenv("SMI_AUX"));

    TH1::SetDefaultSumw2();

    //camera_theta2_cut_lower = input_theta2_cut_lower;
    //camera_theta2_cut_upper = input_theta2_cut_upper;
    //sprintf(theta2_cut_tag, "_Theta2%dto%d", int(camera_theta2_cut_lower), int(camera_theta2_cut_upper));
    sprintf(theta2_cut_tag, "");
    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;
    sprintf(elev_cut_tag, "_TelElev%dto%d%s", int(TelElev_lower), int(TelElev_upper), Azim_region.c_str());
    MSCW_cut_blind = MSCW_cut_moderate;
    MSCL_cut_blind = MSCL_cut_moderate;

    pair<double,double> source_ra_dec = GetSourceRaDec(TString(target));
    mean_tele_point_ra = source_ra_dec.first;
    mean_tele_point_dec = source_ra_dec.second;
    pair<double,double> mean_tele_point_l_b = ConvertRaDecToGalactic(mean_tele_point_ra, mean_tele_point_dec);
    mean_tele_point_l = mean_tele_point_l_b.first;
    mean_tele_point_b = mean_tele_point_l_b.second;

    map_center_x = mean_tele_point_ra;
    map_center_y = mean_tele_point_dec;
    if (UseGalacticCoord)
    {
        map_center_x = mean_tele_point_l;
        map_center_y = mean_tele_point_b;
    }

    TString ONOFF_tag;
    if (isON) 
    {
        source_theta_cut = 0.;
        ONOFF_tag = "ON";
    }
    else
    {
        ONOFF_tag = "OFF";
    }
    char char_SignalStrength[50];
    sprintf(char_SignalStrength, "%i", GammaModel);
    ONOFF_tag += TString("_Model")+TString(char_SignalStrength);


    SizeSecondMax_Cut = 600.;
    if (TString(target).Contains("V4")) SizeSecondMax_Cut = 400.;
    if (TString(target).Contains("V5")) SizeSecondMax_Cut = 400.;

    std::cout << "Get a list of target observation runs" << std::endl;

    std::istringstream ss_input;
    ss_input.str(target);
    std::string source_strip("");
    std::string source_off_strip("");
    std::string segment;
    std::vector<std::string> seglist;
    while(std::getline(ss_input, segment, '_'))
    {
        seglist.push_back(segment);
    }
    for (int s=0;s<seglist.size()-1;s++)
    {
        if (s>0)
        {
            source_strip.append("_");
            source_off_strip.append("_");
        }
        source_strip.append(seglist.at(s));
        source_off_strip.append(seglist.at(s));
    }
    source_off_strip.append("_OFF");

    vector<pair<string,int>> Data_runlist_init = GetRunList(source_strip);
    std::cout << "Getting runlist elevation..." << std::endl;
    vector<double> Data_runlist_init_elev = GetRunElevationList(Data_runlist_init);
    std::cout << "Getting runlist RaDec..." << std::endl;
    vector<double> Data_runlist_init_radec = GetRunRaDecList(Data_runlist_init);
    std::cout << "Sorting list by elevation..." << std::endl;
    std::cout << "Data_runlist_init.size() = " << Data_runlist_init.size() << std::endl;
    std::cout << "Data_runlist_init_elev.size() = " << Data_runlist_init_elev.size() << std::endl;
    SortingList(&Data_runlist_init,&Data_runlist_init_elev);
    //SortingList(&Data_runlist_init,&Data_runlist_init_radec);

    std::cout << "Getting background list..." << std::endl;
    vector<pair<string,int>> Dark_runlist_init;
    if (!UseDBOnly)
    {
        if (TString(target).Contains("V6")) Dark_runlist_init = GetRunList("OffRunsV6");
        if (TString(target).Contains("V5")) Dark_runlist_init = GetRunList("OffRunsV5");
        //Dark_runlist_init = GetRunList(source_off_strip);
    }
    else
    {
        if (TString(target).Contains("V6")) Dark_runlist_init = GetDBRunList("V6", tel_elev_lower_input, tel_elev_upper_input);
        if (TString(target).Contains("V5")) Dark_runlist_init = GetDBRunList("V5", tel_elev_lower_input, tel_elev_upper_input);
    }
    
    std::cout << "initial Dark_runlist size = " << Dark_runlist_init.size() << std::endl;
    TTree RunListTree("RunListTree","ON data runn list tree");
    if (!doImposter)
    {
        SelectONRunList(&RunListTree, Data_runlist_init, Dark_runlist_init, tel_elev_lower_input, tel_elev_upper_input, source_ra_dec, source_strip);
    }
    else
    {
        //SelectONRunList(&RunListTree, Data_runlist_init, Dark_runlist_init, tel_elev_lower_input, tel_elev_upper_input, source_ra_dec);
        int iteration = 1;
        if (TString(target).Contains("Imposter1")) iteration = 1;
        if (TString(target).Contains("Imposter2")) iteration = 2;
        if (TString(target).Contains("Imposter3")) iteration = 3;
        if (TString(target).Contains("Imposter4")) iteration = 4;
        if (TString(target).Contains("Imposter5")) iteration = 5;
        if (TString(target).Contains("Imposter6")) iteration = 6;
        TString root_file_name = "Netflix_RunList_"+TString(target)+"_"+TString(output_file_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+TString(mjd_cut_tag)+"_"+ONOFF_tag+"_G0.root";
        //vector<pair<string,int>> Data_runlist_new = GetNewRunListFromTree(&RunListTree);
        //RunListTree.Reset();
        //std::cout << "Empty the tree before searching for imposters, RunListTree entries = " << RunListTree.GetEntries() << std::endl;
        SelectImposterRunList(&RunListTree, root_file_name, Data_runlist_init, Dark_runlist_init, tel_elev_lower_input, tel_elev_upper_input, iteration, source_ra_dec, source_strip);
    }

    RunListTree.SetBranchStatus("*", 0);
    RunListTree.SetBranchStatus("ON_runnumber", 1);
    RunListTree.SetBranchStatus("ON_pointing_RA", 1);
    RunListTree.SetBranchStatus("ON_pointing_Dec", 1);
    RunListTree.SetBranchStatus("ON_exposure_hour", 1);
    RunListTree.SetBranchStatus("OFF_runnumber", 1);
    RunListTree.SetBranchStatus("OFF_exposure_hour", 1);
    int ON_runnumber_ptr;
    double ON_exposure_hour_ptr;
    double ON_pointing_RA_ptr;
    double ON_pointing_Dec_ptr;
    vector<int>* OFF_runnumber_ptr = new std::vector<int>(10);
    vector<double>* OFF_exposure_hour_ptr = new std::vector<double>(10);
    RunListTree.SetBranchAddress("ON_runnumber",&ON_runnumber_ptr);
    RunListTree.SetBranchAddress("ON_exposure_hour",&ON_exposure_hour_ptr);
    RunListTree.SetBranchAddress("ON_pointing_RA",&ON_pointing_RA_ptr);
    RunListTree.SetBranchAddress("ON_pointing_Dec",&ON_pointing_Dec_ptr);
    RunListTree.SetBranchAddress("OFF_runnumber",&OFF_runnumber_ptr);
    RunListTree.SetBranchAddress("OFF_exposure_hour",&OFF_exposure_hour_ptr);
    int ON_runnumber;
    double ON_exposure_hour;
    double ON_pointing_RA;
    double ON_pointing_Dec;
    vector<int> OFF_runnumber;
    vector<double> OFF_exposure_hour;
    TTree RunListTree_subgroup("RunListTree_subgroup","ON data run list tree subgroup");
    RunListTree_subgroup.Branch("ON_runnumber",&ON_runnumber,"ON_runnumber/I");
    RunListTree_subgroup.Branch("ON_exposure_hour",&ON_exposure_hour,"ON_exposure_hour/D");
    RunListTree_subgroup.Branch("ON_pointing_RA",&ON_pointing_RA,"ON_pointing_RA/D");
    RunListTree_subgroup.Branch("ON_pointing_Dec",&ON_pointing_Dec,"ON_pointing_Dec/D");
    RunListTree_subgroup.Branch("OFF_runnumber","std::vector<int>",&OFF_runnumber);
    RunListTree_subgroup.Branch("OFF_exposure_hour","std::vector<double>",&OFF_exposure_hour);
    std::cout << __LINE__ << std::endl;
    std::cout << "RunListTree.GetEntries() = " << RunListTree.GetEntries() << std::endl;
    int usable_runs = 0;
    int saved_files = 0;
    double exposure_hour_limit = exposure_limit;
    for (int x_idx=0;x_idx<Skymap_nzones_x;x_idx++)
    {
        for (int y_idx=0;y_idx<Skymap_nzones_y;y_idx++)
        {
            int group_index = 0;
            double exposure_hour_sum = 0.;
            for (int on_run=0;on_run<RunListTree.GetEntries();on_run++)
            {
                sprintf(group_tag, "_G%d", group_index);
                sprintf(map_x_tag, "_X%d", x_idx);
                sprintf(map_y_tag, "_Y%d", y_idx);
                double map_x_bin_upper = double(x_idx+1)*2.0*Skymap_size_x/double(Skymap_nzones_x)-Skymap_size_x+map_center_x;
                double map_x_bin_lower = double(x_idx)*2.0*Skymap_size_x/double(Skymap_nzones_x)-Skymap_size_x+map_center_x;
                double map_y_bin_upper = double(y_idx+1)*2.0*Skymap_size_y/double(Skymap_nzones_y)-Skymap_size_y+map_center_y;
                double map_y_bin_lower = double(y_idx)*2.0*Skymap_size_y/double(Skymap_nzones_y)-Skymap_size_y+map_center_y;

                RunListTree.GetEntry(on_run);
                ON_runnumber = ON_runnumber_ptr;
                ON_exposure_hour = ON_exposure_hour_ptr;
                ON_pointing_RA = ON_pointing_RA_ptr;
                ON_pointing_Dec = ON_pointing_Dec_ptr;
                OFF_runnumber = *OFF_runnumber_ptr;
                OFF_exposure_hour = *OFF_exposure_hour_ptr;

                double run_pointing_x = ON_pointing_RA;
                double run_pointing_y = ON_pointing_Dec;
                pair<double,double> run_pointing_l_b = ConvertRaDecToGalactic(run_pointing_x,run_pointing_y);
                if (UseGalacticCoord)
                {
                    run_pointing_x = run_pointing_l_b.first;
                    run_pointing_y = run_pointing_l_b.second;
                }
                if (Skymap_nzones_x>1 || Skymap_nzones_y>1)
                {
                    if (run_pointing_x>map_x_bin_upper) continue;
                    if (run_pointing_x<map_x_bin_lower) continue;
                    if (run_pointing_y>map_y_bin_upper) continue;
                    if (run_pointing_y<map_y_bin_lower) continue;
                }

                RunListTree_subgroup.Fill();

                if (OFF_runnumber.at(0)!=0)
                {
                    usable_runs += 1;
                }

                exposure_hour_sum += ON_exposure_hour;
                
                if (exposure_hour_sum>exposure_hour_limit || on_run==RunListTree.GetEntries()-1)
                {
                    TFile OutputFile(TString(SMI_OUTPUT)+"/Netflix_RunList_"+TString(target)+"_"+TString(output_file_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+TString(mjd_cut_tag)+"_"+ONOFF_tag+group_tag+map_x_tag+map_y_tag+".root","recreate");
                    RunListTree_subgroup.Write();
                    OutputFile.Close();
                    RunListTree_subgroup.Reset();
                    exposure_hour_sum = 0.;
                    group_index += 1;
                    saved_files += 1;
                }

            }
        }
    }

    std::cout << "initial runs = " << Data_runlist_init.size() << std::endl;
    std::cout << "usable runs = " << usable_runs << std::endl;
    std::cout << "saved files = " << saved_files << std::endl;
    std::cout << "Done." << std::endl;

}
