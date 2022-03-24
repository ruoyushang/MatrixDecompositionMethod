
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

//#include "../../EventDisplay/VEvndispRunParameter.h"
#include "/home/rshang/MatrixDecompositionMethod/EventDisplay/VEvndispRunParameter.h"

#include "GetRunList.h"
#include "NetflixParameters.h"

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

    ifstream myfile (SMI_AUX+"/diagnostics_20220111.txt");
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
                    if (nth_delimiter==103 && std::stoi(acc_runnumber,nullptr,10)==run_number) 
                    {
                        NSB = std::stod(acc_nsb,&sz);
                        if (std::stoi(acc_version,nullptr,10)!=2) NSB = 0.;
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
                else if (nth_delimiter==103)
                {
                    acc_nsb += line[i];
                }
            }
            nth_line += 1;
        }
        myfile.close();
    }
    else std::cout << "Unable to open file diagnostics.txt" << std::endl; 

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
    ifstream myfile (SMI_AUX+"/category_allruns.txt");
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
                if (acc_runtype=="science" || acc_runtype=="NULL" || acc_runtype=="None")
                {
                    runtype = 0;
                }
                else if (acc_runtype=="reducedhv")
                {
                    runtype = 1;
                    //runtype = 0;
                    if (doPrint)
                    {
                        std::cout << "Run " << run_number << " rejected. Run Type " << acc_runtype << std::endl;
                    }
                }
                else
                {
                    runtype = 2;
                    if (doPrint)
                    {
                        std::cout << "Run " << run_number << " rejected. Run Type " << acc_runtype << std::endl;
                    }
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

    if (usable_time==0.)
    {
        char run_number_char[50];
        sprintf(run_number_char, "%i", run_number);
        TFile*  input_file = TFile::Open(file_name.c_str());
        TTree* pointing_tree = nullptr;
        pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number_char)+"/stereo/pointingDataReduced"));
        pointing_tree->SetBranchStatus("*",0);
        pointing_tree->SetBranchStatus("Time",1);
        pointing_tree->SetBranchAddress("Time",&Time);
        double total_entries = (double)pointing_tree->GetEntries();
        pointing_tree->GetEntry(0);
        int time_start = Time;
        pointing_tree->GetEntry(total_entries-1);
        int time_end = Time;
        input_file->Close();
        usable_time = time_end-time_start;
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

    int run_mjd = 0.;

    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());
    TTree* pointing_tree = nullptr;
    pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/pointingDataReduced"));
    pointing_tree->SetBranchStatus("*",0);
    pointing_tree->SetBranchStatus("MJD",1);
    pointing_tree->SetBranchAddress("MJD",&MJD_UInt_t);
    double total_entries = (double)pointing_tree->GetEntries();
    pointing_tree->GetEntry(0);
    run_mjd = MJD_UInt_t;
    input_file->Close();
    //std::cout << "root file MJD = " << run_mjd << std::endl;

    //string line;
    //char delimiter = ' ';
    //string acc_runnumber = "";
    //string acc_version = "";
    //string acc_nsb = "";
    //int nth_line = 0;
    //int nth_delimiter = 0;
    //std::string::size_type sz;

    //ifstream myfile (SMI_AUX+"/diagnostics_20210705.txt");
    //if (myfile.is_open())
    //{
    //    while ( getline(myfile,line) )
    //    {
    //        acc_runnumber = "";
    //        acc_version = "";
    //        acc_nsb = "";
    //        nth_delimiter = 0;
    //        for(int i = 0; i < line.size(); i++)
    //        {
    //            if (nth_line<84) continue;
    //            if(line[i] == delimiter)
    //            {
    //                if (nth_delimiter==4 && std::stoi(acc_runnumber,nullptr,10)==run) 
    //                {
    //                    run_mjd = std::stoi(acc_nsb,&sz);
    //                    if (std::stoi(acc_version,nullptr,10)!=2) run_mjd = 0;
    //                }
    //                nth_delimiter += 1;
    //            }
    //            else if (nth_delimiter==0)
    //            {
    //                acc_runnumber += line[i];
    //            }
    //            else if (nth_delimiter==1)
    //            {
    //                acc_version += line[i];
    //            }
    //            else if (nth_delimiter==4)
    //            {
    //                acc_nsb += line[i];
    //            }
    //        }
    //        nth_line += 1;
    //    }
    //    myfile.close();
    //}
    //else std::cout << "Unable to open file diagnostics.txt" << std::endl; 
    //std::cout << "diag file MJD = " << run_mjd << std::endl;

    return run_mjd;
}
bool PointingSelection(string file_name,int run, double Elev_cut_lower, double Elev_cut_upper)
{
    if (run>100000) return true;
    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());

    TTree* pointing_tree = nullptr;
    pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/pointingDataReduced"));
    pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
    pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
    pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
    pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
    double total_entries = (double)pointing_tree->GetEntries();
    pointing_tree->GetEntry(int(total_entries/2.));
    int MJD = GetRunMJD(file_name,run);
    if (TelElevation<Elev_cut_lower) 
    {
        input_file->Close();
        return false;
    }
    if (TelElevation>Elev_cut_upper)
    {
        input_file->Close();
        return false;
    }
    if (Azim_region=="North")
    {
        if (!(TelAzimuth<45. || TelAzimuth>315.))
        {
            input_file->Close();
            return false;
        }
    }
    if (Azim_region=="South")
    {
        if (!(TelAzimuth>180.-45. && TelAzimuth<180.+45.))
        {
            input_file->Close();
            return false;
        }
    }
    if (Azim_region=="East")
    {
        if (!(TelAzimuth>90.-45. && TelAzimuth<90.+45.))
        {
            input_file->Close();
            return false;
        }
    }
    if (Azim_region=="West")
    {
        if (!(TelAzimuth>270.-45. && TelAzimuth<270.+45.))
        {
            input_file->Close();
            return false;
        }
    }
    if (TString(target).Contains("SS433"))
    {
        double TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
        double TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
        double delta_RA = TelRAJ2000_tmp - 287.9565;
        double delta_Dec = TelDecJ2000_tmp - 4.9827;
        double distance = pow(delta_RA*delta_RA + delta_Dec*delta_Dec,0.5);
        //if (distance>2.0)
        //{
        //    input_file->Close();
        //    return false;
        //}
    }
    else if (TString(target).Contains("MGRO_J2031"))
    {
        double TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
        double TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
        double delta_RA = TelRAJ2000_tmp - 308.041666667;
        double delta_Dec = TelDecJ2000_tmp - 41.4594444444;
        double distance = pow(delta_RA*delta_RA + delta_Dec*delta_Dec,0.5);
        if (distance>2.0)
        {
            input_file->Close();
            return false;
        }
    }
    else if (TString(target).Contains("MGRO_J1908"))
    {
        double TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
        double TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
        double delta_RA = TelRAJ2000_tmp - 287.188333333;
        double delta_Dec = TelDecJ2000_tmp - 6.16233333333;
        double distance = pow(delta_RA*delta_RA + delta_Dec*delta_Dec,0.5);
        if (distance>2.0)
        {
            input_file->Close();
            return false;
        }
    }
    else if (TString(target).Contains("WComae"))
    {
        double TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
        double TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
        //double delta_RA = TelRAJ2000_tmp - 185.382;
        //double delta_Dec = TelDecJ2000_tmp - 28.233;
        double delta_RA = TelRAJ2000_tmp - (185.360+184.616)/2.;
        double delta_Dec = TelDecJ2000_tmp - (30.191+30.130)/2.;
        double distance = pow(delta_RA*delta_RA + delta_Dec*delta_Dec,0.5);
        if (distance>2.0)
        {
            input_file->Close();
            return false;
        }
        //if (MJD>57250)
        //{
        //    input_file->Close();
        //    return false;
        //}
    }
    int MJD_start = 59246; // 2021-02-01
    int MJD_end = 59327; // 2021-04-23
    double X = ( double(MJD_start) - double(MJD) )/162.5;
    double Phase=X-floor(X);
    if (TString(target).Contains("SS433Half1"))
    {
        if (Phase>0.5)
        {
            input_file->Close();
            return false;
        }
    }
    if (TString(target).Contains("SS433Half2"))
    {
        if (Phase<=0.5)
        {
            input_file->Close();
            return false;
        }
    }
    if (TString(target).Contains("MGRO_J1908Half1"))
    {
        if (Phase>0.5)
        {
            input_file->Close();
            return false;
        }
    }
    if (TString(target).Contains("MGRO_J1908Half2"))
    {
        if (Phase<=0.5)
        {
            input_file->Close();
            return false;
        }
    }
    if (TString(target).Contains("SS433Quad1"))
    {
        if (!(Phase<0.25))
        {
            input_file->Close();
            return false;
        }
    }
    if (TString(target).Contains("SS433Quad2"))
    {
        if (!(Phase>=0.25 && Phase<0.5))
        {
            input_file->Close();
            return false;
        }
    }
    if (TString(target).Contains("SS433Quad3"))
    {
        if (!(Phase>=0.5 && Phase<0.75))
        {
            input_file->Close();
            return false;
        }
    }
    if (TString(target).Contains("SS433Quad4"))
    {
        if (!(Phase>=0.75))
        {
            input_file->Close();
            return false;
        }
    }

    input_file->Close();
    return true;
}

pair<double,double> GetRunRaDec(string file_name, int run)
{

    double TelRAJ2000_tmp = 0.;
    double TelDecJ2000_tmp = 0.;

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
    //std::cout << "root file RA = " << TelRAJ2000_tmp << " Dec = " << TelDecJ2000_tmp << std::endl;

    //string line;
    //char delimiter = ' ';
    //string acc_runnumber = "";
    //string acc_version = "";
    //string acc_elev = "";
    //string acc_az = "";
    //int nth_line = 0;
    //int nth_delimiter = 0;
    //std::string::size_type sz;

    //ifstream myfile (SMI_AUX+"/diagnostics_20210705.txt");
    //if (myfile.is_open())
    //{
    //    while ( getline(myfile,line) )
    //    {
    //        acc_runnumber = "";
    //        acc_version = "";
    //        acc_elev = "";
    //        acc_az = "";
    //        nth_delimiter = 0;
    //        for(int i = 0; i < line.size(); i++)
    //        {
    //            if (nth_line<84) continue;
    //            if(line[i] == delimiter)
    //            {
    //                if (nth_delimiter==45 && std::stoi(acc_runnumber,nullptr,10)==run) 
    //                {
    //                    TelRAJ2000_tmp = std::stod(acc_elev,&sz);
    //                    if (std::stoi(acc_version,nullptr,10)!=2) TelRAJ2000_tmp = 0.;
    //                }
    //                if (nth_delimiter==46 && std::stoi(acc_runnumber,nullptr,10)==run) 
    //                {
    //                    TelDecJ2000_tmp = std::stod(acc_az,&sz);
    //                    if (std::stoi(acc_version,nullptr,10)!=2) TelDecJ2000_tmp = 0.;
    //                }
    //                nth_delimiter += 1;
    //            }
    //            else if (nth_delimiter==0)
    //            {
    //                acc_runnumber += line[i];
    //            }
    //            else if (nth_delimiter==1)
    //            {
    //                acc_version += line[i];
    //            }
    //            else if (nth_delimiter==45)
    //            {
    //                acc_elev += line[i];
    //            }
    //            else if (nth_delimiter==46)
    //            {
    //                acc_az += line[i];
    //            }
    //        }
    //        nth_line += 1;
    //    }
    //    myfile.close();
    //}
    //else std::cout << "Unable to open file diagnostics.txt" << std::endl; 
    //std::cout << "diag file RA = " << TelRAJ2000_tmp << " Dec = " << TelDecJ2000_tmp << std::endl;

    return std::make_pair(TelRAJ2000_tmp,TelDecJ2000_tmp);
}

pair<double,double> GetRunElevAzim(string file_name, int run)
{
    double TelElevation_avg = 0.;
    double TelAzimuth_avg = 0.;

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
    for (int entry=0;entry<total_entries;entry++)
    {
        pointing_tree->GetEntry(entry);
        TelElevation_avg += TelElevation;
        TelAzimuth_avg += TelAzimuth;
    }
    TelElevation_avg = TelElevation_avg/total_entries;
    TelAzimuth_avg = TelAzimuth_avg/total_entries;
    input_file->Close();

    //std::cout << "root file elev = " << TelElevation_avg << " azim = " << TelAzimuth_avg << std::endl;

    //string line;
    //char delimiter = ' ';
    //string acc_runnumber = "";
    //string acc_version = "";
    //string acc_elev = "";
    //string acc_az = "";
    //int nth_line = 0;
    //int nth_delimiter = 0;
    //std::string::size_type sz;

    //ifstream myfile (SMI_AUX+"/diagnostics_20210705.txt");
    //if (myfile.is_open())
    //{
    //    while ( getline(myfile,line) )
    //    {
    //        acc_runnumber = "";
    //        acc_version = "";
    //        acc_elev = "";
    //        acc_az = "";
    //        nth_delimiter = 0;
    //        for(int i = 0; i < line.size(); i++)
    //        {
    //            if (nth_line<84) continue;
    //            if(line[i] == delimiter)
    //            {
    //                if (nth_delimiter==7 && std::stoi(acc_runnumber,nullptr,10)==run) 
    //                {
    //                    TelElevation_avg = std::stod(acc_elev,&sz);
    //                    if (std::stoi(acc_version,nullptr,10)!=2) TelElevation_avg = 0.;
    //                }
    //                if (nth_delimiter==8 && std::stoi(acc_runnumber,nullptr,10)==run) 
    //                {
    //                    TelAzimuth_avg = std::stod(acc_az,&sz);
    //                    if (std::stoi(acc_version,nullptr,10)!=2) TelAzimuth_avg = 0.;
    //                }
    //                nth_delimiter += 1;
    //            }
    //            else if (nth_delimiter==0)
    //            {
    //                acc_runnumber += line[i];
    //            }
    //            else if (nth_delimiter==1)
    //            {
    //                acc_version += line[i];
    //            }
    //            else if (nth_delimiter==7)
    //            {
    //                acc_elev += line[i];
    //            }
    //            else if (nth_delimiter==8)
    //            {
    //                acc_az += line[i];
    //            }
    //        }
    //        nth_line += 1;
    //    }
    //    myfile.close();
    //}
    //else std::cout << "Unable to open file diagnostics.txt" << std::endl; 
    //std::cout << "diag file elev = " << TelElevation_avg << " azim = " << TelAzimuth_avg << std::endl;

    return std::make_pair(TelElevation_avg,TelAzimuth_avg);
}
bool MJDSelection(string file_name,int run, int MJD_start_cut, int MJD_end_cut)
{
    if (MJD_start_cut==0 && MJD_end_cut==0) return true;
    if (run>100000) return true;
    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());

    TTree* pointing_tree = nullptr;
    pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/pointingDataReduced"));
    pointing_tree->SetBranchAddress("MJD",&MJD_UInt_t);
    double total_entries = (double)pointing_tree->GetEntries();
    pointing_tree->GetEntry(0);
    if (MJD_UInt_t<MJD_start_cut || MJD_UInt_t>MJD_end_cut) 
    {
        input_file->Close();
        return false;
    }
    input_file->Close();
    return true;
}
void SortingList(vector<pair<string,int>>* list, vector<double>* list_pointing)
{
    pair<string,int> temp_list;
    double temp_list_pointing;
    for (int run1=0;run1<list->size()-1;run1++)
    {
        for (int run2=run1+1;run2<list->size();run2++)
        {
            if (list_pointing->at(run1)<list_pointing->at(run2))
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

int FindAMatchedRun(int ON_runnumber, pair<double,double> ON_pointing, double ON_NSB, double ON_L3Rate, double ON_MJD, vector<int> OFF_runnumber, vector<pair<double,double>> OFF_pointing, vector<double> OFF_NSB, vector<double> OFF_L3Rate, vector<double> OFF_MJD, vector<int> exclusion_list) 
{
    int matched_runnumber = 0;
    double match_chi2 = 1e10;
    double threshold_dElev = 2.0;
    double threshold_dAirmass = 0.1;
    double threshold_dNSB = 0.2;
    double threshold_dL3Rate = 0.3;
    double threshold_dMJD = 3.*365.;
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
        double delta_elev = abs(ON_pointing.first-OFF_pointing[off_run].first);
        double delta_airmass = abs(1./sin(ON_pointing.first*M_PI/180.)-1./sin(OFF_pointing[off_run].first*M_PI/180.));

        double delta_azim = abs(ON_pointing.second-OFF_pointing[off_run].second);
        if (delta_azim>180.) delta_azim = 360.-delta_azim;

        if (delta_airmass>threshold_dAirmass) continue;
        if (delta_nsb>threshold_dNSB) continue;
        //if (delta_l3rate/ON_L3Rate>threshold_dL3Rate) continue;
        //if (delta_elev>threshold_dElev) continue;
        //if (delta_mjd>threshold_dMJD) continue;
        //
        //double chi2 = pow(delta_nsb,2); // matrix method doesn't make good prediction when dNSB is large.
        //double chi2 = pow(delta_airmass,2);
        //double chi2 = pow(delta_mjd,2);
        //double chi2 = pow(delta_l3rate,2);
        double chi2 = pow(delta_azim,2);

        if (chi2<match_chi2)
        {
            match_chi2 = chi2;
            matched_runnumber = OFF_runnumber[off_run];
        }
    }
    return matched_runnumber;
}

void SelectImposterRunList(TTree * RunListTree, vector<pair<string,int>> Data_runlist, vector<pair<string,int>> OFF_runlist_input, double Elev_cut_lower, double Elev_cut_upper, int MJD_start_cut, int MJD_end_cut, int iteration)
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
        if (!MJDSelection(filename,int(Data_runlist[run].second),MJD_start_cut,MJD_end_cut)) 
        {
            std::cout << int(Data_runlist[run].second) << " MJD rejected." << std::endl;
            continue;
        }
        double L3_rate = GetRunL3Rate(Data_runlist[run].second);
        if (L3_rate<150.)
        {
            std::cout << int(Data_runlist[run].second) << " L3 rate rejected." << std::endl;
            continue;
        }

        new_list.push_back(Data_runlist[run].second);
        new_list_pointing.push_back(GetRunElevAzim(filename,int(Data_runlist[run].second)));
        ON_NSB.push_back(NSB_thisrun);
        double exposure_thisrun = GetRunUsableTime(filename,Data_runlist[run].second);
        ON_time.push_back(exposure_thisrun);
        ON_L3Rate.push_back(GetRunL3Rate(Data_runlist[run].second));

    }

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
        double run_elevation = GetRunElevAzim(OFF_filename,int(OFF_runlist_input[off_run].second)).first;
        if (run_elevation<Elev_cut_lower-5.) continue;
        if (run_elevation>Elev_cut_upper+5.) continue;
        if (RunTypeCategory(OFF_runlist_input[off_run].second,true)!=0) 
        {
            std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " category rejected." << std::endl;
            continue;
        }
        double L3_rate = GetRunL3Rate(OFF_runlist_input[off_run].second);
        if (L3_rate<150.)
        {
            std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " L3 rate rejected." << std::endl;
            continue;
        }
        OFF_runlist.push_back(OFF_runlist_input[off_run].second);
        OFF_pointing.push_back(GetRunElevAzim(OFF_filename,int(OFF_runlist_input[off_run].second)));
        OFF_pointing_radec.push_back(GetRunRaDec(OFF_filename,int(OFF_runlist_input[off_run].second)));
        double NSB_thisrun = GetRunPedestalVar(int(OFF_runlist_input[off_run].second));
        OFF_NSB.push_back(NSB_thisrun);
        double exposure_thisrun = GetRunUsableTime(OFF_filename,OFF_runlist_input[off_run].second);
        OFF_time.push_back(exposure_thisrun);
        OFF_L3Rate.push_back(GetRunL3Rate(OFF_runlist_input[off_run].second));
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
    vector<int> imposter_list;
    vector<pair<double,double>> imposter_list_radec;
    for (int trial=0;trial<iteration;trial++)
    {
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

            int matched_runnumber = FindAMatchedRun(ON_runnumber, ON_pointing, on_run_NSB, ON_L3Rate, on_run_MJD, OFF_runlist, OFF_pointing, OFF_NSB, OFF_L3Rate, OFF_MJD, exclusion_list); 
            if (matched_runnumber==0)
            {
                std::cout << "ON run " << ON_runnumber << " failed to find an imposter." << std::endl;
            }
            else
            {
                std::cout << "ON run " << ON_runnumber << " found an imposter " << matched_runnumber << std::endl;
                if (trial==iteration-1)
                {
                    imposter_list.push_back(matched_runnumber);
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


    for (int on_run=0;on_run<imposter_list.size();on_run++)
    {

        ON_runnumber = imposter_list.at(on_run);
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
        while (total_off_time<1.0*ON_exposure_hour && continue_to_find_match)
        {

            int matched_runnumber = FindAMatchedRun(ON_runnumber, ON_pointing, on_run_NSB, ON_L3Rate, on_run_MJD, OFF_runlist, OFF_pointing, OFF_NSB, OFF_L3Rate, OFF_MJD, exclusion_list); 
            if (matched_runnumber==0)
            {
                std::cout << "Imposter run " << ON_runnumber << " failed to find a match." << std::endl;
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
            }

        }

        RunListTree->Fill();
    }

}

void SelectONRunList(TTree * RunListTree, vector<pair<string,int>> Data_runlist, vector<pair<string,int>> OFF_runlist_input, double Elev_cut_lower, double Elev_cut_upper, int MJD_start_cut, int MJD_end_cut)
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
        if (!MJDSelection(filename,int(Data_runlist[run].second),MJD_start_cut,MJD_end_cut)) 
        {
            std::cout << int(Data_runlist[run].second) << " MJD rejected." << std::endl;
            continue;
        }
        double L3_rate = GetRunL3Rate(Data_runlist[run].second);
        if (L3_rate<150.)
        {
            std::cout << int(Data_runlist[run].second) << " L3 rate rejected." << std::endl;
            continue;
        }

        new_list.push_back(Data_runlist[run].second);
        new_list_pointing.push_back(GetRunElevAzim(filename,int(Data_runlist[run].second)));
        ON_NSB.push_back(NSB_thisrun);
        double exposure_thisrun = GetRunUsableTime(filename,Data_runlist[run].second);
        ON_time.push_back(exposure_thisrun);
        ON_L3Rate.push_back(GetRunL3Rate(Data_runlist[run].second));

    }

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
        double run_elevation = GetRunElevAzim(OFF_filename,int(OFF_runlist_input[off_run].second)).first;
        if (run_elevation<Elev_cut_lower-5.) continue;
        if (run_elevation>Elev_cut_upper+5.) continue;
        if (RunTypeCategory(OFF_runlist_input[off_run].second,true)!=0) 
        {
            std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " category rejected." << std::endl;
            continue;
        }
        double L3_rate = GetRunL3Rate(OFF_runlist_input[off_run].second);
        if (L3_rate<150.)
        {
            std::cout << "OFF run " << int(OFF_runlist_input[off_run].second) << " L3 rate rejected." << std::endl;
            continue;
        }
        OFF_runlist.push_back(OFF_runlist_input[off_run].second);
        OFF_pointing.push_back(GetRunElevAzim(OFF_filename,int(OFF_runlist_input[off_run].second)));
        OFF_pointing_radec.push_back(GetRunRaDec(OFF_filename,int(OFF_runlist_input[off_run].second)));
        double NSB_thisrun = GetRunPedestalVar(int(OFF_runlist_input[off_run].second));
        OFF_NSB.push_back(NSB_thisrun);
        double exposure_thisrun = GetRunUsableTime(OFF_filename,OFF_runlist_input[off_run].second);
        OFF_time.push_back(exposure_thisrun);
        OFF_L3Rate.push_back(GetRunL3Rate(OFF_runlist_input[off_run].second));
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
        while (total_off_time<2.0*ON_exposure_hour && continue_to_find_match)
        {

            int matched_runnumber = FindAMatchedRun(ON_runnumber, ON_pointing, on_run_NSB, ON_L3Rate, on_run_MJD, OFF_runlist, OFF_pointing, OFF_NSB, OFF_L3Rate, OFF_MJD, exclusion_list); 
            if (matched_runnumber==0)
            {
                std::cout << "ON run " << ON_runnumber << " failed to find a match." << std::endl;
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
            }

        }

        RunListTree->Fill();
    }

}

void PrepareRunList(string target_data, double tel_elev_lower_input, double tel_elev_upper_input, int MJD_start_cut, int MJD_end_cut, double input_theta2_cut_lower, double input_theta2_cut_upper, bool isON, bool doImposter, int GammaModel)
{

    SMI_INPUT = string(std::getenv("SMI_INPUT"));
    SMI_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    SMI_DIR = string(std::getenv("SMI_DIR"));
    SMI_AUX = string(std::getenv("SMI_AUX"));

    TH1::SetDefaultSumw2();

    if (MJD_start_cut!=0 || MJD_end_cut!=0)
    {
        sprintf(mjd_cut_tag, "_MJD%dto%d", MJD_start_cut, MJD_end_cut);
    }
    camera_theta2_cut_lower = input_theta2_cut_lower;
    camera_theta2_cut_upper = input_theta2_cut_upper;
    sprintf(theta2_cut_tag, "_Theta2%dto%d", int(camera_theta2_cut_lower), int(camera_theta2_cut_upper));
    sprintf(target, "%s", target_data.c_str());
    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;
    sprintf(elev_cut_tag, "_TelElev%dto%d%s", int(TelElev_lower), int(TelElev_upper), Azim_region.c_str());
    MSCW_cut_blind = MSCW_cut_moderate;
    MSCL_cut_blind = MSCL_cut_moderate;

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
    char char_SignalStrength[50];
    sprintf(char_SignalStrength, "%i", GammaModel);
    ONOFF_tag += TString("_Model")+TString(char_SignalStrength);


    SizeSecondMax_Cut = 600.;
    if (TString(target).Contains("V4")) SizeSecondMax_Cut = 400.;
    if (TString(target).Contains("V5")) SizeSecondMax_Cut = 400.;

    std::cout << "Get a list of target observation runs" << std::endl;
    vector<pair<string,int>> Data_runlist_init = GetRunList(target);
    vector<pair<string,int>> Dark_runlist_init = GetRunList("OffRunsV6");
    if (TString(target).Contains("V5")) Dark_runlist_init = GetRunList("OffRunsV5");
    if (TString(target).Contains("V4")) Dark_runlist_init = GetRunList("OffRunsV4");
    std::cout << "initial Dark_runlist size = " << Dark_runlist_init.size() << std::endl;
    TTree RunListTree("RunListTree","ON data runn list tree");
    if (!doImposter)
    {
        SelectONRunList(&RunListTree, Data_runlist_init, Dark_runlist_init, tel_elev_lower_input, tel_elev_upper_input, MJD_start_cut, MJD_end_cut);
    }
    else
    {
        int iteration = 1;
        if (TString(target).Contains("Imposter1")) iteration = 1;
        if (TString(target).Contains("Imposter2")) iteration = 2;
        if (TString(target).Contains("Imposter3")) iteration = 3;
        if (TString(target).Contains("Imposter4")) iteration = 4;
        if (TString(target).Contains("Imposter5")) iteration = 5;
        SelectImposterRunList(&RunListTree, Data_runlist_init, Dark_runlist_init, tel_elev_lower_input, tel_elev_upper_input, MJD_start_cut, MJD_end_cut, iteration);
    }


    TFile OutputFile(TString(SMI_OUTPUT)+"/Netflix_RunList_"+TString(target)+"_"+TString(output_file_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+TString(mjd_cut_tag)+"_"+ONOFF_tag+".root","recreate");

    RunListTree.Write();
    OutputFile.Close();

    std::cout << "initial runs = " << Data_runlist_init.size() << std::endl;
    std::cout << "selected runs = " << RunListTree.GetEntries() << std::endl;
    std::cout << "Done." << std::endl;

}
