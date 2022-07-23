
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
double ra_sky_imposter = 0;
double dec_sky_imposter = 0;
double ra_sky_dark = 0;
double dec_sky_dark = 0;
double exposure_hours = 0.;
double shower_count = 0.;
double mean_elev = 0.;
double mean_azim = 0.;
double mean_nsb = 0.;
double exposure_hours_usable = 0.;
double exposure_hours_ref = 0.;
int MJD_Start = 2147483647;
int MJD_End = 0;
double map_x_bin_upper = 0.;
double map_x_bin_lower = 0.;
double map_y_bin_upper = 0.;
double map_y_bin_lower = 0.;

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

string SMI_INPUT;
string SMI_OUTPUT;
string SMI_DIR;
string SMI_AUX;

double GetCrabFlux(double energy_gev)
{
    double flux = 3.75*pow(10,-7)*pow(energy_gev/1000.,-2.467-0.16*log(energy_gev/1000.));
    return flux;
}
double GetModelFlux(double energy_gev, double amp, double index)
{
    double flux = amp*pow(10,-7)*pow(energy_gev/1000.,-1.*index);
    return flux;
}
double FindBrightStarWeight(double ra, double dec)
{
    double weight = 1.;
    for (int star=0;star<FaintStars_Data.size();star++)
    {
        double star_ra = FaintStars_Data.at(star).at(0);
        double star_dec = FaintStars_Data.at(star).at(1);
        double star_brightness = FaintStars_Data.at(star).at(3);
        double radius = pow((ra-star_ra)*(ra-star_ra)+(dec-star_dec)*(dec-star_dec),0.5);
        if (star_brightness>faint_brightness_cut) continue;
        if (radius>bright_star_radius_cut) continue;
        if (star_brightness<3.0)
        {
            weight = 1.08;
        }
        else if (star_brightness<4.0)
        {
            weight = 1.10;
        }
        else
        {
            weight = 1.04;
        }
    }
    return weight;
}
bool CoincideWithBrightStars(double ra, double dec)
{
    bool isCoincident = false;
    for (int star=0;star<BrightStars_Data.size();star++)
    {
        double star_ra = BrightStars_Data.at(star).at(0);
        double star_dec = BrightStars_Data.at(star).at(1);
        double star_brightness = BrightStars_Data.at(star).at(3);
        double radius = pow((ra-star_ra)*(ra-star_ra)+(dec-star_dec)*(dec-star_dec),0.5);
        if (star_brightness>brightness_cut) continue;
        if (radius>bright_star_radius_cut) continue;
        isCoincident = true;
    }
    return isCoincident;
}
bool CoincideWithGammaSources(double ra, double dec)
{
    bool isCoincident = false;
    for (int star=0;star<GammaSource_Data.size();star++)
    {
        double star_ra = GammaSource_Data.at(star).at(0);
        double star_dec = GammaSource_Data.at(star).at(1);
        double radius = pow((ra-star_ra)*(ra-star_ra)+(dec-star_dec)*(dec-star_dec),0.5);
        if (radius>0.2) continue;
        isCoincident = true;
    }
    return isCoincident;
}
pair<double,double> ConvertRaDecToGalactic(double Ra, double Dec)
{
    double delta = Dec*M_PI/180.;
    double delta_G = 27.12825*M_PI/180.;
    double alpha = Ra*M_PI/180.;
    double alpha_G = 192.85948*M_PI/180.;
    double l_NCP = 122.93192*M_PI/180.;
    double sin_b = sin(delta)*sin(delta_G)+cos(delta)*cos(delta_G)*cos(alpha-alpha_G);
    double cos_b = cos(asin(sin_b));
    double sin_l_NCP_m_l = cos(delta)*sin(alpha-alpha_G)/cos_b;
    double cos_l_NCP_m_l = (cos(delta_G)*sin(delta)-sin(delta_G)*cos(delta)*cos(alpha-alpha_G))/cos_b;
    double b = (asin(sin_b))*180./M_PI;
    double l = (l_NCP-atan2(sin_l_NCP_m_l,cos_l_NCP_m_l))*180./M_PI;
    double b_round = floor(b*pow(10,2))/pow(10,2);
    double l_round = floor(l*pow(10,2))/pow(10,2);
    return std::make_pair(l_round,b_round);
}
pair<double,double> GetSourceRaDec(TString source_name)
{
    double Source_RA = 0.;
    double Source_Dec = 0.;
    if (source_name.Contains("SgrA"))
    {
            Source_RA = 266.415;
                Source_Dec = -29.006;
    }
    if (source_name.Contains("Geminga"))
    {
            Source_RA = 98.476;
                Source_Dec = 17.770;
    }
    if (source_name.Contains("PKS1441"))
    {
            Source_RA = 220.987;
                Source_Dec = 25.029;
    }
    if (source_name.Contains("MS1221"))
    {
            Source_RA = 186.101;
                Source_Dec = 24.607;
    }
    if (source_name.Contains("Tycho"))
    {
            Source_RA = 6.340;
                Source_Dec = 64.130;
    }
    if (source_name.Contains("CTA1"))
    {
            Source_RA = 1.608;
                Source_Dec = 72.984;
    }
    if (source_name.Contains("S3_1227"))
    {
            Source_RA = 187.559;
                Source_Dec = 25.302;
    }
    if (source_name.Contains("SNR_G150p3Plus04p5"))
    {
            Source_RA = 67.82;
                Source_Dec = 55.89;
    }
    if (source_name.Contains("Crab"))
    {
            Source_RA = 83.633;
                Source_Dec = 22.014;
    }
    if (source_name.Contains("Mrk421"))
    {
            Source_RA = 166.079;
                Source_Dec = 38.195;
    }
    if (source_name.Contains("H1426"))
    {
            Source_RA = 217.136;
                Source_Dec = 42.673;
    }
    if (source_name.Contains("1ES0229"))
    {
            Source_RA = 38.222;
                Source_Dec = 20.273;
    }
    if (source_name.Contains("PKS1424"))
    {
            Source_RA = 216.750;
                Source_Dec = 23.783;
    }
    if (source_name.Contains("3C264"))
    {
            Source_RA = 176.271;
                Source_Dec = 19.606;
    }
    if (source_name.Contains("OJ287"))
    {
            Source_RA = 133.705;
                Source_Dec = 20.100;
    }
    if (source_name.Contains("RBS0413"))
    {
            Source_RA = 49.946;
                Source_Dec = 18.762;
    }
    if (source_name.Contains("PG1553"))
    {
            Source_RA = 238.936;
                Source_Dec = 11.195;
    }
    if (source_name.Contains("Segue1"))
    {
            Source_RA = 151.767;
                Source_Dec = 16.082;
    }
    if (source_name.Contains("Coma"))
    {
            Source_RA = 194.953;
                Source_Dec = 27.981;
    }
    if (source_name.Contains("1ES1011"))
    {
            Source_RA = 153.767;
                Source_Dec = 49.434;
    }
    if (source_name.Contains("NGC1275"))
    {
            Source_RA = 49.950;
                Source_Dec = 41.512;
    }
    if (source_name.Contains("1ES0647"))
    {
            Source_RA = 102.694;
                Source_Dec = 25.050;
    }
    if (source_name.Contains("TriII"))
    {
            Source_RA = 34.32861263958333;
                Source_Dec = 36.0;
    }
    if (source_name.Contains("1ES1440"))
    {
            Source_RA = 220.701;
                Source_Dec = 12.011;
    }
    if (source_name.Contains("1ES1741"))
    {
            Source_RA = 266.005;
                Source_Dec = 19.546;
    }
    if (source_name.Contains("IC443HotSpot"))
    {
            //Source_RA = 94.511;
            //    Source_Dec = 22.660;
            Source_RA = 94.213;
                Source_Dec = 22.503;
    }
    if (source_name.Contains("RGBJ0710"))
    {
            Source_RA = 107.610;
                Source_Dec = 59.150;
    }
    if (source_name.Contains("CasA"))
    {
            Source_RA = 350.808;
                Source_Dec = 58.807;
    }
    if (source_name.Contains("M82"))
    {
            Source_RA = 148.970;
                Source_Dec = 69.679;
    }
    if (source_name.Contains("M87"))
    {
            Source_RA = 187.70593076;
                Source_Dec = 12.3911232939;
    }
    if (source_name.Contains("Boomerang"))
    {
            Source_RA = 337.183333333;
                Source_Dec = 61.1666666667;
    }
    if (source_name.Contains("BLLac"))
    {
            Source_RA = 330.680416667;
                Source_Dec = 42.2777777778;
    }
    if (source_name.Contains("UrsaMajorII"))
    {
            Source_RA = 132.875;
                Source_Dec = 63.13;
    }
    if (source_name.Contains("LHAASO_J2108"))
    {
            Source_RA = 317.15;
                Source_Dec = 51.95;
    }
    if (source_name.Contains("LHAASO_J0341"))
    {
            Source_RA = 55.34;
                Source_Dec = 52.97;
    }
    if (source_name.Contains("LHAASO_J1929"))
    {
            Source_RA = 292.25;
                Source_Dec = 17.75;
    }
    if (source_name.Contains("LHAASO_J1843"))
    {
            Source_RA = 280.75;
                Source_Dec = -3.63;
    }
    if (source_name.Contains("LHAASO_J1956"))
    {
            Source_RA = 299.05;
                Source_Dec = 28.75;
    }
    if (source_name.Contains("Perseus"))
    {
            Source_RA = 52.9;
                Source_Dec = 30.9;
    }
    if (source_name.Contains("PSRB0355plus54"))
    {
            Source_RA = 59.72083333333333;
                Source_Dec = 54.22027777777778;
    }
    if (source_name.Contains("1ES0414"))
    {
            Source_RA = 64.2206666667;
                Source_Dec = 1.089;
    }
    if (source_name.Contains("3C273"))
    {
            Source_RA = 187.277915345;
                Source_Dec = 2.05238856846;
    }
    if (source_name.Contains("1ES0502"))
    {
            Source_RA = 76.9839421535;
                Source_Dec = 67.6234172932;
    }
    if (source_name.Contains("Draco"))
    {
            Source_RA = 260.059729167;
                Source_Dec = 57.9212194444;
    }
    if (source_name.Contains("GammaCygni"))
    {
            Source_RA = 305.557091;
                Source_Dec = 40.256679166666665;
    }
    if (source_name.Contains("G079"))
    {
            Source_RA = 308.119;
                Source_Dec = 40.328;
    }
    if (source_name.Contains("WComae"))
    {
            Source_RA = 185.360;
                Source_Dec = 30.191;
            //Source_RA = (185.360+184.616)/2.;
            //    Source_Dec = (30.191+30.130)/2.;
    }
    if (source_name.Contains("1ES1218"))
    {
            Source_RA = 185.360;
                Source_Dec = 30.191;
    }
    if (source_name.Contains("MGRO_J2019"))
    {
            Source_RA = 304.854166667;
                Source_Dec = 36.8038888889;
    }
    if (source_name.Contains("MGRO_J1908"))
    {
            Source_RA = 286.975;
                Source_Dec = 6.03777777778;
            //Source_RA = 287.05;
            //    Source_Dec = 6.39;
    }
    if (source_name.Contains("SS433"))
    {
            Source_RA = 287.9565;
                Source_Dec = 4.98272222222;
    }
    if (source_name.Contains("MAGIC_J1857"))
    {
            Source_RA = 284.303;
                Source_Dec = 2.729;
    }
    if (source_name.Contains("MGRO_J2031"))
    {
            Source_RA = 308.041666667;
                Source_Dec = 41.4594444444;
    }
    if (source_name.Contains("HESS_J1825"))
    {
            Source_RA = 276.37;
                Source_Dec = -13.83;
    }
    if (source_name.Contains("Cygnus"))
    {
            Source_RA = 304.646;
                Source_Dec = 36.833;
    }
    if (source_name.Contains("2HWC_J1953"))
    {
            Source_RA = 298.260;
                Source_Dec = 29.480;
    }
    if (source_name.Contains("2HWC_J1930"))
    {
            Source_RA = 292.150;
                Source_Dec = 17.780;
    }
    if (source_name=="Proton")
    {
            Source_RA = 0.;
                Source_Dec = 0.;
    }
    if (source_name=="Proton_NSB750")
    {
            Source_RA = 0.;
                Source_Dec = 0.;
    }
    return std::make_pair(Source_RA,Source_Dec);
}

TObject* getEffAreaHistogram( TFile* fAnasumDataFile, int runnumber, double offset)
{
  double iSlizeY = -9999;
  //iSlizeY = offset;
  string dirname = "energyHistograms";
  string hisname = "herecEffectiveArea_on";
  if( !fAnasumDataFile )
  {
    return 0;
  }
  
  char dx[600];
  if( runnumber > 1 )
  {
    sprintf( dx, "run_%d/stereo/%s", runnumber, dirname.c_str() );
  }
  else
  {
    if( runnumber == 0 )
    {
      sprintf( dx, "total/stereo/%s", dirname.c_str() );
    }
    else if( runnumber == 1 )
    {
      sprintf( dx, "total_%d/stereo/%s", runnumber, dirname.c_str() );
    }
    else
    {
      sprintf( dx, "total_%d/stereo/%s", -1 * runnumber, dirname.c_str() );
    }
  }
  
  fAnasumDataFile->cd( dx );
  TDirectory* iDir = gDirectory;
  if( !iDir )
  {
    return 0;
  }
  
  TObject* h = ( TObject* )iDir->Get( hisname.c_str() );
  
  if( h && iSlizeY < -9998. )
  {
    return h->Clone();
  }
  else if( h )
  {
    string iClassName = h->ClassName();
    if( iClassName.find( "TH2" ) != string::npos )
    {
      TH2* i_h2 = ( TH2* )h;
      string iN = hisname + "px";
      TH1* i_h = ( TH1* )i_h2->ProjectionX( iN.c_str(), i_h2->GetYaxis()->FindBin( iSlizeY ), i_h2->GetYaxis()->FindBin( iSlizeY ) );
      return i_h->Clone();
    }
  }
  
  
  return 0;
}
bool DarkFoV() {
    if (R2off>9.) return false;
    if (R2off<camera_theta2_cut_lower) return false;
    if (R2off>camera_theta2_cut_upper) return false;
    //if (Yoff<camera_theta2_cut_lower) return false;
    //if (Yoff>camera_theta2_cut_upper) return false;
    //if (CoincideWithBrightStars(ra_sky_dark,dec_sky_dark)) return false;
    if (CoincideWithGammaSources(ra_sky_dark,dec_sky_dark)) return false;
    return true;
}
bool MCFoV() {
    //double angular_dist = pow(pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2),0.5);
    //if (angular_dist>1.5) return false;
    return true;
}
bool FoV(bool doImposter) {

    if (R2off>9.) return false;
    if (R2off<camera_theta2_cut_lower) return false;
    if (R2off>camera_theta2_cut_upper) return false;
    
    double x = ra_sky-mean_tele_point_ra;
    double y = dec_sky-mean_tele_point_dec;
    if (abs(x)>Skymap_size) return false;
    if (abs(y)>Skymap_size) return false;
    if (source_theta2_cut>(x*x+y*y)) return false;
    if (x<map_x_bin_lower) return false;
    if (x>map_x_bin_upper) return false;
    if (y<map_y_bin_lower) return false;
    if (y>map_y_bin_upper) return false;
    //if (CoincideWithBrightStars(ra_sky,dec_sky)) return false;
    //if (CoincideWithGammaSources(ra_sky,dec_sky)) return false;
    if (doImposter)
    {
        if (CoincideWithGammaSources(ra_sky_imposter,dec_sky_imposter)) return false;
    } 
    
    //vector<double> ss433_ra;
    //vector<double> ss433_dec;
    //ss433_ra.push_back(287.9565);
    //ss433_dec.push_back(4.9827);
    //ss433_ra.push_back(288.404);
    //ss433_dec.push_back(4.930);
    //ss433_ra.push_back(288.58);
    //ss433_dec.push_back(4.91);
    //ss433_ra.push_back(287.42);
    //ss433_dec.push_back(5.04);
    //for (int star=0;star<ss433_ra.size();star++)
    //{
    //    double star_ra = ss433_ra.at(star);
    //    double star_dec = ss433_dec.at(star);
    //    double radius = pow((ra_sky-star_ra)*(ra_sky-star_ra)+(dec_sky-star_dec)*(dec_sky-star_dec),0.5);
    //    if (radius<0.3) return false;
    //}

    return true;
}
bool SourceFoV() {

    if (R2off>9.) return false;
    if (R2off<camera_theta2_cut_lower) return false;
    if (R2off>camera_theta2_cut_upper) return false;
    
    double x = ra_sky-mean_tele_point_ra;
    double y = dec_sky-mean_tele_point_dec;
    if ((x*x+y*y)>0.3*0.3) return false;
    return true;
}
bool SourceRingFoV() {

    if (R2off>9.) return false;
    if (R2off<camera_theta2_cut_lower) return false;
    if (R2off>camera_theta2_cut_upper) return false;
    
    double x = ra_sky-mean_tele_point_ra;
    double y = dec_sky-mean_tele_point_dec;
    if ((x*x+y*y)<0.3*0.3) return false;
    if ((x*x+y*y)>2.*0.3*0.3) return false;
    return true;
}
bool RoIFoV(int which_roi) {
    double x = ra_sky-roi_ra.at(which_roi);
    double y = dec_sky-roi_dec.at(which_roi);
    double radius = pow(x*x+y*y,0.5);
    if (radius<roi_radius_inner.at(which_roi)) return false;
    if (radius>roi_radius_outer.at(which_roi)) return false;
    return true;
}
double GetShowerDepth(double shower_height, double elevation)
{
    double atm_depth = 25.; //km
    double earth_radius = 6371.; //km
    double zenith_rad = M_PI/180.*(90.-elevation);
    double shower_depth = pow(pow(earth_radius+atm_depth,2)-pow(earth_radius*sin(zenith_rad),2),0.5)-earth_radius*cos(zenith_rad)-shower_height;
    return shower_depth;
}
double RescaleMSCW(double MSCW_input, double R2off_input, double MSCW_shift)
{
    //double MSCW_new = MSCW_input-MSCW_shift*min(R2off_input,1.0);
    double MSCW_new = MSCW_input*(1.+MSCW_shift);
    return MSCW_new;
}
bool SelectNImages()
{
    if (NImages<2) return false;
    return true;
}
bool ApplyTimeCuts(double event_time, vector<pair<double,double>> timecut)
{
    bool pass_cut = true;
    for (int cut=0;cut<timecut.size();cut++)
    {
        double timecut_start = timecut.at(cut).first;
        double timecut_end = timecut.at(cut).second;
        //std::cout << "event_time " << event_time << ", timecut_start " << timecut_start << ", timecut_end " << timecut_end << std::endl;
        if (event_time>timecut_start && event_time<timecut_end)
        {
            pass_cut = false;
        }
    }
    return pass_cut;
}
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
        if (!pointing_tree)
        {
            return 0.;
        }
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
    if (!pointing_tree)
    {
        return 0;
    }
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
    if (!pointing_tree)
    {
        return false;
    }
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
    if (!pointing_tree)
    {
        return std::make_pair(0.,0.);
    }
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

void SetEventDisplayTreeBranch(TTree* Data_tree)
{
    Data_tree->SetBranchStatus("*",0);
    Data_tree->SetBranchStatus("Xoff",1);
    Data_tree->SetBranchStatus("Yoff",1);
    Data_tree->SetBranchStatus("theta2",1);
    Data_tree->SetBranchStatus("Xoff_derot",1);
    Data_tree->SetBranchStatus("Yoff_derot",1);
    Data_tree->SetBranchStatus("ErecS",1);
    Data_tree->SetBranchStatus("EChi2S",1);
    Data_tree->SetBranchStatus("MSCW",1);
    Data_tree->SetBranchStatus("MSCL",1);
    Data_tree->SetBranchStatus("NImages",1);
    Data_tree->SetBranchStatus("Xcore",1);
    Data_tree->SetBranchStatus("Ycore",1);
    Data_tree->SetBranchStatus("SizeSecondMax",1);
    Data_tree->SetBranchStatus("EmissionHeight",1);
    Data_tree->SetBranchStatus("Time",1);
    Data_tree->SetBranchStatus("Shower_Ze",1);
    Data_tree->SetBranchStatus("Shower_Az",1);
    Data_tree->SetBranchStatus("MJD",1);
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
}

void GetRunCosmicRayAcceptance(string file_name, int run, TH2D* hist_input)
{
    TH2D hist_acc = TH2D("hist_acc","",8,-2,2,8,-2,2);
    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());
    TTree* event_tree = nullptr;
    event_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/data_on"));
    SetEventDisplayTreeBranch(event_tree);
    for (int entry=0;entry<event_tree->GetEntries();entry++) 
    {
        ErecS = 0;
        NImages = 0;
        SizeSecondMax = 0;
        MSCW = 0;
        MSCL = 0;
        Xoff = 0;
        Yoff = 0;
        event_tree->GetEntry(entry);
        if (!SelectNImages()) continue;
        if (SizeSecondMax<SizeSecondMax_Cut) continue;
        if (MSCW<0.7 && MSCL<0.7) continue;
        hist_acc.Fill(Xoff,Yoff);
    }
    input_file->Close();
    hist_input->Add(&hist_acc);
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
    if (!pointing_tree)
    {
        return std::make_pair(0.,0.);
    }
    pointing_tree->SetBranchStatus("*",0);
    pointing_tree->SetBranchStatus("TelElevation",1);
    pointing_tree->SetBranchStatus("TelAzimuth",1);
    pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
    pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
    double total_entries = (double)pointing_tree->GetEntries();
    //for (int entry=0;entry<total_entries;entry++)
    //{
    //    pointing_tree->GetEntry(entry);
    //    TelElevation_avg += TelElevation;
    //    TelAzimuth_avg += TelAzimuth;
    //}
    //TelElevation_avg = TelElevation_avg/total_entries;
    //TelAzimuth_avg = TelAzimuth_avg/total_entries;
    pointing_tree->GetEntry(int(total_entries/2.));
    TelElevation_avg = TelElevation;
    TelAzimuth_avg = TelAzimuth;
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
    if (!pointing_tree)
    {
        return false;
    }
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
vector<pair<string,int>> SelectONRunList(vector<pair<string,int>> Data_runlist, double Elev_cut_lower, double Elev_cut_upper, int MJD_start_cut, int MJD_end_cut)
{
    std::cout << "initial runs = " << Data_runlist.size() << std::endl;
    vector<pair<string,int>> new_list;
    vector<double> new_list_pointing;
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

        new_list.push_back(std::make_pair(Data_runlist[run].first,Data_runlist[run].second));
        new_list_pointing.push_back(GetRunElevAzim(filename,int(Data_runlist[run].second)).first);
    }
    std::cout << "selected runs = " << new_list.size() << std::endl;
    if (new_list.size()>1) SortingList(&new_list, &new_list_pointing);
    return new_list;
}

double GetHistogramChi2(TH2D* hist_1, TH2D* hist_2)
{
    double chi2 = 0.;
    double dof = double(hist_1->GetNbinsX())*double(hist_1->GetNbinsY());
    for (int binx=0;binx<hist_1->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist_1->GetNbinsY();biny++)
        {
            double hist_1_count = hist_1->GetBinContent(binx+1,biny+1);
            double hist_2_count = hist_2->GetBinContent(binx+1,biny+1);
            double sigma = max(1.,pow(hist_1_count,0.5));
            chi2 += pow((hist_1_count-hist_2_count)/sigma,2);
        }
    }
    return chi2/dof;
}

vector<vector<vector<pair<string,int>>>> SelectDarkRunList(vector<pair<string,int>> ON_runlist, vector<pair<string,int>> OFF_runlist_input, double tel_elev_lower, double tel_elev_upper, bool nsb_reweight)
{

    int nsb_bins = 1;
    if (nsb_reweight) nsb_bins = 20;
    TH2D Hist_OnData_ElevNSB = TH2D("Hist_OnData_ElevNSB","",nsb_bins,0,10,18,0,90);
    TH2D Hist_OffData_ElevNSB = TH2D("Hist_OffData_ElevNSB","",nsb_bins,0,10,18,0,90);

    std::cout << "Load ON run info" << std::endl;
    vector<pair<double,double>> ON_pointing;
    vector<pair<double,double>> ON_pointing_radec;
    vector<double> ON_count;
    vector<double> ON_time;
    vector<vector<double>> Dark_count;
    vector<double> ON_NSB;
    vector<double> ON_MJD;
    vector<double> ON_L3Rate;
    for (int on_run=0;on_run<ON_runlist.size();on_run++)
    {
        char ON_runnumber[50];
        char ON_observation[50];
        sprintf(ON_runnumber, "%i", int(ON_runlist[on_run].second));
        sprintf(ON_observation, "%s", ON_runlist[on_run].first.c_str());
        string ON_filename;
        ON_filename = TString(SMI_INPUT+"/"+string(ON_runnumber)+".anasum.root");
        if (TString(ON_observation).Contains("Proton")) ON_pointing.push_back(std::make_pair(70,0));
        else ON_pointing.push_back(GetRunElevAzim(ON_filename,int(ON_runlist[on_run].second)));
        if (TString(ON_observation).Contains("Proton")) ON_pointing_radec.push_back(std::make_pair(0,0));
        else ON_pointing_radec.push_back(GetRunRaDec(ON_filename,int(ON_runlist[on_run].second)));
        double NSB_thisrun = GetRunPedestalVar(int(ON_runlist[on_run].second));
        ON_NSB.push_back(NSB_thisrun);
        double exposure_thisrun = GetRunUsableTime(ON_filename,ON_runlist[on_run].second);
        ON_time.push_back(exposure_thisrun);
        ON_MJD.push_back(GetRunMJD(ON_filename,int(ON_runlist[on_run].second)));
        ON_L3Rate.push_back(GetRunL3Rate(ON_runlist[on_run].second));
        ON_count.push_back(exposure_thisrun*GetRunL3Rate(ON_runlist[on_run].second));
        if (MJD<MJD_Start) MJD_Start = MJD;
        if (MJD>MJD_End) MJD_End = MJD;
        Hist_OnData_ElevNSB.Fill(ON_NSB[ON_NSB.size()-1],ON_pointing[ON_pointing.size()-1].first,exposure_thisrun);
    }

    std::cout << "Load OFF run info" << std::endl;
    vector<pair<string,int>> OFF_runlist;
    vector<pair<double,double>> OFF_pointing;
    vector<pair<double,double>> OFF_pointing_radec;
    vector<double> OFF_time;
    vector<double> OFF_MJD;
    vector<double> OFF_L3Rate;
    vector<double> OFF_NSB;
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
        if (run_elevation<tel_elev_lower-5.) continue;
        if (run_elevation>tel_elev_upper+5.) continue;
        OFF_runlist.push_back(OFF_runlist_input[off_run]);
        if (TString(OFF_observation).Contains("Proton")) OFF_pointing.push_back(std::make_pair(70,0));
        else OFF_pointing.push_back(GetRunElevAzim(OFF_filename,int(OFF_runlist_input[off_run].second)));
        if (TString(OFF_observation).Contains("Proton")) OFF_pointing_radec.push_back(std::make_pair(0,0));
        else OFF_pointing_radec.push_back(GetRunRaDec(OFF_filename,int(OFF_runlist_input[off_run].second)));
        double NSB_thisrun = GetRunPedestalVar(int(OFF_runlist_input[off_run].second));
        OFF_NSB.push_back(NSB_thisrun);
        double exposure_thisrun = GetRunUsableTime(OFF_filename,OFF_runlist_input[off_run].second);
        OFF_time.push_back(exposure_thisrun);
        OFF_MJD.push_back(GetRunMJD(OFF_filename,int(OFF_runlist_input[off_run].second)));
        OFF_L3Rate.push_back(GetRunL3Rate(OFF_runlist_input[off_run].second));
        Hist_OffData_ElevNSB.Fill(OFF_NSB[OFF_NSB.size()-1],OFF_pointing[OFF_pointing.size()-1].first,exposure_thisrun);
    }

    std::cout << "Select matched runs" << std::endl;
    vector<vector<vector<pair<string,int>>>> new_list;
    for (int on_run=0;on_run<ON_runlist.size();on_run++)
    {
        vector<vector<pair<string,int>>> the_samples;
        vector<double> Dark_count_thisrun;
        for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
        {
            vector<pair<string,int>> the_runs;
            the_samples.push_back(the_runs);
            Dark_count_thisrun.push_back(0.);
        }
        new_list.push_back(the_samples);
        Dark_count.push_back(Dark_count_thisrun);
    }
    vector<pair<double,double>> ON_pointing_radec_new;
    TH2D hist_on_acc = TH2D("hist_on_acc","",8,-2,2,8,-2,2);
    TH2D hist_off_acc = TH2D("hist_off_acc","",8,-2,2,8,-2,2);
    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
    {
        std::cout << "Complete " << nth_sample << "/" << n_dark_samples << std::endl;
        double total_dNSB_chi2 = 0.;
        double total_dElev_chi2 = 0.;
        double total_dMJD_chi2 = 0.;
        double total_dL3Rate_chi2 = 0.;
        double dNSB_chi2 = 0.;
        double dElev_chi2 = 0.;
        double dMJD_chi2 = 0.;
        double dL3Rate_chi2 = 0.;
        for (int on_run=0;on_run<ON_runlist.size();on_run++)
        {
            hist_on_acc.Reset();
            char ON_runnumber[50];
            sprintf(ON_runnumber, "%i", int(ON_runlist[on_run].second));
            string ON_filename;
            ON_filename = TString(SMI_INPUT+"/"+string(ON_runnumber)+".anasum.root");
            //GetRunCosmicRayAcceptance(ON_filename, int(ON_runlist[on_run].second), &hist_on_acc);
            std::cout << "finding matches for " << on_run << "/" << ON_runlist.size() << " runs..." << std::endl;
            std::cout << "ON_pointing[on_run].first = " << ON_pointing[on_run].first << std::endl;
            std::cout << "ON_count[on_run] = " << ON_count[on_run] << std::endl;
            std::cout << "ON_time[on_run] = " << ON_time[on_run] << std::endl;
            std::cout << "ON_L3Rate[on_run] = " << ON_L3Rate[on_run] << std::endl;
            std::cout << "ON_NSB[on_run] = " << ON_NSB[on_run] << std::endl;
            int number_of_search = 0;
            double accumulated_count = 0.;
            double offset_NSB = 0.;
            double offset_Elev = 0.;
            double threshold_dNSB = 0.5;
            double threshold_dElev = 2.0;
            double threshold_dAzim = 45.;
            double threshold_dMJD = 3.*365.;
            double threshold_dL3Rate = 0.3;
            double threshold_dTime = 15.*60.;
            while (accumulated_count<2.0*ON_count[on_run])
            {
                pair<string,int> best_match;
                pair<double,double> best_pointing;
                double best_chi2 = 1e10;
                int best_off_run = 0;
                double best_count = 0.;
                bool found_match = false;
                for (int off_run=0;off_run<OFF_runlist.size();off_run++)
                {

                    if (RunTypeCategory(OFF_runlist[off_run].second,false)!=0) continue;
                    double diff_ra = ON_pointing_radec[on_run].first-OFF_pointing_radec[off_run].first;
                    double diff_dec = ON_pointing_radec[on_run].second-OFF_pointing_radec[off_run].second;
                    if ((diff_ra*diff_ra+diff_dec*diff_dec)<10.*10.) continue;

                    //M87
                    double M87_RA = 187.70593076;
                    double M87_Dec = 12.3911232939;
                    diff_ra = M87_RA-OFF_pointing_radec[off_run].first;
                    diff_dec = M87_Dec-OFF_pointing_radec[off_run].second;
                    if ((diff_ra*diff_ra+diff_dec*diff_dec)<10.*10.) continue;

                    if (ON_runlist[on_run].first.find("Proton")==std::string::npos)
                    {
                        if (ON_runlist[on_run].first.compare(OFF_runlist[off_run].first) == 0) continue;
                    }
                    bool already_used_run = false;
                    for (int on_run2=0;on_run2<ON_runlist.size();on_run2++)
                    {
                        if (int(ON_runlist[on_run2].second)==int(OFF_runlist[off_run].second)) already_used_run = true; // this OFF run is in ON runlist
                    }
                    for (int nth_other_sample=0;nth_other_sample<n_dark_samples;nth_other_sample++)
                    {
                        for (int new_run=0;new_run<new_list.at(on_run).at(nth_other_sample).size();new_run++)
                        {
                            if (int(new_list.at(on_run).at(nth_other_sample)[new_run].second)==int(OFF_runlist[off_run].second)) already_used_run = true;
                        }
                    }
                    if (already_used_run) continue;

                    double delta_azimuth = abs(ON_pointing[on_run].second-OFF_pointing[off_run].second);
                    if (delta_azimuth>180.) delta_azimuth = 360.-delta_azimuth;

                    double chi2 = 0.;
                    found_match = true;
                    //chi2 = pow(ON_pointing[on_run].first-OFF_pointing[off_run].first+offset_Elev,2);
                    chi2 = pow(delta_azimuth,2);
                    //chi2 = pow(ON_L3Rate[on_run]-OFF_L3Rate[off_run],2);
                    //chi2 = pow(ON_NSB[on_run]-OFF_NSB[off_run]-offset_NSB,2);
                    if (pow(ON_pointing[on_run].first-OFF_pointing[off_run].first+offset_Elev,2)>threshold_dElev*threshold_dElev)
                    {
                        found_match = false;
                        continue;
                    }
                    //if (pow(delta_azimuth,2)>threshold_dAzim*threshold_dAzim)
                    //{
                    //    found_match = false;
                    //    continue;
                    //}
                    if (pow(ON_NSB[on_run]-OFF_NSB[off_run]-offset_NSB,2)>threshold_dNSB*threshold_dNSB)
                    {
                        found_match = false;
                        continue;
                    }
                    if (pow(ON_MJD[on_run]-OFF_MJD[off_run],2)>threshold_dMJD*threshold_dMJD)
                    {
                        found_match = false;
                        continue;
                    }
                    if (pow((ON_L3Rate[on_run]-OFF_L3Rate[off_run])/ON_L3Rate[on_run],2)>threshold_dL3Rate*threshold_dL3Rate)
                    {
                        found_match = false;
                        continue;
                    }
                    if (pow(ON_time[on_run]-OFF_time[off_run],2)>threshold_dTime*threshold_dTime)
                    {
                        found_match = false;
                        continue;
                    }

                    hist_off_acc.Reset();
                    char OFF_runnumber[50];
                    sprintf(OFF_runnumber, "%i", int(OFF_runlist[off_run].second));
                    string OFF_filename;
                    OFF_filename = TString(SMI_INPUT+"/"+string(OFF_runnumber)+".anasum.root");
                    //GetRunCosmicRayAcceptance(OFF_filename, int(OFF_runlist[off_run].second), &hist_off_acc);
                    //std::cout << "hist_on_acc.Integral() = " << hist_on_acc.Integral() << std::endl;
                    //std::cout << "hist_off_acc.Integral() = " << hist_off_acc.Integral() << std::endl;
                    //double on_off_scale = hist_on_acc.Integral()/hist_off_acc.Integral();
                    //hist_off_acc.Scale(on_off_scale);
                    //double chi2_acc = GetHistogramChi2(&hist_on_acc, &hist_off_acc);
                    //std::cout << "chi2_acc = " << chi2_acc << std::endl;
                    //chi2 = chi2_acc;

                    if (found_match && best_chi2>chi2)
                    {
                        best_chi2 = chi2;
                        best_match = OFF_runlist[off_run];
                        best_pointing = OFF_pointing[off_run];
                        best_off_run = off_run;
                        best_count = OFF_time[off_run]*OFF_L3Rate[off_run];
                        dElev_chi2 = pow(ON_pointing[on_run].first-OFF_pointing[off_run].first,2);
                        dNSB_chi2 = pow(ON_NSB[on_run]-OFF_NSB[off_run],2);
                        dMJD_chi2 = pow(ON_MJD[on_run]-OFF_MJD[off_run],2);
                        dL3Rate_chi2 = pow(ON_L3Rate[on_run]-OFF_L3Rate[off_run],2);
                    }
                    //if (best_chi2<1.0) break;
                }
                if (best_chi2<1e10) 
                {
                    std::cout << on_run << "/" << ON_runlist.size() << ", found a match " << best_match.first << " " << best_match.second << ", best_chi2 = " << best_chi2 << std::endl;
                    //std::cout << "OFF_pointing[off_run].first = " << OFF_pointing[best_off_run].first << std::endl;
                    //std::cout << "OFF_time[off_run] = " << OFF_time[best_off_run] << std::endl;
                    //std::cout << "OFF_L3Rate[off_run] = " << OFF_L3Rate[best_off_run] << std::endl;
                    //std::cout << "OFF_NSB[off_run] = " << OFF_NSB[best_off_run] << std::endl;
                    new_list.at(on_run).at(nth_sample).push_back(best_match);
                    ON_pointing_radec_new.push_back(ON_pointing_radec[on_run]);
                    n_good_matches += 1;
                    accumulated_count += best_count;
                    total_dNSB_chi2 += dNSB_chi2;
                    total_dElev_chi2 += dElev_chi2;
                    total_dMJD_chi2 += dMJD_chi2;
                    total_dL3Rate_chi2 += dL3Rate_chi2;
                }
                else
                {
                    // searched whole OFF list and found no match.
                    number_of_search += 1;
                    std::cout << on_run << "/" << ON_runlist.size() << ", number_of_search = " << number_of_search << std::endl;
                    if (number_of_search>=1)
                    {
                        std::cout << "couldn't find a matched run, relax time cuts." << std::endl;
                        threshold_dMJD += 10.0*365.;
                        threshold_dTime += 5.0*60.;
                    }
                    if (number_of_search>=2)
                    {
                        std::cout << "couldn't find a matched run, relax elevation and NSB cuts." << std::endl;
                        threshold_dNSB += 0.5;
                        threshold_dElev += 2.;
                    }
                    if (number_of_search>=3)
                    {
                        std::cout << "couldn't find a matched run for " << int(ON_runlist[on_run].second) << ", break." << std::endl;
                        std::cout << "ON run elevation " << ON_pointing[on_run].first << ", azimuth " << ON_pointing[on_run].second << ", NSB " << ON_NSB[on_run] << std::endl; 
                        break;
                        std::cout << "couldn't find a matched run, relax all cuts." << std::endl;
                        threshold_dNSB += 10.;
                        threshold_dElev += 10.;
                        threshold_dMJD += 10.0*365.;
                        threshold_dL3Rate += 0.3;
                        threshold_dTime += 5.0*60.;
                    }
                }
            }
            Dark_count.at(on_run).at(nth_sample) = accumulated_count;
        }
        total_dNSB_chi2 = pow(total_dNSB_chi2/double(OFF_runlist.size()),0.5);
        total_dElev_chi2 = pow(total_dElev_chi2/double(OFF_runlist.size()),0.5);
        total_dMJD_chi2 = pow(total_dMJD_chi2/double(OFF_runlist.size()),0.5);
        total_dL3Rate_chi2 = pow(total_dL3Rate_chi2/double(OFF_runlist.size()),0.5);
        std::cout << "total_dNSB_chi2 = " << total_dNSB_chi2 << std::endl;
        std::cout << "total_dElev_chi2 = " << total_dElev_chi2 << std::endl;
        std::cout << "total_dMJD_chi2 = " << total_dMJD_chi2 << std::endl;
        std::cout << "total_dL3Rate_chi2 = " << total_dL3Rate_chi2 << std::endl;
    }
    for (int on_run=0;on_run<ON_runlist.size();on_run++)
    {
        vector<double> Dark_weight_thisrun;
        for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
        {
            double weight = 0.;
            if (Dark_count.at(on_run).at(nth_sample)>0.)
            {
                weight = ON_count.at(on_run)/Dark_count.at(on_run).at(nth_sample);
            }
            Dark_weight_thisrun.push_back(weight);
        }
        Dark_weight.push_back(Dark_weight_thisrun);
    }

    return new_list;

}


void GetGammaSources()
{
    std::ifstream astro_file(SMI_AUX+"/TeVCat_RaDec.txt");
    std::string line;
    // Read one line at a time into the variable line:
    while(std::getline(astro_file, line))
    {
        if (line.empty()) continue;
        std::vector<double>   lineData;
        std::stringstream  lineStream(line);
        double value;
        // Read an integer at a time from the line
        while(lineStream >> value)
        {
            // Add the integers from a line to a 1D array (vector)
            lineData.push_back(value);
        }
        if (lineData.size()!=2) continue;
        double star_ra = lineData.at(0);
        double star_dec = lineData.at(1);
        GammaSource_Data.push_back(lineData);
    }
    std::cout << "I found " << GammaSource_Data.size() << " gamma-ray sources." << std::endl;
}

void GetBrightStars()
{
    std::ifstream astro_file(SMI_AUX+"/Hipparcos_MAG8_1997.dat");
    std::string line;
    // Read one line at a time into the variable line:
    while(std::getline(astro_file, line))
    {
        if (line.find("#")!=std::string::npos) continue;
        if (line.find("*")!=std::string::npos) continue;
        if (line.empty()) continue;
        std::vector<double>   lineData;
        std::stringstream  lineStream(line);
        double value;
        // Read an integer at a time from the line
        while(lineStream >> value)
        {
            // Add the integers from a line to a 1D array (vector)
            lineData.push_back(value);
        }
        if (lineData.size()!=5) continue;
        double star_ra = lineData.at(0);
        double star_dec = lineData.at(1);
        double star_brightness = lineData.at(3)+lineData.at(4);
        if (pow((mean_tele_point_ra-star_ra)*(mean_tele_point_ra-star_ra)+(mean_tele_point_dec-star_dec)*(mean_tele_point_dec-star_dec),0.5)>Skymap_size) continue;
        if (star_brightness<brightness_cut)
        {
            BrightStars_Data.push_back(lineData);
        }
        if (star_brightness<faint_brightness_cut)
        {
            FaintStars_Data.push_back(lineData);
        }
    }
    std::cout << "I found " << BrightStars_Data.size() << " bright stars" << std::endl;
}

bool SignalSelectionTheta2()
{
    if (MSCW>MSCW_cut_blind) return false;
    if (MSCW<-1.*MSCW_cut_blind) return false;
    if (MSCL>MSCL_cut_blind) return false;
    if (MSCL<-1.*MSCL_cut_blind) return false;
    return true;
}
bool ControlSelectionTheta2()
{
    if (SignalSelectionTheta2()) return false;
    if (MSCW<MSCW_cut_blind && MSCL<MSCL_cut_blind) return false;
    //if (MSCW<MSCW_cut_blind || MSCL<MSCL_cut_blind) return false;
    if (MSCW<MSCW_cut_blind) return false;
    double boundary = 0.8;
    if (MSCW>boundary*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind) return false;

    //if (MSCW<MSCW_cut_blind) return false;
    //if (MSCL>MSCL_cut_blind) return false;
    //double boundary = 1.0;
    //if (MSCL>boundary*gamma_hadron_dim_ratio_l[0]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind) return false;
    //if (MSCW>boundary*gamma_hadron_dim_ratio_w[0]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind) return false;
    return true;
}

void Smooth2DMap(TH2D* hist, double smooth_size)
{

    TH2D Hist_Unit = TH2D("Hist_Unit","",Skymap_nbins/3,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins/3,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size);
    TH2D Hist_Smooth = TH2D("Hist_Smooth","",Skymap_nbins/3,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins/3,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size);
    double bin_size = Hist_Unit.GetXaxis()->GetBinCenter(2)-Hist_Unit.GetXaxis()->GetBinCenter(1);
    int nbin_smooth = int(2*smooth_size/bin_size) + 1;
    for (int bx1=1;bx1<=Hist_Unit.GetNbinsX();bx1++)
    {
        for (int by1=1;by1<=Hist_Unit.GetNbinsY();by1++)
        {
            double bin_content = 0;
            double bin_error = 0;
            double bin_norm = 0;
            double locationx1 = Hist_Unit.GetXaxis()->GetBinCenter(bx1);
            double locationy1 = Hist_Unit.GetYaxis()->GetBinCenter(by1);
            for (int bx2=bx1-nbin_smooth;bx2<=bx1+nbin_smooth;bx2++)
            {
                for (int by2=by1-nbin_smooth;by2<=by1+nbin_smooth;by2++)
                {
                    if (bx2<1) continue;
                    if (by2<1) continue;
                    if (bx2>Hist_Unit.GetNbinsX()) continue;
                    if (by2>Hist_Unit.GetNbinsY()) continue;
                    double locationx2 = Hist_Unit.GetXaxis()->GetBinCenter(bx2);
                    double locationy2 = Hist_Unit.GetYaxis()->GetBinCenter(by2);
                    double distance = pow(pow(locationx1-locationx2,2)+pow(locationy1-locationy2,2),0.5);
                    bin_content += TMath::Gaus(distance,0,smooth_size)*hist->GetBinContent(bx2,by2);
                    bin_norm += TMath::Gaus(distance,0,smooth_size);
                }
            }
            Hist_Smooth.SetBinContent(bx1,by1,bin_content);
            Hist_Unit.SetBinContent(bx1,by1,bin_norm);
        }
    }
    Hist_Smooth.Divide(&Hist_Unit);
    for (int bx1=1;bx1<=Hist_Unit.GetNbinsX();bx1++)
    {
        for (int by1=1;by1<=Hist_Unit.GetNbinsY();by1++)
        {
            hist->SetBinContent(bx1,by1,Hist_Smooth.GetBinContent(bx1,by1));
        }
    }

}

void PrepareDarkData_SubGroup(string target_data, double tel_elev_lower_input, double tel_elev_upper_input, int MJD_start_cut, int MJD_end_cut, double input_theta2_cut_lower, double input_theta2_cut_upper, bool isON, bool doImposter, int GammaModel, int group_index, int map_x_index, int map_y_index)
{

    SMI_INPUT = string(std::getenv("SMI_INPUT"));
    SMI_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    SMI_DIR = string(std::getenv("SMI_DIR"));
    SMI_AUX = string(std::getenv("SMI_AUX"));

    TH1::SetDefaultSumw2();

    sprintf(group_tag, "_G%d", group_index);
    sprintf(map_x_tag, "_X%d", map_x_index);
    sprintf(map_y_tag, "_Y%d", map_y_index);
    map_x_bin_upper = double(map_x_index+1)*2.0*Skymap_size/double(Skymap_normalization_nbins)-Skymap_size;
    map_x_bin_lower = double(map_x_index)*2.0*Skymap_size/double(Skymap_normalization_nbins)-Skymap_size;
    map_y_bin_upper = double(map_y_index+1)*2.0*Skymap_size/double(Skymap_normalization_nbins)-Skymap_size;
    map_y_bin_lower = double(map_y_index)*2.0*Skymap_size/double(Skymap_normalization_nbins)-Skymap_size;

    roi_name.clear();
    roi_ra.clear();
    roi_dec.clear();
    roi_radius_inner.clear();
    roi_radius_outer.clear();
    BrightStars_Data.clear();
    FaintStars_Data.clear();
    GammaSource_Data.clear();
    Dark_weight.clear();
    exposure_hours = 0.;
    exposure_hours_usable = 0.;

    std::cout << "Get systematic error histograms" << std::endl;
    vector<TH1D> Hist_NormSystErr;
    for (int elev=0;elev<N_elev_bins;elev++) 
    {
        char elev_tag[50];
        sprintf(elev_tag,"TelElev%ito%i",int(elev_bins[elev]),int(elev_bins[elev+1]));
        Hist_NormSystErr.push_back(TH1D("Hist_NormSystErr_Temp_"+TString(elev_tag),"",N_energy_bins,energy_bins));
    }
    vector<vector<TH2D>> Hist_ShapeSystErr;
    vector<vector<TH1D>> Hist_ShapeSystErr_1D;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_tag[50];
        sprintf(e_tag,"ErecS%ito%i",int(energy_bins[e]),int(energy_bins[e+1]));
        vector<TH2D> Hist_ShapeSystErr_ThisBin;
        vector<TH1D> Hist_ShapeSystErr_1D_ThisBin;
        for (int xybin=0;xybin<N_integration_radii;xybin++) 
        {
            char xybin_tag[50];
            sprintf(xybin_tag,"Bin%i",xybin);
            Hist_ShapeSystErr_ThisBin.push_back(TH2D("Hist_ShapeSystErr_Temp_"+TString(e_tag)+TString("_")+TString(xybin_tag),"",36,-3,3,36,-3,3));
            Hist_ShapeSystErr_1D_ThisBin.push_back(TH1D("Hist_ShapeSystErr_1D_Temp_"+TString(e_tag)+TString("_")+TString(xybin_tag),"",5,0.,2.));
        }
        Hist_ShapeSystErr.push_back(Hist_ShapeSystErr_ThisBin);
        Hist_ShapeSystErr_1D.push_back(Hist_ShapeSystErr_1D_ThisBin);
    }

    TFile InputSystFile(TString(SMI_AUX)+"/SystErrors.root");
    TString hist_name;
    for (int elev=0;elev<N_elev_bins;elev++) 
    {
        char elev_tag[50];
        sprintf(elev_tag,"TelElev%ito%i",int(elev_bins[elev]),int(elev_bins[elev+1]));
        hist_name  = "Hist_NormSystErr_"+TString(elev_tag);
        Hist_NormSystErr.at(elev).Add( (TH1D*)InputSystFile.Get(hist_name) );
        std::cout << "Get hist " << hist_name << std::endl;
        Hist_NormSystErr.at(elev).Print("All");
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_tag[50];
        sprintf(e_tag,"ErecS%ito%i",int(energy_bins[e]),int(energy_bins[e+1]));
        for (int xybin=0;xybin<N_integration_radii;xybin++) 
        {
            char xybin_tag[50];
            sprintf(xybin_tag,"Bin%i",xybin);
            hist_name  = "Hist_ShapeSystErr_"+TString(e_tag)+TString("_")+TString(xybin_tag);
            std::cout << "Get hist " << hist_name << std::endl;
            Hist_ShapeSystErr.at(e).at(xybin).Add( (TH2D*)InputSystFile.Get(hist_name) );
            hist_name  = "Hist_ShapeSystErr_1D_"+TString(e_tag)+TString("_")+TString(xybin_tag);
            std::cout << "Get hist " << hist_name << std::endl;
            Hist_ShapeSystErr_1D.at(e).at(xybin).Add( (TH1D*)InputSystFile.Get(hist_name) );
        }
    }


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

    SizeSecondMax_Cut = 600.;
    if (TString(target).Contains("V4")) SizeSecondMax_Cut = 400.;
    if (TString(target).Contains("V5")) SizeSecondMax_Cut = 400.;


    mean_tele_point_ra = 0.;
    mean_tele_point_dec = 0.;
    pair<double,double> source_ra_dec = GetSourceRaDec(TString(target));
    mean_tele_point_ra = source_ra_dec.first;
    mean_tele_point_dec = source_ra_dec.second;
    pair<double,double> mean_tele_point_l_b = ConvertRaDecToGalactic(mean_tele_point_ra, mean_tele_point_dec);
    mean_tele_point_l = mean_tele_point_l_b.first;
    mean_tele_point_b = mean_tele_point_l_b.second;

    GetBrightStars();
    GetGammaSources();

    if (!isON) 
    {
        roi_name.push_back("Central region");
        roi_ra.push_back(mean_tele_point_ra);
        roi_dec.push_back(mean_tele_point_dec);
        roi_radius_inner.push_back(0.0);
        roi_radius_outer.push_back(1.0);

        roi_name.push_back("Control region");
        roi_ra.push_back(mean_tele_point_ra);
        roi_dec.push_back(mean_tele_point_dec);
        roi_radius_inner.push_back(1.0);
        roi_radius_outer.push_back(10.);
    }
    else
    {
        roi_name.push_back("Central region");
        roi_ra.push_back(mean_tele_point_ra);
        roi_dec.push_back(mean_tele_point_dec);
        roi_radius_inner.push_back(0.0);
        roi_radius_outer.push_back(1.0);

        roi_name.push_back("Control region");
        roi_ra.push_back(mean_tele_point_ra);
        roi_dec.push_back(mean_tele_point_dec);
        roi_radius_inner.push_back(1.0);
        roi_radius_outer.push_back(10.);

        if (TString(target).Contains("MGRO_J1908")) 
        {
            roi_name.push_back("VHE region");
            roi_ra.push_back(286.84);
            roi_dec.push_back(6.22);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.5);

            roi_name.push_back("PSR region");
            roi_ra.push_back(286.98);
            roi_dec.push_back(6.04);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.15);

            roi_name.push_back("PSR tail 1");
            roi_ra.push_back(287.09);
            roi_dec.push_back(6.32);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.15);

            roi_name.push_back("PSR tail 2");
            roi_ra.push_back(287.2);
            roi_dec.push_back(6.6);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.15);

            roi_name.push_back("G40.5-0.5");
            roi_ra.push_back(286.786);
            roi_dec.push_back(6.498);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.2);

            roi_name.push_back("VERITAS (HAWC region)");
            roi_ra.push_back(287.05);
            roi_dec.push_back(6.39);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(1.5);

            roi_name.push_back("Hot spot north");
            roi_ra.push_back(286.8);
            roi_dec.push_back(7.1);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);

            roi_name.push_back("Hot spot west");
            roi_ra.push_back(286.2);
            roi_dec.push_back(6.4);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);

            roi_name.push_back("Hot spot east");
            roi_ra.push_back(288.1);
            roi_dec.push_back(6.4);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);

            double ring_center_ra = 286.91;
            double ring_center_dec = 6.32;
            roi_name.push_back("Ring 1");
            roi_ra.push_back(ring_center_ra);
            roi_dec.push_back(ring_center_dec);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.5);

            roi_name.push_back("Ring 2");
            roi_ra.push_back(ring_center_ra);
            roi_dec.push_back(ring_center_dec);
            roi_radius_inner.push_back(0.5);
            roi_radius_outer.push_back(1.0);

            roi_name.push_back("Ring 3");
            roi_ra.push_back(ring_center_ra);
            roi_dec.push_back(ring_center_dec);
            roi_radius_inner.push_back(1.0);
            roi_radius_outer.push_back(1.5);

        }
        else if (TString(target).Contains("SS433")) 
        {
            // Use RoI definition from here: https://veritas.sao.arizona.edu/wiki/index.php/SS433_Observations#EventDisplay_v4.80a_analysis_on_SS433_.28Payel.29
            
            roi_name.push_back("SS 433");
            roi_ra.push_back(287.9565);
            roi_dec.push_back(4.9827);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.15);

            roi_name.push_back("SS 433 e1");
            roi_ra.push_back(288.404);
            roi_dec.push_back(4.930);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.13);

            roi_name.push_back("SS 433 e1 max");
            roi_ra.push_back(288.207);
            roi_dec.push_back(4.951);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.13);

            roi_name.push_back("SS 433 e2");
            roi_ra.push_back(288.58);
            roi_dec.push_back(4.91);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.32);

            roi_name.push_back("SS 433 w1");
            roi_ra.push_back(287.42);
            roi_dec.push_back(5.04);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.14);
        }
        else if (TString(target).Contains("MAGIC_J1857")) 
        {
            roi_name.push_back("PSR J1856+0245");
            roi_ra.push_back(284.211666667);
            roi_dec.push_back(2.76394444444);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);

            roi_name.push_back("HESS J1858+020");
            roi_ra.push_back(284.583333333);
            roi_dec.push_back(2.09);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);
        }
        else if (TString(target).Contains("MGRO_J2031")) 
        {
            roi_name.push_back("MGRO J2031+41");
            roi_ra.push_back(mean_tele_point_ra);
            roi_dec.push_back(mean_tele_point_dec);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(1.0);
            roi_name.push_back("PSR J2032+4127");
            roi_ra.push_back(308.041666667);
            roi_dec.push_back(41.4594444444);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(1.0);
        }
        else if (TString(target).Contains("Crab")) 
        {
            if (TString(target).Contains("Offset"))
            {
                roi_name.push_back("Crab int. r=0.3 deg");
                roi_ra.push_back(mean_tele_point_ra);
                roi_dec.push_back(mean_tele_point_dec);
                roi_radius_inner.push_back(0.);
                roi_radius_outer.push_back(0.3);
                roi_name.push_back("Crab int. r=0.5 deg");
                roi_ra.push_back(mean_tele_point_ra);
                roi_dec.push_back(mean_tele_point_dec);
                roi_radius_inner.push_back(0.);
                roi_radius_outer.push_back(0.5);
                roi_name.push_back("Crab int. r=1.0 deg");
                roi_ra.push_back(mean_tele_point_ra);
                roi_dec.push_back(mean_tele_point_dec);
                roi_radius_inner.push_back(0.);
                roi_radius_outer.push_back(1.0);
                roi_name.push_back("Crab int. r=1.5 deg");
                roi_ra.push_back(mean_tele_point_ra);
                roi_dec.push_back(mean_tele_point_dec);
                roi_radius_inner.push_back(0.);
                roi_radius_outer.push_back(1.5);
            } 
            else
            {
                roi_name.push_back("Crab");
                roi_ra.push_back(mean_tele_point_ra);
                roi_dec.push_back(mean_tele_point_dec);
                roi_radius_inner.push_back(0.);
                roi_radius_outer.push_back(0.3);
            }
        }
        else if (TString(target).Contains("Mrk421")) 
        {
            roi_name.push_back("Mrk421");
            roi_ra.push_back(mean_tele_point_ra);
            roi_dec.push_back(mean_tele_point_dec);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.5);
        }
        else if (TString(target).Contains("Geminga")) 
        {
            roi_name.push_back("Geminga Pulsar");
            roi_ra.push_back(98.476);
            roi_dec.push_back(17.770);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(2.0);
        }
        else if (TString(target).Contains("Cygnus")) 
        {
            roi_name.push_back("Cygnus");
            roi_ra.push_back(mean_tele_point_ra);
            roi_dec.push_back(mean_tele_point_dec);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.5);
        }
        else if (TString(target).Contains("IC443")) 
        {
            roi_name.push_back("IC 443");
            roi_ra.push_back(mean_tele_point_ra);
            roi_dec.push_back(mean_tele_point_dec);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.5);
        }
        else if (TString(target).Contains("1ES1218")) 
        {
            roi_name.push_back("W Comae");
            roi_ra.push_back(185.382);
            roi_dec.push_back(28.233);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);

            roi_name.push_back("1ES 1218+304");
            roi_ra.push_back(185.360);
            roi_dec.push_back(30.191);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);

            roi_name.push_back("1ES 1215+303");
            roi_ra.push_back(184.616);
            roi_dec.push_back(30.130);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);
        }
        else if (TString(target).Contains("WComae")) 
        {
            roi_name.push_back("W Comae");
            roi_ra.push_back(185.382);
            roi_dec.push_back(28.233);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);

            roi_name.push_back("1ES 1218+304");
            roi_ra.push_back(185.360);
            roi_dec.push_back(30.191);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);

            roi_name.push_back("1ES 1215+303");
            roi_ra.push_back(184.616);
            roi_dec.push_back(30.130);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);
        }
        else if (TString(target).Contains("Boomerang")) 
        {
            roi_name.push_back("Boomerang");
            roi_ra.push_back(337.001);
            roi_dec.push_back(60.978);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.5);

            roi_name.push_back("Unknown");
            roi_ra.push_back(340.278);
            roi_dec.push_back(61.112);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);
        }
        else if (TString(target).Contains("GammaCygni")) 
        {
            roi_name.push_back("PSR J2021+4026");
            roi_ra.push_back(305.377066667);
            roi_dec.push_back(40.4481944444);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.24);

            roi_name.push_back("VER J2019+407");
            roi_ra.push_back(304.95);
            roi_dec.push_back(40.9);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.24);
        }
        else if (TString(target).Contains("LHAASO_J1929"))
        {
            roi_name.push_back("J1929+1745");
            roi_ra.push_back(mean_tele_point_ra);
            roi_dec.push_back(mean_tele_point_dec);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.6);
        }
        else if (TString(target).Contains("LHAASO_J1843"))
        {
            roi_name.push_back("J1843-0338");
            roi_ra.push_back(mean_tele_point_ra);
            roi_dec.push_back(mean_tele_point_dec);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.6);
        }
        else if (TString(target).Contains("LHAASO_J1956"))
        {
            roi_name.push_back("J1956+2845");
            roi_ra.push_back(mean_tele_point_ra);
            roi_dec.push_back(mean_tele_point_dec);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.6);
        }
        else
        {
            roi_name.push_back("Center");
            roi_ra.push_back(mean_tele_point_ra);
            roi_dec.push_back(mean_tele_point_dec);
            roi_radius_inner.push_back(0.);
            roi_radius_outer.push_back(0.3);
        }
    }


    //for (int star=0;star<FaintStars_Data.size();star++)
    //{
    //    roi_name.push_back("b-mag "+ std::to_string(FaintStars_Data.at(star).at(3)));
    //    roi_ra.push_back(FaintStars_Data.at(star).at(0));
    //    roi_dec.push_back(FaintStars_Data.at(star).at(1));
    //    roi_radius_inner.push_back(0.);
    //    roi_radius_outer.push_back(bright_star_radius_cut);
    //}

    shower_count = 0.;
    mean_elev = 0.;
    mean_azim = 0.;
    mean_nsb = 0.;
    TH1D Hist_Elev = TH1D("Hist_Elev","",N_elev_bins,elev_bins);
    TH1D Hist_MJD = TH1D("Hist_MJD","",N_MJD_bins,MJD_bins);
    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    TH1D Hist_ErecS_fine = TH1D("Hist_ErecS_fine","",N_energy_fine_bins,energy_fine_bins);
    TH1D Hist_EffArea = TH1D("Hist_EffArea","",N_energy_fine_bins,energy_fine_bins);
    TH1D Hist_Dark_NSB = TH1D("Hist_Dark_NSB","",20,4,14);
    TH2D Hist_Dark_ShowerDirection = TH2D("Hist_Dark_ShowerDirection","",180,0,360,90,0,90);
    TH2D Hist_Dark_ElevNSB = TH2D("Hist_Dark_ElevNSB","",20,0,10,18,0,90);
    TH2D Hist_Dark_ElevAzim = TH2D("Hist_Dark_ElevAzim","",18,0,360,18,0,90);
    TH1D Hist_Data_NSB = TH1D("Hist_Data_NSB","",20,4,14);
    TH2D Hist_Data_ShowerDirection = TH2D("Hist_Data_ShowerDirection","",180,0,360,90,0,90);
    TH2D Hist_Data_ElevNSB = TH2D("Hist_Data_ElevNSB","",20,0,10,18,0,90);
    TH2D Hist_Data_ElevAzim = TH2D("Hist_Data_ElevAzim","",18,0,360,18,0,90);
    TH2D Hist_Data_Skymap = TH2D("Hist_Data_Skymap","",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size);
    TH2D Hist_Data_Elev_Skymap = TH2D("Hist_Data_Elev_Skymap","",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size);
    TH2D Hist_Data_Azim_Skymap = TH2D("Hist_Data_Azim_Skymap","",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size);
    TH2D Hist_Data_NSB_Skymap = TH2D("Hist_Data_NSB_Skymap","",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size);

    vector<TH1D> Hist_OnData_Incl_CR_Zenith;
    vector<TH1D> Hist_OnDark_Incl_CR_Zenith;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));

        MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-(-1.*MSCW_cut_blind))+MSCW_cut_blind;
        MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-(-1.*MSCL_cut_blind))+MSCL_cut_blind;
        MSCW_plot_lower = -0.5*gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-(-1.*MSCW_cut_blind))-MSCW_cut_blind;
        MSCL_plot_lower = -0.5*gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-(-1.*MSCL_cut_blind))-MSCL_cut_blind;
        N_bins_for_deconv = N_bins_for_deconv_func_E[e];

        Hist_OnData_Incl_CR_Zenith.push_back(TH1D("Hist_OnData_Incl_CR_Zenith_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
        Hist_OnDark_Incl_CR_Zenith.push_back(TH1D("Hist_OnDark_Incl_CR_Zenith_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
    }

    vector<TH1D> Hist_OnData_Correct_R2off;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        int R2off_bins = 18;
        if (energy_bins[e]>=650.) R2off_bins = 9;
        if (energy_bins[e]>=4000.) R2off_bins = 1;
        Hist_OnData_Correct_R2off.push_back(TH1D("Hist_Stage1_OnData_Correct_R2off_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",R2off_bins,0,9));
    }


    vector<TH2D> Hist_OnData_MSCLW;
    vector<TH2D> Hist_OnData_Point_MSCLW;
    vector<TH2D> Hist_OnData_Ring_MSCLW;
    vector<TH1D> Hist_OnData_SR_Energy;
    vector<TH1D> Hist_OnData_CR_Energy;
    vector<TH1D> Hist_OnDark_SR_Energy;
    vector<TH1D> Hist_OnData_SR_Skymap_Theta2;
    vector<TH1D> Hist_OnData_CR_Skymap_Theta2;
    vector<TH1D> Hist_NormSyst_Skymap_Theta2;
    vector<TH1D> Hist_ShapeSyst_Skymap_Theta2;
    vector<TH1D> Hist_OnData_SR_Yoff;
    vector<TH1D> Hist_OnData_SR_Xoff;
    vector<TH2D> Hist_OnData_SR_XYoff;
    vector<TH1D> Hist_OnData_CR_Yoff;
    vector<TH1D> Hist_OnData_CR_Xoff;
    vector<TH1D> Hist_OnData_SR_R2off;
    vector<TH1D> Hist_OnDark_SR_R2off;
    vector<TH1D> Hist_OnData_CR_R2off;
    vector<TH1D> Hist_OnData_ISR_R2off;
    vector<TH2D> Hist_OnData_CR_XYoff;
    vector<TH2D> Hist_OnDark_SR_XYoff;
    vector<TH1D> Hist_OnData_CR_Yoff_Raw;
    vector<TH2D> Hist_Photon_Exp_Skymap;
    vector<TH2D> Hist_Photon_Raw_Skymap;
    vector<TH2D> Hist_OnData_EffArea_Skymap;
    vector<TH2D> Hist_OnData_ISR_Skymap;
    vector<TH2D> Hist_OnData_SR_Skymap;
    vector<TH2D> Hist_OnData_CR_Skymap;
    vector<TH2D> Hist_NormSyst_Skymap;
    vector<vector<TH2D>> Hist_ShapeSyst_Skymap;
    vector<TH2D> Hist_OnData_Expo_Skymap;
    vector<TH2D> Hist_OnDark_SR_Skymap;
    vector<TH2D> Hist_OnData_SR_Skymap_Galactic;
    vector<TH2D> Hist_OnData_CR_Skymap_Galactic;
    vector<TH1D> Hist_OnData_SR_Height;
    vector<TH1D> Hist_OnData_CR_Height;
    vector<TH1D> Hist_OnData_SR_Depth;
    vector<TH1D> Hist_OnData_CR_Depth;
    vector<TH1D> Hist_OnData_SR_Zenith;
    vector<TH1D> Hist_OnData_CR_Zenith;
    vector<TH1D> Hist_OnData_SR_Rcore;
    vector<TH1D> Hist_OnData_CR_Rcore;
    vector<TH1D> Hist_OnData_SR_Energy_CamCenter;
    vector<TH1D> Hist_OnData_CR_Energy_CamCenter;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));

        int XYoff_bins = 36;

        MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-(-1.*MSCW_cut_blind))+MSCW_cut_blind;
        MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-(-1.*MSCL_cut_blind))+MSCL_cut_blind;
        MSCW_plot_lower = -0.5*gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-(-1.*MSCW_cut_blind))-MSCW_cut_blind;
        MSCL_plot_lower = -0.5*gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-(-1.*MSCL_cut_blind))-MSCL_cut_blind;
        N_bins_for_deconv = N_bins_for_deconv_func_E[e];

        Hist_OnData_MSCLW.push_back(TH2D("Hist_Stage1_OnData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnData_Point_MSCLW.push_back(TH2D("Hist_Stage1_OnData_Point_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnData_Ring_MSCLW.push_back(TH2D("Hist_Stage1_OnData_Ring_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnData_SR_Skymap_Theta2.push_back(TH1D("Hist_Stage1_OnData_SR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_CR_Skymap_Theta2.push_back(TH1D("Hist_Stage1_OnData_CR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_NormSyst_Skymap_Theta2.push_back(TH1D("Hist_Stage1_NormSyst_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_ShapeSyst_Skymap_Theta2.push_back(TH1D("Hist_Stage1_ShapeSyst_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_OnData_SR_Yoff.push_back(TH1D("Hist_Stage1_OnData_SR_Yoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",30,-3,3));
        Hist_OnData_SR_Xoff.push_back(TH1D("Hist_Stage1_OnData_SR_Xoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",30,-3,3));
        Hist_OnData_SR_XYoff.push_back(TH2D("Hist_Stage1_OnData_SR_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",XYoff_bins,-3,3,XYoff_bins,-3,3));
        Hist_OnData_CR_Yoff.push_back(TH1D("Hist_Stage1_OnData_CR_Yoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",30,-3,3));
        Hist_OnData_CR_Xoff.push_back(TH1D("Hist_Stage1_OnData_CR_Xoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",30,-3,3));
        Hist_OnData_SR_R2off.push_back(TH1D("Hist_Stage1_OnData_SR_R2off_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",36,0,9));
        Hist_OnDark_SR_R2off.push_back(TH1D("Hist_Stage1_OnDark_SR_R2off_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",36,0,9));
        Hist_OnData_CR_R2off.push_back(TH1D("Hist_Stage1_OnData_CR_R2off_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",36,0,9));
        Hist_OnData_ISR_R2off.push_back(TH1D("Hist_Stage1_OnData_ISR_R2off_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",9,0,9));
        Hist_OnData_CR_XYoff.push_back(TH2D("Hist_Stage1_OnData_CR_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",XYoff_bins,-3,3,XYoff_bins,-3,3));
        Hist_OnDark_SR_XYoff.push_back(TH2D("Hist_Stage1_OnDark_SR_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",XYoff_bins,-3,3,XYoff_bins,-3,3));
        Hist_OnData_CR_Yoff_Raw.push_back(TH1D("Hist_Stage1_OnData_CR_Yoff_Raw_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",30,-3,3));
        Hist_Photon_Exp_Skymap.push_back(TH2D("Hist_Stage1_Photon_Exp_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_Photon_Raw_Skymap.push_back(TH2D("Hist_Stage1_Photon_Raw_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_EffArea_Skymap.push_back(TH2D("Hist_Stage1_OnData_EffArea_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_ISR_Skymap.push_back(TH2D("Hist_Stage1_OnData_ISR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_SR_Skymap.push_back(TH2D("Hist_Stage1_OnData_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_CR_Skymap.push_back(TH2D("Hist_Stage1_OnData_CR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_NormSyst_Skymap.push_back(TH2D("Hist_Stage1_NormSyst_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));

        vector<TH2D> Hist_ShapeSyst_Skymap_ThisBin;
        for (int xybin=0;xybin<N_integration_radii;xybin++) 
        {
            char xybin_tag[50];
            sprintf(xybin_tag,"Bin%i",xybin);
            Hist_ShapeSyst_Skymap_ThisBin.push_back(TH2D("Hist_Stage1_ShapeSyst_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up)+TString("_")+TString(xybin_tag),"",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        }
        Hist_ShapeSyst_Skymap.push_back(Hist_ShapeSyst_Skymap_ThisBin);

        Hist_OnData_Expo_Skymap.push_back(TH2D("Hist_Stage1_OnData_Expo_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins/2,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins/2,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnDark_SR_Skymap.push_back(TH2D("Hist_Stage1_OnDark_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins,mean_tele_point_ra-Skymap_size,mean_tele_point_ra+Skymap_size,Skymap_nbins,mean_tele_point_dec-Skymap_size,mean_tele_point_dec+Skymap_size));
        Hist_OnData_SR_Energy.push_back(TH1D("Hist_Stage1_OnData_SR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_CR_Energy.push_back(TH1D("Hist_Stage1_OnData_CR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnDark_SR_Energy.push_back(TH1D("Hist_Stage1_OnDark_SR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_SR_Height.push_back(TH1D("Hist_Stage1_OnData_SR_Height_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",25,0,25.));
        Hist_OnData_CR_Height.push_back(TH1D("Hist_Stage1_OnData_CR_Height_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",25,0,25.));
        Hist_OnData_SR_Depth.push_back(TH1D("Hist_Stage1_OnData_SR_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",25,0,25.));
        Hist_OnData_CR_Depth.push_back(TH1D("Hist_Stage1_OnData_CR_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",25,0,25.));
        Hist_OnData_SR_Zenith.push_back(TH1D("Hist_Stage1_OnData_SR_Zenith_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
        Hist_OnData_CR_Zenith.push_back(TH1D("Hist_Stage1_OnData_CR_Zenith_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",45,0,90));
        Hist_OnData_SR_Rcore.push_back(TH1D("Hist_Stage1_OnData_SR_Rcore_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",25,0,500.));
        Hist_OnData_CR_Rcore.push_back(TH1D("Hist_Stage1_OnData_CR_Rcore_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",25,0,500.));
        Hist_OnData_SR_Energy_CamCenter.push_back(TH1D("Hist_Stage1_OnData_SR_Energy_CamCenter_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_OnData_CR_Energy_CamCenter.push_back(TH1D("Hist_Stage1_OnData_CR_Energy_CamCenter_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        pair<double,double> tele_point_l_b = ConvertRaDecToGalactic(mean_tele_point_ra, mean_tele_point_dec);
        mean_tele_point_l = tele_point_l_b.first;
        mean_tele_point_b = tele_point_l_b.second;
        Hist_OnData_SR_Skymap_Galactic.push_back(TH2D("Hist_Stage1_OnData_SR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,Skymap_nbins,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
        Hist_OnData_CR_Skymap_Galactic.push_back(TH2D("Hist_Stage1_OnData_CR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins,tele_point_l_b.first-Skymap_size,tele_point_l_b.first+Skymap_size,Skymap_nbins,tele_point_l_b.second-Skymap_size,tele_point_l_b.second+Skymap_size));
    }


    vector<vector<TH2D>> Hist_OnDark_MSCLW;
    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
    {
        char sample_tag[50];
        sprintf(sample_tag, "%i", nth_sample);
        vector<TH2D> Hist_OnDark_OneSample_MSCLW;
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-(-1.*MSCW_cut_blind))+MSCW_cut_blind;
            MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-(-1.*MSCL_cut_blind))+MSCL_cut_blind;
            MSCW_plot_lower = -0.5*gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-(-1.*MSCW_cut_blind))-MSCW_cut_blind;
            MSCL_plot_lower = -0.5*gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-(-1.*MSCL_cut_blind))-MSCL_cut_blind;
            N_bins_for_deconv = N_bins_for_deconv_func_E[e];

            Hist_OnDark_OneSample_MSCLW.push_back(TH2D("Hist_Stage1_OnDark_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        }
        Hist_OnDark_MSCLW.push_back(Hist_OnDark_OneSample_MSCLW);
    }

    vector<TH2D> Hist_SRDark_XYoff;
    vector<TH2D> Hist_CRDark_XYoff;
    vector<TH2D> Hist_SRCRDarkRatio_XYoff;
    vector<TH2D> Hist_SRCRDarkRatio_XYoff_Smooth;
    vector<TH2D> Hist_SRDark_RaDec;
    vector<TH2D> Hist_CRDark_RaDec;
    vector<TH2D> Hist_SRCRDarkRatio_RaDec;
    vector<TH2D> Hist_SRCRDarkRatio_RaDec_Smooth;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        int Xoff_bins = 8;
        int Yoff_bins = 8;
        if (e>2)
        {
            Xoff_bins = 4;
            Yoff_bins = 4;
        }
        Hist_SRDark_XYoff.push_back(TH2D("Hist_SRDark_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Xoff_bins,-2,2,Yoff_bins,-2,2));
        Hist_CRDark_XYoff.push_back(TH2D("Hist_CRDark_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Xoff_bins,-2,2,Yoff_bins,-2,2));
        Hist_SRCRDarkRatio_XYoff.push_back(TH2D("Hist_SRCRDarkRatio_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Xoff_bins,-2,2,Yoff_bins,-2,2));
        Hist_SRCRDarkRatio_XYoff_Smooth.push_back(TH2D("Hist_SRCRDarkRatio_XYoff_Smooth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",32,-2,2,32,-2,2));
        int RaDec_bins = 8;
        if (e>=2)
        {
            RaDec_bins = 4;
        }
        Hist_SRDark_RaDec.push_back(TH2D("Hist_SRDark_RaDec_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",2*RaDec_bins,mean_tele_point_ra-2*Skymap_size,mean_tele_point_ra+2*Skymap_size,2*RaDec_bins,mean_tele_point_dec-2*Skymap_size,mean_tele_point_dec+2*Skymap_size));
        Hist_CRDark_RaDec.push_back(TH2D("Hist_CRDark_RaDec_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",2*RaDec_bins,mean_tele_point_ra-2*Skymap_size,mean_tele_point_ra+2*Skymap_size,2*RaDec_bins,mean_tele_point_dec-2*Skymap_size,mean_tele_point_dec+2*Skymap_size));
        Hist_SRCRDarkRatio_RaDec.push_back(TH2D("Hist_SRCRDarkRatio_RaDec_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",2*RaDec_bins,mean_tele_point_ra-2*Skymap_size,mean_tele_point_ra+2*Skymap_size,2*RaDec_bins,mean_tele_point_dec-2*Skymap_size,mean_tele_point_dec+2*Skymap_size));
        Hist_SRCRDarkRatio_RaDec_Smooth.push_back(TH2D("Hist_SRCRDarkRatio_RaDec_Smooth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",2*Skymap_nbins,mean_tele_point_ra-2*Skymap_size,mean_tele_point_ra+2*Skymap_size,2*Skymap_nbins,mean_tele_point_dec-2*Skymap_size,mean_tele_point_dec+2*Skymap_size));
    }
    vector<TH1D> Hist_SRDark_R2off;
    vector<TH1D> Hist_CRDark_R2off;
    vector<TH1D> Hist_SRDark_Energy;
    vector<TH1D> Hist_CRDark_Energy;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        Hist_SRDark_R2off.push_back(TH1D("Hist_SRDark_R2off_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",18,0,9));
        Hist_CRDark_R2off.push_back(TH1D("Hist_CRDark_R2off_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",18,0,9));
        Hist_SRDark_Energy.push_back(TH1D("Hist_SRDark_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_CRDark_Energy.push_back(TH1D("Hist_CRDark_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
    }


    vector<vector<TH1D>> Hist_OnData_SR_RoI_Energy;
    vector<vector<TH1D>> Hist_OnData_CR_RoI_Energy;
    vector<vector<TH1D>> Hist_NormSyst_RoI_Energy;
    vector<vector<TH1D>> Hist_ShapeSyst_RoI_Energy;
    vector<vector<TH1D>> Hist_OnData_SR_Skymap_RoI_Theta2;
    vector<vector<TH1D>> Hist_OnData_CR_Skymap_RoI_Theta2;
    vector<vector<TH1D>> Hist_NormSyst_Skymap_RoI_Theta2;
    vector<vector<TH1D>> Hist_ShapeSyst_Skymap_RoI_Theta2;
    vector<vector<TH1D>> Hist_OnData_SR_Skymap_RoI_X;
    vector<vector<TH1D>> Hist_OnData_CR_Skymap_RoI_X;
    vector<vector<TH1D>> Hist_OnData_SR_Skymap_RoI_Y;
    vector<vector<TH1D>> Hist_OnData_CR_Skymap_RoI_Y;
    vector<vector<TH1D>> Hist_OnData_SR_RoI_MJD;
    vector<vector<TH1D>> Hist_OnData_CR_RoI_MJD;
    for (int nth_roi=0;nth_roi<roi_ra.size();nth_roi++)
    {
        char roi_tag[50];
        sprintf(roi_tag, "%i", nth_roi);
        vector<TH1D> Hist_OnData_OneRoI_SR_RoI_Energy;
        vector<TH1D> Hist_OnData_OneRoI_CR_RoI_Energy;
        vector<TH1D> Hist_OneNormSyst_RoI_Energy;
        vector<TH1D> Hist_OneShapeSyst_RoI_Energy;
        vector<TH1D> Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2;
        vector<TH1D> Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2;
        vector<TH1D> Hist_OneNormSyst_Skymap_RoI_Theta2;
        vector<TH1D> Hist_OneShapeSyst_Skymap_RoI_Theta2;
        vector<TH1D> Hist_OnData_OneRoI_SR_Skymap_RoI_X;
        vector<TH1D> Hist_OnData_OneRoI_CR_Skymap_RoI_X;
        vector<TH1D> Hist_OnData_OneRoI_SR_Skymap_RoI_Y;
        vector<TH1D> Hist_OnData_OneRoI_CR_Skymap_RoI_Y;
        vector<TH1D> Hist_OnData_OneRoI_SR_RoI_MJD;
        vector<TH1D> Hist_OnData_OneRoI_CR_RoI_MJD;
        double roi_range = 2.*roi_radius_outer.at(nth_roi);
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));

            MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-(-1.*MSCW_cut_blind))+MSCW_cut_blind;
            MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-(-1.*MSCL_cut_blind))+MSCL_cut_blind;
            MSCW_plot_lower = -0.5*gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-(-1.*MSCW_cut_blind))-MSCW_cut_blind;
            MSCL_plot_lower = -0.5*gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-(-1.*MSCL_cut_blind))-MSCL_cut_blind;
            N_bins_for_deconv = N_bins_for_deconv_func_E[e];

            Hist_OnData_OneRoI_SR_RoI_Energy.push_back(TH1D("Hist_Stage1_OnData_SR_RoI_Energy_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
            Hist_OnData_OneRoI_CR_RoI_Energy.push_back(TH1D("Hist_Stage1_OnData_CR_RoI_Energy_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
            Hist_OneNormSyst_RoI_Energy.push_back(TH1D("Hist_Stage1_NormSyst_RoI_Energy_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
            Hist_OneShapeSyst_RoI_Energy.push_back(TH1D("Hist_Stage1_ShapeSyst_RoI_Energy_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
            Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2.push_back(TH1D("Hist_Stage1_OnData_SR_Skymap_RoI_Theta2_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,1.0));
            Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2.push_back(TH1D("Hist_Stage1_OnData_CR_Skymap_RoI_Theta2_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,1.0));
            Hist_OneNormSyst_Skymap_RoI_Theta2.push_back(TH1D("Hist_Stage1_NormSyst_Skymap_RoI_Theta2_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,1.0));
            Hist_OneShapeSyst_Skymap_RoI_Theta2.push_back(TH1D("Hist_Stage1_ShapeSyst_Skymap_RoI_Theta2_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,1.0));
            Hist_OnData_OneRoI_SR_Skymap_RoI_X.push_back(TH1D("Hist_Stage1_OnData_SR_Skymap_RoI_X_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",60,-roi_range,roi_range));
            Hist_OnData_OneRoI_CR_Skymap_RoI_X.push_back(TH1D("Hist_Stage1_OnData_CR_Skymap_RoI_X_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",60,-roi_range,roi_range));
            Hist_OnData_OneRoI_SR_Skymap_RoI_Y.push_back(TH1D("Hist_Stage1_OnData_SR_Skymap_RoI_Y_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",60,-roi_range,roi_range));
            Hist_OnData_OneRoI_CR_Skymap_RoI_Y.push_back(TH1D("Hist_Stage1_OnData_CR_Skymap_RoI_Y_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",60,-roi_range,roi_range));
            Hist_OnData_OneRoI_SR_RoI_MJD.push_back(TH1D("Hist_Stage1_OnData_SR_RoI_MJD_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",800,56200-4000,56200+4000));
            Hist_OnData_OneRoI_CR_RoI_MJD.push_back(TH1D("Hist_Stage1_OnData_CR_RoI_MJD_V"+TString(roi_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",800,56200-4000,56200+4000));
        }
        Hist_OnData_SR_RoI_Energy.push_back(Hist_OnData_OneRoI_SR_RoI_Energy);
        Hist_OnData_CR_RoI_Energy.push_back(Hist_OnData_OneRoI_CR_RoI_Energy);
        Hist_NormSyst_RoI_Energy.push_back(Hist_OneNormSyst_RoI_Energy);
        Hist_ShapeSyst_RoI_Energy.push_back(Hist_OneShapeSyst_RoI_Energy);
        Hist_OnData_SR_Skymap_RoI_Theta2.push_back(Hist_OnData_OneRoI_SR_Skymap_RoI_Theta2);
        Hist_OnData_CR_Skymap_RoI_Theta2.push_back(Hist_OnData_OneRoI_CR_Skymap_RoI_Theta2);
        Hist_NormSyst_Skymap_RoI_Theta2.push_back(Hist_OneNormSyst_Skymap_RoI_Theta2);
        Hist_ShapeSyst_Skymap_RoI_Theta2.push_back(Hist_OneShapeSyst_Skymap_RoI_Theta2);
        Hist_OnData_SR_Skymap_RoI_X.push_back(Hist_OnData_OneRoI_SR_Skymap_RoI_X);
        Hist_OnData_CR_Skymap_RoI_X.push_back(Hist_OnData_OneRoI_CR_Skymap_RoI_X);
        Hist_OnData_SR_Skymap_RoI_Y.push_back(Hist_OnData_OneRoI_SR_Skymap_RoI_Y);
        Hist_OnData_CR_Skymap_RoI_Y.push_back(Hist_OnData_OneRoI_CR_Skymap_RoI_Y);
        Hist_OnData_SR_RoI_MJD.push_back(Hist_OnData_OneRoI_SR_RoI_MJD);
        Hist_OnData_CR_RoI_MJD.push_back(Hist_OnData_OneRoI_CR_RoI_MJD);
    }
    
    vector<TH1D> Hist_Source_Theta2;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        Hist_Source_Theta2.push_back(TH1D("Hist_Source_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",5,0.,5.*source_theta2_cut));
    }


    std::cout << "Prepare dark run samples..." << std::endl;
    std::cout << "Reading file " << TString(SMI_OUTPUT)+"/Netflix_RunList_"+TString(target)+"_"+TString(output_file_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+TString(mjd_cut_tag)+"_"+ONOFF_tag+group_tag+".root" << std::endl;

    TFile InputRunListFile(TString(SMI_OUTPUT)+"/Netflix_RunList_"+TString(target)+"_"+TString(output_file_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+TString(mjd_cut_tag)+"_"+ONOFF_tag+group_tag+".root");
    TTree* RunListTree_ptr = nullptr;
    RunListTree_ptr = (TTree*) InputRunListFile.Get("RunListTree_subgroup");
    RunListTree_ptr->SetBranchStatus("*", 0);
    RunListTree_ptr->SetBranchStatus("ON_runnumber", 1);
    RunListTree_ptr->SetBranchStatus("ON_pointing_RA", 1);
    RunListTree_ptr->SetBranchStatus("ON_pointing_Dec", 1);
    RunListTree_ptr->SetBranchStatus("ON_exposure_hour", 1);
    RunListTree_ptr->SetBranchStatus("OFF_runnumber", 1);
    RunListTree_ptr->SetBranchStatus("OFF_exposure_hour", 1);
    int ON_runnumber_ptr;
    RunListTree_ptr->SetBranchAddress("ON_runnumber",&ON_runnumber_ptr);
    double ON_exposure_hour_ptr;
    RunListTree_ptr->SetBranchAddress("ON_exposure_hour",&ON_exposure_hour_ptr);
    double ON_pointing_RA_ptr;
    RunListTree_ptr->SetBranchAddress("ON_pointing_RA",&ON_pointing_RA_ptr);
    double ON_pointing_Dec_ptr;
    RunListTree_ptr->SetBranchAddress("ON_pointing_Dec",&ON_pointing_Dec_ptr);
    vector<int>* OFF_runnumber_ptr = new std::vector<int>(10);
    RunListTree_ptr->SetBranchAddress("OFF_runnumber",&OFF_runnumber_ptr);
    vector<double>* OFF_exposure_hour_ptr = new std::vector<double>(10);
    RunListTree_ptr->SetBranchAddress("OFF_exposure_hour",&OFF_exposure_hour_ptr);

    vector<pair<string,int>> Data_runlist;
    vector<double> Data_exposure_hour;
    vector<pair<double,double>> Data_runlist_RaDec;
    vector<vector<vector<pair<string,int>>>> Dark_runlist;
    vector<vector<vector<double>>> Dark_exposure_hour;
    for (int on_run=0;on_run<RunListTree_ptr->GetEntries();on_run++)
    {
        RunListTree_ptr->GetEntry(on_run);
        if (OFF_runnumber_ptr->size()==0) continue;
        if (OFF_runnumber_ptr->at(0)==0) continue;
        std::cout << "ON run " << ON_runnumber_ptr << ", number of OFF runs = " << OFF_runnumber_ptr->size() << std::endl;
        pair<string,int> on_run_pair = std::make_pair("",ON_runnumber_ptr);
        Data_runlist.push_back(on_run_pair);
        pair<double,double> on_run_radec = std::make_pair(ON_pointing_RA_ptr,ON_pointing_Dec_ptr);
        std::cout << "ON run RA = " << on_run_radec.first << ", Dec = " << on_run_radec.second << std::endl;
        Data_runlist_RaDec.push_back(on_run_radec);
        Data_exposure_hour.push_back(ON_exposure_hour_ptr);
        vector<vector<pair<string,int>>> off_run_vtr;
        vector<pair<string,int>> off_run_vtr2;
        vector<vector<double>> off_expo_vtr;
        vector<double> off_expo_vtr2;
        for (int off_run=0;off_run<OFF_runnumber_ptr->size();off_run++)
        {
            pair<string,int> off_run_pair = std::make_pair("",OFF_runnumber_ptr->at(off_run));
            off_run_vtr2.push_back(off_run_pair);
            off_expo_vtr2.push_back(OFF_exposure_hour_ptr->at(off_run));
        }
        off_run_vtr.push_back(off_run_vtr2);
        off_expo_vtr.push_back(off_expo_vtr2);
        Dark_runlist.push_back(off_run_vtr);
        Dark_exposure_hour.push_back(off_expo_vtr);
    }

    for (int run=0;run<Data_runlist.size();run++)
    {

        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        string filename;
        filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

        //std::cout << "Get telescope pointing RA and Dec..." << std::endl;
        pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
        //if (!TString(target).Contains("Proton")) tele_point_ra_dec = GetRunRaDec(filename,int(Data_runlist[run].second));
        tele_point_ra_dec = Data_runlist_RaDec[run];
        run_tele_point_ra = tele_point_ra_dec.first;
        run_tele_point_dec = tele_point_ra_dec.second;
        pair<double,double> tele_point_ra_dec_imposter = std::make_pair(0,0);
        if (doImposter) 
        {
            tele_point_ra_dec_imposter = GetRunRaDec(filename,int(Data_runlist[run].second));
        }

        for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
        {
            for (int off_run=0;off_run<Dark_runlist.at(run).at(nth_sample).size();off_run++)
            {

                std::cout << "Prepare off run ..." << int(Dark_runlist.at(run).at(nth_sample)[off_run].second) << " for on run " << run_number << std::endl;
                if (Dark_runlist.at(run).at(nth_sample)[off_run].second==0) continue;
                char run_number[50];
                char Dark_observation[50];
                sprintf(run_number, "%i", int(Dark_runlist.at(run).at(nth_sample)[off_run].second));
                sprintf(Dark_observation, "%s", Dark_runlist.at(run).at(nth_sample)[off_run].first.c_str());
                string filename_dark;
                filename_dark = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

                TFile*  dark_input_file = TFile::Open(filename_dark.c_str());
                TH1* i_hEffAreaP = ( TH1* )getEffAreaHistogram(dark_input_file,Dark_runlist.at(run).at(nth_sample)[off_run].second, 0.5);
                TString root_file = "run_"+string(run_number)+"/stereo/data_on";
                TTree* Dark_tree = (TTree*) dark_input_file->Get(root_file);
                if (!Dark_tree)
                {
                    std::cout << "TTree does not exist: " << root_file << std::endl;
                    continue;
                }
                SetEventDisplayTreeBranch(Dark_tree);
                vector<pair<double,double>> timecut_thisrun = GetRunTimecuts(Dark_runlist.at(run).at(nth_sample)[off_run].second);

                pair<double,double> dark_tele_point_ra_dec = std::make_pair(0,0);
                dark_tele_point_ra_dec = GetRunRaDec(filename_dark,int(Dark_runlist.at(run).at(nth_sample)[off_run].second));
                double mean_dark_tele_point_ra = 0.;
                double mean_dark_tele_point_dec = 0.;
                pair<double,double> dark_source_ra_dec = GetSourceRaDec(TString(Dark_observation));
                mean_dark_tele_point_ra = source_ra_dec.first;
                mean_dark_tele_point_dec = source_ra_dec.second;

                Dark_tree->GetEntry(0);
                double time_0 = Time;
                Dark_tree->GetEntry(Dark_tree->GetEntries()-1);
                double time_1 = Time;

                double NSB_thisrun = GetRunPedestalVar(int(Dark_runlist.at(run).at(nth_sample)[off_run].second));
                Hist_Dark_NSB.Fill(NSB_thisrun);
                double tele_elev_off = GetRunElevAzim(filename_dark,int(Dark_runlist.at(run).at(nth_sample)[off_run].second)).first;
                double tele_azim_off = GetRunElevAzim(filename_dark,int(Dark_runlist.at(run).at(nth_sample)[off_run].second)).second;

                for (int e=0;e<N_energy_bins;e++) 
                {
                    Hist_Source_Theta2.at(e).Reset();
                }

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
                    if (doImposter)
                    {
                        ra_sky_imposter = tele_point_ra_dec_imposter.first+Xoff_derot;
                        dec_sky_imposter = tele_point_ra_dec_imposter.second+Yoff_derot;
                    }
                    //if (doRaster)
                    //{
                    //    double delta_phi = 2*M_PI*double(entry)/double(Dark_tree->GetEntries());
                    //    double delta_r = 1.0;
                    //    ra_sky += delta_r*cos(delta_phi);
                    //    dec_sky += delta_r*sin(delta_phi);
                    //}
                    ra_sky_dark = dark_tele_point_ra_dec.first+Xoff_derot;
                    dec_sky_dark = dark_tele_point_ra_dec.second+Yoff_derot;
                    // redefine theta2
                    theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
                    theta2_dark = pow(ra_sky_dark-mean_dark_tele_point_ra,2)+pow(dec_sky_dark-mean_dark_tele_point_dec,2);
                    pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
                    int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
                    int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
                    int elevation = Hist_Elev.FindBin(tele_elev_off)-1;
                    if (energy<0) continue;
                    if (energy>=N_energy_bins) continue;
                    if (elevation<0) continue;
                    if (elevation>=N_elev_bins) continue;
                    if (!SelectNImages()) continue;
                    if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
                    if (SizeSecondMax<SizeSecondMax_Cut) continue;
                    //if (EmissionHeight<6.) continue;
                    double shower_depth = GetShowerDepth(EmissionHeight,tele_elev_off);
                    //if (shower_depth>4.) continue;
                    if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
                    //if (pow(Xcore*Xcore+Ycore*Ycore,0.5)<100) continue;
                    //if (R2off>4.) continue;
                    MSCW = RescaleMSCW(MSCW, R2off, MSCW_rescale[energy]);
                    MSCL = RescaleMSCW(MSCL, R2off, MSCL_rescale[energy]);
                    double run_weight = Data_exposure_hour[run]/Dark_exposure_hour.at(run).at(nth_sample)[off_run];
                    double weight = run_weight;
                    if (DarkFoV())
                    {
                        if (SignalSelectionTheta2())
                        {
                            Hist_Source_Theta2.at(energy).Fill(theta2_dark);
                        }
                    }
                }

                vector<double> source_weight;
                for (int e=0;e<N_energy_bins;e++) 
                {
                    int nbins = Hist_Source_Theta2.at(e).GetNbinsX();
                    double background_counts = 0.;
                    for (int bin=2;bin<=nbins;bin++) // first bin is the signal region
                    {
                        background_counts += Hist_Source_Theta2.at(e).GetBinContent(bin);
                    }
                    double background_avg = background_counts/double(nbins-1);
                    if (Hist_Source_Theta2.at(e).GetBinContent(1)>0.)
                    {
                        source_weight.push_back(background_avg/Hist_Source_Theta2.at(e).GetBinContent(1));
                    }
                    else
                    {
                        source_weight.push_back(0.);
                    }
                }

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
                    if (doImposter)
                    {
                        ra_sky_imposter = tele_point_ra_dec_imposter.first+Xoff_derot;
                        dec_sky_imposter = tele_point_ra_dec_imposter.second+Yoff_derot;
                    }
                    //if (doRaster)
                    //{
                    //    double delta_phi = 2*M_PI*double(entry)/double(Dark_tree->GetEntries());
                    //    double delta_r = 1.0;
                    //    ra_sky += delta_r*cos(delta_phi);
                    //    dec_sky += delta_r*sin(delta_phi);
                    //}
                    ra_sky_dark = dark_tele_point_ra_dec.first+Xoff_derot;
                    dec_sky_dark = dark_tele_point_ra_dec.second+Yoff_derot;
                    // redefine theta2
                    theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
                    theta2_dark = pow(ra_sky_dark-mean_dark_tele_point_ra,2)+pow(dec_sky_dark-mean_dark_tele_point_dec,2);
                    pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
                    int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
                    int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
                    int elevation = Hist_Elev.FindBin(tele_elev_off)-1;
                    if (energy<0) continue;
                    if (energy>=N_energy_bins) continue;
                    if (elevation<0) continue;
                    if (elevation>=N_elev_bins) continue;
                    if (!SelectNImages()) continue;
                    if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
                    if (SizeSecondMax<SizeSecondMax_Cut) continue;
                    //if (EmissionHeight<6.) continue;
                    double shower_depth = GetShowerDepth(EmissionHeight,tele_elev_off);
                    //if (shower_depth>4.) continue;
                    if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
                    //if (pow(Xcore*Xcore+Ycore*Ycore,0.5)<100) continue;
                    //if (R2off>4.) continue;
                    MSCW = RescaleMSCW(MSCW, R2off, MSCW_rescale[energy]);
                    MSCL = RescaleMSCW(MSCL, R2off, MSCL_rescale[energy]);
                    Hist_Dark_ShowerDirection.Fill(Shower_Az,Shower_Ze);
                    Hist_Dark_ElevNSB.Fill(NSB_thisrun,tele_elev_off);
                    Hist_Dark_ElevAzim.Fill(NSB_thisrun,tele_azim_off);
                    double run_weight = Data_exposure_hour[run]/Dark_exposure_hour.at(run).at(nth_sample)[off_run];
                    double weight = run_weight;
                    //if (theta2_dark<source_theta2_cut && SignalSelectionTheta2()) weight = run_weight*source_weight.at(energy);

                    if (DarkFoV())
                    {
                        if (SignalSelectionTheta2())
                        {
                            Hist_SRDark_XYoff.at(energy).Fill(Xoff,Yoff,weight);
                            Hist_SRDark_RaDec.at(energy).Fill(ra_sky,dec_sky,weight);
                            Hist_OnDark_SR_XYoff.at(energy).Fill(Xoff,Yoff,weight);
                            Hist_OnDark_SR_R2off.at(energy).Fill(R2off,weight);
                            Hist_SRDark_Energy.at(energy).Fill(ErecS*1000.,weight);
                            Hist_SRDark_R2off.at(energy).Fill(R2off,weight);
                            if (FoV(false))
                            {
                                Hist_OnDark_SR_Skymap.at(energy).Fill(ra_sky,dec_sky,weight);
                                Hist_OnDark_SR_Energy.at(energy).Fill(ErecS*1000.,weight);
                            }
                        }
                        else if (ControlSelectionTheta2())
                        {
                            Hist_CRDark_XYoff.at(energy).Fill(Xoff,Yoff,weight);
                            Hist_CRDark_RaDec.at(energy).Fill(ra_sky,dec_sky,weight);
                            Hist_CRDark_Energy.at(energy).Fill(ErecS*1000.,weight);
                            Hist_CRDark_R2off.at(energy).Fill(R2off,weight);
                        }

                        if (theta2_dark>source_theta2_cut)
                        {
                            if (FoV(false))
                            {
                                Hist_OnDark_MSCLW.at(nth_sample).at(energy).Fill(MSCL,MSCW,weight);
                            }
                        }
                    }

                }
                dark_input_file->Close();
            }
        }
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        Hist_SRCRDarkRatio_RaDec.at(e).Reset();
        Hist_SRCRDarkRatio_RaDec.at(e).Add(&Hist_SRDark_RaDec.at(e));
        Hist_SRCRDarkRatio_RaDec.at(e).Divide(&Hist_CRDark_RaDec.at(e));
        for (int binx=1;binx<=Hist_SRCRDarkRatio_RaDec_Smooth.at(e).GetNbinsX();binx++)
        {
            for (int biny=1;biny<=Hist_SRCRDarkRatio_RaDec_Smooth.at(e).GetNbinsY();biny++)
            {
                double ra_sky = Hist_SRCRDarkRatio_RaDec_Smooth.at(e).GetXaxis()->GetBinCenter(binx);
                double dec_sky = Hist_SRCRDarkRatio_RaDec_Smooth.at(e).GetYaxis()->GetBinCenter(biny);
                int big_bin_ra = Hist_SRCRDarkRatio_RaDec.at(e).GetXaxis()->FindBin(ra_sky);
                int big_bin_dec = Hist_SRCRDarkRatio_RaDec.at(e).GetYaxis()->FindBin(dec_sky);
                double ratio_content = Hist_SRCRDarkRatio_RaDec.at(e).GetBinContent(big_bin_ra,big_bin_dec);
                Hist_SRCRDarkRatio_RaDec_Smooth.at(e).SetBinContent(binx,biny,ratio_content);
            }
        }
        double smooth_size = 0.5;
        Smooth2DMap(&Hist_SRCRDarkRatio_RaDec_Smooth.at(e), smooth_size);
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        Hist_SRCRDarkRatio_XYoff.at(e).Reset();
        Hist_SRCRDarkRatio_XYoff.at(e).Add(&Hist_SRDark_XYoff.at(e));
        Hist_SRCRDarkRatio_XYoff.at(e).Divide(&Hist_CRDark_XYoff.at(e));
        for (int binx=1;binx<=Hist_SRCRDarkRatio_XYoff_Smooth.at(e).GetNbinsX();binx++)
        {
            for (int biny=1;biny<=Hist_SRCRDarkRatio_XYoff_Smooth.at(e).GetNbinsY();biny++)
            {
                double x_off = Hist_SRCRDarkRatio_XYoff_Smooth.at(e).GetXaxis()->GetBinCenter(binx);
                double y_off = Hist_SRCRDarkRatio_XYoff_Smooth.at(e).GetYaxis()->GetBinCenter(biny);
                int big_bin_x_off = Hist_SRCRDarkRatio_XYoff.at(e).GetXaxis()->FindBin(x_off);
                int big_bin_y_off = Hist_SRCRDarkRatio_XYoff.at(e).GetYaxis()->FindBin(y_off);
                double ratio_content = Hist_SRCRDarkRatio_XYoff.at(e).GetBinContent(big_bin_x_off,big_bin_y_off);
                Hist_SRCRDarkRatio_XYoff_Smooth.at(e).SetBinContent(binx,biny,ratio_content);
            }
        }
        double smooth_size = 0.5;
        Smooth2DMap(&Hist_SRCRDarkRatio_XYoff_Smooth.at(e), smooth_size);
    }

    std::cout << "Build acceptance function from cosmic rays." << std::endl;
    vector<double> Data_runlist_exposure;
    for (int run=0;run<Data_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        string filename;
        filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

        std::cout << "Get telescope pointing RA and Dec for run " << run_number << std::endl;
        pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
        //if (!TString(target).Contains("Proton")) tele_point_ra_dec = GetRunRaDec(filename,int(Data_runlist[run].second));
        tele_point_ra_dec = Data_runlist_RaDec[run];
        run_tele_point_ra = tele_point_ra_dec.first;
        run_tele_point_dec = tele_point_ra_dec.second;
        pair<double,double> tele_point_ra_dec_imposter = std::make_pair(0,0);
        if (doImposter) 
        {
            tele_point_ra_dec_imposter = GetRunRaDec(filename,int(Data_runlist[run].second));
        }

        TFile*  input_file = TFile::Open(filename.c_str());
        TString root_file = "run_"+string(run_number)+"/stereo/data_on";
        TTree* Data_tree = (TTree*) input_file->Get(root_file);
        if (!Data_tree)
        {
            std::cout << "TTree does not exist: " << root_file << std::endl;
            continue;
        }
        SetEventDisplayTreeBranch(Data_tree);
        std::cout << "Get time cuts for run " << run_number << std::endl;
        vector<pair<double,double>> timecut_thisrun = GetRunTimecuts(Data_runlist[run].second);
        std::cout << "Get elev. and azim. for run " << run_number << std::endl;
        double tele_elev = GetRunElevAzim(filename,int(Data_runlist[run].second)).first;

        Data_tree->GetEntry(0);
        double time_0 = Time;
        Data_tree->GetEntry(Data_tree->GetEntries()-1);
        double time_1 = Time;
        exposure_hours += (time_1-time_0)/3600.;
        std::cout << "Get usable time for run " << run_number << std::endl;
        double exposure_thisrun = GetRunUsableTime(filename,Data_runlist[run].second)/3600.;
        exposure_hours_usable += exposure_thisrun;
        Data_runlist_exposure.push_back(exposure_thisrun);

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
            if (doImposter)
            {
                ra_sky_imposter = tele_point_ra_dec_imposter.first+Xoff_derot;
                dec_sky_imposter = tele_point_ra_dec_imposter.second+Yoff_derot;
            }
            //if (doRaster)
            //{
            //    double delta_phi = 2*M_PI*double(entry)/double(Data_tree->GetEntries());
            //    double delta_r = 1.0;
            //    ra_sky += delta_r*cos(delta_phi);
            //    dec_sky += delta_r*sin(delta_phi);
            //}
            // redefine theta2
            theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
            pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            int elevation = Hist_Elev.FindBin(tele_elev)-1;
            int year = Hist_MJD.FindBin(MJD)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            if (elevation<0) continue;
            if (elevation>=N_elev_bins) continue;
            if (!SelectNImages()) continue;
            if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
            if (SizeSecondMax<SizeSecondMax_Cut) continue;
            //if (EmissionHeight<6.) continue;
            double shower_depth = GetShowerDepth(EmissionHeight,tele_elev);
            //if (shower_depth>4.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            //if (pow(Xcore*Xcore+Ycore*Ycore,0.5)<100) continue;
            //if (TString(target).Contains("Crab") && theta2<0.3) continue;
            //if (TString(target).Contains("Mrk421") && theta2<0.3) continue;
            //if (R2off>4.) continue;
            MSCW = RescaleMSCW(MSCW, R2off, MSCW_rescale[energy]);
            MSCL = RescaleMSCW(MSCL, R2off, MSCL_rescale[energy]);
            
            int bin_energy = Hist_CRDark_Energy.at(energy).FindBin(ErecS*1000.);
            int big_bin_xoff = Hist_SRCRDarkRatio_XYoff_Smooth.at(energy).GetXaxis()->FindBin(Xoff);
            int big_bin_yoff = Hist_SRCRDarkRatio_XYoff_Smooth.at(energy).GetYaxis()->FindBin(Yoff);
            int big_bin_ra = Hist_SRCRDarkRatio_RaDec_Smooth.at(energy).GetXaxis()->FindBin(ra_sky);
            int big_bin_dec = Hist_SRCRDarkRatio_RaDec_Smooth.at(energy).GetYaxis()->FindBin(dec_sky);
            double data_cr_content_energy = Hist_SRDark_Energy.at(energy).GetBinContent(bin_energy);
            double dark_cr_content_energy = Hist_CRDark_Energy.at(energy).GetBinContent(bin_energy);
            double xyoff_weight = Hist_SRCRDarkRatio_XYoff_Smooth.at(energy).GetBinContent(big_bin_xoff,big_bin_yoff);
            double radec_weight = Hist_SRCRDarkRatio_RaDec_Smooth.at(energy).GetBinContent(big_bin_ra,big_bin_dec);
            double energy_weight = 1.;
            if (dark_cr_content_energy>0.)
            {
                energy_weight = data_cr_content_energy/dark_cr_content_energy;
            }

            if (ControlSelectionTheta2())
            {
                if (R2off<camera_theta2_cut_upper)
                {
                    Hist_OnData_Correct_R2off.at(energy).Fill(R2off,xyoff_weight);
                }
            }

        }
        input_file->Close();
    }

    std::cout << "Build exposure map from cosmic rays." << std::endl;
    for (int run=0;run<Data_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        string filename;
        filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

        //std::cout << "Get telescope pointing RA and Dec..." << std::endl;
        pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
        //if (!TString(target).Contains("Proton")) tele_point_ra_dec = GetRunRaDec(filename,int(Data_runlist[run].second));
        tele_point_ra_dec = Data_runlist_RaDec[run];
        run_tele_point_ra = tele_point_ra_dec.first;
        run_tele_point_dec = tele_point_ra_dec.second;
        pair<double,double> tele_point_ra_dec_imposter = std::make_pair(0,0);
        if (doImposter) 
        {
            tele_point_ra_dec_imposter = GetRunRaDec(filename,int(Data_runlist[run].second));
        }

        TFile*  input_file = TFile::Open(filename.c_str());
        TString root_file = "run_"+string(run_number)+"/stereo/data_on";
        TTree* Data_tree = (TTree*) input_file->Get(root_file);
        if (!Data_tree)
        {
            std::cout << "TTree does not exist: " << root_file << std::endl;
            continue;
        }
        SetEventDisplayTreeBranch(Data_tree);
        vector<pair<double,double>> timecut_thisrun = GetRunTimecuts(Data_runlist[run].second);
        double tele_elev = GetRunElevAzim(filename,int(Data_runlist[run].second)).first;

        Data_tree->GetEntry(0);
        double time_0 = Time;
        Data_tree->GetEntry(Data_tree->GetEntries()-1);
        double time_1 = Time;
        double exposure_thisrun = GetRunUsableTime(filename,Data_runlist[run].second)/3600.;

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
            if (doImposter)
            {
                ra_sky_imposter = tele_point_ra_dec_imposter.first+Xoff_derot;
                dec_sky_imposter = tele_point_ra_dec_imposter.second+Yoff_derot;
            }
            //if (doRaster)
            //{
            //    double delta_phi = 2*M_PI*double(entry)/double(Data_tree->GetEntries());
            //    double delta_r = 1.0;
            //    ra_sky += delta_r*cos(delta_phi);
            //    dec_sky += delta_r*sin(delta_phi);
            //}
            // redefine theta2
            theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
            pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            int elevation = Hist_Elev.FindBin(tele_elev)-1;
            int year = Hist_MJD.FindBin(MJD)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            if (elevation<0) continue;
            if (elevation>=N_elev_bins) continue;
            if (!SelectNImages()) continue;
            if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
            if (SizeSecondMax<SizeSecondMax_Cut) continue;
            //if (EmissionHeight<6.) continue;
            double shower_depth = GetShowerDepth(EmissionHeight,tele_elev);
            //if (shower_depth>4.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            //if (pow(Xcore*Xcore+Ycore*Ycore,0.5)<100) continue;
            //if (TString(target).Contains("Crab") && theta2<0.3) continue;
            //if (TString(target).Contains("Mrk421") && theta2<0.3) continue;
            //if (R2off>4.) continue;
            MSCW = RescaleMSCW(MSCW, R2off, MSCW_rescale[energy]);
            MSCL = RescaleMSCW(MSCL, R2off, MSCL_rescale[energy]);
            
            int bin_energy = Hist_CRDark_Energy.at(energy).FindBin(ErecS*1000.);
            int big_bin_xoff = Hist_SRCRDarkRatio_XYoff_Smooth.at(energy).GetXaxis()->FindBin(Xoff);
            int big_bin_yoff = Hist_SRCRDarkRatio_XYoff_Smooth.at(energy).GetYaxis()->FindBin(Yoff);
            int big_bin_ra = Hist_SRCRDarkRatio_RaDec_Smooth.at(energy).GetXaxis()->FindBin(ra_sky);
            int big_bin_dec = Hist_SRCRDarkRatio_RaDec_Smooth.at(energy).GetYaxis()->FindBin(dec_sky);
            double data_cr_content_energy = Hist_SRDark_Energy.at(energy).GetBinContent(bin_energy);
            double dark_cr_content_energy = Hist_CRDark_Energy.at(energy).GetBinContent(bin_energy);
            double xyoff_weight = Hist_SRCRDarkRatio_XYoff_Smooth.at(energy).GetBinContent(big_bin_xoff,big_bin_yoff);
            double radec_weight = Hist_SRCRDarkRatio_RaDec_Smooth.at(energy).GetBinContent(big_bin_ra,big_bin_dec);
            double energy_weight = 1.;
            if (dark_cr_content_energy>0.)
            {
                energy_weight = data_cr_content_energy/dark_cr_content_energy;
            }

            if (ControlSelectionTheta2())
            {
                if (R2off<camera_theta2_cut_upper)
                {
                    Hist_OnData_Expo_Skymap.at(energy).Fill(ra_sky,dec_sky,xyoff_weight);
                }
            }

        }
        input_file->Close();
    }

    std::cout << "Build templates from cosmic rays." << std::endl;
    for (int run=0;run<Data_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        string filename;
        filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

        //std::cout << "Get telescope pointing RA and Dec..." << std::endl;
        pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
        //if (!TString(target).Contains("Proton")) tele_point_ra_dec = GetRunRaDec(filename,int(Data_runlist[run].second));
        tele_point_ra_dec = Data_runlist_RaDec[run];
        run_tele_point_ra = tele_point_ra_dec.first;
        run_tele_point_dec = tele_point_ra_dec.second;
        pair<double,double> tele_point_ra_dec_imposter = std::make_pair(0,0);
        if (doImposter) 
        {
            tele_point_ra_dec_imposter = GetRunRaDec(filename,int(Data_runlist[run].second));
        }

        TFile*  input_file = TFile::Open(filename.c_str());
        TString root_file = "run_"+string(run_number)+"/stereo/data_on";
        TTree* Data_tree = (TTree*) input_file->Get(root_file);
        if (!Data_tree)
        {
            std::cout << "TTree does not exist: " << root_file << std::endl;
            continue;
        }
        SetEventDisplayTreeBranch(Data_tree);
        vector<pair<double,double>> timecut_thisrun = GetRunTimecuts(Data_runlist[run].second);
        double tele_elev = GetRunElevAzim(filename,int(Data_runlist[run].second)).first;
        Data_tree->GetEntry(0);
        double time_0 = Time;
        Data_tree->GetEntry(Data_tree->GetEntries()-1);
        double time_1 = Time;

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
            if (doImposter)
            {
                ra_sky_imposter = tele_point_ra_dec_imposter.first+Xoff_derot;
                dec_sky_imposter = tele_point_ra_dec_imposter.second+Yoff_derot;
            }
            //if (doRaster)
            //{
            //    double delta_phi = 2*M_PI*double(entry)/double(Data_tree->GetEntries());
            //    double delta_r = 1.0;
            //    ra_sky += delta_r*cos(delta_phi);
            //    dec_sky += delta_r*sin(delta_phi);
            //}
            // redefine theta2
            theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
            pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            int elevation = Hist_Elev.FindBin(tele_elev)-1;
            int year = Hist_MJD.FindBin(MJD)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            if (elevation<0) continue;
            if (elevation>=N_elev_bins) continue;
            if (!SelectNImages()) continue;
            if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
            if (SizeSecondMax<SizeSecondMax_Cut) continue;
            //if (EmissionHeight<6.) continue;
            double shower_depth = GetShowerDepth(EmissionHeight,tele_elev);
            //if (shower_depth>4.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            //if (pow(Xcore*Xcore+Ycore*Ycore,0.5)<100) continue;
            //if (TString(target).Contains("Crab") && theta2<0.3) continue;
            //if (TString(target).Contains("Mrk421") && theta2<0.3) continue;
            //if (R2off>4.) continue;
            MSCW = RescaleMSCW(MSCW, R2off, MSCW_rescale[energy]);
            MSCL = RescaleMSCW(MSCL, R2off, MSCL_rescale[energy]);

            double norm_syst_err = Hist_NormSystErr.at(elevation).GetBinContent(energy+1); 
            double shape_syst_err[N_integration_radii] = {0.,0.,0.,0.,0.};
            for (int xybin=0;xybin<N_integration_radii;xybin++) 
            {
                //int xoff_bin = Hist_ShapeSystErr.at(energy).at(xybin).GetXaxis()->FindBin(Xoff);
                //int yoff_bin = Hist_ShapeSystErr.at(energy).at(xybin).GetYaxis()->FindBin(Yoff);
                //shape_syst_err[xybin] = Hist_ShapeSystErr.at(energy).at(xybin).GetBinContent(xoff_bin,yoff_bin);
                int roff_bin = Hist_ShapeSystErr_1D.at(energy).at(xybin).GetXaxis()->FindBin(pow(R2off,0.5));
                shape_syst_err[xybin] = Hist_ShapeSystErr_1D.at(energy).at(xybin).GetBinContent(roff_bin);
            }
            
            int bin_energy = Hist_CRDark_Energy.at(energy).FindBin(ErecS*1000.);
            int big_bin_xoff = Hist_SRCRDarkRatio_XYoff_Smooth.at(energy).GetXaxis()->FindBin(Xoff);
            int big_bin_yoff = Hist_SRCRDarkRatio_XYoff_Smooth.at(energy).GetYaxis()->FindBin(Yoff);
            int big_bin_ra = Hist_SRCRDarkRatio_RaDec_Smooth.at(energy).GetXaxis()->FindBin(ra_sky);
            int big_bin_dec = Hist_SRCRDarkRatio_RaDec_Smooth.at(energy).GetYaxis()->FindBin(dec_sky);
            double data_cr_content_energy = Hist_SRDark_Energy.at(energy).GetBinContent(bin_energy);
            double dark_cr_content_energy = Hist_CRDark_Energy.at(energy).GetBinContent(bin_energy);
            double xyoff_weight = Hist_SRCRDarkRatio_XYoff_Smooth.at(energy).GetBinContent(big_bin_xoff,big_bin_yoff);
            double radec_weight = Hist_SRCRDarkRatio_RaDec_Smooth.at(energy).GetBinContent(big_bin_ra,big_bin_dec);
            double energy_weight = 1.;
            if (dark_cr_content_energy>0.)
            {
                energy_weight = data_cr_content_energy/dark_cr_content_energy;
            }
            int bin_ra = Hist_OnData_Expo_Skymap.at(energy).GetXaxis()->FindBin(ra_sky);
            int bin_dec = Hist_OnData_Expo_Skymap.at(energy).GetYaxis()->FindBin(dec_sky);

            if (FoV(doImposter))
            {
                if (ControlSelectionTheta2())
                {
                    if (R2off<1.0*1.0)
                    {
                        Hist_OnData_CR_Energy_CamCenter.at(energy).Fill(ErecS*1000.,energy_weight);
                    }
                    Hist_OnData_CR_Skymap_Theta2.at(energy).Fill(theta2,xyoff_weight);
                    Hist_NormSyst_Skymap_Theta2.at(energy).Fill(theta2,xyoff_weight*norm_syst_err);
                    Hist_ShapeSyst_Skymap_Theta2.at(energy).Fill(theta2,xyoff_weight*shape_syst_err[0]);
                    Hist_OnData_CR_Yoff.at(energy).Fill(Yoff,xyoff_weight);
                    Hist_OnData_CR_Xoff.at(energy).Fill(Xoff,xyoff_weight);
                    Hist_OnData_CR_XYoff.at(energy).Fill(Xoff,Yoff,xyoff_weight);
                    Hist_OnData_CR_R2off.at(energy).Fill(R2off,xyoff_weight);
                    Hist_OnData_CR_Yoff_Raw.at(energy).Fill(Yoff,1.);
                    Hist_OnData_CR_Skymap.at(energy).Fill(ra_sky,dec_sky,xyoff_weight);
                    Hist_NormSyst_Skymap.at(energy).Fill(ra_sky,dec_sky,xyoff_weight*norm_syst_err);
                    for (int xybin=0;xybin<N_integration_radii;xybin++) 
                    {
                        Hist_ShapeSyst_Skymap.at(energy).at(xybin).Fill(ra_sky,dec_sky,xyoff_weight*shape_syst_err[xybin]);
                    }
                    Hist_OnData_CR_Skymap_Galactic.at(energy).Fill(evt_l_b.first,evt_l_b.second,xyoff_weight);
                    Hist_OnData_CR_Energy.at(energy).Fill(ErecS*1000.,xyoff_weight);
                    Hist_OnData_CR_Zenith.at(energy).Fill(Shower_Ze,xyoff_weight);
                    //Hist_OnData_CR_Height.at(energy).Fill(EmissionHeight,energy_weight);
                    //Hist_OnData_CR_Depth.at(energy).Fill(shower_depth,energy_weight);
                    //Hist_OnData_CR_Rcore.at(energy).Fill(pow(Xcore*Xcore+Ycore*Ycore,0.5),energy_weight);
                    for (int nth_roi=0;nth_roi<roi_ra.size();nth_roi++)
                    {
                        if (RoIFoV(nth_roi)) 
                        {
                            Hist_OnData_CR_RoI_Energy.at(nth_roi).at(energy).Fill(ErecS*1000.,xyoff_weight);
                            Hist_NormSyst_RoI_Energy.at(nth_roi).at(energy).Fill(ErecS*1000.,xyoff_weight*norm_syst_err);
                            double roi_bin_size = roi_radius_outer.at(nth_roi);
                            int xybin_up = N_integration_radii-1;
                            double bin_size_up = integration_radii[xybin_up];
                            double bin_size_low = integration_radii[xybin_up-1];
                            double shape_syst_err_intpl = shape_syst_err[xybin_up-1];
                            for (int xybin=1;xybin<N_integration_radii;xybin++) 
                            {
                                double current_bin_size = integration_radii[xybin];
                                if (current_bin_size>roi_bin_size)
                                {
                                    xybin_up = xybin;
                                    bin_size_up = integration_radii[xybin];
                                    bin_size_low = integration_radii[xybin-1];
                                    shape_syst_err_intpl = shape_syst_err[xybin_up-1] + (shape_syst_err[xybin_up]-shape_syst_err[xybin_up-1])/(bin_size_up-bin_size_low)*(roi_bin_size-bin_size_low);
                                    break;
                                }
                            }
                            Hist_ShapeSyst_RoI_Energy.at(nth_roi).at(energy).Fill(ErecS*1000.,xyoff_weight*shape_syst_err_intpl);
                            Hist_OnData_CR_RoI_MJD.at(nth_roi).at(energy).Fill(MJD,xyoff_weight);
                        }
                        theta2_roi = pow(ra_sky-roi_ra.at(nth_roi),2)+pow(dec_sky-roi_dec.at(nth_roi),2);
                        Hist_OnData_CR_Skymap_RoI_Theta2.at(nth_roi).at(energy).Fill(theta2_roi,xyoff_weight);
                        Hist_NormSyst_Skymap_RoI_Theta2.at(nth_roi).at(energy).Fill(theta2_roi,xyoff_weight*norm_syst_err);
                        Hist_ShapeSyst_Skymap_RoI_Theta2.at(nth_roi).at(energy).Fill(theta2_roi,xyoff_weight*shape_syst_err[0]);
                        proj_x_roi = 1.*(ra_sky-roi_ra.at(nth_roi))+0.*(dec_sky-roi_dec.at(nth_roi));
                        proj_y_roi = 0.*(ra_sky-roi_ra.at(nth_roi))+1.*(dec_sky-roi_dec.at(nth_roi));
                        if (roi_name.at(nth_roi)=="Geminga Pulsar")
                        {
                            proj_x_roi = cos(atan2(97.,138.))*(ra_sky-roi_ra.at(nth_roi))+sin(atan2(97.,138.))*(dec_sky-roi_dec.at(nth_roi));
                            proj_y_roi = -sin(atan2(97.,138.))*(ra_sky-roi_ra.at(nth_roi))+cos(atan2(97.,138.))*(dec_sky-roi_dec.at(nth_roi));
                        }
                        if (abs(proj_y_roi)<roi_radius_outer.at(nth_roi))
                        {
                            Hist_OnData_CR_Skymap_RoI_X.at(nth_roi).at(energy).Fill(proj_x_roi,xyoff_weight);
                        }
                        if (abs(proj_x_roi)<roi_radius_outer.at(nth_roi))
                        {
                            Hist_OnData_CR_Skymap_RoI_Y.at(nth_roi).at(energy).Fill(proj_y_roi,xyoff_weight);
                        }
                    }
                }
            }

        }
        input_file->Close();
    }


    std::cout << "Prepare ON run samples..." << std::endl;
    int final_runs = 0;
    for (int run=0;run<Data_runlist.size();run++)
    {
        std::cout << "Prepare run ..." << int(Data_runlist[run].second) << std::endl;
        final_runs += 1;

        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        string filename;
        filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

        //std::cout << "Get telescope pointing RA and Dec..." << std::endl;
        pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
        //if (!TString(target).Contains("Proton")) tele_point_ra_dec = GetRunRaDec(filename,int(Data_runlist[run].second));
        tele_point_ra_dec = Data_runlist_RaDec[run];
        run_tele_point_ra = tele_point_ra_dec.first;
        run_tele_point_dec = tele_point_ra_dec.second;
        pair<double,double> tele_point_ra_dec_imposter = std::make_pair(0,0);
        if (doImposter) 
        {
            tele_point_ra_dec_imposter = GetRunRaDec(filename,int(Data_runlist[run].second));
        }

        TFile*  input_file = TFile::Open(filename.c_str());
        TString root_file = "run_"+string(run_number)+"/stereo/data_on";
        TTree* Data_tree = (TTree*) input_file->Get(root_file);
        if (!Data_tree)
        {
            std::cout << "TTree does not exist: " << root_file << std::endl;
            continue;
        }
        SetEventDisplayTreeBranch(Data_tree);
        vector<pair<double,double>> timecut_thisrun = GetRunTimecuts(Data_runlist[run].second);

        double NSB_thisrun = GetRunPedestalVar(int(Data_runlist[run].second));
        Hist_Data_NSB.Fill(NSB_thisrun);
        double tele_elev = GetRunElevAzim(filename,int(Data_runlist[run].second)).first;
        double tele_azim = GetRunElevAzim(filename,int(Data_runlist[run].second)).second;

        // Get effective area and livetime and determine the cosmic electron counts for this run.
        //std::cout << "Get effective area and livetime..." << std::endl;
        double wobble_angle = 0.5;
        if (TString(target).Contains("Offset"))
        {
            wobble_angle = 1.0;
        }
        TH1* i_hEffAreaP = ( TH1* )getEffAreaHistogram(input_file,Data_runlist[run].second, wobble_angle);

        Data_tree->GetEntry(0);
        double time_0 = Time;
        Data_tree->GetEntry(Data_tree->GetEntries()-1);
        double time_1 = Time;
        double exposure_thisrun = GetRunUsableTime(filename,Data_runlist[run].second)/3600.;

        for (int e=0;e<N_energy_fine_bins;e++) 
        {
            double eff_area = i_hEffAreaP->GetBinContent( i_hEffAreaP->FindBin( log10(0.5*(energy_fine_bins[e]+energy_fine_bins[e+1])/1000.)));
            Hist_EffArea.SetBinContent(e+1,Hist_EffArea.GetBinContent(e+1)+eff_area*(3600.*exposure_thisrun));
        }
        for (int e=0;e<N_energy_bins;e++) 
        {
            double eff_area = i_hEffAreaP->GetBinContent( i_hEffAreaP->FindBin( log10(0.5*(energy_bins[e]+energy_bins[e+1])/1000.)));
            double binx_size = Hist_Photon_Exp_Skymap.at(e).GetXaxis()->GetBinCenter(2)-Hist_Photon_Exp_Skymap.at(e).GetXaxis()->GetBinCenter(1);
            double biny_size = Hist_Photon_Exp_Skymap.at(e).GetYaxis()->GetBinCenter(2)-Hist_Photon_Exp_Skymap.at(e).GetYaxis()->GetBinCenter(1);
            for (int binx=0;binx<Hist_Photon_Exp_Skymap.at(e).GetNbinsX();binx++)
            {
                for (int biny=0;biny<Hist_Photon_Exp_Skymap.at(e).GetNbinsY();biny++)
                {
                    double bin_ra = Hist_Photon_Exp_Skymap.at(e).GetXaxis()->GetBinCenter(binx+1);
                    double bin_dec = Hist_Photon_Exp_Skymap.at(e).GetYaxis()->GetBinCenter(biny+1);
                    double source_extension = 1.0;
                    double offset_ra = 0.;
                    double offset_dec = 0.;
                    double distance = 0.;
                    double radial_density = 0.;
                    double PercentCrab = 0.;
                    double gamma_flux = 0.;
                    double expected_photons = 0.;
                    double expected_photons_local = 0.;
                    double old_content = Hist_Photon_Exp_Skymap.at(e).GetBinContent(binx+1,biny+1);
                    if (GammaModel==0) // no MC signal
                    {
                        expected_photons_local = 0.;
                    }
                    else if (GammaModel==1) // 2 percent Crab, 1 Guassian, RMS = 1 degrees
                    {
                        offset_ra = 0.;
                        offset_dec = 0.;
                        source_extension = 1.;
                        distance = pow(pow(mean_tele_point_ra+offset_ra-bin_ra,2)+pow(mean_tele_point_dec+offset_dec-bin_dec,2),0.5);
                        radial_density = exp(-0.5*distance*distance/(source_extension*source_extension))/(source_extension*source_extension*2.*M_PI);
                        PercentCrab = 2.;
                        gamma_flux = PercentCrab/100.*GetCrabFlux((energy_bins[e+1]+energy_bins[e])/2.);
                        expected_photons = gamma_flux*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
                        expected_photons_local = expected_photons*radial_density*(binx_size*biny_size);
                    }
                    else if (GammaModel==2) // 5 percent Crab, 1 Guassian, RMS = 1 degrees
                    {
                        offset_ra = 0.;
                        offset_dec = 0.;
                        source_extension = 1.;
                        distance = pow(pow(mean_tele_point_ra+offset_ra-bin_ra,2)+pow(mean_tele_point_dec+offset_dec-bin_dec,2),0.5);
                        radial_density = exp(-0.5*distance*distance/(source_extension*source_extension))/(source_extension*source_extension*2.*M_PI);
                        PercentCrab = 5.;
                        gamma_flux = PercentCrab/100.*GetCrabFlux((energy_bins[e+1]+energy_bins[e])/2.);
                        expected_photons = gamma_flux*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
                        expected_photons_local = expected_photons*radial_density*(binx_size*biny_size);
                    }
                    else if (GammaModel==3) // 10 percent Crab, 1 Guassian, RMS = 1 degrees
                    {
                        offset_ra = 0.;
                        offset_dec = 0.;
                        source_extension = 1.;
                        distance = pow(pow(mean_tele_point_ra+offset_ra-bin_ra,2)+pow(mean_tele_point_dec+offset_dec-bin_dec,2),0.5);
                        radial_density = exp(-0.5*distance*distance/(source_extension*source_extension))/(source_extension*source_extension*2.*M_PI);
                        PercentCrab = 10.;
                        gamma_flux = PercentCrab/100.*GetCrabFlux((energy_bins[e+1]+energy_bins[e])/2.);
                        expected_photons = gamma_flux*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
                        expected_photons_local = expected_photons*radial_density*(binx_size*biny_size);
                    }
                    else if (GammaModel==4) // 2 Guassians, RMS = 0.5 degrees
                    {
                        offset_ra = 1.;
                        offset_dec = 1.;
                        source_extension = 0.5;
                        gamma_flux = GetModelFlux((energy_bins[e+1]+energy_bins[e])/2.,0.2,2.5);
                        distance = pow(pow(mean_tele_point_ra+offset_ra-bin_ra,2)+pow(mean_tele_point_dec+offset_dec-bin_dec,2),0.5);
                        radial_density = exp(-0.5*distance*distance/(source_extension*source_extension))/(source_extension*source_extension*2.*M_PI);
                        expected_photons = gamma_flux*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
                        expected_photons_local += expected_photons*radial_density*(binx_size*biny_size);
                        offset_ra = -1.;
                        offset_dec = -1.;
                        source_extension = 0.5;
                        gamma_flux = GetModelFlux((energy_bins[e+1]+energy_bins[e])/2.,0.1,3.5);
                        distance = pow(pow(mean_tele_point_ra+offset_ra-bin_ra,2)+pow(mean_tele_point_dec+offset_dec-bin_dec,2),0.5);
                        radial_density = exp(-0.5*distance*distance/(source_extension*source_extension))/(source_extension*source_extension*2.*M_PI);
                        expected_photons = gamma_flux*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
                        expected_photons_local += expected_photons*radial_density*(binx_size*biny_size);
                    }
                    else if (GammaModel==5) // 1 Guassian, RMS = 1.0 degrees
                    {
                        offset_ra = 0.;
                        offset_dec = 0.;
                        source_extension = 1.0;
                        gamma_flux = GetModelFlux((energy_bins[e+1]+energy_bins[e])/2.,0.2,3.0);
                        distance = pow(pow(mean_tele_point_ra+offset_ra-bin_ra,2)+pow(mean_tele_point_dec+offset_dec-bin_dec,2),0.5);
                        radial_density = exp(-0.5*distance*distance/(source_extension*source_extension))/(source_extension*source_extension*2.*M_PI);
                        expected_photons = gamma_flux*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
                        expected_photons_local += expected_photons*radial_density*(binx_size*biny_size);
                    }
                    else if (GammaModel==6) // 1 Guassian, Geminga model
                    {
                        offset_ra = 0.;
                        offset_dec = 0.;
                        source_extension = 2.0;
                        gamma_flux = GetModelFlux((energy_bins[e+1]+energy_bins[e])/2.,0.35,2.2);
                        distance = pow(pow(mean_tele_point_ra+offset_ra-bin_ra,2)+pow(mean_tele_point_dec+offset_dec-bin_dec,2),0.5);
                        radial_density = exp(-0.5*distance*distance/(source_extension*source_extension))/(source_extension*source_extension*2.*M_PI);
                        expected_photons = gamma_flux*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
                        expected_photons_local += expected_photons*radial_density*(binx_size*biny_size);
                    }
                    Hist_Photon_Exp_Skymap.at(e).SetBinContent(binx+1,biny+1,expected_photons_local+old_content);
                }
            }
        }

        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_Source_Theta2.at(e).Reset();
        }
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
            if (doImposter)
            {
                ra_sky_imposter = tele_point_ra_dec_imposter.first+Xoff_derot;
                dec_sky_imposter = tele_point_ra_dec_imposter.second+Yoff_derot;
            }
            //if (doRaster)
            //{
            //    double delta_phi = 2*M_PI*double(entry)/double(Data_tree->GetEntries());
            //    double delta_r = 1.0;
            //    ra_sky += delta_r*cos(delta_phi);
            //    dec_sky += delta_r*sin(delta_phi);
            //}
            // redefine theta2
            theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
            pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            int elevation = Hist_Elev.FindBin(tele_elev)-1;
            int year = Hist_MJD.FindBin(MJD)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            if (elevation<0) continue;
            if (elevation>=N_elev_bins) continue;
            if (!SelectNImages()) continue;
            if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
            if (SizeSecondMax<SizeSecondMax_Cut) continue;
            //if (EmissionHeight<6.) continue;
            double shower_depth = GetShowerDepth(EmissionHeight,tele_elev);
            //if (shower_depth>4.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            //if (pow(Xcore*Xcore+Ycore*Ycore,0.5)<100) continue;
            MSCW = RescaleMSCW(MSCW, R2off, MSCW_rescale[energy]);
            MSCL = RescaleMSCW(MSCL, R2off, MSCL_rescale[energy]);
            if (FoV(doImposter))
            {
                if (SignalSelectionTheta2())
                {
                    Hist_Source_Theta2.at(energy).Fill(theta2);
                }
            }
        }
        vector<double> source_weight;
        for (int e=0;e<N_energy_bins;e++) 
        {
            int nbins = Hist_Source_Theta2.at(e).GetNbinsX();
            double background_counts = 0.;
            for (int bin=2;bin<=nbins;bin++) // first bin is the signal region
            {
                background_counts += Hist_Source_Theta2.at(e).GetBinContent(bin);
            }
            double background_avg = background_counts/double(nbins-1);
            if (Hist_Source_Theta2.at(e).GetBinContent(1)>0.)
            {
                source_weight.push_back(background_avg/Hist_Source_Theta2.at(e).GetBinContent(1));
            }
            else
            {
                source_weight.push_back(0.);
            }
        }

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
            if (doImposter)
            {
                ra_sky_imposter = tele_point_ra_dec_imposter.first+Xoff_derot;
                dec_sky_imposter = tele_point_ra_dec_imposter.second+Yoff_derot;
            }
            //if (doRaster)
            //{
            //    double delta_phi = 2*M_PI*double(entry)/double(Data_tree->GetEntries());
            //    double delta_r = 1.0;
            //    ra_sky += delta_r*cos(delta_phi);
            //    dec_sky += delta_r*sin(delta_phi);
            //}
            // redefine theta2
            theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
            pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            int elevation = Hist_Elev.FindBin(tele_elev)-1;
            int year = Hist_MJD.FindBin(MJD)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            if (elevation<0) continue;
            if (elevation>=N_elev_bins) continue;
            if (!SelectNImages()) continue;
            if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
            if (SizeSecondMax<SizeSecondMax_Cut) continue;
            //if (EmissionHeight<6.) continue;
            double shower_depth = GetShowerDepth(EmissionHeight,tele_elev);
            //if (shower_depth>4.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            //if (pow(Xcore*Xcore+Ycore*Ycore,0.5)<100) continue;
            //if (TString(target).Contains("Crab") && theta2<0.3) continue;
            //if (TString(target).Contains("Mrk421") && theta2<0.3) continue;
            //if (R2off>4.) continue;
            MSCW = RescaleMSCW(MSCW, R2off, MSCW_rescale[energy]);
            MSCL = RescaleMSCW(MSCL, R2off, MSCL_rescale[energy]);
            double weight = 1.;
            //if (theta2<source_theta2_cut && SignalSelectionTheta2()) weight = source_weight.at(energy);
            int bin_ra = Hist_OnData_Expo_Skymap.at(energy).GetXaxis()->FindBin(ra_sky);
            int bin_dec = Hist_OnData_Expo_Skymap.at(energy).GetYaxis()->FindBin(dec_sky);

            double eff_area = i_hEffAreaP->GetBinContent( i_hEffAreaP->FindBin( log10(0.5*(energy_bins[energy]+energy_bins[energy+1])/1000.)));
            //double eff_area_weight = (eff_area*exposure_thisrun*3600.*(energy_bins[energy+1]-energy_bins[energy])/1000.);
            double eff_area_weight = (eff_area*(energy_bins[energy+1]-energy_bins[energy])/1000.);

            shower_count += 1.;
            mean_elev += tele_elev;
            double shifted_azim = tele_azim+90.;
            if (shifted_azim>360.) shifted_azim = shifted_azim-360.;
            mean_azim += shifted_azim;
            mean_nsb += NSB_thisrun;
            Hist_Data_ShowerDirection.Fill(Shower_Az,Shower_Ze);
            Hist_Data_ElevNSB.Fill(NSB_thisrun,tele_elev);
            Hist_Data_ElevAzim.Fill(tele_azim,tele_elev);
            Hist_Data_Skymap.Fill(ra_sky,dec_sky);
            Hist_Data_Elev_Skymap.Fill(ra_sky,dec_sky,tele_elev);
            Hist_Data_Azim_Skymap.Fill(ra_sky,dec_sky,tele_azim);
            Hist_Data_NSB_Skymap.Fill(ra_sky,dec_sky,NSB_thisrun);
            if (FoV(doImposter) || Data_runlist[run].first.find("Proton")!=std::string::npos)
            {
                Hist_OnData_MSCLW.at(energy).Fill(MSCL,MSCW,weight);
                if (SourceFoV())
                {
                    Hist_OnData_Point_MSCLW.at(energy).Fill(MSCL,MSCW,weight);
                }
                else if (SourceRingFoV())
                {
                    Hist_OnData_Ring_MSCLW.at(energy).Fill(MSCL,MSCW,weight);
                }
            }
            if (!SignalSelectionTheta2())
            {
                if (FoV(doImposter))
                {
                    Hist_OnData_EffArea_Skymap.at(energy).Fill(ra_sky,dec_sky,eff_area_weight);
                    Hist_OnData_ISR_Skymap.at(energy).Fill(ra_sky,dec_sky,weight);
                    Hist_OnData_ISR_R2off.at(energy).Fill(R2off,weight);
                }
            }
            if (SignalSelectionTheta2())
            {
                if (FoV(doImposter))
                {
                    if (R2off<1.0*1.0)
                    {
                        Hist_OnData_SR_Energy_CamCenter.at(energy).Fill(ErecS*1000.,weight);
                    }
                    Hist_OnData_SR_Yoff.at(energy).Fill(Yoff,weight);
                    Hist_OnData_SR_Xoff.at(energy).Fill(Xoff,weight);
                    Hist_OnData_SR_XYoff.at(energy).Fill(Xoff,Yoff,weight);
                    Hist_OnData_SR_R2off.at(energy).Fill(R2off,weight);
                    Hist_OnData_SR_Skymap_Theta2.at(energy).Fill(theta2,weight);
                    Hist_OnData_SR_Skymap.at(energy).Fill(ra_sky,dec_sky,weight);
                    Hist_OnData_SR_Skymap_Galactic.at(energy).Fill(evt_l_b.first,evt_l_b.second,weight);
                    Hist_OnData_SR_Energy.at(energy).Fill(ErecS*1000.,weight);
                    Hist_OnData_SR_Zenith.at(energy).Fill(Shower_Ze,weight);
                    if (SourceFoV())
                    {
                        Hist_OnData_SR_Height.at(energy).Fill(EmissionHeight,weight);
                        Hist_OnData_SR_Depth.at(energy).Fill(shower_depth,weight);
                        Hist_OnData_SR_Rcore.at(energy).Fill(pow(Xcore*Xcore+Ycore*Ycore,0.5),weight);
                    }
                    else if (SourceRingFoV())
                    {
                        Hist_OnData_CR_Height.at(energy).Fill(EmissionHeight,weight);
                        Hist_OnData_CR_Depth.at(energy).Fill(shower_depth,weight);
                        Hist_OnData_CR_Rcore.at(energy).Fill(pow(Xcore*Xcore+Ycore*Ycore,0.5),weight);
                    }
                    for (int nth_roi=0;nth_roi<roi_ra.size();nth_roi++)
                    {
                        if (RoIFoV(nth_roi)) 
                        {
                            Hist_OnData_SR_RoI_Energy.at(nth_roi).at(energy).Fill(ErecS*1000.,weight);
                            Hist_OnData_SR_RoI_MJD.at(nth_roi).at(energy).Fill(MJD,weight);
                        }
                        theta2_roi = pow(ra_sky-roi_ra.at(nth_roi),2)+pow(dec_sky-roi_dec.at(nth_roi),2);
                        Hist_OnData_SR_Skymap_RoI_Theta2.at(nth_roi).at(energy).Fill(theta2_roi,weight);
                        proj_x_roi = 1.*(ra_sky-roi_ra.at(nth_roi))+0.*(dec_sky-roi_dec.at(nth_roi));
                        proj_y_roi = 0.*(ra_sky-roi_ra.at(nth_roi))+1.*(dec_sky-roi_dec.at(nth_roi));
                        if (roi_name.at(nth_roi)=="Geminga Pulsar")
                        {
                            proj_x_roi = cos(atan2(97.,138.))*(ra_sky-roi_ra.at(nth_roi))+sin(atan2(97.,138.))*(dec_sky-roi_dec.at(nth_roi));
                            proj_y_roi = -sin(atan2(97.,138.))*(ra_sky-roi_ra.at(nth_roi))+cos(atan2(97.,138.))*(dec_sky-roi_dec.at(nth_roi));
                        }
                        if (abs(proj_y_roi)<roi_radius_outer.at(nth_roi))
                        {
                            Hist_OnData_SR_Skymap_RoI_X.at(nth_roi).at(energy).Fill(proj_x_roi,weight);
                        }
                        if (abs(proj_x_roi)<roi_radius_outer.at(nth_roi))
                        {
                            Hist_OnData_SR_Skymap_RoI_Y.at(nth_roi).at(energy).Fill(proj_y_roi,weight);
                        }
                    }
                }
            }
        }
        input_file->Close();
    }

    for (int e=0;e<N_energy_fine_bins;e++) 
    {
        Hist_EffArea.SetBinContent(e+1,Hist_EffArea.GetBinContent(e+1)/(3600.*exposure_hours_usable));
    }

    vector<int> Data_runlist_number;
    vector<string> Data_runlist_name;
    vector<int> Dark_runlist_number;
    vector<string> Dark_runlist_name;
    vector<double> Data_runlist_elev;
    vector<double> Data_runlist_L3Rate;
    vector<double> Data_runlist_NSB;
    vector<int> Data_runlist_MJD;
    for (int run=0;run<Data_runlist.size();run++)
    {
        Data_runlist_name.push_back(Data_runlist[run].first);
        Data_runlist_number.push_back(Data_runlist[run].second);
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        string filename;
        filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");
        Data_runlist_elev.push_back(GetRunElevAzim(filename,int(Data_runlist[run].second)).first);
        Data_runlist_L3Rate.push_back(GetRunL3Rate(int(Data_runlist[run].second)));
        Data_runlist_NSB.push_back(GetRunPedestalVar(int(Data_runlist[run].second)));
        Data_runlist_MJD.push_back(GetRunMJD(filename,int(Data_runlist[run].second)));
        for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
        {
            for (int off_run=0;off_run<Dark_runlist.at(run).at(nth_sample).size();off_run++)
            {
                Dark_runlist_name.push_back(Dark_runlist.at(run).at(nth_sample)[off_run].first);
                Dark_runlist_number.push_back(Dark_runlist.at(run).at(nth_sample)[off_run].second);
            }
        }
    }

    std::cout << "prepare photon template" << std::endl;
    vector<pair<string,int>> PhotonMC_runlist = GetRunList("Photon");
    vector<pair<string,int>> PhotonData_runlist = GetRunList("CrabV5");
    //PhotonData_runlist = SelectONRunList(PhotonData_runlist,TelElev_lower,TelElev_upper,0,0);
    
    vector<TH2D> Hist_GammaMC_MSCLW;
    vector<TH2D> Hist_GammaData_MSCLW;
    vector<TH2D> Hist_GammaDataON_MSCLW;
    vector<TH2D> Hist_GammaDataOFF_MSCLW;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));

        MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-(-1.*MSCW_cut_blind))+MSCW_cut_blind;
        MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-(-1.*MSCL_cut_blind))+MSCL_cut_blind;
        MSCW_plot_lower = -0.5*gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-(-1.*MSCW_cut_blind))-MSCW_cut_blind;
        MSCL_plot_lower = -0.5*gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-(-1.*MSCL_cut_blind))-MSCL_cut_blind;
        N_bins_for_deconv = N_bins_for_deconv_func_E[e];

        Hist_GammaMC_MSCLW.push_back(TH2D("Hist_Stage1_GammaMC_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_GammaData_MSCLW.push_back(TH2D("Hist_Stage1_GammaData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_GammaDataON_MSCLW.push_back(TH2D("Hist_Stage1_GammaDataON_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_GammaDataOFF_MSCLW.push_back(TH2D("Hist_Stage1_GammaDataOFF_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
    }

    // Get Gamma ray MC template
    //for (int run=0;run<PhotonMC_runlist.size();run++)
    //{
    //    char run_number[50];
    //    char Data_observation[50];
    //    sprintf(run_number, "%i", int(PhotonMC_runlist[run].second));
    //    sprintf(Data_observation, "%s", PhotonMC_runlist[run].first.c_str());
    //    string filename;
    //    filename = TString("$SMI_INPUT/"+string(run_number)+".anasum.root");
    //    TFile*  input_file = TFile::Open(filename.c_str());
    //    TString root_file = "run_"+string(run_number)+"/stereo/data_on";
    //    TTree* Data_tree = (TTree*) input_file->Get(root_file);
    //    Data_tree->SetBranchAddress("Xoff",&Xoff);
    //    Data_tree->SetBranchAddress("Yoff",&Yoff);
    //    Data_tree->SetBranchAddress("theta2",&theta2);
    //    Data_tree->SetBranchAddress("ra",&ra_sky);
    //    Data_tree->SetBranchAddress("dec",&dec_sky);
    //    Data_tree->SetBranchAddress("ErecS",&ErecS);
    //    Data_tree->SetBranchAddress("EChi2S",&EChi2S);
    //    Data_tree->SetBranchAddress("MSCW",&MSCW);
    //    Data_tree->SetBranchAddress("MSCL",&MSCL);
    //    Data_tree->SetBranchAddress("NImages",&NImages);
    //    Data_tree->SetBranchAddress("Xcore",&Xcore);
    //    Data_tree->SetBranchAddress("Ycore",&Ycore);
    //    Data_tree->SetBranchAddress("SizeSecondMax",&SizeSecondMax);
    //    Data_tree->SetBranchAddress("Time",&Time);
    //    Data_tree->SetBranchAddress("Shower_Ze",&Shower_Ze);
    //    Data_tree->SetBranchAddress("Shower_Az",&Shower_Az);

    //    for (int entry=0;entry<Data_tree->GetEntries();entry++) 
    //    {
    //        ErecS = 0;
    //        EChi2S = 0;
    //        NImages = 0;
    //        Xcore = 0;
    //        Ycore = 0;
    //        SizeSecondMax = 0;
    //        MSCW = 0;
    //        MSCL = 0;
    //        R2off = 0;
    //        Data_tree->GetEntry(entry);
    //        R2off = Xoff*Xoff+Yoff*Yoff;
    //        ra_sky = mean_tele_point_ra+ra_sky;
    //        dec_sky = mean_tele_point_dec+dec_sky;
    //        pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
    //        int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
    //        int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
    //        if (energy<0) continue;
    //        if (energy>=N_energy_bins) continue;
    //        if (!SelectNImages()) continue;
    //        MSCW = RescaleMSCW(MSCW, R2off, MSCW_rescale[energy]);
    //        MSCL = RescaleMSCW(MSCL, R2off, MSCL_rescale[energy]);
    //        Hist_GammaMC_MSCLW.at(energy).Fill(MSCL,MSCW);
    //        if (MCFoV() && FoV())
    //        {
    //            Hist_Photon_Raw_Skymap.at(energy).Fill(ra_sky,dec_sky);
    //        }
    //    }

    //    for (int entry=0;entry<Data_tree->GetEntries();entry++) 
    //    {
    //        ErecS = 0;
    //        EChi2S = 0;
    //        NImages = 0;
    //        Xcore = 0;
    //        Ycore = 0;
    //        SizeSecondMax = 0;
    //        MSCW = 0;
    //        MSCL = 0;
    //        R2off = 0;
    //        Data_tree->GetEntry(entry);
    //        R2off = Xoff*Xoff+Yoff*Yoff;
    //        ra_sky = mean_tele_point_ra+ra_sky;
    //        dec_sky = mean_tele_point_dec+dec_sky;
    //        pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
    //        int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
    //        int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
    //        if (energy<0) continue;
    //        if (energy>=N_energy_bins) continue;
    //        if (!SelectNImages()) continue;
    //        MSCW = RescaleMSCW(MSCW, R2off, MSCW_rescale[energy]);
    //        MSCL = RescaleMSCW(MSCL, R2off, MSCL_rescale[energy]);

    //        int binx_sky = Hist_Photon_Exp_Skymap.at(energy).GetXaxis()->FindBin(ra_sky);
    //        int biny_sky = Hist_Photon_Exp_Skymap.at(energy).GetYaxis()->FindBin(dec_sky);
    //        double expected_photons = Hist_Photon_Exp_Skymap.at(energy).GetBinContent(binx_sky,biny_sky);
    //        double raw_photons = Hist_Photon_Raw_Skymap.at(energy).GetBinContent(binx_sky,biny_sky);
    //        double gamma_weight = 0.;
    //        if (raw_photons>0.)
    //        {
    //            gamma_weight = expected_photons/raw_photons;
    //        }
    //        if (MCFoV() && FoV())
    //        {
    //            Hist_OnData_MSCLW.at(energy).Fill(MSCL,MSCW,gamma_weight);
    //        }
    //        if (SignalSelectionTheta2())
    //        {
    //            if (MCFoV() && FoV())
    //            {
    //                if (R2off<1.0*1.0)
    //                {
    //                    Hist_OnData_SR_Energy_CamCenter.at(energy).Fill(ErecS*1000.,gamma_weight);
    //                }
    //                Hist_OnData_SR_Skymap_Theta2.at(energy).Fill(theta2,gamma_weight);
    //                Hist_OnData_SR_Skymap.at(energy).Fill(ra_sky,dec_sky,gamma_weight);
    //                Hist_OnData_SR_Skymap_Galactic.at(energy).Fill(evt_l_b.first,evt_l_b.second,gamma_weight);
    //                Hist_OnData_SR_Energy.at(energy).Fill(ErecS*1000.,gamma_weight);
    //                for (int nth_roi=0;nth_roi<roi_ra.size();nth_roi++)
    //                {
    //                    if (nth_roi>0 || !isON)
    //                    {
    //                        if (RoIFoV(nth_roi)) 
    //                        {
    //                            Hist_OnData_SR_RoI_Energy.at(nth_roi).at(energy).Fill(ErecS*1000.,gamma_weight);
    //                        }
    //                        theta2_roi = pow(ra_sky-roi_ra.at(nth_roi),2)+pow(dec_sky-roi_dec.at(nth_roi),2);
    //                        Hist_OnData_SR_Skymap_RoI_Theta2.at(nth_roi).at(energy).Fill(theta2_roi,gamma_weight);
    //                    }
    //                    else
    //                    {
    //                        bool isRoI = false;
    //                        for (int nth_roi2=1;nth_roi2<roi_ra.size();nth_roi2++)
    //                        {
    //                            if (RoIFoV(nth_roi2)) isRoI = true; 
    //                        }
    //                        if (!isRoI)
    //                        {
    //                            Hist_OnData_SR_RoI_Energy.at(nth_roi).at(energy).Fill(ErecS*1000.,gamma_weight);
    //                        }
    //                        theta2_roi = pow(ra_sky-roi_ra.at(nth_roi),2)+pow(dec_sky-roi_dec.at(nth_roi),2);
    //                        Hist_OnData_SR_Skymap_RoI_Theta2.at(nth_roi).at(energy).Fill(theta2_roi,gamma_weight);
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}


    // Get Gamma ray data template
    //std::cout << "PhotonData_runlist.size() = " << PhotonData_runlist.size() << std::endl;
    //for (int run=0;run<PhotonData_runlist.size();run++)
    //{
    //    char run_number[50];
    //    char Data_observation[50];
    //    sprintf(run_number, "%i", int(PhotonData_runlist[run].second));
    //    sprintf(Data_observation, "%s", PhotonData_runlist[run].first.c_str());
    //    string filename;
    //    filename = TString("$SMI_INPUT/"+string(run_number)+".anasum.root");

    //    TFile*  input_file = TFile::Open(filename.c_str());
    //    TString root_file = "run_"+string(run_number)+"/stereo/data_on";
    //    TTree* Data_tree = (TTree*) input_file->Get(root_file);
    //    SetEventDisplayTreeBranch(Data_tree);

    //    for (int entry=0;entry<Data_tree->GetEntries();entry++) 
    //    {
    //        ErecS = 0;
    //        EChi2S = 0;
    //        NImages = 0;
    //        Xcore = 0;
    //        Ycore = 0;
    //        SizeSecondMax = 0;
    //        MSCW = 0;
    //        MSCL = 0;
    //        R2off = 0;
    //        Data_tree->GetEntry(entry);
    //        R2off = Xoff*Xoff+Yoff*Yoff;
    //        Phioff = atan2(Yoff,Xoff)+M_PI;
    //        int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
    //        int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
    //        if (energy<0) continue;
    //        if (energy>=N_energy_bins) continue;
    //        if (!SelectNImages()) continue;
    //        if (SizeSecondMax<400.) continue;
    //        if (EmissionHeight<6.) continue;
    //        if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
    //        MSCW = RescaleMSCW(MSCW, R2off, MSCW_rescale[energy]);
    //        MSCL = RescaleMSCW(MSCL, R2off, MSCL_rescale[energy]);
    //        if (theta2<0.3)
    //        {
    //            Hist_GammaDataON_MSCLW.at(energy).Fill(MSCL,MSCW);
    //        }
    //        else if (theta2>0.3 && theta2<0.5)
    //        {
    //            Hist_GammaDataOFF_MSCLW.at(energy).Fill(MSCL,MSCW);
    //        }
    //    }
    //    input_file->Close();
    //}
    //for (int e=0;e<N_energy_bins;e++) 
    //{

    //    MSCW_plot_upper = gamma_hadron_dim_ratio_w[e]*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
    //    MSCL_plot_upper = gamma_hadron_dim_ratio_l[e]*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;
    //    N_bins_for_deconv = N_bins_for_deconv_func_E[e];

    //    int binx_blind = Hist_GammaDataON_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind)-1;
    //    int binx_lower = Hist_GammaDataON_MSCLW.at(e).GetXaxis()->FindBin(MSCL_plot_lower);
    //    int biny_blind = Hist_GammaDataON_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind)-1;
    //    int biny_lower = Hist_GammaDataON_MSCLW.at(e).GetYaxis()->FindBin(MSCW_plot_lower);
    //    double GammaDataON_SR_Integral = (double) Hist_GammaDataON_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
    //    double GammaDataOFF_SR_Integral = (double) Hist_GammaDataOFF_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
    //    double GammaDataON_all_Integral = (double) Hist_GammaDataON_MSCLW.at(e).Integral();
    //    double GammaDataOFF_all_Integral = (double) Hist_GammaDataOFF_MSCLW.at(e).Integral();
    //    double GammaDataON_CR_Integral = GammaDataON_all_Integral-GammaDataON_SR_Integral;
    //    double GammaDataOFF_CR_Integral = GammaDataOFF_all_Integral-GammaDataOFF_SR_Integral;
    //    double scale = GammaDataON_CR_Integral/GammaDataOFF_CR_Integral;
    //    Hist_GammaDataOFF_MSCLW.at(e).Scale(scale);
    //    Hist_GammaData_MSCLW.at(e).Add(&Hist_GammaDataON_MSCLW.at(e));
    //    Hist_GammaData_MSCLW.at(e).Add(&Hist_GammaDataOFF_MSCLW.at(e),-1.);
    //    for (int binx=0;binx<Hist_GammaData_MSCLW.at(e).GetNbinsX();binx++)
    //    {
    //        for (int biny=0;biny<Hist_GammaData_MSCLW.at(e).GetNbinsY();biny++)
    //        {
    //            double old_content = Hist_GammaData_MSCLW.at(e).GetBinContent(binx+1,biny+1);
    //            double old_error = Hist_GammaData_MSCLW.at(e).GetBinError(binx+1,biny+1);
    //            if (old_content<0)
    //            {
    //                Hist_GammaData_MSCLW.at(e).SetBinContent(binx+1,biny+1,0);
    //                Hist_GammaData_MSCLW.at(e).SetBinError(binx+1,biny+1,0);
    //            }
    //        }
    //    }
    //    std::cout << "Hist_GammaData_MSCLW.at(e).Integral() = " << Hist_GammaData_MSCLW.at(e).Integral() << std::endl;
    //}

    mean_elev = mean_elev/shower_count;
    mean_azim = mean_azim/shower_count;
    mean_azim = mean_azim-90.;  // north will be around +/- 0 and sourth will be around +/- 180
    mean_nsb = mean_nsb/shower_count;

    TFile OutputFile(TString(SMI_OUTPUT)+"/Netflix_"+TString(target)+"_"+TString(output_file_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+TString(mjd_cut_tag)+"_"+ONOFF_tag+group_tag+map_x_tag+map_y_tag+".root","recreate");

    TTree InfoTree("InfoTree","info tree");
    InfoTree.Branch("N_bins_for_deconv",&N_bins_for_deconv,"N_bins_for_deconv/I");
    InfoTree.Branch("MSCW_cut_blind",&MSCW_cut_blind,"MSCW_cut_blind/D");
    InfoTree.Branch("MSCL_cut_blind",&MSCL_cut_blind,"MSCL_cut_blind/D");
    InfoTree.Branch("MSCW_plot_lower",&MSCW_plot_lower,"MSCW_plot_lower/D");
    InfoTree.Branch("MSCL_plot_lower",&MSCL_plot_lower,"MSCL_plot_lower/D");
    InfoTree.Branch("mean_elev",&mean_elev,"mean_elev/D");
    InfoTree.Branch("mean_azim",&mean_azim,"mean_azim/D");
    InfoTree.Branch("mean_nsb",&mean_nsb,"mean_nsb/D");
    InfoTree.Branch("exposure_hours",&exposure_hours,"exposure_hours/D");
    InfoTree.Branch("exposure_hours_usable",&exposure_hours_usable,"exposure_hours_usable/D");
    InfoTree.Branch("exposure_hours_ref",&exposure_hours_ref,"exposure_hours_ref/D");
    InfoTree.Branch("MJD_Start",&MJD_Start,"MJD_Start/I");
    InfoTree.Branch("MJD_End",&MJD_End,"MJD_End/I");
    InfoTree.Branch("mean_tele_point_ra",&mean_tele_point_ra,"mean_tele_point_ra/D");
    InfoTree.Branch("mean_tele_point_dec",&mean_tele_point_dec,"mean_tele_point_dec/D");
    InfoTree.Branch("mean_tele_point_l",&mean_tele_point_l,"mean_tele_point_l/D");
    InfoTree.Branch("mean_tele_point_b",&mean_tele_point_b,"mean_tele_point_b/D");
    InfoTree.Branch("Data_runlist_name","std::vector<std::string>",&Data_runlist_name);
    InfoTree.Branch("Data_runlist_number","std::vector<int>",&Data_runlist_number);
    InfoTree.Branch("Dark_runlist_name","std::vector<std::string>",&Dark_runlist_name);
    InfoTree.Branch("Dark_runlist_number","std::vector<int>",&Dark_runlist_number);
    InfoTree.Branch("Data_runlist_MJD","std::vector<int>",&Data_runlist_MJD);
    InfoTree.Branch("Data_runlist_elev","std::vector<double>",&Data_runlist_elev);
    InfoTree.Branch("Data_runlist_NSB","std::vector<double>",&Data_runlist_NSB);
    InfoTree.Branch("Data_runlist_L3Rate","std::vector<double>",&Data_runlist_L3Rate);
    InfoTree.Branch("Data_runlist_exposure","std::vector<double>",&Data_runlist_exposure);
    InfoTree.Branch("roi_name","std::vector<std::string>",&roi_name);
    InfoTree.Branch("roi_ra","std::vector<double>",&roi_ra);
    InfoTree.Branch("roi_dec","std::vector<double>",&roi_dec);
    InfoTree.Branch("roi_radius_inner","std::vector<double>",&roi_radius_inner);
    InfoTree.Branch("roi_radius_outer","std::vector<double>",&roi_radius_outer);
    InfoTree.Branch("Skymap_size",&Skymap_size,"Skymap_size/D");
    InfoTree.Branch("Skymap_nbins",&Skymap_nbins,"Skymap_nbins/I");
    InfoTree.Fill();
    InfoTree.Write();

    TTree StarTree("StarTree","star tree");
    double star_ra;
    double star_dec;
    double star_brightness;
    StarTree.Branch("star_ra",&star_ra,"star_ra/D");
    StarTree.Branch("star_dec",&star_dec,"star_dec/D");
    StarTree.Branch("star_brightness",&star_brightness,"star_brightness/D");
    for (int star=0;star<BrightStars_Data.size();star++)
    {
        star_ra = BrightStars_Data.at(star).at(0);
        star_dec = BrightStars_Data.at(star).at(1);
        star_brightness = BrightStars_Data.at(star).at(3);
        StarTree.Fill();
    }
    StarTree.Write();
    TTree FaintStarTree("FaintStarTree","faint star tree");
    double faint_star_ra;
    double faint_star_dec;
    double faint_star_brightness;
    FaintStarTree.Branch("faint_star_ra",&faint_star_ra,"faint_star_ra/D");
    FaintStarTree.Branch("faint_star_dec",&faint_star_dec,"faint_star_dec/D");
    FaintStarTree.Branch("faint_star_brightness",&faint_star_brightness,"faint_star_brightness/D");
    for (int star=0;star<FaintStars_Data.size();star++)
    {
        faint_star_ra = FaintStars_Data.at(star).at(0);
        faint_star_dec = FaintStars_Data.at(star).at(1);
        faint_star_brightness = FaintStars_Data.at(star).at(3);
        FaintStarTree.Fill();
    }
    FaintStarTree.Write();

    Hist_Dark_ElevNSB.Write();
    Hist_Dark_ElevAzim.Write();
    Hist_Data_ElevNSB.Write();
    Hist_Data_ElevAzim.Write();
    Hist_Data_Skymap.Write();
    Hist_Data_Elev_Skymap.Write();
    Hist_Data_Azim_Skymap.Write();
    Hist_Data_NSB_Skymap.Write();
    Hist_EffArea.Write();
    for (int e=0;e<N_energy_bins;e++) 
    {
        Hist_GammaMC_MSCLW.at(e).Write();
        Hist_GammaData_MSCLW.at(e).Write();
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        std::cout << "Hist_OnData_MSCLW.at(e).Integral() = " << Hist_OnData_MSCLW.at(e).Integral() << std::endl;
        Hist_OnData_MSCLW.at(e).Write();
        Hist_OnData_Point_MSCLW.at(e).Write();
        Hist_OnData_Ring_MSCLW.at(e).Write();
        Hist_OnData_SR_Skymap_Theta2.at(e).Write();
        Hist_OnData_CR_Skymap_Theta2.at(e).Write();
        Hist_NormSyst_Skymap_Theta2.at(e).Write();
        Hist_ShapeSyst_Skymap_Theta2.at(e).Write();
        Hist_OnData_SR_Yoff.at(e).Write();
        Hist_OnData_SR_Xoff.at(e).Write();
        Hist_OnData_SR_XYoff.at(e).Write();
        Hist_OnData_CR_Yoff.at(e).Write();
        Hist_OnData_CR_Xoff.at(e).Write();
        Hist_OnData_SR_R2off.at(e).Write();
        Hist_OnDark_SR_R2off.at(e).Write();
        Hist_OnData_CR_R2off.at(e).Write();
        Hist_OnData_ISR_R2off.at(e).Write();
        Hist_OnData_CR_XYoff.at(e).Write();
        Hist_OnDark_SR_XYoff.at(e).Write();
        Hist_OnData_CR_Yoff_Raw.at(e).Write();
        Hist_OnData_EffArea_Skymap.at(e).Write();
        Hist_OnData_ISR_Skymap.at(e).Write();
        Hist_OnData_SR_Skymap.at(e).Write();
        Hist_OnData_CR_Skymap.at(e).Write();
        Hist_NormSyst_Skymap.at(e).Write();
        for (int xybin=0;xybin<N_integration_radii;xybin++) 
        {
            Hist_ShapeSyst_Skymap.at(e).at(xybin).Write();
        }
        Hist_OnDark_SR_Skymap.at(e).Write();
        Hist_OnData_SR_Skymap_Galactic.at(e).Write();
        Hist_OnData_CR_Skymap_Galactic.at(e).Write();
        Hist_OnData_SR_Energy.at(e).Write();
        Hist_OnData_CR_Energy.at(e).Write();
        Hist_OnDark_SR_Energy.at(e).Write();
        Hist_OnData_SR_Height.at(e).Write();
        Hist_OnData_CR_Height.at(e).Write();
        Hist_OnData_SR_Depth.at(e).Write();
        Hist_OnData_CR_Depth.at(e).Write();
        Hist_OnData_SR_Rcore.at(e).Write();
        Hist_OnData_CR_Rcore.at(e).Write();
        Hist_OnData_SR_Zenith.at(e).Write();
        Hist_OnData_CR_Zenith.at(e).Write();
        Hist_OnData_SR_Energy_CamCenter.at(e).Write();
        Hist_OnData_CR_Energy_CamCenter.at(e).Write();
    }
    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
    {
        for (int e=0;e<N_energy_bins;e++) 
        {
            std::cout << "Hist_OnDark_MSCLW.at(e).Integral() = " << Hist_OnDark_MSCLW.at(nth_sample).at(e).Integral() << std::endl;
            Hist_OnDark_MSCLW.at(nth_sample).at(e).Write();
        }
    }
    for (int nth_roi=0;nth_roi<roi_ra.size();nth_roi++)
    {
        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OnData_SR_RoI_Energy.at(nth_roi).at(e).Write();
            Hist_OnData_CR_RoI_Energy.at(nth_roi).at(e).Write();
            Hist_NormSyst_RoI_Energy.at(nth_roi).at(e).Write();
            Hist_ShapeSyst_RoI_Energy.at(nth_roi).at(e).Write();
            Hist_OnData_SR_RoI_MJD.at(nth_roi).at(e).Write();
            Hist_OnData_CR_RoI_MJD.at(nth_roi).at(e).Write();
            Hist_OnData_SR_Skymap_RoI_Theta2.at(nth_roi).at(e).Write();
            Hist_OnData_CR_Skymap_RoI_Theta2.at(nth_roi).at(e).Write();
            Hist_NormSyst_Skymap_RoI_Theta2.at(nth_roi).at(e).Write();
            Hist_ShapeSyst_Skymap_RoI_Theta2.at(nth_roi).at(e).Write();
            Hist_OnData_SR_Skymap_RoI_X.at(nth_roi).at(e).Write();
            Hist_OnData_CR_Skymap_RoI_X.at(nth_roi).at(e).Write();
            Hist_OnData_SR_Skymap_RoI_Y.at(nth_roi).at(e).Write();
            Hist_OnData_CR_Skymap_RoI_Y.at(nth_roi).at(e).Write();
        }
    }
    
    OutputFile.Close();

    std::cout << "selected runs = " << Data_runlist.size() << std::endl;
    std::cout << "final runs = " << final_runs << std::endl;
    std::cout << "Done." << std::endl;
}

void PrepareDarkData(string target_data, double tel_elev_lower_input, double tel_elev_upper_input, int MJD_start_cut, int MJD_end_cut, double input_theta2_cut_lower, double input_theta2_cut_upper, bool isON, bool doImposter, int GammaModel)
{

    SMI_INPUT = string(std::getenv("SMI_INPUT"));
    SMI_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    SMI_DIR = string(std::getenv("SMI_DIR"));
    SMI_AUX = string(std::getenv("SMI_AUX"));

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

    bool file_exists = true;
    int group_index = 0;
    int n_groups = 0;
    while (file_exists)
    {
        sprintf(group_tag, "_G%d", group_index);
        std::cout << "Reading file " << TString(SMI_OUTPUT)+"/Netflix_RunList_"+TString(target)+"_"+TString(output_file_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+TString(mjd_cut_tag)+"_"+ONOFF_tag+group_tag+".root" << std::endl;
        if (gSystem->AccessPathName(TString(SMI_OUTPUT)+"/Netflix_RunList_"+TString(target)+"_"+TString(output_file_tag)+TString(elev_cut_tag)+TString(theta2_cut_tag)+TString(mjd_cut_tag)+"_"+ONOFF_tag+group_tag+".root"))
        {
            std::cout << "file does not exist." << std::endl;
            file_exists = false;
        }
        else
        {
            std::cout << "file exists." << std::endl;
            n_groups += 1;
        }
        group_index += 1;
    }
    std::cout << "n_groups = " << n_groups << std::endl;

    for (int g_idx=0;g_idx<n_groups;g_idx++)
    {
        std::cout << "===============================================================================" << std::endl;
        std::cout << "Prepare sub-group " << g_idx+1 << "/" << n_groups << std::endl;
        std::cout << "===============================================================================" << std::endl;
        for (int x_idx=0;x_idx<Skymap_normalization_nbins;x_idx++)
        {
            for (int y_idx=0;y_idx<Skymap_normalization_nbins;y_idx++)
            {
                PrepareDarkData_SubGroup(target_data, tel_elev_lower_input, tel_elev_upper_input, MJD_start_cut, MJD_end_cut, input_theta2_cut_lower, input_theta2_cut_upper, isON, doImposter, GammaModel, g_idx, x_idx, y_idx);
            }
        }
    }

}
