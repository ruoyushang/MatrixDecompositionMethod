
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "NetflixParameters.h"

pair<double,double> GetSourceRaDec(TString source_name)
{
    double Source_RA = 0.;
    double Source_Dec = 0.;
    if (source_name.Contains("PSR_J2238_p5903"))
    {
            Source_RA = 339.6177917;
            Source_Dec = 59.0613333;
    }
    if (source_name.Contains("PSR_J1930_p1852"))
    {
            Source_RA = 292.62554;
            Source_Dec = 18.87058;
    }
    if (source_name.Contains("PSR_J2022_p3842"))
    {
            Source_RA = 305.59037;
            Source_Dec = 38.704117;
    }
    if (source_name.Contains("PSR_J0023_p09"))
    {
            Source_RA = 5.82032291;
            Source_Dec = 9.38996121;
    }
    if (source_name.Contains("PSR_J0248_p6021"))
    {
            Source_RA = 42.077571;
            Source_Dec = 60.359644;
    }
    if (source_name.Contains("PSR_J0633_p0632"))
    {
            Source_RA = 98.43421;
            Source_Dec = 6.5430;
    }
    if (source_name.Contains("LSI_p61_303"))
    {
            Source_RA = 40.1416667;
            Source_Dec = 61.2569444;
    }
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
            //Source_Dec = 22.660;
            Source_RA = 94.213;
            Source_Dec = 22.503;
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
    if (source_name.Contains("UrsaMinor"))
    {
            Source_RA = 227.2854167;
            Source_Dec = 67.2225000;
    }
    if (source_name.Contains("RGB_J0710_p591"))
    {
            Source_RA = 107.6100000;
            Source_Dec = 59.1500000;
    }
    if (source_name.Contains("LHAASO_J2032"))
    {
            Source_RA = 308.0500000;
            Source_Dec = 41.0500000;
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
            //Source_RA = 185.382;
            //Source_Dec = 28.233;
            Source_RA = (185.360+184.616)/2.;
            Source_Dec = (30.191+30.130)/2.;
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
            //Source_RA = 286.975;
            //Source_Dec = 6.269;
            Source_RA = 286.975;
            Source_Dec = 6.03777777778;
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

vector<std::pair<string,int>> GetRunListFromFile(string source)
{
    string line;
    char delimiter = ',';
    string acc_runnumber = "";
    string acc_source = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;
    vector<std::pair<string,int>> list;

    string SMI_DIR;
    SMI_DIR = string(std::getenv("SMI_DIR"));
    ifstream myfile (SMI_DIR+"/runlist_backup/RunList_"+source+".txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            acc_source = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==1)
                {
                    acc_runnumber += line[i];
                }
                else if (nth_delimiter==0)
                {
                    acc_source += line[i];
                }
                if (i==line.size()-1)
                {
                    double runnumber = std::stod(acc_runnumber,&sz);
                    list.push_back(std::make_pair(acc_source,runnumber));
                }
            }
            nth_line += 1;
        }
        myfile.close();
    }
    else cout << "Unable to open file"; 

    return list;
}

vector<std::pair<string,int>> GetRunList(string source) {
        vector<std::pair<string,int>> list;
        vector<std::pair<string,int>> list_temp;
        if (UseDL3Tree)
        {
            if (source.find("CrabV6") != std::string::npos) {
                list = GetRunListFromFile("CrabV6_DL3");
            }
            if (source.find("Crab_Elev70_V6") != std::string::npos) {
                list = GetRunListFromFile("CrabV6_DL3");
            }
            if (source.find("Crab_Elev50_V6") != std::string::npos) {
                list = GetRunListFromFile("CrabV6_DL3");
            }
            if (source.find("Crab_Elev30_V6") != std::string::npos) {
                list = GetRunListFromFile("CrabV6_DL3");
            }
        }
        else
        {
            if (source.find("CrabV6") != std::string::npos) {
                list = GetRunListFromFile("CrabV6");
            }
            if (source.find("Crab_Elev70_V6") != std::string::npos) {
                list = GetRunListFromFile("CrabV6");
            }
            if (source.find("Crab_Elev50_V6") != std::string::npos) {
                list = GetRunListFromFile("CrabV6");
            }
            if (source.find("Crab_Elev30_V6") != std::string::npos) {
                list = GetRunListFromFile("CrabV6");
            }
        }
        if (source.find("Segue1V6") != std::string::npos)
        {
            list = GetRunListFromFile("Segue1V6");
        }
        if (source.find("Ton599") != std::string::npos) {
            list = GetRunListFromFile("Ton599");
        }
        if (source.find("3C264V6") != std::string::npos) {
            list = GetRunListFromFile("3C264V6");
        }
        if (source.find("PKS1424V6") != std::string::npos) {
            list = GetRunListFromFile("PKS1424V6");
        }
        if (source.find("PG1553V6") != std::string::npos) {
            list = GetRunListFromFile("PG1553V6");
        }
        if (source.find("IC443HotSpotV6") != std::string::npos) {
            list = GetRunListFromFile("IC443HotSpotV6");
        }
        if (source.find("H1426V6") != std::string::npos) {
            list = GetRunListFromFile("H1426V6");
        }
        if (source.find("CrabRHVV6") != std::string::npos) {
            list = GetRunListFromFile("CrabRHVV6");
        }
        if (source.find("Crab_Offset_1p0_V6") != std::string::npos) {
            list = GetRunListFromFile("Crab_Offset_1p0_V6");
        }
        if (source.find("Crab_Offset_1p5_V6") != std::string::npos) {
            list = GetRunListFromFile("Crab_Offset_1p5_V6");
        }
        if (source.find("SNR_G150p3Plus04p5_V6") != std::string::npos) {
            list = GetRunListFromFile("SNR_G150p3Plus04p5_V6");
        }
        if (source.find("SNR_G150p3Plus04p5_Jamie") != std::string::npos) {
            list = GetRunListFromFile("SNR_G150p3Plus04p5_Jamie");
        }
        if (source.find("S3_1227_V6") != std::string::npos)
        {
            list = GetRunListFromFile("S3_1227_V6");
        }
        if (source.find("1ES0229V6") != std::string::npos)
        {
            list = GetRunListFromFile("1ES0229V6");
        }
        if (source.find("SgrAV6") != std::string::npos)
        {
            list = GetRunListFromFile("SgrAV6");
        }
        if (source.find("GemingaV6") != std::string::npos)
        {
            list = GetRunListFromFile("GemingaV6");
        }
        if (source.find("CasAV6") != std::string::npos)
        {
            list = GetRunListFromFile("CasAV6");
        }
        if (source=="Photon")
        {
            list = GetRunListFromFile("Photon");
        }
        if (source=="Proton_NSB200")
        {
            list = GetRunListFromFile("Proton_NSB200");
        }
        if (source.find("G079") != std::string::npos)
        {
            list = GetRunListFromFile("G079");
        }
        if (source.find("Mrk421") != std::string::npos)
        {
            list = GetRunListFromFile("Mrk421");
        }
        if (source.find("WComaeV6") != std::string::npos)
        {
            list = GetRunListFromFile("WComaeV6");
        }
        if (source.find("1ES1218V6") != std::string::npos)
        {
            list = GetRunListFromFile("1ES1218V6");
        }
        if (source.find("OJ287V6") != std::string::npos)
        {
            list = GetRunListFromFile("OJ287V6");
        }
        if (source.find("1ES1011V6") != std::string::npos)
        {
            list = GetRunListFromFile("1ES1011V6");
        }
        if (source.find("NGC1275V6") != std::string::npos)
        {
            list = GetRunListFromFile("NGC1275V6");
        }
        if (source.find("1ES0647V6") != std::string::npos)
        {
            list = GetRunListFromFile("1ES0647V6");
        }
        if (source.find("TriIIV6") != std::string::npos)
        {
            list = GetRunListFromFile("TriIIV6");
        }
        if (source.find("1ES1440V6") != std::string::npos)
        {
            list = GetRunListFromFile("1ES1440V6");
        }
        if (source.find("1ES1741V6") != std::string::npos)
        {
            list = GetRunListFromFile("1ES1741V6");
        }
        if (source.find("MGRO_J2019_V6") != std::string::npos)
        {
            list = GetRunListFromFile("MGRO_J2019_V6");
        }
        if (source.find("MGRO_J1908_V6") != std::string::npos)
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source.find("SS433_V6") != std::string::npos)
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source.find("SS433Half1_V6") != std::string::npos)
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source.find("SS433Half2_V6") != std::string::npos)
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source.find("MAGIC_J1857_V6") != std::string::npos)
        {
            //list = GetRunListFromFile("MGRO_J1908_V6");
            list = GetRunListFromFile("HESS_J1857_V6");
        }
        if (source.find("MGRO_J2031_V6") != std::string::npos)
        {
            list = GetRunListFromFile("MGRO_J2031_V6");
        }
        if (source.find("CygnusV6") != std::string::npos)
        {
            list = GetRunListFromFile("CygnusV6");
        }
        if (source.find("RBS0413V6") != std::string::npos)
        {
            list = GetRunListFromFile("RBS0413V6");
        }
        if (source.find("PG1553V6") != std::string::npos)
        {
            list = GetRunListFromFile("PG1553V6");
        }
        if (source.find("PKS1441V6") != std::string::npos)
        {
            list = GetRunListFromFile("PKS1441V6");
        }
        if (source.find("MS1221V6") != std::string::npos)
        {
            list = GetRunListFromFile("MS1221V6");
        }
        if (source.find("TychoV6") != std::string::npos)
        {
            list = GetRunListFromFile("TychoV6");
        }
        if (source.find("2HWC_J1953V6") != std::string::npos)
        {
            list = GetRunListFromFile("2HWC_J1953V6");
        }
        if (source.find("2HWC_J1930V6") != std::string::npos)
        {
            list = GetRunListFromFile("2HWC_J1930V6");
        }
        if (source.find("ComaV6") != std::string::npos) {
            list = GetRunListFromFile("ComaV6");
        }
        if (source.find("BoomerangV6") != std::string::npos)
        {
            list = GetRunListFromFile("BoomerangV6");
        }
        if (source.find("BLLacV6") != std::string::npos)
        {
            list = GetRunListFromFile("BLLacV6");
        }
        if (source.find("M82V6") != std::string::npos)
        {
            list = GetRunListFromFile("M82V6");
        }
        if (source.find("M87V6") != std::string::npos)
        {
            list = GetRunListFromFile("M87V6");
        }
        if (source.find("GammaCygniV6") != std::string::npos)
        {
            list = GetRunListFromFile("GammaCygniV6");
        }
        if (source.find("HESS_J1825_V6") != std::string::npos)
        {
            list = GetRunListFromFile("HESS_J1825_V6");
        }
        if (source.find("3C273V6") != std::string::npos)
        {
            list = GetRunListFromFile("3C273V6");
        }
        if (source.find("1ES0502V6") != std::string::npos)
        {
            list = GetRunListFromFile("1ES0502V6");
        }
        if (source.find("DracoV6") != std::string::npos)
        {
            list = GetRunListFromFile("DracoV6");
        }
        if (source.find("PSR_J2238_p5903_V6") != std::string::npos)
        {
            list = GetRunListFromFile("PSR_J2238_p5903_V6");
        }
        if (source.find("PSR_J0248_p6021_V6") != std::string::npos)
        {
            list = GetRunListFromFile("PSR_J0248_p6021_V6");
        }
        if (source.find("PSR_J0633_p0632_V6") != std::string::npos)
        {
            list = GetRunListFromFile("PSR_J0633_p0632_V6");
        }
        if (source.find("LSI_p61_303_V6") != std::string::npos)
        {
            list = GetRunListFromFile("PSR_J0248_p6021_V6");
        }
        if (source.find("LHAASO_J2032_V6") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J2032_V6");
        }
        if (source.find("LHAASO_J2032_Baseline_V6") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J2032_Baseline_V6");
        }
        if (source.find("LHAASO_J2032_Fall2017_V6") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J2032_Fall2017_V6");
        }
        if (source.find("LHAASO_J2108_V6") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J2108_V6");
        }
        if (source.find("LHAASO_J0341_V6") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J0341_V6");
        }
        if (source.find("LHAASO_J1929_V6") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J1929_V6");
        }
        if (source.find("LHAASO_J1843_V6") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J1843_V6");
        }
        if (source.find("LHAASO_J1956_V6") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J1956_V6");
        }
        if (source.find("Perseus_V6") != std::string::npos)
        {
            list = GetRunListFromFile("Perseus_V6");
        }
        if (source.find("PSRB0355plus54_V6") != std::string::npos)
        {
            list = GetRunListFromFile("PSRB0355plus54_V6");
        }
        if (source.find("UrsaMajorIIV6") != std::string::npos)
        {
            list = GetRunListFromFile("UrsaMajorIIV6");
        }
        if (source.find("UrsaMinorV6") != std::string::npos)
        {
            list = GetRunListFromFile("UrsaMinorV6");
        }
        if (source.find("RGB_J0710_p591_V6") != std::string::npos)
        {
            list = GetRunListFromFile("RGB_J0710_p591_V6");
        }
        if (source.find("PSR_J0023_p09_V6") != std::string::npos)
        {
            list = GetRunListFromFile("PSR_J0023_p09_V6");
        }
        if (source.find("PSR_J1930_p1852_V6") != std::string::npos)
        {
            list = GetRunListFromFile("PSR_J1930_p1852_V6");
        }
        if (source.find("PSR_J2022_p3842_V6") != std::string::npos)
        {
            list = GetRunListFromFile("PSR_J2022_p3842_V6");
        }
        if (source=="OffRunsV6")
        {
            if (RHVData)
            {
                list_temp = GetRunListFromFile("CrabRHV_OFF_V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("CrabRHV_OFF_042720222_V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
            }
            else
            {
                if (UseDL3Tree)
                {
                    list_temp = GetRunListFromFile("CasA_Imposter_V6");
                    list.insert(list.end(), list_temp.begin(), list_temp.end());
                    list_temp = GetRunListFromFile("HESS_J1825_Imposter_V6_DQM");
                    list.insert(list.end(), list_temp.begin(), list_temp.end());
                    list_temp = GetRunListFromFile("HESS_J1825_OFF_V6");
                    list.insert(list.end(), list_temp.begin(), list_temp.end());
                    list_temp = GetRunListFromFile("MGRO_J1908_Imposter_V6_DL3");
                    list.insert(list.end(), list_temp.begin(), list_temp.end());
                    list_temp = GetRunListFromFile("UrsaMinorV6");
                    list.insert(list.end(), list_temp.begin(), list_temp.end());
                    list_temp = GetRunListFromFile("RGB_J0710_p591_V6");
                    list.insert(list.end(), list_temp.begin(), list_temp.end());
                    list_temp = GetRunListFromFile("LHAASO_J2032_Imposter_V6");
                    list.insert(list.end(), list_temp.begin(), list_temp.end());
                }
                list_temp = GetRunListFromFile("LowElevationDarkV6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("PKS1424V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("UrsaMajorIIV6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("PG1553V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("3C273V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("1ES0502V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("DracoV6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("BLLacV6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("M82V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("Segue1V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("TriIIV6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("1ES1011V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("3C264V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("H1426V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("1ES0229V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("OJ287V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("1ES0647V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("NGC1275V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                if (!UseDL3Tree)
                {
                    list_temp = GetRunListFromFile("MGRO_J1908_OFF_V6");
                    list.insert(list.end(), list_temp.begin(), list_temp.end());
                    list_temp = GetRunListFromFile("MGRO_J1908_Imposter_V6");
                    list.insert(list.end(), list_temp.begin(), list_temp.end());
                    list_temp = GetRunListFromFile("MGRO_J1908_Imposter_V6_2022June");
                    list.insert(list.end(), list_temp.begin(), list_temp.end());
                }
            }
            
            //list_temp = GetRunListFromFile("RBS0413V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("IC443_MatchedDark_V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("WComae_MatchedDark_V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("MGRO_J1908_MatchedDark_V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("SgrA_MatchedDark_V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("Cygnus_MatchedDark_V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("Geminga_MatchedDark_V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("CasA_MatchedDark_V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("H1426_MatchedDark_V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
        }
        if (source.find("UrsaMinorV5") != std::string::npos)
        {
            list = GetRunListFromFile("UrsaMinorV5");
        }
        if (source.find("RGB_J0710_p591_V5") != std::string::npos)
        {
            list = GetRunListFromFile("RGB_J0710_p591_V5");
        }
        if (source.find("Mrk421V5") != std::string::npos)
        {
            list = GetRunListFromFile("Mrk421V5");
        }
        if (source.find("CTA1V5") != std::string::npos) {
            list = GetRunListFromFile("CTA1V5");
        }
        if (source.find("MGRO_J2019_V5") != std::string::npos) {
            list = GetRunListFromFile("MGRO_J2019_V5");
        }
        if (source.find("MGRO_J1908_V5") != std::string::npos) {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source.find("SS433_V5") != std::string::npos) {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source.find("SS433Half1_V5") != std::string::npos) {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source.find("SS433Half2_V5") != std::string::npos) {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source.find("MAGIC_J1857_V5") != std::string::npos) {
            //list = GetRunListFromFile("MGRO_J1908_V5");
            list = GetRunListFromFile("HESS_J1857_V5");
        }
        if (source.find("IC443HotSpotV5") != std::string::npos)
        {
            list = GetRunListFromFile("IC443HotSpotV5");
        }
        if (source.find("Segue1V5") != std::string::npos)
        {
            list = GetRunListFromFile("Segue1V5");
        }
        if (source.find("GemingaV5") != std::string::npos)
        {
            list = GetRunListFromFile("GemingaV5");
        }
        if (UseDL3Tree)
        {
            if (source.find("CrabV5") != std::string::npos) {
                list = GetRunListFromFile("CrabV5_DL3");
            }
            if (source.find("Crab_Elev70_V5") != std::string::npos)
            {
                list = GetRunListFromFile("CrabV5_DL3");
            }
            if (source.find("Crab_Elev50_V5") != std::string::npos)
            {
                list = GetRunListFromFile("CrabV5_DL3");
            }
            if (source.find("Crab_Elev30_V5") != std::string::npos)
            {
                list = GetRunListFromFile("CrabV5_DL3");
            }
        }
        else
        {
            if (source.find("CrabV5") != std::string::npos)
            {
                list = GetRunListFromFile("CrabV5");
            }
            if (source.find("Crab_Elev70_V5") != std::string::npos)
            {
                list = GetRunListFromFile("CrabV5");
            }
            if (source.find("Crab_Elev50_V5") != std::string::npos)
            {
                list = GetRunListFromFile("CrabV5");
            }
            if (source.find("Crab_Elev30_V5") != std::string::npos)
            {
                list = GetRunListFromFile("CrabV5");
            }
        }
        if (source.find("WComaeV5") != std::string::npos)
        {
            list = GetRunListFromFile("WComaeV5");
        }
        if (source.find("CygnusV5") != std::string::npos)
        {
            list = GetRunListFromFile("CygnusV5");
        }
        if (source.find("MGRO_J2031_V5") != std::string::npos)
        {
            list = GetRunListFromFile("MGRO_J2031_V5");
        }
        if (source.find("PG1553V5") != std::string::npos)
        {
            list = GetRunListFromFile("PG1553V5");
        }
        if (source.find("RBS0413V5") != std::string::npos)
        {
            list = GetRunListFromFile("RBS0413V5");
        }
        if (source.find("1ES0229V5") != std::string::npos)
        {
            list = GetRunListFromFile("1ES0229V5");
        }
        if (source.find("1ES0414V5") != std::string::npos)
        {
            list = GetRunListFromFile("1ES0414V5");
        }
        if (source.find("PKS1424V5") != std::string::npos)
        {
            list = GetRunListFromFile("PKS1424V5");
        }
        if (source.find("M82V5") != std::string::npos)
        {
            list = GetRunListFromFile("M82V5");
        }
        if (source.find("M87V5") != std::string::npos)
        {
            list = GetRunListFromFile("M87V5");
        }
        if (source.find("BoomerangV5") != std::string::npos)
        {
            list = GetRunListFromFile("BoomerangV5");
        }
        if (source.find("BLLacV5") != std::string::npos)
        {
            list = GetRunListFromFile("BLLacV5");
        }
        if (source.find("GammaCygniV5") != std::string::npos)
        {
            list = GetRunListFromFile("GammaCygniV5");
        }
        if (source.find("TychoV5") != std::string::npos)
        {
            list = GetRunListFromFile("TychoV5");
        }
        if (source.find("3C273V5") != std::string::npos)
        {
            list = GetRunListFromFile("3C273V5");
        }
        if (source.find("1ES0502V5") != std::string::npos)
        {
            list = GetRunListFromFile("1ES0502V5");
        }
        if (source.find("DracoV5") != std::string::npos)
        {
            list = GetRunListFromFile("DracoV5");
        }
        if (source.find("LHAASO_J2032_V5") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J2032_V5");
        }
        if (source.find("LHAASO_J2032_Baseline_V5") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J2032_Baseline_V5");
        }
        if (source=="OffRunsV5")
        {
            list_temp = GetRunListFromFile("BLLacV5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("PKS1424V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("1ES0229V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("1ES0414V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("RBS0413V5");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("PG1553V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("Segue1V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("ComaV5");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("M82V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("3C273V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("1ES0502V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("DracoV5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            if (!UseDL3Tree)
            {
                list_temp = GetRunListFromFile("MGRO_J1908_OFF_V5");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("MGRO_J1908_Imposter_V5");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("MGRO_J1908_Imposter_V5_2022June");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
            }
            if (UseDL3Tree)
            {
                list_temp = GetRunListFromFile("MGRO_J1908_Imposter_V5_DL3");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("UrsaMinorV5");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("RGB_J0710_p591_V5");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("LHAASO_J2032_Imposter_V5");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
            }
        }
        if (source.find("WComaeV4") != std::string::npos)
        {
            list = GetRunListFromFile("WComaeV4");
        }
        if (source.find("MGRO_J1908_V4") != std::string::npos) {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source.find("SS433_V4") != std::string::npos) {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source.find("SS433Half1_V4") != std::string::npos) {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source.find("SS433Half2_V4") != std::string::npos) {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source.find("MAGIC_J1857_V4") != std::string::npos) {
            //list = GetRunListFromFile("MGRO_J1908_V4");
            list = GetRunListFromFile("HESS_J1857_V4");
        }
        if (source.find("MGRO_J2031_V4") != std::string::npos) {
            list = GetRunListFromFile("MGRO_J2031_V4");
        }
        if (source.find("IC443HotSpotV4") != std::string::npos)
        {
            list = GetRunListFromFile("IC443HotSpotV4");
        }
        if (source.find("ComaV4") != std::string::npos) {
            list = GetRunListFromFile("ComaV4");
        }
        if (source.find("M82V4") != std::string::npos)
        {
            list = GetRunListFromFile("M82V4");
        }
        if (source.find("BoomerangV4") != std::string::npos)
        {
            list = GetRunListFromFile("BoomerangV4");
        }
        if (source.find("GammaCygniV4") != std::string::npos)
        {
            list = GetRunListFromFile("GammaCygniV4");
        }
        if (source.find("CrabV4") != std::string::npos)
        {
            list = GetRunListFromFile("CrabV4");
        }
        if (source.find("TychoV4") != std::string::npos)
        {
            list = GetRunListFromFile("TychoV4");
        }
        if (source=="OffRunsV4")
        {
            list_temp = GetRunListFromFile("WComae_MatchedDark_V4");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("MGRO_J1908_MatchedDark_V4");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("IC443_MatchedDark_V4");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("ComaV4");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("M82V4");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
        }

        return list;
}
