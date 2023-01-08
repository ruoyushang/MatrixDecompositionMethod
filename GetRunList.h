
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
    if (source_name.Contains("GalacticPlane_All_l30"))
    {
        Source_RA = 281.522;
        Source_Dec = -2.609;
    }
    if (source_name.Contains("GalacticPlane_All_l40"))
    {
        Source_RA = 286.095;
        Source_Dec = 6.287;
    }
    if (source_name.Contains("GalacticPlane_All_l50"))
    {
        Source_RA = 290.829;
        Source_Dec = 15.142;
    }
    if (source_name.Contains("GalacticPlane_All_l60"))
    {
        Source_RA = 295.976;
        Source_Dec = 23.890;
    }
    if (source_name.Contains("GalacticPlane_All_l70"))
    {
        Source_RA = 301.866;
        Source_Dec = 32.442;
    }
    if (source_name.Contains("GalacticPlane_All_l80"))
    {
        Source_RA = 308.971;
        Source_Dec = 40.664;
    }
    if (source_name.Contains("GalacticPlane_All_l90"))
    {
        Source_RA = 318.004;
        Source_Dec = 48.330;
    }
    if (source_name.Contains("GalacticPlane_All_l100"))
    {
        Source_RA = 330.004;
        Source_Dec = 55.050;
    }
    if (source_name.Contains("GalacticPlane_All_l110"))
    {
        Source_RA = 346.131;
        Source_Dec = 60.160;
    }
    if (source_name.Contains("GalacticPlane_All_l120"))
    {
        Source_RA = 6.451;
        Source_Dec = 62.726;
    }
    if (source_name.Contains("GalacticPlane_All_l130"))
    {
        Source_RA = 28.071;
        Source_Dec = 62.034;
    }
    if (source_name.Contains("GalacticPlane_All_l140"))
    {
        Source_RA = 46.813;
        Source_Dec = 58.298;
    }
    if (source_name.Contains("GalacticPlane_All_l150"))
    {
        Source_RA = 61.117;
        Source_Dec = 52.420;
    }
    if (source_name.Contains("GalacticPlane_All_l160"))
    {
        Source_RA = 71.744;
        Source_Dec = 45.246;
    }
    if (source_name.Contains("GalacticPlane_All_l170"))
    {
        Source_RA = 79.873;
        Source_Dec = 37.315;
    }
    if (source_name.Contains("GalacticPlane_All_l180"))
    {
        Source_RA = 86.405;
        Source_Dec = 28.936;
    }
    if (source_name.Contains("GalacticPlane_All_l190"))
    {
        Source_RA = 91.940;
        Source_Dec = 20.290;
    }
    if (source_name.Contains("GalacticPlane_All_l200"))
    {
        Source_RA = 96.882;
        Source_Dec = 11.489;
    }
    if (source_name.Contains("GalacticPlane_All_l210"))
    {
        Source_RA = 101.522;
        Source_Dec = 2.609;
    }
    if (source_name.Contains("GalacticPlane_All_l220"))
    {
        Source_RA = 106.095;
        Source_Dec = -6.287;
    }
    if (source_name.Contains("SNR_G109_m01"))
    {
            Source_RA = 345.40;
            Source_Dec = 58.88;
    }
    if (source_name.Contains("SNR_G111_m02"))
    {
            Source_RA = 350.85;
            Source_Dec = 58.82;
    }
    if (source_name.Contains("PSR_J1946_p2052"))
    {
            Source_RA = 296.56;
            Source_Dec = 20.87;
    }
    if (source_name.Contains("PSR_B1937_p21"))
    {
            Source_RA = 294.91;
            Source_Dec = 21.58;
    }
    if (source_name.Contains("PSR_J1846_m0258"))
    {
            Source_RA = 281.60;
            Source_Dec = -2.98;
    }
    if (source_name.Contains("PSR_J0007_p7303"))
    {
            Source_RA = 1.76;
            Source_Dec = 73.05;
    }
    if (source_name.Contains("PSR_J0023_p0923"))
    {
            Source_RA = 5.82;
            Source_Dec = 9.39;
    }
    if (source_name.Contains("PSR_B0355_p54"))
    {
            Source_RA = 59.72;
            Source_Dec = 54.22;
    }
    if (source_name.Contains("PSR_J0030_p0451"))
    {
            Source_RA = 7.61;
            Source_Dec = 4.86;
    }
    if (source_name.Contains("PSR_J0218_p4232"))
    {
            Source_RA = 34.53;
            Source_Dec = 42.54;
    }
    if (source_name.Contains("PSR_J0205_p6449"))
    {
            Source_RA = 31.41;
            Source_Dec = 64.83;
    }
    if (source_name.Contains("PSR_J0357_p3205"))
    {
            Source_RA = 59.47;
            Source_Dec = 32.09;
    }
    if (source_name.Contains("PSR_J0633_p0632"))
    {
            //Source_RA = 98.4342083;
            //Source_Dec = 6.5430278;
            Source_RA = 98.2533333;
            Source_Dec = 5.7941667;
    }
    if (source_name.Contains("PSR_B0611_p22"))
    {
            Source_RA = 93.57;
            Source_Dec = 22.50;
    }
    if (source_name.Contains("PSR_J1841_m0345"))
    {
            Source_RA = 280.41;
            Source_Dec = -3.81;
    }
    if (source_name.Contains("PSR_J1849_m0003"))
    {
            Source_RA = 282.35;
            Source_Dec = -0.05;
    }
    if (source_name.Contains("PSR_J1856_p0245"))
    {
            Source_RA = 284.21;
            Source_Dec = 2.76;
    }
    if (source_name.Contains("PSR_J1913_p0904"))
    {
            Source_RA = 288.34;
            Source_Dec = 9.08;
    }
    if (source_name.Contains("PSR_J1938_p2213"))
    {
            Source_RA = 294.56;
            Source_Dec = 22.22;
    }
    if (source_name.Contains("PSR_J2021_p3651"))
    {
            Source_RA = 305.27;
            Source_Dec = 36.85;
    }
    if (source_name.Contains("PSR_J2021_p4026"))
    {
            Source_RA = 305.37;
            Source_Dec = 40.45;
    }
    if (source_name.Contains("PSR_J2032_p4127"))
    {
            Source_RA = 308.05;
            Source_Dec = 41.46;
    }
    if (source_name.Contains("PSR_J2238_p5903"))
    {
            Source_RA = 339.62;
            Source_Dec = 59.06;
    }
    if (source_name.Contains("PSR_J0517_p2212"))
    {
            Source_RA = 79.32;
            Source_Dec = 22.21;
    }
    if (source_name.Contains("PSR_J0751_p1807"))
    {
            Source_RA = 117.79;
            Source_Dec = 18.13;
    }
    if (source_name.Contains("PSR_J1023_p0038"))
    {
            Source_RA = 155.95;
            Source_Dec = 0.64;
    }
    if (source_name.Contains("PSR_J1024_m0719"))
    {
            Source_RA = 156.16;
            Source_Dec = -7.32;
    }
    if (source_name.Contains("PSR_J2022_p3842"))
    {
            Source_RA = 305.59;
            Source_Dec = 38.70;
    }
    if (source_name.Contains("PSR_J2229_p6114"))
    {
            Source_RA = 337.27;
            Source_Dec = 61.24;
    }
    if (source_name.Contains("PSR_B2127_p11"))
    {
            Source_RA = 322.49;
            Source_Dec = 12.17;
    }
    if (source_name.Contains("SNR_G073_p00"))
    {
            Source_RA = 303.56;
            Source_Dec = 36.20;
    }
    if (source_name.Contains("V_V725_Tau"))
    {
            Source_RA = 84.728;
            Source_Dec = 26.316;
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

    std::cout << "GetRunList source = " << source << std::endl;

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
        if (source.find("Crab_Elev60_V6") != std::string::npos) {
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
    if (source.find("CrabRHV_V6") != std::string::npos) {
        list = GetRunListFromFile("CrabRHV_V6");
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
        list = GetRunListFromFile("SS433_V6");
    }
    if (source.find("SS433_extra_V6") != std::string::npos)
    {
        list = GetRunListFromFile("SS433_extra_V6");
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
        //list = GetRunListFromFile("TychoV6");
        list_temp = GetRunListFromFile("GalacticPlane_All_l110_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l130_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
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
    if (source.find("V_V725_Tau_V6") != std::string::npos)
    {
        list = GetRunListFromFile("V_V725_Tau_V6");
    }
    if (source.find("GalacticPlane_All")!=std::string::npos && source.find("V6")!=std::string::npos)
    {
        list_temp = GetRunListFromFile("GalacticPlane_All_l30_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l50_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l70_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l90_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l110_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l130_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l150_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l170_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l190_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l210_V6");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
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
    if (source.find("PSR_J2032_p4127_Baseline_V6") != std::string::npos)
    {
        list = GetRunListFromFile("LHAASO_J2032_Baseline_V6");
    }
    if (source.find("PSR_J2032_p4127_Fall2017_V6") != std::string::npos)
    {
        list = GetRunListFromFile("LHAASO_J2032_Fall2017_V6");
    }
    if (source.find("PSR_")!=std::string::npos && source.find("_V6")!=std::string::npos)
    {
        list_temp = GetRunListFromFile(source);
        list.insert(list.end(), list_temp.begin(), list_temp.end());
    }
    if (source.find("SNR_")!=std::string::npos && source.find("_V6")!=std::string::npos)
    {
        list_temp = GetRunListFromFile(source);
        list.insert(list.end(), list_temp.begin(), list_temp.end());
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
    if (source=="OffRunsV6")
    {
        if (RHVData)
        {
            list_temp = GetRunListFromFile("CrabRHV_Imposter_V6");
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
                list_temp = GetRunListFromFile("GalacticPlane_All_l210_Imposter_V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("GalacticPlane_All_l190_Imposter_V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("GalacticPlane_All_l180_Imposter_V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("GalacticPlane_All_l150_Imposter_V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("GalacticPlane_All_l140_Imposter_V6");
                list.insert(list.end(), list_temp.begin(), list_temp.end());
                list_temp = GetRunListFromFile("GalacticPlane_All_l120_Imposter_V6");
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
    if (source.find("GalacticPlane_All")!=std::string::npos && source.find("V5")!=std::string::npos)
    {
        list_temp = GetRunListFromFile("GalacticPlane_All_l30_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l50_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l70_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l90_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l110_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l130_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l150_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l170_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l190_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l210_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
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
        list = GetRunListFromFile("SS433_V5");
    }
    if (source.find("SS433_extra_V5") != std::string::npos) {
        list = GetRunListFromFile("SS433_extra_V5");
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
        if (source.find("Crab_Elev60_V5") != std::string::npos)
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
        //list = GetRunListFromFile("TychoV5");
        list_temp = GetRunListFromFile("GalacticPlane_All_l110_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
        list_temp = GetRunListFromFile("GalacticPlane_All_l130_V5");
        list.insert(list.end(), list_temp.begin(), list_temp.end());
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
    if (source.find("PSR_J2032_p4127_Baseline_V5") != std::string::npos)
    {
        list = GetRunListFromFile("LHAASO_J2032_Baseline_V5");
    }
    if (source.find("PSR_")!=std::string::npos && source.find("_V5")!=std::string::npos)
    {
        list_temp = GetRunListFromFile(source);
        list.insert(list.end(), list_temp.begin(), list_temp.end());
    }
    if (source.find("SNR_")!=std::string::npos && source.find("_V5")!=std::string::npos)
    {
        list_temp = GetRunListFromFile(source);
        list.insert(list.end(), list_temp.begin(), list_temp.end());
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
            list_temp = GetRunListFromFile("GalacticPlane_All_l210_Imposter_V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("GalacticPlane_All_l190_Imposter_V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("GalacticPlane_All_l180_Imposter_V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("GalacticPlane_All_l140_Imposter_V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("GalacticPlane_All_l120_Imposter_V5");
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
