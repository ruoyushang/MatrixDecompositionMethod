
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

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

    ifstream myfile ("/home/rshang/EventDisplay/NewBkgMethodExtendedSource/RunList_"+source+".txt");
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
        if (source=="Segue1V6")
        {
            list = GetRunListFromFile("Segue1V6");
        }
        if (source=="CrabV4")
        {
            list = GetRunListFromFile("CrabV4");
        }
        if (source=="Ton599") {
            list = GetRunListFromFile("Ton599");
        }
        if (source=="3C264V6") {
            list = GetRunListFromFile("3C264V6");
        }
        if (source=="PKS1424V6") {
            list = GetRunListFromFile("PKS1424V6");
        }
        if (source=="IC443HotSpotV6") {
            list = GetRunListFromFile("IC443HotSpotV6");
        }
        if (source=="H1426V6") {
            list = GetRunListFromFile("H1426V6");
        }
        if (source=="CrabV6") {
            list = GetRunListFromFile("CrabV6");
        }
        if (source=="SNRG150p3_V6") {
            list = GetRunListFromFile("SNRG150p3_V6");
        }
        if (source=="S3_1227_V6")
        {
            list = GetRunListFromFile("S3_1227_V6");
        }
        if (source=="1ES0229V6")
        {
            list = GetRunListFromFile("1ES0229V6");
        }
        if (source=="SgrAV6")
        {
            list = GetRunListFromFile("SgrAV6");
        }
        if (source=="GemingaV6")
        {
            list = GetRunListFromFile("GemingaV6");
        }
        if (source=="CasA")
        {
            list = GetRunListFromFile("CasA");
        }
        if (source=="Photon")
        {
            list = GetRunListFromFile("Photon");
        }
        if (source=="Proton_NSB200")
        {
            list = GetRunListFromFile("Proton_NSB200");
        }
        if (source=="G079")
        {
            list = GetRunListFromFile("G079");
        }
        if (source=="RGBJ0710")
        {
            list = GetRunListFromFile("RGBJ0710");
        }
        if (source=="Mrk421")
        {
            list = GetRunListFromFile("Mrk421");
        }
        if (source=="WComaeV6")
        {
            list = GetRunListFromFile("WComaeV6");
        }
        if (source=="1ES1218V6")
        {
            list = GetRunListFromFile("1ES1218V6");
        }
        if (source=="OJ287V6")
        {
            list = GetRunListFromFile("OJ287V6");
        }
        if (source=="1ES1011V6")
        {
            list = GetRunListFromFile("1ES1011V6");
        }
        if (source=="NGC1275V6")
        {
            list = GetRunListFromFile("NGC1275V6");
        }
        if (source=="1ES0647V6")
        {
            list = GetRunListFromFile("1ES0647V6");
        }
        if (source=="1ES1440V6")
        {
            list = GetRunListFromFile("1ES1440V6");
        }
        if (source=="1ES1741V6")
        {
            list = GetRunListFromFile("1ES1741V6");
        }
        if (source=="MGRO_J1908_V6")
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source=="SS433_V6")
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source=="SS433Half1_V6")
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source=="SS433Half2_V6")
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source=="SS433Quad1_V6")
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source=="SS433Quad2_V6")
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source=="SS433Quad3_V6")
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source=="SS433Quad4_V6")
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source=="MAGIC_J1857_V6")
        {
            list = GetRunListFromFile("MGRO_J1908_V6");
        }
        if (source=="MGRO_J2031_V6")
        {
            list = GetRunListFromFile("MGRO_J2031_V6");
        }
        if (source=="CygnusV6")
        {
            list = GetRunListFromFile("CygnusV6");
        }
        if (source=="RBS0413V6")
        {
            list = GetRunListFromFile("RBS0413V6");
        }
        if (source=="PG1553V6")
        {
            list = GetRunListFromFile("PG1553V6");
        }
        if (source=="PKS1441V6")
        {
            list = GetRunListFromFile("PKS1441V6");
        }
        if (source=="MS1221V6")
        {
            list = GetRunListFromFile("MS1221V6");
        }
        if (source=="TychoV6")
        {
            list = GetRunListFromFile("TychoV6");
        }
        if (source=="2HWC_J1953V6")
        {
            list = GetRunListFromFile("2HWC_J1953V6");
        }
        if (source=="2HWC_J1930V6")
        {
            list = GetRunListFromFile("2HWC_J1930V6");
        }
        if (source=="ComaV6") {
            list = GetRunListFromFile("ComaV6");
        }
        if (source=="BoomerangV6")
        {
            list = GetRunListFromFile("BoomerangV6");
        }
        if (source=="BLLacV6")
        {
            list = GetRunListFromFile("BLLacV6");
        }
        if (source=="GammaCygniV6")
        {
            list = GetRunListFromFile("GammaCygniV6");
        }
        if (source=="Everything")
        {
            list_temp = GetRunListFromFile("BLLacV6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("NGC1275V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("Segue1V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("1ES0647V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("RBS0413V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("1ES1011V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("3C264V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("PKS1424V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("H1426V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("1ES0229V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("OJ287V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("IC443HotSpotV6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("IC443_MatchedDark_V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("WComaeV6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("WComae_MatchedDark_V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("MGRO_J1908_MatchedDark_V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("SgrA_MatchedDark_V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("Cygnus_MatchedDark_V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("Geminga_MatchedDark_V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
        }
        if (source=="Mrk421V5")
        {
            list = GetRunListFromFile("Mrk421V5");
        }
        if (source=="CTA1V5") {
            list = GetRunListFromFile("CTA1V5");
        }
        if (source=="MGRO_J1908_V5") {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source=="SS433_V5") {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source=="SS433Half1_V5") {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source=="SS433Half2_V5") {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source=="SS433Quad1_V5") {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source=="SS433Quad2_V5") {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source=="SS433Quad3_V5") {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source=="SS433Quad4_V5") {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source=="MAGIC_J1857_V5") {
            list = GetRunListFromFile("MGRO_J1908_V5");
        }
        if (source=="IC443HotSpotV5")
        {
            list = GetRunListFromFile("IC443HotSpotV5");
        }
        if (source=="Segue1V5")
        {
            list = GetRunListFromFile("Segue1V5");
        }
        if (source=="GemingaV5")
        {
            list = GetRunListFromFile("GemingaV5");
        }
        if (source=="CrabV5")
        {
            list = GetRunListFromFile("CrabV5");
        }
        if (source=="WComaeV5")
        {
            list = GetRunListFromFile("WComaeV5");
        }
        if (source=="CygnusV5")
        {
            list = GetRunListFromFile("CygnusV5");
        }
        if (source=="MGRO_J2031_V5")
        {
            list = GetRunListFromFile("MGRO_J2031_V5");
        }
        if (source=="PG1553V5")
        {
            list = GetRunListFromFile("PG1553V5");
        }
        if (source=="RBS0413V5")
        {
            list = GetRunListFromFile("RBS0413V5");
        }
        if (source=="1ES0229V5")
        {
            list = GetRunListFromFile("1ES0229V5");
        }
        if (source=="RGBJ0710V5")
        {
            list = GetRunListFromFile("RGBJ0710V5");
        }
        if (source=="M82V5")
        {
            list = GetRunListFromFile("M82V5");
        }
        if (source=="BoomerangV5")
        {
            list = GetRunListFromFile("BoomerangV5");
        }
        if (source=="BLLacV5")
        {
            list = GetRunListFromFile("BLLacV5");
        }
        if (source=="GammaCygniV5")
        {
            list = GetRunListFromFile("GammaCygniV5");
        }
        if (source=="EverythingV5")
        {
            list_temp = GetRunListFromFile("BLLacV5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("RGBJ0710V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("1ES0229V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("RBS0413V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("PG1553V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("Segue1V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("IC443HotSpotV5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("WComaeV5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("WComae_MatchedDark_V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("MGRO_J1908_MatchedDark_V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("IC443_MatchedDark_V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("Cygnus_MatchedDark_V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("ComaV5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("M82V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
        }
        if (source=="WComaeV4")
        {
            list = GetRunListFromFile("WComaeV4");
        }
        if (source=="MGRO_J1908_V4") {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source=="SS433_V4") {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source=="SS433Half1_V4") {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source=="SS433Half2_V4") {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source=="SS433Quad1_V4") {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source=="SS433Quad2_V4") {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source=="SS433Quad3_V4") {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source=="SS433Quad4_V4") {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source=="MAGIC_J1857_V4") {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source=="MGRO_J2031_V4") {
            list = GetRunListFromFile("MGRO_J2031_V4");
        }
        if (source=="IC443HotSpotV4")
        {
            list = GetRunListFromFile("IC443HotSpotV4");
        }
        if (source=="ComaV4") {
            list = GetRunListFromFile("ComaV4");
        }
        if (source=="M82V4")
        {
            list = GetRunListFromFile("M82V4");
        }
        if (source=="BoomerangV4")
        {
            list = GetRunListFromFile("BoomerangV4");
        }
        if (source=="GammaCygniV4")
        {
            list = GetRunListFromFile("GammaCygniV4");
        }
        if (source=="EverythingV4")
        {
            list_temp = GetRunListFromFile("IC443HotSpotV4");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("WComaeV4");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
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
