
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
        if (source=="ComaV4") {
            list = GetRunListFromFile("ComaV4");
        }
        if (source=="ComaV6") {
            list = GetRunListFromFile("ComaV6");
        }
        if (source=="Ton599") {
            list = GetRunListFromFile("Ton599");
        }
        if (source=="3C264") {
            list = GetRunListFromFile("3C264");
        }
        if (source=="PKS1424") {
            list = GetRunListFromFile("PKS1424");
        }
        if (source=="IC443HotSpotV6") {
            list = GetRunListFromFile("IC443HotSpotV6");
        }
        if (source=="H1426") {
            list = GetRunListFromFile("H1426");
        }
        if (source=="Crab") {
            list = GetRunListFromFile("Crab");
        }
        if (source=="S3_1227_V6")
        {
            list = GetRunListFromFile("S3_1227_V6");
        }
        if (source=="1ES0229")
        {
            list = GetRunListFromFile("1ES0229");
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
        if (source=="M82")
        {
            list = GetRunListFromFile("M82");
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
        if (source=="MGRO_J1908_V6_new")
        {
            list = GetRunListFromFile("MGRO_J1908_V6_new");
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
        if (source=="2HWC_J1953V6_new")
        {
            list = GetRunListFromFile("2HWC_J1953V6_new");
        }
        if (source=="2HWC_J1930V6")
        {
            list = GetRunListFromFile("2HWC_J1930V6");
        }
        if (source=="Everything")
        {
            list_temp = GetRunListFromFile("IC443HotSpotV6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("IC443_MatchedDark_V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("WComaeV6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
        }
        if (source=="CTA1V5") {
            list = GetRunListFromFile("CTA1V5");
        }
        if (source=="MGRO_J1908_V5") {
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
        if (source=="EverythingV5")
        {
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
        }
        if (source=="WComaeV4")
        {
            list = GetRunListFromFile("WComaeV4");
        }
        if (source=="MGRO_J1908_V4") {
            list = GetRunListFromFile("MGRO_J1908_V4");
        }
        if (source=="IC443HotSpotV4")
        {
            list = GetRunListFromFile("IC443HotSpotV4");
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
        }

        return list;
}
