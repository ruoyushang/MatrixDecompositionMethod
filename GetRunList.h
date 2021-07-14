
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
        if (source.find("CrabV6") != std::string::npos) {
            list = GetRunListFromFile("CrabV6");
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
        if (source.find("RGBJ0710") != std::string::npos)
        {
            list = GetRunListFromFile("RGBJ0710");
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
        if (source.find("1ES1440V6") != std::string::npos)
        {
            list = GetRunListFromFile("1ES1440V6");
        }
        if (source.find("1ES1741V6") != std::string::npos)
        {
            list = GetRunListFromFile("1ES1741V6");
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
        if (source.find("LHAASO_J2108_V6") != std::string::npos)
        {
            list = GetRunListFromFile("LHAASO_J2108_V6");
        }
        if (source=="Everything")
        {
            list_temp = GetRunListFromFile("PKS1424V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("PG1553V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("SNR_G150p3Plus04p5_V6");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            //list_temp = GetRunListFromFile("3C273V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
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
            list_temp = GetRunListFromFile("1ES0647V6");
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

            //list_temp = GetRunListFromFile("NGC1275V6");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
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
        if (source.find("Mrk421V5") != std::string::npos)
        {
            list = GetRunListFromFile("Mrk421V5");
        }
        if (source.find("CTA1V5") != std::string::npos) {
            list = GetRunListFromFile("CTA1V5");
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
        if (source.find("CrabV5") != std::string::npos)
        {
            list = GetRunListFromFile("CrabV5");
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
        if (source.find("RGBJ0710V5") != std::string::npos)
        {
            list = GetRunListFromFile("RGBJ0710V5");
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
        if (source=="EverythingV5")
        {
            list_temp = GetRunListFromFile("BLLacV5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("RGBJ0710V5");
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
            //list_temp = GetRunListFromFile("3C273V5");
            //list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("1ES0502V5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
            list_temp = GetRunListFromFile("DracoV5");
            list.insert(list.end(), list_temp.begin(), list_temp.end());
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
        if (source=="EverythingV4")
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
