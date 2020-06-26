
#include "MakePrediction.C"

int n_iterations_syst = 10;

vector<string> GetSourceList()
{
    vector<string> target_list;

    target_list.push_back("RGBJ0710V5");
    target_list.push_back("1ES0229V5");
    target_list.push_back("RBS0413V5");
    target_list.push_back("PG1553V5");
    target_list.push_back("Segue1V5");

    target_list.push_back("PKS1424V6");
    target_list.push_back("1ES0229V6");
    target_list.push_back("1ES0647V6");
    target_list.push_back("RBS0413V6");
    target_list.push_back("1ES1011V6");
    target_list.push_back("3C264V6");
    target_list.push_back("H1426V6");
    target_list.push_back("OJ287V6");
    target_list.push_back("Segue1V6");

    return target_list;
}
void SystematicsAnalysis(double tel_elev_lower_input, double tel_elev_upper_input, int MJD_start_cut, int MJD_end_cut, bool isON)
{

    TH1::SetDefaultSumw2();

    vector<string> target_list = GetSourceList();
    TString ONOFF_tag = "OFF";
    if (MJD_start_cut!=0 || MJD_end_cut!=0)
    {
        sprintf(mjd_cut_tag, "_MJD%dto%d", MJD_start_cut, MJD_end_cut);
    }
    int rank_variation = NumberOfEigenvectors;
    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;
    MSCW_cut_blind = MSCW_cut_moderate;
    MSCL_cut_blind = MSCL_cut_moderate;
    MSCW_plot_upper = gamma_hadron_dim_ratio*(MSCW_cut_blind-MSCW_plot_lower)+MSCW_cut_blind;
    MSCL_plot_upper = gamma_hadron_dim_ratio*(MSCL_cut_blind-MSCL_plot_lower)+MSCL_cut_blind;

    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);

    vector<vector<TH2D>> Hist_OneGroup_Dark_MSCLW;
    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
    {
        char sample_tag[50];
        sprintf(sample_tag, "%i", nth_sample);
        vector<TH2D> Hist_OneSample_Dark_MSCLW;
        for (int e=0;e<N_energy_bins;e++) 
        {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));
            Hist_OneSample_Dark_MSCLW.push_back(TH2D("Hist_OneSample_Dark_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        }
        Hist_OneGroup_Dark_MSCLW.push_back(Hist_OneSample_Dark_MSCLW);
    }
    vector<vector<TH2D>> Hist_OnData_MSCLW;
    vector<vector<vector<TH2D>>> Hist_OnBkgd_MSCLW;
    vector<vector<TH2D>> Hist_OnDark_MSCLW;
    vector<TH2D> Hist_OneGroup_Data_MSCLW;
    vector<TH2D> Hist_OneGroup_Bkgd_MSCLW;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        Hist_OneGroup_Data_MSCLW.push_back(TH2D("Hist_OneGroup_Data_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OneGroup_Bkgd_MSCLW.push_back(TH2D("Hist_OneGroup_Bkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
    }

    vector<string>* Data_runlist_name_ptr = new std::vector<string>(10);
    vector<int>* Data_runlist_number_ptr = new std::vector<int>(10);
    vector<int>* Data_runlist_MJD_ptr = new std::vector<int>(10);
    vector<double>* Data_runlist_elev_ptr = new std::vector<double>(10);
    vector<double>* Data_runlist_exposure_ptr = new std::vector<double>(10);
    vector<string>* roi_name_ptr = new std::vector<string>(10);

    char group_tag[50];
    vector<int> group_size;
    for (int e=0;e<N_energy_bins;e++) 
    {
        group_size.push_back(0);
    }

    exposure_hours = 0.;
    for (int e=0;e<N_energy_bins;e++) 
    {
        std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "energy = " << energy_bins[e] << std::endl;
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));

        int nth_group = 0;
        vector<TH2D> Hist_NthGroup_Data_MSCLW;
        vector<TH2D> Hist_NthGroup_Dark_MSCLW;
        vector<vector<TH2D>> Hist_NthGroup_Bkgd_MSCLW;
        for (int source=0;source<target_list.size();source++)
        {
            sprintf(target, "%s", target_list.at(source).c_str());
            TFile InputDataFile("../Netflix_"+TString(target)+"_"+TString(output_file_tag)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_"+ONOFF_tag+".root");
            TTree* InfoTree_ptr = nullptr;
            InfoTree_ptr = (TTree*) InputDataFile.Get("InfoTree");
            InfoTree_ptr->SetBranchAddress("Data_runlist_name",&Data_runlist_name_ptr);
            InfoTree_ptr->SetBranchAddress("Data_runlist_number",&Data_runlist_number_ptr);
            InfoTree_ptr->SetBranchAddress("Data_runlist_MJD",&Data_runlist_MJD_ptr);
            InfoTree_ptr->SetBranchAddress("Data_runlist_elev",&Data_runlist_elev_ptr);
            InfoTree_ptr->SetBranchAddress("Data_runlist_exposure",&Data_runlist_exposure_ptr);
            InfoTree_ptr->SetBranchAddress("roi_name",&roi_name_ptr);
            InfoTree_ptr->GetEntry(0);
            int FirstRun = 0;
            int LastRun = Data_runlist_name_ptr->size();

            TString hist_name;

            for (int on_run=FirstRun;on_run<LastRun;on_run++)
            {

                bool use_this_run = true;
                if (MJD_start_cut!=0 || MJD_end_cut!=0)
                {
                    if (MJD_start_cut>Data_runlist_MJD_ptr->at(on_run)) use_this_run = false;
                    if (MJD_end_cut<Data_runlist_MJD_ptr->at(on_run)) use_this_run = false;
                    if (use_this_run) std::cout << Data_runlist_MJD_ptr->at(on_run) << std::endl;
                }
                if (use_this_run)
                {
                    if (e==0)
                    {
                        exposure_hours += Data_runlist_exposure_ptr->at(on_run);
                    }
                    if (MJD_Start>Data_runlist_MJD_ptr->at(on_run)) MJD_Start = Data_runlist_MJD_ptr->at(on_run);
                    if (MJD_End<Data_runlist_MJD_ptr->at(on_run)) MJD_End = Data_runlist_MJD_ptr->at(on_run);
                    char sample_tag[50];
                    sprintf(sample_tag, "%i", on_run);

                    int binx_lower = Hist_OneGroup_Data_MSCLW.at(0).GetXaxis()->FindBin(MSCL_cut_lower);
                    binx_blind_global = Hist_OneGroup_Data_MSCLW.at(0).GetXaxis()->FindBin(MSCL_cut_blind)-1;
                    int binx_upper = Hist_OneGroup_Data_MSCLW.at(0).GetXaxis()->FindBin(1.)-1;
                    int biny_lower = Hist_OneGroup_Data_MSCLW.at(0).GetYaxis()->FindBin(MSCW_cut_lower);
                    biny_blind_global = Hist_OneGroup_Data_MSCLW.at(0).GetYaxis()->FindBin(MSCW_cut_blind)-1;
                    int biny_upper = Hist_OneGroup_Data_MSCLW.at(0).GetYaxis()->FindBin(1.)-1;

                    hist_name  = "Hist_OnData_MSCLW_V"+TString(sample_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                    Hist_OneGroup_Data_MSCLW.at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
                    {
                        char sample2_tag[50];
                        sprintf(sample2_tag, "%i", nth_sample);
                        hist_name  = "Hist_OnDark_MSCLW_R"+TString(sample_tag)+"_V"+TString(sample2_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up);
                        Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e).Add( (TH2D*)InputDataFile.Get(hist_name) );
                    }

                    group_size.at(e) += 1;
                }

                if (group_size.at(e)==100)
                {

                    sprintf(group_tag, "%i", int(nth_group));
                    Hist_NthGroup_Data_MSCLW.push_back(TH2D("Hist_NthGroup_Data_MSCLW_V"+TString(group_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                    Hist_NthGroup_Dark_MSCLW.push_back(TH2D("Hist_NthGroup_Dark_MSCLW_V"+TString(group_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));

                    char iteration_tag[50];
                    vector<TH2D> Hist_NthIteration_Bkgd_MSCLW;
                    for (int nth_iteration=0; nth_iteration<n_iterations_syst; nth_iteration++)
                    {
                        sprintf(iteration_tag, "%i", int(nth_iteration));
                        Hist_NthIteration_Bkgd_MSCLW.push_back(TH2D("Hist_NthGroup_Bkgd_MSCLW_I"+TString(iteration_tag)+"_V"+TString(group_tag)+"_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                    }
                    Hist_NthGroup_Bkgd_MSCLW.push_back(Hist_NthIteration_Bkgd_MSCLW);

                    mtx_data = fillMatrix(&Hist_OneGroup_Data_MSCLW.at(e));
                    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
                    Hist_NthGroup_Data_MSCLW.at(nth_group).Add(&Hist_OneGroup_Data_MSCLW.at(e));
                    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
                    {
                        NormalizaDarkMatrix(&Hist_OneGroup_Data_MSCLW.at(e), &Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e));
                        mtx_dark = fillMatrix(&Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e));
                        eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
                        Hist_NthGroup_Dark_MSCLW.at(nth_group).Add(&Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e),1./double(n_dark_samples));
                    }
                    for (int nth_iteration=0; nth_iteration<n_iterations_syst; nth_iteration++)
                    {
                        for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
                        {
                            NormalizaDarkMatrix(&Hist_OneGroup_Data_MSCLW.at(e), &Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e));
                            mtx_dark = fillMatrix(&Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e));
                            eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
                            LeastSquareSolutionMethod(rank_variation,nth_iteration+1);
                            TH2D Hist_Temp_Bkgd = TH2D("Hist_Temp_Bkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
                            fill2DHistogram(&Hist_Temp_Bkgd,mtx_data_bkgd);
                            Hist_OneGroup_Bkgd_MSCLW.at(e).Add(&Hist_Temp_Bkgd,1./double(n_dark_samples));
                        }
                        Hist_NthGroup_Bkgd_MSCLW.at(nth_group).at(nth_iteration).Add(&Hist_OneGroup_Bkgd_MSCLW.at(e));
                        Hist_OneGroup_Bkgd_MSCLW.at(e).Reset();
                    }
                    for (int nth_sample=0;nth_sample<n_dark_samples;nth_sample++)
                    {
                        Hist_OneGroup_Dark_MSCLW.at(nth_sample).at(e).Reset();
                    }
                    Hist_OneGroup_Data_MSCLW.at(e).Reset();
                    group_size.at(e) = 0;
                    nth_group += 1;
                }
            }
        }
        Hist_OnData_MSCLW.push_back(Hist_NthGroup_Data_MSCLW);
        Hist_OnDark_MSCLW.push_back(Hist_NthGroup_Dark_MSCLW);
        Hist_OnBkgd_MSCLW.push_back(Hist_NthGroup_Bkgd_MSCLW);

    }

    TFile OutputFile("../Netflix_Syst_"+TString(output_file_tag)+"_"+TString(output_file2_tag)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+TString(mjd_cut_tag)+"_"+ONOFF_tag+".root","recreate");

    int number_groups = Hist_OnData_MSCLW.at(0).size();
    TTree InfoTree("InfoTree","info tree");
    InfoTree.Branch("N_bins_for_deconv",&N_bins_for_deconv,"N_bins_for_deconv/I");
    InfoTree.Branch("MSCW_cut_blind",&MSCW_cut_blind,"MSCW_cut_blind/D");
    InfoTree.Branch("MSCL_cut_blind",&MSCL_cut_blind,"MSCL_cut_blind/D");
    InfoTree.Branch("exposure_hours",&exposure_hours,"exposure_hours/D");
    InfoTree.Branch("number_groups",&number_groups,"number_groups/I");
    InfoTree.Fill();
    InfoTree.Write();

    for (int e=0;e<N_energy_bins;e++)
    {
        for (int nth_group=0;nth_group<Hist_OnData_MSCLW.at(e).size();nth_group++)
        {
            Hist_OnData_MSCLW.at(e).at(nth_group).Write();
            for (int nth_iteration=0; nth_iteration<n_iterations_syst; nth_iteration++)
            {
                Hist_OnBkgd_MSCLW.at(e).at(nth_group).at(nth_iteration).Write();
            }
            Hist_OnDark_MSCLW.at(e).at(nth_group).Write();
        }
    }

}
