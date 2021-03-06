
int N_bins_for_deconv = 16; // 8 should be the lowest bin number

bool UseTruncatedONData = false;
bool solution_w_regularizations = false;
bool solution_w_constraints = true;
bool UseReplacedONData = false;
bool UseMinChi2 = false;

//char output_file2_tag[50] = "mdm_nominal";
//bool TruncateNoise = true;
//int NumberOfRealEigenvectors = 3;

//char output_file2_tag[50] = "mdm_1vec";
//bool TruncateNoise = true;
//int NumberOfRealEigenvectors = 1;

//char output_file2_tag[50] = "mdm_2vec";
//bool TruncateNoise = true;
//int NumberOfRealEigenvectors = 2;

char output_file2_tag[50] = "mdm_full";
bool TruncateNoise = false;
int NumberOfRealEigenvectors = 3;


int NumberOfEigenvectors_Stable = 3;
double PercentCrab = 0.;
//double MSCW_cut_moderate = 0.5;
//double MSCL_cut_moderate = 0.6;
double MSCW_cut_moderate = 0.9;
double MSCL_cut_moderate = 1.0;
double MSCW_cut_loose = 0.9;
double MSCL_cut_loose = 1.0;
double camera_theta2_cut_lower = 0.;
double camera_theta2_cut_upper = 1.;
double source_theta2_cut = 0.2;

int n_iterations = 3;
int n_dark_samples = 1;
const int N_energy_bins = 6;
double energy_bins[N_energy_bins+1] = {pow(10,2.0),pow(10,2.33),pow(10,2.66),pow(10,3.0),pow(10,3.33),pow(10,3.66),pow(10,4.0)};
const int N_energy_fine_bins = 20;
double energy_fine_bins[N_energy_fine_bins+1] = {pow(10,2.0),pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9),pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8),pow(10,3.9),pow(10,4.0)};
double gamma_flux[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double gamma_count[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double raw_gamma_count[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
const int N_elev_bins = 12;
double elev_bins[N_elev_bins+1] = {25,30,35,40,45,50,55,60,65,70,75,80,85};

//int group_size_limit[N_energy_bins] = {50000,50000,50000,50000,50000};
//int group_size_limit[N_energy_bins] = {100000,100000,100000,100000,100000};
int group_size_limit[N_energy_bins] = {50000000,50000000,50000000,50000000,50000000,50000000};

double gamma_hadron_dim_ratio_w[N_energy_bins] = {1.,1.,1.,1.,1.,1.};
double gamma_hadron_dim_ratio_l[N_energy_bins] = {1.,1.,1.,1.,1.,1.};
char output_file_tag[50] = "16bins";
int N_bins_for_deconv_func_E[N_energy_bins] = {16,16,16,16,16,16};
//char output_file_tag[50] = "24bins";
//int N_bins_for_deconv_func_E[N_energy_bins] = {24,24,24,24,24,24};
//char output_file_tag[50] = "32bins";
//int N_bins_for_deconv_func_E[N_energy_bins] = {32,32,32,32,32,32};
int N_max_ranks_func_E[N_energy_bins] = {3,3,3,3,3,3};
double MSCW_cut_lower = -1.0;
double MSCL_cut_lower = -1.0;
double MSCW_cut_blind = 1.0;
double MSCL_cut_blind = 1.0;
double MSCW_plot_upper = 3.;
double MSCL_plot_upper = 3.;
double MSCW_plot_lower = -1.;
double MSCL_plot_lower = -1.;

double Skymap_size = 4.;

double brightness_cut = 1.0;
double faint_brightness_cut = 6.0;
double bright_star_radius_cut = 0.3;
