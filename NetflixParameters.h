
int N_bins_for_deconv = 16; // 8 should be the lowest bin number

const int N_energy_bins = 6;
double Log10_alpha[N_energy_bins] = {-6.2,-6.3,-6.0,-4.6,-4.6,-3.0};
int N_bins_for_deconv_func_E[N_energy_bins] = {N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv};
int RankTruncation[N_energy_bins] = {3,3,2,2,1,1};
double energy_bins[N_energy_bins+1] = {pow(10,2.0),pow(10,2.33),pow(10,2.66),pow(10,3.0),pow(10,3.33),pow(10,3.66),pow(10,4.0)};
double gamma_hadron_dim_ratio_w[N_energy_bins] = {1.,1.,1.,1.,1.,1.};
double gamma_hadron_dim_ratio_l[N_energy_bins] = {1.,1.,1.,1.,1.,1.};

//const int N_energy_bins = 8;
//double Log10_alpha[N_energy_bins] = {-5.9,-5.9,-5.7,-5.7,-5.7,-4.6,-4.6,-3.0};
//int N_bins_for_deconv_func_E[N_energy_bins] = {N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv};
//int RankTruncation[N_energy_bins] = {3,3,3,2,2,2,1,1};
//double energy_bins[N_energy_bins+1] = {pow(10,2.0),pow(10,2.16),pow(10,2.33),pow(10,2.50),pow(10,2.66),pow(10,3.0),pow(10,3.33),pow(10,3.66),pow(10,4.0)};
//double gamma_hadron_dim_ratio_w[N_energy_bins] = {1.,1.,1.,1.,1.,1.,1.,1.};
//double gamma_hadron_dim_ratio_l[N_energy_bins] = {1.,1.,1.,1.,1.,1.,1.,1.};

bool EigenDecomposition = false;

double MSCW_plot_lower = -0.5;
double MSCL_plot_lower = -0.5;

char output_file_tag[50] = "tight";
double MSCW_cut_moderate = 0.5;
double MSCL_cut_moderate = 0.6;
//char output_file_tag[50] = "loose";
//double MSCW_cut_moderate = 0.9;
//double MSCL_cut_moderate = 1.0;

bool UseTruncatedONData = false;
bool UseMinChi2 = false;

bool UsePerturbation = false;
bool TruncateNoise = false;
bool UseReplacedONData = false;
bool UseRegularization = false;

char output_file2_tag[50] = "mdm_default";
int RegularizationType = 0;
//char output_file2_tag[50] = "mdm_rank3";
//int RegularizationType = 1;
//char output_file2_tag[50] = "mdm_rank5";
//int RegularizationType = 2;
//char output_file2_tag[50] = "mdm_cutoff";
//int RegularizationType = 3;
//char output_file2_tag[50] = "mdm_tikhonov";
//int RegularizationType = 4;
//char output_file2_tag[50] = "mdm_cutoff_eigen";
//int RegularizationType = 5;

int NumberOfRealEigenvectors = 4;
int NumberOfEigenvectors_Stable = 3;
double PercentCrab = 0.;
double MSCW_cut_loose = 0.9;
double MSCL_cut_loose = 1.0;
double camera_theta2_cut_lower = 0.;
double camera_theta2_cut_upper = 1.;
double source_theta2_cut = 0.2;

int n_dark_samples = 1;
const int N_energy_fine_bins = 20;
double energy_fine_bins[N_energy_fine_bins+1] = {pow(10,2.0),pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9),pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8),pow(10,3.9),pow(10,4.0)};
double gamma_flux[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double gamma_count[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double raw_gamma_count[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
const int N_elev_bins = 12;
double elev_bins[N_elev_bins+1] = {25,30,35,40,45,50,55,60,65,70,75,80,85};

double MSCW_cut_blind = 1.0;
double MSCL_cut_blind = 1.0;
double MSCW_plot_upper = 3.;
double MSCL_plot_upper = 3.;

double Skymap_size = 3.;
int Skymap_nbins = 120;

double brightness_cut = 1.0;
//double faint_brightness_cut = 5.0;
double faint_brightness_cut = 6.0;
double bright_star_radius_cut = 0.3;
