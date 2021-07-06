
//int N_bins_for_deconv = 20; // 8 should be the lowest bin number
int N_bins_for_deconv = 16; // 8 should be the lowest bin number
//int N_bins_for_deconv = 12; // 8 should be the lowest bin number
//int N_bins_for_deconv = 8; // 8 should be the lowest bin number

const int N_energy_bins = 6;
double Log10_alpha[N_energy_bins] = {0.,0.,0.,0.,0.,0.};
double Syst_LRR[N_energy_bins] = {0.013,0.011,0.023,0.029,0.061,0.210};
double Syst_RFoV[N_energy_bins] = {0.007,0.016,0.028,0.050,0.125,0.420};
int N_bins_for_deconv_func_E[N_energy_bins] = {N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv};
int RankTruncation[N_energy_bins] = {3,3,2,2,1,1};
double energy_bins[N_energy_bins+1] = {pow(10,2.0),pow(10,2.33),pow(10,2.66),pow(10,3.0),pow(10,3.33),pow(10,3.66),pow(10,4.0)};
double gamma_hadron_dim_ratio_w[N_energy_bins] = {1.,1.,1.,1.,1.,1.};
double gamma_hadron_dim_ratio_l[N_energy_bins] = {1.,1.,1.,1.,1.,1.};

bool EigenDecomposition = false;

double MSCW_plot_lower = -0.6;
double MSCL_plot_lower = -0.6;

char output_file_tag[50] = "tight";
double MSCW_cut_moderate = 0.5;
double MSCL_cut_moderate = 0.6;
double MSCW_cut_buffer = 0.5;
double MSCL_cut_buffer = 0.6;

bool UseTruncatedONData = false;
bool UseMinChi2 = false;

bool UsePerturbation = false;
bool TruncateNoise = false;
bool UseReplacedONData = false;
bool UseRegularization = false;

char output_file2_tag[50] = "mdm_default";
int RegularizationType = 0;
int WeightingType = 0;

int NumberOfRealEigenvectors = 4;
int NumberOfEigenvectors_Stable = 3;
double MSCW_cut_loose = 0.9;
double MSCL_cut_loose = 1.0;
double camera_theta2_cut_lower = 0.;
double camera_theta2_cut_upper = 1.;
double source_theta2_cut = 0.09;

int n_dark_samples = 1;
//const int N_energy_fine_bins = 10;
//double energy_fine_bins[N_energy_fine_bins+1] = {pow(10,2.0),pow(10,2.2),pow(10,2.4),pow(10,2.6),pow(10,2.8),pow(10,3.0),pow(10,3.2),pow(10,3.4),pow(10,3.6),pow(10,3.8),pow(10,4.0)};
const int N_energy_fine_bins = 20;
double energy_fine_bins[N_energy_fine_bins+1] = {pow(10,2.0),pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9),pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8),pow(10,3.9),pow(10,4.0)};
const int N_elev_bins = 8;
double elev_bins[N_elev_bins+1] = {45,50,55,60,65,70,75,80,85};
//const int N_MJD_bins = 16; // Aug 31st of each year
//double MJD_bins[N_MJD_bins+1] = {53613,53978,54343,54709,55074,55439,55804,56170,56535,56900,57265,57631,57996,58361,58726,59092,59457};
const int N_MJD_bins = 4; // Aug 31st of every 4 years
double MJD_bins[N_MJD_bins+1] = {53613,55074,56535,57996,59457};

double MSCW_cut_blind = 1.0;
double MSCL_cut_blind = 1.0;
double MSCW_chi2_upper = 3.;
double MSCL_chi2_upper = 3.;
double MSCW_plot_upper = 3.;
double MSCL_plot_upper = 3.;

double Skymap_size = 3.;
int Skymap_nbins = 120;

double brightness_cut = 1.0;
//double faint_brightness_cut = 5.0;
double faint_brightness_cut = 6.0;
double bright_star_radius_cut = 0.3;
