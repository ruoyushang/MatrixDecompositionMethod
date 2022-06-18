
const int N_energy_bins = 4;

//int N_bins_for_deconv = 16; // 8 should be the lowest bin number
//double Log10_alpha[N_energy_bins] = {-4.5,-4.5,-4.5,-4.5};
int N_bins_for_deconv = 8; // 8 should be the lowest bin number
double Log10_alpha[N_energy_bins] = {-4.5,-4.5,-4.5,-4.5};
//int N_bins_for_deconv = 4; // 8 should be the lowest bin number
//double Log10_alpha[N_energy_bins] = {-3.5,-3.5,-3.5,-3.5};
//int N_bins_for_deconv = 2; // 8 should be the lowest bin number
//double Log10_alpha[N_energy_bins] = {-3.5,-3.5,-3.5,-3.5};

//const int N_energy_bins = 6;
//double Log10_alpha[N_energy_bins] = {-5.,-3.,-3.,-3.8,-5.,-5.};
//int N_bins_for_deconv_func_E[N_energy_bins] = {N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv};
//int RankTruncation[N_energy_bins] = {3,3,2,2,1,1};
//double energy_bins[N_energy_bins+1] = {100.,251.,631.,1585.,3981.,10000.,25118.};
//double gamma_hadron_dim_ratio_w[N_energy_bins] = {1.,1.,1.,1.,1.,1.};
//double gamma_hadron_dim_ratio_l[N_energy_bins] = {1.,1.,1.,1.,1.,1.};

int N_bins_for_deconv_func_E[N_energy_bins] = {N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv};
int RankTruncation[N_energy_bins] = {3,3,2,2};
double energy_bins[N_energy_bins+1] = {100.,316.,1000.,3162.,10000.};
double gamma_hadron_dim_ratio_w[N_energy_bins] = {1.,1.,1.,1.};
double gamma_hadron_dim_ratio_l[N_energy_bins] = {1.,1.,1.,1.};

bool UseDBOnly = false;
bool RHVData = false;
bool EigenDecomposition = false;

bool UseRegularization = true;
//bool UseRegularization = false;

bool AcceptanceCorrection = false;
bool ExposureCorrection = false;

double MSCW_plot_lower = -0.6;
double MSCL_plot_lower = -0.6;

char output_file_tag[50] = "tight";
double MSCW_cut_moderate = 0.6;
double MSCL_cut_moderate = 0.6;
//double MSCW_rescale[N_energy_bins] = {0.0,0.05,0.1,0.15,0.2,0.25};
//double MSCL_rescale[N_energy_bins] = {0.0,0.05,0.1,0.15,0.2,0.25};
double MSCW_rescale[N_energy_bins] = {0.0,0.0,0.0,0.0};
double MSCL_rescale[N_energy_bins] = {0.0,0.0,0.0,0.0};

bool UseTruncatedONData = false;
bool UseMinChi2 = false;

bool UsePerturbation = false;
bool TruncateNoise = false;
bool UseReplacedONData = false;

char output_file2_tag[50] = "mdm_default";
int RegularizationType = 0;
int WeightingType = 0;

string Azim_region = "";
//string Azim_region = "North";
//string Azim_region = "South";
//string Azim_region = "West";
//string Azim_region = "East";

int NumberOfRealEigenvectors = 4;
int NumberOfEigenvectors_Stable = 3;
double MSCW_cut_loose = 0.9;
double MSCL_cut_loose = 1.0;
double camera_theta2_cut_lower = 0.;
double camera_theta2_cut_upper = 1.;
double source_theta2_cut = 0.09;

int n_dark_samples = 1;
//const int N_energy_fine_bins = 6;
//double energy_fine_bins[N_energy_fine_bins+1] = {pow(10,2.0),pow(10,2.4),pow(10,2.8),pow(10,3.2),pow(10,3.6),pow(10,4.0),pow(10,4.4)};
const int N_energy_fine_bins = 4;
double energy_fine_bins[N_energy_fine_bins+1] = {pow(10,2.0),pow(10,2.5),pow(10,3.0),pow(10,3.5),pow(10,4.0)};
const int N_elev_bins = 5;
double elev_bins[N_elev_bins+1] = {35,45,55,65,75,85};
const int N_XY_bins = 5;
int XY_bins[N_XY_bins] = {36,18,9,3,1};
const int N_integration_radii = 5;
double integration_radii[N_integration_radii] = {0.1,0.5,1.0,1.5,2.0};
//const int N_MJD_bins = 16; // Aug 31st of each year
//double MJD_bins[N_MJD_bins+1] = {53613,53978,54343,54709,55074,55439,55804,56170,56535,56900,57265,57631,57996,58361,58726,59092,59457};
const int N_MJD_bins = 4; // Aug 31st of every 4 years
double MJD_bins[N_MJD_bins+1] = {53613,55074,56535,57996,59457};
const int N_NSB_bins = 5;
double NSB_bins[N_NSB_bins+1] = {3,4,5,6,7,8};
const int N_azim_bins = 5;
double azim_bins[N_azim_bins+1] = {0,45,135,225,315,360};

double MSCW_cut_blind = 1.0;
double MSCL_cut_blind = 1.0;
double MSCW_chi2_upper = 3.;
double MSCL_chi2_upper = 3.;
double MSCW_plot_upper = 3.;
double MSCL_plot_upper = 3.;

double Skymap_size = 2.;
int Skymap_nbins = 60;
int Skymap_normalization_nbins = 1;

double brightness_cut = 5.0;
//double faint_brightness_cut = 5.0;
double faint_brightness_cut = 7.0;
double bright_star_radius_cut = 0.3;
