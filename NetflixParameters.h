
const int N_energy_bins = 10;
double energy_bins[N_energy_bins+1] = {100.,200.,251.,316.,398.,501.,794.,1259.,1995.,5011.,12589.};
//double energy_bins[N_energy_bins+1] = {100.,200.,251.,316.,398.,631.,1000.,1995.,5011.,12589.};
//double energy_bins[N_energy_bins+1] = {100.,200.,316.,501.,1000.,1995.,5011.,12589.};
//double energy_bins[N_energy_bins+1] = {100.,200.,398.,794.,1585.,3162.,6310.,12589.};
const int N_energy_fine_bins = 10;
double energy_fine_bins[N_energy_fine_bins+1] = {100.,200.,251.,316.,398.,501.,794.,1259.,1995.,5011.,12589.};
//double energy_fine_bins[N_energy_fine_bins+1] = {100.,200.,251.,316.,398.,631.,1000.,1995.,5011.,12589.};
//double energy_fine_bins[N_energy_fine_bins+1] = {100.,200.,316.,501.,1000.,1995.,5011.,12589.};
//double energy_fine_bins[N_energy_fine_bins+1] = {100.,200.,398.,794.,1585.,3162.,6310.,12589.};

int N_bins_for_deconv = 12;
int matrix_rank[N_energy_bins] = {1,2,2,2,2,2,2,1,1,1};

int MatchingSelection = 0; // default
//int MatchingSelection = 1; // free elevation
//int MatchingSelection = 2; // free azimuth
//int MatchingSelection = 3; // free NSB
//int MatchingSelection = 4; // free MJD
//
double MatchRun_dElev = 0.1;
double MatchRun_dAzim = 45.;
double MatchRun_dNSB = 10.0;

double Log10_alpha_LE = 10.;
double Log10_alpha_ME = 10.;
double Log10_alpha_HE = 10.;
double Log10_alpha[N_energy_bins] = {};
double norm_scale[N_energy_bins] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//double norm_scale[N_energy_bins] = {0.9863, 0.9962, 0.9982, 0.9977, 0.9983, 0.9994, 0.9931, 0.9810, 0.9731};
//double norm_scale[N_energy_bins] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//double norm_scale[N_energy_bins] = {1.0717, 1.0040, 1.0037, 1.0036, 0.9908, 0.9958, 0.9793};

double exposure_limit = 5.; // default

double Log10_beta[N_energy_bins] =  {};
double optimiz_alpha_lower[N_energy_bins] = {};
double optimiz_alpha_upper[N_energy_bins] = {};
double optimiz_beta_lower[N_energy_bins] = {};
double optimiz_beta_upper[N_energy_bins] = {};
double gamma_hadron_dim_ratio_w = 1.;
double gamma_hadron_dim_ratio_l = 1.;
double gamma_hadron_low_end = 0.;
double MSCW_rescale[N_energy_bins] = {};
double MSCL_rescale[N_energy_bins] = {};


bool UseDL3Tree = true;
//bool RHVData = true;
bool RHVData = false;


double Skymap_large_size_x = 5.;
double Skymap_large_size_y = 5.;
bool UseGalacticCoord = false;
double Skymap_size_x = 2.5;
double Skymap_size_y = 2.5;
int Skymap_nbins_x = 100;
int Skymap_nbins_y = 100;
int Skymap_nzones_x = 1;
int Skymap_nzones_y= 1;

bool UseDBOnly = false;

bool UseRegularization = true;
//bool UseRegularization = false;

int AcceptanceCorrection = 1;  // RaDec weight 
//int AcceptanceCorrection = 2;  // XYoff weight // nominal
//int AcceptanceCorrection = 3;  // Roff weight

bool ExposureCorrection = false;

bool doRaster = false;

int nbins_fitting = 4;

char output_file_tag[50] = "tight";
double MSCW_cut_moderate = 0.5;
double MSCL_cut_moderate = 0.7;
int min_NImages = 3;
double EmissionHeight_cut = 6.;
double max_Rcore = 250.;

int MJD_start_cut = 0;
int MJD_end_cut = 0;

//double ring_radius_inner = 0.5;
double ring_radius_inner = 1.0;
//double ring_radius_inner = 1.5;

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
double camera_theta2_cut_upper = 4.0;
double source_theta_cut = 0.2;

int n_dark_samples = 1;
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

double MSCW_plot_upper = 3.;
double MSCL_plot_upper = 3.;
double MSCW_plot_lower = -0.5;
double MSCL_plot_lower = -0.5;

double MSCW_cut_blind = 1.0;
double MSCL_cut_blind = 1.0;
double MSCW_chi2_upper = 3.;
double MSCL_chi2_upper = 3.;

double brightness_cut = 6.0;
double faint_brightness_cut = 7.0;
double bright_star_radius_cut = 0.25;
//double bright_star_radius_cut = 0.;
