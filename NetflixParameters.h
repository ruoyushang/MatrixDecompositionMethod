
const int N_energy_bins = 6;

//int N_bins_for_deconv = 4;
//int N_bins_for_deconv = 6;
//int N_bins_for_deconv = 8; 
int N_bins_for_deconv = 12; // nominal
//int N_bins_for_deconv = 16; 
//int N_bins_for_deconv = 20; 

//double elbow_ratio[N_energy_bins] = {1.0,1.0,1.0,1.0,1.0,1.0}; 
//double elbow_ratio[N_energy_bins] = {2.0,2.0,2.0,2.0,2.0,2.0};
double elbow_ratio[N_energy_bins] = {4.0,4.0,4.0,4.0,4.0,4.0}; // nominal
//double elbow_ratio[N_energy_bins] = {8.0,8.0,8.0,8.0,8.0,8.0};

int MatchingSelection = 0; // default
//int MatchingSelection = 1; // free elevation
//int MatchingSelection = 2; // free azimuth
//int MatchingSelection = 3; // free NSB
//int MatchingSelection = 4; // free MJD

double Log10_alpha_single = 1.0;
double Log10_alpha[N_energy_bins] = {Log10_alpha_single,Log10_alpha_single,Log10_alpha_single,Log10_alpha_single,Log10_alpha_single,Log10_alpha_single};


double Log10_beta[N_energy_bins] =  {0.0,0.0,0.0,0.0,0.0,0.0};
double optimiz_alpha_lower[N_energy_bins] = {-1.5,-1.5,-1.5,-1.5,-1.5,-1.5};
double optimiz_alpha_upper[N_energy_bins] = {1.5,1.5,1.5,1.5,1.5,1.5};
double optimiz_beta_lower[N_energy_bins] = {-1.5,-1.5,-1.5,-1.5,-1.5,-1.5};
double optimiz_beta_upper[N_energy_bins] = {1.5,1.5,1.5,1.5,1.5,1.5};
double energy_bins[N_energy_bins+1] = {200.,398.,794.,1585.,3162.,6310.,12589.};
double gamma_hadron_dim_ratio_w = 1.;
double gamma_hadron_dim_ratio_l = 1.;
double gamma_hadron_low_end = 0.;
double MSCW_rescale[N_energy_bins] = {0.0,0.0,0.0,0.0,0.0,0.0};
double MSCL_rescale[N_energy_bins] = {0.0,0.0,0.0,0.0,0.0,0.0};
const int N_energy_fine_bins = 6;
double energy_fine_bins[N_energy_fine_bins+1] = {200.,398.,794.,1585.,3162.,6310.,12589.};


bool UseDL3Tree = true;
//bool RHVData = true;
bool RHVData = false;


//bool UseGalacticCoord = true;
//double Skymap_size_x = 6.;
//int Skymap_nbins_x = 60;
//double Skymap_size_y = 6.;
//int Skymap_nbins_y = 60;
//int Skymap_nzones_x = 3;
//int Skymap_nzones_y= 3;
bool UseGalacticCoord = false;
double Skymap_size_x = 2.;
int Skymap_nbins_x = 45;
double Skymap_size_y = 2.;
int Skymap_nbins_y = 45;
int Skymap_nzones_x = 1;
int Skymap_nzones_y= 1;
//bool UseGalacticCoord = false;
//double Skymap_size_x = 3.;
//int Skymap_nbins_x = 68;
//double Skymap_size_y = 3.;
//int Skymap_nbins_y = 68;
//int Skymap_nzones_x = 2;
//int Skymap_nzones_y= 2;

bool UseDBOnly = false;
//bool UseDBOnly = true;
double exposure_limit = 5.; // default
//double exposure_limit = 1000.;

bool UseRegularization = true;
//bool UseRegularization = false;

//int AcceptanceCorrection = 1;  // RaDec weight 
int AcceptanceCorrection = 2;  // XYoff weight // nominal
//int AcceptanceCorrection = 3;  // Roff weight

bool ExposureCorrection = false;

int nbins_fitting = 4;

char output_file_tag[50] = "tight";
double MSCW_cut_moderate = 0.6;
double MSCL_cut_moderate = 0.6;
//double EmissionHeight_cut = 0.;
double EmissionHeight_cut = 6.;

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
double camera_theta2_cut_upper = 3.;
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
double bright_star_radius_cut = 0.3;
