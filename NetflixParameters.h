
int N_bins_for_deconv = 16; // 8 should be the lowest bin number
//int N_bins_for_deconv = 8; // 8 should be the lowest bin number
//int N_bins_for_deconv = 4; // 8 should be the lowest bin number
//int N_bins_for_deconv = 2; // 8 should be the lowest bin number

const int N_energy_bins = 7;
double Log10_alpha[N_energy_bins] = {-5.,-3.,-3.,-3.8,-5.,-5.,-5.};
int N_bins_for_deconv_func_E[N_energy_bins] = {N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv,N_bins_for_deconv};
int RankTruncation[N_energy_bins] = {3,3,2,2,1,1,1};
double energy_bins[N_energy_bins+1] = {100.,214.,457.,1000.,2140.,4570.,10000.,25000.};
double gamma_hadron_dim_ratio_w[N_energy_bins] = {1.,1.,1.,1.,1.,1.,1.};
double gamma_hadron_dim_ratio_l[N_energy_bins] = {1.,1.,1.,1.,1.,1.,1.};

bool EigenDecomposition = false;

bool UseRegularization = true;
//bool UseRegularization = false;

//bool AcceptanceCorrection = true;
bool AcceptanceCorrection = false;
//bool ExposureCorrection = true;
bool ExposureCorrection = false;

double MSCW_plot_lower = -0.6;
double MSCL_plot_lower = -0.6;

char output_file_tag[50] = "tight";
double MSCW_cut_moderate = 0.6;
double MSCL_cut_moderate = 0.6;
double MSCW_rescale[N_energy_bins] = {1.,1.,1.1,1.2,1.5,2.,2.4};
double MSCL_rescale[N_energy_bins] = {1.,1.,1.1,1.2,1.5,2.,2.4};

bool UseTruncatedONData = false;
bool UseMinChi2 = false;

bool UsePerturbation = false;
bool TruncateNoise = false;
bool UseReplacedONData = false;

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
const int N_energy_fine_bins = 12;
double energy_fine_bins[N_energy_fine_bins+1] = {pow(10,2.0),pow(10,2.2),pow(10,2.4),pow(10,2.6),pow(10,2.8),pow(10,3.0),pow(10,3.2),pow(10,3.4),pow(10,3.6),pow(10,3.8),pow(10,4.0),pow(10,4.2),pow(10,4.4)};
const int N_elev_bins = 5;
double elev_bins[N_elev_bins+1] = {40,50,60,70,80,90};
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

double Skymap_size = 3.;
int Skymap_nbins = 120;

double brightness_cut = 1.0;
//double faint_brightness_cut = 5.0;
double faint_brightness_cut = 6.0;
double bright_star_radius_cut = 0.3;
