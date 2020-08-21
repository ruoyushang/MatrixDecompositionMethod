
//char output_file_tag[50] = "8bins";
//int N_bins_for_deconv = 8; // 8 should be the lowest bin number
char output_file_tag[50] = "16bins";
int N_bins_for_deconv = 16; // 8 should be the lowest bin number

char output_file2_tag[50] = "constrained";
bool solution_w_constraints = true;
//char output_file2_tag[50] = "unconstrained";
//bool solution_w_constraints = false;

int NumberOfEigenvectors = 3;
int NumberOfEigenvectors_Stable = 3;
double PercentCrab = 0.;
double MSCW_cut_moderate = 0.5;
double MSCL_cut_moderate = 0.6;
double MSCW_cut_loose = 0.9;
double MSCL_cut_loose = 1.0;
double camera_theta2_cut_lower = 0.;
double camera_theta2_cut_upper = 1.;
double source_theta2_cut = 0.2;

bool DoSequential = true;
int n_iterations = 5;
int n_dark_samples = 3;
const int N_energy_bins = 5;
double energy_bins[N_energy_bins+1] = {pow(10,2.0),pow(10,2.3),pow(10,2.6),pow(10,3.0),pow(10,3.4),pow(10,4.0)};
const int N_energy_fine_bins = 20;
double energy_fine_bins[N_energy_fine_bins+1] = {pow(10,2.0),pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9),pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8),pow(10,3.9),pow(10,4.0)};
double gamma_flux[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double gamma_count[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double raw_gamma_count[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
const int N_elev_bins = 12;
double elev_bins[N_elev_bins+1] = {25,30,35,40,45,50,55,60,65,70,75,80,85};

//int group_size_limit[N_energy_bins] = {50000,50000,50000,50000,50000};
int group_size_limit[N_energy_bins] = {50000000,50000000,50000000,50000000,50000000};

double gamma_hadron_dim_ratio_w[N_energy_bins] = {1.,1.,1.,3.,7.};
double gamma_hadron_dim_ratio_l[N_energy_bins] = {1.,1.,1.,1.,1.};
double MSCW_cut_lower = -1.0;
double MSCL_cut_lower = -1.0;
double MSCW_cut_blind = 1.0;
double MSCL_cut_blind = 1.0;
double MSCW_plot_upper = 3.;
double MSCL_plot_upper = 3.;
double MSCW_plot_lower = -1.;
double MSCL_plot_lower = -1.;

double Skymap_size = 4.;

double brightness_cut = 5.0;
double faint_brightness_cut = 6.0;
double bright_star_radius_cut = 0.3;
