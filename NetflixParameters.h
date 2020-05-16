
char output_file_tag[50] = "16bins";
int NumberOfEigenvectors = 3;
double PercentCrab = 0.;
double tel_elev_lower_input = 45.;
double tel_elev_upper_input = 85.;
double MSCW_cut_moderate = 0.35;
double MSCL_cut_moderate = 0.5;
double MSCW_cut_loose = 0.5;
double MSCL_cut_loose = 0.7;
double camera_theta2_cut = 9.;
int N_bins_for_deconv = 16; // 8 should be the lowest bin number

int n_control_samples = 6;
const int N_energy_bins = 1;
double energy_bins[N_energy_bins+1] = {pow(10,2.0),pow(10,4.0)};
const int N_energy_fine_bins = 20;
double energy_fine_bins[N_energy_fine_bins+1] = {pow(10,2.0),pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9),pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8),pow(10,3.9),pow(10,4.0)};
double gamma_flux[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double gamma_count[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double raw_gamma_count[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

double gamma_hadron_dim_ratio = 1.;
double MSCW_cut_lower = -1.0;
double MSCL_cut_lower = -1.0;
double MSCW_cut_blind = 1.0;
double MSCL_cut_blind = 1.0;
double MSCW_plot_upper = 3.;
double MSCL_plot_upper = 3.;
double MSCW_plot_lower = -1.;
double MSCL_plot_lower = -1.;
