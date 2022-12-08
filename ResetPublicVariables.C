

void ResetPublicVariables()
{

    string MY_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    std::cout << "MY_OUTPUT = " << MY_OUTPUT << std::endl;

    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_db_search")
    {
        UseDBOnly = true;
        exposure_limit = 1000.;
    }

    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_4x4")
    {
        N_bins_for_deconv = 4;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_6x6")
    {
        N_bins_for_deconv = 6;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_8x8")
    {
        N_bins_for_deconv = 8;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_12x12")
    {
        N_bins_for_deconv = 12;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_16x16")
    {
        N_bins_for_deconv = 16;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_20x20")
    {
        N_bins_for_deconv = 20;
    }

    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_elbow1p0")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            elbow_ratio[eb] = 1.0; 
        }
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_elbow2p0")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            elbow_ratio[eb] = 2.0; 
        }
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_elbow4p0")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            elbow_ratio[eb] = 4.0; 
        }
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_elbow8p0")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            elbow_ratio[eb] = 8.0; 
        }
    }

    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_FreeElev")
    {
        MatchingSelection = 1;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_FreeAzim")
    {
        MatchingSelection = 2;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_FreeNSB")
    {
        MatchingSelection = 3;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_FreeMJD")
    {
        MatchingSelection = 4;
    }

    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_LogAlpha_m1p5")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            Log10_alpha[eb] = -1.5; 
        }
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_LogAlpha_m1p0")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            Log10_alpha[eb] = -1.0; 
        }
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_LogAlpha_m0p5")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            Log10_alpha[eb] = -0.5; 
        }
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_LogAlpha_p0p0")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            Log10_alpha[eb] = 0.0; 
        }
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_LogAlpha_p0p5")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            Log10_alpha[eb] = 0.5; 
        }
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_LogAlpha_p1p0")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            Log10_alpha[eb] = 1.0; 
        }
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_LogAlpha_p1p5")
    {
        for (int eb=0;eb<N_energy_bins;eb++)
        {
            Log10_alpha[eb] = 1.5; 
        }
    }

}
