

void ResetPublicVariables(TString target_name)
{

    //if (target_name.Contains("Crab"))
    //{
    //    doRaster = true;
    //}

    string MY_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    std::cout << "MY_OUTPUT = " << MY_OUTPUT << std::endl;

    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_db_search")
    {
        UseDBOnly = true;
        exposure_limit = 1000.;
    }

    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_tight")
    {
        MSCW_cut_moderate = 0.4;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_medium")
    {
        MSCW_cut_moderate = 0.5;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_loose")
    {
        MSCW_cut_moderate = 0.6;
        MSCL_cut_moderate = 0.8;
    }

    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_1hrs")
    {
        exposure_limit = 1.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_2hrs")
    {
        exposure_limit = 2.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_5hrs")
    {
        exposure_limit = 5.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_10hrs")
    {
        exposure_limit = 10.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_20hrs")
    {
        exposure_limit = 20.;
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

    //if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_elbow1p2")
    //{
    //    for (int eb=0;eb<N_energy_bins;eb++)
    //    {
    //        elbow_ratio[eb] = 1.2; 
    //    }
    //}
    //if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_elbow1p5")
    //{
    //    for (int eb=0;eb<N_energy_bins;eb++)
    //    {
    //        elbow_ratio[eb] = 1.5; 
    //    }
    //}
    //if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_elbow2p0")
    //{
    //    for (int eb=0;eb<N_energy_bins;eb++)
    //    {
    //        elbow_ratio[eb] = 2.0; 
    //    }
    //}
    //if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_elbow2p5")
    //{
    //    for (int eb=0;eb<N_energy_bins;eb++)
    //    {
    //        elbow_ratio[eb] = 2.5; 
    //    }
    //}
    //if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_elbow3p0")
    //{
    //    for (int eb=0;eb<N_energy_bins;eb++)
    //    {
    //        elbow_ratio[eb] = 3.0; 
    //    }
    //}

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

    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dEl_0p2")
    {
        MatchRun_dElev = 0.2;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_05")
    {
        MatchRun_dAzim = 5.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_10")
    {
        MatchRun_dAzim = 10.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_20")
    {
        MatchRun_dAzim = 20.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_40")
    {
        MatchRun_dAzim = 40.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_80")
    {
        MatchRun_dAzim = 80.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_180")
    {
        MatchRun_dAzim = 1000.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_150")
    {
        MatchRun_dAzim = 150.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_120")
    {
        MatchRun_dAzim = 120.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_90")
    {
        MatchRun_dAzim = 90.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_60")
    {
        MatchRun_dAzim = 60.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dAz_30")
    {
        MatchRun_dAzim = 30.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dNSB_0p5")
    {
        MatchRun_dNSB = 0.5;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dNSB_1p0")
    {
        MatchRun_dNSB = 1.0;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dNSB_1p5")
    {
        MatchRun_dNSB = 1.5;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dNSB_2p0")
    {
        MatchRun_dNSB = 2.0;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dNSB_2p5")
    {
        MatchRun_dNSB = 2.5;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dNSB_3p0")
    {
        MatchRun_dNSB = 3.0;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dNSB_3p5")
    {
        MatchRun_dNSB = 3.5;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_dNSB_4p0")
    {
        MatchRun_dNSB = 4.0;
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
