/// \file
/// \ingroup tutorial_FITS
/// \notebook -draw
/// Open a FITS file and retrieve the first plane of the image array
/// as a TImage object
///
/// \macro_image
/// \macro_code
/// \macro_output
///
/// \author Claudi Martinez
//
///https://root.cern.ch/doc/master/FITS__tutorial1_8C_source.html

void FITS_tutorial1()
{
   // Here we open a FITS file that contains only the primary HDU, consisting on an image.
   // The object you will see is a snapshot of the NGC7662 nebula,
   // which was taken by the author on November 2009 in Barcelona.

   // Open primary HDU from file
   //TString dir = gROOT->GetTutorialDir();
   //std::cout << "dir = " << dir << std::endl;
   //TFITSHDU hdu(dir + "/fitsio/sample1.fits");
   TString dir = "/Users/rshang/Veritas/HAWC_result/data/";
   std::cout << "dir = " << dir << std::endl;
   //TFITSHDU hdu(dir + "skv19487154506678_1420MHz_GammaCygni.fits");
   //string output_file = "skv19487154506678_1420MHz_GammaCygni.txt";
   //double CRVAL1     = 305.557091;
   //double CRVAL2     = 40.256679;
   //TFITSHDU hdu(dir + "skv19486561528738_1420MHz_J1908.fits");
   //string output_file = "skv19486561528738_1420MHz_J1908.txt";
   //double CRVAL1     = 286.975;
   //double CRVAL2     = 6.268611;
   //TFITSHDU hdu(dir + "skv19487414377965_1420MHz_J1857.fits");
   //string output_file = "skv19487414377965_1420MHz_J1857.txt";
   //double CRVAL1     = 284.2958333333333;
   //double CRVAL2     = 2.6666666666666665;
   //TFITSHDU hdu(dir + "skv22327805333650_SwiftBat20keV_J1908.fits");
   //string output_file = "skv22327805333650_SwiftBat20keV_J1908.txt";
   //double CRVAL1     = 287.18;
   //double CRVAL2     = 6.18;
   //TFITSHDU hdu(dir + "skv22338613766121_SwiftBat14keV_J1908.fits");
   //string output_file = "skv22338613766121_SwiftBat14keV_J1908.txt";
   //double CRVAL1     = 287.18;
   //double CRVAL2     = 6.18;
   TFITSHDU hdu(dir + "skv23774123153055_fermi5_J1908.fits");
   string output_file = "skv23774123153055_fermi5_J1908.txt";
   double CRVAL1     = 287.18;
   double CRVAL2     = 6.18;

   // Dump the HDUs within the FITS file
   // and also their metadata
   // printf("Press ENTER to see summary of all data stored in the file:"); getchar();

   hdu.Print("F+");

   // Here we get the exposure time.
   printf("Exposure time = %s\n", hdu.GetKeywordValue("EXPTIME").Data());

   // Read the primary array as a matrix, selecting only layer 0.
   // This function may be useful to do image processing, e.g. custom filtering

   int NAXIS1     = 300;
   int NAXIS2     = 300;

   double skymap_size_1 = 8.;
   double skymap_size_2 = 8.;

   std::unique_ptr<TMatrixD> mat(hdu.ReadAsMatrix(0));
   //mat->Print();
   const Double_t *mtx_data = mat->GetMatrixArray();

   // Read the primary array as an image, selecting only layer 0.
   TImage * im = (TImage *)hdu.ReadAsImage(0);

   // Read the primary array as a histogram. Depending on array dimensions, the returned
   // histogram will be 1D, 2D or 3D.
   TH2 * hist = (TH2 *)hdu.ReadAsHistogram();

   double skymap_x_upper = CRVAL1+0.5*skymap_size_1;
   double skymap_x_lower = CRVAL1-0.5*skymap_size_1;
   double skymap_y_upper = CRVAL2+0.5*skymap_size_2;
   double skymap_y_lower = CRVAL2-0.5*skymap_size_2;
   TH2D * hist_skymap = new TH2D("hist_skymap","",hist->GetNbinsX(),skymap_x_lower,skymap_x_upper,hist->GetNbinsY(),skymap_y_lower,skymap_y_upper);
   for (int binx=1;binx<=hist_skymap->GetNbinsX();binx++)
   {
       for (int biny=1;biny<=hist_skymap->GetNbinsY();biny++)
       {
           int binx_reflect = hist_skymap->GetNbinsX()-binx+1;
           //hist_skymap->SetBinContent(binx_reflect,biny,hist->GetBinContent(binx,biny));
           int index = (binx-1) + (biny-1)*hist_skymap->GetNbinsX();
           hist_skymap->SetBinContent(binx_reflect,biny,max(0.,mtx_data[index]));
       }
   }
   //hist_skymap->Rebin2D(1,1);
   std::cout << "hist_skymap->GetMaximum() = " << hist_skymap->GetMaximum() << std::endl;
   std::cout << "hist_skymap->GetMinimum() = " << hist_skymap->GetMinimum() << std::endl;

   ofstream myfile;
   myfile.open (output_file);
   for (int binx=1;binx<=hist_skymap->GetNbinsX();binx++)
   {
       for (int biny=1;biny<=hist_skymap->GetNbinsY();biny++)
       {
           double RA = hist_skymap->GetXaxis()->GetBinCenter(binx);
           double Dec = hist_skymap->GetYaxis()->GetBinCenter(biny);
           double content = hist_skymap->GetBinContent(binx,biny);
           myfile << RA << " " << Dec << " " << content << "\n";
       }
   }
   myfile.close();

   //auto c = new TCanvas("c1", "FITS tutorial #1", 1400, 800);
   //c->Divide(2, 1);
   //c->cd(1);
   //im->Draw();
   //c->cd(2);
   ////gPad->SetLogz();
   //hist->Draw("COL");
   auto c = new TCanvas("c1", "FITS tutorial #1", 800, 800);
   c->cd();
   hist_skymap->Draw("COL4Z");

}
