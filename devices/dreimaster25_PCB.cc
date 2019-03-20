// Simon Spannagel (DESY) January 2016

#include "TCanvas.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TF1.h"

#include "assembly.h"
#include "propagate.h"
#include "materials.h"
#include "constants.h"
#include "log.h"
#include "math.h"

using namespace std;
using namespace gblsim;
using namespace unilog;
#define energies 14

int main(int argc, char* argv[]) {

  /*
   * dreimaster resolution simulation with 25x100 um pixel sensors and r4s chip at the DESY TB21 beam line
   * Three planes with 20mm spacing, intrinsic sensor resolution FIX
   */

  
  Log::ReportingLevel() = Log::FromString("DEBUG");

  for (int i = 1; i < argc; i++) {
    // Setting verbosity:
    if (std::string(argv[i]) == "-v") { 
      Log::ReportingLevel() = Log::FromString(std::string(argv[++i]));
      continue;
    } 
  }

  TFile * out = TFile::Open("/home/zoiirene/Programming/resolution/Output/dreimaster25PCB-resolution.root","RECREATE");
  gDirectory->pwd();

  TCanvas *c1 = new TCanvas("c1","resolution",700,700);
  TProfile *resolution = new TProfile("resolution"," ",100,0,7);

  //----------------------------------------------------------------------------
  // Preparation of the dreimaster and beam properties:

  // dreimaster planes consist of 150um silicon + 700um chip plus 1.56mm PCB plus 50um copper:
  double sensor = 850e-3 / X0_Si + 1.56/X0_PCB + 50e-3/X0_Cu;
  // or 150um silicon + 700um chip if there is cut out:
  //double sensor = 850e-3 / X0_Si; // + 1.56/X0_PCB + 50e-3/X0_Cu;
  cout << " sensor thcikness / x0 " << sensor << endl;
  // The intrinsic resolution one hit resolution
  double RES = 2.68e-3 ; // / sqrt(3./2.); //updated from M-Scan 2832-2840 
  
  // sensor DUT(sensor) sensor
  //    |       |       |
  //    |       |       |
  //    |<----->|       |
  //       DIST      
  
  // Distance between 3M planes in mm:
  double DIST = 20;

  double BEAM[energies];
  double INVBEAM2[energies];
  double INVBEAM2error[energies];
  double RESOLUTION[energies];
  double resError = 0.03; 
  double RESOLUTION2[energies];
  double RESOLUTION2error[energies];
  //----------------------------------------------------------------------------
  // Build the trajectory through the telescope device:

  // Build a vector of all telescope planes:
  std::vector<plane> dreimaster;
  double position = 0;
  
  // Upstream 3M plane:
  for(int i = 0; i < 1; i++) {
    dreimaster.push_back(plane(position,sensor,true,RES));
    }

  // Downstream 3M plane:
  //  position = 2*DIST + 2*DUT_DIST;
  position = DIST*2;
  for(int i = 0; i < 1; i++) {
    dreimaster.push_back(plane(position,sensor,true,RES));
    }


    // Prepare the DUT (no measurement, just scatterer
    plane dut(DIST, sensor, false); 

    // Duplicate the planes vector and add the current DUT:
    std::vector<plane> planes = dreimaster;
    planes.push_back(dut);
    for (int i = 0; i< energies ; i++)
      {
	BEAM[i] = 1.2 +0.4*i; //GeV
	if(i==energies - 1) BEAM[i] = 100; //GeV
	cout << BEAM[i] << endl;
	INVBEAM2[i] = 1./(BEAM[i]*BEAM[i]);
	INVBEAM2error[i] = 0.;
	// Build the telescope:
	telescope mytel(planes, BEAM[i]);
	
	// Get the resolution at plane-vector position (x):
	LOG(logRESULT) << "Track resolution at BEAM energy " << BEAM[i] << "GeV: " << mytel.getResolution(1);
	RESOLUTION[i] = mytel.getResolution(1);
	RESOLUTION2[i] = mytel.getResolution(1)*mytel.getResolution(1);
	RESOLUTION2error[i] = 2*mytel.getResolution(1)*resError;
	resolution->Fill(BEAM[i],RESOLUTION[i],1);
	LOG(logRESULT) << "Track resolution at 1/(BEAM energy)^2 " << 1./(BEAM[i]*BEAM[i]) << "GeV: " << mytel.getResolution(1); //*mytel.getResolution(1);
	//resolutionsquare->Fill(INVBEAM2[i],mytel.getResolution(1)*mytel.getResolution(1),1);
      }


    cout << " plotting canvas 1" <<endl;
    c1->cd();
    resolution->SetTitle("3M Track Resolution  at DUT;beam momentum scan ;resolution at DUT #left[#mum#right]");
    resolution->GetYaxis()->SetRangeUser(0.,16.);
    resolution->SetMarkerStyle(20);
    resolution->SetLineColor(kRed+1);
    resolution->SetLineWidth(2);
    resolution->SetMarkerColor(kRed+1);
    resolution->SetMinimum(0.);  
    resolution->Draw();
    c1->Write();

    cout << " plotting canvas 2" <<endl;


    TCanvas *c2 = new TCanvas("c2","resolutionsquare",700,700);
    TGraph *resolutionsquare = new TGraph(energies,INVBEAM2,RESOLUTION2);//"resolutionsquare"," ",100,0.,1.);
    
    c2->cd();
    resolutionsquare->SetTitle("3M Track Resolution^{2}  at DUT;beam momentum scan ;resolution at DUT #left[#mum^{2}#right]");
    resolutionsquare->GetYaxis()->SetRangeUser(0.,300.);
    resolutionsquare->GetXaxis()->SetRangeUser(0.,1.);
    resolutionsquare->SetMarkerStyle(20);
    //    resolutionsquare->SetLineColor(kRed+1);

    
    resolutionsquare->SetLineWidth(2);
    resolutionsquare->SetMarkerColor(kRed+1);
    resolutionsquare->SetMinimum(0.);  
    resolutionsquare->Draw();
    c2->Write();


    //trying adding errors from data
    TCanvas *c2e = new TCanvas("c2e","resolutionsquareWE",700,700);
    TGraphErrors *resolutionsquareWE = new TGraphErrors(energies,INVBEAM2,RESOLUTION2,INVBEAM2error,RESOLUTION2error);//"resolutionsquareWE"," ",100,0.,1.);

    c2e->cd();
    resolutionsquareWE->SetTitle("3M Track Resolution^{2}  at DUT;beam momentum scan ;resolution at DUT #left[#mum^{2}#right]");
    resolutionsquareWE->GetYaxis()->SetRangeUser(0.,300.);
    resolutionsquareWE->GetXaxis()->SetRangeUser(0.,1.);
    resolutionsquareWE->SetMarkerStyle(20);
    //    resolutionsquareWE->SetLineColor(kRed+1);


    TF1 *fit32 = new TF1("fit32","pol1", 0., 0.8);
    fit32->SetLineColor(kBlue);
    //  fit1->SetParameter(0,3);
    //fit1->SetParameter(1,15);
    fit32->SetParName(0,"#sigma_{hit}^{2}");
    fit32->SetParName(1,"#sigma_{MS}^{2}");

    resolutionsquareWE->Fit("fit32","R");


    ostringstream strfit2[2];
    TString ss_fit2[2];
    ostringstream strfit_err2[2];
    TString ss_fit_err2[2];

    for(int i=0; i<2;i++)
      {
	strfit2[i] << setprecision(3) << sqrt(fit32->GetParameter(i));
	ss_fit2[i]=strfit2[i].str();
	strfit_err2[i] << setprecision(1) << 0.5*fit32->GetParError(i)/sqrt(fit32->GetParameter(i));
	ss_fit_err2[i]=strfit_err2[i].str();
	

      }

    cout << "#sigma_{hit} = (" << ss_fit2[0] << " #pm " << ss_fit_err2[0] << ") #mum" << endl;
    cout << "#sigma_{MS} = (" << ss_fit2[1] << " #pm " << ss_fit_err2[1] << ") #mum*GeV" << endl;
    
    TLatex Tl_22;
    Tl_22.SetTextAlign(12);
    Tl_22.SetTextSize(0.05);
    Tl_22.DrawLatexNDC(0.2,0.8,"#sigma_{hit} = ("+ss_fit2[0]+" #pm "+ss_fit_err2[0]+") #mum");

    Tl_22.SetTextSize(0.04);
    Tl_22.DrawLatexNDC(0.2,0.72,"#sigma_{MS} = ("+ss_fit2[1]+" #pm "+ss_fit_err2[1]+") #mum*GeV");

    TPaveText pt(0.3,0.82,0.5,0.89,"NDC");
    pt.SetTextSize(0.05);
    pt.AddText("#sigma_{hit} = ("+ss_fit2[0]+" #pm "+ss_fit_err2[0]+") #mum");
    pt.SetBorderSize(0);
    pt.SetFillColor(kWhite);
    pt.Draw("same");

    

    
    resolutionsquareWE->SetLineWidth(2);
    resolutionsquareWE->SetMarkerColor(kRed+1);
    resolutionsquareWE->SetMinimum(0.);
    resolutionsquareWE->Draw();
    c2e->Write();
    


  
  // Write result to file
  out->Write();
  return 0;
}
