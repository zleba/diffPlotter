//********************************************************************
//     
//     fnlo-tk-h1diffpdf.cc
//     Program to read fastNLO v2 tables and derive
//     QCD cross sections using PDFs e.g. from LHAPDF
//     
//********************************************************************
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include "fastnlotk/fastNLODiffReader.h"
#include "fastNLODiffAlphas.h"

#include "TH1D.h"


// Function prototype for flexible-scale function 
double Function_Mu(double s1, double s2 );

double Function_q(double q, double pt );
double Function_pt(double q, double pt );


vector<double> getXsec(string str, DPDFset &dpdf, int imem)
{
  fastNLODiffAlphas fnlodiff(str);
  fnlodiff.setDPDF(dpdf);
  fnlodiff.SetLHAPDFMember(imem);

  //  If you want to receive your cross section in
  //   pb/GeV or in pb. Here we choose pb/GeV
  fnlodiff.SetUnits(fastNLO::kPublicationUnits);
   fnlodiff.SetFit(fastNLODiffAlphas::Fit19);
    fnlodiff.SettIntegratedRange(-1.);
    fnlodiff.SetProtonE(920.);
  


	fnlodiff.SetContributionON(fastNLO::kFixedOrder,0,true);
	fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,true);
	//fnlodiff.SetContributionON(fastNLO::kFixedOrder,2,false);
	fnlodiff.SetContributionON(fastNLO::kFixedOrder,2,true);


	//cout << "The file name is " << argv[1] << endl;


  // Set the xpom integration interval and method
  // -------- Boris LRG dijets
  //fnlodiff.SetXPomLinSlicing( 5, xpomBins[i],  xpomBins[i+1]); // e.g. Boris LRG dijets
  fnlodiff.SetXPomLinSlicing( 30, 0.000 ,  0.030 ); // e.g. Radek VFPS dijets

  fnlodiff.SetExternalFuncForMuR (&Function_Mu);
  fnlodiff.SetExternalFuncForMuF (&Function_Mu);
  fnlodiff.SettIntegratedRange(-1.);

  vector<double>  xs = fnlodiff.GetDiffCrossSection();
  fnlodiff.PrintCrossSections();
  return xs;

}










//__________________________________________________________________________________________________________________________________

int main(int argc, char** argv){
	

  // namespaces
  using namespace std;
  using namespace say;		// namespace for 'speaker.h'-verbosity levels
  using namespace fastNLO;	// namespace for fastNLO constants

	SetGlobalVerbosity(ERROR);


  //Overall had-corr 4%
  vector<double> xpomBins = {0.0025, 0.0050, 0.0079, 0.0126, 0.0199, 0.0300};
  //vector<double> hadCorr = {0.0025, 0.0050, 0.0079, 0.0126, 0.0199, 0.0300};

vector<double> hadCorr = {
1.33391458703146,
1.18524590163938,
1.0896731753159,
1.02595802443358,
0.9834812571787 };




  //for(int i = 0; i < xpomBins.size() - 1; ++i) {
  
    vector<string> files = {
"../tables/nnlojet/FPS/nnlo/H1-LQall-8.diff_FPS_ptj1.newver.tab",
"../tables/nnlojet/FPS/nnlo/H1-LQall-8.diff_FPS_deltaeta.newver.tab",
"../tables/nnlojet/FPS/nnlo/H1-LQall-8.diff_FPS_q2.newver.tab",
"../tables/nnlojet/FPS/nnlo/H1-LQall-8.diff_FPS_y.newver.tab",
"../tables/nnlojet/FPS/nnlo/H1-LQall-8.diff_FPS_xi2zIP.newver.tab",
"../tables/nnlojet/VFPS/nnlo/H1-LQall-8.diff_VFPS_etaavg.newver.tab",
"../tables/nnlojet/VFPS/nnlo/H1-LQall-8.diff_VFPS_etadel.newver.tab",
"../tables/nnlojet/VFPS/nnlo/H1-LQall-8.diff_VFPS_ptavg_12.newver.tab",
"../tables/nnlojet/VFPS/nnlo/H1-LQall-8.diff_VFPS_ptj1.newver.tab",
"../tables/nnlojet/VFPS/nnlo/H1-LQall-8.diff_VFPS_q2.newver.tab",
"../tables/nnlojet/VFPS/nnlo/H1-LQall-8.diff_VFPS_xi2zIP.newver.tab",
"../tables/nnlojet/VFPS/nnlo/H1-LQall-8.diff_VFPS_y.newver.tab",
"../tables/nnlojet/VFPS/nnlo/H1-LQall-8.diff_VFPS_yMx.newver.tab",
    };


    DPDFset dpdf   = DPDFset("H1_DPDF_2006B_NLO_pom", "H1_DPDF_2006B_NLO_reg");

    vector<vector<vector<double>>> res(files.size());
    for(int k = 0; k < res.size(); ++k)
        res[k].resize(31);

#pragma omp parallel for
    for(int k = 0; k < res[0].size(); ++k) {
        for(int i = 0; i < files.size(); ++i) {
            res[i][k] = getXsec(files[i], dpdf, k);
        }
    }


    /*
    for(auto r : res) {
        for(auto k : r)
            cout << k <<" ";
        cout << endl;
    }
    */
    return 0;
  //say::SetGlobalVerbosity(say::DEBUG);
  fastNLODiffAlphas fnlodiff("../tables/nnlojet/FPS/nnlo/H1-LQall-8.diff_FPS_ptj1.newver.tab");
  
  //  If you want to receive your cross section in
  //   pb/GeV or in pb. Here we choose pb/GeV
  fnlodiff.SetUnits(fastNLO::kPublicationUnits);
   fnlodiff.SetFit(fastNLODiffAlphas::ABCDE);
    fnlodiff.SettIntegratedRange(-1.);
    fnlodiff.SetProtonE(920.);
  
	fnlodiff.SetContributionON(fastNLO::kFixedOrder,0,true);
	fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,true);
	//fnlodiff.SetContributionON(fastNLO::kFixedOrder,2,false);
	fnlodiff.SetContributionON(fastNLO::kFixedOrder,2,true);


	//cout << "The file name is " << argv[1] << endl;


  // Set the xpom integration interval and method
  // -------- Boris LRG dijets
  //fnlodiff.SetXPomLinSlicing( 5, xpomBins[i],  xpomBins[i+1]); // e.g. Boris LRG dijets
  fnlodiff.SetXPomLinSlicing( 30, 0.000 ,  0.030 ); // e.g. Radek VFPS dijets

  fnlodiff.SetExternalFuncForMuR (&Function_Mu);
  fnlodiff.SetExternalFuncForMuF (&Function_Mu);
  fnlodiff.SettIntegratedRange(-1.);
  fnlodiff.PrintCrossSections();

  /*
  vector<double> bins = { 
      0.0032,
      0.0063,
      0.0126,
      0.0251,
      0.0501,
      0.1000,
      0.1995,
      0.3981};

  TH1D *hBeta = new TH1D("hBeta", "Beta", bins.size()-1, bins.data());

  int nBins = 120;
  double step = 0.03 / nBins;


  for(int i = 0; i < nBins; ++i) {
    double xpLow = i*step;
    double xpHi  = (i+1)*step;
    
    double xpavg = (xpLow + xpHi) / 2;
   
    fnlodiff.SetXPomLinSlicing( 1, xpLow,  xpHi); // e.g. Boris LRG dijets

    vector<double>  xs = fnlodiff.GetDiffCrossSection();
    for(int k = 0; k < xs.size(); ++k) {
        cout << "Helenka " << k <<" "<< fnlodiff.GetObsBinLoBound(k,0) << endl;
        double xBj = (fnlodiff.GetObsBinLoBound(k,0) + fnlodiff.GetObsBinUpBound(k,0))/2;

        double w = fnlodiff.GetObsBinUpBound(k,0) - fnlodiff.GetObsBinLoBound(k,0);

        double beta = xBj / xpavg;

        hBeta->Fill(beta, w*xs[k]);
    }
  }

  hBeta->Scale(1.0, "width");

  for(int i = 1; i <= hBeta->GetNbinsX(); ++i)
     cout << i <<" "<< hBeta->GetBinContent(i) << endl;
  
  return 0;
*/


  //fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,false);
  
    //fnlodiff.SetExternalFuncForMuR (&Function_Mu);
    //fnlodiff.SetExternalFuncForMuF (&Function_Mu);


  //   fnlodiff.SetAlphasMz(0.1180);

  //  const int nStep = 1;
  

//    double xpom[nStep] = {0.0299028233407105 };
//    double dxpom[nStep] = {0.000315879371143547};

 
//  fnlodiff.SetXPomSlicing( nStep, xpom, dxpom );

  /*
  // -------- Radek VFPS dijets
  fnlodiff.SetExternalFuncForMuR (&Function_Mu);
  fnlodiff.SetExternalFuncForMuF (&Function_Mu);
  fnlodiff.SettIntegratedRange(-0.6);
  fnlodiff.SetXPomLinSlicing( 30, 0.010 ,  0.024 ); // e.g. Radek VFPS dijets
  
  // Radek's slicing
//   int nStep = 56;
//   double xpom[56] = { 0.01025, 0.01075, 0.01125, 0.01175, 0.01225, 0.01275, 0.01325, 0.01375, 0.01425, 0.01475, 0.01525, 0.01575, 0.01625, 0.01675, 0.01725, 0.01775, 0.01825, 0.01875, 0.01925, 0.01975, 0.02025, 0.02075, 0.02125, 0.02175, 0.02225, 0.02275,  0.02325, 0.02375,  };
//   double dxpom[] = {0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005 , 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005};
//   fnlodiff.SetXPomSlicing( nStep, xpom, dxpom );
  
  // switch NLO contribution off
  //fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,false);
     
  // make your scale definition (see above)
  //fnlodiff.SetScaleFactorsMuRMuF(1.0,1.0);
  */
  
  // calculate and access the cross section
  vector<double>  xs = fnlodiff.GetDiffCrossSection();
  fnlodiff.PrintCrossSections();
  cout << "Print out finished" << endl;

  return 0;
  //fnlodiff.PrintCrossSectionsWithReference();

  cout<<"double[5] = {"
      <<xs[0]<<", "
      <<xs[1]<<", "
      <<xs[2]<<", "
      <<xs[3]<<", "
      <<xs[4]<<"};"<<endl;
     

  double xSec = 0;
  for(unsigned i = 0; i < xs.size(); ++i) {
      double s = abs( fnlodiff.GetObsBinLoBound(i,0) - fnlodiff.GetObsBinUpBound(i,0) );
      xSec += s * xs[i];
  }
  //cout << "Total is " << xpomBins[i]<<" "<<  xSec / (xpomBins[i+1] - xpomBins[i]) * hadCorr[i] / 1.2  << endl;

  //}


  return 0;

  /*
  
  const int nbins = fnlodiff.GetNObsBin();
  const int npdfall = 30;
  // const int npdf = 30;
  vector<vector<double> > xsPDFerr;
  xsPDFerr.push_back(xs);
  for ( int i = 1; i<=npdfall ; i++ ) {
     fnlodiff.SetLHAPDFMember(i);
     //fnlodiff.CalcCrossSection();
     vector<double>  xsErr = fnlodiff.GetDiffCrossSection();
     xsPDFerr.push_back(xsErr);
     for ( int b = 0 ; b<nbins ; b++ ) {
	cout<<"  "<<xsErr[b]/xs[b]*100.<<endl;
     }
     cout<<endl;
  }

  
  vector<double> PDFerrUp, PDFerrDn, PDFerrMasterUp, PDFerrMasterDn;
//   cout<<" bin     \t\tDPDF error"<<endl;
//   for ( int i = 1; i<=npdf ; i++ ) {
//      cout<<"\n-------- PDFset "<<i<<" --------"<<endl;
//       for ( int b = 0 ; b<nbins ; b++ ) {
// 	 //cout<<"  "<<b<<"\t"<<xsPDFerr[i][b]<<"\t\t"<<xsPDFerr[i][b]/xsPDFerr[0][b]<<endl;
// 	 cout<<"  "<<b<<"\t\t"<<xsPDFerr[i][b]/xsPDFerr[0][b]<<endl;
//       }    
//   }  
  //   cout<<"\n===================================="<<endl<<endl;
  //   for ( int i = npdf+1; i<=npdfall ; i++ ) {
  //      cout<<"\n-------- PDFset "<<i<<" --------"<<endl;
  //       for ( int b = 0 ; b<nbins ; b++ ) {
  // 	 cout<<"  "<<b<<"\t"<<xsPDFerr[i][b]<<"\t\t"<<xsPDFerr[i][b]/xsPDFerr[0][b]<<endl;
  //       }    
  //   }  

  //zpom 0.2 0.4 646.5 646.1 646.3 643.6 648.4 643.3 648.2 663.8 628.1 641.4 622 603 659 644.5 647.5 682 612.2 613.6 676 643.9 647.4 643.6 648 583.8 667.8 653.7 637.1 650.6 640.1 610.7 681 
  cout<<endl;
  for ( unsigned int b = 0 ; b<xs.size() ; b++ ) {
      cout<<argv[1];
      //     cout<<" "<<fnlodiff.GetLoBin(b,0)<<" "<<fnlodiff.GetUpBin(b,0);
      cout<<" "<<fnlodiff.GetObsBinLoBound(b,0)<<" "<<fnlodiff.GetObsBinUpBound(b,0);
      for ( int i = 0; i<=npdfall ; i++ ) {
          cout<<" "<<xsPDFerr[i][b];
      }
      cout<<endl;
  }
  cout<<endl;
  */
  
  return 0;
}



//__________________________________________________________________________________________________________________________________

double Function_q(double q, double pt ){
  // --- fastNLO user: This is an example function
  //     to demonstrate how you might perform the
  //     definition of the scales using a 
    return sqrt(q*q);
}
double Function_pt(double q, double pt ){
  // --- fastNLO user: This is an example function
  //     to demonstrate how you might perform the
  //     definition of the scales using a 
    return 1.1*sqrt(pt*pt);
}



double Function_Mu(double s1, double s2 ){
  // --- fastNLO user: This is an example function
  //     to demonstrate how you might perform the
  //     definition of the scales using a 
  //     'flexible-scale'-table
   //double mu = s1*exp(0.3*s2);

   //double mu = sqrt(s1*s1/4. + s2*s2);
   double mu = sqrt(s1*s1 + s2*s2);
   return mu;
}

//__________________________________________________________________________________________________________________________________
