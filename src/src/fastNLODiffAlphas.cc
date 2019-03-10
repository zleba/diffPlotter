#include "fastNLODiffAlphas.h"

fastNLODiffAlphas::fastNLODiffAlphas(string filename) : fastNLODiffReader(filename), fAlphasMz(0.1184) ,fierr(0),  ftint(-1.) {

}


//______________________________________________________________________________


// double fastNLODiffAlphas::EvolveAlphas(double Q) const {
//    // --- fastNLO user:
//    // Implementation of Alpha_s evolution as function of the
//    // factorization scale [and alphas(Mz)].
//    //

//    static const int NF=4; // from h12006B_wrapper.h
//    static const double b0 = (11. - 2./3.*NF);  // The beta coefficients of the QCD beta function
//    static const  double b1 = (51. - 19./3.*NF);

//    //      double lmd = 0.399; // according to matthias
//    static const double lmd = 0.3395; // according to matthias
//    double t = log(Q/lmd);
//    double asMz = 1.0/(b0*t);

//    return asMz*(1.0-b1/b0*asMz*log(2.0*t)) *TWOPI ;
// }


//______________________________________________________________________________


bool fastNLODiffAlphas::InitPDF() {
   // --- fastNLO user:
   //  Initalize PDF parameters if necessary
   //
   // nothing todo!
   return true;
}


//______________________________________________________________________________



vector<double> fastNLODiffAlphas::GetDiffXFX(double xpom, double zpom, double muf) const {
   //
   //  GetDiffXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  xpom, zpom and factorisation scale.
   //

   //vector < double > xfxFitB(13), xfxFitJets(13), xfxZEUS(13);
	vector <double> xfx(13, 0.0);
   //diffpdf_(&xpom,&zpom,&muf,&xfx[0]);
   int ierr = fierr;
   double tmax = ftint;
   //cout <<"RADEK "<< xpom <<" "<< zpom <<" "<< muf <<" "<< tmax << endl;
   //cout << "Tmax is " << tmax << endl;

   dpdf.zfzQ2xp (ierr, zpom, muf*muf, xpom, 0, abs(tmax), xfx);

   //cout << "Is only quark " << onlyQuarks << endl;
   if(onlyQuarks) xfx[6] = 0;
	//double sumB = accumulate(xfxFitB.begin(),xfxFitB.end(), 0.0);
	//double sumJ = accumulate(xfxFitJets.begin(),xfxFitJets.end(), 0.0);
   //cout << "Fit output " << xfxFitB[6] <<" "<< xfxFitJets[6] <<" "<< xfxZEUS[6]<< endl;
   //cout << "Fit output " << sumB <<" "<< sumJ << endl;
   //logger.info<<"xpom="<<xpom<<"\tzpom="<<zpom<<"\tmuf="<<muf<<"\tgluon = "<<xfx[6]<<endl;
   return xfx;
}


//______________________________________________________________________________



//______________________________________________________________________________
// Alpha_s krams
//______________________________________________________________________________
double fastNLODiffAlphas::GetAlphasMz() const {
   return fAlphasMz;
};

void fastNLODiffAlphas::SetMz(double Mz) {
   Alphas::SetMz(Mz);
}

void fastNLODiffAlphas::SetNFlavor(int nflavor) {
   if (nflavor == 0) {
      Alphas::SetFlavorMatchingOn(true);
      Alphas::SetNf(6);
      logger.warn["SetNFlavor"]<<"GRV evolution of alpha_s is implemented for Nf=5 only.\n";
      logger.warn["SetNFlavor"]<<"You chose a variable Nf with Nfmax=6, i.e. results for Nf other than 5 presumably are wrong!\n";
   } else if (nflavor == 5) {
      Alphas::SetNf(nflavor);
   } else {
      logger.error["SetNFlavor"]<<"GRV evolution of alpha_s is implemented for Nf=5 only.\n";
      exit(1);
   }
}

void fastNLODiffAlphas::SetNLoop(int nloop) {
   Alphas::SetNLoop(nloop);
}

void fastNLODiffAlphas::SetAlphasMz(double AlphasMz , bool ReCalcCrossSection) {
   logger.debug["SetAlphasMz"]<<"Setting alpha_s(Mz)="<<AlphasMz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the alpha_s value at M_Z
   //
   fAlphasMz    = AlphasMz;             // new alpha_s value
   if (ReCalcCrossSection) CalcCrossSection();
}


double fastNLODiffAlphas::EvolveAlphas(double Q) const {
   // --------------------------------
   // code like nlojet++ ctep-pdf.h
   // --------------------------------
//     static const double alphasMZ = 0.1179;
//     static const double b0  = 1.2202;
//     static const double b1  = 0.4897;
//     static const double Mz     = 91.70;
//     double L = log(Q/Mz);
//     L = (b0 + alphasMZ*b1)*L;
//     return alphasMZ/(1.0 + alphasMZ*L);
   // --------------------------------
   // code like h12006B_wrapper.h
   // --------------------------------
//     static const int NF=5; // from h12006B_wrapper.h
//     static const double b0 = (11. - 2./3.*NF);  // The beta coefficients of the QCD beta function
//     static const  double b1 = (51. - 19./3.*NF);
//     const double pi = 3.14159265; 
//     //      double lmd = 0.399; // according to matthias
//     static const double lmd = 0.3395; // according to matthias
//     double t = log(Q/lmd);
//     double asMz = 1.0/(b0*t);
//     return asMz*(1.0-b1/b0*asMz*log(2.0*t)) *2.*pi ;
   // --------------------------------
   // Code like Boris has implemented in pdf-cteq_err.h
   // --------------------------------
//      double mr2=Q*Q;
//      double lm = 0.3395; // according to matthias
//      double lmd = lm*lm;
//      unsigned int nf = 5;
//      double pi = 3.14159265;
//      double b0 = 11.0 - 2.0/3.0*nf;
//      double t = log(mr2/lmd);
//      double b1 = (51.0 - 19.0/3.0*nf);
//      double as_nlo = 2./(b0*t)*(1. - 2*b1/(b0*b0)*log(t)/t)*2.*pi;
//      return as_nlo; 
   // --------------------------------
   // Code like Boris IS using in nlo-process_i1f0.h 
   // Now he is using lm = 0.226 to 0.228 !!!
   // --------------------------------
//     double mr2=Q*Q;
//     double lm = 0.228; // 
//     double lmd = lm*lm;
//     int nf = 5;
//     double b0 = (33.0 - 2.0*nf)/6.;
//     double t = std::log(mr2/lmd);
//     double b1 = (153.-19.*nf)/6.;
//     double as = 1./(b0*t)*(1. - b1*std::log(t)/(t*std::pow(b0, 2)));
//     as*=2*3.14159265; // two pi for fastNLO
//     return as;
   // --------------------------------
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   return Alphas::CalcAlphasMu(Q , fAlphasMz);

}





void fastNLODiffAlphas::SetGRVtoPDG2012_2loop() {
   logger.info["SetGrVtoPDF2012"]<<"Resetting to GRV Alphas::Alphas evolution."<<endl;
   Alphas::SetMz(91.1876); // PDG 2012
   Alphas::SetNf(5);
   Alphas::SetNLoop(2);
   Alphas::SetFlavorMatchingOn(false);
   if (logger.info.GetSpeak()) {
      logger.info<<"Calling Alphas::PrintInfo()."<<endl;
      logger.info<<"Alpha_s(Mz) value is taken from fastNLODiffAlphas, instead of Alphas::Alphas."<<endl;
      Alphas::PrintInfo();
   }
}
