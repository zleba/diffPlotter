// Author: Daniel Britzger
// DESY, 02/04/2012

#ifndef fASTNLODIFFALPHAS
#define fASTNLODIFFALPHAS


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLODiffUSER                                                     //
//                                                                      //
//  FastNLODiffReader is a standalone code for reading                  //
//  diffractive FastNLO tables of version 2.0 for DIS processes         //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include "fastnlotk/fastNLODiffReader.h"
#include "fastnlotk/Alphas.h"
#include <numeric>
#include <cassert>
#include "DPDFset.h"

using namespace std;

//void LoadGrid(const string& Label, int n=0);
//void Set_tmin(double t);

/*
extern "C" {
   //   void diffpdf_(double* xpom, double*  zpom, double*  Q2, double *pdfs);
   void diffpdferr_(double* xpom, double*  zpom, double*  Q2, int* ipdf, int* ierr, double *pdfs, double *tmax);
   void strpriv_(double* X, double *MUF, double *xpom, double *tcut, double *XPQ);
   void initialisedpdf_(int *iSet);
   void getdpdf_(double *xPom, double *tcut, double *z, double *Qsq, const char *PomOrReg, double *DPDF);
}
*/


class fastNLODiffAlphas : public fastNLODiffReader {

public:
	enum Fits{ FitB, FitA, FitJets, zeusSJ, Fit19, ABCDE, MRW};

    bool onlyQuarks = false;
	Fits fit;

   fastNLODiffAlphas(string filename);
   ~fastNLODiffAlphas(void) {
      ;
   };

   // ---- Alphas vars ---- //
   // Setters
   void SetMz(double Mz);
   void SetNFlavor(int nflavor);
   void SetNLoop(int nloop);
   void SetAlphasMz(double AlphasMz , bool ReCalcCrossSection = false);
   // Getters
   double GetAlphasMz() const;
   void SetGRVtoPDG2012_2loop();

   void SetLHAPDFFilename(const char* sth) const {;}
   void SetLHAPDFMember(const int a) {fierr=a;}
   void SetFit(Fits fit_) {
        fit = fit_;
        //if(fit == zeusSJ) {LoadGrid("zeuspdf/grids/zeusD_SJ"); Set_tmin(ftint);} 
        //if(fit == MRW) {int iset = 1; initialisedpdf_(&iset); }
   }
   void IncludeOnlyQuarks(bool st =true) {onlyQuarks = st;}

   int GetFitID() const { return fifit;}
   void SetFitID( const int ifit) { fifit = ifit; } // set fit id: fitA=1, fitB=2

   void SettIntegratedRange(double tmax) { ftint = tmax;   }
   double GettIntegratedRange() const {return  ftint;}

   void setDPDF(DPDFset &dpdf_) { dpdf = dpdf_; }

protected:

   // inherited functions
   double EvolveAlphas(double Q) const ;
   bool InitPDF();
   vector<double> GetDiffXFX(double xpom, double zpom, double muf) const ;

   // ---- Alphas vars ---- //
   double fAlphasMz;

   // ----- diff pdf ----- /
   int fierr;
   int fifit;
   double ftint;
   vector<DPDFset> dpdfsLha;
   DPDFset dpdf;

public:
   // Dummy functions for fitting code
   int GetNPDFMembers() const { 
      if ( fifit == 2 ) return 30;
      else if ( fifit == 1 ) return 32;}


};

//______________________________________________________________________________


	void ZeusDpdf3(double xP, double zP, double QQ, double f[7], int xpow=1);



//______________________________________________________________________________


//vector<DPDFset> fastNLODiffAlphas::dpdfsLha = {};


fastNLODiffAlphas::fastNLODiffAlphas(string filename) : fastNLODiffReader(filename), fAlphasMz(0.1184) ,fierr(0), fifit(2), ftint(-1.) {
	fit = FitB;	

    if(dpdfsLha.size() == 0) {
        dpdfsLha.resize(6);
        //auto pdfs  = mkPDFs("H1_DPDF_2006B_NLO_pom");

        /*
        dpdfsLha[FitB]    = DPDFset("H1_DPDF_2006B_NLO_pom", "H1_DPDF_2006B_NLO_reg");
        cout << "Bint OK" << endl;
        dpdfsLha[FitA]    = DPDFset("H1_DPDF_2006A_NLO_pom", "H1_DPDF_2006A_NLO_reg");
        cout << "Aint OK" << endl;
        dpdfsLha[FitJets] = DPDFset("H1_DPDF_2007Jets_NLO_pom", "H1_DPDF_2007Jets_NLO_reg");
        cout << "Jetsint OK" << endl;
        dpdfsLha[zeusSJ]  = DPDFset("ZEUS_DPDF_2009SJ_NLO_pom");
        cout << "ZEUSint OK" << endl;
        dpdfsLha[Fit19]  = DPDFset("lhaTest_pom", "lhaTest_reg");
        cout << "Fit19 OK" << endl;
        dpdfsLha[ABCDE]  = DPDFset("ABCDE_pom", "ABCDE_reg");
        cout << "ABCDE OK" << endl;
        */
    }
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
	vector <double> xfx(13, 1.0);
   //diffpdf_(&xpom,&zpom,&muf,&xfx[0]);
   int ifit = fifit;
   int ierr = fierr;
   double tmax = ftint;
   //cout <<"RADEK "<< xpom <<" "<< zpom <<" "<< muf <<" "<< tmax << endl;
   //cout << "Tmax is " << tmax << endl;
	if(fit == FitB) {
		ifit = 2;
		//diffpdferr_(&xpom,&zpom,&muf,&ifit,&ierr,&xfx[0],&tmax);
        dpdfsLha[FitB].zfzQ2xp (ierr, zpom, muf*muf, xpom, 0, abs(tmax), xfx);
	}
	else if(fit == FitA) {
		ifit = 1;
		//diffpdferr_(&xpom,&zpom,&muf,&ifit,&ierr,&xfxOld[0],&tmax);
        dpdfsLha[FitA].zfzQ2xp (ierr, zpom, muf*muf, xpom, 0, abs(tmax), xfx);
	}
	else if(fit == FitJets) {
		//strpriv_(&zpom, &muf, &xpom, &tmax, &xfxOld[0]);
        dpdfsLha[FitJets].zfzQ2xp (ierr, zpom, muf*muf, xpom, 0, abs(tmax), xfx);
    }
	else if(fit == Fit19) {
		//strpriv_(&zpom, &muf, &xpom, &tmax, &xfxOld[0]);
        //dpdfsLha[Fit19].zfzQ2xp (ierr, zpom, muf*muf, xpom, 0, abs(tmax), xfx);

        dpdf.zfzQ2xp (ierr, zpom, muf*muf, xpom, 0, abs(tmax), xfx);
        ierr = ierr;
    }
	else if(fit == ABCDE) {
		//strpriv_(&zpom, &muf, &xpom, &tmax, &xfxOld[0]);
        dpdfsLha[ABCDE].zfzQ2xp (ierr, zpom, muf*muf, xpom, 0, abs(tmax), xfx);
    }

    /*
	else if(fit == zeusSJ) {
		double q2 = min(10000.0, max(muf*muf, 1.8));
		//cout << "ZEUS q2 " << q2 << endl;
		ZeusDpdf3(xpom, zpom, q2, &xfx[6], 1);
		for(int i = 1; i <=6; ++i)
			xfx[6-i] = xfx[6+i];
		for(int i = 0; i < 13; ++i)
			xfx[i] *= 1.2;
	}


    else if(fit == MRW) {
        double q2 = muf*muf;
	    vector <double> DPDFPom(13), DPDFReg(13);
        getdpdf_(&xpom,&tmax,&zpom,&q2,"Pom",&DPDFPom[0]);
        getdpdf_(&xpom,&tmax,&zpom,&q2,"Reg",&DPDFReg[0]);

		for(int i = 0; i < 13; ++i) {
            xfx[i] = 1.2*(DPDFPom[i] + DPDFReg[i]);
        }
    }
    */

	else 
		assert(0 && "Unknown Fit");

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



#endif
