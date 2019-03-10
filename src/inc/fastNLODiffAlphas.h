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


class fastNLODiffAlphas : public fastNLODiffReader {

public:

    bool onlyQuarks = false;

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
   void IncludeOnlyQuarks(bool st =true) {onlyQuarks = st;}


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
   double ftint;
   DPDFset dpdf;

public:
   // Dummy functions for fitting code
   int GetNPDFMembers() const { 
       return dpdf.size();
   }

};




#endif
