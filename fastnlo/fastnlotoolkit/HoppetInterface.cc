#include "fastnlotk/HoppetInterface.h"
#include "fastnlotk/fastNLOReader.h"
#include "hoppet_v1.h"
#include <LHAPDF/LHAPDF.h>
#include "fastnlotk/speaker.h"
#include <cstdlib>


// Default values (PDG)
bool HoppetInterface::IsInitialized = false;
fastNLOReader *HoppetInterface::fnlo = NULL;
double HoppetInterface::fAlphasMz = PDG_ASMZ;
double HoppetInterface::fMz = PDG_MZ;
int HoppetInterface::fnFlavor = -1;
int HoppetInterface::fnLoop = 1;
double HoppetInterface::QMass[6] = {PDG_MD, PDG_MU, PDG_MS, PDG_MC, PDG_MB, PDG_MT};

void HoppetInterface::InitHoppet(fastNLOReader& lfnlo) {

   if (!IsInitialized) {
      StartHoppet();
      fnlo = &lfnlo;
   }

   //If fnFlavor smaller than 1 use VFNS (NNPDF reports nf=-1)
   hoppetsetpolemassvfn_(QMass[3], QMass[4], 10000000000.0);
   say::info["InitHoppet"] << "Using variable-flavour number scheme with the the given masses."
                          << " M_Top is set to 10000000000.0 GeV, so effectively nf_max = 5." <<  std::endl;
   // Carry out evolution
   hoppetEvolve(fAlphasMz, fMz, 2, 1.0, &HoppetInterface::LHAsub, 2.00001);
   // Fills the HOPPET PDF represenation using PDF provided by LHAPDF.
   //hoppetAssign(&evolvepdf_);
}

void HoppetInterface::StartHoppet(){
   double ymax = 12.0;
   double dy = 0.1;
   int order = -6;
   double dlnlnQ = dy/4.0;
   double Qmin = 1.0;
   double Qmax = 28000;
   int fnLoop = 2.;
   hoppetStartExtended( ymax, dy, Qmin, Qmax, dlnlnQ, fnLoop, order, factscheme_MSbar);
   IsInitialized = true;
}

void HoppetInterface::LHAsub(const double & x, const double & Q, double * pdf) {
   //
   //Provides PDF for Hoppet
   for (int i=0; i<13; i++)
   {
   pdf[i] = HoppetInterface::fnlo->GetXFX(x,Q)[i];
   }
}


std::vector<double> HoppetInterface::GetSpl(double x, double Q){
   if (!IsInitialized) {
      say::error["GetSpl"] << "Hoppet not correctly initialized!" <<  std::endl;
      std::exit(1);
   }

   //! Returns the splitting functions
   static std::vector<double> xfx(13);
   hoppetEvalSplit(x, Q, 1, 5, &xfx[0]);
   return xfx;
}

std::vector<double> HoppetInterface::GetXFX(double x, double Q){
   //! Returns PDF
   if (!IsInitialized) {
      say::error["GetSpl"] << "Hoppet not correctly initialized!" <<  std::endl;
      std::exit(1);
   }
   static std::vector<double> xfx(13);
   hoppetEval(x, Q, &xfx[0]);
   return xfx;
}

double HoppetInterface::EvolveAlphas(double Q) {
   if (!IsInitialized) {
      say::error["EvolveAlphas"] << "Hoppet not correctly initialized!" <<  std::endl;
      std::exit(1);
   }
   return hoppetAlphaS(Q);
}
