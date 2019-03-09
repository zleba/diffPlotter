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
#include "TApplication.h"

#include <fstream>
#include <map>
#include <set>
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TFrame.h"
#include <functional>
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TGaxis.h"
#include <cassert>
#include <algorithm>
#include <cstring>
#include <memory>
#include "RemoveOverlaps.h"
#include "plottingHelper.h"

using namespace PlottingHelper;

//TString defFile = "H1-LQall-8c.diff";
TString defFile = "newver";
//TString tablesDir = "/home/radek/moje/daniel/tables/";
TString tablesDir = "../tables/"; //  /home/radek/moje/daniel/tables/";

const map<TString, TString> analMap = {
    {"FPS",      "#splitline{ H1 FPS}{(HERA #Iota#Iota)}"},
    {"VFPS",     "#splitline{H1 VFPS}{(HERA #Iota#Iota)}"},
    {"LRG",      "#splitline{ H1 LRG}{(HERA #Iota#Iota)}"},
    {"LRG_H1",   "#splitline{ H1 LRG}{(HERA #Iota)}"},
    {"LRGH1820", "#splitline{ H1 LRG}{(300 GeV)}"},
    {"LRGZEUS",  "#splitline{ZEUS LRG}{  (HERA #Iota)}"},
};
const map<TString, TString> analMapInline = {
    {"FPS",      "H1 FPS (HERA #Iota#Iota)"},
    {"VFPS",     "H1 VFPS (HERA #Iota#Iota)"},
    {"LRG",      "H1 LRG (HERA #Iota#Iota)"},
    {"LRG_H1",   "H1 LRG (HERA #Iota)"},
    {"LRGH1820", "H1 LRG (300 GeV)"},
    {"LRGZEUS",  "ZEUS LRG (HERA #Iota)"},
};


//const vector<TString> arrQ[]={
	//"4 < Q^{2} < 6 GeV^{2}", "6 < Q^{2} < 10 GeV^{2}", "10 < Q^{2} < 18 GeV^{2}",
    //"18 < Q^{2} < 34 GeV^{2}", "34 < Q^{2} < 100 GeV^{2}",""};
const map<TString,TString> names2D {
    {"ptjet1_q2_4_6",   "4 < Q^{2} < 6 GeV^{2}"},
    {"ptjet1_q2_6_10",  "6 < Q^{2} < 10 GeV^{2}"},
    {"ptjet1_q2_10_18", "10 < Q^{2} < 18 GeV^{2}"},
    {"ptjet1_q2_18_34", "18 < Q^{2} < 34 GeV^{2}"},
    {"ptjet1_q2_34_100","34 < Q^{2} < 100 GeV^{2}"},

    {"zpom_q2_4_10",   "4 < Q^{2} < 10 GeV^{2}"},
    {"zpom_q2_10_20",  "10 < Q^{2} < 20 GeV^{2}"},
    {"zpom_q2_20_40",  "20 < Q^{2} < 40 GeV^{2}"},
    {"zpom_q2_40_100", "40 < Q^{2} < 100 GeV^{2}"},


    {"zpom_ptjet1_5_6p5",  "5 < p_{T}^{*jet1} < 6.5 GeV"},
    {"zpom_ptjet1_6p5_8",  "6.5 < p_{T}^{*jet1} < 8 GeV"},
    {"zpom_ptjet1_8_16",   "8 < p_{T}^{*jet1} < 16 GeV"},

    {"zpom_q2_5_12",   "4 < Q^{2} < 12 GeV^{2}"},
    {"zpom_q2_12_25",  "12 < Q^{2} < 25 GeV^{2}"},
    {"zpom_q2_25_50",  "25 < Q^{2} < 50 GeV^{2}"},
    {"zpom_q2_50_100", "50 < Q^{2} < 100 GeV^{2}"},
};

int GetNpads(TCanvas *can)
{
    int i;
    for(i = 1; i < 100; ++i) {
        if(!can->GetPad(i)) break;
    }
    return i - 1;
}



double getChi2(TH1D *hData, TGraphAsymmErrors *grAll, TH1D *hModel);

double funQ2pPt2(double q, double pt) {return sqrt(q*q+pt*pt);}
double funQ2(double q, double pt) {return sqrt(q*q);}
double funPt2(double q, double pt) {return sqrt(pt*pt);}
double funForthQ2pPt2(double q, double pt) {return sqrt(q*q/4+pt*pt);}
double funSqrtQ4pPt4(double q, double pt) {return sqrt(sqrt(q*q*q*q+pt*pt*pt*pt));}
double funZEUSmuR(double q, double pt) { return pt; }
double funZEUSmuF(double q, double pt) { return q; }

double ScoreFitB[2]={}, ScoreFitA[2]={}, ScoreJets[2]={}, ScoreSJ[2]={};
double ScoreFitBup[2]={}, ScoreFitBdn[2]={};
double ScoreJetsup[2]={}, ScoreJetsdn[2]={};


#define SF TString::Format 


// Function prototype for flexible-scale function 
double Function_Mu(double s1, double s2 );

class Histogram;

struct HistoErr;

void PlotTotal(TCanvas *can, Histogram *h1, Histogram *h2=0, Histogram *h3=0, Histogram *h4=0, Histogram *h5=0, Histogram *h6=0);

void PlotTotalDPDFs(TCanvas *can, Histogram *h1, Histogram *h2=0, Histogram *h3=0, Histogram *h4=0, Histogram *h5=0, Histogram *h6=0);

double GetRandom(vector<double> &fastBinsLo, vector<double> &fastBinsHi, TVectorD &vec, int i);


void plotXi12(TString fileName, const char *xpomN, const char *zpomN);


struct Setting
{
	double renFactor, factFactor;
	function<double(double,double)> renFunc, factFunc;
	TString scaleTag;
	TString pdfName;
	int pdfVecId;
    bool onlyQuarks = false;


	TString getFileName() {
		map<double,char> scaleMap;
		scaleMap[0.5]='d';
		scaleMap[1]='c';
		scaleMap[2]='u';

		return pdfName+"-"+SF("%d",pdfVecId)+"-"+scaleTag+"-"+scaleMap[renFactor]+scaleMap[factFactor];
	}
	void CreateFitB(int pdfID, TString ScaleTag, double CoefQ2, double CoefPt2, double sclFact=1) {
        /*
        if(ScaleTag != "SqrtQ4pPt4")
            renFunc=factFunc=[=](double q, double pt) {return  sqrt(CoefQ2*q*q + CoefPt2*pt*pt);};
        else
            renFunc=factFunc=[=](double q, double pt) {return  sqrt(sqrt(CoefQ2*q*q*q*q + CoefPt2*pt*pt*pt*pt));};
        */

        if(ScaleTag == "SqrtQ4pPt4")
            renFunc=factFunc=[=](double q, double pt) {return  sqrt(sqrt(CoefQ2*q*q*q*q + CoefPt2*pt*pt*pt*pt));};
        else if(ScaleTag == "ZEUSscale") {
            renFunc =[=](double q, double pt) {return  sqrt(CoefPt2*pt*pt);};
            factFunc=[=](double q, double pt) {return  sqrt(CoefQ2*q*q);};
        }
        else
            renFunc=factFunc=[=](double q, double pt) {return  sqrt(CoefQ2*q*q + CoefPt2*pt*pt);};


		renFactor = factFactor = sclFact;
		scaleTag = ScaleTag;
		pdfName = "FitB";
		pdfVecId = pdfID;
	}


	void CreateFitA(int pdfID, TString ScaleTag, double CoefQ2, double CoefPt2, double sclFact=1) {
		renFunc=factFunc=[=](double q, double pt) {return  sqrt(CoefQ2*q*q + CoefPt2*pt*pt);};
		renFactor = factFactor = sclFact;
		scaleTag = ScaleTag;
		pdfName = "FitA";
		pdfVecId = pdfID;
	}
	void CreateFitJets(int pdfID, TString ScaleTag, double CoefQ2, double CoefPt2, double sclFact=1) {
		renFunc=factFunc=[=](double q, double pt) {return  sqrt(CoefQ2*q*q + CoefPt2*pt*pt);};
		renFactor = factFactor = sclFact;
		scaleTag = ScaleTag;
		pdfName = "FitJets";
		pdfVecId = pdfID;
	}
	void CreateZeusSJ(int pdfID, TString ScaleTag, double CoefQ2, double CoefPt2, double sclFact=1) {
		renFunc=factFunc=[=](double q, double pt) {return  sqrt(CoefQ2*q*q + CoefPt2*pt*pt);};
		renFactor = factFactor = sclFact;
		scaleTag = ScaleTag;
		pdfName = "zeusSJ";
		pdfVecId = pdfID;
	}

	void CreateFit19(int pdfID, TString ScaleTag, double CoefQ2, double CoefPt2, double sclFact=1) {
		renFunc=factFunc=[=](double q, double pt) {return  sqrt(CoefQ2*q*q + CoefPt2*pt*pt);};
		renFactor = factFactor = sclFact;
		scaleTag = ScaleTag;
		pdfName = "Fit19";
		pdfVecId = pdfID;
	}


	void CreateMRW(int pdfID, TString ScaleTag, double CoefQ2, double CoefPt2, double sclFact=1) {
		renFunc=factFunc=[=](double q, double pt) {return  sqrt(CoefQ2*q*q + CoefPt2*pt*pt);};
		renFactor = factFactor = sclFact;
		scaleTag = ScaleTag;
		pdfName = "MRW";
		pdfVecId = pdfID;
	}
};




//__________________________________________________________________________________________________________________________________
pair<double,double> GetSystTot(vector<double> *sysVec, int binId) ;

struct Theory {
	int iOrder; //1--LO, 2--NLO, 3--NNLO
	TString file_name;
	vector<double> xsc;
	Setting setting;

	Theory & operator+=(const Theory &th) {
		assert(xsc.size() == th.xsc.size() );
		//cout << "Old order " << iOrder << endl;
		//cout << "New order " << th.iOrder << endl;
		//assert(iOrder == th.iOrder); TODO
		for(unsigned i = 0; i < xsc.size(); ++i) {
			xsc[i]     += th.xsc[i];
		}
		return *this;
	}

	void ReverseBinning() {
		auto reverseAll = [](vector<double> &v) { reverse(begin(v), end(v) ); };
		reverseAll(xsc);
	}

	void Print() const {
		for(double Xsc : xsc) 
			cout << Xsc << " ";
		cout << endl;
	}
};


class Histogram {

public:
	void loadData(TString File_name, TString Var_name);
	Theory LoadTheory(TString theor_file, Setting setting);

	//void LoadNLO() {   nlo = LoadTheory(nlo_file); }
	//void LoadNNLO() { nnlo = LoadTheory(nnlo_file); }

	void LoadTheories() {
		for(TString file : theor_files) 
			for(unsigned i = 0; i < settings.size(); ++i) {
				TString setTag = settings[i].getFileName();
				if(!file.Contains("H1-LQall-8c") && !file.Contains(defFile) && i>0) continue;
                //if(!file.Contains(
				theories[file+":"+setTag] = LoadTheory(file, settings[i]);
			}
	}

	vector<double> LoadHistogramsNloNnloData(TString refName, HistoErr &grData, HistoErr &grNlo, HistoErr &grNnlo, bool sysInside = false);
	void LoadHistogramsFitJetsSJ(HistoErr &grNloJets, HistoErr &grNnloJets,
											  HistoErr &grNloSJ, HistoErr &grNnloSJ);
	HistoErr LoadHistograms(TString name, TString file, TString histName, TString refName);

	void plotDPDFStudies();

	vector<double> GetBinning();
	int getNbins() {return nBins;}

	void plotNLOvsNNLO();
	void plotNLOvsNNLOratio(TString plotStyle, double minY, double maxY);
	void plotNLOvsNNLOratioAbs(TString plotStyle, double minY, double maxY, double minYabs, double maxYabs, double Factors);
	void plotScaleChoices();
	void plotScaleStudies();

    void  CalculateGluonFraction();
	map<double, vector<vector<double>>>  calculateScaleDependence(TString scaleTag, TString tag);

	TString getTag() { return theor_path(theor_path.First('/')+1, theor_path.Last('/')-theor_path.First('/')-1 ); }

	friend void PlotComparison(TCanvas *can,  TString plotStyle,  vector<Histogram*> histos);
	friend void PlotComparisonWithAbs(TCanvas *can,  TString plotStyle,  vector<Histogram*> histos);
	friend void PlotSingle(TCanvas *can,  TString plotStyle,  Histogram* hist);
	friend void plotScaleDependences(TString scaleTag, Histogram *h1, Histogram *h2, Histogram *h3, Histogram *h4, Histogram *h5);
	friend void PlotTotal(TCanvas *can, Histogram *h1, Histogram *h2, Histogram *h3, Histogram *h4, Histogram *h5, Histogram *h6);
	friend void plotXi12(TString fileName, const char *xpomN, const char *zpomN);
	void plotScaleDependence(int binId, TString scaleTag);

	void IncludeSettings(vector<Setting> settings_) {settings=settings_;}

	//Histo Style
	struct {
		bool isLogY{false}, isLogX{false};
		double legX{0.5}, legY{0.6};

	} style;

private:

	TString data_file;
	TString theor_path;
	//TString nlo_file;
	//TString nnlo_file;
	set<TString> theor_files;
	TString var_name;
	TString xTitle;
	TString yTitle;
	TString  Title;

	int nBins;
	static int npdfall;

	double dissFactor;

	vector<double> xMin;
	vector<double> xMax;

	vector<double> data;
	vector<double> dataStatErr;
	vector<double> dataSystErr;
	vector<double> radCorr;
	vector<double> hadrCorr;

	//Theory nlo, nnlo;
	map<TString, Theory> theories;
	//vector<double> nlo;
	//vector<double> nloUp;
	//vector<double> nloDown;
	//vector< vector<double> > nloSyst;


	vector<Setting> settings;

	Theory CalcNLO(TString theor_file, Setting setting);
	//double GetSystTot(int binId) const;


	void ConvertToTotal(Theory &th, vector<double> &fastBinsLo, vector<double> &fastBinsHi);//For total and xpom
	void ConvertToW(Theory &th, vector<double> &fastBinsLo, vector<double> &fastBinsHi);//For total and xpom
	Theory CalculateRawNLO( fastNLODiffAlphas  &fnlodiff,  vector<double> &fastBinsLo, vector<double> &fastBinsHi  );//without HadCorr
	Theory CalculateXpomNLO( fastNLODiffAlphas  &fnlodiff, vector<double> &fastBinsLo, vector<double> &fastBinsHi  );//Xpom Calc
	Theory CalculateVarNLO( fastNLODiffAlphas  &fnlodiff,  vector<double> &fastBinsLo, vector<double> &fastBinsHi,
	                     std::function<double(double,double,double)>  func); //xpom,q2,tableVar


	TMatrixD CreateInterpolationMatrix(vector<double> &fastBinsLo, vector<double> &fastBinsHi);

};

int Histogram::npdfall = 30;

void Histogram::loadData(TString File_name, TString Var_name)
{

	data_file = File_name;
	var_name = Var_name;

	ifstream file;
	file.open( data_file.Data() );

	string str;
	TString strR;
	do {
		getline(file, str);
		if(!file.good()) {
			cout << "Histogram data for \"" << var_name << "\" in "<< data_file<<" not exist" << endl;
			exit(1);
		}
		strR = str;
	} while( !strR.EqualTo(var_name) );

	char cStr[1000];

	//Read theor path
	getline(file, str);
	theor_path = str;
	theor_path += "/";

	//Read NLO file name
	while(1) {
		getline(file, str);

		if(str.find(':') == string::npos)
			break;

		strcpy(cStr, str.c_str());
		TString start = strtok(cStr, ":");
		start = ","+start+",";
		TString TableFile  = strtok(NULL, ":");


		if(start.Contains(",nlo,") ) {
			theor_files.insert("nlo:" +TableFile);
			//nlo_file  = "nlo:" +TableFile;
		}
		if(start.Contains(",nnlo,") ) {
			theor_files.insert("nnlo:" +TableFile);
			//nnlo_file = "nnlo:"+TableFile;
		}
	}

	//nlo_file = str;

	//Read NNLO file name
	//getline(file, str);
	//nnlo_file = str;

	//Read titles
	//getline(file, str);
	strcpy(cStr, str.c_str());
	//xTitle = strR( 0, strR.First(',') );
	//yTitle = strR( strR.First(',')+1, 1000 );
	xTitle = strtok(cStr, ",");
	yTitle = strtok(NULL, ",");
	Title  = strtok(NULL, ",");

	//Skip line
	getline(file, str);
	
	//Read table
	while(1){
		double xcntr, xmin, xmax, xsec, tot, stat, syst, hadr, hadrP, hadrM, rad;
		double systUnc, systCor;
		bool status;
		if( theor_path.Contains("/VFPS/") )
			status= static_cast<bool>( file >> xcntr>> xmin >> xmax >> xsec >> stat >> syst >> rad >> hadr );
		else if( theor_path.Contains("/LRG/") ) 	
			status= static_cast<bool>( file >> xcntr>> xmin >> xmax >> xsec >> stat >> syst >> hadr >> hadrP>> hadrM >> rad );
		else if( theor_path.Contains("/LRG_H1/") ) {	
			status= static_cast<bool>( file >> xmin >> xmax >> xsec >> tot >> stat >>
			                         systUnc >> systCor >> hadr >> hadrP  );
			double syst1 = sqrt(tot*tot - stat*stat);
			//double syst2 = hypot(systCor, systUnc);
			syst = syst1;
		}
		else if( theor_path.Contains("/LRGH1820/") ) {	
			status= static_cast<bool>( file >> xmin >> xmax >> xsec >> stat >> systCor >> tot >> hadr >> hadrP  );
			double syst1 = sqrt(tot*tot - stat*stat);
			syst = syst1;
		}
		else if( theor_path.Contains("/FPS/") ) {
			status= static_cast<bool>( file >> xmin >> xmax >> xsec >> tot >> stat >> syst >> hadr );
			stat = xsec * stat/100.;
			syst = xsec * syst/100.;
			tot  = xsec * tot/100.;
			//cout << "It is FPS " << xmin <<" "<< xmax<< endl;
		}
		else if( theor_path.Contains("/LRGZEUS/") ) {	
			double sysUp, sysDown, corrUp, corrDown, Delta;
			status= static_cast<bool>( file >> xmin >> xmax >> xsec >> stat >> sysUp >> sysDown >> corrUp>> corrDown >> Delta >> hadr  );
			syst = sqrt(sysUp*sysUp + sysDown*sysDown + corrUp*corrUp + corrDown*corrDown)/sqrt(2);
			//hadr = 1;
			if(Var_name == "xpom") { //Correction for bad units in ZEUS data
                double C = log(xmax/xmin) / (xmax-xmin);
				xsec *= C;
				stat *= C;
				syst *= C;
			}
		}
		else {
			cout << "Wrong file type " << theor_path << endl;
			exit(1);
		}
		
		if(!status) break;

		xMin.push_back(xmin);
		xMax.push_back(xmax);
		data.push_back(xsec);
		dataStatErr.push_back(stat);
		dataSystErr.push_back(syst);
		hadrCorr.push_back(hadr);

	}
	nBins = xMin.size();
	if( nBins == 0 ) {
		cout <<"filesPath "<< theor_path << endl;
		//cout <<"nloFile "<< nlo_file << endl;
		//cout <<"nnloFile "<< nnlo_file << endl;
		cout <<"xTitle "<< xTitle << endl;
		cout <<"yTitle "<< yTitle << endl;
		cout <<"Title "<< Title << endl;
		cout <<"My size " <<  xMin.size() << endl;
		exit(1);
	}
	
	file.close();
	/*
	cout << Var_name << endl;
	cout << theor_path << " " << nlo_file << endl;
	for(int i = 0; i < xMin.size(); ++i)
		cout << i <<" "<< xMin[i]<<" "<<xMax[i] << endl;
	exit(1);
	*/

}

double Function_MuDISENTren(double s1, double s2 ){
	// --- fastNLO user: This is an example function
	//     to demonstrate how you might perform the
	//     definition of the scales using a 
	//     'flexible-scale'-table
	//double mu = s1*exp(0.3*s2);
	//q2/4 + meanPt^2
	double mu = s2;
	return mu;
}

double Function_MuDISENTfac(double s1, double s2 ){
	// --- fastNLO user: This is an example function
	//     to demonstrate how you might perform the
	//     definition of the scales using a 
	//     'flexible-scale'-table
	//double mu = s1*exp(0.3*s2);
	//q2/4 + meanPt^2
	double mu = 6.2;
	return mu;
}


Theory Histogram::CalculateRawNLO( fastNLODiffAlphas  &fnlodiff,  vector<double> &fastBinsLo, vector<double> &fastBinsHi  )
{
	Theory th;

	//  If you want to receive your cross section in
	//   pb/GeV or in pb. Here we choose pb/GeV
	fnlodiff.SetUnits(fastNLO::kPublicationUnits);


	//fnlodiff.SetExternalFuncForMuR (&Function_Mu);
	//fnlodiff.SetExternalFuncForMuF (&Function_Mu);
	//fnlodiff.SetExternalFuncForMuR (&Function_MuDISENTren);
	//fnlodiff.SetExternalFuncForMuF (&Function_MuDISENTfac);



	fastBinsLo = fnlodiff.GetObsBinsLoBounds(0);
	fastBinsHi = fnlodiff.GetObsBinsUpBounds(0);

	// calculate and access the cross section
	//static int i = 0;
	//fnlodiff.SetXPomLogSlicing( 5, pow(10,xMin[i]) ,  pow(10,xMax[i]) ); 
	//++i;
	//fnlodiff.FillPDFCache();

	//fnlodiff.SetScaleFactorsMuRMuF(1.0, 1.0);
	th.xsc = fnlodiff.GetDiffCrossSection();
	
	/*
	cout << "KRISTINA start" << endl;
	for(int k = 0; k < th.xsc.size(); ++k)
		cout << k << " "<< th.xsc[k] << endl;
	cout << "KRISTINA end" << endl;

	return th;
	*/
	
	//XsUncertainty  nloUnc = fnlodiff.GetScaleUncertainty( fastNLO::EScaleUncertaintyStyle::kSymmetricTwoPoint );

	//fnlodiff.SetScaleFactorsMuRMuF(2., 2.);
	//th.xscUp=fnlodiff.GetDiffCrossSection();

	//fnlodiff.SetScaleFactorsMuRMuF(1./2., 1./2.);
	//th.xscDown=fnlodiff.GetDiffCrossSection();

	/*
		for(int i = 0; i < nBins; ++i) {
		cout << i <<" "<< nlo[i]<<" "<< nloUnc.xs[i]/nlo[i]<<" "<< fastBinsHi[i]-fastBinsLo[i] <<" "<< nloUnc.dxsl[i] << " "<<  nloUnc.dxsu[i] <<  endl;;
		}
		exit(1);
	*/

	/*
	//DPDF errors
	fnlodiff.SetScaleFactorsMuRMuF(1.0, 1.0);
	th.xscSyst.resize(npdfall+1);
	th.xscSyst[0] = th.xsc;
	for ( int i = 1; i<=npdfall ; i++ ) {
		fnlodiff.SetLHAPDFMember(i);
		th.xscSyst[i] = fnlodiff.GetDiffCrossSection();
	}
	*/


	//Convert histogram to total
	if(var_name.EqualTo("total") ||  var_name.EqualTo("xpom") ||  var_name.EqualTo("logxpom")  ) {
		//cout <<"RADEKnow "<< th.xsc[0] << " "<<th.xsc[1]<<" "<< th.xsc[2]<< endl;
		ConvertToTotal(th, fastBinsLo, fastBinsHi);
		//cout << "Cross section now RADEK " << th.xsc[0] << endl;
	}
	if(var_name == "w") {
		ConvertToW(th, fastBinsLo, fastBinsHi);
	}

	return th;

}

/// ahoj
///This method calculates NLO or NNLO differential cross section for given table file theor_file
///NLO and NNLO is specified by ":" tag, e.g. nlo:filename 
/// test

Theory Histogram::CalcNLO(TString theor_file, Setting setting)
{
    cout <<"DANIEL " <<  theor_file<< " "<< setting.getFileName() << endl;
	double normFactor = 1.;
	map<TString, double> dDiff = {
	{"ptjet1_q2_4_6", 2},
	{"ptjet1_q2_6_10", 4},
	{"ptjet1_q2_10_18",8},
	{"ptjet1_q2_18_34",16},
	{"ptjet1_q2_34_100", 66},
	{"zpom_q2_4_10", 6},
	{"zpom_q2_10_20", 10},
	{"zpom_q2_20_40", 20},
	{"zpom_q2_40_100", 60},

	{"zpom_ptjet1_5_6p5", 1.5},
	{"zpom_ptjet1_6p5_8", 1.5},
	{"zpom_ptjet1_8_16",  8},
	{"zpom_q2_5_12",   7},
	{"zpom_q2_12_25", 13},
	{"zpom_q2_25_50", 25},
	{"zpom_q2_50_100",50},
	};

	if( dDiff.count(var_name)>0)
		normFactor = 1./dDiff[var_name];


	//say::SetGlobalVerbosity(say::DEBUG);
	TString TheoryFile = tablesDir +theor_path+ theor_file;
	//TheoryFile = "../tables/nnlojet/FPS/nlo/H1-LQall-7.diff_FPS_q2.abs.tab";
	TheoryFile.ReplaceAll(":", "/");
	cout <<"Theory file is " <<  TheoryFile << endl;

	vector<fastNLODiffAlphas> fnlodiffs;

	fnlodiffs.push_back( fastNLODiffAlphas(TheoryFile.Data()) );

	if(var_name == "ptJ" || var_name == "etaJ") {
		TString TheoryFileJet1 = TheoryFile;
		TheoryFileJet1.ReplaceAll("ptj2", "ptj1");
		TheoryFileJet1.ReplaceAll("etaj2", "etaj1");
		cout << "Add theor file is " << TheoryFileJet1 << endl;
		fnlodiffs.push_back( fastNLODiffAlphas(TheoryFileJet1.Data()) );
	}


	//fastNLODiffAlphas fnlodiff( TheoryFile.Data() );
	int Order = 0;

	for(auto & fnlodiff : fnlodiffs) {
		fnlodiff.SetMuRFunctionalForm(fastNLO::kQuadraticMean);
		fnlodiff.SetMuFFunctionalForm(fastNLO::kQuadraticMean);

		if(theor_file.Contains("nnlo:")) {
			fnlodiff.SetContributionON(fastNLO::kFixedOrder,0,true);
			fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,true);
			fnlodiff.SetContributionON(fastNLO::kFixedOrder,2,true);
			cout << "Is nnlo" << endl;
			Order = 3;
		}
		else if(theor_file.Contains("nlo:")) {
			fnlodiff.SetContributionON(fastNLO::kFixedOrder,0,true);
			fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,true);
			fnlodiff.SetContributionON(fastNLO::kFixedOrder,2,false);
			cout << "Is nlo" << endl;
			Order = 2;
		}
		else {
			fnlodiff.SetContributionON(fastNLO::kFixedOrder,0,true);
			fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,false);
			fnlodiff.SetContributionON(fastNLO::kFixedOrder,2,false);
			cout << "Is lo" << endl;
			Order = 1;
		}
		//set other parameters according to Setting


        

		if(setting.scaleTag == "Q2pPt2") {
			fnlodiff.SetExternalFuncForMuR(funQ2pPt2); 
			fnlodiff.SetExternalFuncForMuF(funQ2pPt2);
		}
		else if(setting.scaleTag == "Q2") {
			fnlodiff.SetExternalFuncForMuR(funQ2);
			fnlodiff.SetExternalFuncForMuF(funQ2);
		}
		else if(setting.scaleTag == "Pt2") {
			fnlodiff.SetExternalFuncForMuR(funPt2);
			fnlodiff.SetExternalFuncForMuF(funPt2);
		}
		else if(setting.scaleTag == "0.25Q2pPt2") {
			fnlodiff.SetExternalFuncForMuR(funForthQ2pPt2);
			fnlodiff.SetExternalFuncForMuF(funForthQ2pPt2);
		}
		else if(setting.scaleTag == "SqrtQ4pPt4") {
			fnlodiff.SetExternalFuncForMuR(funSqrtQ4pPt4);
			fnlodiff.SetExternalFuncForMuF(funSqrtQ4pPt4);
		}
		else if(setting.scaleTag == "ZEUSscale") {
			fnlodiff.SetExternalFuncForMuR(funZEUSmuR);
			fnlodiff.SetExternalFuncForMuF(funZEUSmuF);
		}
		else {
			cout << "Unknown function tag" << endl;
			assert(0);
		}



		fnlodiff.SetScaleFactorsMuRMuF(setting.renFactor, setting.factFactor);


		// Set the xpom integration interval and method
		// -------- Boris LRG dijets
		if( theor_path.Contains("/VFPS/") ) {
			fnlodiff.SetXPomLinSlicing( 20, 0.010 ,  .024 ); // VFPS range
			fnlodiff.SettIntegratedRange(-0.6);
			fnlodiff.SetProtonE(920.);
		}
		else if( theor_path.Contains("/FPS/") ) {
			fnlodiff.SetXPomLinSlicing( 20, 0.00 ,  0.10 ); // FPS range
			fnlodiff.SettIntegratedRange(-1.);
			fnlodiff.SetProtonE(920.);
		}
		else if( theor_path.Contains("/LRGH1820/") ) {
			fnlodiff.SetXPomLinSlicing( 20, 0. ,  0.03 ); // LRG range
			fnlodiff.SettIntegratedRange(-1.);
			fnlodiff.SetProtonE(820.);
			normFactor *= pow(920./820.,2);//TODO check
		}
		else {
			fnlodiff.SetXPomLinSlicing( 20, 0. ,  0.03 ); // LRG range
			fnlodiff.SettIntegratedRange(-1.);
			fnlodiff.SetProtonE(920.);
		}


		if(setting.pdfName == "FitB")
			fnlodiff.SetFit(fastNLODiffAlphas::FitB);
		else if(setting.pdfName == "FitA")
			fnlodiff.SetFit(fastNLODiffAlphas::FitA);
		else if(setting.pdfName == "FitJets")
			fnlodiff.SetFit(fastNLODiffAlphas::FitJets);
		else if(setting.pdfName == "Fit19")
			fnlodiff.SetFit(fastNLODiffAlphas::Fit19);
		else if(setting.pdfName == "zeusSJ")
			fnlodiff.SetFit(fastNLODiffAlphas::zeusSJ);
		else if(setting.pdfName == "MRW")
			fnlodiff.SetFit(fastNLODiffAlphas::MRW);
		else
			assert(0 && "Wrong fit name");

		fnlodiff.SetLHAPDFMember(setting.pdfVecId);

        fnlodiff.IncludeOnlyQuarks(setting.onlyQuarks);







	
	}

	//reference to first entry
	auto & fnlodiff = fnlodiffs[0];

	
	//cout << "long " << sizeof(int) << " "<< sizeof(long)<<" "<<sizeof(long long) << endl;

	//cout<<"DANIEL " << fastNLO::kExtern<<endl;
	//exit(1);

	/*
	//TEST PART
	//  If you want to receive your cross section in
	//   pb/GeV or in pb. Here we choose pb/GeV
	fnlodiff.SetUnits(fastNLO::kPublicationUnits);


	fnlodiff.SetExternalFuncForMuR (&Function_Mu);
	fnlodiff.SetExternalFuncForMuF (&Function_Mu);
	//fnlodiff.SetExternalFuncForMuR (&Function_MuDISENTren);
	//fnlodiff.SetExternalFuncForMuF (&Function_MuDISENTfac);


	//fastBinsLo = fnlodiff.GetObsBinsLoBounds(0);
	//fastBinsHi = fnlodiff.GetObsBinsUpBounds(0);

	// calculate and access the cross section

	fnlodiff.SetScaleFactorsMuRMuF(1.0, 1.0);
	vector<double> myXsc = fnlodiff.GetDiffCrossSection();

	for(auto x : myXsc)
		cout << x << "   ";
	cout << endl;

	cout << "Helenka " << dissFactor << " "<< normFactor <<" "<< hadrCorr[0]<< endl;
	exit(1);
	*/

	//TEST PART end




	//fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,false);

	vector<double> fastBinsLo, fastBinsHi ;

	auto MX = [] (double xpom, double Q2, double y) {
	              const double s = 4*27.6*920.;
	              return  sqrt( max(0.0,y * s * xpom - Q2)); };

	auto ZP = [] (double xpom, double Q2, double xi12) {
	              return xi12/xpom; };

	auto BETA= [] (double xpom, double Q2, double xBjor) {
	              return xBjor/xpom; };

	Theory th;
	th.iOrder = Order; //setting order of theory
	th.file_name = theor_file;

	//cout << "RADEKHERE " << __LINE__ << endl;
	if(var_name.Contains("xpom") )
		th = CalculateXpomNLO( fnlodiff,  fastBinsLo,  fastBinsHi );
	else if(var_name == "mx" )
		th = CalculateVarNLO( fnlodiff,  fastBinsLo,  fastBinsHi, MX );
	else if(var_name.BeginsWith("zpom") )
		th = CalculateVarNLO( fnlodiff,  fastBinsLo,  fastBinsHi, ZP );
	else if(var_name == "beta" )
		th = CalculateVarNLO( fnlodiff,  fastBinsLo,  fastBinsHi, BETA );
	else if(var_name == "ptJ"|| var_name == "etaJ") {
		assert(fnlodiffs.size() == 2);
		th  = CalculateRawNLO( fnlodiffs[0],  fastBinsLo,  fastBinsHi  );
		th += CalculateRawNLO( fnlodiffs[1],  fastBinsLo,  fastBinsHi  );
		if(var_name == "etaJ") { //Correct for -sign convention
			th.ReverseBinning();
			reverse(begin(fastBinsLo), end(fastBinsLo) );
			reverse(begin(fastBinsHi), end(fastBinsHi) );
			for(unsigned i = 0; i < fastBinsLo.size(); ++i) {
				fastBinsLo[i] *= -1;
				fastBinsHi[i] *= -1;
			}
			swap(fastBinsLo, fastBinsHi);
		}
	}
	else
		th = CalculateRawNLO( fnlodiff,  fastBinsLo,  fastBinsHi  );



	//fnlodiff.SetMuFFunctionalForm(kScale2);
	//fnlodiff.SetMuRFunctionalForm(kScale2);

	for(unsigned i = 0; i < th.xsc.size(); ++i) {
		th.xsc[i] *= dissFactor * normFactor* hadrCorr[i];

		//th.xscUp[i] *= dissFactor * normFactor * hadrCorr[i];
		//th.xscDown[i] *= dissFactor * normFactor * hadrCorr[i];

		//for(int j=0; j <= npdfall; ++j)
			//th.xscSyst[j][i] *= dissFactor * normFactor * hadrCorr[i];

	}



	try {
		if( static_cast<int>(fastBinsLo.size()) != nBins ||
		    static_cast<int>(fastBinsHi.size()) != nBins ) {
			throw 20;
		}
		for(int i = 0; i < nBins; ++i) {
			if( abs( fastBinsLo[i] - xMin[i] ) > 1e-4 )
				throw 21;
			if( abs( fastBinsHi[i] - xMax[i] ) > 1e-4 )
				throw 21;
			
		}
	}
	catch ( int e ) {
		cout << "Bins numbers do no corresponds for variable " << var_name << " "<<  data_file << endl;
		cout << fastBinsLo.size() << " "<< nBins << endl;
		cout << "Bins in data:" << endl;
		for(int i = 0; i < nBins; ++i)
			cout << xMin[i]<<","<<xMax[i]<<":";
		cout <<endl<< "Bins in Theory:" << endl;
		for(unsigned i = 0; i < fastBinsLo.size(); ++i)
			cout << fastBinsLo[i]<<","<<fastBinsHi[i]<<":";

		cout << "Exception " << e << endl;
		exit(1);
	}

	th.setting = setting;

	return th;

}



Theory Histogram::LoadTheory(TString theor_file, Setting setting)
{
	if( theor_path.Contains("/VFPS/") ) {
		dissFactor = 1/1.2;
	}
	else if( theor_path.Contains("/FPS/") ) {
		dissFactor = 1/1.2;
	}
	else if( theor_path.Contains("/LRGZEUS/") ) {
		dissFactor = 1/1.2;
	}
	else {
		dissFactor = 1.;
	}

	TString TextFile = theor_file(theor_file.First(':')+1, 100000); 
	TextFile = TextFile(0, TextFile.Last('.') ); 
	
	TString orderTag = theor_file(0,theor_file.First(':'));
	TString rawFileN = theor_file(theor_file.First(':')+1, 10000);


	TString fileName = TString("../data/Theory/") + theor_path + rawFileN;//  theor_file; 
	fileName.ReplaceAll(":", "/");
	TString setTag = setting.getFileName();
	fileName = fileName(0, fileName.Last('/') )+ "/" +var_name+"/"+orderTag+"/"+ TextFile + "-" +setTag+"-"+ var_name + ".txt";
	TString fileDirectory=fileName(0, fileName.Last('/') );

	//cout << fileDirectory << endl;
	//cout << fileName << endl;
	//exit(1);
    //cout << "KAREL " << theor_file << endl;

	Theory th;

	std::ifstream nloXsec(fileName);
	if( nloXsec.good() ) {
		//Load theory
		th.xsc.resize(nBins);
		//th.xscUp.resize(nBins);
		//th.xscDown.resize(nBins);
		//th.xscSyst.resize(npdfall+1);
		//for(auto &a : th.xscSyst)
			//a.resize(nBins);

		try {
			for(int i = 0; i < nBins; ++i)
				nloXsec >> th.xsc[i];
			/*
			for(int i = 0; i < nBins; ++i)
				nloXsec >> th.xscUp[i];
			for(int i = 0; i < nBins; ++i)
				nloXsec >> th.xscDown[i];
			for(int j = 0; j <= npdfall; ++j) {
				for(int i = 0; i < nBins; ++i)
					nloXsec >> th.xscSyst[j][i];
			}
			*/
		}
		catch (std::ifstream::failure e) {
			std::cerr << "Exception opening/reading/closing file\n";
		}

		nloXsec.close();
		th.setting = setting;
	}
	else {
		th = CalcNLO(theor_file, setting);


		//Save theory
		std::ofstream nloXsecNew(fileName);
		if(!nloXsecNew.good()) {
			int st = system( (TString("mkdir -p ")+fileDirectory).Data() );
			if(st != 0) {
				cout << "Troubles with file creating " << endl;
				cout << "Status " << st << endl;
				exit(1);
			}
			nloXsecNew.open(fileName, std::ofstream::out);
			if(!nloXsecNew.good()) assert(0);
		}

		for(int i = 0; i < nBins; ++i)
			nloXsecNew << th.xsc[i] << " ";
		/*
		nloXsecNew << endl;
		for(int i = 0; i < nBins; ++i)
			nloXsecNew << th.xscUp[i] << " ";
		nloXsecNew << endl;
		for(int i = 0; i < nBins; ++i)
			nloXsecNew << th.xscDown[i] << " ";
		nloXsecNew << endl;

		for(int j = 0; j <= npdfall; ++j) {
			for(int i = 0; i < nBins; ++i)
				nloXsecNew << th.xscSyst[j][i] << " ";
			nloXsecNew << endl;
		}
		*/

		nloXsecNew.close();
		cout <<"Theory saved to " << fileName << endl;
	}

	return th;

	//fastNLODiffAlphas fnlodiff( (TString("../tables/") + nlo_file).Data() );

}

pair<double,double> GetSystTot(vector<double> *sysVec, int binId) 
{
	double sum2Up = 0;
	double sum2Dn = 0;

	for(unsigned i = 1; i <= 30; ++i) {
		double d = sysVec[i][binId] - sysVec[0][binId];

        sum2Up += pow(max(0.0, d), 2);
        sum2Dn += pow(max(0.0,-d), 2);

		//sum2 += d*d;
	}
	//return sqrt(sum2/2.);// TODO, dummy now
    return make_pair(sqrt(sum2Up), sqrt(sum2Dn));
}


vector<double> Histogram::GetBinning()
{
	int nBins = xMin.size();
	vector<double> bins(nBins+1);
	for(int i = 0; i < nBins; ++i)
		bins[i] = xMin[i];
	bins[nBins] = xMax[nBins-1];
	return bins;
};

vector<double> Hist2Vec(TH1D *h)
{
    vector<double> vals;
    for(int i = 1; i <= h->GetNbinsX(); ++i)
        vals.push_back(h->GetBinContent(i));
    return vals;
}

struct HistoErr {
	TH1D *h, *hRatio;

	TGraphAsymmErrors *gr, *grAll;
	TGraphAsymmErrors *grRatio, *grAllRatio;

	TString name;
	vector<double> hBins;
	int nBins;

	void Init(TString Name, vector<double> &xMin, vector<double> &xMax) {
		nBins = xMin.size();
		name = Name;

		hBins.resize(nBins+1);
		for(int i = 0; i < nBins; ++i)
			hBins[i] = xMin[i];
		hBins[nBins] = xMax[nBins-1];

		int hash = rand();

		h = new TH1D( SF("h%s%d",name.Data(), hash), name, nBins, hBins.data());

		hRatio = new TH1D( SF("h%sRatio%d",name.Data(),hash), SF("%sRatio",name.Data()), nBins, hBins.data());

		gr      = new TGraphAsymmErrors(nBins);
		grAll   = new TGraphAsymmErrors(nBins);
		grRatio = new TGraphAsymmErrors(nBins);
		grAllRatio = new TGraphAsymmErrors(nBins);
	}

	void SetBin(int i, double val, double valRef,
	            double err0Down, double err0Up, double err1Down, double err1Up) {
		
		double valX = h->GetBinCenter(i);
		double binSize = h->GetBinWidth(i) / 2.;
		if(name.Contains("Data")) binSize = 0;

		h->SetBinContent(i, val);
		h->SetBinError(i, 0);

		gr->SetPoint(i-1, valX, val);
		//cout << "Name " << name << " : "<< binSize << endl;
		gr->SetPointError(i-1, binSize, binSize, err0Down, err0Up);

		grAll->SetPoint(i-1, valX, val);
		double totErUp   = hypot(err0Up, err1Up);
		double totErDown = hypot(err0Down, err1Down);
		grAll->SetPointError(i-1, binSize, binSize,  totErDown, totErUp );

		hRatio->SetBinContent(i, val/valRef );
		hRatio->SetBinError(i, 0.0 );

		grRatio->SetPoint(i-1, valX, val/valRef );
		grRatio->SetPointError(i-1, binSize, binSize, err0Down/valRef,  err0Up/valRef );

		grAllRatio->SetPoint(i-1, valX, val/valRef );
		grAllRatio->SetPointError(i-1, binSize, binSize,  totErDown/valRef, totErUp/valRef );


	}
	double GetMaxY(int i) {
		return h->GetBinContent(i+1) + grAll->GetErrorYhigh(i);
	}
	double GetMinY(int i) {
		return h->GetBinContent(i+1) - grAll->GetErrorYlow(i);
	}

	double GetMaxYratio(int i) {
		return hRatio->GetBinContent(i+1) + grAllRatio->GetErrorYhigh(i);
	}
	double GetMinYratio(int i) {
		return hRatio->GetBinContent(i+1) - grAllRatio->GetErrorYlow(i);
	}


};

HistoErr Merge(TString name, vector<HistoErr> histos)
{
	HistoErr sum;
	vector<double> xMin, xMax;
	for(unsigned i = 0; i < histos.size(); ++i) {
		xMin.push_back(i+0.5);
		xMax.push_back(i+1.5);
	}

	sum.Init(name, xMin, xMax);
	for(unsigned i = 0; i < histos.size(); ++i) {

		sum.h->SetBinContent(i+1, histos[i].h->GetBinContent(1));
		sum.h->SetBinError(i+1,   histos[i].h->GetBinError(1));
		sum.hRatio->SetBinContent(i+1, histos[i].hRatio->GetBinContent(1));
		sum.hRatio->SetBinError(i+1,   histos[i].hRatio->GetBinError(1));

		double x, y, yRat;
		histos[i].gr->GetPoint(0, x, y);
		histos[i].grRatio->GetPoint(0, x, yRat);

		sum.gr->SetPoint(i, i+1., y);
		sum.grAll->SetPoint(i, i+1., y);
		sum.grRatio->SetPoint(i, i+1., yRat);
		sum.grAllRatio->SetPoint(i, i+1., yRat);


		auto CopyErr = [](int i, TGraphAsymmErrors *grIn, TGraphAsymmErrors *grOut) {
			double xh = grIn->GetErrorXhigh(0);
			double xl = grIn->GetErrorXlow(0);
			double yh = grIn->GetErrorYhigh(0);
			double yl = grIn->GetErrorYlow(0);
			grOut->SetPointError(i, xl, xh, yl, yh);
		};
		CopyErr(i, histos[i].gr, sum.gr);
		CopyErr(i, histos[i].grAll, sum.grAll);
		CopyErr(i, histos[i].grRatio, sum.grRatio);
		CopyErr(i, histos[i].grAllRatio, sum.grAllRatio);


	}
	return sum;

}


vector<double> Histogram::LoadHistogramsNloNnloData(TString refName, HistoErr &grData, HistoErr &grNlo, HistoErr &grNnlo, bool sysInside)
{

	TString refTag  = refName(0, refName.First(':')+1);
	TString refEnd  = refName(refName.First(':')+1, 10000);

	vector<double>  refXsc;

	vector<double>  Nlo,  NloUp, NloDown;
	vector<double>  NNlo,  NNloUp, NNloDown;
	vector<double> NloSyst[32];
	vector<double> NNloSyst[32];

	//vector<TString> NloNames, NNloNames;
    int isFound = 0;
	for(auto & th : theories) {
		TString name = th.first;
        //cout << "Pusinka " << name << endl;
		if( name.Contains(defFile) ) {
			if(name.BeginsWith(refTag) && name.EndsWith(refEnd))
				refXsc = th.second.xsc;

			if( name.BeginsWith("nlo:") ) {
                ++isFound;
				//cout << "RADEK " << name << endl;
				if(name.EndsWith("Fit19-0-Q2pPt2-cc"))
					Nlo = th.second.xsc ;
				if(name.EndsWith("Fit19-0-Q2pPt2-uu"))
					NloUp= th.second.xsc;
				if(name.EndsWith("Fit19-0-Q2pPt2-dd"))
					NloDown =  th.second.xsc;

				for(int s = 0; s <= 30; ++s) 
					if(name.EndsWith(SF("Fit19-%d-Q2pPt2-cc",s)))
						NloSyst[s] = th.second.xsc ;

			}
			if( name.BeginsWith("nnlo:")  ) {
                ++isFound;
				if(name.EndsWith("Fit19-0-Q2pPt2-cc"))
					NNlo= th.second.xsc ;
				if(name.EndsWith("Fit19-0-Q2pPt2-uu"))
					NNloUp =  th.second.xsc;
				if(name.EndsWith("Fit19-0-Q2pPt2-dd"))
					NNloDown =  th.second.xsc;

				for(int s = 0; s <= 30; ++s) 
					if(name.EndsWith(SF("Fit19-%d-Q2pPt2-cc",s)))
						NNloSyst[s] = th.second.xsc ;

			}
		}
	}
    if(isFound<2) {
        cout << "Nothing found for "<< refName << endl;
        exit(0);
    }
    //cout << "Test "<<NloUp.size()<<" "<< Nlo.size()<<" "<<NloDown.size() << endl;


	grData.Init("Data", xMin, xMax);
	grNlo.Init("NLO", xMin, xMax);
	grNnlo.Init("NNLO", xMin, xMax);

	double dispMax = -1e30;
	double dispMin = +1e30;
	double ratMax = -100;
	double ratMin = +100;

	for(unsigned i = 0; i < xMin.size(); ++i) {

		grData.SetBin(i+1, data[i], refXsc[i], dataStatErr[i], dataStatErr[i],
		                                   dataSystErr[i], dataSystErr[i]);

		//LOAD NLO
		//double nloSyst = 0.;//GetSystTot(NloSyst,  i);

        double nloSystUp=0, nloSystDn = 0;
		//tie(nloSystUp,nloSystDn) = GetSystTot(NloSyst, i); //Comment it for plotting

		double nloErUp    = max({0., NloUp[i] - Nlo[i], NloDown[i] - Nlo[i]});
		double nloErDown  =-min({0., NloUp[i] - Nlo[i], NloDown[i] - Nlo[i]});

		if(sysInside == false)
			grNlo.SetBin(i+1, Nlo[i], refXsc[i], nloErDown, nloErUp,
			                                   nloSystDn, nloSystUp);
		else
			grNlo.SetBin(i+1, Nlo[i], refXsc[i], nloSystDn, nloSystUp,
			                               nloErDown, nloErUp);

		//LOAD NNLO
        double nnloSystUp=0, nnloSystDn = 0;
		//tie(nnloSystUp,nnloSystDn) = GetSystTot(NNloSyst, i); //without systematics
		double nnloErUp    = max({0., NNloUp[i] - NNlo[i], NNloDown[i] - NNlo[i]});
		double nnloErDown  =-min({0., NNloUp[i] - NNlo[i], NNloDown[i] - NNlo[i]});

		if(sysInside == false)
			grNnlo.SetBin(i+1, NNlo[i], refXsc[i], nnloErDown, nnloErUp,
			                                       nnloSystDn, nnloSystUp);
		else
			grNnlo.SetBin(i+1, NNlo[i], refXsc[i], nnloSystDn, nnloSystUp,
			                                       nnloErDown, nnloErUp);


		dispMax = max({dispMax,grData.GetMaxY(i), grNlo.GetMaxY(i), grNnlo.GetMaxY(i) });
		dispMin = min({dispMin,grData.GetMinY(i), grNlo.GetMinY(i), grNnlo.GetMinY(i) });

		ratMax = max({ratMax,grData.GetMaxYratio(i), grNlo.GetMaxYratio(i), grNnlo.GetMaxYratio(i) });
		ratMin = min({ratMin,grData.GetMinYratio(i), grNlo.GetMinYratio(i), grNnlo.GetMinYratio(i) });

	}
	return {dispMin, dispMax, ratMin, ratMax};

}

HistoErr Histogram::LoadHistograms(TString name, TString file, TString histName, TString refName)
{
	vector<double> ref, xSec;

	TString histTag = histName(0, histName.First(':')+1);
	TString refTag  = refName(0, refName.First(':')+1);

	TString histEnd = histName(histName.First(':')+1, 10000);
	TString refEnd  = refName(refName.First(':')+1, 10000);
	//cout << "Helenka " << histName << endl;
	//cout << "RADECEK " << refTag << "         " << refEnd << endl;

	for(auto & th : theories) {
		TString name = th.first;
		if( name.Contains(file) ) {
			if( name.BeginsWith(refTag) && name.EndsWith(refEnd) ) 
				ref = th.second.xsc;
			if( name.BeginsWith(histTag) && name.EndsWith(histEnd) ) 
				xSec = th.second.xsc;
		}
	}
	assert(ref.size() > 0);
	assert(xSec.size() > 0);

	HistoErr hist;
	hist.Init(name, xMin, xMax);

	for(unsigned i = 0; i < xMin.size(); ++i) 
		hist.SetBin(i+1, xSec[i], ref[i], 0, 0, 0, 0);

	return hist;
}


void Histogram::LoadHistogramsFitJetsSJ(HistoErr &grNloJets, HistoErr &grNnloJets,
                                                  HistoErr &grNloSJ, HistoErr &grNnloSJ)
{

	vector<double> Nlo;
	vector<double>  NloJets, NloSJ;
	vector<double>  NNloJets, NNloSJ;

	for(auto & th : theories) {
		TString name = th.first;

		if( name.Contains(defFile) ) {
			if( name.BeginsWith("nlo:") ) {
				//cout << "RADEK " << name << endl;
				if(name.EndsWith("FitJets-0-Q2pPt2-cc"))
					NloJets = th.second.xsc ;
				if(name.EndsWith("zeusSJ-0-Q2pPt2-cc"))
					NloSJ = th.second.xsc ;
				if(name.EndsWith("FitB-0-Q2pPt2-cc"))
					Nlo = th.second.xsc ;
			}
			if( name.BeginsWith("nnlo:")  ) {

				if(name.EndsWith("FitJets-0-Q2pPt2-cc"))
					NNloJets= th.second.xsc ;
				if(name.EndsWith("zeusSJ-0-Q2pPt2-cc"))
					NNloSJ= th.second.xsc ;
			}
		}
	}

	grNloSJ.Init("nloSJ", xMin, xMax);
	grNnloSJ.Init("nnloSJ", xMin, xMax);
	grNloJets.Init("nloJets", xMin, xMax);
	grNnloJets.Init("nnloJets", xMin, xMax);

	//double dispMax = -1e30;
	//double dispMin = +1e30;
	//double ratMax = -100;
	//double ratMin = +100;

	for(unsigned i = 0; i < xMin.size(); ++i) {


		grNloSJ.SetBin(i+1, NloSJ[i], Nlo[i], 0, 0, 0, 0);
		grNnloSJ.SetBin(i+1, NNloSJ[i], Nlo[i], 0, 0, 0, 0);
		grNloJets.SetBin(i+1, NloJets[i], Nlo[i], 0, 0, 0, 0);
		grNnloJets.SetBin(i+1, NNloJets[i], Nlo[i], 0, 0, 0, 0);

		//dispMax = max({dispMax,grData.GetMaxY(i), grNlo.GetMaxY(i), grNnlo.GetMaxY(i) });
		//dispMin = min({dispMin,grData.GetMinY(i), grNlo.GetMinY(i), grNnlo.GetMinY(i) });

		//ratMax = max({ratMax,grData.GetMaxYratio(i), grNlo.GetMaxYratio(i), grNnlo.GetMaxYratio(i) });
		//ratMin = min({ratMin,grData.GetMinYratio(i), grNlo.GetMinYratio(i), grNnlo.GetMinYratio(i) });

	}
	//return {dispMin, dispMax, ratMin, ratMax};

}












void CopyStyle(TGraphAsymmErrors *grOrg, TGraphAsymmErrors *grRat)
{
	grRat->SetLineColor( grOrg->GetLineColor() );
	grRat->SetLineWidth( grOrg->GetLineWidth() );
	grRat->SetFillColor( grOrg->GetFillColor());
	grRat->SetFillStyle( grOrg->GetFillStyle() );

	grRat->SetMarkerStyle(grOrg->GetMarkerStyle());
	grRat->SetMarkerSize(grOrg->GetMarkerSize());
}

void CopyStyle(TH1D *hOrg, TH1D *hRat)
{
	hRat->SetLineColor(hOrg->GetLineColor());
	hRat->SetLineStyle(hOrg->GetLineStyle());
	hRat->SetLineWidth(hOrg->GetLineWidth());
	hRat->SetMarkerColor(hOrg->GetMarkerColor());
}




#if 0
void Histogram::plotNLOvsNNLO()
{


	HistoErr grData, grNlo, grNnlo;
	vector<double> range = LoadHistogramsNloNnloData("nlo:FitB-0-Q2pPt2-cc", grData, grNlo, grNnlo);

	/*
	HistoErr grNloSJ = LoadHistograms("sjNLO", "H1-LQall-8c.diff",
	                         "nlo:zeusSJ-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloSJ = LoadHistograms("sjNNLO", "H1-LQall-8c.diff",
	                         "nnlo:zeusSJ-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNloJets = LoadHistograms("NLOjets", "H1-LQall-8c.diff",
	                         "nlo:FitJets-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloJets = LoadHistograms("NNLOjets", "H1-LQall-8c.diff",
	                         "nnlo:FitJets-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	*/


	double leftMargin = 0.15;
	double rightMargin = 0.03;
	double bottomMargin = 0.35;
	double fontS = 0.08;

	double pointSize = 0.8;

	TVirtualPad * currpad = 0;

	currpad = gPad;//canvas->cd(1);

	//currpad -> Divide(1,3,0.,0.,0); // rozdelim prvni pad na tri podsebou
	//currpad -> Divide(1,2,0.,0.,0);

	// a ty se zacnou znovu cislovat od jedne

	//Setting style NLO
	grNlo.grAll->SetLineColor(kWhite);
	grNlo.grAll->SetFillColor(kYellow);
	CopyStyle(grNlo.grAll, grNlo.grAllRatio);

	grNlo.gr->SetLineColor(kWhite);
	grNlo.gr->SetFillColor(kOrange);
	CopyStyle(grNlo.gr, grNlo.grRatio);

	grNlo.h->SetLineColor(kWhite);
	CopyStyle(grNlo.h, grNlo.hRatio);


	//Setting style NNLO
	grNnlo.grAll->SetLineColor(kGreen);
	grNnlo.grAll->SetFillColor(kGreen+3);
	grNnlo.grAll->SetFillStyle(3004);
	CopyStyle(grNnlo.grAll, grNnlo.grAllRatio);

	grNnlo.gr->SetLineColor(kGreen);
	grNnlo.gr->SetFillColor(kGreen);
	grNnlo.gr->SetFillStyle(3004);
	CopyStyle(grNnlo.gr, grNnlo.grRatio);

	grNnlo.h->SetLineColor(kGreen);
	grNnlo.h->SetLineStyle(1);
	CopyStyle(grNnlo.h, grNnlo.hRatio);


	//Setting style Data
	grData.gr->SetMarkerStyle(20);
	grData.gr->SetMarkerSize(pointSize);
	CopyStyle(grData.gr, grData.grRatio);

	grData.grAll->SetMarkerStyle(20);
	grData.grAll->SetMarkerSize(pointSize);
	CopyStyle(grData.grAll, grData.grAllRatio);



	//currpad -> cd(1);// a prepnu se do prvniho padu 
	int hash = rand();
	TPad *upPad = new TPad( SF("upPad%d",hash), "upPad", 0, 0.5, 1, 1 );
	upPad->SetBottomMargin(0);
	upPad->SetLeftMargin(leftMargin);
	upPad->SetRightMargin(rightMargin);
	upPad->Draw();
	upPad->cd();

	if(style.isLogY) gPad->SetLogy();
	if(style.isLogX) gPad->SetLogx();

	//gPad -> SetRightMargin(0.02);

	grData.h->SetLineColor(kBlack);
	grData.h->Draw("AXIS");
	grData.h->GetYaxis()->SetTitle(yTitle);
	grData.h->GetYaxis()->SetTitleSize(fontS);
	grData.h->GetYaxis()->SetTitleOffset(0.85);
	grData.h->GetYaxis()->SetLabelSize(fontS);
	if(!style.isLogY) {
		grData.h->SetMinimum(0);
		grData.h->SetMaximum(range[1] * 1.1);
	}
	else {
		grData.h->SetMinimum(range[0] * 0.8);
		grData.h->SetMaximum(range[1] * 1.2);
	}


	grNlo.grAll->Draw("e2 same");
	grNlo.gr->Draw("e2 same");

	grNlo.h->Draw("same");

	grNnlo.grAll->Draw("e2 same");
	grNnlo.gr->Draw("e2 same");

	grNnlo.h->Draw("same");

	//grNloJets.h->Draw("same");
	//grNloSJ.h->Draw("same");
	//grNnloJets.h->SetLineColor(kBlue);
	//grNnloJets.h->Draw("same");
	//grNnloSJ.h->SetLineColor(kRed);
	//grNnloSJ.h->Draw("same");

	gStyle->SetEndErrorSize(4);
	grData.gr->Draw("e same");
	grData.grAll->Draw("pz same");



	//Plot Legend
	TLegend *leg = new TLegend( style.legX+0.5-0.5, style.legY+0.60-0.6,
	                            style.legX+0.95-0.5,style.legY+0.88-0.6 );
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	TString tag =  theor_path(theor_path.First('/')+1, theor_path.Last('/')-theor_path.First('/')-1 );
	leg->SetHeader(Title);
	leg->AddEntry(grData.gr,    SF("H1 %s data", tag.Data()), "ep");

	leg->AddEntry(grNlo.grAll,     SF("NLO x %.2f",dissFactor) , "fl");
	leg->AddEntry(grNnlo.grAll,   SF("NNLO  x %.2f", dissFactor) , "fl");


	leg->Draw();

	gPad->Update();
	gPad->RedrawAxis();
	TFrame *frU = (TFrame *) gPad->FindObject("TFrame");
	frU->SetFillStyle(0);
	frU->Draw();
	gPad->Update();


	currpad->cd();
	TPad *downPad = new TPad( SF("downPad%d",hash), "downPad", 0, 0.0, 1, 0.5 );
	downPad->SetTopMargin(0);
	downPad->SetLeftMargin(leftMargin);
	downPad->SetRightMargin(rightMargin);
	downPad->SetBottomMargin(bottomMargin);
	downPad->Draw();
	downPad->cd();
	//gPad -> SetRightMargin(0.02);
	if(style.isLogX) gPad->SetLogx();

	//double Max = max( hDataRatio->GetMaximum(), hNLORatio->GetMaximum() );
	//grData.hRatio->SetMaximum(Max);

	grData.hRatio->GetYaxis()->SetNdivisions(204);

	grData.hRatio->Draw("AXIS");
	grData.hRatio->SetMinimum( max(0.,range[2] - 0.1*range[3]) );
	grData.hRatio->SetMaximum( 1.1*range[3] );


	grData.hRatio->GetXaxis()->SetTitle(xTitle);
	grData.hRatio->GetXaxis()->SetTitleSize(fontS * 5./5);

	grData.hRatio->GetYaxis()->SetLabelSize(fontS * 5./5);
	grData.hRatio->GetXaxis()->SetLabelSize(fontS * 5./5);

	grData.hRatio->GetYaxis()->SetTitle("Data/NLO");
	grData.hRatio->GetYaxis()->SetTitleSize(fontS * 5./5);
	grData.hRatio->GetYaxis()->SetTitleOffset(0.84);

	grNlo.grAllRatio->Draw("e2 same");

	grNlo.grRatio->Draw("e2 same");


	grNlo.hRatio->Draw("same");


	grNnlo.grAllRatio->Draw("e2 same");
	grNnlo.grRatio->Draw("e2 same");

	grNnlo.hRatio->Draw("same");


	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");

	gPad->Update();
	gPad->RedrawAxis();
	TFrame *frD = (TFrame *) gPad->FindObject("TFrame");
	frD->SetFillStyle(0);
	frD->Draw();
	gPad->Update();
}
#endif

#if 0
void Histogram::plotScaleStudies()
{


	HistoErr grData, grNlo, grNnlo;
	vector<double> range = LoadHistogramsNloNnloData("nnlo:FitB-0-Q2pPt2-cc", grData, grNlo, grNnlo);


	HistoErr grNNloQ2  = LoadHistograms("NNLOq2", "H1-LQall-8c.diff",
	                         "nnlo:FitB-0-Q2-cc", "nnlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloPt = LoadHistograms("NNLOpt", "H1-LQall-8c.diff",
	                         "nnlo:FitB-0-Pt2-cc", "nnlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloQ2Pt = LoadHistograms("NNLOq2pt", "H1-LQall-8c.diff",
	                         "nnlo:FitB-0-0.25Q2pPt2-cc", "nnlo:FitB-0-Q2pPt2-cc");



	double leftMargin = 0.15;
	double rightMargin = 0.03;
	double bottomMargin = 0.35;
	double fontS = 0.08;

	double pointSize = 0.8;

	TVirtualPad * currpad = 0;

	currpad = gPad;//canvas->cd(1);

	//currpad -> Divide(1,3,0.,0.,0); // rozdelim prvni pad na tri podsebou
	//currpad -> Divide(1,2,0.,0.,0);

	// a ty se zacnou znovu cislovat od jedne

	//Setting style NNLO
	grNnlo.grAll->SetLineColor(kGreen);
	grNnlo.grAll->SetFillColor(kGreen+3);
	grNnlo.grAll->SetFillStyle(3004);
	CopyStyle(grNnlo.grAll, grNnlo.grAllRatio);

	grNnlo.gr->SetLineColor(kGreen);
	grNnlo.gr->SetFillColor(kGreen);
	grNnlo.gr->SetFillStyle(3004);
	CopyStyle(grNnlo.gr, grNnlo.grRatio);

	grNnlo.h->SetLineColor(kGreen);
	grNnlo.h->SetLineStyle(1);
	CopyStyle(grNnlo.h, grNnlo.hRatio);

	grNNloQ2.h->SetLineColor(kRed); CopyStyle(grNNloQ2.h, grNNloQ2.hRatio);
	grNnloPt.h->SetLineColor(kBlack); CopyStyle(grNnloPt.h, grNnloPt.hRatio);
	grNnloQ2Pt.h->SetLineColor(kBlue);CopyStyle(grNnloQ2Pt.h, grNnloQ2Pt.hRatio);




	//Setting style Data
	grData.gr->SetMarkerStyle(20);
	grData.gr->SetMarkerSize(pointSize);
	CopyStyle(grData.gr, grData.grRatio);

	grData.grAll->SetMarkerStyle(20);
	grData.grAll->SetMarkerSize(pointSize);
	CopyStyle(grData.grAll, grData.grAllRatio);



	//currpad -> cd(1);// a prepnu se do prvniho padu 
	int hash = rand();
	TPad *upPad = new TPad( SF("upPad%d",hash), "upPad", 0, 0.5, 1, 1 );
	upPad->SetBottomMargin(0);
	upPad->SetLeftMargin(leftMargin);
	upPad->SetRightMargin(rightMargin);
	upPad->Draw();
	upPad->cd();

	if(style.isLogY) gPad->SetLogy();
	if(style.isLogX) gPad->SetLogx();

	//gPad -> SetRightMargin(0.02);

	grData.h->SetLineColor(kBlack);
	grData.h->Draw("AXIS");
	grData.h->GetYaxis()->SetTitle(yTitle);
	grData.h->GetYaxis()->SetTitleSize(fontS);
	grData.h->GetYaxis()->SetTitleOffset(0.85);
	grData.h->GetYaxis()->SetLabelSize(fontS);
	if(!style.isLogY) {
		grData.h->SetMinimum(0);
		grData.h->SetMaximum(range[1] * 1.1);
	}
	else {
		grData.h->SetMinimum(range[0] * 0.8);
		grData.h->SetMaximum(range[1] * 1.2);
	}



	grNnlo.grAll->Draw("e2 same");
	grNnlo.gr->Draw("e2 same");
	grNnlo.h->Draw("same");

	grNNloQ2.h->Draw("same");
	grNnloPt.h->Draw("same");
	grNnloQ2Pt.h->Draw("same");

	gStyle->SetEndErrorSize(4);
	grData.gr->Draw("e same");
	grData.grAll->Draw("pz same");



	//Plot Legend
	TLegend *leg = new TLegend( style.legX+0.5-0.5, style.legY+0.60-0.6,
	                            style.legX+0.95-0.5,style.legY+0.88-0.6 );
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	TString tag =  theor_path(theor_path.First('/')+1, theor_path.Last('/')-theor_path.First('/')-1 );
	leg->SetHeader(Title);
	leg->AddEntry(grData.gr,    SF("H1 %s data", tag.Data()), "ep");
	leg->AddEntry(grNnlo.grAll, SF("NNLO  x %.2f", dissFactor) , "fl");
	leg->AddEntry(grNNloQ2.h,   SF("NNLO (#mu^{2}=Q^{2}) x %.2f", dissFactor) , "l");
	leg->AddEntry(grNnloPt.h,   SF("NNLO (#mu^{2}=p_{T}^{2}) x %.2f", dissFactor) , "l");
	leg->AddEntry(grNnloQ2Pt.h, SF("NNLO (#mu^{2}=#frac{Q^{2}}{4}+p_{T}^{2}) x %.2f", dissFactor) , "l");


	leg->Draw();

	gPad->Update();
	gPad->RedrawAxis();
	TFrame *frU = (TFrame *) gPad->FindObject("TFrame");
	frU->SetFillStyle(0);
	frU->Draw();
	gPad->Update();


	currpad->cd();
	TPad *downPad = new TPad( SF("downPad%d",hash), "downPad", 0, 0.0, 1, 0.5 );
	downPad->SetTopMargin(0);
	downPad->SetLeftMargin(leftMargin);
	downPad->SetRightMargin(rightMargin);
	downPad->SetBottomMargin(bottomMargin);
	downPad->Draw();
	downPad->cd();
	//gPad -> SetRightMargin(0.02);
	if(style.isLogX) gPad->SetLogx();

	//double Max = max( hDataRatio->GetMaximum(), hNLORatio->GetMaximum() );
	//grData.hRatio->SetMaximum(Max);

	grData.hRatio->GetYaxis()->SetNdivisions(204);

	grData.hRatio->Draw("AXIS");
	grData.hRatio->SetMinimum( max(0.,range[2] - 0.1*range[3]) );
	grData.hRatio->SetMaximum( 1.1*range[3] );


	grData.hRatio->GetXaxis()->SetTitle(xTitle);
	grData.hRatio->GetXaxis()->SetTitleSize(fontS * 5./5);

	grData.hRatio->GetYaxis()->SetLabelSize(fontS * 5./5);
	grData.hRatio->GetXaxis()->SetLabelSize(fontS * 5./5);

	grData.hRatio->GetYaxis()->SetTitle("Data/NNLO");
	grData.hRatio->GetYaxis()->SetTitleSize(fontS * 5./5);
	grData.hRatio->GetYaxis()->SetTitleOffset(0.84);

	grNnlo.grAllRatio->Draw("e2 same");
	grNnlo.grRatio->Draw("e2 same");
	grNnlo.hRatio->Draw("same");


	grNNloQ2.hRatio->Draw("same");
	grNnloPt.hRatio->Draw("same");
	grNnloQ2Pt.hRatio->Draw("same");


	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");

	gPad->Update();
	gPad->RedrawAxis();
	TFrame *frD = (TFrame *) gPad->FindObject("TFrame");
	frD->SetFillStyle(0);
	frD->Draw();
	gPad->Update();
}
#endif

struct myCanvas {

	TVirtualPad * currpad;

	double leftMargin = 0.15;
	double rightMargin = 0.03;
	double bottomMargin = 0.35;
	double fontS = 0.08;

	double pointSize = 0.8;

	TString xTitle, yTitleUp, yTitleDown;

	TPad *upPad, *downPad;
	vector<double> range;

	void CreateFrame() {
		currpad = gPad;
		//currpad -> cd(1);// a prepnu se do prvniho padu 
		int hash = rand();
		upPad = new TPad( SF("upPad%d",hash), "upPad", 0, 0.5, 1, 1 );
		upPad->SetBottomMargin(0);
		upPad->SetLeftMargin(leftMargin);
		upPad->SetRightMargin(rightMargin);
		upPad->Draw();

		currpad->cd();
		downPad = new TPad( SF("downPad%d",hash), "downPad", 0, 0.0, 1, 0.5 );
		downPad->SetTopMargin(0);
		downPad->SetLeftMargin(leftMargin);
		downPad->SetRightMargin(rightMargin);
		downPad->SetBottomMargin(bottomMargin);
		downPad->Draw();

		currpad->cd();
	}

	static void SetNLOstyle(HistoErr &grNlo) {
		//Setting style NLO
		int green = kGreen - 9;
		green = kGreen - 3;

		grNlo.grAll->SetLineColor(green);
		grNlo.grAll->SetFillColor(green);
		grNlo.grAll->SetFillStyle(3013);//3144
		CopyStyle(grNlo.grAll, grNlo.grAllRatio);

		grNlo.gr->SetLineColor(green);
		grNlo.gr->SetFillColor(green); //or -6
		grNlo.gr->SetFillStyle(3013); //3144
		CopyStyle(grNlo.gr, grNlo.grRatio);

		//grNlo.h->SetLineColor(kWhite);
		grNlo.h->SetLineColor(kTeal+3);
		grNlo.h->SetMarkerColor(kTeal+3);
		grNlo.h->SetLineWidth(2);
		CopyStyle(grNlo.h, grNlo.hRatio);
	}
	static void SetNNLOstyle(HistoErr &grNnlo) {
		//Setting style NNLO
		double allTr = 0.2;
		double inTr = 0.5;
		allTr = 0.35;
		inTr = 0.45;
		grNnlo.grAll->SetLineWidth(0);
		grNnlo.grAll->SetLineColorAlpha(   kBlue-7, allTr /* kAzure-4*/);
		grNnlo.grAll->SetFillColorAlpha(kBlue-7, allTr);
		grNnlo.grAll->SetFillStyle(1001);
		CopyStyle(grNnlo.grAll, grNnlo.grAllRatio);

		grNnlo.gr->SetLineWidth(0);
		grNnlo.gr->SetLineColorAlpha( kBlue-7, inTr);
		grNnlo.gr->SetFillColorAlpha(kBlue-7, inTr);
		grNnlo.gr->SetFillStyle(1001);
		CopyStyle(grNnlo.gr, grNnlo.grRatio);

		grNnlo.h->SetLineColor(kBlue+2);
		grNnlo.h->SetMarkerColor(kBlue+2);
		grNnlo.h->SetLineStyle(1);
		grNnlo.h->SetLineWidth(2);
		CopyStyle(grNnlo.h, grNnlo.hRatio);
	}

	void SetDataStyle(HistoErr &grData, bool isLogX, bool isLogY, double Max = -1.) {
		//Setting style Data
		grData.gr->SetMarkerStyle(20);
		grData.gr->SetMarkerSize(pointSize);
		CopyStyle(grData.gr, grData.grRatio);

		grData.grAll->SetMarkerStyle(20);
		grData.grAll->SetMarkerSize(pointSize);
		CopyStyle(grData.grAll, grData.grAllRatio);
		

		upPad->cd();

		if(isLogY) gPad->SetLogy();
		if(isLogX) gPad->SetLogx();


		grData.h->SetLineColor(kBlack);
		grData.h->Draw("AXIS");
		grData.h->GetYaxis()->SetTitle(yTitleUp);
		grData.h->GetYaxis()->SetTitleSize(fontS);
		grData.h->GetYaxis()->SetTitleOffset(0.85);
		grData.h->GetYaxis()->SetLabelSize(fontS);
		grData.h->GetYaxis()->SetMoreLogLabels();
		if(isLogY) {
			grData.h->SetMinimum(0);
			if(Max == -1.)
				grData.h->SetMaximum(range[1] * 1.1);
			else
				grData.h->SetMaximum(Max);
		}
		else {
			grData.h->SetMinimum(range[0] * 0.8);
			grData.h->SetMaximum(range[1] * 1.2);
		}


		downPad->cd();

		if(isLogX) gPad->SetLogx();
		grData.hRatio->GetYaxis()->SetNdivisions(204);

		grData.hRatio->SetMinimum( max(0.,range[2] - 0.1*range[3]) );
		grData.hRatio->SetMaximum( 1.1*range[3] );


		grData.hRatio->GetXaxis()->SetTitle(xTitle);
		grData.hRatio->GetXaxis()->SetTitleSize(fontS * 5./5);

		grData.hRatio->GetYaxis()->SetLabelSize(fontS * 5./5);
		grData.hRatio->GetXaxis()->SetLabelSize(fontS * 5./5);
		grData.hRatio->GetXaxis()->SetMoreLogLabels();

		grData.hRatio->GetYaxis()->SetTitle(yTitleDown);
		grData.hRatio->GetYaxis()->SetTitleSize(fontS * 5./5);
		grData.hRatio->GetYaxis()->SetTitleOffset(0.84);
		grData.hRatio->Draw("AXIS");

		currpad->cd();
	}
	void MakeNoTickX(HistoErr &grData) {
		upPad->cd();
		grData.h->GetXaxis()->SetTickLength(0.0);
	}

	static void UpdataFrame() {
		gPad->Update();
		gPad->RedrawAxis();
		TFrame *frU = (TFrame *) gPad->FindObject("TFrame");
		frU->SetFillStyle(0);
		frU->Draw();
		gPad->Update();
	}




};

void Histogram::plotDPDFStudies()
{
	myCanvas myCan;


	HistoErr grData, grNlo, grNnlo;
	vector<double> range = LoadHistogramsNloNnloData("nnlo:FitB-0-Q2pPt2-cc", grData, grNlo, grNnlo, true);
	myCan.range = range;
	myCan.yTitleUp =  yTitle;
	myCan.yTitleDown = "Data/NNLO";
	myCan.xTitle = xTitle;

	HistoErr grNnloSJ = LoadHistograms("sjNNLO", defFile,
	                         "nnlo:zeusSJ-0-Q2pPt2-cc", "nnlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloJets = LoadHistograms("NNLOjets", defFile,
	                         "nnlo:FitJets-0-Q2pPt2-cc", "nnlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloFitA = LoadHistograms("NNLOfitA", defFile,
	                         "nnlo:FitA-0-Q2pPt2-cc", "nnlo:FitB-0-Q2pPt2-cc");

	myCan.CreateFrame();

	//Setting NNLO style
	myCanvas::SetNNLOstyle(grNnlo);

	grNnloSJ.h->SetLineColor(kRed); CopyStyle(grNnloSJ.h, grNnloSJ.hRatio);
	grNnloJets.h->SetLineColor(kBlack); CopyStyle(grNnloJets.h, grNnloJets.hRatio);
	grNnloFitA.h->SetLineColor(kBlue); CopyStyle(grNnloFitA.h, grNnloFitA.hRatio);


	//Setting Data style and frame
	myCan.SetDataStyle(grData, style.isLogX, style.isLogY);

	////////////////////////////////////////////
	//Up Frame
	////////////////////////////////////////////

	myCan.upPad->cd();


	if(!style.isLogY) {
		grData.h->SetMinimum(0);
		grData.h->SetMaximum(range[1] * 1.1);
	}
	else {
		grData.h->SetMinimum(range[0] * 0.8);
		grData.h->SetMaximum(range[1] * 1.2);
	}

	//grData.h->Draw("AXIS");



	grNnlo.grAll->Draw("e2 same");
	grNnlo.gr->Draw("e2 same");
	grNnlo.h->Draw("same");

	grNnloSJ.h->Draw("same");
	grNnloJets.h->Draw("same");
	grNnloFitA.h->Draw("same");

	gStyle->SetEndErrorSize(4);
	grData.gr->Draw("e same");
	grData.grAll->Draw("pz same");



	//Plot Legend
	TLegend *leg = new TLegend( style.legX+0.5-0.5, style.legY+0.60-0.6,
	                            style.legX+0.95-0.5,style.legY+0.88-0.6 );
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	TString tag =  theor_path(theor_path.First('/')+1, theor_path.Last('/')-theor_path.First('/')-1 );
	leg->SetHeader(Title);
	leg->AddEntry(grData.gr,    SF("H1 %s data", tag.Data()), "ep");
	leg->AddEntry(grNnlo.grAll, SF("NNLO (Fit B) x %.2f", dissFactor) , "fl");
	leg->AddEntry(grNnloSJ.h,   SF("NNLO (ZEUS SJ) x %.2f", dissFactor) , "l");
	leg->AddEntry(grNnloJets.h, SF("NNLO (Fit Jets) x %.2f", dissFactor) , "l");
	leg->AddEntry(grNnloFitA.h, SF("NNLO (Fit A) x %.2f", dissFactor) , "l");


	leg->Draw();

	myCan.UpdataFrame();

	////////////////////////////////////////////
	//Down Frame
	////////////////////////////////////////////

	myCan.downPad->cd();


	grNnlo.grAllRatio->Draw("e2 same");
	grNnlo.grRatio->Draw("e2 same");
	grNnlo.hRatio->Draw("same");


	grNnloSJ.hRatio->Draw("same");
	grNnloJets.hRatio->Draw("same");
	grNnloFitA.hRatio->Draw("same");


	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");

	myCan.UpdataFrame();
}


void Histogram::plotScaleStudies()
{
	myCanvas myCan;


	HistoErr grData, grNlo, grNnlo;
	vector<double> range = LoadHistogramsNloNnloData("nnlo:FitB-0-Q2pPt2-cc", grData, grNlo, grNnlo, true);
	myCan.range = range;
	myCan.yTitleUp =  yTitle;
	myCan.yTitleDown = "Data/NNLO";
	myCan.xTitle = xTitle;


	HistoErr grNNloQ2  = LoadHistograms("NNLOq2", defFile,
	                         "nnlo:FitB-0-Q2-cc", "nnlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloPt = LoadHistograms("NNLOpt", defFile,
	                         "nnlo:FitB-0-Pt2-cc", "nnlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloQ2Pt = LoadHistograms("NNLOq2pt", defFile,
	                         "nnlo:FitB-0-0.25Q2pPt2-cc", "nnlo:FitB-0-Q2pPt2-cc");


	myCan.CreateFrame();

	//Setting NNLO style
	myCanvas::SetNNLOstyle(grNnlo);


	grNNloQ2.h->SetLineColor(kRed); CopyStyle(grNNloQ2.h, grNNloQ2.hRatio);
	grNnloPt.h->SetLineColor(kBlack); CopyStyle(grNnloPt.h, grNnloPt.hRatio);
	grNnloQ2Pt.h->SetLineColor(kBlue);CopyStyle(grNnloQ2Pt.h, grNnloQ2Pt.hRatio);


	//Setting Data style and frame
	myCan.SetDataStyle(grData, style.isLogX, style.isLogY);

	////////////////////////////////////////////
	//Up Frame
	////////////////////////////////////////////

	myCan.upPad->cd();


	if(!style.isLogY) {
		grData.h->SetMinimum(0);
		grData.h->SetMaximum(range[1] * 1.1);
	}
	else {
		grData.h->SetMinimum(range[0] * 0.8);
		grData.h->SetMaximum(range[1] * 1.2);
	}

	//grData.h->Draw("AXIS");



	grNnlo.grAll->Draw("e2 same");
	grNnlo.gr->Draw("e2 same");
	grNnlo.h->Draw("same");

	grNNloQ2.h->Draw("same");
	grNnloPt.h->Draw("same");
	grNnloQ2Pt.h->Draw("same");

	gStyle->SetEndErrorSize(4);
	grData.gr->Draw("e same");
	grData.grAll->Draw("pz same");



	//Plot Legend

	TLegend *leg = new TLegend( style.legX+0.5-0.5, style.legY+0.60-0.6,
	                            style.legX+0.95-0.5,style.legY+0.88-0.6 );
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	TString tag =  theor_path(theor_path.First('/')+1, theor_path.Last('/')-theor_path.First('/')-1 );
	leg->SetHeader(Title);
	leg->AddEntry(grData.gr,    SF("H1 %s data", tag.Data()), "ep");
	leg->AddEntry(grNnlo.grAll, SF("NNLO  x %.2f", dissFactor) , "fl");
	leg->AddEntry(grNNloQ2.h,   SF("NNLO (#mu^{2}=Q^{2}) x %.2f", dissFactor) , "l");
	leg->AddEntry(grNnloPt.h,   SF("NNLO (#mu^{2}=p_{T}^{2}) x %.2f", dissFactor) , "l");
	leg->AddEntry(grNnloQ2Pt.h, SF("NNLO (#mu^{2}=#frac{Q^{2}}{4}+p_{T}^{2}) x %.2f", dissFactor) , "l");

	leg->Draw();

	myCan.UpdataFrame();

	////////////////////////////////////////////
	//Down Frame
	////////////////////////////////////////////

	myCan.downPad->cd();


	grNnlo.grAllRatio->Draw("e2 same");
	grNnlo.grRatio->Draw("e2 same");
	grNnlo.hRatio->Draw("same");


	grNNloQ2.hRatio->Draw("same");
	grNnloPt.hRatio->Draw("same");
	grNnloQ2Pt.hRatio->Draw("same");


	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");

	myCan.UpdataFrame();
}





void Histogram::plotNLOvsNNLO()
{
	myCanvas myCan;


	HistoErr grData, grNlo, grNnlo;
	vector<double> range = LoadHistogramsNloNnloData("nlo:FitB-0-Q2pPt2-cc", grData, grNlo, grNnlo, false);
	myCan.range = range;
	myCan.yTitleUp =  yTitle;
	myCan.yTitleDown = "Data/NLO";
	myCan.xTitle = xTitle;


	myCan.CreateFrame();

	//Setting NNLO style
	myCanvas::SetNNLOstyle(grNnlo);
	myCanvas::SetNLOstyle(grNlo);


	//Setting Data style and frame
	myCan.SetDataStyle(grData, style.isLogX, style.isLogY);

	////////////////////////////////////////////
	//Up Frame
	////////////////////////////////////////////

	myCan.upPad->cd();


	if(!style.isLogY) {
		grData.h->SetMinimum(0);
		grData.h->SetMaximum(range[1] * 1.1);
	}
	else {
		grData.h->SetMinimum(range[0] * 0.8);
		grData.h->SetMaximum(range[1] * 1.2);
	}


	grNlo.grAll->Draw("e2 same");
	grNlo.gr->Draw("e2 same");
	grNlo.h->Draw("same");


	grNnlo.grAll->Draw("e2 same");
	grNnlo.gr->Draw("e2 same");
	grNnlo.h->Draw("same");


	gStyle->SetEndErrorSize(4);
	grData.gr->Draw("e same");
	grData.grAll->Draw("pz same");



	//Plot Legend

	TLegend *leg = new TLegend( style.legX+0.5-0.5, style.legY+0.60-0.6,
	                            style.legX+0.95-0.5,style.legY+0.88-0.6 );
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	TString tag =  theor_path(theor_path.First('/')+1, theor_path.Last('/')-theor_path.First('/')-1 );
	leg->SetHeader(Title);
	leg->AddEntry(grData.gr,    SF("H1 %s data", tag.Data()), "ep");

	leg->AddEntry(grNlo.grAll,     SF("NLO x %.2f",dissFactor) , "fl");
	leg->AddEntry(grNnlo.grAll,   SF("NNLO  x %.2f", dissFactor) , "fl");

	leg->Draw();

	myCan.UpdataFrame();

	////////////////////////////////////////////
	//Down Frame
	////////////////////////////////////////////

	myCan.downPad->cd();


	grNlo.grAllRatio->Draw("e2 same");
	grNlo.grRatio->Draw("e2 same");
	grNlo.hRatio->Draw("same");


	grNnlo.grAllRatio->Draw("e2 same");
	grNnlo.grRatio->Draw("e2 same");
	grNnlo.hRatio->Draw("same");


	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");

	myCan.UpdataFrame();
}


void Histogram::plotNLOvsNNLOratio(TString plotStyle, double minY, double maxY)
{
	//myCanvas myCan;


	HistoErr grData, grNlo, grNnlo;
	vector<double> range = LoadHistogramsNloNnloData("nlo:FitB-0-Q2pPt2-cc", grData, grNlo, grNnlo, plotStyle == "DPDF");
	//myCan.range = range;
	//myCan.yTitleUp =  yTitle;
	//myCan.yTitleDown = "Data/NLO";
	//myCan.xTitle = xTitle;


	HistoErr grNnloSJ = LoadHistograms("sjNNLO", defFile,
	                         "nnlo:zeusSJ-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloJets = LoadHistograms("NNLOjets", defFile,
	                         "nnlo:FitJets-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloFitA = LoadHistograms("NNLOfitA", defFile,
	                         "nnlo:FitA-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloMRW = LoadHistograms("NNLOMRW", defFile,
	                         "nnlo:MRW-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");

	HistoErr grNnloQ2, grNnloPt, grNnloQ2Pt;
	if(plotStyle == "Scale") {
		grNnloQ2  = LoadHistograms("NNLOq2", defFile,
											 "nnlo:FitB-0-Q2-cc", "nlo:FitB-0-Q2pPt2-cc");
		grNnloPt = LoadHistograms("NNLOpt", defFile,
											 "nnlo:FitB-0-Pt2-cc", "nlo:FitB-0-Q2pPt2-cc");
		grNnloQ2Pt = LoadHistograms("NNLOq2pt", defFile,
											 "nnlo:FitB-0-0.25Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	}


	//myCan.CreateFrame();

	//Setting NNLO style
	myCanvas::SetNNLOstyle(grNnlo);
	myCanvas::SetNLOstyle(grNlo);


	//Setting Data style and frame
	//myCan.SetDataStyle(grData, style.isLogX, style.isLogY);

	double pointSize = 0.8;
	gStyle->SetEndErrorSize(4);
	grData.grRatio->SetMarkerStyle(20);
	grData.grRatio->SetMarkerSize(pointSize);

	grData.grAllRatio->SetMarkerStyle(20);
	grData.grAllRatio->SetMarkerSize(pointSize);




	////////////////////////////////////////////
	//Down Frame
	////////////////////////////////////////////
	double tickSize = 0.02;

	gPad->SetLeftMargin(0);
	gPad->SetRightMargin(0);

	gPad->SetBottomMargin(0.25);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks(1,1);

	//cout << "RADEK paht " << theor_path << endl;
	//grData.hRatio->SetTitle("test");
	//grNlo.grAllRatio->SetTitle("test");
	if(style.isLogX) gPad->SetLogx();
	if(style.isLogX)
		grData.hRatio->GetXaxis()->SetLabelOffset(-0.060);
	else
		grData.hRatio->GetXaxis()->SetLabelOffset(-0.020);

	grData.hRatio->GetXaxis()->SetTitleOffset(0.7);
	grData.hRatio->GetXaxis()->SetTitle(xTitle);
	grData.hRatio->GetXaxis()->SetTitleSize(0.10);
	grData.hRatio->GetXaxis()->SetLabelSize(0.10);
	grData.hRatio->GetXaxis()->CenterTitle();

	grData.hRatio->GetYaxis()->SetNdivisions(505);
	grData.hRatio->GetYaxis()->SetTickLength(3.8*tickSize);
	grData.hRatio->GetXaxis()->SetTickLength(tickSize);

	grData.hRatio->GetXaxis()->SetNdivisions(505);

	if(style.isLogX)
		grData.hRatio->GetXaxis()->SetMoreLogLabels();



	if(minY < 0 || maxY < 0) {
		grData.hRatio->SetMinimum( max(0.,range[2] - 0.1*range[3]) );
		grData.hRatio->SetMaximum( 1.1*range[3] );
	}
	else {
		grData.hRatio->SetMinimum(minY);
		grData.hRatio->SetMaximum(maxY);
	}

	grData.hRatio->Draw("AXIS");

	//RemoveOverlaps(gPad, grData.hRatio->GetXaxis(), true, true);
	//grNlo.grAllRatio->Draw("e2 same");

	gStyle->SetHatchesSpacing(4);
	gStyle->SetHatchesLineWidth(1.8);

	if(plotStyle == "NLOvsNNLO") grNlo.grRatio->Draw("e2 same");
	grNlo.hRatio->Draw("same");


	grNnlo.grAllRatio->Draw("e2 same");
	grNnlo.grRatio->Draw("e2 same");
	grNnlo.hRatio->Draw("same");



	if(plotStyle == "DPDF") {
		grNnloSJ.hRatio->SetLineColor(kRed); grNnloSJ.hRatio->SetLineWidth(2);
		grNnloJets.hRatio->SetLineColor(kMagenta); grNnloJets.hRatio->SetLineWidth(2);
		grNnloFitA.hRatio->SetLineColor(kBlack); grNnloFitA.hRatio->SetLineWidth(2);
		grNnloMRW.hRatio->SetLineColor(kOrange); grNnloMRW.hRatio->SetLineWidth(2);


		grNnloSJ.hRatio->Draw("same");
		grNnloJets.hRatio->Draw("same");
		grNnloFitA.hRatio->Draw("same");
		grNnloMRW.hRatio->Draw("same");
	}
	else if(plotStyle == "Scale") {
		grNnloQ2.hRatio->SetLineColor(kRed); grNnloQ2.hRatio->SetLineWidth(2);
		grNnloPt.hRatio->SetLineColor(kBlack); grNnloPt.hRatio->SetLineWidth(2);
		grNnloQ2Pt.hRatio->SetLineColor(kCyan); grNnloQ2Pt.hRatio->SetLineWidth(2);

		grNnloQ2.hRatio->Draw("same");
		grNnloPt.hRatio->Draw("same");
		grNnloQ2Pt.hRatio->Draw("same");
	}

	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");








	double chiNLO  = getChi2(grData.h, grData.gr, grNlo.h);
	double chiNNLO = getChi2(grData.h, grData.gr, grNnlo.h);

	//cout <<"HOPE "<< chiNLO << " "<< chiNNLO << endl;
	/*
	TLegend *leg = new TLegend(0.2, 0.27, 0.8, 0.32);
	leg->SetBorderSize(0);
	if(chiNLO < chiNNLO)
		leg->SetFillColorAlpha( kGreen - 9, 0.5);
	else
		leg->SetFillColorAlpha(kBlue - 7, 0.5);

	leg->SetTextSize(0.08);
	leg->SetHeader( SF("#chi^{2}_{nlo}=%2.1f, #chi^{2}_{nnlo}=%2.1f", chiNLO, chiNNLO) );

	leg->Draw();
	*/
	//leg->AddEntry(grData.gr,    SF("H1 %s data", tag.Data()), "ep");
	//leg->AddEntry(grNlo.grAll,     SF("NLO x %.2f",dissFactor) , "fl");
	//leg->AddEntry(grNnlo.grAll,   SF("NNLO  x %.2f", dissFactor) , "fl");

	//leg->Draw();



	HistoErr grNloSJ = LoadHistograms("sjNLO", defFile,
	                         "nlo:zeusSJ-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNloJets = LoadHistograms("NLOjets", defFile,
	                         "nlo:FitJets-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNloFitA = LoadHistograms("NLOfitA", defFile,
	                         "nlo:FitA-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");

	HistoErr grNloFitBup = LoadHistograms("NLOfitBup", defFile,
	                         "nlo:FitB-0-Q2pPt2-uu", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloFitBup = LoadHistograms("NNLOfitBup", defFile,
	                         "nnlo:FitB-0-Q2pPt2-uu", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNloFitBdn = LoadHistograms("NLOfitBdn", defFile,
	                         "nlo:FitB-0-Q2pPt2-dd", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloFitBdn = LoadHistograms("NNLOfitBdn", defFile,
	                         "nnlo:FitB-0-Q2pPt2-dd", "nlo:FitB-0-Q2pPt2-cc");

	HistoErr grNloJetsup = LoadHistograms("NLOJetsup", defFile,
	                         "nlo:FitJets-0-Q2pPt2-uu", "nlo:FitJets-0-Q2pPt2-cc");
	HistoErr grNnloJetsup = LoadHistograms("NNLOJetsup", defFile,
	                         "nnlo:FitJets-0-Q2pPt2-uu", "nlo:FitJets-0-Q2pPt2-cc");
	HistoErr grNloJetsdn = LoadHistograms("NLOJetsdn", defFile,
	                         "nlo:FitJets-0-Q2pPt2-dd", "nlo:FitJets-0-Q2pPt2-cc");
	HistoErr grNnloJetsdn = LoadHistograms("NNLOJetsdn", defFile,
	                         "nnlo:FitJets-0-Q2pPt2-dd", "nlo:FitJets-0-Q2pPt2-cc");








	if(chiNLO < chiNNLO)
		++ScoreFitB[0];
	else
		++ScoreFitB[1];

	chiNLO  = getChi2(grData.h, grData.gr, grNloSJ.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloSJ.h);

	if(chiNLO < chiNNLO)
		++ScoreSJ[0];
	else
		++ScoreSJ[1];

	chiNLO  = getChi2(grData.h, grData.gr, grNloJets.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloJets.h);

	if(chiNLO < chiNNLO)
		++ScoreJets[0];
	else
		++ScoreJets[1];

	chiNLO  = getChi2(grData.h, grData.gr, grNloFitA.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloFitA.h);

	if(chiNLO < chiNNLO)
		++ScoreFitA[0];
	else
		++ScoreFitA[1];


	//Fit Bup
	chiNLO  = getChi2(grData.h, grData.gr, grNloFitBup.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloFitBup.h);

	if(chiNLO < chiNNLO)
		++ScoreFitBup[0];
	else
		++ScoreFitBup[1];

	//Fit Bdn
	chiNLO  = getChi2(grData.h, grData.gr, grNloFitBdn.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloFitBdn.h);

	if(chiNLO < chiNNLO)
		++ScoreFitBdn[0];
	else
		++ScoreFitBdn[1];



	//Fit Jetsup
	chiNLO  = getChi2(grData.h, grData.gr, grNloJetsup.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloJetsup.h);

	if(chiNLO < chiNNLO)
		++ScoreJetsup[0];
	else
		++ScoreJetsup[1];

	//Fit Jetsdn
	chiNLO  = getChi2(grData.h, grData.gr, grNloJetsdn.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloJetsdn.h);

	if(chiNLO < chiNNLO)
		++ScoreJetsdn[0];
	else
		++ScoreJetsdn[1];










	myCanvas::UpdataFrame();
}

double getBestChi2(TString anal, TString var, vector<double> thVec);

void Histogram::plotNLOvsNNLOratioAbs(TString plotStyle, double minY, double maxY, double minYabs, double maxYabs, double Factor)
{
	//myCanvas myCan;


	HistoErr grData, grNlo, grNnlo;
	vector<double> range = LoadHistogramsNloNnloData("nlo:FitB-0-Q2pPt2-cc", grData, grNlo, grNnlo, plotStyle == "DPDF");
	//myCan.range = range;
	//myCan.yTitleUp =  yTitle;
	//myCan.yTitleDown = "Data/NLO";
	//myCan.xTitle = xTitle;


	HistoErr grNnloSJ = LoadHistograms("sjNNLO", defFile,
	                         "nnlo:zeusSJ-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloJets = LoadHistograms("NNLOjets", defFile,
	                         "nnlo:FitJets-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloFitA = LoadHistograms("NNLOfitA", defFile,
	                         "nnlo:FitA-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloMRW = LoadHistograms("NNLOMRW", defFile,
	                         "nnlo:MRW-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloFit19 = LoadHistograms("NNLOFit19", defFile,
	                         "nnlo:Fit19-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");

	HistoErr grNnloQ2, grNnloPt, grNnloQ2Pt;
	if(plotStyle == "Scale") {
		grNnloQ2  = LoadHistograms("NNLOq2", defFile,
											 "nnlo:FitB-0-Q2-cc", "nlo:FitB-0-Q2pPt2-cc");
		grNnloPt = LoadHistograms("NNLOpt", defFile,
											 "nnlo:FitB-0-Pt2-cc", "nlo:FitB-0-Q2pPt2-cc");
		grNnloQ2Pt = LoadHistograms("NNLOq2pt", defFile,
											 "nnlo:FitB-0-0.25Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	}


	//myCan.CreateFrame();

	//Setting NNLO style
	myCanvas::SetNNLOstyle(grNnlo);
	myCanvas::SetNLOstyle(grNlo);


	//Setting Data style and frame
	//myCan.SetDataStyle(grData, style.isLogX, style.isLogY);

	double pointSize = 0.8;
	gStyle->SetEndErrorSize(4);
    auto decorateData = [&](TGraphAsymmErrors *gr){gr->SetMarkerStyle(20);gr->SetMarkerSize(pointSize);};
    decorateData(grData.grRatio);
    decorateData(grData.grAllRatio);
    decorateData(grData.gr);
    decorateData(grData.grAll);

	//grData.grRatio->SetMarkerStyle(20);
	//grData.grRatio->SetMarkerSize(pointSize);
	//grData.grAllRatio->SetMarkerStyle(20);
	//grData.grAllRatio->SetMarkerSize(pointSize);




	////////////////////////////////////////////
	//Down Frame
	////////////////////////////////////////////

	gPad->SetTicks(1,1);

	if(style.isLogX) gPad->SetLogx();

    //grData.hRatio->SetName(SF("%s_%d", grData.hRatio->GetName(), rand()));
	grData.hRatio->Draw("AXIS");


    /*
    TH1D *h1 = new TH1D(SF("h1%d",rand()), "h1", 20, -5, 5);
    for(int i = 1; i <= h1->GetNbinsX(); ++i) {
        h1->SetBinContent(i, 1);
        h1->SetBinError(i, 0.4);
    }
    TH1D *h2 = (TH1D*) h1->Clone(SF("%d",rand()));

	gStyle->SetHatchesSpacing(8);
    h1->SetFillStyle(3144);
    h1->SetFillColor(kOrange);

    h1->Draw("e2" );
    SetFTO({12}, {5}, {1.4, 2.6, 0.3, 3.7});
    gPad->Update();
    */

    
    auto Decorator = [&](TString xTitle) {
        SetFTO({12}, {5}, {1.4 + (var_name == "beta")*0.2, 2.6, 0.3, 3.7});
        GetXaxis()->SetNdivisions(505);
        GetYaxis()->SetNdivisions(505);
        GetXaxis()->SetTitle(xTitle);
        GetXaxis()->CenterTitle();
        GetYaxis()->CenterTitle();
    };

    Decorator(xTitle);
    GetYaxis()->SetTitle("#sigma/#sigma_{NLO}");

	if(style.isLogX && var_name != "beta") GetXaxis()->SetMoreLogLabels();

	if(minY < 0 || maxY < 0) {
		GetFrame()->SetMinimum( max(0.,range[2] - 0.1*range[3]) );
		GetFrame()->SetMaximum( 1.1*range[3] );
	}
	else {
		GetFrame()->SetMinimum(minY);
		GetFrame()->SetMaximum(maxY);
	}

	RemoveOverlaps(gPad, grData.hRatio->GetXaxis(), true, true);
	//grNlo.grAllRatio->Draw("e2 same");


	gStyle->SetHatchesSpacing(4);
	gStyle->SetHatchesLineWidth(1.8);

	if(plotStyle == "NLOvsNNLO") grNlo.grRatio->Draw("e2 same");
    

	grNlo.hRatio->Draw("same ][");


	grNnlo.grAllRatio->Draw("e2 same");
	grNnlo.grRatio->Draw("e2 same");
	grNnlo.hRatio->Draw("same ][");




	if(plotStyle == "DPDF") {
		grNnloSJ.hRatio->SetLineColor(kRed); grNnloSJ.hRatio->SetLineWidth(2);
		grNnloJets.hRatio->SetLineColor(kMagenta); grNnloJets.hRatio->SetLineWidth(2);
		grNnloFitA.hRatio->SetLineColor(kBlack); grNnloFitA.hRatio->SetLineWidth(2);
		grNnloMRW.hRatio->SetLineColor(kOrange); grNnloMRW.hRatio->SetLineWidth(2);


		grNnloSJ.hRatio->Draw("same");
		grNnloJets.hRatio->Draw("same");
		grNnloFitA.hRatio->Draw("same");
		grNnloMRW.hRatio->Draw("same");
	}
	else if(plotStyle == "Scale") {
		grNnloQ2.hRatio->SetLineColor(kRed); grNnloQ2.hRatio->SetLineWidth(2);
		grNnloPt.hRatio->SetLineColor(kBlack); grNnloPt.hRatio->SetLineWidth(2);
		grNnloQ2Pt.hRatio->SetLineColor(kCyan); grNnloQ2Pt.hRatio->SetLineWidth(2);

		grNnloQ2.hRatio->Draw("same ][");
		grNnloPt.hRatio->Draw("same ][");
		grNnloQ2Pt.hRatio->Draw("same ][");
	}


	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");


	myCanvas::UpdataFrame();


	////////////////////////////////////////////
	//Up Frame
	////////////////////////////////////////////

    TCanvas *cAll = gPad->GetCanvas();
    cout << "Katerina " << gPad->GetNumber() << endl;
    cAll->cd(gPad->GetNumber() - GetNpads(cAll)/2 );
	gPad->SetTicks(1,1);

    //h2->DrawCopy("e2" );
    //return;


	if(style.isLogY) gPad->SetLogy();
	if(style.isLogX) gPad->SetLogx();

    auto scaleGraph = [](TGraphAsymmErrors *gr, double scl) {
        for(int i = 0; i < gr->GetN(); ++i) {
            double x, y;
            gr->GetPoint(i, x, y);
            double eL = gr->GetErrorYlow(i);
            double eH = gr->GetErrorYhigh(i);

            eL *= scl;
            eH *= scl;
            y  *= scl;
            
            gr->SetPoint(i, x, y);
            gr->SetPointEYlow(i,  eL);
            gr->SetPointEYhigh(i,  eH);
        }
    };

    grData.h->Scale(Factor);
    grNlo.h->Scale(Factor);
    grNnlo.h->Scale(Factor);

    scaleGraph(grData.gr, Factor);
    scaleGraph(grData.grAll, Factor);
    scaleGraph(grNlo.gr, Factor);
    scaleGraph(grNnlo.gr, Factor);
    scaleGraph(grNnlo.grAll, Factor);


	grData.h->Draw("AXIS");

    //return;

	if(plotStyle == "NLOvsNNLO") grNlo.gr->Draw("e2 same");
	grNlo.h->Draw("same ][");


	grNnlo.grAll->Draw("e2 same");
	grNnlo.gr->Draw("e2 same");
	grNnlo.h->Draw("same ][");


	if(plotStyle == "DPDF") {
        grNnloSJ.h->Scale(Factor);
        grNnloJets.h->Scale(Factor);
        grNnloFitA.h->Scale(Factor);
        grNnloMRW.h->Scale(Factor);

		grNnloSJ.h->SetLineColor(kRed); grNnloSJ.h->SetLineWidth(2);
		grNnloJets.h->SetLineColor(kMagenta); grNnloJets.h->SetLineWidth(2);
		grNnloFitA.h->SetLineColor(kBlack); grNnloFitA.h->SetLineWidth(2);
		grNnloMRW.h->SetLineColor(kOrange); grNnloMRW.h->SetLineWidth(2);


		grNnloSJ.h->Draw("same");
		grNnloJets.h->Draw("same");
		grNnloFitA.h->Draw("same");
		grNnloMRW.h->Draw("same");
	}
	else if(plotStyle == "Scale") {

		grNnloQ2.h->Scale(Factor);
		grNnloPt.h->Scale(Factor);
		grNnloQ2Pt.h->Scale(Factor);


		grNnloQ2.h->SetLineColor(kRed); grNnloQ2.h->SetLineWidth(2);
		grNnloPt.h->SetLineColor(kBlack); grNnloPt.h->SetLineWidth(2);
		grNnloQ2Pt.h->SetLineColor(kCyan); grNnloQ2Pt.h->SetLineWidth(2);

		grNnloQ2.h->Draw("same ][");
		grNnloPt.h->Draw("same ][");
		grNnloQ2Pt.h->Draw("same ][");
	}







	grData.gr->Draw("p same");
	grData.grAll->Draw("pz same");



    Decorator("");
    GetFrame()->SetMinimum(minYabs);
    GetFrame()->SetMaximum(maxYabs);
    GetYaxis()->SetTitle(yTitle);




	//cout <<"HOPE "<< chiNLO << " "<< chiNNLO << endl;
	/*
	TLegend *leg = new TLegend(0.2, 0.27, 0.8, 0.32);
	leg->SetBorderSize(0);
	if(chiNLO < chiNNLO)
		leg->SetFillColorAlpha( kGreen - 9, 0.5);
	else
		leg->SetFillColorAlpha(kBlue - 7, 0.5);

	leg->SetTextSize(0.08);
	leg->SetHeader( SF("#chi^{2}_{nlo}=%2.1f, #chi^{2}_{nnlo}=%2.1f", chiNLO, chiNNLO) );

	leg->Draw();
	*/
	//leg->AddEntry(grData.gr,    SF("H1 %s data", tag.Data()), "ep");
	//leg->AddEntry(grNlo.grAll,     SF("NLO x %.2f",dissFactor) , "fl");
	//leg->AddEntry(grNnlo.grAll,   SF("NNLO  x %.2f", dissFactor) , "fl");

	//leg->Draw();


    /*

	HistoErr grNloSJ = LoadHistograms("sjNLO", defFile,
	                         "nlo:zeusSJ-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNloJets = LoadHistograms("NLOjets", defFile,
	                         "nlo:FitJets-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNloFitA = LoadHistograms("NLOfitA", defFile,
	                         "nlo:FitA-0-Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");

	HistoErr grNloFitBup = LoadHistograms("NLOfitBup", defFile,
	                         "nlo:FitB-0-Q2pPt2-uu", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloFitBup = LoadHistograms("NNLOfitBup", defFile,
	                         "nnlo:FitB-0-Q2pPt2-uu", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNloFitBdn = LoadHistograms("NLOfitBdn", defFile,
	                         "nlo:FitB-0-Q2pPt2-dd", "nlo:FitB-0-Q2pPt2-cc");
	HistoErr grNnloFitBdn = LoadHistograms("NNLOfitBdn", defFile,
	                         "nnlo:FitB-0-Q2pPt2-dd", "nlo:FitB-0-Q2pPt2-cc");

	HistoErr grNloJetsup = LoadHistograms("NLOJetsup", defFile,
	                         "nlo:FitJets-0-Q2pPt2-uu", "nlo:FitJets-0-Q2pPt2-cc");
	HistoErr grNnloJetsup = LoadHistograms("NNLOJetsup", defFile,
	                         "nnlo:FitJets-0-Q2pPt2-uu", "nlo:FitJets-0-Q2pPt2-cc");
	HistoErr grNloJetsdn = LoadHistograms("NLOJetsdn", defFile,
	                         "nlo:FitJets-0-Q2pPt2-dd", "nlo:FitJets-0-Q2pPt2-cc");
	HistoErr grNnloJetsdn = LoadHistograms("NNLOJetsdn", defFile,
	                         "nnlo:FitJets-0-Q2pPt2-dd", "nlo:FitJets-0-Q2pPt2-cc");

    */


    #if 1
    vector<TString> fitNames = {"FitB", "FitA", "FitJets", "zeusSJ", "MRW"};
    cout << "Chi2 " << getTag() <<" "<< var_name <<" "<< xTitle<<" ";
    for(const TString &Sn : fitNames) {
        const char *n = Sn.Data();

        HistoErr GrNlo = LoadHistograms(SF("NLO%scntMy",n), defFile,
                                 SF("nlo:%s-0-Q2pPt2-cc",n), SF("nlo:%s-0-Q2pPt2-cc",n));
        HistoErr GrNnlo= LoadHistograms(SF("NNLO%scntMy",n), defFile,
                                 SF("nnlo:%s-0-Q2pPt2-cc",n), SF("nlo:%s-0-Q2pPt2-cc",n));


        HistoErr GrNloUp = LoadHistograms(SF("NLO%supMy",n), defFile,
                                 SF("nlo:%s-0-Q2pPt2-uu",n), SF("nlo:%s-0-Q2pPt2-cc",n));
        HistoErr GrNnloUp= LoadHistograms(SF("NNLO%supMy",n), defFile,
                                 SF("nnlo:%s-0-Q2pPt2-uu",n), SF("nlo:%s-0-Q2pPt2-cc",n));

        HistoErr GrNloDn = LoadHistograms(SF("NLO%sdnMy",n), defFile,
                                 SF("nlo:%s-0-Q2pPt2-dd",n), SF("nlo:%s-0-Q2pPt2-cc",n));
        HistoErr GrNnloDn= LoadHistograms(SF("NNLO%sdnMy",n), defFile,
                                 SF("nnlo:%s-0-Q2pPt2-dd",n), SF("nlo:%s-0-Q2pPt2-cc",n));

        //FitCentral
        //double chiNLO  = getChi2(grData.h, grData.grAll, GrNlo.h);
        //double chiNNLO = getChi2(grData.h, grData.grAll, GrNnlo.h);
        //cout <<  chiNLO << " "<< chiNNLO << "   ";

        //My chi2
        double chiNLOmy  = getBestChi2(getTag(), var_name, Hist2Vec(GrNlo.h));
        double chiNNLOmy = getBestChi2(getTag(), var_name, Hist2Vec(GrNnlo.h));
        cout <<  chiNLOmy << " "<< chiNNLOmy << "   ";
        //cout << endl<<"RADEKold  " << chiNLO << " "<< chiNNLO << "   " << endl;;
        //cout << "RADEKchi2 " << chiNLOmy << " " <<  chiNNLOmy << endl;

        //Fit Bup
        //chiNLO  = getChi2(grData.h, grData.grAll, GrNloUp.h);
        //chiNNLO = getChi2(grData.h, grData.grAll, GrNnloUp.h);
        //cout <<  chiNLO << " "<< chiNNLO << "   ";

        chiNLOmy  = getBestChi2(getTag(), var_name, Hist2Vec(GrNloUp.h));
        chiNNLOmy = getBestChi2(getTag(), var_name, Hist2Vec(GrNnloUp.h));
        cout <<  chiNLOmy << " "<< chiNNLOmy << "   ";


        //Fit Bdn
        //chiNLO  = getChi2(grData.h, grData.grAll, GrNloDn.h);
        //chiNNLO = getChi2(grData.h, grData.grAll, GrNnloDn.h);
        //cout <<  chiNLO << " "<< chiNNLO <<"     ";

        chiNLOmy  = getBestChi2(getTag(), var_name, Hist2Vec(GrNloDn.h));
        chiNNLOmy = getBestChi2(getTag(), var_name, Hist2Vec(GrNnloDn.h));
        cout <<  chiNLOmy << " "<< chiNNLOmy << "   ";

    }
    cout << endl;
    #endif







    /*
    //FitB
	double chiNLO  = getChi2(grData.h, grData.gr, grNlo.h);
	double chiNNLO = getChi2(grData.h, grData.gr, grNnlo.h);
    cout << "Chi2 " << getTag() <<" "<< var_name <<" "<< xTitle<<" "<<  chiNLO << " "<< chiNNLO << "   ";


	//Fit Bup
	chiNLO  = getChi2(grData.h, grData.gr, grNloFitBup.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloFitBup.h);
    cout <<  chiNLO << " "<< chiNNLO << "   ";

	//Fit Bdn
	chiNLO  = getChi2(grData.h, grData.gr, grNloFitBdn.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloFitBdn.h);
    cout <<  chiNLO << " "<< chiNNLO <<"   ";

	//FitA
	chiNLO  = getChi2(grData.h, grData.gr, grNloFitA.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloFitA.h);
    cout <<  chiNLO << " "<< chiNNLO <<"   ";


    //FitJets
	chiNLO  = getChi2(grData.h, grData.gr, grNloJets.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloJets.h);
    cout <<  chiNLO << " "<< chiNNLO <<"   ";

    //FitSJ
	chiNLO  = getChi2(grData.h, grData.gr, grNloSJ.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloSJ.h);
    cout <<  chiNLO << " "<< chiNNLO << endl;
    */

    /*
	//Fit Jetsup
	chiNLO  = getChi2(grData.h, grData.gr, grNloJetsup.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloJetsup.h);


	//Fit Jetsdn
	chiNLO  = getChi2(grData.h, grData.gr, grNloJetsdn.h);
	chiNNLO = getChi2(grData.h, grData.gr, grNnloJetsdn.h);
    */


	myCanvas::UpdataFrame();
}




void PlotTotal(TCanvas *can, Histogram *h1, Histogram *h2, Histogram *h3, Histogram *h4, Histogram *h5, Histogram *h6)
{

	vector<Histogram> histos = {*h1,*h2,*h3,*h4,*h5,*h6};


	HistoErr grData, grNlo, grNnlo;
	HistoErr grNnloSJ, grNnloJets, grNnloFitA;

	vector<HistoErr> grDataVec(6), grNloVec(6), grNnloVec(6);
	vector<HistoErr> grNnloSJVec(6), grNnloJetsVec(6), grNnloFitAVec(6);

	vector<double> range = { 1e30, -1e30, 1e30, -1e30 };
	for(unsigned i = 0; i < histos.size(); ++i) {
		vector<double> rangeTemp = histos[i].LoadHistogramsNloNnloData("nlo:FitB-0-Q2pPt2-cc", grDataVec[i], grNloVec[i], grNnloVec[i], false);

		range[0] = min(range[0], rangeTemp[0]);
		range[2] = min(range[2], rangeTemp[2]);
		range[1] = max(range[1], rangeTemp[1]);
		range[3] = max(range[3], rangeTemp[3]);
	}
	grData = Merge("Data", grDataVec);
	grNlo  = Merge("NLO", grNloVec);
	grNnlo = Merge("NNLO", grNnloVec);

    //cout << "ZEUSerr " << grNlo.grAll->GetErrorYhigh(6) << " "<< grNlo.grAll->GetErrorYlow(6) << endl;

	TString Names[] = { "FPS",   "VFPS", "LRG", "LRG_H1", "LRGH1820", "LRGZEUS" } ;

	//map<TString, TString> analMap;
	//analMap["FPS"] = "#splitline{H1 FPS}{(HERA #Iota#Iota)}";
	//analMap["VFPS"] = "#splitline{H1 VFPS}{(HERA #Iota#Iota)}";
	//analMap["LRG"] = "#splitline{H1 LRG}{(HERA #Iota#Iota)}}";
	//analMap["LRG_H1"] = "#splitline{H1 LRG}{(HERA #Iota)}}";
	//analMap["LRGH1820"] = "#splitline{H1 LRG}{    (300 GeV)}";
	//analMap["LRGZEUS"] = "#splitline{ZEUS LRG}{(HERA #Iota)}";




	for(int i = 0; i < 6; ++i) {
		grData.hRatio->GetXaxis()->SetBinLabel(i+1, analMap.at(Names[i]));
	}
	grData.hRatio->LabelsOption("h");
	grData.hRatio->LabelsDeflate();


	myCanvas myCan;


	myCan.range = range;
	myCan.yTitleUp =  "#sigma_{tot} [pb]";
	myCan.yTitleDown = "#sigma/#sigma_{NLO}";
	myCan.xTitle = "";
	myCan.fontS  = 0.065;



	//My start
	myCan.CreateFrame();

	//Setting NNLO style
	myCanvas::SetNNLOstyle(grNnlo);

	//Setting NLO style
	myCanvas::SetNLOstyle(grNlo);




	//Setting Data style and frame
	myCan.SetDataStyle(grData, false, true);
	grData.h->GetYaxis()->SetNoExponent();

	//myCan.MakeNoTickX(grData);
	grData.h->GetXaxis()->SetTickLength(0.0);
	grData.hRatio->GetXaxis()->SetLabelOffset(0.03);
	grData.hRatio->GetXaxis()->SetLabelSize(myCan.fontS * 1.4722 * 0.9);
	//cout << "RADEK " << grData.hRatio->GetLabelFont("X") << " "<< grData.hRatio->GetLabelFont("Y") << endl;


	double pointSize = 1.2;
	gStyle->SetEndErrorSize(8.);
	grData.grAllRatio->SetMarkerStyle(20);
	grData.grAllRatio->SetMarkerSize(pointSize);
	grData.grAll->SetMarkerStyle(20);
	grData.grAll->SetMarkerSize(pointSize);



	//exit(0);

	////////////////////////////////////////////
	//Up Frame
	////////////////////////////////////////////

	myCan.upPad->cd();


	if(!true) {
		grData.h->SetMinimum(0);
		grData.h->SetMaximum(range[1] * 1.1);
	}
	else {
		grData.h->SetMinimum(range[0] * 0.8);
		grData.h->SetMaximum(range[1] * 1.2);
	}

	//grData.h->Draw("AXIS");

	grNlo.h->SetBinError(1, 1e-6);
	grNnlo.h->SetBinError(1, 1e-6);


	grNlo.grAll->Draw("e2 same");
	grNlo.gr->Draw("e2 same");
	grNlo.h->Draw("same");


	grNnlo.grAll->Draw("e2 same");
	grNnlo.gr->Draw("e2 same");
	grNnlo.h->Draw("same");


    //Print out of NNLO and NLO
    ofstream tableTot("tableTot.txt");
    for(int i = 0; i < grNlo.h->GetNbinsX(); ++i) {
        //Inside is Scale Err
        double NLOsclDn = grNlo.gr->GetErrorYlow(i);
        double NLOsclUp = grNlo.gr->GetErrorYhigh(i);

        double NNLOsclDn = grNnlo.gr->GetErrorYlow(i);
        double NNLOsclUp = grNnlo.gr->GetErrorYhigh(i);

        double NLOallDn = grNlo.grAll->GetErrorYlow(i);
        double NLOallUp = grNlo.grAll->GetErrorYhigh(i);

        double NNLOallDn = grNnlo.grAll->GetErrorYlow(i);
        double NNLOallUp = grNnlo.grAll->GetErrorYhigh(i);

        double NLOpdfDn = sqrt(pow(NLOallDn,2) - pow(NLOsclDn,2));
        double NLOpdfUp = sqrt(pow(NLOallUp,2) - pow(NLOsclUp,2));

        double NNLOpdfDn = sqrt(pow(NNLOallDn,2) - pow(NNLOsclDn,2));
        double NNLOpdfUp = sqrt(pow(NNLOallUp,2) - pow(NNLOsclUp,2));

        double nlo  = grNlo.h->GetBinContent(i+1);
        double nnlo = grNnlo.h->GetBinContent(i+1);
        double data = grData.h->GetBinContent(i+1);
        //cout <<"OLA " <<  i << grData.h->GetXaxis()->GetBinLabel(i) << " "<<grData.h->GetBinContent(i)<<" : " << grNlo.h->GetBinContent(i) << " "<< grNnlo.h->GetBinContent(i) << endl;
        tableTot << Names[i] << " $"<< nlo <<"^{+" << NLOsclUp <<"}_{-" << NLOsclDn <<"}$ & " 
                       << "$" << nnlo <<"^{+" << NNLOsclUp <<"}_{-" << NNLOsclDn <<"}$ &"
                       << "$^+{"<<NLOpdfUp <<"}_{" << NLOpdfDn <<"}$ &"
                       << "$^+{"<<NNLOpdfUp <<"}_{" << NNLOpdfDn <<"}$ " << endl;
    }
    tableTot.close();


	gStyle->SetEndErrorSize(4);
	grData.gr->Draw("e same");
	grData.grAll->Draw("pz same");


    gPad->SetTicks(1,1);

	//Plot Legend
	double w = 0.27;
	/*
	TLegend *leg1 = new TLegend( 0.35, 0.60, 0.35+w,0.88);
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	leg1->SetTextSize(myCan.fontS);
	leg1->SetHeader("Total cross sections");
	*/

	TGraphAsymmErrors *grTemp = new TGraphAsymmErrors;
	grTemp->SetLineWidth(0);
	grTemp->SetLineColorAlpha(kBlue-7, 0.0);
	grTemp->SetFillColorAlpha(kBlue-7, 0.6425); //0.7 org


	
	TLegend *leg2 = new TLegend( 0.30, 0.60, 0.30+2*w,0.88);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	leg2->SetTextSize(myCan.fontS);
	leg2-> SetNColumns(2);

	leg2->AddEntry((TObject*)0, "Total dijet cross sections", "");
	leg2->AddEntry((TObject*)0, "", "");

	leg2->AddEntry(grData.gr, "HERA Data", "ep");
	leg2->AddEntry(grNnlo.hRatio,    "NNLO (H1 Fit B)", "l");

	leg2->AddEntry(grNlo.hRatio,    "NLO (H1 Fit B)", "l");
	leg2->AddEntry(grTemp,   "NNLO scale unc.", "f");

	leg2->AddEntry(grNlo.grRatio  ,   "NLO scale unc.", "f");
	leg2->AddEntry(grNnlo.grAllRatio,"NNLO scale+DPDF unc.", "f");



	//leg1->Draw();
	leg2->Draw();

	myCan.UpdataFrame();
    DrawLatexUp(1.5, "#scale[2.0]{#font[72]{NNLOJET}}", -1, "l"); //RADEK

	////////////////////////////////////////////
	//Down Frame
	////////////////////////////////////////////

	myCan.downPad->cd();

    gPad->SetTicks(1,1);

	grNlo.grAllRatio->Draw("e2 same");
	grNlo.grRatio->Draw("e2 same");
	grNlo.hRatio->Draw("same");



	grNnlo.grAllRatio->Draw("e2 same");
	grNnlo.grRatio->Draw("e2 same");

	grNnlo.hRatio->SetBinError(1, 1e-6);
	grNnlo.hRatio->Draw("same");


	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");

	myCan.UpdataFrame();


	can->SaveAs( can->GetTitle() );
	can->Clear();

}





void PlotTotalDPDFs(TCanvas *can, Histogram *h1, Histogram *h2, Histogram *h3, Histogram *h4, Histogram *h5, Histogram *h6)
{

	vector<Histogram> histos = {*h1,*h2,*h3,*h4,*h5,*h6};


	HistoErr grData, grNlo, grNnlo;
	HistoErr grNnloSJ, grNnloJets, grNnloFitA, grNnloMRW, grNnloFitB;

	vector<HistoErr> grDataVec(6), grNloVec(6), grNnloVec(6);
	vector<HistoErr> grNnloSJVec(6), grNnloJetsVec(6), grNnloFitAVec(6), grNnloMRWVec(6), grNnloFitBVec(6);

	vector<double> range = { 1e30, -1e30, 1e30, -1e30 };
	for(unsigned i = 0; i < histos.size(); ++i) {
		vector<double> rangeTemp = histos[i].LoadHistogramsNloNnloData("nlo:Fit19-0-Q2pPt2-cc", grDataVec[i], grNloVec[i], grNnloVec[i], true);

		grNnloSJVec[i] = histos[i].LoadHistograms("sjNNLO", defFile,
										 "nnlo:zeusSJ-0-Q2pPt2-cc", "nlo:Fit19-0-Q2pPt2-cc");
		grNnloJetsVec[i] = histos[i].LoadHistograms("NNLOjets", defFile,
										 "nnlo:FitJets-0-Q2pPt2-cc", "nlo:Fit19-0-Q2pPt2-cc");
		grNnloFitAVec[i] = histos[i].LoadHistograms("NNLOfitA", defFile,
										 "nnlo:FitA-0-Q2pPt2-cc", "nlo:Fit19-0-Q2pPt2-cc");
		grNnloMRWVec[i] = histos[i].LoadHistograms("NNLOMRW", defFile,
										 "nnlo:MRW-0-Q2pPt2-cc", "nlo:Fit19-0-Q2pPt2-cc");

		grNnloFitBVec[i] = histos[i].LoadHistograms("NNLOFitB", defFile,
										 "nnlo:FitB-0-Q2pPt2-cc", "nlo:Fit19-0-Q2pPt2-cc");

		range[0] = min(range[0], rangeTemp[0]);
		range[2] = min(range[2], rangeTemp[2]);
		range[1] = max(range[1], rangeTemp[1]);
		range[3] = max(range[3], rangeTemp[3]);
	}
	grData = Merge("Data", grDataVec);
	grNlo  = Merge("NLO", grNloVec);
	grNnlo = Merge("NNLO", grNnloVec);

	grNnloSJ = Merge("NnloSJ", grNnloSJVec);
	grNnloJets = Merge("NnloJets", grNnloJetsVec);
	grNnloFitA = Merge("NnloFitA", grNnloFitAVec);
	grNnloMRW   = Merge("NnloMRW", grNnloMRWVec);
	grNnloFitB = Merge("NnloFitB", grNnloFitBVec);


	TString Names[] = { "FPS",   "VFPS", "LRG", "LRG_H1", "LRGH1820", "LRGZEUS" } ;

	//map<TString, TString> analMap;
	//analMap["FPS"] = "#splitline{#splitline{      H1}{(HERA #Iota#Iota)}}{    FPS}";
	//analMap["VFPS"] = "#splitline{#splitline{      H1}{(HERA #Iota#Iota)}}{   VFPS}";
	//analMap["LRG"] = "#splitline{#splitline{      H1}{(HERA #Iota#Iota)}}{    LRG}";
	//analMap["LRG_H1"] = "#splitline{#splitline{      H1}{(HERA #Iota)}}{    LRG}";
	//analMap["LRGH1820"] = "#splitline{#splitline{         H1}{    (HERA #Iota)}}{LRG, 820GeV}";
	//analMap["LRGZEUS"] = "#splitline{#splitline{   ZEUS}{(HERA #Iota)}}{    LRG}";


	for(int i = 0; i < 6; ++i) {
		grData.hRatio->GetXaxis()->SetBinLabel(i+1, analMap.at(Names[i]));
	}
	grData.hRatio->LabelsOption("h");
	grData.hRatio->LabelsDeflate();


	myCanvas myCan;


	myCan.range = range;
	myCan.yTitleUp =  "#sigma_{tot} [pb]";
	myCan.yTitleDown = "#sigma/#sigma_{NLO}";
	myCan.xTitle = "";
	myCan.fontS  = 0.065;



	//My start
	myCan.CreateFrame();

	//Setting NNLO style
	myCanvas::SetNNLOstyle(grNnlo);

	//Setting NNLO style
	myCanvas::SetNLOstyle(grNlo);

	grNnloSJ.h->SetLineColor(kRed); grNnloSJ.h->SetMarkerColor(kRed);  grNnloSJ.h->SetLineWidth(2);
	CopyStyle(grNnloSJ.h, grNnloSJ.hRatio);
	grNnloJets.h->SetLineColor(kMagenta); grNnloJets.h->SetMarkerColor(kMagenta); grNnloJets.h->SetLineWidth(2);
	CopyStyle(grNnloJets.h, grNnloJets.hRatio);
	grNnloFitA.h->SetLineColor(kBlack); grNnloFitA.h->SetMarkerColor(kBlack); grNnloFitA.h->SetLineWidth(2);
	CopyStyle(grNnloFitA.h, grNnloFitA.hRatio);
	//grNnloMRW.h->SetLineColor(kOrange); grNnloMRW.h->SetMarkerColor(kOrange); grNnloMRW.h->SetLineWidth(2);
	//CopyStyle(grNnloMRW.h, grNnloMRW.hRatio);

	grNnloFitB.h->SetLineColor(kOrange); grNnloFitB.h->SetMarkerColor(kOrange); grNnloFitB.h->SetLineWidth(2);
	CopyStyle(grNnloFitB.h, grNnloFitB.hRatio);
	//myCanvas::SetNNLOstyle(grNnloFitB);


	//Setting Data style and frame
	gStyle->SetEndErrorSize(16.);
	myCan.SetDataStyle(grData, false, true, 800);
	grData.h->GetYaxis()->SetNoExponent();

	grData.h->GetXaxis()->SetTickLength(0.0);
	grData.hRatio->GetXaxis()->SetLabelOffset(0.03);
	grData.hRatio->GetXaxis()->SetLabelSize(myCan.fontS * 1.4722 * 0.9);


	double pointSize = 1.2;
	grData.grAllRatio->SetMarkerStyle(20);
	grData.grAllRatio->SetMarkerSize(pointSize);
	grData.grAll->SetMarkerStyle(20);
	grData.grAll->SetMarkerSize(pointSize);





	////////////////////////////////////////////
	//Up Frame
	////////////////////////////////////////////

	myCan.upPad->cd();
    gPad->SetTicks(1,1);


	if(!true) {
		grData.h->SetMinimum(0);
		grData.h->SetMaximum(range[1] * 1.1);
	}
	else {
		grData.h->SetMinimum(range[0] * 0.8);
		grData.h->SetMaximum(range[1] * 2.2);
	}

	//grData.h->Draw("AXIS");


	grNlo.h->SetBinError(1, 1e-6);
	grNnlo.h->SetBinError(1, 1e-6);
	grNnloSJ.h->SetBinError(1, 1e-6);
	grNnloJets.h->SetBinError(1, 1e-6);
	//grNnloFitA.h->SetBinError(1, 1e-6);
	grNnloMRW.h->SetBinError(1, 1e-6);
	grNnloFitB.h->SetBinError(1, 1e-6);






	grNlo.h->Draw("same");

	grNnlo.grAll->Draw("e2 same");
	grNnlo.gr->Draw("e2 same");
	grNnlo.h->Draw("same");




	grNnloSJ.h->Draw("same");
	grNnloJets.h->Draw("same");
	//grNnloFitA.h->Draw("same");
	//grNnloMRW.h->Draw("same");
	grNnloFitB.h->Draw("same");
    //for(int i = 1; i < grNnloMRW.h->GetNbinsX(); ++i) {
        //cout << "LULULU " << grNnlo.h->GetBinContent(i) <<" "<< grNnloMRW.h->GetBinContent(i) << endl;
    //}

	gStyle->SetEndErrorSize(4);
	grData.gr->Draw("e same");
	grData.grAll->Draw("pz same");



	//Plot Legend
	/*
	TLegend *leg = new TLegend( 0.5, 0.60, 0.95,0.88);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetHeader("Total cross sections");
	leg->AddEntry(grData.gr, "Data", "ep");
	leg->AddEntry(grNnlo.grAll, "NNLO (H1 2006 Fit B)" , "fl");
	leg->AddEntry(grNnloFitA.h, "NNLO (H1 2006 Fit A)" , "l");
	leg->AddEntry(grNnloSJ.h,   "NNLO (ZEUS SJ)" , "l");
	leg->AddEntry(grNnloJets.h, "NNLO (H1 2007 Fit Jets)" , "l");
	*/

	TGraphAsymmErrors *grTemp = new TGraphAsymmErrors;
	grTemp->SetLineWidth(0);
	grTemp->SetLineColorAlpha(kBlue-7, 0.0);
	grTemp->SetFillColorAlpha(kBlue-7, 0.6425); //0.7 org



	double w = 0.32;
	TLegend *leg = new TLegend( 0.29, 0.50, 0.29+2*w,0.88);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.95*myCan.fontS);
	leg-> SetNColumns(2);

	leg->AddEntry((TObject*)0, "Total dijet cross sec.", "");
	//leg->AddEntry((TObject*)0, "", "");
	leg->AddEntry(grData.gr, "HERA Data", "ep");

	leg->AddEntry(grNnloFitB.h, "NNLO (H1 Fit B)" , "l");
	leg->AddEntry(grNlo.hRatio,    "NLO (H1 Fit 19)", "l");

	leg->AddEntry(grNnloJets.h, "NNLO (H1 Fit Jets)" , "l");
	leg->AddEntry(grNnlo.hRatio,"NNLO (H1 Fit 19)", "l");

	leg->AddEntry(grNnloSJ.h,   "NNLO (ZEUS SJ)" , "l");
	leg->AddEntry(grTemp /*grNnlo.grRatio*/,   "DPDF unc.", "f");

	//leg->AddEntry(grNnloMRW.h,   "NNLO (MRW)" , "l");
	leg->AddEntry((TObject*)nullptr,   "" , "");
	leg->AddEntry(grNnlo.grAllRatio,"DPDF+scale unc.", "f");


	leg->Draw();

	myCan.UpdataFrame();
    DrawLatexUp(1.5, "#scale[2.0]{#font[72]{NNLOJET}}", -1, "l"); //RADEK

	////////////////////////////////////////////
	//Down Frame
	////////////////////////////////////////////

	myCan.downPad->cd();
    gPad->SetTicks(1,1);

	grNlo.hRatio->Draw("same");

	grNnlo.grAllRatio->Draw("e2 same");
	grNnlo.grRatio->Draw("e2 same");
	grNnlo.hRatio->SetBinError(1, 1e-6);
	grNnlo.hRatio->Draw("same");


	grNnloSJ.hRatio->SetBinError(1, 1e-6);
	grNnloJets.hRatio->SetBinError(1, 1e-6);
	grNnloFitA.hRatio->SetBinError(1, 1e-6);
	grNnloMRW.hRatio->SetBinError(1, 1e-6);
	grNnloFitB.hRatio->SetBinError(1, 1e-6);


	grNnloSJ.hRatio->Draw("same");
	grNnloJets.hRatio->Draw("same");
	//grNnloFitA.hRatio->Draw("same");
	//grNnloMRW.hRatio->Draw("same");
	grNnloFitB.hRatio->Draw("same");


	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");

	myCan.UpdataFrame();


	can->SaveAs( can->GetTitle() );
	can->Clear();

}


void plotScaleStudiesTotal(TCanvas *can, Histogram *h1, Histogram *h2, Histogram *h3, Histogram *h4, Histogram *h5, Histogram *h6)
{

	vector<Histogram> histos = {*h1,*h2,*h3,*h4,*h5,*h6};


	HistoErr grData, grNlo, grNnlo;
	HistoErr grNnloQ2, grNnloPt, grNnloQ2Pt;

	vector<HistoErr> grDataVec(6), grNloVec(6), grNnloVec(6);
	vector<HistoErr> grNnloQ2Vec(6), grNnloPtVec(6), grNnloQ2PtVec(6);

	vector<double> range = { 1e30, -1e30, 1e30, -1e30 };
	for(unsigned i = 0; i < histos.size(); ++i) {
		vector<double> rangeTemp = histos[i].LoadHistogramsNloNnloData("nlo:FitB-0-Q2pPt2-cc", grDataVec[i], grNloVec[i], grNnloVec[i], false);

		grNnloQ2Vec[i]  = histos[i].LoadHistograms("NNLOq2", defFile,
										 "nnlo:FitB-0-Q2-cc", "nlo:FitB-0-Q2pPt2-cc");
		grNnloPtVec[i] = histos[i].LoadHistograms("NNLOpt", defFile,
										 "nnlo:FitB-0-Pt2-cc", "nlo:FitB-0-Q2pPt2-cc");
		grNnloQ2PtVec[i] = histos[i].LoadHistograms("NNLOq2pt", defFile,
										 "nnlo:FitB-0-0.25Q2pPt2-cc", "nlo:FitB-0-Q2pPt2-cc");


		range[0] = min(range[0], rangeTemp[0]);
		range[2] = min(range[2], rangeTemp[2]);
		range[1] = max(range[1], rangeTemp[1]);
		range[3] = max(range[3], rangeTemp[3]);
	}
	grData = Merge("Data", grDataVec);
	grNlo  = Merge("NLO", grNloVec);
	grNnlo = Merge("NNLO", grNnloVec);

	grNnloQ2 = Merge("NnloQ2", grNnloQ2Vec);
	grNnloPt = Merge("NnloPt", grNnloPtVec);
	grNnloQ2Pt = Merge("NnloQ2Pt", grNnloQ2PtVec);


	TString Names[] = { "FPS",   "VFPS", "LRG", "LRG_H1", "LRGH1820", "LRGZEUS" } ;

	//map<TString, TString> analMap;
	//analMap["FPS"] = "#splitline{#splitline{      H1}{(HERA #Iota#Iota)}}{    FPS}";
	//analMap["VFPS"] = "#splitline{#splitline{      H1}{(HERA #Iota#Iota)}}{   VFPS}";
	//analMap["LRG"] = "#splitline{#splitline{      H1}{(HERA #Iota#Iota)}}{    LRG}";
	//analMap["LRG_H1"] = "#splitline{#splitline{      H1}{(HERA #Iota)}}{    LRG}";
	//analMap["LRGH1820"] = "#splitline{#splitline{         H1}{    (HERA #Iota)}}{LRG, 820GeV}";
	//analMap["LRGZEUS"] = "#splitline{#splitline{   ZEUS}{(HERA #Iota)}}{    LRG}";


	for(int i = 0; i < 6; ++i) {
		grData.hRatio->GetXaxis()->SetBinLabel(i+1, analMap.at(Names[i]));
	}
	grData.hRatio->LabelsOption("h");
	grData.hRatio->LabelsDeflate();


	myCanvas myCan;


	myCan.range = range;
	myCan.yTitleUp =  "#sigma_{tot} [pb]";
	myCan.yTitleDown = "#sigma/#sigma_{NLO}";
	myCan.xTitle = "";
	myCan.fontS  = 0.065;



	//My start
	myCan.CreateFrame();

	//Setting NNLO style
	myCanvas::SetNNLOstyle(grNnlo);
	myCanvas::SetNLOstyle(grNlo);

	grNnloQ2.h->SetLineColor(kRed); grNnloQ2.h->SetMarkerColor(kRed); grNnloQ2.h->SetLineWidth(2);
	CopyStyle(grNnloQ2.h, grNnloQ2.hRatio);
	grNnloPt.h->SetLineColor(kBlack); grNnloPt.h->SetMarkerColor(kBlack);  grNnloPt.h->SetLineWidth(2);
	CopyStyle(grNnloPt.h, grNnloPt.hRatio);
	grNnloQ2Pt.h->SetLineColor(kCyan); grNnloQ2Pt.h->SetMarkerColor(kCyan); grNnloQ2Pt.h->SetLineWidth(2);
	CopyStyle(grNnloQ2Pt.h, grNnloQ2Pt.hRatio);


	//Setting Data style and frame
	myCan.SetDataStyle(grData, false, true);
	grData.h->GetYaxis()->SetNoExponent();

	grData.h->GetXaxis()->SetTickLength(0.0);
	grData.hRatio->GetXaxis()->SetLabelOffset(0.03);
	grData.hRatio->GetXaxis()->SetLabelSize(myCan.fontS * 1.4722 * 0.9);


	double pointSize = 1.2;
	gStyle->SetEndErrorSize(8.);
	grData.grAllRatio->SetMarkerStyle(20);
	grData.grAllRatio->SetMarkerSize(pointSize);
	grData.grAll->SetMarkerStyle(20);
	grData.grAll->SetMarkerSize(pointSize);




	////////////////////////////////////////////
	//Up Frame
	////////////////////////////////////////////

	myCan.upPad->cd();
    gPad->SetTicks(1,1);


	if(!true) {
		grData.h->SetMinimum(0);
		grData.h->SetMaximum(range[1] * 1.1);
	}
	else {
		grData.h->SetMinimum(range[0] * 0.8);
		grData.h->SetMaximum(range[1] * 2.2);
	}

	grNlo.h->SetBinError(1, 1e-6);
	grNnlo.h->SetBinError(1, 1e-6);

	grNlo.h->Draw("same");

	grNnlo.grAll->Draw("e2 same");
	grNnlo.gr->Draw("e2 same");
	grNnlo.h->Draw("same");


	grNnloQ2.h->SetBinError(1, 1e-6);
	grNnloPt.h->SetBinError(1, 1e-6);
	grNnloQ2Pt.h->SetBinError(1, 1e-6);

	grNnloQ2.h->Draw("same");
	grNnloPt.h->Draw("same");
	grNnloQ2Pt.h->Draw("same");

	gStyle->SetEndErrorSize(4);
	grData.gr->Draw("e same");
	grData.grAll->Draw("pz same");



	//Plot Legend
	/*
	TLegend *leg = new TLegend( 0.5, 0.60, 0.95,0.88);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetHeader("Total cross sections");
	leg->AddEntry(grData.gr, "Data", "ep");
	leg->AddEntry(grNnlo.grAll, "NNLO (#mu^{2}=Q^{2}+p_{T}^2)", "fl");
	leg->AddEntry(grNnloQ2.h,   "NNLO (#mu^{2}=Q^{2})" , "l");
	leg->AddEntry(grNnloPt.h,   "NNLO (#mu^{2}=p_{T}^{2})"  , "l");
	leg->AddEntry(grNnloQ2Pt.h, "NNLO (#mu^{2}=#frac{Q^{2}}{4}+p_{T}^{2})" , "l");
	leg->Draw();
	*/
	TGraphAsymmErrors *grTemp = new TGraphAsymmErrors;
	grTemp->SetLineWidth(0);
	grTemp->SetLineColorAlpha(kBlue-7, 0.0);
	grTemp->SetFillColorAlpha(kBlue-7, 0.6425); //0.7 org


	double w = 0.30;
	TLegend *leg = new TLegend( 0.30, 0.47, 0.30+2*w,0.88);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(myCan.fontS*0.95);
	leg-> SetNColumns(2);

	leg->AddEntry((TObject*)0, "Total dijet cross sections", "");
	leg->AddEntry((TObject*)0, "", "");

	leg->AddEntry(grData.gr, "HERA Data", "ep");
	leg->AddEntry(grNlo.hRatio,    "NLO (#mu^{2}=Q^{2}+p_{T}^{2})", "l");

	leg->AddEntry(grNnloQ2Pt.h,   "NNLO (#mu^{2}=#frac{Q^{2}}{4}+p_{T}^{2})" , "l");
	leg->AddEntry(grNnlo.hRatio,"NNLO (#mu^{2}=Q^{2}+p_{T}^{2})", "l");

	leg->AddEntry(grNnloPt.h, "NNLO (#mu^{2}=p_{T}^{2})" , "l");
	leg->AddEntry(grTemp /*grNnlo.grRatio*/,   "Scale unc.", "f");

	leg->AddEntry(grNnloQ2.h, "NNLO (#mu^{2}=Q^{2})" , "l");
	leg->AddEntry(grNnlo.grAllRatio,"Scale+DPDF unc.", "f");


	leg->Draw();









	myCan.UpdataFrame();
    DrawLatexUp(1.5, "#scale[2.0]{#font[72]{NNLOJET}}", -1, "l"); //RADEK

	////////////////////////////////////////////
	//Down Frame
	////////////////////////////////////////////

	myCan.downPad->cd();
    gPad->SetTicks(1,1);


	grNlo.hRatio->SetBinError(1, 1e-6);
	grNnlo.hRatio->SetBinError(1, 1e-6);

	grNlo.hRatio->Draw("same");

	grNnlo.grAllRatio->Draw("e2 same");
	grNnlo.grRatio->Draw("e2 same");
	grNnlo.hRatio->Draw("same");

	grNnloQ2.hRatio->SetBinError(1, 1e-6);
	grNnloPt.hRatio->SetBinError(1, 1e-6);
	grNnloQ2Pt.hRatio->SetBinError(1, 1e-6);



	grNnloQ2.hRatio->Draw("same");
	grNnloPt.hRatio->Draw("same");
	grNnloQ2Pt.hRatio->Draw("same");


	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");

	myCan.UpdataFrame();


	can->SaveAs( can->GetTitle() );
	can->Clear();

}


















void Histogram::plotScaleChoices()
{

	cout << "Plotting : \""<< var_name<<"\" " << nBins << " "<< xMin.size() << endl;

	Double_t *hBins = new Double_t[nBins+1];
	for(int i = 0; i < nBins; ++i)
		hBins[i] = xMin[i];
	hBins[nBins] = xMax[nBins-1];

	int hash = rand();


	//Count NLO

	vector<double>  Nlo,  NloUp, NloDown;
	vector<double>  NNlo,  NNloUp, NNloDown;
	vector<double> NloA[4], NNloA[4];
	vector<TString> NloNames, NNloNames;
	Theory *nloNow;
	for(auto & th : theories) {
		TString name = th.first;

		if( name.Contains(defFile) ) {
			if( name.BeginsWith("nlo:") ) {
				//cout << "RADEK " << name << endl;
				if(name.EndsWith("FitB-0-Q2pPt2-cc"))
					Nlo  = th.second.xsc ;
				if(name.EndsWith("FitB-0-Q2pPt2-uu"))
					NloUp= th.second.xsc;
				if(name.EndsWith("FitB-0-Q2pPt2-dd"))
					NloDown =  th.second.xsc;

				if(name.EndsWith("FitB-0-Q2-cc"))
					NloA[0] =  th.second.xsc;
				if(name.EndsWith("FitB-0-Pt2-cc"))
					NloA[1] =  th.second.xsc;
				if(name.EndsWith("FitB-0-0.25Q2pPt2-cc"))
					NloA[2] =  th.second.xsc;

				NloNames.insert(NloNames.begin(),  name);
				nloNow = &th.second;
			}
			if( name.BeginsWith("nnlo:")  ) {

				if(name.EndsWith("FitB-0-Q2pPt2-cc"))
					NNlo= th.second.xsc ;
				if(name.EndsWith("FitB-0-Q2pPt2-uu"))
					NNloUp =  th.second.xsc;
				if(name.EndsWith("FitB-0-Q2pPt2-dd"))
					NNloDown =  th.second.xsc;

				if(name.EndsWith("FitB-0-Q2-cc"))
					NNloA[0] =  th.second.xsc;
				if(name.EndsWith("FitB-0-Pt2-cc"))
					NNloA[1] =  th.second.xsc;
				if(name.EndsWith("FitB-0-0.25Q2pPt2-cc"))
					NNloA[2] =  th.second.xsc;


				NNloNames.insert(NNloNames.begin(), name);
			}
		}
	}
	//exit(0);

	//vector<double> &Nlo     = nloNow->xsc;
	//vector<double> &NloUp   = nloNow->xscUp;
	//vector<double> &NloDown = nloNow->xscDown;

	//vector<double> &NNlo     = nnloNow->xsc;





	TH1D *hData = new TH1D( SF("hData%d",hash), "data", nBins, hBins);

	TH1D *hNLO, *hNLORatio, *hNNLO, *hNNLORatio;
	TH1D *hNLOa[4], *hNNLOa[4];
	TH1D *hNLOaRatio[4], *hNNLOaRatio[4];

	hNLO      = new TH1D( SF("hNLO%d",hash), "NLO", nBins, hBins);
	hNLORatio = new TH1D( SF("hNLORatio%d",hash), "NLORatio", nBins, hBins);
	hNNLO     = new TH1D( SF("hNNLO%d",hash), "NNLO", nBins, hBins);
	hNNLORatio= new TH1D( SF("hNNLORatio%d",hash), "NNLORatio", nBins, hBins);

	for(int i = 0; i < 3; ++i) {
		hNLOa[i]      = new TH1D( SF("hNLOa%d_%d",i,hash), "NLOa", nBins, hBins);
		hNNLOa[i]     = new TH1D( SF("hNNLOa%d_%d",i,hash), "NNLOa", nBins, hBins);
		hNLOaRatio[i]      = new TH1D( SF("hNLOaRatio%d_%d",i,hash), "NLOaRat", nBins, hBins);
		hNNLOaRatio[i]     = new TH1D( SF("hNNLOaRatio%d_%d",i,hash), "NNLOaRat", nBins, hBins);
	}

	/*
	TH1D *hNLO = new TH1D( SF("hNLO%d",hash), "NLO", nBins, hBins);
	TH1D *hNLORatio = new TH1D( SF("hNLORatio%d",hash), "NLORatio", nBins, hBins);

	TH1D *hNNLO = new TH1D( SF("hNNLO%d",hash), "NNLO", nBins, hBins);
	TH1D *hNNLORatio = new TH1D( SF("hNNLORatio%d",hash), "NNLORatio", nBins, hBins);
	*/


	TH1D *hDataRatio = new TH1D( SF("hDataRatio%d",hash), "DataRatio", nBins, hBins);

	TGraphAsymmErrors *gData      = new TGraphAsymmErrors(nBins);
	TGraphAsymmErrors *gDataAll   = new TGraphAsymmErrors(nBins);
	TGraphAsymmErrors *gDataRatio = new TGraphAsymmErrors(nBins);
	TGraphAsymmErrors *gDataAllRatio = new TGraphAsymmErrors(nBins);

	TGraphAsymmErrors *gNLO       = new TGraphAsymmErrors(nBins);
	TGraphAsymmErrors *gNLOscale  = new TGraphAsymmErrors(nBins);
	TGraphAsymmErrors *gNLORatio  = new TGraphAsymmErrors(nBins);
	TGraphAsymmErrors *gNLOscaleRatio=new TGraphAsymmErrors(nBins);


	TGraphAsymmErrors *gNNLO       = new TGraphAsymmErrors(nBins);
	TGraphAsymmErrors *gNNLOscale  = new TGraphAsymmErrors(nBins);
	TGraphAsymmErrors *gNNLORatio  = new TGraphAsymmErrors(nBins);
	TGraphAsymmErrors *gNNLOscaleRatio=new TGraphAsymmErrors(nBins);

	double dispMax = -1e30;
	double dispMin = +1e30;
	double ratMax = -100;
	double ratMin = +100;


	for(int i = 0; i < nBins; ++i) {
		hData->SetBinContent(i+1, data[i]);
		hData->SetBinError(i+1, dataStatErr[i]);

		gData->SetPoint(i, hData->GetBinCenter(i+1), data[i] );
		gData->SetPointError(i, 0, 0, dataStatErr[i], dataStatErr[i] );

		gDataAll->SetPoint(i, hData->GetBinCenter(i+1), data[i] );
		double totEr = hypot(dataStatErr[i], dataSystErr[i]);
		gDataAll->SetPointError(i, 0, 0,  totEr, totEr );

		hNLO->SetBinContent(i+1, Nlo[i]);

		hNLORatio->SetBinContent(i+1, Nlo[i]/Nlo[i] );
		hNLORatio->SetBinError(i+1, 0 );

		hNNLO->SetBinContent(i+1, NNlo[i]);
		hNNLORatio->SetBinContent(i+1, NNlo[i]/Nlo[i] );
		hNNLORatio->SetBinError(i+1, 0 );


		//Alternative scales
		for(int k = 0; k < 3; ++k) {
			hNLOa[k]->SetBinContent(i+1, NloA[k][i]);
			hNLOa[k]->SetBinError(i+1, 0);
			hNNLOa[k]->SetBinContent(i+1, NNloA[k][i]);
			hNNLOa[k]->SetBinError(i+1, 0);
			hNLOaRatio[k]->SetBinContent(i+1, NloA[k][i]/Nlo[i]);
			hNLOaRatio[k]->SetBinError(i+1, 0);
			hNNLOaRatio[k]->SetBinContent(i+1, NNloA[k][i]/Nlo[i]);
			hNNLOaRatio[k]->SetBinError(i+1, 0);
		}





		hDataRatio->SetBinContent(i+1, data[i]/Nlo[i] );
		hDataRatio->SetBinError(i+1, dataStatErr[i]/Nlo[i] );

		gDataRatio->SetPoint(i, hData->GetBinCenter(i+1), data[i]/Nlo[i] );
		gDataRatio->SetPointError(i, 0, 0, dataStatErr[i]/Nlo[i],  dataStatErr[i]/Nlo[i]  );

		gDataAllRatio->SetPoint(i, hData->GetBinCenter(i+1), data[i]/Nlo[i] );
		gDataAllRatio->SetPointError(i, 0, 0,  totEr/Nlo[i], totEr/Nlo[i] );

		//LOAD NLO

		double nloSyst = 0;//nloNow->GetSystTot(i);
		double nloH    = max({0., NloUp[i] - Nlo[i], NloDown[i] - Nlo[i]});
		double nloL    =-min({0., NloUp[i] - Nlo[i], NloDown[i] - Nlo[i]});

		double nloTotU = hypot(nloSyst,nloL);
		double nloTotD = hypot(nloSyst,nloH);

		gNLOscale->SetPoint(i, hData->GetBinCenter(i+1), Nlo[i] );
		gNLOscale->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nloL, nloH );

		gNLO->SetPoint(i, hData->GetBinCenter(i+1), Nlo[i] );
		gNLO->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2,   nloTotU, nloTotD );


		gNLORatio->SetPoint(i, hData->GetBinCenter(i+1), Nlo[i]/Nlo[i] );
		gNLORatio->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nloTotU/Nlo[i], nloTotD/Nlo[i] );

		gNLOscaleRatio->SetPoint(i, hData->GetBinCenter(i+1), Nlo[i]/Nlo[i] );
		gNLOscaleRatio->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nloL/Nlo[i], nloH/Nlo[i] );

		//LOAD NNLO

		double nnloSyst = 0;//nloNow->GetSystTot(i);
		double nnloH    = max({0., NNloUp[i] - NNlo[i], NNloDown[i] - NNlo[i]});
		double nnloL    =-min({0., NNloUp[i] - NNlo[i], NNloDown[i] - NNlo[i]});

		double nnloTotU = hypot(nnloSyst,nnloL);
		double nnloTotD = hypot(nnloSyst,nnloH);

		gNNLOscale->SetPoint(i, hData->GetBinCenter(i+1), NNlo[i] );
		gNNLOscale->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nnloL, nnloH );

		gNNLO->SetPoint(i, hData->GetBinCenter(i+1), NNlo[i] );
		gNNLO->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2,   nnloTotU, nnloTotD );


		gNNLORatio->SetPoint(i, hData->GetBinCenter(i+1), NNlo[i]/Nlo[i] );
		gNNLORatio->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nnloTotU/Nlo[i], nnloTotD/Nlo[i] );

		gNNLOscaleRatio->SetPoint(i, hData->GetBinCenter(i+1), NNlo[i]/Nlo[i] );
		gNNLOscaleRatio->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nnloL/Nlo[i], nnloH/Nlo[i] );



		dispMax = max(dispMax, nloTotD + Nlo[i] );
		dispMax = max(dispMax, nnloTotD + NNlo[i] );
		dispMax = max(dispMax, data[i] + totEr );

		dispMin = min(dispMin, -nloTotU + Nlo[i] );
		dispMin = min(dispMin, -nnloTotU + NNlo[i] );
		dispMin = min(dispMin, data[i] - totEr );


		
		ratMax  = max(ratMax, Nlo[i]/Nlo[i]  + nloTotD/Nlo[i]);
		ratMax  = max(ratMax, NNlo[i]/Nlo[i] + nnloTotD/Nlo[i]);
		ratMax  = max(ratMax,  data[i]/Nlo[i] + totEr/Nlo[i] );

		ratMin  = min(ratMin, Nlo[i]/Nlo[i]  - nloTotU/Nlo[i]);
		ratMin  = min(ratMin, NNlo[i]/Nlo[i] - nnloTotU/Nlo[i]);
		ratMin  = min(ratMin, data[i]/Nlo[i] - totEr/Nlo[i] );


	}

	double leftMargin = 0.15;
	double rightMargin = 0.03;
	double bottomMargin = 0.35;
	double fontS = 0.08;

	double pointSize = 0.8;

	TVirtualPad * currpad = 0;

	currpad = gPad;//canvas->cd(1);

	//currpad -> Divide(1,3,0.,0.,0); // rozdelim prvni pad na tri podsebou
	//currpad -> Divide(1,2,0.,0.,0);

	// a ty se zacnou znovu cislovat od jedne

	//currpad -> cd(1);// a prepnu se do prvniho padu 
	TPad *upPad = new TPad( SF("upPad%d",hash), "upPad", 0, 0.5, 1, 1 );
	upPad->SetBottomMargin(0);
	upPad->SetLeftMargin(leftMargin);
	upPad->SetRightMargin(rightMargin);
	upPad->Draw();
	upPad->cd();

	if(style.isLogY) gPad->SetLogy();
	if(style.isLogX) gPad->SetLogx();

	//gPad -> SetRightMargin(0.02);

	hData->SetLineColor(kBlack);
	hData->Draw("AXIS");
	hData->GetYaxis()->SetTitle(yTitle);
	hData->GetYaxis()->SetTitleSize(fontS);
	hData->GetYaxis()->SetTitleOffset(0.85);
	hData->GetYaxis()->SetLabelSize(fontS);
	if(!style.isLogY) {
		hData->SetMinimum(0);
		hData->SetMaximum(dispMax * 1.1);
	}
	else {
		hData->SetMinimum(dispMin * 0.8);
		hData->SetMaximum(dispMax * 1.2);
	}


	gNLO->SetLineColor(kWhite);
	gNLO->SetFillColor(kOrange);
	gNLO->Draw("e2 same");

	gNLOscale->SetLineColor(kWhite);
	gNLOscale->SetFillColor(kOrange);
	gNLOscale->Draw("e2 same");

	hNLO->SetLineColor(kWhite);
	hNLO->Draw("same");



	gNNLO->SetLineColor(kGreen);
	gNNLO->SetFillColor(kGreen);
	gNNLO->SetFillStyle(3004);
	gNNLO->Draw("e2 same");


	gNNLOscale->SetLineColor(kGreen);
	gNNLOscale->SetFillColor(kGreen);
	gNNLOscale->SetFillStyle(3004);
	gNNLOscale->Draw("e2 same");

	hNNLO->SetLineColor(kGreen);
	hNNLO->SetLineStyle(1);
	hNNLO->Draw("same");


	hNLOa[0]->SetLineColor(kRed);
	hNLOa[0]->SetLineStyle(1);
	hNLOa[0]->Draw("same");
	hNLOa[1]->SetLineColor(kRed);
	hNLOa[1]->SetLineStyle(2);
	hNLOa[1]->Draw("same");
	hNLOa[2]->SetLineColor(kRed);
	hNLOa[2]->SetLineStyle(3);
	hNLOa[2]->Draw("same");

	hNNLOa[0]->SetLineColor(kBlue);
	hNNLOa[0]->SetLineStyle(1);
	hNNLOa[0]->Draw("same");
	hNNLOa[1]->SetLineColor(kBlue);
	hNNLOa[1]->SetLineStyle(2);
	hNNLOa[1]->Draw("same");
	hNNLOa[2]->SetLineColor(kBlue);
	hNNLOa[2]->SetLineStyle(3);
	hNNLOa[2]->Draw("same");






	gStyle->SetEndErrorSize(4);
	gData->SetMarkerStyle(20);
	gData->SetMarkerSize(pointSize);
	gData->Draw("e same");

	gDataAll->SetMarkerStyle(20);
	gDataAll->SetMarkerSize(pointSize);
	gDataAll->Draw("pz same");


	//Plot Legend
	TLegend *leg = new TLegend( style.legX+0.5-0.5, style.legY+0.60-0.6,
	                            style.legX+0.95-0.5,style.legY+0.88-0.6 );
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	TString tag =  theor_path(theor_path.First('/')+1, theor_path.Last('/')-theor_path.First('/')-1 );
	leg->SetHeader(Title);
	leg->AddEntry(gData,    SF("H1 %s data", tag.Data()), "ep");

	leg->AddEntry(gNLO,     SF("NLO x %.2f (Q^{2}+p_{T}^{2})",dissFactor) , "fl");
	leg->AddEntry(hNLOa[0], SF("NLO x %.2f (Q^{2})",dissFactor) , "l");
	leg->AddEntry(hNLOa[1], SF("NLO x %.2f (p_{T}^{2})",dissFactor) , "l");
	leg->AddEntry(hNLOa[2], SF("NLO x %.2f (Q^{2}/4+p_{T}^{2})",dissFactor) , "l");


	leg->AddEntry(gNNLOscale,SF("NNLO x %.2f (Q^{2}+p_{T}^{2})", dissFactor) , "fl");
	leg->AddEntry(hNNLOa[0], SF("NNLO x %.2f (Q^{2})",dissFactor) , "l");
	leg->AddEntry(hNNLOa[1], SF("NNLO x %.2f (p_{T}^{2})",dissFactor) , "l");
	leg->AddEntry(hNNLOa[2], SF("NNLO x %.2f (Q^{2}/4+p_{T}^{2})",dissFactor) , "l");




	/*
	for(unsigned i = 1; i < Nlo.size(); ++i) {
		TString name = NloNames[i](0,NloNames[i].Last('.') );
		name = name( name.Last('.')+1, 10000);
		leg->AddEntry(hNLO[i], SF("NLO %s x %.2f",name.Data(), dissFactor) , "l");
	}

	//leg->AddEntry(hNNLO[0], SF("NNLO H12006 Fit-B x %.2f",dissFactor) , "l");
	for(unsigned i = 0; i < Nlo.size(); ++i) {
		TString name = NNloNames[i](0,NNloNames[i].Last('.') );
		name = name( name.Last('.')+1, 10000);
		leg->AddEntry(hNNLO[i], SF("NNLO %s x %.2f",name.Data(), dissFactor) , "l");
	}
	*/

	leg->Draw();

	gPad->Update();
	gPad->RedrawAxis();
	TFrame *frU = (TFrame *) gPad->FindObject("TFrame");
	frU->SetFillStyle(0);
	frU->Draw();
	gPad->Update();


	currpad->cd();
	TPad *downPad = new TPad( SF("downPad%d",hash), "downPad", 0, 0.0, 1, 0.5 );
	downPad->SetTopMargin(0);
	downPad->SetLeftMargin(leftMargin);
	downPad->SetRightMargin(rightMargin);
	downPad->SetBottomMargin(bottomMargin);
	downPad->Draw();
	downPad->cd();
	//gPad -> SetRightMargin(0.02);
	if(style.isLogX) gPad->SetLogx();

	double Max = max( hDataRatio->GetMaximum(), hNLORatio->GetMaximum() );
	hDataRatio->SetMaximum(Max);

	hDataRatio->GetYaxis()->SetNdivisions(204);

	hDataRatio->Draw("AXIS");
	hDataRatio->SetMinimum( max(0.,ratMin - 0.1*ratMax) );
	hDataRatio->SetMaximum( 1.1*ratMax );


	hDataRatio->GetXaxis()->SetTitle(xTitle);
	hDataRatio->GetXaxis()->SetTitleSize(fontS * 5./5);

	hDataRatio->GetYaxis()->SetLabelSize(fontS * 5./5);
	hDataRatio->GetXaxis()->SetLabelSize(fontS * 5./5);

	hDataRatio->GetYaxis()->SetTitle("NLO/Data");
	hDataRatio->GetYaxis()->SetTitleSize(fontS * 5./5);
	hDataRatio->GetYaxis()->SetTitleOffset(0.84);

	gNLORatio->SetFillColor(kOrange);
	gNLORatio->Draw("e2 same");

	gNLOscaleRatio->SetFillColor(kOrange);
	gNLOscaleRatio->Draw("e2 same");


	hNLORatio->SetLineColor(kWhite);
	hNLORatio->Draw("same");

	hNNLORatio->SetLineColor(kBlack);
	hNNLORatio->SetLineStyle(2);
	hNNLORatio->Draw("same");


	gNNLORatio->SetLineColor(kGreen);
	gNNLORatio->SetFillColor(kGreen);
	gNNLORatio->SetFillStyle(3004);
	gNNLORatio->Draw("e2 same");


	gNNLOscaleRatio->SetLineColor(kGreen);
	gNNLOscaleRatio->SetFillColor(kGreen);
	gNNLOscaleRatio->SetFillStyle(3004);
	gNNLOscaleRatio->Draw("e2 same");

	hNNLORatio->SetLineColor(kGreen);
	hNNLORatio->SetLineStyle(1);
	hNNLORatio->Draw("same");

	hNLOaRatio[0]->SetLineColor(kRed);
	hNLOaRatio[0]->SetLineStyle(1);
	hNLOaRatio[0]->Draw("same");
	hNLOaRatio[1]->SetLineColor(kRed);
	hNLOaRatio[1]->SetLineStyle(2);
	hNLOaRatio[1]->Draw("same");
	hNLOaRatio[2]->SetLineColor(kRed);
	hNLOaRatio[2]->SetLineStyle(3);
	hNLOaRatio[2]->Draw("same");

	hNNLOaRatio[0]->SetLineColor(kBlue);
	hNNLOaRatio[0]->SetLineStyle(1);
	hNNLOaRatio[0]->Draw("same");
	hNNLOaRatio[1]->SetLineColor(kBlue);
	hNNLOaRatio[1]->SetLineStyle(2);
	hNNLOaRatio[1]->Draw("same");
	hNNLOaRatio[2]->SetLineColor(kBlue);
	hNNLOaRatio[2]->SetLineStyle(3);
	hNNLOaRatio[2]->Draw("same");










	gDataRatio->SetMarkerStyle(20);
	gDataRatio->SetMarkerSize(pointSize);
	gDataRatio->Draw("p same");

	gDataAllRatio->Draw("pz same");

	gPad->Update();
	gPad->RedrawAxis();
	TFrame *frD = (TFrame *) gPad->FindObject("TFrame");
	frD->SetFillStyle(0);
	frD->Draw();
	gPad->Update();
}


void  Histogram::CalculateGluonFraction()
{
    TString tag = "newver";

    TString currFileNNLO, currFileNLO, currFileLO;
    for(TString f: theor_files) {
        if(f.Contains(tag) && f.BeginsWith("nnlo:") ) currFileNNLO = f;
        if(f.Contains(tag) && f.BeginsWith("nlo:")  ) currFileNLO = f;
    }

    currFileLO = currFileNLO;
    currFileLO.ReplaceAll("nlo:", "lo:");

    auto GetTheories = [&](bool onlyQ) -> vector<Theory> {
        Setting setting;
        //setting.CreateFitB(0, "Q2pPt2", 1, 1, 1.0); 
        setting.CreateZeusSJ(0, "ZEUSscale", 1, 1, 1.0); 
        setting.onlyQuarks = onlyQ;

        Theory NNLO = CalcNLO(currFileNNLO, setting);
        Theory NLO  = CalcNLO(currFileNLO, setting);
        Theory LO   = CalcNLO(currFileNLO, setting);
        return {LO, NLO, NNLO};
    };

    auto thFull = GetTheories(false);

    for(unsigned i = 0; i < thFull[1].xsc.size(); ++i) 
        cout <<i <<" " << thFull[1].xsc[i] << endl;;

    /*
    auto thQ    = GetTheories(true);

   

    cout << "Fraction output " << endl;
    for(int order = 0; order < 3; ++order) {
        cout << "Order is " << order << endl;
        for(unsigned i = 0; i < thFull[order].xsc.size(); ++i) 
            cout <<i <<" " << thFull[order].xsc[i] <<" "<<  thQ[order].xsc[i] << " "<< thQ[order].xsc[i]/thFull[order].xsc[i]<< endl;;
    }
    */


}


//ren, fact, both
map<double, vector<vector<double>>>  Histogram::calculateScaleDependence(TString scaleTag, TString tag)
{

	//TString fileTag = isRen ? "renScaleDep"+tag+".txt" : "factScaleDep"+tag+".txt";
	TString fileTag = scaleTag+ "ScaleDep"+tag+".txt";

	TString fileName = TString("../data/Theory/") + theor_path + var_name +"/" + fileTag;//  theor_file; 
	cout << fileName << endl;


	map<double, vector<vector<double>>> myMap;

	//fileName = "myTest.txt";

	fstream file;
	file.open(fileName.Data(), fstream::in );

	if(file.good()) {
		int nBins;
		double scale;
		while(file >> nBins >> scale) {
			vector<vector<double>> temp;
			for(int i = 0; i < nBins; ++i) {
				vector<double> vals(9);
				for(int j = 0; j < 9; ++j)
					file >> vals[j];
				temp.push_back(vals);
			}
			myMap[scale] = temp;
		}
		file.close();
	}
	//Create file
	else {
		file.open( fileName.Data(), fstream::out );

		TString currFileNNLO, currFileNLO, currFileLO;
		for(TString f: theor_files) {
			if(f.Contains(tag) && f.BeginsWith("nnlo:") ) currFileNNLO = f;
			if(f.Contains(tag) && f.BeginsWith("nlo:")  ) currFileNLO = f;
		}

		currFileLO = currFileNLO;
		currFileLO.ReplaceAll("nlo:", "lo:");


		double minScal = 0.1, maxScal = 10;

		int nSteps = 16;
		double step = pow(maxScal/minScal, 1./nSteps);

		int scaleId = 0;
		for(double scale = minScal; scaleId <= nSteps; scale *= step, ++scaleId) {

			auto GetTheories = [&](double renScale, double factScale) -> vector<Theory> {
				Setting setting;
				setting.CreateFitB(0, "Q2pPt2", 1, 1, 1.0); 
				setting.factFactor = factScale;
				setting.renFactor  = renScale;

				Theory NNLO = CalcNLO(currFileNNLO, setting);
				Theory NLO  = CalcNLO(currFileNLO, setting);
				Theory LO   = CalcNLO(currFileLO, setting);
				return {LO, NLO, NNLO};
			};

			vector<Theory> th, thUp, thDown;

			if(scaleTag == "ren") {
				th     = GetTheories(scale, 1.0);
				thUp   = GetTheories(scale, 2.0);
				thDown = GetTheories(scale, 0.5);
			}
			else if(scaleTag == "fact") {
				th     = GetTheories(1.0,scale);
				thUp   = GetTheories(2.0,scale);
				thDown = GetTheories(0.5,scale);
			}
            else if(scaleTag == "both") {
				th     = GetTheories(scale,scale);
				thUp   = GetTheories(scale,2*scale);
				thDown = GetTheories(scale,0.5*scale);
            }
			vector<vector<double>> temp;

			file << th[0].xsc.size() <<" "<< scale << endl;
			for(unsigned i = 0; i < th[0].xsc.size(); ++i) {
				file << th[0].xsc[i] <<" "<<  thUp[0].xsc[i] << " "<< thDown[0].xsc[i] <<" ";
				file << th[1].xsc[i] <<" "<<  thUp[1].xsc[i] << " "<< thDown[1].xsc[i] <<" ";
				file << th[2].xsc[i] <<" "<<  thUp[2].xsc[i] << " "<< thDown[2].xsc[i] << endl;
				temp.push_back({th[0].xsc[i], thUp[0].xsc[i], thDown[0].xsc[i],
				                th[1].xsc[i], thUp[1].xsc[i], thDown[1].xsc[i],
				                th[2].xsc[i], thUp[2].xsc[i], thDown[2].xsc[i]});
			}
			myMap[scale] = temp;
		}
		file.close();

	}

	return myMap;

}

void plotScaleDependences(TString scaleTag, Histogram *h1, Histogram *h2, Histogram *h3, Histogram *h4, Histogram *h5)
{
	Histogram *histos[] = {h1, h2, h3, h4, h5};

	int hash = rand();
	TVirtualPad *can = gPad;//new TCanvas(TString::Format("ccc%d", hash), "Scale Dependence", 400, 600);

	int nBins = h1->getNbins();
	can->Divide(nBins, 5, 1e-5, 1e-5);

	for(int i = 0; i < 5; ++i) {
		for(int j = 0; j < nBins; ++j) {
			can->cd(3*i + j + 1);
			cout << "My bin " << j+1 << endl;
			histos[i]->plotScaleDependence(j+1, scaleTag);

			TLegend *leg = new TLegend( 0.2, 0.2, 0.6, 0.4 );
			leg->SetHeader(histos[i]->Title);
			leg->SetBorderSize(0);
			leg->SetFillStyle(0);
			leg->Draw();
		}
	}

	can->cd();
	/*
	if(isRen)
		can->SaveAs("MyNewHopeRen.pdf");
	else 
		can->SaveAs("MyNewHopeFac.pdf");
	*/

}

void Histogram::plotScaleDependence(int binId, TString scaleTag)
{
	TGraphAsymmErrors *grNNLO;
	TGraphAsymmErrors *grNLO;
	TGraphAsymmErrors *grLO;

	TGraphAsymmErrors *grDataStat;
	TGraphAsymmErrors *grDataTot;


	grNNLO = new TGraphAsymmErrors;
	grNLO = new TGraphAsymmErrors;
	grLO = new TGraphAsymmErrors;
	grDataStat = new TGraphAsymmErrors;
	grDataTot = new TGraphAsymmErrors;
	TH1D *hData;


    gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.16);

	map<double, vector<vector<double>>> myMap = calculateScaleDependence(scaleTag, defFile);

	int nBins = myMap.begin()->second.size();
	assert(binId <= nBins);

	{
		int i = binId-1;
		grDataStat->SetPoint(0, 1, data[i]);
		grDataStat->SetPointError(0, 0.94, 10 , dataStatErr[i], dataStatErr[i]);

		grDataTot->SetPoint(0, 1, data[i]);
		grDataTot->SetPointError(0, 0.94, 10, hypot(dataStatErr[i],dataSystErr[i]), hypot(dataStatErr[i],dataSystErr[i]) );

		


	}

	//grNNLO->SetFillColorAlpha(kBlue-7, 0.2);
	//grNNLO->SetFillStyle(1001);
	//grNLO->SetFillColorAlpha(kBlue-7, 0.2);
	//grNLO->SetFillStyle(1001);
	grLO->SetFillColorAlpha(kRed-7, 0.2);
	grLO->SetFillStyle(3144);
	grLO->SetLineColor(kBlack);
	//grNLO->SetLineColor(kBlack);
	//grNNLO->SetLineColor(kBlack);


	int green = kGreen - 3;
	//grNLO->SetLineColor(green);
	grNLO->SetFillColorAlpha(green, 0.2);
	grNLO->SetFillStyle(3144);


	double allTr = 0.35;
	double inTr = 0.45;
	grNNLO->SetLineWidth(0);
	grNNLO->SetLineColorAlpha(   kBlue-7, allTr /* kAzure-4*/);
	grNNLO->SetFillColorAlpha(kBlue-7, allTr);
	grNNLO->SetFillStyle(1001);






	//Set scale graphs
	int idScale = 0;
	for(const auto &el : myMap) {
		double scale = el.first;
		vector<vector<double>> vec = el.second;
		assert( vec.size() == nBins);
		{
			int i = binId-1;
			grLO->SetPoint(idScale, scale, vec[i][0]);
			grNLO->SetPoint(idScale, scale, vec[i][3]);
			grNNLO->SetPoint(idScale, scale, vec[i][6]);

			double errUp, errDown;
			errUp   = max({0., vec[i][1]-vec[i][0], vec[i][2]-vec[i][0]});
			errDown =-min({0., vec[i][1]-vec[i][0], vec[i][2]-vec[i][0]});
			grLO->SetPointError(idScale,0,0, errDown, errUp);

			errUp   = max({0., vec[i][4]-vec[i][3], vec[i][5]-vec[i][3]});
			errDown =-min({0., vec[i][4]-vec[i][3], vec[i][5]-vec[i][3]});
			grNLO->SetPointError(idScale,0,0, errDown, errUp);

			errUp   = max({0., vec[i][7]-vec[i][6], vec[i][8]-vec[i][6]});
			errDown =-min({0., vec[i][7]-vec[i][6], vec[i][8]-vec[i][6]});
			grNNLO->SetPointError(idScale,0,0, errDown, errUp);

		}
		++idScale;
	}


	//TCanvas *can = new TCanvas("ccc", "Scale Dependence");

	{
		int i = binId -1;
		gPad->SetLogx();
		gPad->SetLogy();

		double x, yMax, yMin;
		grNLO->GetPoint(0, x, yMax);
		grLO->GetPoint(grLO->GetN()-1, x, yMin);

		TH1 *hTemp;
		if(scaleTag == "ren") {
			//hTemp = gPad->DrawFrame(0.1, 0.7*yMin, 10, 1.3*yMax);
			hTemp = gPad->DrawFrame(0.1, 15.1489, 10, 326.664);
			cout << "HelenkaR " << 0.7*yMin << " "<< 1.3*yMax << endl;
		}
		else if(var_name == "total")
			hTemp = gPad->DrawFrame(0.2, 15.1489, 5, 326.664);
		else
			hTemp = gPad->DrawFrame(0.2, 0.7*yMin, 5, 1.3*yMax);

		//TH1 *hTemp = gPad->DrawFrame(0.1, 0.7*yMin, 10, 1.3*yMax);

		if(scaleTag == "ren") 
			hTemp->GetXaxis()->SetTitle("#mu_{R}/#sqrt{Q^{2}+p_{T}^{2}}");
		else if(scaleTag == "fact")
			hTemp->GetXaxis()->SetTitle("#mu_{F}/#sqrt{Q^{2}+p_{T}^{2}}");
        else
			hTemp->GetXaxis()->SetTitle("#mu_{R,F}/#sqrt{Q^{2}+p_{T}^{2}}");

		hTemp->GetXaxis()->SetTitleOffset(1.4);
		hTemp->GetYaxis()->SetTitleOffset(1.5);

		hTemp->GetXaxis()->SetMoreLogLabels();
		hTemp->GetYaxis()->SetMoreLogLabels();
		hTemp->GetXaxis()->SetNoExponent();
		hTemp->GetYaxis()->SetNoExponent();


		hTemp->GetYaxis()->SetTitle( yTitle );
		cout << "RADECEK bin " << i << " " << xMin[i]<<" "<< xMin[i] << endl;
		if(xTitle != "dummy")
			hTemp->SetTitle(SF("%g < %s < %g", xMin[i], xTitle.Data(), xMax[i]) );
		else
			hTemp->SetTitle("");

		//grNNLO->SetMaximum(1.3*yMax);
		//grNNLO->SetMinimum(0.7*yMin);

		grNNLO->SetLineStyle(1);
		//grNLO->SetLineStyle(2);
		//grLO->SetLineStyle(3);

		grNNLO->SetLineColor(kBlue+2);
		grNNLO->SetLineStyle(1);
		grNNLO->SetLineWidth(2);

		grNLO->SetLineColor(kTeal+3);
		grNLO->SetLineWidth(2);

		grLO->SetLineColor(kRed);
		grLO->SetLineWidth(2);


		//grNNLO->GetXaxis()->SetRangeUser(1./4, 4);
		//grNNLO->Draw("a c3");
		grNNLO->Draw("c3 same");

		grNLO->Draw("c3 same");
		grLO->Draw("c3 same");


		grDataTot->SetFillStyle(3004);
		grDataTot->SetLineColor(kWhite);
		grDataTot->SetFillColor(kBlack);
		grDataTot->Draw("e2 same");
		grDataStat->SetFillStyle(3004);
		grDataStat->SetLineColor(kBlack);
		grDataStat->SetFillColor(kBlack);
		grDataStat->SetLineWidth(2);
		grDataStat->Draw("Le2 same");


		hData = new TH1D(TString::Format("hData%d", rand()),"hist", 1, GetXaxis()->GetXmin(), GetXaxis()->GetXmax());
		hData->SetBinContent(1, data[0]);


		hData->SetLineColor(kBlack);
		hData->SetLineWidth(2);
		hData->Draw("same ][");

		//grNNLO->Draw("c4 same");
		//grNLO->Draw("c4 same");
		//grLO->Draw("c4 same");

        SetFonts(17);


		//grNNLO->SetLineStyle(1);
		//grNLO->SetLineStyle(2);
		//grLO->SetLineStyle(3);
        TLegend *leg;
        //if(isRen)
        //leg = new TLegend( 0.5, 0.65, 0.9, 0.85 );
        //else


        leg = new TLegend( 0.52-0.05, 0.67-0.01, 0.92-0.05, 0.87+0.02 );
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->SetTextSize( hTemp->GetXaxis()->GetTitleSize() /*0.03427*/);

		//char lab = (scaleTag != "ren")  ? 'R' : 'F';
		char lab = (scaleTag == "fact")  ? 'R' : 'F';
		leg->AddEntry(grNNLO,  SF("NNLO + (#mu_{%c} unc.)", lab), "lf");
		leg->AddEntry(grNLO,   SF("NLO + (#mu_{%c} unc.)", lab), "lf");
		leg->AddEntry(grLO,    SF("LO + (#mu_{%c} unc.)", lab), "lf");
		leg->AddEntry(grDataStat, "H1 LRG Data" , "lf");
		leg->Draw();
        DrawLatexUp(1.1, "#scale[1.2]{#font[72]{NNLOJET}}", -1, "l"); //RADEK
	}

	//can->SaveAs("ScaleDependenceRadecek.pdf");

}

/*
void Histogram::plotScaleDependence()
{

	TString currFileNNLO, currFileNLO, currFileLO;
	for(TString file : theor_files) {
		if(file.Contains("H1-LQall-8c") && file.BeginsWith("nnlo:") ) currFileNNLO = file;
		if(file.Contains("H1-LQall-8c") && file.BeginsWith("nlo:")  ) currFileNLO = file;
	}

	currFileLO = currFileNLO;
	currFileLO.ReplaceAll("nlo:", "lo:");


	TGraph *grNNLO[20];
	TGraph *grNLO[20];
	TGraph *grLO[20];
	TGraph *grRat[20];
	TGraphAsymmErrors *grDataStat[20];
	TGraphAsymmErrors *grDataTot[20];
	for(int i = 0; i < 20; ++i) {
		grNNLO[i] = new TGraph;
		grNLO[i] = new TGraph;
		grLO[i] = new TGraph;
		grRat[i] = new TGraph;
		grDataStat[i] = new TGraphAsymmErrors;
		grDataTot[i] = new TGraphAsymmErrors;
	}


	int scaleId = 0;
	for(double scale = 1./pow(1.4,7); scale < pow(1.4,8); scale *= 1.4) {

		Setting setting;
		setting.CreateFitB(0, "Q2pPt2", 1, 1, 1.0); //Scale up
		setting.factFactor = 1.0;
		setting.renFactor  = scale;

		Theory NNLO = CalcNLO(currFileNNLO, setting);
		Theory NLO  = CalcNLO(currFileNLO, setting);
		Theory LO   = CalcNLO(currFileLO, setting);
		cout << scale <<" ";
		for(unsigned i = 0; i < NNLO.xsc.size(); ++i) {
			grNNLO[i]->SetPoint(scaleId, scale, NNLO.xsc[i]);
			grNLO[i] ->SetPoint(scaleId, scale, NLO.xsc[i]);
			grLO[i] ->SetPoint(scaleId, scale, LO.xsc[i]);
			grRat[i]->SetPoint(scaleId, scale, NLO.xsc[i]/LO.xsc[i]);


			//cout << NNLO.xsc[i] <<" ";
		}
		//cout << endl;
		++scaleId;
	}

	for(int i = 0; i < 4; ++i) {
		grDataStat[i]->SetPoint(0, 1, data[i]);
		grDataStat[i]->SetPointError(0, 0.94, 10 , dataStatErr[i], dataStatErr[i]);

		grDataTot[i]->SetPoint(0, 1, data[i]);
		grDataTot[i]->SetPointError(0, 0.94, 10, hypot(dataStatErr[i],dataSystErr[i]), hypot(dataStatErr[i],dataSystErr[i]) );
	}



	#define SF TString::Format 

	TCanvas *can = new TCanvas("ccc", "Scale Dependence");
	can->Divide(2,2);
	for(int i = 0; i < 4 ; ++i) {
		can->cd(i+1);
		gPad->SetLogx();
		gPad->SetLogy();

		double x, yMax, yMin;
		grNLO[i]->GetPoint(0, x, yMax);
		grLO[i]->GetPoint(grLO[i]->GetN()-1, x, yMin);

		grNNLO[i]->GetXaxis()->SetTitle("#mu_{r}/#sqrt{Q^{2}+p_{T}^{2}}");
		grNNLO[i]->GetYaxis()->SetTitle( yTitle );
		grNNLO[i]->SetTitle(SF("%g < %s < %g", xMin[i], xTitle.Data(), xMax[i]) );

		grNNLO[i]->SetMaximum(1.2*yMax);
		grNNLO[i]->SetMinimum(0.8*yMin);

		grNNLO[i]->SetLineStyle(1);
		grNLO[i]->SetLineStyle(2);
		grLO[i]->SetLineStyle(3);

		grNNLO[i]->Draw("acp");
		grNLO[i]->Draw("cp same");
		grLO[i]->Draw("cp same");


		grDataTot[i]->SetFillStyle(3004);
		grDataTot[i]->SetLineColor(kWhite);
		grDataTot[i]->SetFillColor(kBlue);
		grDataTot[i]->Draw("e2 same");
		grDataStat[i]->SetFillStyle(3004);
		grDataStat[i]->SetLineColor(kRed);
		grDataStat[i]->SetFillColor(kRed);
		grDataStat[i]->Draw("Le2 same");

		grNNLO[i]->Draw("cp same");
		grNLO[i]->Draw("cp same");
		grLO[i]->Draw("cp same");

		grNNLO[i]->SetLineStyle(1);
		grNLO[i]->SetLineStyle(2);
		grLO[i]->SetLineStyle(3);

		TLegend *leg = new TLegend( 0.5, 0.7, 0.9, 0.9 );
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);

		leg->AddEntry(grNNLO[i],  "NNLO", "l");
		leg->AddEntry(grNLO[i],   "NLO" , "l");
		leg->AddEntry(grLO[i],    "LO" , "l");
		leg->AddEntry(grDataStat[i], "Data" , "lf");
		leg->Draw();
	}


	for(int k = 0; k < 3; ++k) {
		cout << "Cross section table for ptRange "<< k << endl;
		for(int i = 0; i < grLO[k]->GetN(); ++i) {
			double x, yLO, yNLO, yNNLO;
			grLO[k]->GetPoint(i, x, yLO);
			grNLO[k]->GetPoint(i, x, yNLO);
			grNNLO[k]->GetPoint(i, x, yNNLO);
			cout << x << " "<< 6*x<<" "<<  yLO<<" "<< yNLO<<" "<< yNNLO << endl;
		}
	}



	can->SaveAs("ScaleDependenceLRGxpom.eps");

	TCanvas *dan = new TCanvas("ddd", "Scale Dependence - Ratio");
	dan->Divide(2,2);
	for(int i = 0; i < 4 ; ++i) {
		dan->cd(i+1);
		gPad->SetLogx();
		grRat[i]->SetMinimum(0);
		grRat[i]->SetMaximum(3.5);
		grRat[i]->SetTitle(SF("%g < %s < %g", xMin[i], xTitle.Data(), xMax[i]));

		grRat[i]->GetXaxis()->SetTitle("#mu_{r}/(Q^{2}+p_{T}^{2})");
		grRat[i]->GetYaxis()->SetTitle("#sigma_{NLO}/#sigma_{LO}");

		grRat[i]->Draw("acp");
	}
	dan->SaveAs("ScaleDependenceRatioFactLRGxpom.eps");

}
*/









void readHistograms(vector<Setting> setting, map<const char *, Histogram> &hist, TString fileName,
                   const char *n1,   const char *n2=0, const char *n3=0, const char *n4=0, const char *n5=0, const char *n6=0,
                   const char *n7=0, const char *n8=0, const char *n9=0, const char *n10=0)
{
	const char *names[] = {n1, n2, n3, n4, n5, n6, n7, n8, n9, n10};
	for(unsigned i = 0; i < sizeof(names)/sizeof(names[0]); ++i) {
		if(names[i] == 0)
			break;
		if(hist.count(names[i]) > 0) {
			cout << "Histogram " << names[i] << " already exists." << endl;
			exit(1);
		}
		hist[names[i]].loadData( fileName, names[i] );
		hist[names[i]].IncludeSettings(setting);
		hist[names[i]].LoadTheories();
		//hist[names[i]].LoadNLO();
		//hist[names[i]].LoadNNLO();
	}
}

void plotXi12(TString fileName, const char *xpomN, const char *zpomN)
{
	Histogram hXpom, hZpom;

	hXpom.loadData( fileName, xpomN );
	hZpom.loadData( fileName, zpomN );

	bool isLog = strcmp(xpomN, "xpom");  //if not equal to xpom

	TH1D *hXpomData = new TH1D("XpomData", "XpomData", hXpom.nBins, hXpom.GetBinning().data() );
	TH1D *hZpomData = new TH1D("ZpomData", "ZpomData", hZpom.nBins, hZpom.GetBinning().data() );

	auto getDataHist = [](TH1D *h, vector<double> &data) {
		for(unsigned i = 0; i < data.size(); ++i)
			h->SetBinContent(i+1, data[i]);
	};
	getDataHist(hXpomData, hXpom.data);
	getDataHist(hZpomData, hZpom.data);

	double xiMax = isLog ? pow(10,hXpom.GetBinning().back() ) : hXpom.GetBinning().back();
	TH1D *hXi12 = new TH1D("Xi12", "Xi12", 20, 0, xiMax );
	for(unsigned i = 0; i < 100000; ++i) {
		double xpom, zpom, xi12;
		xpom = hXpomData->GetRandom();
		zpom = hZpomData->GetRandom();
		
		if( isLog ) //if not equal to xpom
			xpom = pow(10,xpom);

		//cout <<"RADEK "<<  xpom << " " << zpom << endl;

		xi12 = xpom*zpom;
		hXi12->Fill(xi12);
	}
	hXi12->SetTitle( fileName );
	hXi12->Draw();
	//hXi12->SaveAs("test.root");

}


void PlotFour(TCanvas *can, TString Type, map<const char *, Histogram> hists, const char *n1, const char *n2=0, const char *n3=0, const char *n4=0 )
{
	can->Divide(2,2,0.00001,0.00001); 

	const char *names[] = {n1, n2, n3, n4};

	for(int i = 0; i < 4 && names[i] != 0; ++i) {
		can->cd(i+1);
		try {
			if(Type == "NLOvsNNLO")
				hists.at(names[i]).plotNLOvsNNLO();
			else if(Type == "Scales")
				hists.at(names[i]).plotScaleStudies();
			else if(Type == "DPDF")
				hists.at(names[i]).plotDPDFStudies();
			else {
				cout << "Unknown plot type"<<endl;
				exit(1);
			}
		}
		catch(const std::out_of_range& oor) {
			cout << "Variable " << names[i]<<" not known" << endl;
			exit(1);
		}
	}
	can->SaveAs( can->GetTitle() );
	can->Clear();
}

double getChi2(TH1D *hData, TGraphAsymmErrors *grAll, TH1D *hModel)
{
	/*
	double dataInt = hData->Integral("width");
	double modelInt = hModel->Integral("width");

	double sum = 0;
	for(int i = 1; i <= hModel->GetNbinsX(); ++i) {
		double err = grAll->GetErrorY(i-1);
		sum += pow(hModel->GetBinContent(i)/modelInt*dataInt - hData->GetBinContent(i),2) / err/err;
	}
	return sum/(hData->GetNbinsX()-1);
	*/

	double sumUp = 0, sumDn = 0;
	for(int i = 1; i <= hModel->GetNbinsX(); ++i) {
		double relEr = grAll->GetErrorY(i-1)/hData->GetBinContent(i);
		sumUp += log(hData->GetBinContent(i)/hModel->GetBinContent(i)) / pow(relEr,2);
		sumDn += 1./pow(relEr,2);
	}
	double norm = exp(sumUp/sumDn);

	double sum = 0;
	for(int i = 1; i <= hModel->GetNbinsX(); ++i) {
		double relEr = grAll->GetErrorY(i-1)/hData->GetBinContent(i);
		sum += pow(log(hData->GetBinContent(i)/hModel->GetBinContent(i)/norm),2) / pow(relEr,2);
	}
	return sum/(hData->GetNbinsX()-1);

}


void PlotComparison(TCanvas *can, TString plotStyle,  vector<Histogram*> histos)
{

	//can->SetCanvasSize(500, 500);
	//can->Update();

	can->Divide(histos.size()+1, 1, 0.00001,0.00001); 


	double yMin, yMax;

	yMin = 0.3;
	if(plotStyle == "NLOvsNNLO") {
		if(histos[1]->var_name == "meaneta")
			yMax = 5.;
		else if(histos[1]->var_name == "y")
			yMax = 3.5;
		else if(histos[0]->var_name == "q2")
			yMax = 3.2;
		else if(histos[0]->var_name == "zpom")
			yMax = 4.;
		else if(histos[0]->var_name == "logxpom")
			yMax = 3.5;
		else if(histos[0]->var_name == "deltaetaStar")
			yMax = 3.6;
		else
			yMax = 3.0;
	}
	else if(plotStyle == "DPDF") {
		if(histos[0]->var_name == "zpom")
			yMax = 5.;
		else
			yMax = 4.;

	}
	else if(plotStyle == "Scale") {
		if(histos[0]->var_name == "zpom")
			yMax = 5.;
		else
			yMax = 3.;
	}

	double tickSize = 0.08;

	for(unsigned i = 0; i < histos.size(); ++i) {
		//if(i < 3)
			can->cd(2+i);
		//else
			//can->cd(3+i);
		if(histos[i])
			histos[i]->plotNLOvsNNLOratio(plotStyle, yMin, yMax);
		else { //plot empty
			TH1D * h = new TH1D( SF("hist%d",rand()), "hist", 10, -5, 5);
			gPad->SetLeftMargin(0);
			gPad->SetRightMargin(0);

			gPad->SetBottomMargin(0.25);
			gPad->SetTopMargin(0.06);

			gPad->SetTicks(1,1);
			h->GetXaxis()->SetLabelSize(0);
			h->GetYaxis()->SetNdivisions(505);
			h->GetYaxis()->SetTickLength(tickSize/1.044);
			h->GetXaxis()->SetTickLength(0);
			h->SetMinimum(yMin);
			h->SetMaximum(yMax);
			//h->GetYaxis()->SetTickLength(0);
			h->Draw("axis");	
		}
	}
	can->cd(1);
	gPad->SetBottomMargin(0.25);
	gPad->SetTopMargin(0.06);
//gPad->GetUymin()
	TGaxis *axis = new TGaxis(1.0, 0.25,
                            1.0, 0.94,
                            yMin,yMax,510,"-R");
    axis->SetLabelSize(0.1);
    axis->SetLabelFont(42);
    axis->SetLabelOffset(0.1);
    axis->SetNdivisions(204);
    axis->SetTickSize(0.4);
    axis->SetTitleSize(0.1);
    axis->SetTitleFont(42);
    axis->SetTitleOffset(1.8);
    axis->CenterTitle();
    axis->SetTitle("#sigma/#sigma_{NLO}");

    //cout <<"Helenka "<<  axis->GetLabelFont() << endl;
    axis->Draw();

	//map<TString, TString> analMap;
	//analMap["FPS"] = "#splitline{H1 FPS}{(HERA #Iota#Iota)}";
	//analMap["VFPS"] = "#splitline{H1 VFPS}{(HERA #Iota#Iota)}";
	//analMap["LRG"] = "#splitline{H1 LRG}{(HERA #Iota#Iota)}";
	//analMap["LRG_H1"] = "#splitline{H1 LRG}{(HERA #Iota)}";
	//analMap["LRGH1820"] = "#splitline{H1 LRG}{(300 GeV)}";
	//analMap["LRGZEUS"] = "#splitline{ZEUS LRG}{(HERA #Iota)}";
	
	map<TString, bool> present;
	for(unsigned i = 0; i < histos.size(); ++i)
		if(histos[i]) present[histos[i]->getTag()] = true;
	vector<TString> missing;
	for(const auto &tag : analMap)
		if(!present.count(tag.first))
			missing.push_back(tag.first);

    TString arrQ[]={
	"4 < Q^{2} < 6 GeV^{2}", "6 < Q^{2} < 10 GeV^{2}", "10 < Q^{2} < 18 GeV^{2}",
    "18 < Q^{2} < 34 GeV^{2}", "34 < Q^{2} < 100 GeV^{2}",""};



	int m = 0;
	for(unsigned i = 0; i < histos.size(); ++i) {
		can->cd(i+2);
		TString tag = histos[i] ? histos[i]->getTag() : missing[m++];
		TLegend *leg = new TLegend( 0.3, 0.60, 0.7,0.88);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.1);
		leg->SetTextAlign(23);
		leg->SetHeader( analMap.at(tag));
		//leg->SetHeader(arrQ[i]);
		leg->Draw();
	}


	HistoErr grData, grNlo, grNnlo;
	vector<double> dummLow = {0.}, dummHi= {1.};
	grNlo.Init("NloName", dummLow, dummHi);
	grNnlo.Init("NnloName", dummLow, dummHi);
	grData.Init("DataTemp", dummLow, dummHi);
	myCanvas::SetNLOstyle(grNlo);
	myCanvas::SetNNLOstyle(grNnlo);
	grData.grRatio->SetMarkerStyle(20);
	grData.grRatio->SetMarkerSize(0.8);


	can->cd(4);
	TLegend *leg1 = new TLegend( 0.1, 0.63, 0.6, 0.78 );
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	leg1->SetTextSize(0.1);

	leg1->AddEntry(grData.grRatio,    "Data", "ep");
	if(plotStyle != "Scale")
		leg1->AddEntry(grNlo.hRatio,     "NLO (H1 Fit B)", "l");
	else
		leg1->AddEntry(grNlo.hRatio,    "NLO (#mu^{2}=Q^{2}+p_{T}^{2})", "l");

	if(plotStyle == "NLOvsNNLO") 
		leg1->AddEntry(grNlo.grAllRatio, "scale unc.", "f");
	else
		leg1->AddEntry((TObject*)0, "", "");


	leg1->Draw();

	can->cd(5);
	TLegend *leg2 = new TLegend( 0.1, 0.63, 0.6, 0.78 );
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	leg2->SetTextSize(0.1);


	TGraphAsymmErrors *grTemp = new TGraphAsymmErrors;
	grTemp->SetLineWidth(0);
	grTemp->SetLineColorAlpha(kBlue-7, 0.0);
	grTemp->SetFillColorAlpha(kBlue-7, 0.6425); //0.7 org

	if(plotStyle != "Scale")
		leg2->AddEntry(grNnlo.hRatio,    "NNLO (H1 Fit B)", "l");
	else
		leg2->AddEntry(grNnlo.hRatio,    "NNLO (#mu^{2}=Q^{2}+p_{T}^{2})", "l");

	if(plotStyle != "DPDF") {
		leg2->AddEntry(grTemp /*grNnlo.grRatio*/,   "scale unc.", "f");
		leg2->AddEntry(grNnlo.grAllRatio,"scale+DPDF unc.", "f");
	}
	else {
		leg2->AddEntry(grTemp /*grNnlo.grRatio*/,   "DPDF unc.", "f");
		leg2->AddEntry(grNnlo.grAllRatio,"DPDF+scale unc.", "f");
	}
	leg2->Draw();

	if(plotStyle == "Scale") {
		TLegend *leg3 = new TLegend( 0.10, 0.63-0.16, 0.60, 0.78-0.16 );
		leg3->SetBorderSize(0);
		leg3->SetFillStyle(0);
		leg3->SetTextSize(0.1);

		TH1D *hDummyFitQ2 = new TH1D("hDummyQ2", "test",  1, -1, 1);
		TH1D *hDummyFitPt2 = new TH1D("hDummyPt2", "test",  1, -1, 1);
		TH1D *hDummyFitQ2Pt2 = new TH1D("hDummyQ2Pt2", "test",  1, -1, 1);

		hDummyFitQ2->SetLineColor(kRed); hDummyFitQ2->SetLineWidth(2);
		hDummyFitPt2->SetLineColor(kBlack); hDummyFitPt2->SetLineWidth(2);
		hDummyFitQ2Pt2->SetLineColor(kCyan); hDummyFitQ2Pt2->SetLineWidth(2);


		leg3->AddEntry(hDummyFitQ2Pt2,  "NNLO (#mu^{2}=#frac{Q^{2}}{4}+p_{T}^{2})", "l");
		leg3->AddEntry(hDummyFitPt2,    "NNLO (#mu^{2}=p_{T}^{2})", "l");
		leg3->AddEntry(hDummyFitQ2,  "NNLO (#mu^{2}=Q^{2})", "l");


		leg3->Draw();

	}

	if(plotStyle == "DPDF") {
		can->cd(6);
		TLegend *leg3 = new TLegend( 0.05, 0.63-0.05, 0.55, 0.78 );
		leg3->SetBorderSize(0);
		leg3->SetFillStyle(0);
		leg3->SetTextSize(0.1);

		TH1D* hDummyFitA {new TH1D(SF("hDummyFitA%d",rand()), "test",  1, -1, 1)};
		TH1D* hDummyFitJets { new TH1D(SF("hDummyFitA%d",rand()), "test",  1, -1, 1)};
		TH1D* hDummyFitSJ { new TH1D(SF("hDummyFitSJ%d",rand()), "test",  1, -1, 1)};
		TH1D* hDummyMRW { new TH1D(SF("hDummyMRW%d",rand()), "test",  1, -1, 1)};

		hDummyFitSJ->SetLineColor(kRed); hDummyFitSJ->SetLineWidth(2);
		hDummyFitJets->SetLineColor(kMagenta); hDummyFitJets->SetLineWidth(2);
		hDummyFitA->SetLineColor(kBlack); hDummyFitA->SetLineWidth(2);
		hDummyMRW->SetLineColor(kOrange); hDummyMRW->SetLineWidth(2);

		leg3->AddEntry(hDummyFitA,  "NNLO (H1 Fit A)", "l");
		leg3->AddEntry(hDummyFitJets,  "NNLO (H1 Fit Jets)", "l");
		leg3->AddEntry(hDummyFitSJ,  "NNLO (ZEUS SJ)", "l");
		leg3->AddEntry(hDummyMRW,  "NNLO (MRW)", "l");

		leg3->Draw();


	}

	can->SaveAs( can->GetTitle() );
	can->Clear();

}



void PlotSingle(TCanvas *can, TString plotStyle,  Histogram *hist)
{
    can->Clear();
    can->cd();


    double lMar = 0.15, rMag = 0.1;
    double padSize = (1 - lMar - rMag) / 6.;
    SetLeftRight(lMar, rMag + 5* padSize);

    DividePad({1}, {1,1});

    double yMin, yMax;
    double yMinAbs, yMaxAbs;


    if(hist->var_name == "beta") {
        yMin    = 0.49, yMax = 2.09;
        yMinAbs = 5.1, yMaxAbs = 1896000;
    }
    else if(hist->var_name == "xgamma") {
        yMin    = 0.0, yMax = 2.49;
        yMinAbs = 3.1, yMaxAbs = 45840;
    }

    can->cd(2);
    hist->plotNLOvsNNLOratioAbs(plotStyle, yMin, yMax, yMinAbs, yMaxAbs, 1);



    can->cd(1);

	HistoErr grData, grNlo, grNnlo;
	vector<double> dummLow = {0.}, dummHi= {1.};
	grNlo.Init("NloName", dummLow, dummHi);
	grNnlo.Init("NnloName", dummLow, dummHi);
	grData.Init("DataTemp", dummLow, dummHi);
	myCanvas::SetNLOstyle(grNlo);
	myCanvas::SetNNLOstyle(grNnlo);
	grData.grRatio->SetMarkerStyle(20);
	grData.grRatio->SetMarkerSize(0.8);



    double sh = 0.06;
    double h = 2*0.18;
    double shY = 0.41;
    vector<double> legSizes = { 0.1+sh, shY, 0.18+sh, shY+h };
	TLegend *leg2 = new TLegend(legSizes[0], legSizes[1], legSizes[2], legSizes[3]);


	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	leg2->SetTextSize(GetXaxis()->GetTitleSize());


	TGraphAsymmErrors *grTemp = new TGraphAsymmErrors;
	grTemp->SetLineWidth(0);
	grTemp->SetLineColorAlpha(kBlue-7, 0.0);
	grTemp->SetFillColorAlpha(kBlue-7, 0.6425); //0.7 org


	leg2->AddEntry(grData.grRatio,    "Data", "ep");
	if(plotStyle != "Scale")
		leg2->AddEntry(grNlo.hRatio,     "NLO (H1 Fit B)", "l");
	else
		leg2->AddEntry(grNlo.hRatio,    "NLO (#mu^{2}=Q^{2}+p_{T}^{2})", "l");

	if(plotStyle == "NLOvsNNLO") 
		leg2->AddEntry(grNlo.grAllRatio, "scale unc.", "f");
	else
		leg2->AddEntry((TObject*)0, "", "");


	if(plotStyle != "Scale")
		leg2->AddEntry(grNnlo.hRatio,    "NNLO (H1 Fit B)", "l");
	else
		leg2->AddEntry(grNnlo.hRatio,    "NNLO (#mu^{2}=Q^{2}+p_{T}^{2})", "l");

	if(plotStyle != "DPDF") {
		leg2->AddEntry(grTemp /*grNnlo.grRatio*/,   "scale unc.", "f");
		leg2->AddEntry(grNnlo.grAllRatio,"scale+DPDF unc.", "f");
	}
	else {
		leg2->AddEntry(grTemp /*grNnlo.grRatio*/,   "DPDF unc.", "f");
		leg2->AddEntry(grNnlo.grAllRatio,"DPDF+scale unc.", "f");
	}
	leg2->Draw();

/*
    can->cd();
    TLine *line = new TLine;
    line->SetLineColor(kBlack);
    line->DrawLineNDC(lMar, 0, lMar, 1);
    line->DrawLineNDC(1-rMag-5*padSize, 0, 1-rMag-5*padSize, 1);
    gPad->Update();
*/

	can->SaveAs( can->GetTitle() );
    can->Clear();
}



void PlotComparisonWithAbs(TCanvas *can, TString plotStyle,  vector<Histogram*> histos)
{
	can->Clear();
    can->cd();

	//can->SetCanvasSize(500, 500);
	//can->Update();

	//can->Divide(histos.size()+1, 1, 0.00001,0.00001); 

    double lMar = 0.15, rMag = 0.1;
    double padSize = (1 - lMar - rMag) / 6.;

    if(histos[0] && histos[0]->var_name == "ptjet1_q2_4_6") {
        SetLeftRight(lMar, rMag + padSize);
        DividePad({1,1,1,1,1}, {1,1});
    }
    else if(histos[0] && (histos[0]->var_name == "zpom_q2_4_10" ||  histos[0]->var_name == "zpom_q2_5_12") ) {
        SetLeftRight(lMar, rMag + 2* padSize);
        DividePad({1,1,1,1}, {1,1});
    }
    else if(histos[0] && histos[0]->var_name ==  "zpom_ptjet1_5_6p5") {
        SetLeftRight(lMar, rMag + 3* padSize);
        DividePad({1,1,1}, {1,1});
    }
    else if(histos[0] && (histos[0]->var_name ==  "mx" || histos[0]->var_name == "ptavg")) {
        SetLeftRight(lMar, rMag + 4* padSize);
        DividePad({1,1}, {1,1});
    }
    else if(histos[0] && (histos[0]->var_name ==  "beta" )) {
        cout << "AHHH " << lMar <<" "<< rMag + 5* padSize << endl;
        SetLeftRight(lMar, rMag + 5* padSize);
        cout << "AHHH " << gPad->GetLeftMargin() <<" "<< gPad->GetRightMargin() << endl;
        DividePad({1}, {1,1});
    }
    else {
        SetLeftRight(lMar, rMag);
        DividePad({1,1,1,1,1,1}, {1,1});
    }

	double yMin, yMax;
	double yMinAbs = 0.101, yMaxAbs = 5000;
    vector<double> Factors = {1,1,1,1,1,5};

	yMin = 0.3;
	if(plotStyle == "NLOvsNNLO" || plotStyle == "Scale" || plotStyle == "DPDF") {
		if(histos[1] &&  histos[1]->var_name == "meaneta") {
			yMax = 2.98;
			yMaxAbs = 130.;
            Factors = {0.3, 2.0, 1, 1, 1, 0.5};
        }
		else if(histos[0]->var_name == "y") {
			yMax = 2.49;
			yMaxAbs = 450;
            Factors = {0.3, 2., 1, 1, 200, 200};
        }
		else if(histos[0]->var_name == "q2") {
			yMax = 2.49;
			yMaxAbs = 465;
			yMinAbs = 0.03;
            if(plotStyle == "Scale")
                Factors = {0.3, 2, 1, 1, 1, 1};
            else
                Factors = {0.3, 2, 0.5, 1, 1, 1};
        }
		else if(histos[0]->var_name == "zpom") {
			yMax = 2.49;
			yMaxAbs = 350;
            Factors = {0.3, 2., 1, 1, 1, 1};
        }
		else if(histos[0]->var_name == "logxpom") {
			yMax = 2.49;
            yMaxAbs = 550;
            Factors = {0.4, 0.1, 1, 1, 2, 0.07};
        }
		else if(histos[0]->var_name == "deltaetaStar") {
			yMax = 2.49;
            yMaxAbs = 170;
            Factors = {0.4, 2, 1, 1, 1, 1};
        }
        else if(histos[0]->var_name == "ptjet1_q2_4_6") {
			yMax = 2.29;
            yMaxAbs = 90;
            yMinAbs = 0.04;
            Factors = {1, 2, 6, 6, 24, 1};
        }
        else if(histos[0]->var_name == "zpom_q2_4_10") {
			yMax = 2.29;
            yMaxAbs = 22;
            Factors = {1, 3, 6, 24, 1, 1};
        }
        else if(histos[0]->var_name == "zpom_ptjet1_5_6p5") {
			yMax = 2.89;
            yMaxAbs = 150;
            Factors = {1, 2, 24, 4, 1, 1};
        }
        else if((plotStyle == "DPDF" || plotStyle == "Scale") && histos[0]->var_name == "ptjet1") {
            yMaxAbs = 1e4;
        }
        else if(histos[0]->var_name == "zpom_q2_5_12") {
			yMax = 2.69;
            yMaxAbs = 20;
            Factors = {1, 3, 6, 24, 1, 1};
        }
        else if(histos[0]->var_name == "mx") {
			yMax = 2.49;
			yMin = 0.0;
            yMaxAbs = 15;
            Factors = {2, 1, 6, 24, 1, 1};
        }
        else if(histos[0]->var_name == "ptavg") {
			yMax = 2.49;
			yMin = 0.0;
        }
		else
			yMax = 2.499;
	}
	else if(plotStyle == "DPDF") {
		if(histos[0]->var_name == "zpom")
			yMax = 5.;
		else
			yMax = 4.;

	}
	else if(plotStyle == "Scale") {
		if(histos[0]->var_name == "zpom")
			yMax = 5.;
		else
			yMax = 2.5;
	}


	double tickSize = 0.08;

	for(unsigned i = 0; i < histos.size(); ++i) {
		//if(i < 3)
			can->cd(1 + GetNpads(can)/2 + i);
            //if (i >= 1) continue;
		//else
			//can->cd(3+i);
		if(histos[i]) {

            if(histos[0] && histos[0]->var_name == "y") {
                histos[0]->yTitle = "d#sigma/dy [pb], d#sigma/dW [pb/GeV]";
                //histos[0]->yTitle = "d#sigma/dy [pb] or d#sigma/dW [pb/GeV]";
                //histos[0]->yTitle = "d#sigma/d{y, W} [pb, pb/GeV]";
                //histos[0]->yTitle = "d#sigma/dy [pb] {d#sigma/dW [pb/GeV]}";
            }

			histos[i]->plotNLOvsNNLOratioAbs(plotStyle, yMin, yMax, yMinAbs, yMaxAbs, Factors[i]);

        }
		else { //plot empty
            if(histos[0] && (histos[0]->var_name == "ptjet1_q2_4_6"  ||
                             histos[0]->var_name == "zpom_q2_4_10"  ||
                             histos[0]->var_name == "zpom_ptjet1_5_6p5"  ||
                             histos[0]->var_name == "zpom_q2_5_12" ||
                             histos[0]->var_name == "mx" ||
                             histos[0]->var_name == "ptavg" ||
                             histos[0]->var_name == "beta" 
                             )
            ) continue;

			TH1D * hDown = new TH1D( SF("hist%d",rand()), "hist", 10, -5, 5);
			//gPad->SetLeftMargin(0);
			//gPad->SetRightMargin(0);

			//gPad->SetBottomMargin(0.25);
			//gPad->SetTopMargin(0.06);
			hDown->Draw("axis");	

            SetFTO({12}, {5}, {1.5, 2.6, 0.3, 3.7});
            GetXaxis()->SetNdivisions(505);
            GetYaxis()->SetNdivisions(505);
            GetYaxis()->CenterTitle();
			gPad->SetTicks(1,1);
			GetXaxis()->SetLabelSize(0);
			GetXaxis()->SetTickLength(0);
			hDown->SetMinimum(yMin);
			hDown->SetMaximum(yMax);
			//h->GetYaxis()->SetTickLength(0);

            if(histos[1] && histos[1]->var_name == "meaneta")
                GetYaxis()->SetTitle("d#sigma/d#sigma_{NLO}");


			can->cd(1+i);

			TH1D * hUp = new TH1D( SF("hist%d",rand()), "hist", 10, -5, 5);
            if(histos[1]->var_name == "q2") gPad->SetLogy();

			hUp->Draw("axis");	
			gPad->SetTicks(1,1);
            SetFTO({12}, {5}, {1.5, 2.6, 0.3, 3.7});
            GetXaxis()->SetNdivisions(505);
            GetYaxis()->SetNdivisions(505);
            GetYaxis()->CenterTitle();
			gPad->SetTicks(1,1);
			GetXaxis()->SetLabelSize(0);
			GetXaxis()->SetTickLength(0);
			hUp->SetMinimum(yMinAbs);
			hUp->SetMaximum(yMaxAbs);
            //if(histos[1]->style.isLogY) gPad->SetLogy();

            if(histos[1] && histos[1]->var_name == "meaneta")
                GetYaxis()->SetTitle("d#sigma/d#LT#eta#GT [pb]");

		}
	}

    /*
    if(histos[0] && histos[0]->var_name == "y") {
        can->cd(5);
        DrawLatexUp (-2, "#spade", -1, "r");
        can->cd(6);
        DrawLatexUp (-2, "#spade", -1, "r");
    }
    */


	
	map<TString, bool> present;
	for(unsigned i = 0; i < histos.size(); ++i)
		if(histos[i]) present[histos[i]->getTag()] = true;
	vector<TString> missing;
	for(const auto &tag : analMap)
		if(!present.count(tag.first))
			missing.push_back(tag.first);

    //Calculate legend positions
    unsigned iLeg1 = 4;
    unsigned iLeg2 = 5;

    if(histos[0]) { 
       if( histos[0]->var_name == "zpom_q2_4_10")           iLeg1 = 3, iLeg2 = 4;
       else if( histos[0]->var_name == "zpom_ptjet1_5_6p5") iLeg1 = 2, iLeg2 = 3;
       else if( histos[0]->var_name == "zpom_q2_5_12" )     iLeg1 = 3, iLeg2 = 4;
       else if( histos[0]->var_name == "q2" )               iLeg1 = 3, iLeg2 = 4;
       else if( histos[0]->var_name == "mx" )               iLeg1 = 1, iLeg2 = 2;
       else if( histos[0]->var_name == "ptavg" )            iLeg1 = 1, iLeg2 = 2;
       else if( histos[0]->var_name == "beta" )             iLeg1 = 1, iLeg2 = 1;
       else if( histos[0]->var_name == "ptjet1" && plotStyle == "DPDF" )           iLeg1 = 3, iLeg2 = 4;
       else if( histos[0]->var_name == "ptjet1" && plotStyle == "Scale" )          iLeg1 = 3, iLeg2 = 4;
    }


	int m = 0;
	for(unsigned i = 0; i < histos.size(); ++i) {

        if ( int(i) >= GetNpads(can)/2) continue;

		can->cd(i+1);

        if(histos[0] &&  names2D.count(histos[0]->var_name) > 0 ) {
            //TString tag = histos[i] ? arrQ[i] : missing[m++];
            if(histos[i] &&  names2D.count(histos[i]->var_name) > 0) {
                TString n = names2D.at(histos[i]->var_name);
                DrawLatexUp (-2, n);
            }
        }
        //else if(histos[0] && histos[0]->var_name == "mx") {
            //if(i >= 2) continue;
            //DrawLatexUp (-2,  analMapInline.at(i ? "LRGZEUS" : "VFPS") );
        //}
        else {
            TString tag = histos[i] ? histos[i]->getTag() : missing[m++];
            DrawLatexUp (-2,  analMapInline.at(tag) );
        }


        if(Factors[i] != 1 && histos[i]) {
            DrawLatexUp (-3.3 - (i==iLeg1-1 || i==iLeg2-1)*4.2,  SF("x %g", Factors[i]) );
        }
	}


    //int corr = 0;
    //if(histos[1] && histos[1]->var_name == "meaneta")
        //corr = 1;

    //DrawLatexUp(can->GetPad(1),  -4, "#scale[1.3]{#font[72]{NNLOJET}}"); //RADEK
    //DrawLatexDown(can->GetPad(GetNpads(can)),  -2, "#scale[1.3]{#font[72]{NNLOJET}}"); //RADEK
    DrawLatexUp(can->GetPad(1),  1.1, "#scale[1.2]{#font[72]{NNLOJET}}", -1, "l"); //RADEK
    //DrawLatexUp(can->GetPad(GetNpads(can)/2),  1.3, "#scale[1.3]{#font[72]{NNLOJET}}", -1, "r"); //RADEK

	HistoErr grData, grNlo, grNnlo;
	vector<double> dummLow = {0.}, dummHi= {1.};
	grNlo.Init("NloName", dummLow, dummHi);
	grNnlo.Init("NnloName", dummLow, dummHi);
	grData.Init("DataTemp", dummLow, dummHi);
	myCanvas::SetNLOstyle(grNlo);
	myCanvas::SetNNLOstyle(grNnlo);
	grData.grRatio->SetMarkerStyle(20);
	grData.grRatio->SetMarkerSize(0.8);



    vector<double> legSizes = { 0.1, 0.50-0.04, 0.6, 0.68-0.04 };
	can->cd(iLeg1);
    double H = 1 - gPad->GetLeftMargin();
	TLegend *leg1 = new TLegend(1 - H*(1- legSizes[0]), legSizes[1], 1- H*(1-legSizes[2]), legSizes[3]);
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	leg1->SetTextSize(GetXaxis()->GetTitleSize());


    TString dataTitle = "Data";

    if(histos[0]) {
        if( histos[0]->var_name ==    "zpom_q2_4_10")        dataTitle = analMapInline.at("LRG");
        else if( histos[0]->var_name == "ptjet1_q2_4_6")     dataTitle = analMapInline.at("LRG");
        else if( histos[0]->var_name == "zpom_ptjet1_5_6p5") dataTitle = "ZEUS LRG";
        else if( histos[0]->var_name == "zpom_q2_5_12" )     dataTitle = "ZEUS LRG";
    }



	leg1->AddEntry(grData.grRatio,  dataTitle, "ep");
	if(plotStyle != "Scale")
		leg1->AddEntry(grNlo.hRatio,     "NLO (H1 Fit B)", "l");
	else
		leg1->AddEntry(grNlo.hRatio,    "NLO (#mu^{2}=Q^{2}+p_{T}^{2})", "l");

	if(plotStyle == "NLOvsNNLO") 
		leg1->AddEntry(grNlo.grAllRatio, "scale unc.", "f");
	else
		leg1->AddEntry((TObject*)0, "", "");


	leg1->Draw();

	can->cd(iLeg2);
    H = 1 - gPad->GetRightMargin();
	TLegend *leg2 = new TLegend( legSizes[0]*H, legSizes[1], legSizes[2]*H, legSizes[3]);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	leg2->SetTextSize(GetXaxis()->GetTitleSize());


	TGraphAsymmErrors *grTemp = new TGraphAsymmErrors;
	grTemp->SetLineWidth(0);
	grTemp->SetLineColorAlpha(kBlue-7, 0.0);
	grTemp->SetFillColorAlpha(kBlue-7, 0.6425); //0.7 org

	if(plotStyle != "Scale")
		leg2->AddEntry(grNnlo.hRatio,    "NNLO (H1 Fit B)", "l");
	else
		leg2->AddEntry(grNnlo.hRatio,    "NNLO (#mu^{2}=Q^{2}+p_{T}^{2})", "l");

	if(plotStyle != "DPDF") {
		leg2->AddEntry(grTemp /*grNnlo.grRatio*/,   "scale unc.", "f");
		leg2->AddEntry(grNnlo.grAllRatio,"scale+DPDF unc.", "f");
	}
	else {
		leg2->AddEntry(grTemp /*grNnlo.grRatio*/,   "DPDF unc.", "f");
		leg2->AddEntry(grNnlo.grAllRatio,"DPDF+scale unc.", "f");
	}
	leg2->Draw();

	if(plotStyle == "Scale") {
		//TLegend *leg3 = new TLegend( 0.10, 0.63-0.16, 0.60, 0.78-0.16 );
        double sh = 0.22;
        if( histos[0]->var_name == "ptjet1") {
            sh = 0;
            can->cd(5);
        }

        TLegend *leg3 = new TLegend( legSizes[0]*H, legSizes[1]-sh-0.10, legSizes[2]*H, legSizes[3] -sh);
		leg3->SetBorderSize(0);
		leg3->SetFillStyle(0);
		leg3->SetTextSize(GetXaxis()->GetTitleSize());

		TH1D *hDummyFitQ2 = new TH1D(SF("hDummyQ2%d",rand()), "test",  1, -1, 1);
		TH1D *hDummyFitPt2 = new TH1D(SF("hDummyPt2%d",rand()), "test",  1, -1, 1);
		TH1D *hDummyFitQ2Pt2 = new TH1D(SF("hDummyQ2Pt2%d",rand()), "test",  1, -1, 1);

		hDummyFitQ2->SetLineColor(kRed); hDummyFitQ2->SetLineWidth(2);
		hDummyFitPt2->SetLineColor(kBlack); hDummyFitPt2->SetLineWidth(2);
		hDummyFitQ2Pt2->SetLineColor(kCyan); hDummyFitQ2Pt2->SetLineWidth(2);


		leg3->AddEntry(hDummyFitQ2Pt2,  "NNLO (#mu^{2}=#frac{Q^{2}}{4}+p_{T}^{2})", "l");
		leg3->AddEntry(hDummyFitPt2,    "NNLO (#mu^{2}=p_{T}^{2})", "l");
		leg3->AddEntry(hDummyFitQ2,  "NNLO (#mu^{2}=Q^{2})", "l");


		leg3->Draw();

	}

	if(plotStyle == "DPDF") {
		can->cd(5);
		//TLegend *leg3 = new TLegend( 0.05, 0.63-0.05, 0.55, 0.78 );
        double sh = histos[0]->var_name != "ptjet1" ? 0.18 : 0;


        double d = (legSizes[3] - legSizes[1]) / 3. *4;
		TLegend *leg3 = new TLegend( legSizes[0]*H, legSizes[3]-sh-d, legSizes[2]*H, legSizes[3] -sh);
		leg3->SetBorderSize(0);
		leg3->SetFillStyle(0);
		leg3->SetTextSize(GetXaxis()->GetTitleSize());

		TH1D *hDummyFitA = new TH1D(SF("hDummyFitA%d", rand()), "test",  1, -1, 1);
		TH1D *hDummyFitJets = new TH1D(SF("hDummyFitA%d", rand()), "test",  1, -1, 1);
		TH1D *hDummyFitSJ = new TH1D(SF("hDummyFitSJ%d", rand()), "test",  1, -1, 1);
		TH1D *hDummyMRW = new TH1D(SF("hDummyMRW%d", rand()), "test",  1, -1, 1);

		hDummyFitSJ->SetLineColor(kRed); hDummyFitSJ->SetLineWidth(2);
		hDummyFitJets->SetLineColor(kMagenta); hDummyFitJets->SetLineWidth(2);
		hDummyFitA->SetLineColor(kBlack); hDummyFitA->SetLineWidth(2);
		hDummyMRW->SetLineColor(kOrange); hDummyMRW->SetLineWidth(2);

		leg3->AddEntry(hDummyFitA,  "NNLO (H1 Fit A)", "l");
		leg3->AddEntry(hDummyFitJets,  "NNLO (H1 Fit Jets)", "l");
		leg3->AddEntry(hDummyFitSJ,  "NNLO (ZEUS SJ)", "l");
		leg3->AddEntry(hDummyMRW,  "NNLO (MRW)", "l");

		leg3->Draw();


	}

	can->SaveAs( can->GetTitle() );
	can->Clear();

}












#if 0
void PlotTotal(TCanvas *can, Histogram *h1, Histogram *h2, Histogram *h3, Histogram *h4, Histogram *h5, Histogram *h6)
{
	Histogram *hist[] = {h1, h2, h3, h4, h5, h6};
	cout << "I am here " << __LINE__ << endl;

	int Nhist;
	for(Nhist=0; Nhist < static_cast<int>(sizeof(hist)/sizeof(hist[0])) && hist[Nhist]; ++Nhist)
		;


	TString *tags = new TString[Nhist];
	for(int i = 0; i < Nhist; ++i) {
		TString name = hist[i]->theor_path;
		tags[i] = name(name.First('/')+1, name.Last('/')-name.First('/')-1 );
	}

	//TString tag =  

	cout << "I am here " << __LINE__ << endl;


	int hash = rand();

	#define SF TString::Format 

	TH1D *hData = new TH1D( SF("hData%d",hash), "data", Nhist, -0.5, Nhist-0.5 );
	TH1D *hNLO  = new TH1D( SF("hNLO%d",hash),  "NLO",  Nhist, -0.5, Nhist-0.5);
	TH1D *hNNLO = new TH1D( SF("hNNLO%d",hash), "NNLO", Nhist, -0.5, Nhist-0.5);

	TH1D *hNLORatio = new TH1D( SF("hNLORatio%d",hash), "NLORatio", Nhist, -0.5, Nhist-0.5);
	TH1D *hNNLORatio = new TH1D( SF("hNNLORatio%d",hash), "NNLORatio", Nhist, -0.5, Nhist-0.5);
	TH1D *hDataRatio = new TH1D( SF("hDataRatio%d",hash), "DataRatio", Nhist, -0.5, Nhist-0.5);

	TGraphAsymmErrors *gData      = new TGraphAsymmErrors(Nhist);
	TGraphAsymmErrors *gDataAll   = new TGraphAsymmErrors(Nhist);
	TGraphAsymmErrors *gDataRatio = new TGraphAsymmErrors(Nhist);
	TGraphAsymmErrors *gDataAllRatio = new TGraphAsymmErrors(Nhist);

	TGraphAsymmErrors *gNLO       = new TGraphAsymmErrors(Nhist);
	TGraphAsymmErrors *gNLOscale   = new TGraphAsymmErrors(Nhist);
	TGraphAsymmErrors *gNLORatio  = new TGraphAsymmErrors(Nhist);
	TGraphAsymmErrors *gNLOscaleRatio=new TGraphAsymmErrors(Nhist);

	TGraphAsymmErrors *gNNLO       = new TGraphAsymmErrors(Nhist);
	TGraphAsymmErrors *gNNLOscale   = new TGraphAsymmErrors(Nhist);
	TGraphAsymmErrors *gNNLORatio  = new TGraphAsymmErrors(Nhist);
	TGraphAsymmErrors *gNNLOscaleRatio=new TGraphAsymmErrors(Nhist);


	cout << "I am here " << __LINE__ << endl;
	for(int i = 0; i < Nhist; ++i) {
		hData->GetXaxis()->SetBinLabel(i+1, tags[i]);
		hDataRatio->GetXaxis()->SetBinLabel(i+1, tags[i]);
	}


	cout << "I am here " << __LINE__ << endl;

	double dispMax = -1e30;
	double dispMin = +1e30;

	double ratMax = -1e30;
	double ratMin = +1e30;

	for(int i = 0; i < Nhist; ++i) {
		double xSec = hist[i]->data[0];
		double xSecStat = hist[i]->dataStatErr[0];
		double xSecSyst = hist[i]->dataSystErr[0];
		double totEr = hypot(xSecStat, xSecSyst );

		//vector<double> &nlo     
		//vector<double> &nloUp   
		//vector<double> &nloDown 


		vector<double> nloA(3), nnloA(3);
		double nlo, nloUp, nloDown;
		double nnlo, nnloUp, nnloDown;
		vector<TString> nloNames, nnloNames;
		vector<double> NloSyst[32];
		vector<double> NNloSyst[32];

		nlo=nloUp=nloDown= -1;
		nnlo=nnloUp=nnloDown= -1;
	
		cout << "I am here " << __LINE__ << endl;
		for(auto & th : hist[i]->theories) {
			TString name = th.first;

			if( name.Contains(defFile) ) {
				if( name.BeginsWith("nlo:") ) {
					//cout << "RADEK " << name << endl;
					if(name.EndsWith("FitB-0-Q2pPt2-cc"))
						nlo  = th.second.xsc[0];
					if(name.EndsWith("FitB-0-Q2pPt2-uu"))
						nloUp= th.second.xsc[0];
					if(name.EndsWith("FitB-0-Q2pPt2-dd"))
						nloDown =  th.second.xsc[0];

					if(name.EndsWith("FitB-0-Q2-cc"))
						nloA[0] =  th.second.xsc[0];
					if(name.EndsWith("FitB-0-Pt2-cc"))
						nloA[1] =  th.second.xsc[0];
					if(name.EndsWith("FitB-0-0.25Q2pPt2-cc"))
						nloA[2] =  th.second.xsc[0];

					for(int s = 0; s <= 30; ++s) 
						if(name.EndsWith(SF("FitB-%d-Q2pPt2-cc",s)))
							NloSyst[s] = th.second.xsc ;


					//NloNames.insert(NloNames.begin(),  name);
					//nloNow = &th.second;
				}
				if( name.BeginsWith("nnlo:")  ) {

					if(name.EndsWith("FitB-0-Q2pPt2-cc"))
						nnlo= th.second.xsc[0];
					if(name.EndsWith("FitB-0-Q2pPt2-uu"))
						nnloUp =  th.second.xsc[0];
					if(name.EndsWith("FitB-0-Q2pPt2-dd"))
						nnloDown =  th.second.xsc[0];

					if(name.EndsWith("FitB-0-Q2-cc"))
						nnloA[0] =  th.second.xsc[0];
					if(name.EndsWith("FitB-0-Pt2-cc"))
						nnloA[1] =  th.second.xsc[0];
					if(name.EndsWith("FitB-0-0.25Q2pPt2-cc"))
						nnloA[2] =  th.second.xsc[0];

					for(int s = 0; s <= 30; ++s) 
						if(name.EndsWith(SF("FitB-%d-Q2pPt2-cc",s)))
							NNloSyst[s] = th.second.xsc ;

					//NNloNames.insert(NNloNames.begin(), name);
				}
			}
		}
		assert(nlo != -1 && nloUp != -1 && nloDown != -1);
		assert(nnlo != -1 && nnloUp != -1 && nnloDown != -1);



		//vector<double> &nlo     = hist[i]->nlo.xsc;
		//vector<double> &nloUp   = hist[i]->nlo.xscUp;
		//vector<double> &nloDown = hist[i]->nlo.xscDown;

		//double NLO = nlo;
		//double NNLO = nnlo;

		hData->SetBinContent(i+1, xSec);
		hData->SetBinError(i+1, xSecStat);

		gData->SetPoint(i, hData->GetBinCenter(i+1), xSec);
		gData->SetPointError(i, 0, 0, xSecStat, xSecStat);

		gDataAll->SetPoint(i, hData->GetBinCenter(i+1), xSec );
		gDataAll->SetPointError(i, 0, 0,  totEr, totEr );


		hNLO->SetBinContent(i+1,  nlo);
		hNNLO->SetBinContent(i+1, nnlo);

		hNLORatio->SetBinContent(i+1, nlo/nlo );
		hNLORatio->SetBinError(i+1, 0 );

		hNNLORatio->SetBinContent(i+1, nnlo/nlo );
		hNNLORatio->SetBinError(i+1, 0 );

		hDataRatio->SetBinContent(i+1, xSec/nlo );
		hDataRatio->SetBinError(i+1, xSecStat/nlo  );

		gDataRatio->SetPoint(i, hData->GetBinCenter(i+1), xSec/nlo );
		gDataRatio->SetPointError(i, 0, 0, xSecStat/nlo,  xSecStat/nlo  );

		gDataAllRatio->SetPoint(i, hData->GetBinCenter(i+1), xSec/nlo );
		gDataAllRatio->SetPointError(i, 0, 0,  totEr/nlo, totEr/nlo );

		/*
		double nloSyst = 0;//hist[i]->GetSystTot(0); TODO
		double nloH    = max(nloUp - NLO, nloDown - NLO);
		double nloL    =-min(nloUp - NLO, nloDown - NLO);

		double nloTotU = hypot(nloSyst,nloL);
		double nloTotD = hypot(nloSyst,nloH);

		gNLOsyst->SetPoint(i, hData->GetBinCenter(i+1), NLO );
		gNLOsyst->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nloSyst, nloSyst );

		gNLO->SetPoint(i, hData->GetBinCenter(i+1), NLO );
		gNLO->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2,   nloTotU, nloTotD );


		gNLORatio->SetPoint(i, hData->GetBinCenter(i+1), NLO/xSec );
		gNLORatio->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nloTotU/xSec, nloTotD/xSec );

		gNLOsystRatio->SetPoint(i, hData->GetBinCenter(i+1), NLO/xSec );
		gNLOsystRatio->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nloSyst/xSec, nloSyst/xSec );

		dispMax = max(dispMax, nloTotD + NLO );
		dispMax = max(dispMax, xSec + totEr );
		*/

		//RADEK new

		//LOAD NLO

        
		auto nloSyst = GetSystTot(NloSyst,  0);
		double nloH    = max({0., nloUp - nlo, nloDown - nlo});
		double nloL    =-min({0., nloUp - nlo, nloDown - nlo});

		double nloTotU = hypot(nloSyst.first, nloL);
		double nloTotD = hypot(nloSyst.second,nloH);

		gNLOscale->SetPoint(i, hData->GetBinCenter(i+1), nlo);
		gNLOscale->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nloL, nloH );

		gNLO->SetPoint(i, hData->GetBinCenter(i+1), nlo);
		gNLO->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2,   nloTotU, nloTotD );


		gNLORatio->SetPoint(i, hData->GetBinCenter(i+1), nlo/nlo);
		gNLORatio->SetPointError(i, hData->GetBinWidth(i+1)/2., hData->GetBinWidth(i+1)/2., nloTotU/nlo, nloTotD/nlo );

		gNLOscaleRatio->SetPoint(i, hData->GetBinCenter(i+1), nlo/nlo );
		gNLOscaleRatio->SetPointError(i, hData->GetBinWidth(i+1)/2., hData->GetBinWidth(i+1)/2., nloL/nlo, nloH/nlo );

		//LOAD NNLO

		auto nnloSyst =  GetSystTot(NNloSyst,  0);
		double nnloH    = max({0., nnloUp - nnlo, nnloDown - nnlo});
		double nnloL    =-min({0., nnloUp - nnlo, nnloDown - nnlo});

		double nnloTotU = hypot(nnloSyst.first,nnloL);
		double nnloTotD = hypot(nnloSyst.second,nnloH);

		gNNLOscale->SetPoint(i, hData->GetBinCenter(i+1), nnlo );
		gNNLOscale->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nnloL, nnloH );

		gNNLO->SetPoint(i, hData->GetBinCenter(i+1), nnlo );
		gNNLO->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2,   nnloTotU, nnloTotD );


		gNNLORatio->SetPoint(i, hData->GetBinCenter(i+1), nnlo/nlo );
		gNNLORatio->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nnloTotU/nlo, nnloTotD/nlo );

		gNNLOscaleRatio->SetPoint(i, hData->GetBinCenter(i+1), nnlo/nlo );
		gNNLOscaleRatio->SetPointError(i, hData->GetBinWidth(i+1)/2, hData->GetBinWidth(i+1)/2, nnloL/nlo, nnloH/nlo );



		dispMax = max(dispMax, nloTotD + nlo );
		dispMax = max(dispMax, nnloTotD + nnlo );
		dispMax = max(dispMax, xSec + totEr );

		dispMin = min(dispMin, -nloTotU + nlo );
		dispMin = min(dispMin, -nnloTotU + nnlo );
		dispMin = min(dispMin, xSec - totEr );


		
		ratMax  = max(ratMax, nlo/nlo  + nloTotD/nlo);
		ratMax  = max(ratMax, nnlo/nlo + nnloTotD/nlo);
		ratMax  = max(ratMax,  xSec/nlo + totEr/nlo );

		ratMin  = min(ratMin, nlo/nlo  - nloTotU/nlo);
		ratMin  = min(ratMin, nnlo/nlo - nnloTotU/nlo);
		ratMin  = min(ratMin, xSec/nlo - totEr/nlo );



	}

	double leftMargin = 0.13;
	double rightMargin = 0.03;
	double bottomMargin = 0.35;

	TVirtualPad * currpad = 0;

	currpad = gPad;//canvas->cd(1);

	//currpad -> Divide(1,3,0.,0.,0); // rozdelim prvni pad na tri podsebou
	//currpad -> Divide(1,2,0.,0.,0);

	// a ty se zacnou znovu cislovat od jedne

	//currpad -> cd(1);// a prepnu se do prvniho padu 
	TPad *upPad = new TPad( SF("upPad%d",hash), "upPad", 0, 0.3, 1, 1 );
	upPad->SetBottomMargin(0);
	upPad->SetLeftMargin(leftMargin);
	upPad->SetRightMargin(rightMargin);
	upPad->Draw();
	upPad->cd();
	//gPad -> SetRightMargin(0.02);
	gPad->SetLogy();

	hData->SetLineColor(kBlack);
	hData->Draw("AXIS");
	hData->GetYaxis()->SetTitle("#sigma_{tot} [pb]");
	hData->GetYaxis()->SetTitleSize(0.06);
	hData->GetYaxis()->SetTitleOffset(1.1);
	hData->GetYaxis()->SetLabelSize(0.06);
	hData->GetYaxis()->SetMoreLogLabels();
	hData->GetYaxis()->SetNoExponent();
	hData->SetMinimum(10);
	hData->SetMaximum(1000);


	gNLO->SetLineColor(kWhite);
	gNLO->SetFillColor(kYellow);
	gNLO->Draw("e2 same");

	gNLOscale->SetLineColor(kWhite);
	gNLOscale->SetFillColor(kOrange);
	gNLOscale->Draw("e2 same");


	hNLO->SetLineColor(kWhite);
	hNLO->Draw("same");

	gNNLO->SetLineColor(kGreen);
	gNNLO->SetFillColor(kGreen);
	gNNLO->SetFillStyle(3004);
	gNNLO->Draw("e2 same");


	gNNLOscale->SetLineColor(kGreen);
	gNNLOscale->SetFillColor(kGreen);
	gNNLOscale->SetFillStyle(3004);
	gNNLOscale->Draw("e2 same");


	hNNLO->SetLineColor(kGreen);
	hNNLO->Draw("same");






	gStyle->SetEndErrorSize(4);
	gData->SetMarkerStyle(20);
	//gData->SetMarkerSize(3);
	gData->Draw("e same");

	gDataAll->SetMarkerStyle(20);
	gDataAll->Draw("pz same");


	//Plot Legend
	TLegend *leg = new TLegend( 0.5, 0.7, 0.9, 0.9 );
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);

	leg->AddEntry(gData,  "Data", "ep");
	leg->AddEntry(gNLO,   "H1 2006 Fit B" , "fl");
	leg->Draw();
	gPad->Update();
	gPad->RedrawAxis();
	TFrame *frU = (TFrame *) gPad->FindObject("TFrame");
	frU->SetFillStyle(0);
	frU->Draw();
	gPad->Update();


	currpad->cd();
	TPad *downPad = new TPad( SF("downPad%d",hash), "downPad", 0, 0.0, 1, 0.3 );
	downPad->SetTopMargin(0);
	downPad->SetLeftMargin(leftMargin);
	downPad->SetRightMargin(rightMargin);
	downPad->SetBottomMargin(bottomMargin);
	downPad->Draw();
	downPad->cd();
	//gPad -> SetRightMargin(0.02);

	double Max = max( hDataRatio->GetMaximum(), hNLORatio->GetMaximum() );
	hDataRatio->SetMaximum(Max);

	hDataRatio->GetYaxis()->SetNdivisions(204);

	hDataRatio->Draw("AXIS");
	hDataRatio->SetMinimum(0.5);
	hDataRatio->SetMaximum(2.0);


	hDataRatio->GetXaxis()->SetTitle("");
	hDataRatio->GetXaxis()->SetTitleSize(0.06 * 7./3);

	hDataRatio->GetYaxis()->SetLabelSize(0.06 * 7./3);
	hDataRatio->GetXaxis()->SetLabelSize(0.06 * 7./3);
	hDataRatio->GetXaxis()->SetLabelOffset(0.05);

	hDataRatio->GetYaxis()->SetTitle("NLO/Data");
	hDataRatio->GetYaxis()->SetTitleSize(0.06 * 7./3);
	hDataRatio->GetYaxis()->SetTitleOffset(1.1);

	gNLORatio->SetFillColor(kYellow);
	gNLORatio->Draw("e2 same");

	gNLOscaleRatio->SetFillColor(kOrange);
	gNLOscaleRatio->Draw("e2 same");


	hNLORatio->SetLineColor(kWhite);
	hNLORatio->Draw("same");


	gNNLORatio->SetLineColor(kGreen);
	gNNLORatio->SetFillColor(kGreen);
	gNNLORatio->SetFillStyle(3004);
	gNNLORatio->Draw("e2 same");


	gNNLOscaleRatio->SetLineColor(kGreen);
	gNNLOscaleRatio->SetFillColor(kGreen);
	gNNLOscaleRatio->SetFillStyle(3004);
	gNNLOscaleRatio->Draw("e2 same");


	hNNLORatio->SetLineColor(kGreen);
	hNNLORatio->Draw("same");







	gDataRatio->SetMarkerStyle(20);
	gDataRatio->Draw("p same");

	gDataAllRatio->Draw("pz same");

	gPad->Update();
	gPad->RedrawAxis();
	TFrame *frD = (TFrame *) gPad->FindObject("TFrame");
	frD->SetFillStyle(0);
	frD->Draw();
	gPad->Update();

	can->SaveAs( can->GetTitle() );
	can->Clear();
	//can->Clear();

}

#endif


int main(int argc, char** argv){

	// namespaces
	using namespace std;
	using namespace say;		// namespace for 'speaker.h'-verbosity levels
	using namespace fastNLO;	// namespace for fastNLO constants

	gStyle->SetOptStat(0);


	map<const char *, Histogram> histFPS, histVFPS, histBoris, histMozer, histSchatzel, histZEUS;
	map<const char *, Histogram>  histBorisScales;

	vector<Setting> Settings;
	Setting sItem;
	for(int i = 0; i <= 30; ++i) {
		sItem.CreateFitB(i, "Q2pPt2", 1, 1);
		Settings.push_back(sItem);
	}
	sItem.CreateFitA(0, "Q2pPt2", 1, 1);
	Settings.push_back(sItem);

	sItem.CreateFitJets(0, "Q2pPt2", 1, 1);
	Settings.push_back(sItem);

	sItem.CreateZeusSJ(0, "Q2pPt2", 1, 1);
	Settings.push_back(sItem);

	sItem.CreateMRW(0, "Q2pPt2", 1, 1);
	Settings.push_back(sItem);

	sItem.CreateFit19(0, "Q2pPt2", 1, 1);
	Settings.push_back(sItem);

    //Scale variations
	sItem.CreateFitA(0, "Q2pPt2", 1, 1, 2);//Scale up
	Settings.push_back(sItem);

	sItem.CreateFitA(0, "Q2pPt2", 1, 1, 0.5);//Scale down
	Settings.push_back(sItem);

	sItem.CreateFitJets(0, "Q2pPt2", 1, 1, 2);//Scale up
	Settings.push_back(sItem);

	sItem.CreateFitJets(0, "Q2pPt2", 1, 1, 0.5);//Scale down
	Settings.push_back(sItem);


	sItem.CreateZeusSJ(0, "Q2pPt2", 1, 1, 2);//Scale up
	Settings.push_back(sItem);

	sItem.CreateZeusSJ(0, "Q2pPt2", 1, 1, 0.5);//Scale down
	Settings.push_back(sItem);


	sItem.CreateMRW(0, "Q2pPt2", 1,1,  2);
	Settings.push_back(sItem);

	sItem.CreateMRW(0, "Q2pPt2", 1,1, 0.5);
	Settings.push_back(sItem);

	sItem.CreateFit19(0, "Q2pPt2", 1,1,  2); //Scale up
	Settings.push_back(sItem);

	sItem.CreateFit19(0, "Q2pPt2", 1,1, 0.5); //Scale up
	Settings.push_back(sItem);




	sItem.CreateFitB(0, "Q2pPt2", 1, 1, 2); //Scale up
	Settings.push_back(sItem);
	sItem.CreateFitB(0, "Q2pPt2", 1, 1, 0.5); //Scale down
	Settings.push_back(sItem);

	sItem.CreateFitB(0, "Q2", 1, 0); //Scale Q2
	Settings.push_back(sItem);
	sItem.CreateFitB(0, "Pt2", 0, 1); //Scale Pt2
	Settings.push_back(sItem);
	sItem.CreateFitB(0, "0.25Q2pPt2", 0.25, 1); //Scale Q2/4+Pt2
	Settings.push_back(sItem);


    /*
    //Scale studies
	vector<Setting> SettingsScales;

	for(int i = 0; i <= 32; ++i) {
		sItem.CreateFitA(i, "Q2pPt2", 1, 1);
		SettingsScales.push_back(sItem);
	}


    //Q2 + pT2
	sItem.CreateFitB(0, "Q2pPt2", 1, 1, 1);
	SettingsScales.push_back(sItem);

	sItem.CreateFitB(0, "Q2pPt2", 1, 1, 2);
	SettingsScales.push_back(sItem);

	sItem.CreateFitB(0, "Q2pPt2", 1, 1, 0.5);
	SettingsScales.push_back(sItem);


    //Q2
	sItem.CreateFitB(0, "Q2", 1, 0, 1); //Scale Q2
	SettingsScales.push_back(sItem);

	sItem.CreateFitB(0, "Q2", 1, 0, 2); //Scale Q2
	SettingsScales.push_back(sItem);

	sItem.CreateFitB(0, "Q2", 1, 0, 0.5); //Scale Q2
	SettingsScales.push_back(sItem);


    //pT2
	sItem.CreateFitB(0, "Pt2", 0, 1, 1); //Scale Pt2
	SettingsScales.push_back(sItem);

	sItem.CreateFitB(0, "Pt2", 0, 1, 2); //Scale Pt2
	SettingsScales.push_back(sItem);

	sItem.CreateFitB(0, "Pt2", 0, 1, 0.5); //Scale Pt2
	SettingsScales.push_back(sItem);

    //Scale Q2/4+Pt2
	sItem.CreateFitB(0, "0.25Q2pPt2", 0.25, 1, 1); //Scale Q2/4+Pt2
	SettingsScales.push_back(sItem);

	sItem.CreateFitB(0, "0.25Q2pPt2", 0.25, 1, 2); //Scale Q2/4+Pt2
	SettingsScales.push_back(sItem);

	sItem.CreateFitB(0, "0.25Q2pPt2", 0.25, 1, 0.5); //Scale Q2/4+Pt2
	SettingsScales.push_back(sItem);


    //Scale sqrt(Q2^2+Pt2^2)
	sItem.CreateFitB(0, "SqrtQ4pPt4", 1, 1, 1); 
	SettingsScales.push_back(sItem);

	sItem.CreateFitB(0, "SqrtQ4pPt4", 1, 1, 2); 
	SettingsScales.push_back(sItem);

	sItem.CreateFitB(0, "SqrtQ4pPt4", 1, 1, 0.5);
	SettingsScales.push_back(sItem);


	readHistograms(SettingsScales, histBorisScales,"../data/Measurement/boris.txt", "total");
    */



	
	//readHistograms(Settings, histVFPS,"../data/Measurement/vfps.txt",  "q2");

	//return 0;

    

	readHistograms(Settings, histVFPS,"../data/Measurement/vfps.txt",
	    "zpom", "q2", "ptjet1", "y", "deltaeta", "meaneta", "xpom", "total", "mx");


	readHistograms(Settings, histFPS, "../data/Measurement/fps_central.txt",
	    "zpom", "q2", "ptjet1", "y", "deltaetaStar", "logxpom", "total");


	readHistograms(Settings, histBoris,"../data/Measurement/boris.txt",
	    "zpom", "q2","y", "ptjet1", "ptjet2","ptavg", "deltaetaStar", "logxpom", "total");

	readHistograms(Settings, histBoris,"../data/Measurement/boris.txt",
	   "ptjet1_q2_4_6", "ptjet1_q2_6_10",  "ptjet1_q2_10_18",
	                                                                           "ptjet1_q2_18_34", "ptjet1_q2_34_100");
	readHistograms(Settings, histBoris,"../data/Measurement/boris.txt",
	   "zpom_q2_4_10", "zpom_q2_10_20", "zpom_q2_20_40", "zpom_q2_40_100");


	readHistograms(Settings, histMozer,"../data/Measurement/mozer.txt",  "zpom", "ptjet1", "y", "deltaetaStar", "logxpom", "total");

	readHistograms(Settings, histSchatzel,"../data/Measurement/schatzel.txt", "w", "ptjet1", "q2", "meaneta", "deltaetaStar", "logxpom", "total" );


	readHistograms(Settings, histZEUS,"../data/Measurement/zeus.txt",  "w", "q2", "mx", "ptJ", "etaJ", "xpom", "zpom", "total", "beta", "xgamma");

	readHistograms(Settings, histZEUS,"../data/Measurement/zeus.txt",  "zpom_ptjet1_5_6p5", "zpom_ptjet1_6p5_8", "zpom_ptjet1_8_16" );
	readHistograms(Settings, histZEUS,"../data/Measurement/zeus.txt",  "zpom_q2_5_12", "zpom_q2_12_25", "zpom_q2_25_50", "zpom_q2_50_100" );


	//TApplication *app = new TApplication(argv[0], &argc, argv);


   //map<double, vector<vector<double>>> myMap = histVFPS["ptjet1"].calculateScaleDependence(true);

    //histVFPS["ptjet1"].calculateScaleDependence("both", "newver");
	//histVFPS["ptjet1"].plotScaleDependence(1);

    //histBoris["total"].CalculateGluonFraction();

    //histZEUS["zpom_q2_50_100"].CalculateGluonFraction();
    //return 0;


    TString dir = "forNote/disTurin/";


	TCanvas * canTotDPDF= new TCanvas("canvasDPDF", dir + "totDPDF.pdf", 1000, 1000); 
	PlotTotalDPDFs(canTotDPDF, &histFPS["total"],   &histVFPS["total"],    &histBoris["total"],
	               &histMozer["total"], &histSchatzel["total"],&histZEUS["total"] );

    return 0;

	TCanvas *canTot = new TCanvas("canTot", "Total Dependence", 400, 400);
	//canTot->Print("total.pdf[" );



	histBoris["total"].plotScaleDependence(1, "fact");
	canTot->Print(dir+"muFdepLRG.pdf" );
	canTot->Clear();


	histBoris["total"].plotScaleDependence(1, "ren");
	canTot->Print(dir+"muRdepLRG.pdf" );
	canTot->Clear();

	histBoris["total"].plotScaleDependence(1, "both");
	canTot->Print(dir+"muSimdepLRG.pdf" );
	canTot->Clear();


	//canTot->Print("total.pdf" );
	//canTot->Print("total.pdf]" );

	//return 0;




	/*
	TCanvas *Can = new TCanvas("myCan", "Scale Dependence", 400, 600);

	Can->Print("radecek.pdf[" );
	plotScaleDependences(true, &histBoris["ptjet1_q2_4_6"], &histBoris["ptjet1_q2_6_10"],
		                   &histBoris["ptjet1_q2_10_18"], &histBoris["ptjet1_q2_18_34"], &histBoris["ptjet1_q2_34_100"]);

	Can->SaveAs("radecek.pdf");
	Can->Clear();
	plotScaleDependences(false, &histBoris["ptjet1_q2_4_6"], &histBoris["ptjet1_q2_6_10"],
		                   &histBoris["ptjet1_q2_10_18"], &histBoris["ptjet1_q2_18_34"], &histBoris["ptjet1_q2_34_100"]);

	Can->SaveAs("radecek.pdf");

	Can->SaveAs("radecek.pdf]");
	*/
	//return 0;


	//histVFPS["ptjet1"].plotScaleDependence();
	//histBoris["ptjet1"].plotScaleDependence();
	//histBoris["ptjet1_q2_34_100"].plotScaleDependence();
	//histBoris["logxpom"].plotScaleDependence();

	TString fileName = "DISjetsDiff.pdf";	
	TCanvas * can= new TCanvas("canvas", fileName, 1000, 400); 
	can->Print(fileName+"[" );

	//PlotFour(can, histVFPS, "q2");

	histVFPS["q2"].style.isLogY = true;
	histVFPS["q2"].style.isLogX = true;

	histFPS["ptjet1"].style.isLogY = true;
	histVFPS["ptjet1"].style.isLogY = true;
	histMozer["ptjet1"].style.isLogY = true;
	histSchatzel["ptjet1"].style.isLogY = true;
	histBoris["ptjet1"].style.isLogY = true;
	histBoris["ptjet2"].style.isLogY = true;
	histZEUS["ptJ"].style.isLogY = true;

	histBoris["ptavg"].style.isLogY = true;
	histFPS["q2"].style.isLogY = true;
	histBoris["q2"].style.isLogY = true;
	histSchatzel["q2"].style.isLogY = true;
	histZEUS["q2"].style.isLogY = true;
	histZEUS["beta"].style.isLogY = false;
	histZEUS["beta"].style.isLogX = true;
	histZEUS["xgamma"].style.isLogY = false;



	histVFPS["q2"].style.isLogX = true;
	histFPS["q2"].style.isLogX = true;
	histBoris["q2"].style.isLogX = true;
	histSchatzel["q2"].style.isLogX = true;
	histZEUS["q2"].style.isLogX = true;





	//LegendPos
	histFPS["y"].style.legX = 0.65;
	histVFPS["y"].style.legX = 0.65;
	histBoris["y"].style.legX = 0.65;
	histBoris["zpom"].style.legX = 0.65;
	histBoris["ptjet2"].style.legX = 0.60;
	histVFPS["deltaeta"].style.legX = 0.65;
	histVFPS["zpom"].style.legX = 0.65;
	histVFPS["xpom"].style.legY = 0.15;
	histVFPS["mx"].style.legX = 0.17;
	histVFPS["meaneta"].style.legX = 0.17;
	histFPS["logxpom"].style.legX = 0.4;

	histVFPS["total"].style.legY = 0.13;
	histFPS["total"].style.legY = 0.13;
	histBoris["total"].style.legY = 0.13;
	histSchatzel["total"].style.legY = 0.13;
	histMozer["total"].style.legY = 0.13;
	histZEUS["total"].style.legY = 0.13;
	histVFPS["total"].style.legX = 0.6;
	histFPS["total"].style.legX = 0.6;
	histBoris["total"].style.legX = 0.6;
	histSchatzel["total"].style.legX = 0.6;
	histMozer["total"].style.legX = 0.6;
	histZEUS["total"].style.legX = 0.6;

	histBoris["zpom_q2_4_10"].style.legX = 0.65;
	histBoris["zpom_q2_10_20"].style.legX = 0.65;
	histBoris["zpom_q2_20_40"].style.legX = 0.65;
	histBoris["zpom_q2_40_100"].style.legX = 0.65;

	histBoris["ptjet1_q2_4_6"].style.isLogY = true;
	histBoris["ptjet1_q2_6_10"].style.isLogY = true;
	histBoris["ptjet1_q2_10_18"].style.isLogY = true;
	histBoris["ptjet1_q2_18_34"].style.isLogY = true;
	histBoris["ptjet1_q2_34_100"].style.isLogY = true;


	histSchatzel["w"].style.legX = 0.6;
	histSchatzel["logxpom"].style.legX = 0.4;
	histMozer["logxpom"].style.legX = 0.4;
	histMozer["y"].style.legX = 0.6;
	histMozer["zpom"].style.legX = 0.6;

	histZEUS["mx"].style.legX = 0.6;
	histZEUS["zpom"].style.legX = 0.6;
	histZEUS["xpom"].style.legY = 0.15;
	histZEUS["beta"].style.legX = 0.60;
	histZEUS["beta"].style.legY = 0.65;

	histZEUS["xgamma"].style.legX = 0.20;
	histZEUS["xgamma"].style.legY = 0.65;


	histZEUS["w"].style.legX = 0.37;
	histZEUS["w"].style.legY = 0.10;


	histZEUS["etaJ"].style.legX = 0.55;
	histZEUS["zpom_ptjet1_5_6p5"].style.legX = 0.65;
	histZEUS["zpom_ptjet1_6p5_8"].style.legX = 0.65;
	histZEUS["zpom_ptjet1_8_16"].style.legX = 0.65;
	histZEUS["zpom_q2_5_12"].style.legX = 0.65;
	histZEUS["zpom_q2_12_25"].style.legX = 0.65;
	histZEUS["zpom_q2_25_50"].style.legX = 0.65;
	histZEUS["zpom_q2_50_100"].style.legX = 0.65;

    //Ratios
	TCanvas * canPtTest= new TCanvas("canvasPtTest", dir+"ptjet1DiffTest.pdf", 1000, 400); 
    //canPtTest->SaveAs(canPtTest->GetTitle() + TString("["));

    canPtTest->SetTitle(dir+"ptjet1.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {&histFPS["ptjet1"], &histVFPS["ptjet1"],  &histBoris["ptjet1"],
	               &histMozer["ptjet1"], &histSchatzel["ptjet1"], &histZEUS["ptJ"] });

    #if 1
    canPtTest->SetTitle(dir+"y.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {&histFPS["y"],   &histVFPS["y"],    &histBoris["y"],
	               &histMozer["y"], &histSchatzel["w"], &histZEUS["w"] });
    canPtTest->SetTitle(dir+"q2.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {&histFPS["q2"],   &histVFPS["q2"],    &histBoris["q2"],
	                nullptr, &histSchatzel["q2"], &histZEUS["q2"] });
    canPtTest->SetTitle(dir+"zpom.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {&histFPS["zpom"],   &histVFPS["zpom"],    &histBoris["zpom"],
	               &histMozer["zpom"], nullptr, &histZEUS["zpom"] });


    canPtTest->SetTitle(dir+"xpom.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {&histFPS["logxpom"],   &histVFPS["xpom"],    &histBoris["logxpom"],
	               &histMozer["logxpom"], &histSchatzel["logxpom"], &histZEUS["xpom"] });
    canPtTest->SetTitle(dir+"deltaeta.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {&histFPS["deltaetaStar"],   &histVFPS["deltaeta"],    &histBoris["deltaetaStar"],
	               &histMozer["deltaetaStar"], &histSchatzel["deltaetaStar"], nullptr });
    canPtTest->SetTitle(dir+"meaneta.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {nullptr,   &histVFPS["meaneta"],    nullptr,
	               nullptr, &histSchatzel["meaneta"], &histZEUS["etaJ"] });


    canPtTest->SetTitle(dir+"ptjet1_q2_LRG.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {  &histBoris["ptjet1_q2_4_6"], &histBoris["ptjet1_q2_6_10"],  &histBoris["ptjet1_q2_10_18"],
	                       &histBoris["ptjet1_q2_18_34"], &histBoris["ptjet1_q2_34_100"], nullptr} );

    canPtTest->SetTitle(dir+"zpom_q2_LRG.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {  
	   &histBoris["zpom_q2_4_10"], &histBoris["zpom_q2_10_20"], &histBoris["zpom_q2_20_40"], &histBoris["zpom_q2_40_100"], nullptr, nullptr});


    canPtTest->SetTitle(dir+"zpom_ptjet1_ZEUS.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {  
	&histZEUS["zpom_ptjet1_5_6p5"], &histZEUS["zpom_ptjet1_6p5_8"], &histZEUS["zpom_ptjet1_8_16"], nullptr, nullptr, nullptr} );
    canPtTest->SetTitle(dir+"zpom_q2_ZEUS.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {&histZEUS["zpom_q2_5_12"], &histZEUS["zpom_q2_12_25"], &histZEUS["zpom_q2_25_50"], &histZEUS["zpom_q2_50_100"],
	 nullptr, nullptr} );

    canPtTest->SetTitle(dir+"q2_Scale.pdf");
	PlotComparisonWithAbs(canPtTest, "Scale", {&histFPS["q2"],   &histVFPS["q2"],    &histBoris["q2"],
	                nullptr, &histSchatzel["q2"], &histZEUS["q2"] });

    canPtTest->SetTitle(dir+"zpom_DPDFs.pdf");
	PlotComparisonWithAbs(canPtTest, "DPDF", {&histFPS["zpom"],   &histVFPS["zpom"],    &histBoris["zpom"],
	               &histMozer["zpom"], nullptr, &histZEUS["zpom"] });

    canPtTest->SetTitle(dir+"mx.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {  &histVFPS["mx"], &histZEUS["mx"],  nullptr,
	                       nullptr, nullptr, nullptr} );
    
    canPtTest->SetTitle(dir+"ptavg2_LRG.pdf");
	PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {  &histBoris["ptavg"], &histBoris["ptjet2"],  nullptr,
	                       nullptr, nullptr, nullptr} );


    canPtTest->SetTitle(dir+"ptjet1_Scale.pdf");
	PlotComparisonWithAbs(canPtTest, "Scale", {&histFPS["ptjet1"], &histVFPS["ptjet1"],  &histBoris["ptjet1"],
	               &histMozer["ptjet1"], &histSchatzel["ptjet1"], &histZEUS["ptJ"] });
    canPtTest->SetTitle(dir+"ptjet1_DPDFs.pdf");
	PlotComparisonWithAbs(canPtTest, "DPDF", {&histFPS["ptjet1"], &histVFPS["ptjet1"],  &histBoris["ptjet1"],
	               &histMozer["ptjet1"], &histSchatzel["ptjet1"], &histZEUS["ptJ"] });

    #endif

    //PlotSingle(canPtTest, "NLOvsNNLO", &histZEUS["beta"] );
    //PlotSingle(canPtTest, "NLOvsNNLO", &histZEUS["xgamma"] );
	//PlotComparisonWithAbs(canPtTest, "NLOvsNNLO", {  &histZEUS["beta"], nullptr,  nullptr,
	                        //nullptr, nullptr, nullptr} );

    
    //PlotFour(canPtTest, "NLOvsNNLO", histZEUS, "beta", "zpom", "ptJ", "q2" );

    //canPtTest->SaveAs(canPtTest->GetTitle() + TString("]"));


	//TCanvas * canNew= new TCanvas("canvasNew", dir+"BetaXgamma.pdf", 1000, 1000); 
    //PlotFour(canNew, "NLOvsNNLO", histZEUS, "beta", "xgamma", "ptJ", "q2" );

    #if 0
    return 0;


	TCanvas * canPt= new TCanvas("canvasPt", dir+"ptjet1Diff.pdf", 1000, 400); 
	PlotComparison(canPt, "NLOvsNNLO", {&histFPS["ptjet1"],   &histVFPS["ptjet1"],    &histBoris["ptjet1"],
	               &histMozer["ptjet1"], &histSchatzel["ptjet1"], &histZEUS["ptJ"] });



	TCanvas * canY= new TCanvas("canvasY", dir+"yDiff.pdf", 1000, 400); 
	PlotComparison(canY, "NLOvsNNLO", {&histFPS["y"],   &histVFPS["y"],    &histBoris["y"],
	               &histMozer["y"], &histSchatzel["w"], &histZEUS["w"] });

	TCanvas * canQ2= new TCanvas("canvasQ2", dir+"q2Diff.pdf", 1000, 400); 
	PlotComparison(canQ2, "NLOvsNNLO", {&histFPS["q2"],   &histVFPS["q2"],    &histBoris["q2"],
	                nullptr, &histSchatzel["q2"], &histZEUS["q2"] });

	TCanvas * canQ2scale= new TCanvas("canvasQ2scale", dir+"q2ScaleDiff.pdf", 1000, 400); 
	PlotComparison(canQ2scale, "Scale", {&histFPS["q2"],   &histVFPS["q2"],    &histBoris["q2"],
	                nullptr, &histSchatzel["q2"], &histZEUS["q2"] });


	TCanvas * canZpom= new TCanvas("canvasZpom", dir+"zpomDiff.pdf", 1000, 400); 
	PlotComparison(canZpom, "NLOvsNNLO", {&histFPS["zpom"],   &histVFPS["zpom"],    &histBoris["zpom"],
	               &histMozer["zpom"], nullptr, &histZEUS["zpom"] });

	TCanvas * canZpomDPDF= new TCanvas("canvasZpomDPDF", dir+"zpomDPDFDiff.pdf", 1000, 400); 
	PlotComparison(canZpomDPDF, "DPDF", {&histFPS["zpom"],   &histVFPS["zpom"],    &histBoris["zpom"],
	               &histMozer["zpom"], nullptr, &histZEUS["zpom"] });



	TCanvas * canXpom= new TCanvas("canvasXpom", dir+"xpomDff.pdf", 1000, 400); 
	PlotComparison(canXpom, "NLOvsNNLO", {&histFPS["logxpom"],   &histVFPS["xpom"],    &histBoris["logxpom"],
	               &histMozer["logxpom"], &histSchatzel["logxpom"], &histZEUS["xpom"] });

	TCanvas * canDeltaEta= new TCanvas("canvasDeltaEta", dir+"DeltaEtaDff.pdf", 1000, 400); 
	PlotComparison(canDeltaEta, "NLOvsNNLO", {&histFPS["deltaetaStar"],   &histVFPS["deltaeta"],    &histBoris["deltaetaStar"],
	               &histMozer["deltaetaStar"], &histSchatzel["deltaetaStar"], nullptr });

	TCanvas * canMeanEta= new TCanvas("canvasMeanEta", dir+"MeanEtaDff.pdf", 1000, 400); 
	PlotComparison(canMeanEta, "NLOvsNNLO", {nullptr,   &histVFPS["meaneta"],    nullptr,
	               nullptr, &histSchatzel["meaneta"], &histZEUS["etaJ"] });


	
	TCanvas * canPtjetQ2= new TCanvas("canvasPtjetQ2", dir+"PtjetQ2.pdf", 1000, 400); 
	PlotComparison(canPtjetQ2, "NLOvsNNLO", {  &histBoris["ptjet1_q2_4_6"], &histBoris["ptjet1_q2_6_10"],  &histBoris["ptjet1_q2_10_18"],
	                       &histBoris["ptjet1_q2_18_34"], &histBoris["ptjet1_q2_34_100"], nullptr} );

	PlotComparison(can, "NLOvsNNLO", {  
	   &histBoris["zpom_q2_4_10"], &histBoris["zpom_q2_10_20"], &histBoris["zpom_q2_20_40"], &histBoris["zpom_q2_40_100"], nullptr, nullptr});


	PlotComparison(can, "NLOvsNNLO", {  
	&histZEUS["zpom_ptjet1_5_6p5"], &histZEUS["zpom_ptjet1_6p5_8"], &histZEUS["zpom_ptjet1_8_16"], nullptr, nullptr, nullptr} );
	PlotComparison(can, "NLOvsNNLO", {&histZEUS["zpom_q2_5_12"], &histZEUS["zpom_q2_12_25"], &histZEUS["zpom_q2_25_50"], &histZEUS["zpom_q2_50_100"],
	 nullptr, nullptr} );

	PlotComparison(can, "NLOvsNNLO", {  
	&histVFPS["mx"], &histBoris["ptjet2"], &histBoris["ptavg"], &histZEUS["mx"], &histZEUS["beta"], nullptr} );

	cout << "LOOOOOOOOK" << endl;
	cout << "Fit          NLO  NNLO" << endl;
	cout << "Fit B:      " << ScoreFitB[0] <<" "<< ScoreFitB[1] << endl;
	cout << "Fit A:      " << ScoreFitA[0] <<" "<< ScoreFitA[1] << endl;
	cout << "Fit Jets:   " << ScoreJets[0] <<" "<< ScoreJets[1] << endl;
	cout << "ZEUS SJ:    " << ScoreSJ[0] <<" "<< ScoreSJ[1] << endl;
	cout << "Fit Bup:    " << ScoreFitBup[0] <<" "<< ScoreFitBup[1] << endl;
	cout << "Fit Bdn:    " <<ScoreFitBdn[0] <<" "<< ScoreFitBdn[1] << endl;
	cout << "Fit JetsUp: " << ScoreJetsup[0] <<" "<< ScoreJetsup[1] << endl;
	cout << "Fit JetsDn: " <<ScoreJetsdn[0] <<" "<< ScoreJetsdn[1] << endl;

	can->SaveAs(fileName+"]");
	cout << endl;


	fileName = "DISjetsTot.pdf";	
	TCanvas * dan= new TCanvas("canvasTot", fileName, 1000, 1000); 
	dan->Print(fileName+"[" );
    #endif




	/*
	vector<TString> types;
	types = {"NLOvsNNLO", "Scales", "DPDF"};

	for(TString type : types) {
		PlotFour(can, type, histFPS, "q2", "ptjet1", "y", "deltaetaStar");
		PlotFour(can, type, histFPS, "zpom", "logxpom", "total");

		PlotFour(can, type, histVFPS, "q2", "ptjet1", "y", "deltaeta");
		PlotFour(can, type, histVFPS, "meaneta", "xpom", "mx", "zpom");
		PlotFour(can, type, histVFPS,  "total");



		PlotFour(can, type, histBoris, "q2", "y", "logxpom", "zpom" );
		PlotFour(can, type, histBoris, "ptjet1", "ptjet2", "ptavg", "deltaetaStar"  );
		PlotFour(can, type, histBoris, "ptjet1_q2_4_6", "ptjet1_q2_6_10", "ptjet1_q2_10_18", "ptjet1_q2_18_34");
		PlotFour(can, type, histBoris, "ptjet1_q2_34_100");
		PlotFour(can, type, histBoris, "zpom_q2_4_10", "zpom_q2_10_20", "zpom_q2_20_40", "zpom_q2_40_100");
		PlotFour(can, type, histBoris, "total" );


		
		PlotFour(can, type, histSchatzel,  "ptjet1", "q2", "deltaetaStar", "meaneta" );
		PlotFour(can, type, histSchatzel,  "w", "logxpom", "total" );


		PlotFour(can, type, histMozer, "y", "deltaetaStar", "logxpom", "ptjet1" );
		PlotFour(can, type, histMozer, "zpom", "total" );

		PlotFour(can, type, histZEUS, "q2", "mx", "xpom", "zpom");
		PlotFour(can, type, histZEUS,  "ptJ", "etaJ", "beta", "w" );
		PlotFour(can, type, histZEUS,  "total" );

		PlotFour(can, type, histZEUS, "zpom_ptjet1_5_6p5", "zpom_ptjet1_6p5_8", "zpom_ptjet1_8_16" );
		PlotFour(can, type, histZEUS, "zpom_q2_5_12", "zpom_q2_12_25", "zpom_q2_25_50", "zpom_q2_50_100" );
	}
	*/




	TCanvas * canTotNLOnsNNLO= new TCanvas("canvasNLOvsNNLO", dir + "totNLOvsNNLO.pdf", 1000, 1000); 
	PlotTotal(canTotNLOnsNNLO, &histFPS["total"],   &histVFPS["total"],    &histBoris["total"],
	               &histMozer["total"], &histSchatzel["total"],&histZEUS["total"] );

	TCanvas * canTotDPDFold= new TCanvas("canvasDPDF", dir + "totDPDF.pdf", 1000, 1000); 
	PlotTotalDPDFs(canTotDPDFold, &histFPS["total"],   &histVFPS["total"],    &histBoris["total"],
	               &histMozer["total"], &histSchatzel["total"],&histZEUS["total"] );

	TCanvas * canTotScales= new TCanvas("canvasScales", dir + "totScales.pdf", 1000, 1000); 
	plotScaleStudiesTotal(canTotScales, &histFPS["total"],   &histVFPS["total"],    &histBoris["total"],
	               &histMozer["total"], &histSchatzel["total"],&histZEUS["total"] );


	//app->Run();
	//app->SetReturnFromRun(true);


	/*
	*/

	/*
	can->Divide(2,2);
	can->cd(1);
	plotXi12("../data/Measurement/fps_central.txt", "logxpom", "zpom");
	can->cd(2);
	plotXi12("../data/Measurement/vfps.txt", "xpom", "zpom");
	can->cd(3);
	plotXi12("../data/Measurement/boris.txt", "logxpom", "zpom");
	can->cd(4);
	plotXi12("../data/Measurement/mozer.txt", "logxpom", "zpom");
	can->SaveAs( can->GetTitle() );
	can->Clear();
	*/

	//dan->SaveAs(fileName+"]");



	return 0;


}

void Histogram::ConvertToW(Theory &th, vector<double> &fastBinsLo, vector<double> &fastBinsHi)
{
	auto &xsc     = th.xsc;

    th.Print();
    cout << "Bin Sizes "<< th.file_name<<" " << fastBinsLo.size() <<" "<< fastBinsHi.size()<<" "<< xMin.size()  << endl;
	assert(fastBinsLo.size() == xMin.size() );

	//Check binning compatibility
	for(unsigned i = 0; i < xMin.size(); ++i) {
		double s = 4*27.6*920;
		if(theor_path.Contains("/LRGH1820/"))
			s = 4*27.6*820;
		double wL= sqrt( fastBinsLo[i]*s);
		double wH= sqrt( fastBinsHi[i]*s);
		if( abs(wL-xMin[i]) > 5 || abs(wH-xMax[i]) > 5 ) {
			cout << "Huge difference for "<< i << endl;
			exit(1);
		}
	}

	for(unsigned i = 0; i < xsc.size(); ++i) {
		double factor = (fastBinsHi[i]-fastBinsLo[i]) / (xMax[i]-xMin[i]);
		xsc[i] *= factor;
	}
	fastBinsLo = xMin;
	fastBinsHi = xMax;

}

void Histogram::ConvertToTotal(Theory &th, vector<double> &fastBinsLo, vector<double> &fastBinsHi)
{
	auto &xsc     = th.xsc;

	for(unsigned i = 0; i < xsc.size(); ++i) {
		double bw = fastBinsHi[i]-fastBinsLo[i];
		if(i == 0) bw -= 1;
		xsc[0] += xsc[i] * bw;
	}

	xsc.resize(1);

	fastBinsLo = xMin;
	fastBinsHi = xMax;

}


Theory Histogram::CalculateVarNLO( fastNLODiffAlphas  &fnlodiff,  vector<double> &fastBinsLo, vector<double> &fastBinsHi,
							std::function<double(double,double,double)>  func) //xpom,q2,tableVar
{

	//  If you want to receive your cross section in
	//   pb/GeV or in pb. Here we choose pb/GeV
	fnlodiff.SetUnits(fastNLO::kPublicationUnits);


	//fnlodiff.SetExternalFuncForMuR (&Function_Mu);
	//fnlodiff.SetExternalFuncForMuF (&Function_Mu);



	fastBinsLo = fnlodiff.GetObsBinsLoBounds(0);
	fastBinsHi = fnlodiff.GetObsBinsUpBounds(0);


	// calculate and access the cross section
	typedef std::map<double,  std::vector < std::map< double, double > > >  array3D;

	array3D xsQ2;
	//std::map<double,  std::vector < std::map< double, double > > > xsQ2 = new array3D[3+npdfall];

	double *bins = new double[nBins+1];
	for(int i = 0; i < nBins; ++i)
		bins[i] = xMin[i];
	bins[nBins] = xMax[nBins-1];

	int hash = rand();

	TH1D *hist = new TH1D(TString::Format("hist%d", hash), "histogram", nBins, bins);

	//fnlodiff.SetLHAPDFMember(0);

	//fnlodiff.SetScaleFactorsMuRMuF(1.0, 1.0);
	vector<double> NloXs=fnlodiff.GetDiffCrossSection();
	xsQ2 = fnlodiff.Get3DCrossSection();

	/*
	double Sum=0;
	for(int i = 0; i <NloXs.size(); ++i) {
		Sum+= NloXs[i] * (fastBinsHi[i]-fastBinsLo[i]);
		cout << i <<" : "<<  fastBinsLo[i]<<" "<<fastBinsHi[i]<< "  <>  "<<   NloXs[i] << endl;
	}
	cout << "Sum is " << Sum << endl;
	exit(0);
	*/



	//For interpretation
	TMatrixD mat = CreateInterpolationMatrix(fastBinsLo, fastBinsHi);
	TVectorD nloVec(fastBinsLo.size() );


	//Fill array of Q2
	vector<double> Q2arr;
	for(const auto &a :   (xsQ2.begin()->second)[0])
		Q2arr.push_back(a.first);

	/*
	for(int i = 0; i <Q2arr.size(); ++i) {
		Sum+= NloXs[i] * (fastBinsHi[i]-fastBinsLo[i]);
		cout << i <<" : "<<  fastBinsLo[i]<<" "<<fastBinsHi[i]<< "  <>  "<<   NloXs[i] << endl;
	}
	*/



	//TH1D *hNLOold = new TH1D("histOld","histOld", 300, 10, 40 );
	//TH1D *hNLOnew = new TH1D("histNew","histNew", 300, 10, 40 );

	//const double s = 4*27.6*920;




	const int Niter = 4000;

    srand(1);


	for( const auto &v : xsQ2 ) {
		double xpom = v.first;
		auto xs2D   = v.second;

		for(unsigned i = 0; i < fastBinsLo.size(); ++i) {
			for(const auto &item : xs2D[i]) {
				double Q2 = item.first;
				double xs = item.second*(fastBinsHi[i]-fastBinsLo[i]) / Niter;
				for(int k = 0; k < Niter; ++k) {
					double varHist = fastBinsLo[i] + rand()/(RAND_MAX+0.) * (fastBinsHi[i]-fastBinsLo[i]);
					double VarCalc = func(xpom, Q2, varHist);
					//cout << "Zpom " << VarCalc <<  endl; 
					hist->Fill(VarCalc, xs);
				}
			}
		}
	}




	/*
	for( const auto &v : xsQ2[0] ) {
		double xpom = v.first;
		auto xs2D   = v.second;

		for(int q = 0; q < Q2arr.size(); ++q) {
			double Q2 = Q2arr[q];
			for(int sys = 0; sys < npdfall+3; ++sys) {
				for(int i = 0; i < fastBinsLo.size(); ++i) {
					nloVec[i] = xsQ2[sys][xpom][i][Q2];
					if( xsQ2[sys][xpom][i].size() != Q2arr.size() )
						cout << "Sizes "<<xsQ2[sys][xpom][i].size() <<" "<< Q2arr.size()  << endl;
				}
				TVectorD nloNew = mat*nloVec;
				
				for(int i = 0; i < fastBinsLo.size(); ++i) {
					double xs = nloVec[i]*(fastBinsHi[i]-fastBinsLo[i]) / Niter;
					for(int k = 0; k < Niter; ++k) {
						double varHist = GetRandom(fastBinsLo, fastBinsHi, nloNew, i);
						double VarCalc = func(xpom, Q2, varHist);
						hist[sys]->Fill(VarCalc, xs);

						//if(sys == 0) {
							//hNLOnew->Fill(VarCalc, xs); 
							//double varHist = fastBinsLo[i] + rand()/(RAND_MAX+0.) * (fastBinsHi[i]-fastBinsLo[i]);
							//double VarCalc = func(xpom, Q2, varHist);
							//hNLOold->Fill(VarCalc, xs); 
						//}
					}
				}
			}
		}

	}
	*/

	/*
	TFile *file = TFile::Open("test.root", "recreate");
	hNLOnew->Write();
	hNLOold->Write();
	file->Write();
	exit(0);
	*/

	Theory th;

	th.xsc.resize(nBins);


	for(int i = 0; i < nBins; ++i) {
		double corr = 1./ (xMax[i]-xMin[i]);
		th.xsc[i]=hist->GetBinContent(i+1) * corr;
	}
	delete hist;


	fastBinsLo = xMin;
	fastBinsHi = xMax;

	return th;
}

Theory Histogram::CalculateXpomNLO( fastNLODiffAlphas  &fnlodiff,  vector<double> &fastBinsLo, vector<double> &fastBinsHi  )
{
	Theory th;

	vector<double>	nloXpom(nBins);

	//cout << "RADEKHERE " << __LINE__ << endl;


	for(int i = 0; i < nBins; ++i) {
	//cout << "RADEKHERE " << __LINE__ << endl;

		if(var_name.EqualTo("xpom")) {
			fnlodiff.SetXPomLinSlicing( 5, xMin[i], xMax[i]); // VFPS range
			fnlodiff.SettIntegratedRange(-0.6);
		}
		else if(var_name.EqualTo("logxpom")) {
			if( theor_path.Contains("/LRG_H1/") )
				fnlodiff.SetXPomLogSlicing( 3, pow(10,xMin[i]) ,  pow(10,xMax[i]) ); 
			else
				fnlodiff.SetXPomLogSlicing( 5, pow(10,xMin[i]) ,  pow(10,xMax[i]) ); 

			cout << "RADEK " << i<<" "<<pow(10,xMin[i]) <<" "<< pow(10,xMax[i]) << endl;
			fnlodiff.SettIntegratedRange(-1.);
	//cout << "RADEKHERE " << __LINE__ << endl;
			
		}
		else {
			cout << "error in xpom " <<var_name <<  endl;
		}

		th = CalculateRawNLO( fnlodiff,  fastBinsLo,  fastBinsHi  );

		cout <<"Helenka " <<  th.xsc[0] << endl;
		double bw = xMax[i] - xMin[i];

		nloXpom[i]   = th.xsc[0] / bw;

	}
	//Save it
	th.xsc = nloXpom;

	return th;

}






//__________________________________________________________________________________________________________________________________


double Function_Mu(double s1, double s2 ){
	// --- fastNLO user: This is an example function
	//     to demonstrate how you might perform the
	//     definition of the scales using a 
	//     'flexible-scale'-table
	//double mu = s1*exp(0.3*s2);
	//q2/4 + meanPt^2
	double mu = sqrt(s1*s1/4. + s2*s2);
	return mu;
}

//__________________________________________________________________________________________________________________________________



TMatrixD Histogram::CreateInterpolationMatrix(vector<double> &fastBinsLo, vector<double> &fastBinsHi)
{
	const int Nb = fastBinsLo.size();

	assert(Nb >2);
	
	vector<double> wb(Nb);
	for(int i = 0; i < Nb; ++i)
		wb[i] = fastBinsHi[i] - fastBinsLo[i];

	TMatrixD mat(Nb,Nb);

	mat(0,0)     = 1;
	mat(Nb-1,Nb-1) = 1;

	for(int i = 1; i < Nb-1; ++i) {
		mat(i,i-1) = 1./4 * wb[i] / ( wb[i-1] + wb[i] );
		mat(i,i+1) = 1./4 * wb[i] / ( wb[i+1] + wb[i] );
		mat(i,i  ) = 1./4 * ( 2 + wb[i-1] / (wb[i-1]+wb[i]) + wb[i+1] / (wb[i+1]+wb[i]) );
	}

	return mat.Invert();

}

double GetRandom(vector<double> &fastBinsLo, vector<double> &fastBinsHi, TVectorD &vec, int i)
{
	int n = fastBinsHi.size();
	double w0=0., w2=0.; //dummy values
	double  w1 = fastBinsHi[i] - fastBinsLo[i];
	if(i>0) w0 = fastBinsHi[i-1] - fastBinsLo[i-1];
	if(i<n-1) w2 = fastBinsHi[i+1] - fastBinsLo[i+1];


	if(i == 0) {
		double h1 = (w1*vec[1]+w2*vec[0]) / ( w1 + w2 );
		double h0 = vec[0] - w1/(w2+w1) * (vec[1]-vec[0]);
		double Min = min( {0.0, h0, h1} );
		h0+= -Min;
		h1+= -Min;
		double Max = max(h0, h1);
		double xi, y, yFunc;
		do {
			xi =  rand()/(RAND_MAX+0.);
			yFunc = (1-xi)*h0 + xi*h1;
			y = Max * rand()/(RAND_MAX+0.);
		} while(y > yFunc);

		return ( fastBinsLo[i] + (fastBinsHi[i]-fastBinsLo[i])*xi );
	}
	else if(i == n-1) {
		double h1 = (w0*vec[i]+w1*vec[i-1]) / ( w1 + w0 );
		double h2 = vec[i] + w1/(w0+w1) * (vec[i]-vec[i-1]);
		double Min = min( {0., h2, h1} );
		h2+= -Min;
		h1+= -Min;
		double Max = max(h2, h1);
		double xi, y, yFunc;
		do {
			xi =  rand()/(RAND_MAX+0.);
			yFunc = (1-xi)*h1 + xi*h2;
			y = Max * rand()/(RAND_MAX+0.);
		} while(y > yFunc);

		return ( fastBinsLo[i] + (fastBinsHi[i]-fastBinsLo[i])*xi );
	}
	else if(i > 0 && i < n-1) {
		double hl = (w0*vec[i]+w1*vec[i-1]) / (w1+w0);
		double hr = (w2*vec[i]+w1*vec[i+1]) / (w1+w2);
		double Min = min({0.0, hl, vec[i], hr} );

		double hc;
		hl += -Min;
		hc  = vec[i]-Min;
		hr += -Min;
		double Max = max({hl, hc, hr});

		double xi, y, yFunc;
		do {
			xi = rand()/(RAND_MAX+0.);
			if(xi < 0.5)
				yFunc = (1-2*xi)*hl + 2*xi*hc;
			else
				yFunc = (2*xi-1)*hr + (2-2*xi)*hc;

			y = Max * rand()/(RAND_MAX+0.);

		} while(y > yFunc);

		return ( fastBinsLo[i] + (fastBinsHi[i]-fastBinsLo[i])*xi );

	}
	else {
		assert(0);
		return 0;
	}

}
