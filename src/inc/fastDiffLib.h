#ifndef fastDiff
#define fastDiff

#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <iostream>
#include <functional>
#include "TString.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "fastNLODiffAlphas.h"

#define SF TString::Format 

struct Setting {
	double renFactor, factFactor;
    std::function<double(double,double)> renFunc, factFunc;
	TString scaleTag;
	TString pdfName;
	int pdfVecId;
    bool onlyQuarks = false;

	TString getFileName() {
        std::map<double,char> scaleMap;
		scaleMap[0.5]='d';
		scaleMap[1]='c';
		scaleMap[2]='u';

		return pdfName+"-"+SF("%d",pdfVecId)+"-"+scaleTag+"-"+scaleMap[renFactor]+scaleMap[factFactor];
	}

	void CreateLHAFit(TString fitName, int pdfID, TString ScaleTag, double CoefQ2, double CoefPt2, double sclFact=1) {
		renFunc=factFunc=[=](double q, double pt) {return  sqrt(CoefQ2*q*q + CoefPt2*pt*pt);};
		renFactor = factFactor = sclFact;
		scaleTag = ScaleTag;
		pdfName = fitName;
		pdfVecId = pdfID;
	}

};

struct Theory {
	int iOrder; //1--LO, 2--NLO, 3--NNLO
	TString file_name;
    std::vector<double> xsc;
	Setting setting;

	Theory & operator+=(const Theory &th) {
        //std::assert(xsc.size() == th.xsc.size() );
		for(unsigned i = 0; i < xsc.size(); ++i) {
			xsc[i]     += th.xsc[i];
		}
		return *this;
	}

	void ReverseBinning() {
		auto reverseAll = [](std::vector<double> &v) { reverse(begin(v), end(v) ); };
		reverseAll(xsc);
	}

	void Print() const {
		for(double Xsc : xsc) 
			std::cout << Xsc << " ";
        std::cout << endl;
	}
};

struct HistoErr;

class Histogram {

public:
	void loadData(TString File_name, TString Var_name);
	Theory LoadTheory(TString theor_file, Setting setting);

	void LoadTheories() {
		for(TString file : theor_files) {
            std::vector<Theory> thVec(settings.size());
#pragma omp parallel for
			for(unsigned i = 0; i < settings.size(); ++i) {
                thVec[i] = LoadTheory(file, settings[i]);
			}

			for(unsigned i = 0; i < settings.size(); ++i) {
				TString setTag = settings[i].getFileName();
				theories[file+":"+setTag] = thVec[i];
			}

        }
	}

    std::vector<double> LoadHistogramsNloNnloData(TString refName, HistoErr &grData, HistoErr &grNlo, HistoErr &grNnlo, bool sysInside = false);

    HistoErr LoadThHistogram(TString refName, bool sysInside);
    HistoErr LoadDataHistogram();

	void LoadHistogramsFitJetsSJ(HistoErr &grNloJets, HistoErr &grNnloJets,
											  HistoErr &grNloSJ, HistoErr &grNnloSJ);
	HistoErr LoadHistograms(TString name, TString file, TString histName, TString refName);

	void plotDPDFStudies();

    std::vector<double> GetBinning();
	int getNbins() {return nBins;}

    void plotNLOvsNNLO(vector<TString> typeNames);
	void plotNLOvsNNLOratio(TString plotStyle, double minY, double maxY);
	void plotNLOvsNNLOratioAbs(TString plotStyle, double minY, double maxY, double minYabs, double maxYabs, double Factors);
	void plotScaleChoices();
	void plotScaleStudies();

    void  CalculateGluonFraction();
    std::map<double, std::vector<std::vector<double>>>  calculateScaleDependence(TString scaleTag, TString tag);

	TString getTag() { return theor_path(theor_path.First('/')+1, theor_path.Last('/')-theor_path.First('/')-1 ); }

	friend void PlotComparison(TCanvas *can,  TString plotStyle,  std::vector<Histogram*> histos);
	friend void PlotComparisonWithAbs(TCanvas *can,  TString plotStyle,  std::vector<Histogram*> histos);
	friend void PlotSingle(TCanvas *can,  TString plotStyle,  Histogram* hist);
	friend void plotScaleDependences(TString scaleTag, Histogram *h1, Histogram *h2, Histogram *h3, Histogram *h4, Histogram *h5);
	friend void PlotTotal(TCanvas *can, Histogram *h1, Histogram *h2, Histogram *h3, Histogram *h4, Histogram *h5, Histogram *h6);
	friend void plotXi12(TString fileName, const char *xpomN, const char *zpomN);
	void plotScaleDependence(int binId, TString scaleTag);

	void IncludeSettings(std::vector<Setting> settings_) {settings=settings_;}

	//Histo Style
	struct {
		bool isLogY{false}, isLogX{false};
		double legX{0.5}, legY{0.6};

	} style;

    static map<TString,DPDFset> *dpdfs;
private:

	TString data_file;
	TString theor_path;
	//TString nlo_file;
	//TString nnlo_file;
    std::set<TString> theor_files;
	TString var_name;
	TString xTitle;
	TString yTitle;
	TString  Title;

	int nBins;

	double dissFactor;

    std::vector<double> xMin;
    std::vector<double> xMax;

    std::vector<double> data;
    std::vector<double> dataStatErr;
    std::vector<double> dataSystErr;
    std::vector<double> radCorr;
    std::vector<double> hadrCorr;

	//Theory nlo, nnlo;
    std::map<TString, Theory> theories;
	//vector<double> nlo;
	//vector<double> nloUp;
	//vector<double> nloDown;
	//vector< vector<double> > nloSyst;


    std::vector<Setting> settings;

	Theory CalcTheory(TString theor_file, Setting setting);
	//double GetSystTot(int binId) const;


	void ConvertToTotal(Theory &th, std::vector<double> &fastBinsLo, std::vector<double> &fastBinsHi);//For total and xpom
	void ConvertToW(Theory &th, std::vector<double> &fastBinsLo, std::vector<double> &fastBinsHi);//For total and xpom
	Theory CalculateRawNLO( fastNLODiffAlphas  &fnlodiff,  std::vector<double> &fastBinsLo, std::vector<double> &fastBinsHi  );//without HadCorr
	Theory CalculateXpomNLO( fastNLODiffAlphas  &fnlodiff, std::vector<double> &fastBinsLo, std::vector<double> &fastBinsHi  );//Xpom Calc
	Theory CalculateVarNLO( fastNLODiffAlphas  &fnlodiff,  std::vector<double> &fastBinsLo, std::vector<double> &fastBinsHi,
	                     std::function<double(double,double,double)>  func); //xpom,q2,tableVar

	TMatrixD CreateInterpolationMatrix(std::vector<double> &fastBinsLo, std::vector<double> &fastBinsHi);

};


void readHistograms(std::vector<Setting> setting, std::map<const char *, Histogram> &hist, TString fileName,
                   const char *n1,   const char *n2=0, const char *n3=0, const char *n4=0, const char *n5=0, const char *n6=0,
                   const char *n7=0, const char *n8=0, const char *n9=0, const char *n10=0);

void PlotFour(TCanvas *can, vector<TString> Type, map<const char *, Histogram> hists, const char *n1, const char *n2=0, const char *n3=0, const char *n4=0 );

#endif
