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

#include "fastDiffLib.h"

using namespace PlottingHelper;

//TString defFile = "H1-LQall-8c.diff";
const TString defFile = "newver";
//TString tablesDir = "/home/radek/moje/daniel/tables/";
const TString tablesDir = "../tables/"; //  /home/radek/moje/daniel/tables/";



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

map<TString,DPDFset> *Histogram::dpdfs = nullptr;

#define SF TString::Format 


// Function prototype for flexible-scale function 
double Function_Mu(double s1, double s2 );

class Histogram;

struct HistoErr;

void PlotTotal(TCanvas *can, Histogram *h1, Histogram *h2=0, Histogram *h3=0, Histogram *h4=0, Histogram *h5=0, Histogram *h6=0);

void PlotTotalDPDFs(TCanvas *can, Histogram *h1, Histogram *h2=0, Histogram *h3=0, Histogram *h4=0, Histogram *h5=0, Histogram *h6=0);

double GetRandom(vector<double> &fastBinsLo, vector<double> &fastBinsHi, TVectorD &vec, int i);


void plotXi12(TString fileName, const char *xpomN, const char *zpomN);




//__________________________________________________________________________________________________________________________________
pair<double,double> GetSystTot(vector<double> *sysVec, int binId) ;




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

Theory Histogram::CalcTheory(TString theor_file, Setting setting)
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
	TheoryFile.ReplaceAll("/nlo/", "/nnlo/");

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

        fnlodiff.setDPDF(dpdfs->at(setting.pdfName));

        /*
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
        */

		fnlodiff.SetLHAPDFMember(setting.pdfVecId);

        fnlodiff.IncludeOnlyQuarks(setting.onlyQuarks);
	}

	//reference to first entry
	auto & fnlodiff = fnlodiffs[0];

	
	//cout << "long " << sizeof(int) << " "<< sizeof(long)<<" "<<sizeof(long long) << endl;

	//cout<<"DANIEL " << fastNLO::kExtern<<endl;
	//exit(1);


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
			*/
		}
		catch (std::ifstream::failure e) {
			std::cerr << "Exception opening/reading/closing file\n";
		}

		nloXsec.close();
		th.setting = setting;
	}
	else {
		th = CalcTheory(theor_file, setting);


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



double getBestChi2(TString anal, TString var, vector<double> thVec);



void readHistograms(vector<Setting> setting, map<const char *, Histogram> &hist, TString fileName,
                   const char *n1,   const char *n2, const char *n3, const char *n4, const char *n5, const char *n6,
                   const char *n7, const char *n8, const char *n9, const char *n10)
{
	const char *names[] = {n1, n2, n3, n4, n5, n6, n7, n8, n9, n10};

	for(unsigned i = 0; i < sizeof(names)/sizeof(names[0]); ++i) {
		if(names[i] == 0)
			continue;
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


double getChi2(TH1D *hData, TGraphAsymmErrors *grAll, TH1D *hModel)
{
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


    cout << "RADEK " << __LINE__ << endl;
	// calculate and access the cross section
	typedef std::map<double,  std::vector < std::map< double, double > > >  array3D;

	array3D xsQ2;
	//std::map<double,  std::vector < std::map< double, double > > > xsQ2 = new array3D[3+npdfall];

	double *bins = new double[nBins+1];
	for(int i = 0; i < nBins; ++i)
		bins[i] = xMin[i];
	bins[nBins] = xMax[nBins-1];

    cout << "RADEK " << __LINE__ << endl;
	int hash = rand();

    cout << "RADEK " << __LINE__ << endl;
	TH1D *hist = new TH1D(TString::Format("hist%d", hash), "histogram", nBins, bins);
    cout << "RADEK " << __LINE__ << endl;

	//fnlodiff.SetLHAPDFMember(0);

	//fnlodiff.SetScaleFactorsMuRMuF(1.0, 1.0);
	vector<double> NloXs=fnlodiff.GetDiffCrossSection();
	xsQ2 = fnlodiff.Get3DCrossSection();

    cout << "RADEK " << __LINE__ << endl;
	/*
	double Sum=0;
	for(int i = 0; i <NloXs.size(); ++i) {
		Sum+= NloXs[i] * (fastBinsHi[i]-fastBinsLo[i]);
		cout << i <<" : "<<  fastBinsLo[i]<<" "<<fastBinsHi[i]<< "  <>  "<<   NloXs[i] << endl;
	}
	cout << "Sum is " << Sum << endl;
	exit(0);
	*/

    cout << "RADEK " << __LINE__ << endl;


	//For interpretation
	TMatrixD mat = CreateInterpolationMatrix(fastBinsLo, fastBinsHi);
	TVectorD nloVec(fastBinsLo.size() );


    cout << "RADEK " << __LINE__ << endl;
	//Fill array of Q2
	vector<double> Q2arr;
	for(const auto &a :   (xsQ2.begin()->second)[0])
		Q2arr.push_back(a.first);

    cout << "RADEK " << __LINE__ << endl;




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
