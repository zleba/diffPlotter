#include "fastDiffLib.h"


int main(int argc, char** argv){

	// namespaces
	using namespace std;
	using namespace say;		// namespace for 'speaker.h'-verbosity levels
	using namespace fastNLO;	// namespace for fastNLO constants


	map<const char *, Histogram> histFPS, histVFPS, histBoris, histMozer, histSchatzel, histZEUS;
	map<const char *, Histogram>  histBorisScales;


    TString pdfName = "H1_DPDF_2006B_NLO";

    Histogram::dpdfs = new map<TString,DPDFset> ({
        { "H1_DPDF_2006B_NLO", DPDFset("H1_DPDF_2006B_NLO_pom", "H1_DPDF_2006B_NLO_reg") }
    });



	vector<Setting> Settings;
	Setting sItem;
    for(int i = 0; i < Histogram::dpdfs->at(pdfName).size(); ++i) {
        sItem.CreateLHAFit(pdfName, i, "Q2pPt2", 1, 1, 1); //Scale up
        Settings.push_back(sItem);
    }
	sItem.CreateLHAFit(pdfName, 0, "Q2pPt2", 1, 1, 2); //Scale up
	Settings.push_back(sItem);
	sItem.CreateLHAFit(pdfName, 0, "Q2pPt2", 1, 1, 0.5); //Scale down
	Settings.push_back(sItem);



	readHistograms(Settings, histVFPS,"../data/Measurement/vfps.txt",
	     "q2", "ptjet1", "y", "deltaeta", "meaneta", "xpom", "total");

    return 0;
	readHistograms(Settings, histVFPS,"../data/Measurement/vfps.txt",
	    "zpom", "q2", "ptjet1", "y", "deltaeta", "meaneta", "xpom", "total", "mx");

    return 0;

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


	return 0;


}
