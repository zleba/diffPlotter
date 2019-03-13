#include "fastDiffLib.h"

vector<Setting> getSettings(TString pdfName)
{
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
    return Settings;
}


int main(int argc, char** argv){


	map<const char *, Histogram> histFPS, histVFPS, histBoris, histMozer, histSchatzel, histZEUS;
	map<const char *, Histogram>  histBorisScales;


    //TString pdfName = "H1_DPDF_2006B_NLO";
    TString pdfName = "lhaTest";

    Histogram::dpdfs = new map<TString,DPDFset> ({
        { "H1_DPDF_2006B_NLO", DPDFset("H1_DPDF_2006B_NLO_pom", "H1_DPDF_2006B_NLO_reg") },
        { "lhaTest", DPDFset("lhaTest_pom", "lhaTest_reg") }
    });

    /*
    vector<Setting> settings1 = getSettings("lhaTest");
    vector<Setting> settings2 = getSettings("H1_DPDF_2006B_NLO");
    vector<Setting> Settings(settings1.size() );
    //Settings.insert( Settings.end(), settings1.begin(), settings1.end() );
    Settings.insert( Settings.end(), settings1.begin(), settings1.end() );
    */


	vector<Setting> Settings;
    pdfName = "lhaTest";
	Setting sItem;
    for(int i = 0; i < Histogram::dpdfs->at(pdfName).size(); ++i) {
        sItem.CreateLHAFit(pdfName, i, "Q2pPt2", 1, 1, 1); //Scale up
        Settings.push_back(sItem);
    }
	sItem.CreateLHAFit(pdfName, 0, "Q2pPt2", 1, 1, 2); //Scale up
	Settings.push_back(sItem);
	sItem.CreateLHAFit(pdfName, 0, "Q2pPt2", 1, 1, 0.5); //Scale down
	Settings.push_back(sItem);

    pdfName =  "H1_DPDF_2006B_NLO";
	//Setting sItem;
    for(int i = 0; i < Histogram::dpdfs->at(pdfName).size(); ++i) {
        sItem.CreateLHAFit(pdfName, i, "Q2pPt2", 1, 1, 1); //Scale up
        Settings.push_back(sItem);
    }
	sItem.CreateLHAFit(pdfName, 0, "Q2pPt2", 1, 1, 2); //Scale up
	Settings.push_back(sItem);
	sItem.CreateLHAFit(pdfName, 0, "Q2pPt2", 1, 1, 0.5); //Scale down
	Settings.push_back(sItem);


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

    //return 0;

    //Style
	histVFPS.at("q2").style.isLogY = true;
	histVFPS.at("q2").style.isLogX = true;

	histFPS.at("ptjet1").style.isLogY = true;
	histVFPS.at("ptjet1").style.isLogY = true;
	histMozer.at("ptjet1").style.isLogY = true;
	histSchatzel.at("ptjet1").style.isLogY = true;
	histBoris.at("ptjet1").style.isLogY = true;
	histBoris.at("ptjet2").style.isLogY = true;
	histZEUS.at("ptJ").style.isLogY = true;

	histBoris.at("ptavg").style.isLogY = true;
	histFPS.at("q2").style.isLogY = true;
	histBoris.at("q2").style.isLogY = true;
	histSchatzel.at("q2").style.isLogY = true;
	histZEUS.at("q2").style.isLogY = true;
	histZEUS.at("beta").style.isLogY = false;
	histZEUS.at("beta").style.isLogX = true;
	histZEUS.at("xgamma").style.isLogY = false;

	histBoris.at("ptjet1_q2_4_6").style.isLogY = true;
	histBoris.at("ptjet1_q2_6_10").style.isLogY = true;
	histBoris.at("ptjet1_q2_10_18").style.isLogY = true;
	histBoris.at("ptjet1_q2_18_34").style.isLogY = true;
	histBoris.at("ptjet1_q2_34_100").style.isLogY = true;


	histVFPS.at("q2").style.isLogX = true;
	histFPS.at("q2").style.isLogX = true;
	histBoris.at("q2").style.isLogX = true;
	histSchatzel.at("q2").style.isLogX = true;
	histZEUS.at("q2").style.isLogX = true;



    //histVFPS["q2"].plotNLOvsNNLO("nlo:H1_DPDF_2006B_NLO-Q2pPt2");

    TCanvas *can = new TCanvas("can", "test.pdf");
    can->SaveAs(can->GetTitle()+TString("["));
    //PlotFour(can, {"nlo:H1_DPDF_2006B_NLO-Q2pPt2", "nnlo:H1_DPDF_2006B_NLO-Q2pPt2"}, histVFPS, "q2", "ptjet1", "y", "deltaeta");

    vector<TString> thVec = {"nnlo:lhaTest-Q2pPt2", "nnlo:H1_DPDF_2006B_NLO-Q2pPt2"};
    PlotFour(can, thVec, histVFPS, "zpom", "q2", "ptjet1", "y");
    PlotFour(can, thVec, histVFPS, "deltaeta", "meaneta", "xpom", "mx");
    PlotFour(can, thVec, histVFPS, "total");

    PlotFour(can, thVec, histFPS,  "zpom", "q2", "ptjet1", "y");
    PlotFour(can, thVec, histFPS,  "deltaetaStar", "logxpom", "total");

    PlotFour(can, thVec, histBoris, "zpom", "q2","y", "logxpom");
    PlotFour(can, thVec, histBoris, "ptjet1", "ptjet2", "ptavg", "deltaetaStar");
    PlotFour(can, thVec, histBoris,  "ptjet1_q2_4_6", "ptjet1_q2_6_10",  "ptjet1_q2_10_18", "ptjet1_q2_18_34");
    PlotFour(can, thVec, histBoris,  "ptjet1_q2_34_100");
    PlotFour(can, thVec, histBoris,  "zpom_q2_4_10", "zpom_q2_10_20", "zpom_q2_20_40", "zpom_q2_40_100");
    PlotFour(can, thVec, histBoris,  "total");


    PlotFour(can, thVec, histMozer,  "zpom", "ptjet1", "y", "deltaetaStar");
    PlotFour(can, thVec, histMozer,  "logxpom", "total");

	PlotFour(can, thVec, histSchatzel,  "ptjet1", "q2", "meaneta", "deltaetaStar");
    PlotFour(can, thVec, histSchatzel,  "w", "logxpom", "total" );

    PlotFour(can, thVec, histZEUS, "q2",  "ptJ", "etaJ", "xpom");

    PlotFour(can, thVec, histZEUS,  "zpom_ptjet1_5_6p5", "zpom_ptjet1_6p5_8", "zpom_ptjet1_8_16" );
    PlotFour(can, thVec, histZEUS,   "zpom_q2_5_12", "zpom_q2_12_25", "zpom_q2_25_50", "zpom_q2_50_100" );

    PlotFour(can, thVec, histZEUS,   "total" );


    can->SaveAs(can->GetTitle()+TString("]"));
    return 0;

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


	return 0;


}
