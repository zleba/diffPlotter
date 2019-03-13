#include "fastDiffLib.h"

#include "TLegend.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"


#include "TStyle.h"
#include "TPad.h"
#include "TVirtualPad.h"

#include <iostream> 
#include <regex> 
#include <string> 
#include <cassert> 

#include "plottingHelper.h"
using namespace PlottingHelper;
using namespace std;

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
        assert(nBins > 0);

        cout << "RADEK " << __LINE__ <<" "<< xMin.size() <<" "<< xMax.size()<< endl;
		hBins.resize(nBins+1);
		for(int i = 0; i < nBins; ++i)
			hBins[i] = xMin[i];
		hBins[nBins] = xMax[nBins-1];
        cout << "RADEK " << __LINE__ << endl;

		int hash = rand();

		h = new TH1D( SF("h%s%d",name.Data(), hash), name, nBins, hBins.data());

        cout << "RADEK " << __LINE__ << endl;
		hRatio = new TH1D( SF("h%sRatio%d",name.Data(),hash), SF("%sRatio",name.Data()), nBins, hBins.data());

        cout << "RADEK " << __LINE__ << endl;
		gr      = new TGraphAsymmErrors(nBins);
		grAll   = new TGraphAsymmErrors(nBins);
		grRatio = new TGraphAsymmErrors(nBins);
		grAllRatio = new TGraphAsymmErrors(nBins);
        cout << "RADEK " << __LINE__ << endl;
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

    double getMax() {
        double Max = -1000;
        for(int i = 0; i < gr->GetN(); ++i) {
            double x, y;
            grAll->GetPoint(i, x, y);
            double er = grAll->GetErrorYhigh(i);
            Max = max(Max, y+er);
        }
        return Max;
    }


};











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




static void SetNLOstyle(HistoErr &grNlo)
{
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
static void SetNNLOstyle(HistoErr &grNnlo)
{
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


void SetDataStyle(HistoErr &grData) {
	double pointSize = 0.6;
    //Setting style Data
    grData.gr->SetMarkerStyle(20);
    grData.gr->SetMarkerSize(pointSize);
    CopyStyle(grData.gr, grData.grRatio);

    grData.grAll->SetMarkerStyle(20);
    grData.grAll->SetMarkerSize(pointSize);
    CopyStyle(grData.grAll, grData.grAllRatio);


    grData.h->SetLineColor(kBlack);
}



pair<double,double> GetSystTot(const vector<vector<double>> &sysVec, int binId) 
{
	double sum2Up = 0;
	double sum2Dn = 0;

	for(unsigned i = 1; i < sysVec.size(); ++i) {
		double d = sysVec[i][binId] - sysVec[0][binId];

        sum2Up += pow(max(0.0, d), 2);
        sum2Dn += pow(max(0.0,-d), 2);

		//sum2 += d*d;
	}
	//return sqrt(sum2/2.);// TODO, dummy now
    return make_pair(sqrt(sum2Up), sqrt(sum2Dn));
}



HistoErr Histogram::LoadDataHistogram()
{
    assert(xMin.size() > 0 && xMax.size() > 0);
    HistoErr grData;
	grData.Init("Data", xMin, xMax);

	for(unsigned i = 0; i < xMin.size(); ++i) {
        double xRef = data[i];
		grData.SetBin(i+1, data[i], xRef, dataStatErr[i], dataStatErr[i],
		                                  dataSystErr[i], dataSystErr[i]);
    }
    return grData;
}


HistoErr Histogram::LoadThHistogram(TString refName, TString whatInside)
{

	TString refTag   = refName(0, refName.First(':')+1);
	TString refEnd   = refName(refName.First(':')+1, 10000);
    TString pdfName  = refEnd(0, refEnd.First('-')+0);
    TString scaleFun = refEnd(refEnd.First('-')+1, 10000);

    cout << "refName: " << refTag << endl;
    cout << "scaleFun: " << scaleFun << endl;
    cout << "pdfName: " << pdfName << endl;


	vector<double>  Th,  ThUp, ThDown;
	vector<vector<double>> ThSyst;

	//vector<TString> NloNames, NNloNames;
    int isFound = 0;
	for(auto & th : theories) {
		TString name = th.first;
        cout << name << endl;

        if(name.BeginsWith(refTag) && name.Contains(pdfName+"-") && name.Contains(scaleFun)) {
                ++isFound;
				if(name.EndsWith("-0-Q2pPt2-cc"))
					Th = th.second.xsc ;
				if(name.EndsWith("-0-Q2pPt2-uu"))
					ThUp =  th.second.xsc;
				if(name.EndsWith("-0-Q2pPt2-dd"))
					ThDown =  th.second.xsc;

                regex r("-[0-9]+-");
                smatch m; 

                string nStr = name.Data();
                regex_search(nStr, m, r);
                assert(m.size() > 0);
                cout <<"Radek " <<  m[0] << endl;
                int errId = stoi(string(m[0]).substr(1, string(m[0]).size()-2));

                if(name.EndsWith(SF("-%d-Q2pPt2-cc",errId))) {
                    if(errId >= ThSyst.size())
                        ThSyst.resize(errId+1);
                    ThSyst[errId] = th.second.xsc;
                }

                /*
				for(int s = 0; s <= 10; ++s) 
					if(name.EndsWith(SF("-%d-Q2pPt2-cc",s)))
						ThSyst[s] = th.second.xsc ;
                */
                ++isFound;
        }
	}
    cout << "I am after" << ThSyst.size() <<  endl;

    //HistoErr err;
    //return err;
    //return {0.0};

    if(isFound<2) {
        cout << "Nothing found for "<< refName << endl;
        exit(0);
    }


	//grData.Init("Data", xMin, xMax);
    HistoErr grTh;
	grTh.Init("Th", xMin, xMax);


	for(unsigned i = 0; i < xMin.size(); ++i) {

		//grData.SetBin(i+1, data[i], refXsc[i], dataStatErr[i], dataStatErr[i],
		                                   //dataSystErr[i], dataSystErr[i]);

        double refXsc = data[i];

		//LOAD NNLO
        double thSystUp=0, thSystDn = 0;
		tie(thSystUp,thSystDn) = GetSystTot(ThSyst, i); //without systematics
		double thErUp    = max({0., ThUp[i] - Th[i], ThDown[i] - Th[i]});
		double thErDown  =-min({0., ThUp[i] - Th[i], ThDown[i] - Th[i]});

		if(whatInside.Contains("scaleInside"))
			grTh.SetBin(i+1, Th[i], refXsc, thErDown, thErUp,
			                                       thSystDn, thSystUp);
		else if(whatInside.Contains("sysInside"))
			grTh.SetBin(i+1, Th[i], refXsc, thSystDn, thSystUp,
			                                       thErDown, thErUp);
        else
            assert(0);



	}
    return grTh;

}





void Histogram::plotNLOvsNNLO(vector<TString> typeNames)
{
    if(xMin.size() == 0 || xMax.size() == 0) {
        cout << "Histogram failed : " << var_name << endl;
        assert(0);
    }

    gStyle->SetOptStat(0);

    vector<HistoErr> grTh;
	//vector<double> range = LoadHistogramsNloNnloData("nlo:FitB-0-Q2pPt2-cc", grData, grNlo, grNnlo, false);

	HistoErr grData = LoadDataHistogram();
    for(auto n : typeNames) {
        grTh.push_back(LoadThHistogram(n, "sysInside"));
        //grTh[0] =  LoadThHistogram("nlo:H1_DPDF_2006B_NLO-Q2pPt2", false);
    }

	//Setting NNLO style

    SetDataStyle(grData);
	SetNNLOstyle(grTh[0]);
    if(grTh.size() >= 2)
        SetNLOstyle(grTh[1]);

    //gPad->SetLeftMargin(0.17);
    SetLeftRight(0.17, 0.04);
    SetTopBottom(0.05, 0.17);
    TVirtualPad *can = gPad;
    //DividePad( {1}, {1.3,1});
    DivideTransparent( {1}, {1.4,0,1});

	//myCan.SetDataStyle(grData, style.isLogX, style.isLogY);


    cout << "Radek " << xTitle <<" "<< yTitle << " " << Title << endl;
    auto Decorator = [&]() {
        SetFTO({12}, {10}, {1.3, 2.2, 0.4, 3.3});
        GetXaxis()->SetNdivisions(303);
        GetYaxis()->SetNdivisions(303);
        GetFrame()->SetTitle("");

    };


	////////////////////////////////////////////
	//Up Frame
	////////////////////////////////////////////
    can->cd(1);

    if(style.isLogY) gPad->SetLogy();
    if(style.isLogX) gPad->SetLogx();
    /*
	if(!style.isLogY) {
		grData.h->SetMinimum(0);
		grData.h->SetMaximum(range[1] * 1.1);
	}
	else {
		grData.h->SetMinimum(range[0] * 0.8);
		grData.h->SetMaximum(range[1] * 1.2);
	}
    */

    grTh[0].h->Draw("hist");

    for(int i = grTh.size()-1; i >= 0; --i) {
        grTh[i].grAll->Draw("e2 same");
        grTh[i].gr->Draw("e2 same");
    }
    for(int i = grTh.size()-1; i >= 0; --i) {
        grTh[i].h->Draw("same ][");
    }


	gStyle->SetEndErrorSize(4);
	grData.gr->Draw("e same");
	grData.grAll->Draw("pz same");


    GetYaxis()->SetTitle(yTitle);
    Decorator();
    GetXaxis()->SetLabelOffset(1000);

    double Max = -100;
    for(auto &th : grTh)
        Max = max(Max, th.getMax());
    GetFrame()->SetMaximum(1.2* max(grData.getMax(), Max));
    if(!style.isLogY) GetFrame()->SetMinimum(0);

	//Plot Legend
	TLegend *leg = newLegend(kPos9);
	TString tag =  theor_path(theor_path.First('/')+1, theor_path.Last('/')-theor_path.First('/')-1 );

    if(names2D.count(var_name) > 0 ) {
        TString n = names2D.at(var_name);
        leg->SetHeader(n);
    }
	leg->AddEntry(grData.gr,    SF("H1 %s data", tag.Data()), "ep");
	leg->AddEntry(grTh[0].grAll,     "H1 Fit2019 NNLO" , "fl");
	leg->AddEntry(grTh[1].grAll,     "H1 Fit2006B NLO" , "fl");
    DrawLegends({leg});
	//leg->Draw();


	////////////////////////////////////////////
	//Down Frame
	////////////////////////////////////////////
    can->cd(2);
    if(style.isLogX) gPad->SetLogx();

    grTh[0].hRatio->Draw("hist");


    for(int i = grTh.size()-1; i >= 0; --i) {
        grTh[i].grAllRatio->Draw("e2 same");
        grTh[i].grRatio->Draw("e2 same");
    }
    for(int i = grTh.size()-1; i >= 0; --i) {
        grTh[i].hRatio->Draw("same ][");
    }


	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");

    GetYaxis()->SetRangeUser(0,2);
    GetXaxis()->SetTitle(xTitle);
    GetYaxis()->SetTitle("#sigma/#sigma_{data}");
    GetYaxis()->CenterTitle();

    Decorator();

    //can->SaveAs("plot.pdf");
}

void PlotFour(TCanvas *can, vector<TString> Type, map<const char *, Histogram> hists, const char *n1, const char *n2, const char *n3, const char *n4)
{
	can->Divide(2,2,0.00001,0.00001); 

	const char *names[] = {n1, n2, n3, n4};

	for(int i = 0; i < 4 && names[i] != 0; ++i) {
		can->cd(i+1);
		try {
            cout << names[i] << endl;
            hists.at(names[i]).plotNLOvsNNLO(Type);
		}
		catch(const std::out_of_range& oor) {
			cout << "Variable " << names[i]<<" not known" << endl;
			exit(1);
		}
	}
	can->SaveAs( can->GetTitle() );
	can->Clear();
}
