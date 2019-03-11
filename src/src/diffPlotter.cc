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

#include "plottingHelper.h"
using namespace PlottingHelper;
using namespace std;



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

/*
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
*/


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
    HistoErr grData;
	grData.Init("Data", xMin, xMax);

	for(unsigned i = 0; i < xMin.size(); ++i) {
        double xRef = data[i];
		grData.SetBin(i+1, data[i], xRef, dataStatErr[i], dataStatErr[i],
		                                  dataSystErr[i], dataSystErr[i]);
    }
    return grData;
}


HistoErr Histogram::LoadThHistogram(TString refName, bool sysInside)
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

    cout << "Radek " << __LINE__ << endl;

	for(unsigned i = 0; i < xMin.size(); ++i) {

		//grData.SetBin(i+1, data[i], refXsc[i], dataStatErr[i], dataStatErr[i],
		                                   //dataSystErr[i], dataSystErr[i]);

    cout << "Radek " << __LINE__ << endl;
        double refXsc = data[i];

		//LOAD NNLO
        double thSystUp=0, thSystDn = 0;
		tie(thSystUp,thSystDn) = GetSystTot(ThSyst, i); //without systematics
		double thErUp    = max({0., ThUp[i] - Th[i], ThDown[i] - Th[i]});
		double thErDown  =-min({0., ThUp[i] - Th[i], ThDown[i] - Th[i]});

    cout << "Radek " << __LINE__ << endl;
		if(sysInside == false)
			grTh.SetBin(i+1, Th[i], refXsc, thErDown, thErUp,
			                                       thSystDn, thSystUp);
		else
			grTh.SetBin(i+1, Th[i], refXsc, thSystDn, thSystUp,
			                                       thErDown, thErUp);
    cout << "Radek " << __LINE__ << endl;


	}
    cout << "Radek " << __LINE__ << endl;
    return grTh;

}





void Histogram::plotNLOvsNNLO(vector<TString> typeNames)
{

    gStyle->SetOptStat(0);

    vector<HistoErr> grTh;
	//vector<double> range = LoadHistogramsNloNnloData("nlo:FitB-0-Q2pPt2-cc", grData, grNlo, grNnlo, false);

	HistoErr grData = LoadDataHistogram();
    for(auto n : typeNames) {
        grTh.push_back(LoadThHistogram(n, false));
        //grTh[0] =  LoadThHistogram("nlo:H1_DPDF_2006B_NLO-Q2pPt2", false);
    }
    cout << "Radek " << __LINE__ << endl;

	//Setting NNLO style
	SetNNLOstyle(grTh[0]);
    if(grTh.size() >= 2)
        SetNLOstyle(grTh[1]);

    cout << "Radek " << __LINE__ << endl;
    gPad->SetLeftMargin(0.17);
    TVirtualPad *can = gPad;
    DividePad( {1}, {1,1});

	//myCan.SetDataStyle(grData, style.isLogX, style.isLogY);


    auto Decorator = []() {
        SetFTO({15}, {10}, {1.1, 2.1, 0.4, 2.3});
        GetXaxis()->SetNdivisions(303);
        GetYaxis()->SetNdivisions(303);
        GetFrame()->SetTitle("");
    };


	////////////////////////////////////////////
	//Up Frame
	////////////////////////////////////////////
    cout << "Radek " << __LINE__ << endl;
    can->cd(1);

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

    cout << "Radek " << __LINE__ << endl;
    grTh[0].h->Draw("hist");

    for(HistoErr  &th : grTh) {
        th.grAll->Draw("e2 same");
        th.gr->Draw("e2 same");
        th.h->Draw("same");
    }

    cout << "Radek " << __LINE__ << endl;

	gStyle->SetEndErrorSize(4);
	grData.gr->Draw("e same");
	grData.grAll->Draw("pz same");


    Decorator();

    GetFrame()->SetMaximum(1.2* max(grData.getMax(), grTh[0].getMax()));

	//Plot Legend

	TLegend *leg = new TLegend( 0.3, 0.5, 0.3, 0.5);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
    /*
	TString tag =  theor_path(theor_path.First('/')+1, theor_path.Last('/')-theor_path.First('/')-1 );
	leg->SetHeader(Title);
	leg->AddEntry(grData.gr,    SF("H1 %s data", tag.Data()), "ep");

	leg->AddEntry(grNlo.grAll,     SF("NLO x %.2f",dissFactor) , "fl");
	leg->AddEntry(grNnlo.grAll,   SF("NNLO  x %.2f", dissFactor) , "fl");
    */

    cout << "Radek " << __LINE__ << endl;
	leg->Draw();


	////////////////////////////////////////////
	//Down Frame
	////////////////////////////////////////////
    can->cd(2);
    grTh[0].hRatio->Draw("hist");

    for(HistoErr  &th : grTh) {
        th.grAllRatio->Draw("e2 same");
        th.grRatio->Draw("e2 same");
        th.hRatio->Draw("same");
    }
    cout << "Radek " << __LINE__ << endl;

	grData.grRatio->Draw("p same");
	grData.grAllRatio->Draw("pz same");

    GetYaxis()->SetRangeUser(0,2);

    Decorator();

    cout << "Radek " << __LINE__ << endl;
    //can->SaveAs("plot.pdf");
}

void PlotFour(TCanvas *can, vector<TString> Type, map<const char *, Histogram> hists, const char *n1, const char *n2, const char *n3, const char *n4)
{
	can->Divide(2,2,0.00001,0.00001); 

	const char *names[] = {n1, n2, n3, n4};

	for(int i = 0; i < 4 && names[i] != 0; ++i) {
		can->cd(i+1);
		try {
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
