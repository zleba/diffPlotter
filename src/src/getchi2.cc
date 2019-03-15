#include <iostream>
#include "TString.h"
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <armadillo>
#include <functional>
#include <map>


using namespace std;


double getChi2(arma::mat V, arma::vec data, arma::vec theor);
double getNorm(arma::mat V, arma::vec data, arma::vec theor);

void readVFPS(ifstream &matFile, arma::mat &errMat, arma::vec &data);
void readFPS(ifstream &matFile,  arma::mat &errMat, arma::vec &data);
void readLRG(ifstream &matFile,  arma::mat &errMat, arma::vec &data);

void readMozer(ifstream &matFile,  arma::mat &errMat, arma::vec &data);
void readSchatzel(ifstream &matFile, arma::mat &errMat, arma::vec &data);
void readNoMat(TString type, ifstream &matFile, arma::mat &errMat, arma::vec &data);

void readZEUS(ifstream &matFile,  arma::mat &errMat, arma::vec &data);


map<TString, std::function<void(ifstream&, arma::mat&, arma::vec&)>> Reader = {
    {"VFPS",     readVFPS},
    {"LRG",      readLRG},
    {"FPS",      readFPS},
    {"LRG_H1",   readMozer},
    {"LRGH1820", readSchatzel},
    {"LRGZEUS",  readZEUS},
};

TString glVarName = "";

//get best chi2 (theor normalization is free)
double getChi2(TString anal, TString var, vector<double> thVec)
{
    if(var.Contains('_')) //2D variables not implemented
        return 0;
    if(var.Contains("total")) //total not implemented
        return 0;

    ifstream matFile;
    matFile.open(TString::Format("../data/Measurement/fullErrors/data/%s.txt", anal.Data()));
    //search for variable var
    string line;
    bool isFound = false;
    while (matFile.is_open() &&  getline (matFile,line) ) {
        stringstream sline(line); 
        TString varCand;
        sline >> varCand;

        if(varCand == var) {
            isFound = true;
            break;
        }
    }
    if(!isFound) {
        cout << "Variable " << var<< " for " << anal << " not found" << endl;
        assert(isFound);
    }

    arma::mat errMat;
    arma::vec data;

    glVarName = var;
    //readVFPS(matFile, errMat, data);
    Reader.at(anal)(matFile, errMat, data);

    arma::vec th(thVec.size());
    for(unsigned i = 0; i < thVec.size(); ++i)
        th(i) = thVec[i];

    //cout << "Chi2      is " << getChi2(errMat, data, th) << endl;
    double nBest = 1;//getNorm(errMat, data, th);

    //cout << "Best norm is " << getNorm(errMat, data, th) << endl;

    /*
    for(double n : {nBest/1.01, nBest, nBest*1.01}) {
        cout << "Chi2      is " <<n <<" "<< getChi2(errMat, data, th*n) << endl;
    }
    */

    //for(double n = 1.01; n < 1.03; n+=0.001) {
    //}
    matFile.close();

    //Returning Best chi2
    return  getChi2(errMat, data, th*nBest) / (thVec.size()-1);
}



void readVFPS(ifstream &matFile, arma::mat &errMat, arma::vec &data)
{
    const int nMax = 20;
    arma::vec vals(nMax), errStat(nMax), errSys(nMax);
    arma::mat corStat(nMax,nMax), corSys(nMax, nMax);

    //read the info from the table
    int iNow = 0;
    string line;
    while (matFile.is_open() && getline (matFile,line) ) {
        if(line.size() < 10)
            break;

        double xmin, xmax;
        int id;
        stringstream linNow(line);
        //         0.014     0.019    2    2210    12.8    -14     14.4    -33     1.014    1.003
        linNow >> xmin >> xmax >> id >> vals(iNow) >> errStat(iNow);
        assert(iNow + 1 == id);

        for(int i = 0; i < iNow; ++i){
            linNow >> corStat(iNow,i);
        }

        linNow >> errSys(iNow);
        for(int i = 0; i < iNow; ++i){
            linNow >> corSys(iNow,i);
        }

        //cout << line << '\n';
        ++iNow;
    }
    assert(iNow >= 1);

    corSys.resize(iNow, iNow);
    corStat.resize(iNow, iNow);
    vals.resize(iNow);
    errStat.resize(iNow);
    errSys.resize(iNow);

    errSys  *= 0.01;
    errStat *= 0.01;
    corSys  *= 0.01;
    corStat *= 0.01;

    //matFile.close();

    //Add diagonal
    for(unsigned i = 0; i < corStat.n_rows; ++i) {
        corStat(i,i) = pow(errStat(i), 2);
        corSys(i,i)  = pow(errSys(i) , 2);
    }

    //calc off-diagonal -- erij / sqrt(erii * erjj)
    for(unsigned i = 0; i < corStat.n_rows; ++i) 
    for(unsigned j = 0; j < i; ++j) {
        corStat(i,j) = corStat(i,j) * errStat(i) * errStat(j);
        corStat(j,i) = corStat(i,j);

        corSys(i,j) = corSys(i,j) * errSys(i) * errSys(j);
        corSys(j,i) = corSys(i,j);
    }

    data = vals;
    errMat = corSys + corStat;

}

void readFPS(ifstream &matFile,  arma::mat &errMat, arma::vec &data)
{
    const int nMax = 20;
    arma::vec vals(nMax), errStat(nMax), errSys(nMax);
    arma::mat corStat(nMax,nMax);
    arma::mat shiftMat(nMax, 11);

    //read the info from the table
    int iNow = 0;
    string line;
    while (matFile.is_open() && getline (matFile,line) ) {
        if(line.size() < 10)
            break;

        double xmin, xmax, errTot, hadr, hadrErr;
        stringstream linNow(line);
        //         0.014     0.019    2    2210    12.8    -14     14.4    -33     1.014    1.003
        //$4 \div 6$ &$8.20$ &$13.2$ & $5.7$ & $11.9$ &$1.0$ & $4.5$ & $-3.9$ &$2.1$ & $1.1$ &$2.8$ & $-5.0$ &$1.9$ & $0.7$ &$1.2$ & $-0.9$ &$0.1$ &$1.05 \pm 0.05$ &$1.05$ 

        linNow >> xmin >> xmax >> vals(iNow) >> errTot >> errStat(iNow) >> errSys(iNow);

        //Read statistical correlations
        for(int i = 0; i < 5; ++i)
            linNow >> corStat(iNow, iNow +1+ i);


        for(unsigned i = 0; i < shiftMat.n_cols; ++i) { //read all shifts
           linNow >> shiftMat(iNow, i); 
        }
        linNow >> hadr >> hadrErr;

        //cout << line << '\n';
        ++iNow;
    }


    corStat.resize(iNow, iNow);
    vals.resize(iNow);
    errStat.resize(iNow);
    errSys.resize(iNow);
    shiftMat.resize(iNow, shiftMat.n_cols);

    errSys   *= 0.01;
    errStat  *= 0.01;
    shiftMat *= 0.01;

    //Calculate sysMatrix
    arma::mat corSys(iNow,iNow, arma::fill::zeros);
    for(unsigned i = 0; i < shiftMat.n_cols; ++i) {
        corSys += shiftMat.col(i) * shiftMat.col(i).t();
    }


    //matFile.close();

    //Add diagonal
    for(unsigned i = 0; i < corStat.n_rows; ++i) {
        corStat(i,i) = pow(errStat(i), 2);
    }

    //calc off-diagonal -- erij / sqrt(erii * erjj)
    for(unsigned i = 0;   i < corStat.n_rows; ++i) 
    for(unsigned j = i+1; j < corStat.n_cols; ++j) {
        corStat(i,j) = corStat(i,j) * errStat(i) * errStat(j);
        corStat(j,i) = corStat(i,j);
    }

    data = vals;
    errMat = corSys + corStat;




}

void readLRG(ifstream &matFile, arma::mat &errMat, arma::vec &data)
{
    const int nMax = 20;
    arma::vec vals(nMax), errStat(nMax), errSys(nMax);
    arma::mat corStat(nMax,nMax);
    arma::mat shiftMat(nMax, 12);

    //read the info from the table
    int iNow = 0;
    string line;
    bool hasCorr = false;
    while (matFile.is_open() && getline (matFile,line) ) {
        if(TString(line.c_str()).BeginsWith("corr")) {
            hasCorr = true;
            break;
        }

        double xmin, xmax, errTot, hadr, hadrErr, rad;
        stringstream linNow(line);
        //         0.014     0.019    2    2210    12.8    -14     14.4    -33     1.014    1.003
        //$4 \div 6$ &$8.20$ &$13.2$ & $5.7$ & $11.9$ &$1.0$ & $4.5$ & $-3.9$ &$2.1$ & $1.1$ &$2.8$ & $-5.0$ &$1.9$ & $0.7$ &$1.2$ & $-0.9$ &$0.1$ &$1.05 \pm 0.05$ &$1.05$ 

        linNow >> xmin >> xmax >> vals(iNow) >> errTot >> errStat(iNow) >> errSys(iNow);

        for(unsigned i = 0; i < shiftMat.n_cols; ++i) { //read all shifts
           linNow >> shiftMat(iNow, i); 
        }
        linNow >> hadr >> hadrErr >> rad;

        //cout << line << '\n';
        ++iNow;
    }
    assert(hasCorr);
    assert(iNow >= 1);

    int i = 0;
    while (matFile.is_open() && getline (matFile,line) ) {
        if(line.size() < 5)
            break;

        double xmin, xmax;
        int id;
        stringstream linNow(line);

        linNow >> xmin >> xmax >> id;
        assert(id - 1 == i);
        for(int j = 0; j < iNow; ++j) {
            linNow >> corStat(id-1, j);
        }
        ++i;
    }
    assert(iNow == i);


    corStat.resize(iNow, iNow);
    vals.resize(iNow);
    errStat.resize(iNow);
    errSys.resize(iNow);
    shiftMat.resize(iNow, shiftMat.n_cols);

    errSys   *= 0.01;
    errStat  *= 0.01;
    corStat  *= 0.01;
    shiftMat *= 0.01;

    //Calculate sysMatrix
    arma::mat corSys(iNow,iNow, arma::fill::zeros);
    for(unsigned i = 0; i < shiftMat.n_cols; ++i) {
        corSys += shiftMat.col(i) * shiftMat.col(i).t();
    }


    //matFile.close();

    //Add diagonal
    for(unsigned i = 0; i < corStat.n_rows; ++i) {
        corStat(i,i) = pow(errStat(i), 2);
    }

    //calc off-diagonal -- erij / sqrt(erii * erjj)
    for(unsigned i = 0; i < corStat.n_rows; ++i) 
    for(unsigned j = 0; j < i; ++j) {
        corStat(i,j) = corStat(i,j) * errStat(i) * errStat(j);
        corStat(j,i) = corStat(i,j);
    }

    data = vals;
    errMat = corSys + corStat;

}

void readMozer(ifstream &matFile, arma::mat &errMat, arma::vec &data)
{
    readNoMat("mozer", matFile, errMat, data);
}

void readSchatzel(ifstream &matFile, arma::mat &errMat, arma::vec &data)
{
    readNoMat("schatzel", matFile, errMat, data);
}

void readNoMat(TString type, ifstream &matFile, arma::mat &errMat, arma::vec &data)
{
    const int nMax = 20;
    arma::vec vals(nMax), errStat(nMax), errSysCorr(nMax), errSysUncorr(nMax);

    //read the info from the table
    int iNow = 0;
    string line;
    //bool hasCorr = false;
    while (matFile.is_open() && getline (matFile,line) ) {
        if(line.size() < 5)
            break;

        double xmin, xmax, errTot, hadr, hadrErr;
        stringstream linNow(line);

        if(type == "mozer") {
            linNow >> xmin >> xmax >> vals(iNow) >> errTot >> errStat(iNow) >> errSysUncorr(iNow) >> errSysCorr(iNow);
            linNow >> hadr >> hadrErr;
        }
        else if(type == "schatzel") {
            linNow >> xmin >> xmax >> vals(iNow) >> errStat(iNow) >> errSysCorr(iNow) >> errTot;
            linNow >> hadr >> hadrErr;
            errSysUncorr(iNow) = sqrt(pow(errTot,2) - pow(errSysCorr(iNow),2) -  pow(errStat(iNow),2) );
        }
        //cout << line << '\n';
        ++iNow;
    }
    assert(iNow >= 1);


    vals.resize(iNow);
    errStat.resize(iNow);
    errSysCorr.resize(iNow);
    errSysUncorr.resize(iNow);


    //Diagonal uncorrelated matrix (relative)
    arma::mat corUncor(iNow,iNow, arma::fill::zeros);
    for(unsigned i = 0; i < vals.n_rows; ++i) {
        corUncor(i,i) = (pow(errStat(i),2) + pow(errSysUncorr(i),2)) / pow(vals(i),2);
    }

    //Correlated error matrix (relative)
    errSysCorr /= vals; //corr shift in relative units
    arma::mat corCorr = errSysCorr * errSysCorr.t();

    data = vals;
    errMat = corUncor + corCorr;

}


void readZEUS(ifstream &matFile, arma::mat &errMat, arma::vec &data)
{
    const int nMax = 20;
    arma::vec vals(nMax), errStat(nMax), errSysCorr(nMax), errSysUncorr(nMax);
    arma::mat shiftMat(nMax, 3);

    //read the info from the table
    int iNow = 0;
    string line;
    //bool hasCorr = false;
    while (matFile.is_open() && getline (matFile,line) ) {
        if(line.size() < 5)
            break;

        double xmin, xmax, hadr;
        double errSyst1, errSyst2;
        stringstream linNow(line);

        //5   8   7.4    0.3    +0.3    -0.5      +0.5   -0.5     0.1               1.04719640806102
        linNow >> xmin >> xmax >> vals(iNow) >> errStat(iNow) >> errSyst1 >> errSyst2 >> shiftMat(iNow, 0) >> shiftMat(iNow, 1) >> shiftMat(iNow, 2);
        linNow >> hadr;

        if(glVarName == "xpom") {
            double C = log(xmax/xmin) / (xmax-xmin);
            vals(iNow) *= C;
            errStat(iNow) *= C;
            errSyst1 *= C;
            errSyst2 *= C;
            shiftMat(iNow, 0) *= C;
            shiftMat(iNow, 1) *= C;
            shiftMat(iNow, 2) *= C;
        }

        errSysUncorr(iNow) = hypot(errSyst1, errSyst2);

        //cout << line << '\n';
        ++iNow;
    }
    assert(iNow >= 1);


    vals.resize(iNow);
    errStat.resize(iNow);
    errSysCorr.resize(iNow);
    errSysUncorr.resize(iNow);
    shiftMat.resize(iNow, shiftMat.n_cols);

    //Diagonal uncorrelated matrix (relative)
    arma::mat corUncor(iNow,iNow, arma::fill::zeros);
    for(unsigned i = 0; i < vals.n_rows; ++i) {
        corUncor(i,i) = (pow(errStat(i),2) + pow(errSysUncorr(i),2)) / pow(vals(i),2);
    }

    //Correlated error matrix (relative)
    errSysCorr /= vals; //corr shift in relative units
    arma::mat corCorr = errSysCorr * errSysCorr.t();

    //Calculate sysMatrix
    arma::mat corSys(iNow,iNow, arma::fill::zeros);
    //JetEnergyScaleUp + Down + diffractiveSelErr
    corSys  = 0.5*(shiftMat.col(0) * shiftMat.col(0).t() + shiftMat.col(1) * shiftMat.col(1).t());
    corSys += shiftMat.col(2) * shiftMat.col(2).t();


    data = vals;
    errMat = corUncor + corSys;

}








//get chi2 using cov-matrix V
double getChi2(arma::mat V, arma::vec data, arma::vec theor)
{
    arma::vec sh = log(theor/data);
    arma::mat res = sh.t() * V.i() * sh;
    assert(res.n_cols == 1 && res.n_rows == 1);
    return res(0,0);
}

//Get the theor-normalization with min chi2
double getNorm(arma::mat V, arma::vec data, arma::vec theor)
{
    arma::vec uni(data.n_rows, arma::fill::ones);
    arma::vec sh = log(theor/data);

    arma::mat num = uni.t()*V.i()*sh + sh.t()*V.i()*uni;
    assert(num.n_cols == 1 && num.n_rows == 1);
    arma::mat den = 2.* uni.t()*V.i()*uni;
    assert(den.n_cols == 1 && den.n_rows == 1);

    return 1./exp(num(0,0) / den(0,0));
}



#if 0
int main()
{
    TString dataFile = "LRGZEUS";

    /*
    if(dataFile == "VFPS") {
        readMat(dataFile, "xpom", {2200, 2200, 2200});
    }
    else if(dataFile == "LRG") {
        readMat(dataFile, "q2", {8, 4, 2, 0.7, 0.1});
    }
    else if(dataFile == "FPS") {
        readMat(dataFile, "q2", {21, 10, 3, 1.1, 0.3});
    }
    else if(dataFile == "LRG_H1") {
        readMat(dataFile, "ptjet1", {21, 14, 7, 2, 0.4});
    }
    else if(dataFile == "LRGH1820") {
        readMat(dataFile, "zpom", {60, 30, 20, 5});
    }
    else if(dataFile == "LRGZEUS") {
        readMat(dataFile, "w", {0.4, 0.4, 0.7, 0.7, 0.8, 0.8});
    }
    */

    return 0;
}
#endif
