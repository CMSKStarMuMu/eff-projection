#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooCFunction1Binding.h"

using namespace RooFit ;
using namespace std ;

int maxOrder = 50;

void test1Deff_product()
{
  // Declare observable x
  RooRealVar x("x","x",-TMath::Pi(),TMath::Pi()) ;

  // To construct a proper p.d.f, the formula expression is explicitly normalized internally by dividing 
  // it by a numeric integral of the expresssion over x in the range [-20,20] 
  //
  RooRealVar alpha("alpha","alpha",30,0.1,60) ;
  RooGenericPdf genpdf("genpdf","genpdf","(1+0.1*abs(x)+sin(sqrt(abs(x*alpha+0.1))))",RooArgSet(x,alpha)) ;

  RooRealVar mean("mean","mean",0.5);
  RooRealVar sigma("sigma","sigma",1.5);
  RooGaussian effPdf("effPdf","effPdf",x,mean,sigma);
  RooProdPdf numPdf("numPdf","numPdf",genpdf,effPdf);

  // Generate a toy dataset from the interpreted p.d.f
  RooDataSet* data = genpdf.generate(x,1e5) ;
  RooDataSet* numData = numPdf.generate(x,2e4) ;
  double avgEff = numData->sumEntries() / data->sumEntries();
  cout<<"Average efficiency = "<<avgEff<<endl;

  // Make a plot of the data and the p.d.f overlaid
  RooPlot* xframe = x.frame(Title("Numerator and denominator distributions")) ;
  data->plotOn(xframe) ;
  numData->plotOn(xframe) ;  

  // Draw frame on a canvas
  TCanvas* c = new TCanvas("rf103_interprfuncs","rf103_interprfuncs",800,800) ;
  gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;

  vector < RooRealVar* > sinFac;
  vector < RooRealVar* > cosFac;
  vector < double > sinproj;
  vector < double > cosproj;

  string sumExpression = "";
  RooArgList facList(x);
  
  for (int iOrder=0; iOrder<=maxOrder; ++iOrder) {

    sinFac.push_back( new RooRealVar(Form("s%i",iOrder),Form("s%i",iOrder),0) );
    cosFac.push_back( new RooRealVar(Form("c%i",iOrder),Form("c%i",iOrder),0) );

    sinproj.push_back(0);
    cosproj.push_back(0);

    if ( iOrder==0 ) {
      sumExpression = Form("c%i*cos(%i*x)",iOrder,iOrder);
      facList.add(*cosFac[iOrder]);
    } else {
      sumExpression = sumExpression + Form(" + s%i*sin(%i*x) + c%i*cos(%i*x)",iOrder,iOrder,iOrder,iOrder);
      facList.add(*sinFac[iOrder]);
      facList.add(*cosFac[iOrder]);
    }

  }

  // TF1* tFunc = new TF1("tFunc",sumExpression.c_str(),-TMath::Pi(),TMath::Pi());
  // RooAbsReal* projectedFunc = bindFunction(tFunc,facList);
  // projectedFunc->SetName("projectedFunc");
  RooGenericPdf projectedPdf ("projectedPdf","projectedPdf",sumExpression.c_str(),facList);

  TH1F* denHist = (TH1F*)data->createHistogram("denHist",x,Binning(100,-TMath::Pi(),TMath::Pi()));
  TH1F* numHist = (TH1F*)numData->createHistogram("numHist",x,Binning(100,-TMath::Pi(),TMath::Pi()));
  denHist->Sumw2();
  numHist->Sumw2();

  for (int iBin=1; iBin<=denHist->GetNbinsX(); ++iBin) {

    for (int iOrder=1; iOrder<=maxOrder; ++iOrder) {
      sinproj[iOrder] += TMath::Sin( iOrder*denHist->GetBinCenter(iBin) ) * numHist->GetBinContent(iBin)/denHist->GetBinContent(iBin) * denHist->GetBinWidth(iBin);
      cosproj[iOrder] += TMath::Cos( iOrder*denHist->GetBinCenter(iBin) ) * numHist->GetBinContent(iBin)/denHist->GetBinContent(iBin) * denHist->GetBinWidth(iBin);
    }

  }
  // for (int iPoint=0; iPoint<data->numEntries(); ++iPoint) {
    
  //   xVal = data->get(iPoint)->getRealValue("x");

  //   cosproj[0] += 1;
    
  //   for (int iOrder=1; iOrder<=maxOrder; ++iOrder) {
  //     sinproj[iOrder] += TMath::Sin( iOrder*xVal );
  //     cosproj[iOrder] += TMath::Cos( iOrder*xVal );
  //   }

  // }
      
  cosFac[0]->setVal(avgEff);

  for (int iOrder=1; iOrder<=maxOrder; ++iOrder) {
    sinFac[iOrder]->setVal(sinproj[iOrder]/TMath::Pi());
    cosFac[iOrder]->setVal(cosproj[iOrder]/TMath::Pi());
  }

  cout<<"Constant coefficients: "<<cosFac[0]->getValV()<<endl;
  cout<<"First order sin coefficients: "<<sinproj[1]<<" -> "<<sinFac[1]->getValV()<<endl;

  TH1F* altEffHist = (TH1F*)numHist->Clone("altEffHist");
  altEffHist->Sumw2();
  altEffHist->Divide(denHist);
  altEffHist->SetLineColor(2);
  double scaleFac = 20. / altEffHist->Integral();
  numHist->Scale( scaleFac );
  TEfficiency* effHist = new TEfficiency(*numHist,*denHist);
  effHist->SetName("effHist");
  effHist->SetTitle(Form("Toy efficiency;;Efficiency * %f",scaleFac));
  RooDataSet* fakeData = new RooDataSet("fakeData","fakeData",x);
  for (int i=0; i<20; ++i) fakeData->add(RooArgSet(x),0.5);

  // Make a plot of the data and the p.d.f overlaid
  RooPlot* xframe1 = x.frame(Title("Projected function")) ;
  fakeData->plotOn(xframe1,Invisible()) ;
  effPdf.plotOn(xframe1) ;  
  projectedPdf.plotOn(xframe1,LineColor(kRed)) ;

  // Draw frame on a canvas
  TCanvas* c1 = new TCanvas("Can_Proj_Func","Projected function",1600,800) ;
  gPad->SetLeftMargin(0.15) ; xframe1->GetYaxis()->SetTitleOffset(1.4) ; 
  effHist->Draw();
  // altEffHist->Draw("EPsame") ;
  xframe1->Draw("same") ;

  

}
