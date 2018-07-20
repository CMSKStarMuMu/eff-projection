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

int maxOrder =5;

void test1Deff_comp()
{

  RooRealVar x("x","x",-TMath::Pi(),TMath::Pi()) ;

  //Generate PDFs for denominator, efficiency and numerator
  RooRealVar alpha("alpha","alpha",30,0.1,60) ;
  RooGenericPdf genpdf("genpdf","genpdf","(1+0.1*abs(x)+0.2*sin(sqrt(abs(x*alpha+0.1))))",RooArgSet(x,alpha)) ;

  RooRealVar mean("mean","mean",0.5);
  RooRealVar sigma("sigma","sigma",1.5);
  RooGaussian effPdf("effPdf","effPdf",x,mean,sigma);

  RooProdPdf numPdf("numPdf","numPdf",genpdf,effPdf);

  // Generate toy datasets
  RooDataSet* data = genpdf.generate(x,1e5) ;
  RooDataSet* numData = numPdf.generate(x,2e4) ;
  double avgEff = numData->sumEntries() / data->sumEntries();
  cout<<"Average efficiency = "<<avgEff<<endl;

  // Plot numerator and denominator datasets
  RooPlot* xframe = x.frame(Title("Numerator and denominator toys distributions")) ;
  data->plotOn(xframe) ;
  numData->plotOn(xframe) ;  

  TCanvas* c = new TCanvas("NumDenCanvas","toy_numerator_denominator",800,800) ;
  gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;

  //Declare and initialise the PDFs and all the needed objects
  vector < RooRealVar* > sinFac;
  vector < RooRealVar* > cosFac;
  vector < double > sinproj;
  vector < double > cosproj;
  vector < RooGenericPdf* > partialPdf;

  vector < RooRealVar* > sinFacAs;
  vector < RooRealVar* > cosFacAs;
  vector < double > sinprojAs;
  vector < double > cosprojAs;
  vector < RooGenericPdf* > partialPdfAs;

  string sumExpression = "";
  string sumExpressionAs = "";
  RooArgList facList(x);
  RooArgList facListAs(x);
  
  for (int iOrder=0; iOrder<=maxOrder; ++iOrder) {

    sinFac.push_back( new RooRealVar(Form("s%i",iOrder),Form("s%i",iOrder),0) );
    cosFac.push_back( new RooRealVar(Form("c%i",iOrder),Form("c%i",iOrder),0) );
    sinFacAs.push_back( new RooRealVar(Form("sa%i",iOrder),Form("sa%i",iOrder),0) );
    cosFacAs.push_back( new RooRealVar(Form("ca%i",iOrder),Form("ca%i",iOrder),0) );

    sinproj.push_back(0);
    cosproj.push_back(0);
    sinprojAs.push_back(0);
    cosprojAs.push_back(0);

    if ( iOrder==0 ) {
      sumExpression = Form("c%i*cos(%i*x)",iOrder,iOrder);
      sumExpressionAs = Form("ca%i*cos(%i*x)",iOrder,iOrder);
      facList.add(*cosFac[iOrder]);
      facListAs.add(*cosFacAs[iOrder]);
    } else {
      sumExpression = sumExpression + Form(" + s%i*sin(%i*x) + c%i*cos(%i*x)",iOrder,iOrder,iOrder,iOrder);
      sumExpressionAs = sumExpressionAs + Form(" + sa%i*sin(%i*x) + ca%i*cos(%i*x)",iOrder,iOrder,iOrder,iOrder);
      facList.add(*sinFac[iOrder]);
      facList.add(*cosFac[iOrder]);
      facListAs.add(*sinFacAs[iOrder]);
      facListAs.add(*cosFacAs[iOrder]);
    }

    if (iOrder<maxOrder) {
      partialPdf.push_back( new RooGenericPdf(Form("partPdf%i",iOrder),Form("partPdf%i",iOrder),sumExpression.c_str(),facList) );
      partialPdfAs.push_back( new RooGenericPdf(Form("partPdfAs%i",iOrder),Form("partPdfAs%i",iOrder),sumExpressionAs.c_str(),facListAs) );
    }

  }

  RooGenericPdf projectedPdf ("projectedPdf","projectedPdf",sumExpression.c_str(),facList);
  RooGenericPdf projectedPdfAs ("projectedPdfAs","projectedPdfAs",sumExpressionAs.c_str(),facListAs);

  //Create binned histos (only for projection method)
  TH1F* denHist = (TH1F*)data->createHistogram("denHist",x,Binning(100,-TMath::Pi(),TMath::Pi()));
  TH1F* numHist = (TH1F*)numData->createHistogram("numHist",x,Binning(100,-TMath::Pi(),TMath::Pi()));
  denHist->Sumw2();
  numHist->Sumw2();

  //Compute and set the coefficients
  for (int iOrder=1; iOrder<=maxOrder; ++iOrder) {
    sinprojAs[iOrder] = numData->sumEntries(Form("sin(%i*x)>0",iOrder))/data->sumEntries(Form("sin(%i*x)>0",iOrder))
      - numData->sumEntries(Form("sin(%i*x)<0",iOrder))/data->sumEntries(Form("sin(%i*x)<0",iOrder));
    cosprojAs[iOrder] = numData->sumEntries(Form("cos(%i*x)>0",iOrder))/data->sumEntries(Form("cos(%i*x)>0",iOrder))
      - numData->sumEntries(Form("cos(%i*x)<0",iOrder))/data->sumEntries(Form("cos(%i*x)<0",iOrder));
    for (int iBin=1; iBin<=denHist->GetNbinsX(); ++iBin) {
      sinproj[iOrder] += TMath::Sin( iOrder*denHist->GetBinCenter(iBin) ) * numHist->GetBinContent(iBin)/denHist->GetBinContent(iBin) * denHist->GetBinWidth(iBin);
      cosproj[iOrder] += TMath::Cos( iOrder*denHist->GetBinCenter(iBin) ) * numHist->GetBinContent(iBin)/denHist->GetBinContent(iBin) * denHist->GetBinWidth(iBin);
    }
  }
      
  cosFac[0]->setVal(avgEff);
  cosFacAs[0]->setVal(avgEff);

  for (int iOrder=1; iOrder<=maxOrder; ++iOrder) {
    sinFac[iOrder]->setVal(sinproj[iOrder]/TMath::Pi());
    cosFac[iOrder]->setVal(cosproj[iOrder]/TMath::Pi());
    sinFacAs[iOrder]->setVal(sinprojAs[iOrder]*TMath::Pi()/4);
    cosFacAs[iOrder]->setVal(cosprojAs[iOrder]*TMath::Pi()/4);
  }

  //Plot the PDFs
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

  RooPlot* xframe1 = x.frame(Title("Projected function")) ;
  fakeData->plotOn(xframe1,Invisible()) ;
  effPdf.plotOn(xframe1) ;  
  projectedPdf.plotOn(xframe1,LineColor(kRed)) ;
  projectedPdfAs.plotOn(xframe1,LineColor(kBlack)) ;

  TCanvas* c1 = new TCanvas("Can_Proj_Func","Projected function",1600,800) ;
  gPad->SetLeftMargin(0.15) ; xframe1->GetYaxis()->SetTitleOffset(1.4) ; 
  effHist->Draw();
  xframe1->Draw("same") ;

  //Kolmogorov-Smirnov test
  auto referTest = (effPdf.generate(x,1e5))->createHistogram("referTest",x,Binning(1e6,-TMath::Pi(),TMath::Pi()));
  auto projTest = (projectedPdf.generate(x,1e5))->createHistogram("projTest",x,Binning(1e6,-TMath::Pi(),TMath::Pi()));
  auto asymTest = (projectedPdfAs.generate(x,1e5))->createHistogram("asymTest",x,Binning(1e6,-TMath::Pi(),TMath::Pi()));
  double KStestProj = projTest->KolmogorovTest(referTest,"M");
  double KStestAsym = asymTest->KolmogorovTest(referTest,"M");
  cout<<"KS test result - projected PDF: "<<KStestProj<<endl;
  cout<<"KS test result - asymmetry PDF: "<<KStestAsym<<endl;
  TGraph* partKS = new TGraph();
  TGraph* partKSAs = new TGraph();
  partKS->SetName("partKS");
  partKSAs->SetName("partKSAs");
  for (int iPDF=1; iPDF<partialPdfAs.size(); ++iPDF) {
    auto partTest = (partialPdf[iPDF]->generate(x,1e5))->createHistogram(Form("part%iTest",iPDF),x,Binning(1e6,-TMath::Pi(),TMath::Pi()));
    auto partTestAs = (partialPdfAs[iPDF]->generate(x,1e5))->createHistogram(Form("part%iTestAs",iPDF),x,Binning(1e6,-TMath::Pi(),TMath::Pi()));
    partKS->SetPoint(iPDF-1,iPDF,partTest->KolmogorovTest(referTest,"M"));
    partKSAs->SetPoint(iPDF-1,iPDF,partTestAs->KolmogorovTest(referTest,"M"));
  }
  partKS->SetPoint(partialPdfAs.size()-1,partialPdfAs.size(),KStestProj);
  partKSAs->SetPoint(partialPdfAs.size()-1,partialPdfAs.size(),KStestAsym);
  partKS->SetMarkerStyle(20);
  partKS->SetMarkerColor(2);
  partKSAs->SetMarkerStyle(20);
  TCanvas* c2 = new TCanvas("Can_KSTest","KS test",800,800);
  partKS->Draw("AP");
  partKSAs->Draw("P");
  

}
