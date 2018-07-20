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

using namespace RooFit ;
using namespace std ;

int maxOrder = 15;

void test1D_product()
{
  // Declare observable x
  RooRealVar x("x","x",-TMath::Pi(),TMath::Pi()) ;

  // To construct a proper p.d.f, the formula expression is explicitly normalized internally by dividing 
  // it by a numeric integral of the expresssion over x in the range [-20,20] 
  //
  RooRealVar alpha("alpha","alpha",30,0.1,60) ;
  RooGenericPdf genpdf("genpdf","genpdf","(1+0.1*abs(x)+sin(sqrt(abs(x*alpha+0.1))))",RooArgSet(x,alpha)) ;

  // Generate a toy dataset from the interpreted p.d.f
  RooDataSet* data = genpdf.generate(x,1e5) ;

  // Fit the interpreted p.d.f to the generated data
  genpdf.fitTo(*data) ;

  // Make a plot of the data and the p.d.f overlaid
  RooPlot* xframe = x.frame(Title("Interpreted expression pdf")) ;
  data->plotOn(xframe) ;
  genpdf.plotOn(xframe) ;  

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

  RooGenericPdf projectedPdf ("projectedPdf","projectedPdf",sumExpression.c_str(),facList);

  double xVal;

  for (int iPoint=0; iPoint<data->numEntries(); ++iPoint) {
    
    xVal = data->get(iPoint)->getRealValue("x");

    cosproj[0] += 1;
    
    for (int iOrder=1; iOrder<=maxOrder; ++iOrder) {
      sinproj[iOrder] += TMath::Sin( iOrder*xVal );
      cosproj[iOrder] += TMath::Cos( iOrder*xVal );
    }

  }
      
  cosFac[0]->setVal(cosproj[0]/TMath::Pi()/2.);

  for (int iOrder=1; iOrder<=maxOrder; ++iOrder) {
    sinFac[iOrder]->setVal(sinproj[iOrder]/TMath::Pi());
    cosFac[iOrder]->setVal(cosproj[iOrder]/TMath::Pi());
  }

  // Make a plot of the data and the p.d.f overlaid
  RooPlot* xframe1 = x.frame(Title("Projected pdf")) ;
  data->plotOn(xframe1) ;
  genpdf.plotOn(xframe1) ;  
  projectedPdf.plotOn(xframe1,LineColor(kRed)) ;

  // Draw frame on a canvas
  TCanvas* c1 = new TCanvas("Can_Proj_Pdf","Projected pdf",800,1600) ;
  gPad->SetLeftMargin(0.15) ; xframe1->GetYaxis()->SetTitleOffset(1.4) ; xframe1->Draw() ;

}
