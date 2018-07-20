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

int maxOrder =8;
int nTot = 1e5;

void test2Deff_comp_alt2()
{

  RooRealVar x("x","x",-TMath::Pi(),TMath::Pi());
  RooRealVar y("y","y",-TMath::Pi(),TMath::Pi());

  //Generate PDFs for denominator, efficiency and numerator
  // RooGenericPdf genpdf("genpdf","genpdf","4+(x^2-2)-(y^2/3-2)",RooArgSet(x,y)) ;
  RooGenericPdf genpdf("genpdf","genpdf","4+(0.8*x)-(0.2*y)",RooArgSet(x,y)) ;

  RooRealVar mean("mean","mean",0.);
  RooRealVar sigma("sigma","sigma",1.2);
  // RooGaussian effYPdf("effYPdf","effYPdf",y,mean,sigma);
  // RooGenericPdf effXPdf("effXPdf","effXPdf","5-2/3.14*(x+abs(x))",x);
  // RooGenericPdf effYPdf("effYPdf","effYPdf","3+sin(y)",y);
  // RooGenericPdf effXPdf("effXPdf","effXPdf","3+cos(x)",x);
  // RooProdPdf effPdf("effPdf","effPdf",effXPdf,effYPdf);

  // RooGenericPdf effPdf("effPdf","effPdf","3+0.1*cos(x)*sin(y)-sin(x)+cos(4*x)*cos(y)",RooArgSet(x,y));
  RooGenericPdf effPdf("effPdf","effPdf","(3-y)*4*exp(-1.0/8*x^2)+(3+y)*(x^2/4+1)",RooArgSet(x,y));

  RooProdPdf numPdf("numPdf","numPdf",genpdf,effPdf);

  // Generate toy datasets
  RooDataSet* data = genpdf.generate(RooArgSet(x,y),nTot) ;
  RooDataSet* numData = numPdf.generate(RooArgSet(x,y),nTot/5) ;
  double avgEff = numData->sumEntries() / data->sumEntries();
  cout<<"Average efficiency = "<<avgEff<<endl;

  // Plot numerator and denominator datasets
  RooPlot* xframe = x.frame(Title("Numerator and denominator x distributions")) ;
  data->plotOn(xframe) ;
  numData->plotOn(xframe) ;  
  RooPlot* yframe = y.frame(Title("Numerator and denominator y distributions")) ;
  data->plotOn(yframe) ;
  numData->plotOn(yframe) ;  

  TCanvas* c = new TCanvas("NumDenCanvas","toy_numerator_denominator",800,800) ;
  c->Divide(2,1);
  c->cd(1);
  gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;
  c->cd(2);
  gPad->SetLeftMargin(0.15) ; yframe->GetYaxis()->SetTitleOffset(1.4) ; yframe->Draw() ;

  //Declare and initialise the PDFs and all the needed objects
  vector < RooRealVar* > factors;
  vector < double > proj;
  vector < RooFormulaVar* > vectPdf;

  vector < RooRealVar* > factorsAs;
  vector < double > projAs;
  vector < RooFormulaVar* > vectPdfAs;

  RooArgList facList;
  RooArgList facListAs;
  RooArgList pdfList;
  RooArgList pdfListAs;
    
  for (int xOrder=0; xOrder<=maxOrder; ++xOrder) for (int yOrder=0; yOrder<=maxOrder; ++yOrder) {

	int iOrder = (yOrder + xOrder*(maxOrder+1))*4;

	factors.push_back( new RooRealVar(Form("s%is%i",xOrder,yOrder),Form("s%is%i",xOrder,yOrder),0) );
	factors.push_back( new RooRealVar(Form("s%ic%i",xOrder,yOrder),Form("s%ic%i",xOrder,yOrder),0) );
	factors.push_back( new RooRealVar(Form("c%is%i",xOrder,yOrder),Form("c%is%i",xOrder,yOrder),0) );
	factors.push_back( new RooRealVar(Form("c%ic%i",xOrder,yOrder),Form("c%ic%i",xOrder,yOrder),0) );

	factorsAs.push_back( new RooRealVar(Form("sa%isa%i",xOrder,yOrder),Form("sa%isa%i",xOrder,yOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("sa%ica%i",xOrder,yOrder),Form("sa%ica%i",xOrder,yOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("ca%isa%i",xOrder,yOrder),Form("ca%isa%i",xOrder,yOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("ca%ica%i",xOrder,yOrder),Form("ca%ica%i",xOrder,yOrder),0) );

	vectPdf  .push_back( new RooFormulaVar( Form("pdf_s%i_s%i",xOrder,yOrder),   Form("pdf_s%i_s%i",xOrder,yOrder), 
						Form("sin(%i*x)*sin(%i*y)",xOrder,yOrder), RooArgList(x,y) ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_s%i_s%i",xOrder,yOrder), Form("pdfas_s%i_s%i",xOrder,yOrder), 
						Form("sin(%i*x)*sin(%i*y)",xOrder,yOrder), RooArgList(x,y) ) );
	vectPdf  .push_back( new RooFormulaVar( Form("pdf_s%i_c%i",xOrder,yOrder),   Form("pdf_s%i_c%i",xOrder,yOrder), 
						Form("sin(%i*x)*cos(%i*y)",xOrder,yOrder), RooArgList(x,y) ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_s%i_c%i",xOrder,yOrder), Form("pdfas_s%i_c%i",xOrder,yOrder), 
						Form("sin(%i*x)*cos(%i*y)",xOrder,yOrder), RooArgList(x,y) ) );
	vectPdf  .push_back( new RooFormulaVar( Form("pdf_c%i_s%i",xOrder,yOrder),   Form("pdf_c%i_s%i",xOrder,yOrder), 
						Form("cos(%i*x)*sin(%i*y)",xOrder,yOrder), RooArgList(x,y) ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_c%i_s%i",xOrder,yOrder), Form("pdfas_c%i_s%i",xOrder,yOrder), 
						Form("cos(%i*x)*sin(%i*y)",xOrder,yOrder), RooArgList(x,y) ) );
	vectPdf  .push_back( new RooFormulaVar( Form("pdf_c%i_c%i",xOrder,yOrder),   Form("pdf_c%i_c%i",xOrder,yOrder), 
						Form("cos(%i*x)*cos(%i*y)",xOrder,yOrder), RooArgList(x,y) ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_c%i_c%i",xOrder,yOrder), Form("pdfas_c%i_c%i",xOrder,yOrder), 
						Form("cos(%i*x)*cos(%i*y)",xOrder,yOrder), RooArgList(x,y) ) );

	for (int i=0; i<4; ++i) {
	  proj.push_back(0);
	  projAs.push_back(0);
	}
	
	if ( xOrder>0 && yOrder>0 ) {
	  pdfList  .add(*vectPdf  [iOrder+0]);
	  pdfListAs.add(*vectPdfAs[iOrder+0]);
	  facList  .add(*factors  [iOrder+0]);
	  facListAs.add(*factorsAs[iOrder+0]);
	}
	if ( xOrder>0 ) {
	  pdfList  .add(*vectPdf  [iOrder+1]);
	  pdfListAs.add(*vectPdfAs[iOrder+1]);
	  facList  .add(*factors  [iOrder+1]);
	  facListAs.add(*factorsAs[iOrder+1]);
	}
	if ( yOrder>0 ) {
	  pdfList  .add(*vectPdf  [iOrder+2]);
	  pdfListAs.add(*vectPdfAs[iOrder+2]);
	  facList  .add(*factors  [iOrder+2]);
	  facListAs.add(*factorsAs[iOrder+2]);
	}
	{
	  pdfList  .add(*vectPdf  [iOrder+3]);
	  pdfListAs.add(*vectPdfAs[iOrder+3]);
	  facList  .add(*factors  [iOrder+3]);
	  facListAs.add(*factorsAs[iOrder+3]);
	}

	// if (iOrder<maxOrder) {
	//   partialPdf.push_back( new RooGenericPdf(Form("partPdf%i",iOrder),Form("partPdf%i",iOrder),sumExpression.c_str(),facList) );
	//   partialPdfAs.push_back( new RooGenericPdf(Form("partPdfAs%i",iOrder),Form("partPdfAs%i",iOrder),sumExpressionAs.c_str(),facListAs) );
	// }

      }

  RooAddition projectedPdf   ("projectedPdf"  , "projectedPdf"  , pdfList  , facList  );
  RooAddition projectedPdfAs ("projectedPdfAs", "projectedPdfAs", pdfListAs, facListAs);
  // RooRealSumPdf projectedPdf   ("projectedPdf"  , "projectedPdf"  , pdfList  , facList  );
  // RooRealSumPdf projectedPdfAs ("projectedPdfAs", "projectedPdfAs", pdfListAs, facListAs);

  //Create binned histos (only for projection method)
  TH2F* denHist = (TH2F*)data   ->createHistogram( "denHist",
						   x,     Binning(100,-TMath::Pi(),TMath::Pi()),
						   YVar(y,Binning(100,-TMath::Pi(),TMath::Pi())) );
  TH2F* numHist = (TH2F*)numData->createHistogram( "numHist",
						   x,     Binning(100,-TMath::Pi(),TMath::Pi()),
						   YVar(y,Binning(100,-TMath::Pi(),TMath::Pi())) );
  denHist->Sumw2();
  numHist->Sumw2();

  //Compute and set the coefficients
  factors  [3]->setVal(avgEff);
  factorsAs[3]->setVal(avgEff);

  double xCent, yCent, fact;

  for (int xOrder=0; xOrder<=maxOrder; ++xOrder) for (int yOrder=0; yOrder<=maxOrder; ++yOrder) {
	
	int iOrder = (yOrder + xOrder*(maxOrder+1))*4;

	projAs[iOrder+0] =
	  ( numData->sumEntries(Form("sin(%i*x)*sin(%i*y)>0",xOrder,yOrder)) /
	    data   ->sumEntries(Form("sin(%i*x)*sin(%i*y)>0",xOrder,yOrder)) ) -
	  ( numData->sumEntries(Form("sin(%i*x)*sin(%i*y)<0",xOrder,yOrder)) / 
	    data   ->sumEntries(Form("sin(%i*x)*sin(%i*y)<0",xOrder,yOrder)) );
	projAs[iOrder+1] =
	  ( numData->sumEntries(Form("sin(%i*x)*cos(%i*y)>0",xOrder,yOrder)) /
	    data   ->sumEntries(Form("sin(%i*x)*cos(%i*y)>0",xOrder,yOrder)) ) -
	  ( numData->sumEntries(Form("sin(%i*x)*cos(%i*y)<0",xOrder,yOrder)) / 
	    data   ->sumEntries(Form("sin(%i*x)*cos(%i*y)<0",xOrder,yOrder)) );
	projAs[iOrder+2] =
	  ( numData->sumEntries(Form("cos(%i*x)*sin(%i*y)>0",xOrder,yOrder)) /
	    data   ->sumEntries(Form("cos(%i*x)*sin(%i*y)>0",xOrder,yOrder)) ) -
	  ( numData->sumEntries(Form("cos(%i*x)*sin(%i*y)<0",xOrder,yOrder)) / 
	    data   ->sumEntries(Form("cos(%i*x)*sin(%i*y)<0",xOrder,yOrder)) );
	projAs[iOrder+3] =
	  ( numData->sumEntries(Form("cos(%i*x)*cos(%i*y)>0",xOrder,yOrder)) /
	    data   ->sumEntries(Form("cos(%i*x)*cos(%i*y)>0",xOrder,yOrder)) ) -
	  ( numData->sumEntries(Form("cos(%i*x)*cos(%i*y)<0",xOrder,yOrder)) / 
	    data   ->sumEntries(Form("cos(%i*x)*cos(%i*y)<0",xOrder,yOrder)) );

	// if (iOrder==40) cout<<Form("cos(%i*x)*cos(%i*y)*sin(%i*z)>0",xOrder,yOrder,zOrder) << " " << projAs[iOrder+6] << " "
	// 		    <<numData->sumEntries(Form("cos(%i*x)*cos(%i*y)*sin(%i*z)>0",xOrder,yOrder,zOrder))<< " "
	// 		    <<data   ->sumEntries(Form("cos(%i*x)*cos(%i*y)*sin(%i*z)>0",xOrder,yOrder,zOrder))<< " "
	// 		    <<numData->sumEntries(Form("cos(%i*x)*cos(%i*y)*sin(%i*z)<0",xOrder,yOrder,zOrder))<< " "
	// 		    <<data   ->sumEntries(Form("cos(%i*x)*cos(%i*y)*sin(%i*z)<0",xOrder,yOrder,zOrder))<<endl;

	for (int xBin=1;xBin<=denHist->GetNbinsX();++xBin) for (int yBin=1;yBin<=denHist->GetNbinsY();++yBin) {
	      xCent = denHist->GetXaxis()->GetBinCenter(xBin);
	      yCent = denHist->GetYaxis()->GetBinCenter(yBin);
	      if (denHist->GetBinContent(xBin,yBin)>0 )
		fact = numHist->GetBinContent(xBin,yBin) / denHist->GetBinContent(xBin,yBin) *
		  denHist->GetXaxis()->GetBinWidth(xBin) *
		  denHist->GetYaxis()->GetBinWidth(yBin);
	      else fact=0;

	      if (xOrder==0) fact = fact/2.;
	      if (yOrder==0) fact = fact/2.;

	      proj[iOrder+0] += TMath::Sin(xOrder*xCent) * TMath::Sin(yOrder*yCent) * fact;
	      proj[iOrder+1] += TMath::Sin(xOrder*xCent) * TMath::Cos(yOrder*yCent) * fact;
	      proj[iOrder+2] += TMath::Cos(xOrder*xCent) * TMath::Sin(yOrder*yCent) * fact;
	      proj[iOrder+3] += TMath::Cos(xOrder*xCent) * TMath::Cos(yOrder*yCent) * fact;
	    }

	if ( iOrder>0 ) for (int i=0; i<4; ++i) {
	    factors  [iOrder+i]->setVal(proj  [iOrder+i]/TMath::Pi()/TMath::Pi());
	    factorsAs[iOrder+i]->setVal(projAs[iOrder+i]*TMath::Pi()/4);
	  }
	cout<<iOrder<<" ss\t"<<proj[iOrder+0]/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+0]*TMath::Pi()/4<<endl;
	cout<<iOrder<<" sc\t"<<proj[iOrder+1]/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+1]*TMath::Pi()/4<<endl;
	cout<<iOrder<<" cs\t"<<proj[iOrder+2]/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+2]*TMath::Pi()/4<<endl;
	cout<<iOrder<<" cc\t"<<proj[iOrder+3]/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+3]*TMath::Pi()/4<<endl;
  
      }

  // for (int xOrder=0; xOrder<=maxOrder; ++xOrder) for (int yOrder=0; yOrder<=maxOrder; ++yOrder) for (int zOrder=0; zOrder<=maxOrder; ++zOrder) {
	
  // 	int iOrder = (zOrder + yOrder*(maxOrder+1) + xOrder*(maxOrder+1)*(maxOrder+1))*8;
	
  // 	if ( xOrder>0 && yOrder>0 && zOrder>0 ) 
  // 	  cout<<iOrder+0<<"\t"<<factors  [iOrder+0]->getValV()<<"\t"<<factorsAs[iOrder+0]->getValV()<<endl;
  // 	if ( xOrder>0 && yOrder>0 ) 
  // 	  cout<<iOrder+1<<"\t"<<factors  [iOrder+1]->getValV()<<"\t"<<factorsAs[iOrder+1]->getValV()<<endl;
  // 	if ( xOrder>0 && zOrder>0 ) 
  // 	  cout<<iOrder+2<<"\t"<<factors  [iOrder+2]->getValV()<<"\t"<<factorsAs[iOrder+2]->getValV()<<endl;
  // 	if ( xOrder>0 ) 
  // 	  cout<<iOrder+3<<"\t"<<factors  [iOrder+3]->getValV()<<"\t"<<factorsAs[iOrder+3]->getValV()<<endl;
  // 	if ( yOrder>0 && zOrder>0 ) 
  // 	  cout<<iOrder+4<<"\t"<<factors  [iOrder+4]->getValV()<<"\t"<<factorsAs[iOrder+4]->getValV()<<endl;
  // 	if ( yOrder>0 ) 
  // 	  cout<<iOrder+5<<"\t"<<factors  [iOrder+5]->getValV()<<"\t"<<factorsAs[iOrder+5]->getValV()<<endl;
  // 	if ( zOrder>0 ) 
  // 	  cout<<iOrder+6<<"\t"<<factors  [iOrder+6]->getValV()<<"\t"<<factorsAs[iOrder+6]->getValV()<<endl;
  // 	  cout<<iOrder+7<<"\t"<<factors  [iOrder+7]->getValV()<<"\t"<<factorsAs[iOrder+7]->getValV()<<endl;

  //     }

  //Plot the efficiency slices
  RooDataSet* fakeData = new RooDataSet("fakeData","fakeData",RooArgSet(x,y));
  for (int i=0; i<126; ++i) fakeData->add(RooArgSet(x,y));

  TCanvas* cx1 = new TCanvas("CanX_Proj_Func","Projected function - x projection",3000,600) ;
  TCanvas* cy1 = new TCanvas("CanY_Proj_Func","Projected function - y projection",3000,600) ;
  cx1->Divide(5,1);
  cy1->Divide(5,1);

  double xborder = 0.3;
  double yborder = 0.3;

  vector <TEfficiency*> effHistsX; 
  vector <TEfficiency*> effHistsY;
  vector <RooPlot*> xframes;
  vector <RooPlot*> yframes;
  // vector <TLegend*> xlegs;
  // vector <TLegend*> ylegs;

  TLegend* leg = new TLegend (0.5,0.7,0.9,0.9);

  for (int i=0; i<5; ++i) {

    double centX = -3 + 1.5*i;
    double centY = -3 + 1.5*i;
    double lowX  = TMath::Max( centX - xborder,  1e-4-TMath::Pi() );
    double lowY  = TMath::Max( centY - yborder,  1e-4-TMath::Pi() );
    double highX = TMath::Min( centX + xborder, -1e-4+TMath::Pi() );
    double highY = TMath::Min( centY + yborder, -1e-4+TMath::Pi() );
    
    auto numProjX = numHist->ProjectionX("numProjX", numHist->GetYaxis()->FindBin(lowY), numHist->GetYaxis()->FindBin(highY),"e");
    auto numProjY = numHist->ProjectionY("numProjY", numHist->GetXaxis()->FindBin(lowX), numHist->GetXaxis()->FindBin(highX),"e");
    auto denProjX = denHist->ProjectionX("denProjX", denHist->GetYaxis()->FindBin(lowY), denHist->GetYaxis()->FindBin(highY),"e");
    auto denProjY = denHist->ProjectionY("denProjY", denHist->GetXaxis()->FindBin(lowX), denHist->GetXaxis()->FindBin(highX),"e");
    
    effHistsX.push_back( new TEfficiency(*numProjX,*denProjX) );
    effHistsX[i]->SetName( Form("effHistX_%i",i) );
    effHistsX[i]->SetTitle(Form("Toy efficiency - slice y=%1.1f;x;Efficiency",centX) );
    
    effHistsY.push_back( new TEfficiency(*numProjY,*denProjY) );
    effHistsY[i]->SetName( Form("effHistY_%i",i) );
    effHistsY[i]->SetTitle(Form("Toy efficiency - slice x=%1.1f;y;Efficiency",centY) );

    x.setVal(centX);
    y.setVal(centY);
    xframes.push_back( x.frame(Title( Form("Projected function - x projection %i",i))) );
    fakeData->     plotOn(xframes[i],Invisible()) ;
    effPdf.        plotOn(xframes[i],Slice(RooArgSet(y))                  ,Name(Form("effPdf_%i"        ,i))) ;  
    projectedPdf.  plotOn(xframes[i],LineColor(kRed),Slice(RooArgSet(y))  ,Name(Form("projectedPdf_%i"  ,i))) ;
    projectedPdfAs.plotOn(xframes[i],LineColor(kBlack),Slice(RooArgSet(y)),Name(Form("projectedPdfAs_%i",i))) ;
    yframes.push_back( y.frame(Title( Form("Projected function - y projection %i",i))) );
    fakeData->     plotOn(yframes[i],Invisible()) ;
    effPdf.        plotOn(yframes[i],Slice(RooArgSet(x)));
    projectedPdf.  plotOn(yframes[i],LineColor(kRed),Slice(RooArgSet(x)));
    projectedPdfAs.plotOn(yframes[i],LineColor(kBlack),Slice(RooArgSet(x)));

    cx1->cd(i+1);
    // gPad->SetLeftMargin(0.15) ; xframes[i]->GetYaxis()->SetTitleOffset(1.4) ; 
    effHistsX[i]->Draw();
    cx1->cd(i+1)->Update(); 
    auto graphx = effHistsX[i]->GetPaintedGraph(); 
    graphx->SetMinimum(0);
    graphx->SetMaximum(1);
    graphx->GetYaxis()->SetTitleOffset(1.4);
    cx1->cd(i+1)->Update();
    xframes[i]->Draw("same") ;

    if (i==0) {
      leg->AddEntry(effHistsX[i]               ,"Binned efficiency","lep");
      leg->AddEntry(Form("effPdf_%i"        ,i),"GEN efficiency","l");
      leg->AddEntry(Form("projectedPdf_%i"  ,i),"Standard projection","l");
      leg->AddEntry(Form("projectedPdfAs_%i",i),"Asymmetry projection","l");
    }

    leg->Draw("same");

    cy1->cd(i+1);
    // gPad->SetLeftMargin(0.15) ; yframes[i]->GetYaxis()->SetTitleOffset(1.4) ; 
    effHistsY[i]->Draw();
    cy1->cd(i+1)->Update(); 
    auto graphy = effHistsY[i]->GetPaintedGraph(); 
    graphy->SetMinimum(0);
    graphy->SetMaximum(1); 
    graphy->GetYaxis()->SetTitleOffset(1.4);
    cy1->cd(i+1)->Update();
    yframes[i]->Draw("same") ;
    leg->Draw("same");

  }    
    
}
