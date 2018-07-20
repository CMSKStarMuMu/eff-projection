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

int maxOrder =3;
int nTot = 1e9;

void test2Deff_comp_alt()
{

  RooRealVar x("x","x",-TMath::Pi(),TMath::Pi());
  RooRealVar y("y","y",-TMath::Pi(),TMath::Pi());

  //Generate PDFs for denominator, efficiency and numerator
  RooGenericPdf genpdf("genpdf","genpdf","4+(x^2-2)-(y^2/3-2)",RooArgSet(x,y)) ;

  RooRealVar mean("mean","mean",0.);
  RooRealVar sigma("sigma","sigma",1.2);
  // RooGaussian effYPdf("effYPdf","effYPdf",y,mean,sigma);
  // RooGenericPdf effXPdf("effXPdf","effXPdf","5-2/3.14*(x+abs(x))",x);
  RooGenericPdf effYPdf("effYPdf","effYPdf","3+sin(y)",y);
  RooGenericPdf effXPdf("effXPdf","effXPdf","3+cos(x)",x);
  RooProdPdf effPdf("effPdf","effPdf",effXPdf,effYPdf);

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

  //Plot the PDFs
  double xborder = 0.5;
  double yborder = 0.5;
  auto numProjX = numHist->ProjectionX("numProjX",
				       numHist->GetYaxis()->FindBin(-1*yborder),
				       numHist->GetYaxis()->FindBin(yborder),"e");
  auto numProjY = numHist->ProjectionY("numProjY",
				       numHist->GetXaxis()->FindBin(-1*xborder),
				       numHist->GetXaxis()->FindBin(xborder),"e");
  auto denProjX = denHist->ProjectionX("denProjX",
				       denHist->GetYaxis()->FindBin(-1*yborder),
				       denHist->GetYaxis()->FindBin(yborder),"e");
  auto denProjY = denHist->ProjectionY("denProjY",
				       denHist->GetXaxis()->FindBin(-1*xborder),
				       denHist->GetXaxis()->FindBin(xborder),"e");
  auto altEffHistX = (TH1F*)numProjX->Clone("altEffHistX");
  auto altEffHistY = (TH1F*)numProjY->Clone("altEffHistY");
  altEffHistX->Sumw2();
  altEffHistY->Sumw2();
  altEffHistX->Divide(denProjX);
  altEffHistY->Divide(denProjY);
  double scaleFacX = 20. / altEffHistX->Integral() /3.14;
  double scaleFacY = 20. / altEffHistY->Integral() /3.14;
  // numProjX->Scale( scaleFacX );
  // numProjY->Scale( scaleFacY );
  TEfficiency* effHistX = new TEfficiency(*numProjX,*denProjX);
  effHistX->SetName("effHistX");
  TEfficiency* effHistY = new TEfficiency(*numProjY,*denProjY);
  effHistY->SetName("effHistY");
  effHistX->SetTitle("Toy efficiency - x projection;;Efficiency");
  effHistY->SetTitle("Toy efficiency - y projection;;Efficiency");
  // effHistX->SetTitle(Form("Toy efficiency - x projection;;Efficiency * %f",scaleFacX));
  // effHistY->SetTitle(Form("Toy efficiency - y projection;;Efficiency * %f",scaleFacY));
  RooDataSet* fakeData = new RooDataSet("fakeData","fakeData",RooArgSet(x,y));
  for (int i=0; i<120; ++i) fakeData->add(RooArgSet(x,y));

  RooPlot* xframe1 = x.frame(Title("Projected function - x projection")) ;
  fakeData->plotOn(xframe1,Invisible()) ;
  effPdf.plotOn(xframe1,Slice(RooArgSet(y))) ;  
  // projectedPdf.plotOn(xframe1,LineColor(kRed)) ;
  // projectedPdfAs.plotOn(xframe1,LineColor(kBlack)) ;
  projectedPdf.plotOn(xframe1,LineColor(kRed),Slice(RooArgSet(y))) ;
  projectedPdfAs.plotOn(xframe1,LineColor(kBlack),Slice(RooArgSet(y))) ;

  RooPlot* yframe1 = y.frame(Title("Projected function - y projection")) ;
  fakeData->plotOn(yframe1,Invisible()) ;
  effPdf.plotOn(yframe1,Slice(RooArgSet(x))) ;  
  // projectedPdf.plotOn(yframe1,LineColor(kRed)) ;
  // projectedPdfAs.plotOn(yframe1,LineColor(kBlack)) ;
  projectedPdf.plotOn(yframe1,LineColor(kRed),Slice(RooArgSet(x))) ;
  projectedPdfAs.plotOn(yframe1,LineColor(kBlack),Slice(RooArgSet(x))) ;

  TCanvas* cx1 = new TCanvas("CanX_Proj_Func","Projected function - x projection",1600,800) ;
  gPad->SetLeftMargin(0.15) ; xframe1->GetYaxis()->SetTitleOffset(1.4) ; 
  effHistX->Draw();
  xframe1->Draw("same") ;

  TCanvas* cy1 = new TCanvas("CanY_Proj_Func","Projected function - y projection",1600,800) ;
  gPad->SetLeftMargin(0.15) ; yframe1->GetYaxis()->SetTitleOffset(1.4) ; 
  effHistY->Draw();
  yframe1->Draw("same") ;

  // projectedPdf.Print("v");
  // projectedPdfAs.Print("v");
}
