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
int nTot = 1e5;

void test3Deff()
{

  RooRealVar tK ("tK" ,"tK" ,-TMath::Pi(),TMath::Pi());
  RooRealVar tL ("tL" ,"tL" ,-TMath::Pi(),TMath::Pi());
  RooRealVar phi("phi","phi",-TMath::Pi(),TMath::Pi());
  RooArgSet vars (tK, tL, phi);

  //Generate PDFs for denominator, efficiency and numerator
  RooGenericPdf genpdf("genpdf","genpdf","4+(0*tK)-(0*tL)+(0*phi)",vars) ;

  // RooRealVar mean("mean","mean",0.);
  // RooRealVar sigma("sigma","sigma",1.2);
  // RooGaussian effYPdf("effYPdf","effYPdf",y,mean,sigma);
  // RooGenericPdf effXPdf("effXPdf","effXPdf","5-2/3.14*(x+abs(x))",x);
  // RooGenericPdf effYPdf("effYPdf","effYPdf","3+sin(y)",y);
  // RooGenericPdf effXPdf("effXPdf","effXPdf","3+cos(x)",x);
  // RooProdPdf effPdf("effPdf","effPdf",effXPdf,effYPdf);

  RooGenericPdf effPdf("effPdf","effPdf","0.2+0.01*cos(tK)*sin(tL)-0.07*sin(tK)*cos(phi)+0.07*cos(3*tK)*cos(tL)*sin(phi)",vars);
  // RooGenericPdf effPdf("effPdf","effPdf","(3-tL)*4*exp(-1.0/8*tK^2)+(3+0.3*tL*phi)*(tK^2/4+1)",vars);

  RooProdPdf numPdf("numPdf","numPdf",genpdf,effPdf);

  // Generate toy datasets
  RooDataSet* data = genpdf.generate(vars,nTot) ;
  RooDataSet* numData = numPdf.generate(vars,nTot/5) ;
  double avgEff = numData->sumEntries() / data->sumEntries();
  cout<<"Average efficiency = "<<avgEff<<endl;

  // Plot numerator and denominator datasets
  RooPlot* xframe = tK.frame(Title("Numerator and denominator theta(K) distributions")) ;
  data->plotOn(xframe) ;
  numData->plotOn(xframe) ;
  RooPlot* yframe = tL.frame(Title("Numerator and denominator theta(L) distributions")) ;
  data->plotOn(yframe) ;
  numData->plotOn(yframe) ;
  RooPlot* zframe = phi.frame(Title("Numerator and denominator phi distributions")) ;
  data->plotOn(zframe) ;
  numData->plotOn(zframe) ;

  TCanvas* c = new TCanvas("NumDenCanvas","toy_numerator_denominator",1200,800) ;
  c->Divide(3,1);
  c->cd(1);
  gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;
  c->cd(2);
  gPad->SetLeftMargin(0.15) ; yframe->GetYaxis()->SetTitleOffset(1.4) ; yframe->Draw() ;
  c->cd(3);
  gPad->SetLeftMargin(0.15) ; zframe->GetYaxis()->SetTitleOffset(1.4) ; zframe->Draw() ;

  //Declare and initialise the PDFs and all the needed objects
  // cout<<"Declaration"<<endl;
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
    
  for (int xOrder=0; xOrder<=maxOrder; ++xOrder) for (int yOrder=0; yOrder<=maxOrder; ++yOrder) for (int zOrder=0; zOrder<=maxOrder; ++zOrder) {

	int iOrder = (zOrder + yOrder*(maxOrder+1) + xOrder*(maxOrder+1)*(maxOrder+1))*8;

	factors.push_back( new RooRealVar(Form("s%is%is%i",xOrder,yOrder,zOrder),Form("s%is%is%i",xOrder,yOrder,zOrder),0) );
	factors.push_back( new RooRealVar(Form("s%is%ic%i",xOrder,yOrder,zOrder),Form("s%is%ic%i",xOrder,yOrder,zOrder),0) );
	factors.push_back( new RooRealVar(Form("s%ic%is%i",xOrder,yOrder,zOrder),Form("s%ic%is%i",xOrder,yOrder,zOrder),0) );
	factors.push_back( new RooRealVar(Form("s%ic%ic%i",xOrder,yOrder,zOrder),Form("s%ic%ic%i",xOrder,yOrder,zOrder),0) );
	factors.push_back( new RooRealVar(Form("c%is%is%i",xOrder,yOrder,zOrder),Form("c%is%is%i",xOrder,yOrder,zOrder),0) );
	factors.push_back( new RooRealVar(Form("c%is%ic%i",xOrder,yOrder,zOrder),Form("c%is%ic%i",xOrder,yOrder,zOrder),0) );
	factors.push_back( new RooRealVar(Form("c%ic%is%i",xOrder,yOrder,zOrder),Form("c%ic%is%i",xOrder,yOrder,zOrder),0) );
	factors.push_back( new RooRealVar(Form("c%ic%ic%i",xOrder,yOrder,zOrder),Form("c%ic%ic%i",xOrder,yOrder,zOrder),0) );

	factorsAs.push_back( new RooRealVar(Form("sa%isa%isa%i",xOrder,yOrder,zOrder),Form("sa%isa%isa%i",xOrder,yOrder,zOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("sa%isa%ica%i",xOrder,yOrder,zOrder),Form("sa%isa%ica%i",xOrder,yOrder,zOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("sa%ica%isa%i",xOrder,yOrder,zOrder),Form("sa%ica%isa%i",xOrder,yOrder,zOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("sa%ica%ica%i",xOrder,yOrder,zOrder),Form("sa%ica%ica%i",xOrder,yOrder,zOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("ca%isa%isa%i",xOrder,yOrder,zOrder),Form("ca%isa%isa%i",xOrder,yOrder,zOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("ca%isa%ica%i",xOrder,yOrder,zOrder),Form("ca%isa%ica%i",xOrder,yOrder,zOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("ca%ica%isa%i",xOrder,yOrder,zOrder),Form("ca%ica%isa%i",xOrder,yOrder,zOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("ca%ica%ica%i",xOrder,yOrder,zOrder),Form("ca%ica%ica%i",xOrder,yOrder,zOrder),0) );

	vectPdf  .push_back( new RooFormulaVar( Form("pdf_s%i_s%i_s%i",xOrder,yOrder,zOrder),   Form("pdf_s%i_s%i_s%i",xOrder,yOrder,zOrder),
						Form("sin(%i*tK)*sin(%i*tL)*sin(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_s%i_s%i_s%i",xOrder,yOrder,zOrder), Form("pdfas_s%i_s%i_s%i",xOrder,yOrder,zOrder),
						Form("sin(%i*tK)*sin(%i*tL)*sin(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdf  .push_back( new RooFormulaVar( Form("pdf_s%i_s%i_c%i",xOrder,yOrder,zOrder),   Form("pdf_s%i_s%i_c%i",xOrder,yOrder,zOrder),
						Form("sin(%i*tK)*sin(%i*tL)*cos(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_s%i_s%i_c%i",xOrder,yOrder,zOrder), Form("pdfas_s%i_s%i_c%i",xOrder,yOrder,zOrder),
						Form("sin(%i*tK)*sin(%i*tL)*cos(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdf  .push_back( new RooFormulaVar( Form("pdf_s%i_c%i_s%i",xOrder,yOrder,zOrder),   Form("pdf_s%i_c%i_s%i",xOrder,yOrder,zOrder),
						Form("sin(%i*tK)*cos(%i*tL)*sin(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_s%i_c%i_s%i",xOrder,yOrder,zOrder), Form("pdfas_s%i_c%i_s%i",xOrder,yOrder,zOrder),
						Form("sin(%i*tK)*cos(%i*tL)*sin(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdf  .push_back( new RooFormulaVar( Form("pdf_s%i_c%i_c%i",xOrder,yOrder,zOrder),   Form("pdf_s%i_c%i_c%i",xOrder,yOrder,zOrder),
						Form("sin(%i*tK)*cos(%i*tL)*cos(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_s%i_c%i_c%i",xOrder,yOrder,zOrder), Form("pdfas_s%i_c%i_c%i",xOrder,yOrder,zOrder),
						Form("sin(%i*tK)*cos(%i*tL)*cos(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdf  .push_back( new RooFormulaVar( Form("pdf_c%i_s%i_s%i",xOrder,yOrder,zOrder),   Form("pdf_c%i_s%i_s%i",xOrder,yOrder,zOrder),
						Form("cos(%i*tK)*sin(%i*tL)*sin(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_c%i_s%i_s%i",xOrder,yOrder,zOrder), Form("pdfas_c%i_s%i_s%i",xOrder,yOrder,zOrder),
						Form("cos(%i*tK)*sin(%i*tL)*sin(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdf  .push_back( new RooFormulaVar( Form("pdf_c%i_s%i_c%i",xOrder,yOrder,zOrder),   Form("pdf_c%i_s%i_c%i",xOrder,yOrder,zOrder),
						Form("cos(%i*tK)*sin(%i*tL)*cos(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_c%i_s%i_c%i",xOrder,yOrder,zOrder), Form("pdfas_c%i_s%i_c%i",xOrder,yOrder,zOrder),
						Form("cos(%i*tK)*sin(%i*tL)*cos(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdf  .push_back( new RooFormulaVar( Form("pdf_c%i_c%i_s%i",xOrder,yOrder,zOrder),   Form("pdf_c%i_c%i_s%i",xOrder,yOrder,zOrder),
						Form("cos(%i*tK)*cos(%i*tL)*sin(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_c%i_c%i_s%i",xOrder,yOrder,zOrder), Form("pdfas_c%i_c%i_s%i",xOrder,yOrder,zOrder),
						Form("cos(%i*tK)*cos(%i*tL)*sin(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdf  .push_back( new RooFormulaVar( Form("pdf_c%i_c%i_c%i",xOrder,yOrder,zOrder),   Form("pdf_c%i_c%i_c%i",xOrder,yOrder,zOrder),
						Form("cos(%i*tK)*cos(%i*tL)*cos(%i*phi)",xOrder,yOrder,zOrder), vars ) );
	vectPdfAs.push_back( new RooFormulaVar( Form("pdfas_c%i_c%i_c%i",xOrder,yOrder,zOrder), Form("pdfas_c%i_c%i_c%i",xOrder,yOrder,zOrder),
						Form("cos(%i*tK)*cos(%i*tL)*cos(%i*phi)",xOrder,yOrder,zOrder), vars ) );

	for (int i=0; i<8; ++i) {
	  proj.push_back(0);
	  projAs.push_back(0);
	}
	
	if ( xOrder>0 && yOrder>0 && zOrder>0 ) {
	  pdfList  .add(*vectPdf  [iOrder+0]);
	  pdfListAs.add(*vectPdfAs[iOrder+0]);
	  facList  .add(*factors  [iOrder+0]);
	  facListAs.add(*factorsAs[iOrder+0]);
	}
	if ( xOrder>0 && yOrder>0 ) {
	  pdfList  .add(*vectPdf  [iOrder+1]);
	  pdfListAs.add(*vectPdfAs[iOrder+1]);
	  facList  .add(*factors  [iOrder+1]);
	  facListAs.add(*factorsAs[iOrder+1]);
	}
	if ( xOrder>0 && zOrder>0 ) {
	  pdfList  .add(*vectPdf  [iOrder+2]);
	  pdfListAs.add(*vectPdfAs[iOrder+2]);
	  facList  .add(*factors  [iOrder+2]);
	  facListAs.add(*factorsAs[iOrder+2]);
	}
	if ( xOrder>0 ) {
	  pdfList  .add(*vectPdf  [iOrder+3]);
	  pdfListAs.add(*vectPdfAs[iOrder+3]);
	  facList  .add(*factors  [iOrder+3]);
	  facListAs.add(*factorsAs[iOrder+3]);
	}
	if ( yOrder>0 && zOrder>0 ) {
	  pdfList  .add(*vectPdf  [iOrder+4]);
	  pdfListAs.add(*vectPdfAs[iOrder+4]);
	  facList  .add(*factors  [iOrder+4]);
	  facListAs.add(*factorsAs[iOrder+4]);
	}
	if ( yOrder>0 ) {
	  pdfList  .add(*vectPdf  [iOrder+5]);
	  pdfListAs.add(*vectPdfAs[iOrder+5]);
	  facList  .add(*factors  [iOrder+5]);
	  facListAs.add(*factorsAs[iOrder+5]);
	}
	if ( zOrder>0 ) {
	  pdfList  .add(*vectPdf  [iOrder+6]);
	  pdfListAs.add(*vectPdfAs[iOrder+6]);
	  facList  .add(*factors  [iOrder+6]);
	  facListAs.add(*factorsAs[iOrder+6]);
	}
	{
	  pdfList  .add(*vectPdf  [iOrder+7]);
	  pdfListAs.add(*vectPdfAs[iOrder+7]);
	  facList  .add(*factors  [iOrder+7]);
	  facListAs.add(*factorsAs[iOrder+7]);
	}

      }


  // cout<<"Function declaration"<<endl;

  RooAddition projectedPdf   ("projectedPdf"  , "projectedPdf"  , pdfList  , facList  );
  RooAddition projectedPdfAs ("projectedPdfAs", "projectedPdfAs", pdfListAs, facListAs);

  //Create binned histos (only for projection method)
  // cout<<"Histo declaration"<<endl;
  TH3F* denHist = (TH3F*)data   ->createHistogram( "denHist",
						   tK,      Binning(20,-TMath::Pi(),TMath::Pi())  ,
						   YVar(tL, Binning(20,-TMath::Pi(),TMath::Pi())) ,
						   ZVar(phi,Binning(20,-TMath::Pi(),TMath::Pi())) );
  TH3F* numHist = (TH3F*)numData->createHistogram( "numHist",
						   tK,      Binning(20,-TMath::Pi(),TMath::Pi())  ,
						   YVar(tL, Binning(20,-TMath::Pi(),TMath::Pi())) ,
						   ZVar(phi,Binning(20,-TMath::Pi(),TMath::Pi())) );
  denHist->Sumw2();
  numHist->Sumw2();

  //Compute and set the coefficients
  // cout<<"Calculation"<<endl;
  factors  [7]->setVal(avgEff);
  factorsAs[7]->setVal(avgEff);

  double xCent, yCent, zCent, fact;

  for (int xOrder=0; xOrder<=maxOrder; ++xOrder) for (int yOrder=0; yOrder<=maxOrder; ++yOrder) for (int zOrder=0; zOrder<=maxOrder; ++zOrder) {
	
	int iOrder = (zOrder + yOrder*(maxOrder+1) + xOrder*(maxOrder+1)*(maxOrder+1))*8;

	// fill asymmetry coefficients
	// cout<<"fill asymmetry coefficients"<<endl;

	int totEvNum = numData->sumEntries();
	int totEvDen =    data->sumEntries();

	int countNum0 = numData->sumEntries( Form("sin(%i*tK)*sin(%i*tL)*sin(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countNum1 = numData->sumEntries( Form("sin(%i*tK)*sin(%i*tL)*cos(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countNum2 = numData->sumEntries( Form("sin(%i*tK)*cos(%i*tL)*sin(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countNum3 = numData->sumEntries( Form("sin(%i*tK)*cos(%i*tL)*cos(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countNum4 = numData->sumEntries( Form("cos(%i*tK)*sin(%i*tL)*sin(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countNum5 = numData->sumEntries( Form("cos(%i*tK)*sin(%i*tL)*cos(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countNum6 = numData->sumEntries( Form("cos(%i*tK)*cos(%i*tL)*sin(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countNum7 = numData->sumEntries( Form("cos(%i*tK)*cos(%i*tL)*cos(%i*phi)>0",xOrder,yOrder,zOrder) );

	int countDen0 = data->sumEntries( Form("sin(%i*tK)*sin(%i*tL)*sin(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countDen1 = data->sumEntries( Form("sin(%i*tK)*sin(%i*tL)*cos(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countDen2 = data->sumEntries( Form("sin(%i*tK)*cos(%i*tL)*sin(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countDen3 = data->sumEntries( Form("sin(%i*tK)*cos(%i*tL)*cos(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countDen4 = data->sumEntries( Form("cos(%i*tK)*sin(%i*tL)*sin(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countDen5 = data->sumEntries( Form("cos(%i*tK)*sin(%i*tL)*cos(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countDen6 = data->sumEntries( Form("cos(%i*tK)*cos(%i*tL)*sin(%i*phi)>0",xOrder,yOrder,zOrder) );
	int countDen7 = data->sumEntries( Form("cos(%i*tK)*cos(%i*tL)*cos(%i*phi)>0",xOrder,yOrder,zOrder) );
	// cout<<"-vars prepared"<<endl;

	if (countDen0) projAs[iOrder+0] = ( 1.0 * countNum0 / countDen0 ) - ( 1.0 * (totEvNum-countNum0) / (totEvDen-countDen0) );
	if (countDen1) projAs[iOrder+1] = ( 1.0 * countNum1 / countDen1 ) - ( 1.0 * (totEvNum-countNum1) / (totEvDen-countDen1) );
	if (countDen2) projAs[iOrder+2] = ( 1.0 * countNum2 / countDen2 ) - ( 1.0 * (totEvNum-countNum2) / (totEvDen-countDen2) );
	if (countDen3) projAs[iOrder+3] = ( 1.0 * countNum3 / countDen3 ) - ( 1.0 * (totEvNum-countNum3) / (totEvDen-countDen3) );
	if (countDen4) projAs[iOrder+4] = ( 1.0 * countNum4 / countDen4 ) - ( 1.0 * (totEvNum-countNum4) / (totEvDen-countDen4) );
	if (countDen5) projAs[iOrder+5] = ( 1.0 * countNum5 / countDen5 ) - ( 1.0 * (totEvNum-countNum5) / (totEvDen-countDen5) );
	if (countDen6) projAs[iOrder+6] = ( 1.0 * countNum6 / countDen6 ) - ( 1.0 * (totEvNum-countNum6) / (totEvDen-countDen6) );
	if (countDen7>0 && totEvDen>countDen7) projAs[iOrder+7] = ( 1.0 * countNum7 / countDen7 ) - ( 1.0 * (totEvNum-countNum7) / (totEvDen-countDen7) );
	// cout<<"-vector filled"<<endl;

	// fill standard coefficients
	// cout<<"fill standard coefficients"<<endl;

	for (int xBin=1;xBin<=denHist->GetNbinsX();++xBin) for (int yBin=1;yBin<=denHist->GetNbinsY();++yBin) for (int zBin=1;zBin<=denHist->GetNbinsZ();++zBin) {

	      xCent = denHist->GetXaxis()->GetBinCenter(xBin);
	      yCent = denHist->GetYaxis()->GetBinCenter(yBin);
	      zCent = denHist->GetZaxis()->GetBinCenter(zBin);

 	      if ( denHist->GetBinContent(xBin,yBin,zBin)>0 )
		fact = numHist->GetBinContent(xBin,yBin,zBin) / denHist->GetBinContent(xBin,yBin,zBin) *
		  denHist->GetXaxis()->GetBinWidth(xBin) *
		  denHist->GetYaxis()->GetBinWidth(yBin) *
		  denHist->GetZaxis()->GetBinWidth(zBin);
	      else fact=0;

	      if (xOrder==0) fact = fact/2.;
	      if (yOrder==0) fact = fact/2.;
	      if (zOrder==0) fact = fact/2.;

	      proj[iOrder+0] += TMath::Sin(xOrder*xCent) * TMath::Sin(yOrder*yCent) * TMath::Sin(zOrder*zCent) * fact;
	      proj[iOrder+1] += TMath::Sin(xOrder*xCent) * TMath::Sin(yOrder*yCent) * TMath::Cos(zOrder*zCent) * fact;
	      proj[iOrder+2] += TMath::Sin(xOrder*xCent) * TMath::Cos(yOrder*yCent) * TMath::Sin(zOrder*zCent) * fact;
	      proj[iOrder+3] += TMath::Sin(xOrder*xCent) * TMath::Cos(yOrder*yCent) * TMath::Cos(zOrder*zCent) * fact;
	      proj[iOrder+4] += TMath::Cos(xOrder*xCent) * TMath::Sin(yOrder*yCent) * TMath::Sin(zOrder*zCent) * fact;
	      proj[iOrder+5] += TMath::Cos(xOrder*xCent) * TMath::Sin(yOrder*yCent) * TMath::Cos(zOrder*zCent) * fact;
	      proj[iOrder+6] += TMath::Cos(xOrder*xCent) * TMath::Cos(yOrder*yCent) * TMath::Sin(zOrder*zCent) * fact;
	      proj[iOrder+7] += TMath::Cos(xOrder*xCent) * TMath::Cos(yOrder*yCent) * TMath::Cos(zOrder*zCent) * fact;

	    }
	// cout<<"-vector filled"<<endl;

	if ( iOrder>0 ) for (int i=0; i<8; ++i) {
	    factors  [iOrder+i]->setVal(proj  [iOrder+i]/TMath::Pi()/TMath::Pi()/TMath::Pi());
	    factorsAs[iOrder+i]->setVal(projAs[iOrder+i]*TMath::Pi()/4);
	  }
	// cout<<"Variables normalized"<<endl;

	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<" sss\t"<<proj[iOrder+0]/TMath::Pi()/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+0]*TMath::Pi()/4<<endl;
	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<" ssc\t"<<proj[iOrder+1]/TMath::Pi()/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+1]*TMath::Pi()/4<<endl;
	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<" scs\t"<<proj[iOrder+2]/TMath::Pi()/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+2]*TMath::Pi()/4<<endl;
	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<" scc\t"<<proj[iOrder+3]/TMath::Pi()/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+3]*TMath::Pi()/4<<endl;
	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<" css\t"<<proj[iOrder+4]/TMath::Pi()/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+4]*TMath::Pi()/4<<endl;
	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<" csc\t"<<proj[iOrder+5]/TMath::Pi()/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+5]*TMath::Pi()/4<<endl;
	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<" ccs\t"<<proj[iOrder+6]/TMath::Pi()/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+6]*TMath::Pi()/4<<endl;
	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<" ccc\t"<<proj[iOrder+7]/TMath::Pi()/TMath::Pi()/TMath::Pi()<<"\t"<<projAs[iOrder+7]*TMath::Pi()/4<<endl;
  
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
  RooDataSet* fakeData = new RooDataSet("fakeData","fakeData",vars);
  for (int i=0; i<790; ++i) fakeData->add(vars);

  TCanvas* cx1 = new TCanvas("CanX_Proj_Func","Projected function - x projection",2000,2000) ;
  TCanvas* cy1 = new TCanvas("CanY_Proj_Func","Projected function - y projection",2000,2000) ;
  TCanvas* cz1 = new TCanvas("CanZ_Proj_Func","Projected function - z projection",2000,2000) ;
  cx1->Divide(5,5);
  cy1->Divide(5,5);
  cz1->Divide(5,5);

  double border = 0.3;

  vector <TEfficiency*> effHistsX; 
  vector <TEfficiency*> effHistsY;
  vector <TEfficiency*> effHistsZ;
  vector <RooPlot*> xframes;
  vector <RooPlot*> yframes;
  vector <RooPlot*> zframes;

  TLegend* leg = new TLegend (0.5,0.7,0.9,0.9);

  for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {

    double centA = -3 + 1.5*i;
    double centB = -3 + 1.5*j;
    double lowA  = TMath::Max( centA - border,  1e-4-TMath::Pi() );
    double lowB  = TMath::Max( centB - border,  1e-4-TMath::Pi() );
    double highA = TMath::Min( centA + border, -1e-4+TMath::Pi() );
    double highB = TMath::Min( centB + border, -1e-4+TMath::Pi() );
    
    auto numProjX = numHist->ProjectionX("numProjX", 
					 numHist->GetYaxis()->FindBin(lowA), numHist->GetYaxis()->FindBin(highA),
					 numHist->GetZaxis()->FindBin(lowB), numHist->GetZaxis()->FindBin(highB),"e");
    auto numProjY = numHist->ProjectionY("numProjY", 
					 numHist->GetXaxis()->FindBin(lowA), numHist->GetXaxis()->FindBin(highA),
					 numHist->GetZaxis()->FindBin(lowB), numHist->GetZaxis()->FindBin(highB),"e");
    auto numProjZ = numHist->ProjectionZ("numProjZ", 
					 numHist->GetXaxis()->FindBin(lowA), numHist->GetXaxis()->FindBin(highA),
					 numHist->GetYaxis()->FindBin(lowB), numHist->GetYaxis()->FindBin(highB),"e");
    auto denProjX = denHist->ProjectionX("denProjX", 
					 denHist->GetYaxis()->FindBin(lowA), denHist->GetYaxis()->FindBin(highA),
					 denHist->GetZaxis()->FindBin(lowB), denHist->GetZaxis()->FindBin(highB),"e");
    auto denProjY = denHist->ProjectionY("denProjY", 
					 denHist->GetXaxis()->FindBin(lowA), denHist->GetXaxis()->FindBin(highA),
					 denHist->GetZaxis()->FindBin(lowB), denHist->GetZaxis()->FindBin(highB),"e");
    auto denProjZ = denHist->ProjectionZ("denProjZ", 
					 denHist->GetXaxis()->FindBin(lowA), denHist->GetXaxis()->FindBin(highA),
					 denHist->GetYaxis()->FindBin(lowB), denHist->GetYaxis()->FindBin(highB),"e");

    
    effHistsX.push_back( new TEfficiency(*numProjX,*denProjX) );
    effHistsX.back()->SetName( Form("effHistX_%i_%i",i,j) );
    effHistsX.back()->SetTitle(Form("Toy efficiency - slice y=%1.1f z=%1.1f;x;Efficiency",centA,centB) );
    
    effHistsY.push_back( new TEfficiency(*numProjY,*denProjY) );
    effHistsY.back()->SetName( Form("effHistY_%i_%i",i,j) );
    effHistsY.back()->SetTitle(Form("Toy efficiency - slice x=%1.1f z=%1.1f;y;Efficiency",centA,centB) );

    effHistsZ.push_back( new TEfficiency(*numProjZ,*denProjZ) );
    effHistsZ.back()->SetName( Form("effHistZ_%i_%i",i,j) );
    effHistsZ.back()->SetTitle(Form("Toy efficiency - slice x=%1.1f y=%1.1f;z;Efficiency",centA,centB) );

    tL .setVal(centA);
    phi.setVal(centB);
    xframes.push_back( tK.frame(Title( Form("Projected function - x projection %i %i",i,j))) );
    fakeData->     plotOn(xframes.back(),Invisible()) ;
    effPdf.        plotOn(xframes.back(),Slice(RooArgSet(tL,phi))                  ,Name(Form("effPdfx_%i_%i"        ,i,j))) ;  
    projectedPdf.  plotOn(xframes.back(),LineColor(kRed),Slice(RooArgSet(tL,phi))  ,Name(Form("projectedPdfx_%i_%i"  ,i,j))) ;
    projectedPdfAs.plotOn(xframes.back(),LineColor(kBlack),Slice(RooArgSet(tL,phi)),Name(Form("projectedPdfAsx_%i_%i",i,j))) ;

    tK .setVal(centA);
    phi.setVal(centB);
    yframes.push_back( tL.frame(Title( Form("Projected function - y projection %i %i",i,j))) );
    fakeData->     plotOn(yframes.back(),Invisible()) ;
    effPdf.        plotOn(yframes.back(),Slice(RooArgSet(tK,phi))                  ,Name(Form("effPdfy_%i_%i"        ,i,j))) ;  
    projectedPdf.  plotOn(yframes.back(),LineColor(kRed),Slice(RooArgSet(tK,phi))  ,Name(Form("projectedPdfy_%i_%i"  ,i,j))) ;
    projectedPdfAs.plotOn(yframes.back(),LineColor(kBlack),Slice(RooArgSet(tK,phi)),Name(Form("projectedPdfAsy_%i_%i",i,j))) ;

    tK.setVal(centA);
    tL.setVal(centB);
    zframes.push_back( phi.frame(Title( Form("Projected function - z projection %i %i",i,j))) );
    fakeData->     plotOn(zframes.back(),Invisible()) ;
    effPdf.        plotOn(zframes.back(),Slice(RooArgSet(tK,tL))                  ,Name(Form("effPdfz_%i_%i"        ,i,j))) ;  
    projectedPdf.  plotOn(zframes.back(),LineColor(kRed),Slice(RooArgSet(tK,tL))  ,Name(Form("projectedPdfz_%i_%i"  ,i,j))) ;
    projectedPdfAs.plotOn(zframes.back(),LineColor(kBlack),Slice(RooArgSet(tK,tL)),Name(Form("projectedPdfAsz_%i_%i",i,j))) ;

    cx1->cd(5*j+i+1);
    effHistsX.back()->Draw();
    cx1->cd(5*j+i+1)->Update(); 
    auto graphx = effHistsX.back()->GetPaintedGraph(); 
    graphx->SetMinimum(0);
    graphx->SetMaximum(1);
    graphx->GetYaxis()->SetTitleOffset(1.4);
    cx1->cd(5*j+i+1)->Update();
    xframes.back()->Draw("same") ;

    if (i+j==0) {
      leg->AddEntry(effHistsX.back()                 ,"Binned efficiency","lep");
      leg->AddEntry(Form("effPdfx_%i_%i"        ,i,j),"GEN efficiency","l");
      leg->AddEntry(Form("projectedPdfx_%i_%i"  ,i,j),"Standard projection","l");
      leg->AddEntry(Form("projectedPdfAsx_%i_%i",i,j),"Asymmetry projection","l");
    }

    leg->Draw("same");

    cy1->cd(5*j+i+1);
    effHistsY.back()->Draw();
    cy1->cd(5*j+i+1)->Update(); 
    auto graphy = effHistsY.back()->GetPaintedGraph(); 
    graphy->SetMinimum(0);
    graphy->SetMaximum(1); 
    graphy->GetYaxis()->SetTitleOffset(1.4);
    cy1->cd(5*j+i+1)->Update();
    yframes.back()->Draw("same") ;
    leg->Draw("same");

    cz1->cd(5*j+i+1);
    effHistsZ.back()->Draw();
    cz1->cd(5*j+i+1)->Update(); 
    auto graphz = effHistsZ.back()->GetPaintedGraph(); 
    graphz->SetMinimum(0);
    graphz->SetMaximum(1); 
    graphz->GetYaxis()->SetTitleOffset(1.4);
    cz1->cd(5*j+i+1)->Update();
    zframes.back()->Draw("same") ;
    leg->Draw("same");

  }    
    
}
