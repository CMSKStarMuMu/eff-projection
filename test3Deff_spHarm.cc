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

string plotName = "3d_test_6o-sh_custFunct";

int maxOrder =6;
int nTot = 1e5;

vector< double > getAsymmcoefficients(vector< RooProduct* > funcs, vector< double > maxVal);
int doubleFactorial(int i);

RooRealVar ctK ("ctK","ctK",-1,1);
RooRealVar ctL ("ctL","ctL",-1,1);
RooRealVar phi ("phi","phi",-TMath::Pi(),TMath::Pi());
RooArgSet vars (ctK, ctL, phi);

void test3Deff_spHarm()
{
  //Generate PDFs for denominator, efficiency and numerator
  RooGenericPdf genpdf("genpdf","genpdf","4+(0*ctK)-(0*ctL)+(0*phi)",vars) ;

  // RooRealVar mean("mean","mean",0.);
  // RooRealVar sigma("sigma","sigma",1.2);
  // RooGaussian effYPdf("effYPdf","effYPdf",y,mean,sigma);
  // RooGenericPdf effXPdf("effXPdf","effXPdf","5-2/3.14*(x+abs(x))",x);
  // RooGenericPdf effYPdf("effYPdf","effYPdf","3+sin(y)",y);
  // RooGenericPdf effXPdf("effXPdf","effXPdf","3+cos(x)",x);
  // RooProdPdf effPdf("effPdf","effPdf",effXPdf,effYPdf);

  // RooGenericPdf effPdf("effPdf","effPdf","0.2+0.05*(1-ctK*ctK)*(1-ctL*ctL)*cos(2*phi)-0.05*sqrt(1-ctL*ctL)*(7*ctK*ctK*ctK-3*ctK)*sqrt(1-ctK*ctK)*sin(phi)+0.05*6*ctK*(1-ctK*ctK)*ctL*(1-ctL*ctL)*sin(2*phi)",vars);
  // RooGenericPdf effPdf("effPdf","effPdf","0.2+0.01*cos(ctK)*sin(ctL)-0.07*sin(ctK)*cos(phi)+0.07*cos(3*ctK)*cos(ctL)*sin(phi)",vars);
  RooGenericPdf effPdf("effPdf","effPdf","(1-ctL)*4*exp(-1.0*ctK^2)+(1+ctL*cos(phi))*(2*ctK^2+1)",vars);

  RooProdPdf numPdf("numPdf","numPdf",genpdf,effPdf);

  // Generate toy datasets
  RooDataSet* data = genpdf.generate(vars,nTot) ;
  RooDataSet* numData = numPdf.generate(vars,nTot/5) ;
  double avgEff = numData->sumEntries() / data->sumEntries();
  cout<<"Average efficiency = "<<avgEff<<endl;

  // Plot numerator and denominator datasets
  RooPlot* xframe = ctK.frame(Title("Numerator and denominator cos(theta_K) distributions")) ;
  data->plotOn(xframe) ;
  numData->plotOn(xframe) ;
  RooPlot* yframe = ctL.frame(Title("Numerator and denominator cos(theta_L) distributions")) ;
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
  vector < RooLegendre* > vectPdfLegCosK;
  vector < RooLegendre* > vectPdfLegCosL;
  vector < RooFormulaVar* > vectPdfPoly;
  vector < RooProduct* > vectPdf;

  vector < RooRealVar* > factorsAs;
  vector < double > projAs;
  vector < RooLegendre* > vectPdfLegCosKAs;
  vector < RooLegendre* > vectPdfLegCosLAs;
  vector < RooFormulaVar* > vectPdfPolyAs;
  vector < RooProduct* > vectPdfAs;

  vector< double > maxVal;

  RooArgList facList;
  RooArgList facListAs;
  RooArgList pdfList;
  RooArgList pdfListAs;

  RooPlot* xframeT = ctK.frame(Title("Numerator and denominator cos(theta_K) distributions")) ;
  RooPlot* yframeT = ctL.frame(Title("Numerator and denominator cos(theta_L) distributions")) ;
  RooPlot* zframeT = phi.frame(Title("Numerator and denominator phi distributions")) ;
    
  for (int xOrder=0; xOrder<=maxOrder; ++xOrder)
    for (int yOrder=0; yOrder<=maxOrder; ++yOrder)
      for (int zOrder=-1*TMath::Min(xOrder,yOrder); zOrder<=TMath::Min(xOrder,yOrder); ++zOrder) {

	factors  .push_back( new RooRealVar(Form( "l%i_k%i_m%i",xOrder,yOrder,zOrder),Form( "l%i_k%i_m%i",xOrder,yOrder,zOrder),0) );
	factorsAs.push_back( new RooRealVar(Form("Al%i_k%i_m%i",xOrder,yOrder,zOrder),Form("Al%i_k%i_m%i",xOrder,yOrder,zOrder),0) );

	RooArgList prodList;
	RooArgList prodListAs;

	if (zOrder>0) {
	  vectPdfPoly  .push_back( new RooFormulaVar( Form("pdfPoly%i_%i_%i"  ,xOrder,yOrder,zOrder), Form("pdfPoly%i_%i_%i"  ,xOrder,yOrder,zOrder),
						      Form("cos(%i*phi)",zOrder), phi ) );
	  vectPdfPolyAs.push_back( new RooFormulaVar( Form("pdfPolyas%i_%i_%i",xOrder,yOrder,zOrder), Form("pdfPolyas%i_%i_%i",xOrder,yOrder,zOrder),
						      Form("cos(%i*phi)",zOrder), phi ) );
	  prodList  .add( *vectPdfPoly  .back() );
	  prodListAs.add( *vectPdfPolyAs.back() );
	}
	if (zOrder<0) {
	  vectPdfPoly  .push_back( new RooFormulaVar( Form("pdfPoly%i_%i_%i"  ,xOrder,yOrder,zOrder), Form("pdfPoly%i_%i_%i"  ,xOrder,yOrder,zOrder),
						      Form("sin(%i*phi)",-1*zOrder), phi ) );
	  vectPdfPolyAs.push_back( new RooFormulaVar( Form("pdfPolyas%i_%i_%i",xOrder,yOrder,zOrder), Form("pdfPolyas%i_%i_%i",xOrder,yOrder,zOrder),
						      Form("sin(%i*phi)",-1*zOrder), phi ) );
	  prodList  .add( *vectPdfPoly  .back() );
	  prodListAs.add( *vectPdfPolyAs.back() );
	}

	vectPdfLegCosK  .push_back( new RooLegendre ( Form("pdfLegctK%i_%i_%i"  ,xOrder,yOrder,zOrder),
						      Form("pdfLegctK%i_%i_%i"  ,xOrder,yOrder,zOrder), ctK, xOrder, abs(zOrder) ) );
	vectPdfLegCosKAs.push_back( new RooLegendre ( Form("pdfLegctKAs%i_%i_%i",xOrder,yOrder,zOrder),
						      Form("pdfLegctKAs%i_%i_%i",xOrder,yOrder,zOrder), ctK, xOrder, abs(zOrder) ) );
	prodList  .add( *vectPdfLegCosK  .back() );
	prodListAs.add( *vectPdfLegCosKAs.back() );

	vectPdfLegCosL  .push_back( new RooLegendre ( Form("pdfLegctL%i_%i_%i"  ,xOrder,yOrder,zOrder),
						      Form("pdfLegctL%i_%i_%i"  ,xOrder,yOrder,zOrder), ctL, yOrder, abs(zOrder) ) );
	vectPdfLegCosLAs.push_back( new RooLegendre ( Form("pdfLegctLAs%i_%i_%i",xOrder,yOrder,zOrder),
						      Form("pdfLegctLAs%i_%i_%i",xOrder,yOrder,zOrder), ctL, yOrder, abs(zOrder) ) );
	prodList  .add( *vectPdfLegCosL  .back() );
	prodListAs.add( *vectPdfLegCosLAs.back() );

	vectPdf  .push_back( new RooProduct ( Form("pdf%i_%i_%i"  ,xOrder,yOrder,zOrder), Form("pdf%i_%i_%i"  ,xOrder,yOrder,zOrder), prodList   ) );
	vectPdfAs.push_back( new RooProduct ( Form("pdfAs%i_%i_%i",xOrder,yOrder,zOrder), Form("pdfAs%i_%i_%i",xOrder,yOrder,zOrder), prodListAs ) );

	// maxVal.push_back( 25 );
	maxVal.push_back( 1.5 * TMath::Power( doubleFactorial(2*xOrder-1), 1.0*abs(zOrder)/xOrder ) *
			  1.5 * TMath::Power( doubleFactorial(2*yOrder-1), 1.0*abs(zOrder)/yOrder ) );
	// cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<"\t"
	//     <<2*xOrder-1<<" "<<doubleFactorial(2*xOrder-1)<<" "<<1.5 * TMath::Power( doubleFactorial(2*xOrder-1), 1.0*abs(zOrder)/xOrder )<<"\t"
	//     <<2*yOrder-1<<" "<<doubleFactorial(2*yOrder-1)<<" "<<1.5 * TMath::Power( doubleFactorial(2*yOrder-1), 1.0*abs(zOrder)/yOrder )<<endl;

	proj  .push_back(0);
	projAs.push_back(0);
	
	pdfList  .add( *vectPdf  .back() );
	pdfListAs.add( *vectPdfAs.back() );
	facList  .add( *factors  .back() );
	facListAs.add( *factorsAs.back() );

	// vectPdfAs.back()->plotOn(xframeT,LineColor(vectPdfAs.size()+10),Slice(RooArgSet(ctL,phi))) ;
	// vectPdfAs.back()->plotOn(yframeT,LineColor(vectPdfAs.size()+10),Slice(RooArgSet(ctK,phi))) ;
	// vectPdfAs.back()->plotOn(zframeT,LineColor(vectPdfAs.size()+10),Slice(RooArgSet(ctK,ctL))) ;
	if ( xOrder==yOrder ) {
	  vectPdfLegCosKAs.back()->plotOn(xframeT,LineColor(vectPdfAs.size()+10));
	  vectPdfLegCosLAs.back()->plotOn(yframeT,LineColor(vectPdfAs.size()+10));
	  if ( zOrder!=0 && abs(zOrder)==xOrder ) vectPdfPolyAs.back()->plotOn(zframeT,LineColor(vectPdfAs.size()+10));
	}
      }

  // Plot elements of the function basis

  TCanvas* cT = new TCanvas("NumDenCanvasT","toy_numerator_denominator",1200,800) ;
  cT->Divide(3,1);
  cT->cd(1);
  gPad->SetLeftMargin(0.15) ; xframeT->GetYaxis()->SetTitleOffset(1.4) ; xframeT->Draw() ;
  cT->cd(2);
  gPad->SetLeftMargin(0.15) ; yframeT->GetYaxis()->SetTitleOffset(1.4) ; yframeT->Draw() ;
  cT->cd(3);
  gPad->SetLeftMargin(0.15) ; zframeT->GetYaxis()->SetTitleOffset(1.4) ; zframeT->Draw() ;

  normFacs = getAsymmcoefficients( vectPdfAs, maxVal );

  // cout<<"Function declaration"<<endl;

  RooAddition projectedPdf   ("projectedPdf"  , "projectedPdf"  , pdfList  , facList  );
  RooAddition projectedPdfAs ("projectedPdfAs", "projectedPdfAs", pdfListAs, facListAs);

  //Create binned histos (only for projection method)
  // cout<<"Histo declaration"<<endl;
  TH3F* denHist = (TH3F*)data   ->createHistogram( "denHist",
						   ctK,     Binning(20,-1,1)  ,
						   YVar(ctL,Binning(20,-1,1)) ,
						   ZVar(phi,Binning(20,-TMath::Pi(),TMath::Pi())) );
  TH3F* numHist = (TH3F*)numData->createHistogram( "numHist",
						   ctK,     Binning(20,-1,1)  ,
						   YVar(ctL,Binning(20,-1,1)) ,
						   ZVar(phi,Binning(20,-TMath::Pi(),TMath::Pi())) );
  denHist->Sumw2();
  numHist->Sumw2();

  //Compute and set the coefficients
  // cout<<"Calculation"<<endl;
  factors  [0]->setVal(avgEff);
  factorsAs[0]->setVal(avgEff);

  double fact;
  int iOrder=-1;

  int evNum = numData->sumEntries();
  int evDen =    data->sumEntries();

  for (int xOrder=0; xOrder<=maxOrder; ++xOrder)
    for (int yOrder=0; yOrder<=maxOrder; ++yOrder)
      for (int zOrder=-1*TMath::Min(xOrder,yOrder); zOrder<=TMath::Min(xOrder,yOrder); ++zOrder) {
	
	++iOrder;

	// fill asymmetry coefficients
	// cout<<"fill asymmetry coefficients"<<endl;

	int countNum = 0;
	int countDen = 0;
	
	int totEvNum = 0;
	int totEvDen = 0;

	for (int i=0; i<evNum; ++i) {
	  const RooArgSet *iPoint = numData->get(i);
	  ctK.setVal( iPoint->getRealValue("ctK") );
	  ctL.setVal( iPoint->getRealValue("ctL") );
	  phi.setVal( iPoint->getRealValue("phi") );
	  if ( vectPdfAs[iOrder]->getVal() >  0 ) ++countNum;
	  if ( vectPdfAs[iOrder]->getVal() != 0 ) ++totEvNum;
	}
	for (int i=0; i<evDen; ++i) {
	  const RooArgSet *iPoint = data->get(i);
	  ctK.setVal( iPoint->getRealValue("ctK") );
	  ctL.setVal( iPoint->getRealValue("ctL") );
	  phi.setVal( iPoint->getRealValue("phi") );
	  if ( vectPdfAs[iOrder]->getVal() >  0 ) ++countDen;
	  if ( vectPdfAs[iOrder]->getVal() != 0 ) ++totEvDen;
	}
	// cout<<"-vars prepared"<<endl;

	// cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<" ->\t"<<countNum<<"/"<<countDen<<"\t"<<totEvNum-countNum<<"/"<<totEvDen-countDen<<endl;
	if (countDen>0 && totEvDen>countDen) projAs[iOrder] = ( 1.0 * countNum / countDen ) - ( 1.0 * (totEvNum-countNum) / (totEvDen-countDen) );
	// cout<<"-vector filled"<<endl;

	// fill standard coefficients
	// cout<<"fill standard coefficients"<<endl;

	for (int xBin=1;xBin<=denHist->GetNbinsX();++xBin) for (int yBin=1;yBin<=denHist->GetNbinsY();++yBin) for (int zBin=1;zBin<=denHist->GetNbinsZ();++zBin) {

	      ctK.setVal( denHist->GetXaxis()->GetBinCenter(xBin) );
	      ctL.setVal( denHist->GetYaxis()->GetBinCenter(yBin) );
	      phi.setVal( denHist->GetZaxis()->GetBinCenter(zBin) );

 	      if ( denHist->GetBinContent(xBin,yBin,zBin)>0 )
		proj[iOrder] += numHist->GetBinContent(xBin,yBin,zBin) / denHist->GetBinContent(xBin,yBin,zBin) *
		  denHist->GetXaxis()->GetBinWidth(xBin) *
		  denHist->GetYaxis()->GetBinWidth(yBin) *
		  denHist->GetZaxis()->GetBinWidth(zBin) *
		  vectPdfAs[iOrder]->getVal( vars );

	    }
	// cout<<"-vector filled"<<endl;

	if (zOrder==0) proj[iOrder] = proj[iOrder]/2.0;

	if (iOrder>0) {
	  factors  [iOrder]->setVal( proj[iOrder]
				     *(2*xOrder+1)*TMath::Factorial(xOrder-abs(zOrder))/2/TMath::Factorial(xOrder+abs(zOrder))
				     *(2*yOrder+1)*TMath::Factorial(yOrder-abs(zOrder))/2/TMath::Factorial(yOrder+abs(zOrder))
				     /TMath::Pi() );
	  factorsAs[iOrder]->setVal( projAs[iOrder] / normFacs[iOrder] );
	  // cout<<"Variables normalized"<<endl;
	}

	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<"\t"<<factors[iOrder]->getValV()<<"\t"<<factorsAs[iOrder]->getValV()<<endl;
  
      }

  //Plot the efficiency slices
  RooDataSet* fakeData = new RooDataSet("fakeData","fakeData",vars);
  for (int i=0; i<251; ++i) fakeData->add(vars);
  RooDataSet* fakeData1 = new RooDataSet("fakeData1","fakeData1",vars);
  for (int i=0; i<80; ++i) fakeData1->add(vars);

  TCanvas* cx1 = new TCanvas("CanX_Proj_Func","Projected function - x projection",2000,2000) ;
  TCanvas* cy1 = new TCanvas("CanY_Proj_Func","Projected function - y projection",2000,2000) ;
  TCanvas* cz1 = new TCanvas("CanZ_Proj_Func","Projected function - z projection",2000,2000) ;
  cx1->Divide(5,5);
  cy1->Divide(5,5);
  cz1->Divide(5,5);

  double border = 0.1;

  vector <TEfficiency*> effHistsX; 
  vector <TEfficiency*> effHistsY;
  vector <TEfficiency*> effHistsZ;
  vector <RooPlot*> xframes;
  vector <RooPlot*> yframes;
  vector <RooPlot*> zframes;

  TLegend* leg = new TLegend (0.5,0.7,0.9,0.9);

  for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {

    double centA = -0.9 + 1.8*i/4;
    double centB = -0.9 + 1.8*j/4;
    double lowA  = TMath::Max( centA - border,  1e-4-1 );
    double lowB  = TMath::Max( centB - border,  1e-4-1 );
    double highA = TMath::Min( centA + border, -1e-4+1 );
    double highB = TMath::Min( centB + border, -1e-4+1 );
    
    auto numProjX = numHist->ProjectionX("numProjX", 
					 numHist->GetYaxis()->FindBin(lowA            ), numHist->GetYaxis()->FindBin(highA            ),
					 numHist->GetZaxis()->FindBin(lowB*TMath::Pi()), numHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
    auto numProjY = numHist->ProjectionY("numProjY", 
					 numHist->GetXaxis()->FindBin(lowA            ), numHist->GetXaxis()->FindBin(highA            ),
					 numHist->GetZaxis()->FindBin(lowB*TMath::Pi()), numHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
    auto numProjZ = numHist->ProjectionZ("numProjZ", 
					 numHist->GetXaxis()->FindBin(lowA            ), numHist->GetXaxis()->FindBin(highA            ),
					 numHist->GetYaxis()->FindBin(lowB            ), numHist->GetYaxis()->FindBin(highB            ),"e");
    auto denProjX = denHist->ProjectionX("denProjX", 
					 denHist->GetYaxis()->FindBin(lowA            ), denHist->GetYaxis()->FindBin(highA            ),
					 denHist->GetZaxis()->FindBin(lowB*TMath::Pi()), denHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
    auto denProjY = denHist->ProjectionY("denProjY", 
					 denHist->GetXaxis()->FindBin(lowA            ), denHist->GetXaxis()->FindBin(highA            ),
					 denHist->GetZaxis()->FindBin(lowB*TMath::Pi()), denHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
    auto denProjZ = denHist->ProjectionZ("denProjZ", 
					 denHist->GetXaxis()->FindBin(lowA            ), denHist->GetXaxis()->FindBin(highA            ),
					 denHist->GetYaxis()->FindBin(lowB            ), denHist->GetYaxis()->FindBin(highB            ),"e");

    
    effHistsX.push_back( new TEfficiency(*numProjX,*denProjX) );
    effHistsX.back()->SetName( Form("effHistX_%i_%i",i,j) );
    effHistsX.back()->SetTitle(Form("Toy efficiency - slice y=%1.2f z=%1.2f;x;Efficiency",centA,centB*TMath::Pi()) );
    
    effHistsY.push_back( new TEfficiency(*numProjY,*denProjY) );
    effHistsY.back()->SetName( Form("effHistY_%i_%i",i,j) );
    effHistsY.back()->SetTitle(Form("Toy efficiency - slice x=%1.2f z=%1.2f;y;Efficiency",centA,centB*TMath::Pi()) );

    effHistsZ.push_back( new TEfficiency(*numProjZ,*denProjZ) );
    effHistsZ.back()->SetName( Form("effHistZ_%i_%i",i,j) );
    effHistsZ.back()->SetTitle(Form("Toy efficiency - slice x=%1.2f y=%1.2f;z;Efficiency",centA,centB) );

    ctL.setVal(centA);
    phi.setVal(centB*TMath::Pi());
    xframes.push_back( ctK.frame(Title( Form("Projected function - x projection %i %i",i,j))) );
    fakeData->     plotOn(xframes.back(),Invisible()) ;
    effPdf.        plotOn(xframes.back(),Slice(RooArgSet(ctL,phi))                  ,Name(Form("effPdfx_%i_%i"        ,i,j))) ;  
    projectedPdf.  plotOn(xframes.back(),LineColor(kRed),Slice(RooArgSet(ctL,phi))  ,Name(Form("projectedPdfx_%i_%i"  ,i,j))) ;
    projectedPdfAs.plotOn(xframes.back(),LineColor(kBlack),Slice(RooArgSet(ctL,phi)),Name(Form("projectedPdfAsx_%i_%i",i,j))) ;

    ctK.setVal(centA);
    phi.setVal(centB*TMath::Pi());
    yframes.push_back( ctL.frame(Title( Form("Projected function - y projection %i %i",i,j))) );
    fakeData->     plotOn(yframes.back(),Invisible()) ;
    effPdf.        plotOn(yframes.back(),Slice(RooArgSet(ctK,phi))                  ,Name(Form("effPdfy_%i_%i"        ,i,j))) ;  
    projectedPdf.  plotOn(yframes.back(),LineColor(kRed),Slice(RooArgSet(ctK,phi))  ,Name(Form("projectedPdfy_%i_%i"  ,i,j))) ;
    projectedPdfAs.plotOn(yframes.back(),LineColor(kBlack),Slice(RooArgSet(ctK,phi)),Name(Form("projectedPdfAsy_%i_%i",i,j))) ;

    ctK.setVal(centA);
    ctL.setVal(centB);
    zframes.push_back( phi.frame(Title( Form("Projected function - z projection %i %i",i,j))) );
    fakeData1->    plotOn(zframes.back(),Invisible()) ;
    effPdf.        plotOn(zframes.back(),Slice(RooArgSet(ctK,ctL))                  ,Name(Form("effPdfz_%i_%i"        ,i,j))) ;  
    projectedPdf.  plotOn(zframes.back(),LineColor(kRed),Slice(RooArgSet(ctK,ctL))  ,Name(Form("projectedPdfz_%i_%i"  ,i,j))) ;
    projectedPdfAs.plotOn(zframes.back(),LineColor(kBlack),Slice(RooArgSet(ctK,ctL)),Name(Form("projectedPdfAsz_%i_%i",i,j))) ;

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

  cx1->SaveAs( (plotName + "_x.pdf").c_str() );
  cy1->SaveAs( (plotName + "_y.pdf").c_str() );
  cz1->SaveAs( (plotName + "_z.pdf").c_str() );
    
}

vector< double > getAsymmcoefficients(vector< RooProduct* > funcs, vector< double > maxVal)
{

  vector< double > asymm;

  vector< int > numP;
  vector< int > numM;
  vector< int > denP;
  vector< int > denM;

  int nFunc = funcs.size();
  if (nFunc<1 || maxVal.size()!=nFunc) return asymm;
  
  int iFunc;
  for (iFunc=0; iFunc<nFunc; ++iFunc) {
    numP.push_back(0);
    numM.push_back(0);
    denP.push_back(0);
    denM.push_back(0);
  }

  double fr, fVal;
  TRandom * RG = RooRandom::randomGenerator();

  for (int i=0; i<1e6; ++i) {

    vars.setRealValue( "ctK", RooRandom::uniform(RG)*2-1 );
    vars.setRealValue( "ctL", RooRandom::uniform(RG)*2-1 );
    vars.setRealValue( "phi", (RooRandom::uniform(RG)*2-1)*TMath::Pi() );
    fr = RooRandom::uniform(RG);

    for (iFunc=0; iFunc<nFunc; ++iFunc) {
      fVal = funcs[iFunc]->getVal();
      if ( fVal > 0 ) {
	++denP[iFunc];
	if ( fr * maxVal[iFunc] < fVal ) ++numP[iFunc];
      } else if ( fVal < 0 ) {
	++denM[iFunc];
	if ( -1*fr * maxVal[iFunc] > fVal ) ++numM[iFunc];
      }
    }

  }

  for (iFunc=0; iFunc<nFunc; ++iFunc) asymm.push_back( denP[iFunc]>0 && denM[iFunc]>0 ?
						       ( 1.0*numP[iFunc]/denP[iFunc] + 1.0*numM[iFunc]/denM[iFunc] ) * maxVal[iFunc] : 
						       1.0 );

  return asymm;
}

int doubleFactorial (int i)
{
  int result = 1;
  for (int n=i; n>1; n=n-2) result = result*n;
  return result;
}
