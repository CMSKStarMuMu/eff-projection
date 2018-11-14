#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"

using namespace RooFit ;
using namespace std ;

// want to produce and save plots of distributions and efficiency?
bool plot = true;

static const int nBins = 9;

TCanvas* cx1 [2*nBins];
TCanvas* cy1 [2*nBins];
TCanvas* cz1 [2*nBins];

void projEff_spHarm_fromHistBin(int q2Bin, bool tagFlag, int maxOrder, int xbins, int ybins, int zbins);

RooRealVar ctK ("ctK","cos(#theta_{K})",-1,1);
RooRealVar ctL ("ctL","cos(#theta_{L})",-1,1);
RooRealVar phi ("phi","#phi",-TMath::Pi(),TMath::Pi());
RooArgSet vars (ctK, ctL, phi);

void projEff_spHarm_fromHist(int q2Bin, int tagFlag, int maxOrder = 5, int xbins=25, int ybins = 0, int zbins = 0)
{
  // q2-bin format: [0-8] for one bin, [-1] for each bin recursively
  // tag format:    [1] correctly-tagged. [0] mistagged, [2] each tag, recursively

  if ( maxOrder < 0 ) return;

  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;

  if ( q2Bin > -1 && q2Bin < nBins ) {
    if (tagFlag < 2 && tagFlag > -1) {
      cout<<"Projecting efficiency for q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
      projEff_spHarm_fromHistBin(q2Bin, (tagFlag==1), maxOrder, xbins, ybins, zbins);
    }
    if (tagFlag == 2) {
      cout<<"Projecting efficiency for q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
      projEff_spHarm_fromHistBin(q2Bin, true,  maxOrder, xbins, ybins, zbins);
      projEff_spHarm_fromHistBin(q2Bin, false, maxOrder, xbins, ybins, zbins);
    }
  }
  if (q2Bin == -1) {
    cout<<"Projecting efficiency for all q2 bins"<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if (tagFlag < 2 && tagFlag > -1) {
	cout<<endl<<"q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
	projEff_spHarm_fromHistBin(q2Bin, (tagFlag==1), maxOrder, xbins, ybins, zbins);
      }
      if (tagFlag == 2) {
	cout<<endl<<"q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
	projEff_spHarm_fromHistBin(q2Bin, true,  maxOrder, xbins, ybins, zbins);
	projEff_spHarm_fromHistBin(q2Bin, false, maxOrder, xbins, ybins, zbins);
      }
    }
  }
  
}

void projEff_spHarm_fromHistBin(int q2Bin, bool tagFlag, int maxOrder, int xbins, int ybins, int zbins)
{

  string shortString = Form(tagFlag?"b%ict_JpsiCh":"b%iwt_JpsiCh",q2Bin);
  string longString  = Form(tagFlag?"Jpsi q2 bin %i correct-tag":"Jpsi q2 bin %i wrong-tag",q2Bin);
  int confIndex = (tagFlag?q2Bin:q2Bin+nBins);

  // Load histograms
  TFile* fin = TFile::Open( ("effHist_"+shortString+Form("_%i_%i_%i.root",xbins,ybins,zbins)).c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: effHist_"+shortString+Form("_%i_%i_%i.root",xbins,ybins,zbins)<<endl;
    return;
  }
  TH3F* denHist = (TH3F*)fin->Get(("denHist"+shortString+"__ctK_ctL_phi").c_str());
  TH3F* numHist = (TH3F*)fin->Get(("numHist"+shortString+"__ctK_ctL_phi").c_str());

  //Declare and initialise the functions and all the needed objects
  vector < RooRealVar* > factors;
  vector < double > proj;
  vector < RooLegendre* > vectFuncLegCosK;
  vector < RooLegendre* > vectFuncLegCosL;
  vector < RooFormulaVar* > vectFuncPoly;
  vector < RooProduct* > vectFunc;

  RooArgList facList;
  RooArgList funList;

  for (int xOrder=0; xOrder<=maxOrder; ++xOrder)
    for (int yOrder=0; yOrder<=maxOrder; ++yOrder)
      for (int zOrder=-1*TMath::Min(xOrder,yOrder); zOrder<=TMath::Min(xOrder,yOrder); ++zOrder) {

	// vector of coefficients for the function basis
	factors.push_back( new RooRealVar( Form("l%i_k%i_m%i",xOrder,yOrder,zOrder),
					   Form("l%i_k%i_m%i",xOrder,yOrder,zOrder), 0 ) );

	RooArgList prodList;

	// phi terms by trigonometric polynomials (degree zOrder)
	if (zOrder>0) {
	  vectFuncPoly.push_back( new RooFormulaVar( Form("funcPoly%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("funcPoly%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("cos(%i*phi)",zOrder), phi ) );
	  prodList.add( *vectFuncPoly.back() );
	}
	if (zOrder<0) {
	  vectFuncPoly.push_back( new RooFormulaVar( Form("funcPoly%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("funcPoly%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("sin(%i*phi)",-1*zOrder), phi ) );
	  prodList.add( *vectFuncPoly.back() );
	}

	// cosTK terms by associated Legendre polynomials (degree l=xOrder m=zOrder)
	vectFuncLegCosK.push_back( new RooLegendre ( Form("funcLegctK%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("funcLegctK%i_%i_%i",xOrder,yOrder,zOrder),
						     ctK, xOrder, abs(zOrder) ) );
	prodList.add( *vectFuncLegCosK.back() );

	// cosTL terms by associated Legendre polynomials (degree l=yOrder m=zOrder)
	vectFuncLegCosL.push_back( new RooLegendre ( Form("funcLegctL%i_%i_%i",xOrder,yOrder,zOrder),
						     Form("funcLegctL%i_%i_%i",xOrder,yOrder,zOrder),
						     ctL, yOrder, abs(zOrder) ) );
	prodList.add( *vectFuncLegCosL.back() );

	// build member of the basis of 3D functions
	vectFunc.push_back( new RooProduct ( Form("func%i_%i_%i",xOrder,yOrder,zOrder),
					     Form("func%i_%i_%i",xOrder,yOrder,zOrder),
					     prodList ) );

	// coefficients values to be filled later
	proj.push_back(0);
	
	// preparation of RooArgList objects
	funList.add( *vectFunc.back() );
	facList.add( *factors .back() );

      }

  cout<<"Number of parameters used: "<<factors.size()<<endl;

  // Sum function
  RooAddition projectedFunc ( "projectedFunc", "projectedFunc", funList, facList );

  //Compute and set the coefficients
  double fact;
  int iOrder=-1;

  TStopwatch t;
  t.Start();

  // loop over the coefficients
  for (int xOrder=0; xOrder<=maxOrder; ++xOrder)
    for (int yOrder=0; yOrder<=maxOrder; ++yOrder)
      for (int zOrder=-1*TMath::Min(xOrder,yOrder); zOrder<=TMath::Min(xOrder,yOrder); ++zOrder) {
	
	++iOrder;

	// project the binned efficiency on the [iOrder] function
	for (int xBin=1; xBin<=denHist->GetNbinsX(); ++xBin)
	  for (int yBin=1; yBin<=denHist->GetNbinsY(); ++yBin)
	    for (int zBin=1; zBin<=denHist->GetNbinsZ(); ++zBin) {
	      
	      ctK.setVal( denHist->GetXaxis()->GetBinCenter(xBin) );
	      ctL.setVal( denHist->GetYaxis()->GetBinCenter(yBin) );
	      phi.setVal( denHist->GetZaxis()->GetBinCenter(zBin) );

	      // contribution of one bin
 	      if ( denHist->GetBinContent(xBin,yBin,zBin)>0 )
		proj[iOrder] += ( numHist->GetBinContent(xBin,yBin,zBin) / denHist->GetBinContent(xBin,yBin,zBin) *
				  denHist->GetXaxis()->GetBinWidth(xBin) *
				  denHist->GetYaxis()->GetBinWidth(yBin) *
				  denHist->GetZaxis()->GetBinWidth(zBin) *
				  vectFunc[iOrder]->getVal( vars ) );

	    }

	// normalization of 0-degree trigonometric polynomial differs by a factor 2
	if (zOrder==0) proj[iOrder] = proj[iOrder]/2.0;
	
	// set coefficient value, normalised
	factors[iOrder]->setVal( proj[iOrder]
				 * (2*xOrder+1)*TMath::Factorial(xOrder-abs(zOrder))/2/TMath::Factorial(xOrder+abs(zOrder)) // associated legendre poly
				 * (2*yOrder+1)*TMath::Factorial(yOrder-abs(zOrder))/2/TMath::Factorial(yOrder+abs(zOrder))
				 / TMath::Pi() // trigonometric polynomial
				 );

	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<"\t"<<iOrder<<"\t"<<factors[iOrder]->getValV()<<endl;

      }

  t.Stop();
  t.Print();    

  // save efficiency function to file
  RooWorkspace *ws = new RooWorkspace("ws","Workspace with efficiency parameterisation");
  ws->import( projectedFunc, Silence() );
  ws->writeToFile( ( Form("effProjection_sh%io_",maxOrder)+shortString+Form("_%i_%i_%i.root",xbins,ybins,zbins)).c_str() );

  if (plot) {
    // Plot 1D slices of the efficiency function and binned efficiency
    vector <TEfficiency*> effHistsX; 
    vector <TEfficiency*> effHistsY;
    vector <TEfficiency*> effHistsZ;
    vector <RooPlot*> xframes;
    vector <RooPlot*> yframes;
    vector <RooPlot*> zframes;
    cx1[confIndex] = new TCanvas(("cx1"+shortString).c_str(),("Projected efficiency - "+longString+" - cos(theta_k) slices").c_str(),1500,1500) ;
    cy1[confIndex] = new TCanvas(("cy1"+shortString).c_str(),("Projected efficiency - "+longString+" - cos(theta_l) slices").c_str(),1500,1500) ;
    cz1[confIndex] = new TCanvas(("cz1"+shortString).c_str(),("Projected efficiency - "+longString+" - phi slices"         ).c_str(),1500,1500) ;
    cx1[confIndex]->Divide(5,5);
    cy1[confIndex]->Divide(5,5);
    cz1[confIndex]->Divide(5,5);

    TLegend* leg = new TLegend (0.35,0.8,0.9,0.9);

    // width of the slices in the hidden variables ("border" is half of it)
    double border = 0.035;
    // variables to be filled with global efficiency maximum
    double maxEffX = 0;
    double maxEffY = 0;
    double maxEffZ = 0;

    // loop over slice grid
    for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {

	// central values and borders of the slices in the hidden variables
	double centA = -0.8 + 1.6*i/4;
	double centB = -0.8 + 1.6*j/4;
	double lowA  = TMath::Max( centA - border,  1e-4-1 );
	double lowB  = TMath::Max( centB - border,  1e-4-1 );
	double highA = TMath::Min( centA + border, -1e-4+1 );
	double highB = TMath::Min( centB + border, -1e-4+1 );

	// slicing num and den distributions    
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

	// producing 1D efficiencies from the slices
	effHistsX.push_back( new TEfficiency(*numProjX,*denProjX) );
	effHistsX.back()->SetName( Form("effHistX_%i_%i",i,j) );
	effHistsX.back()->SetTitle( ("Efficiency - "+longString+Form(" - slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Efficiency",centA,centB*TMath::Pi())).c_str() );
    
	effHistsY.push_back( new TEfficiency(*numProjY,*denProjY) );
	effHistsY.back()->SetName( Form("effHistY_%i_%i",i,j) );
	effHistsY.back()->SetTitle( ("Efficiency - "+longString+Form(" - slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Efficiency",centA,centB*TMath::Pi())).c_str() );

	effHistsZ.push_back( new TEfficiency(*numProjZ,*denProjZ) );
	effHistsZ.back()->SetName( Form("effHistZ_%i_%i",i,j) );
	effHistsZ.back()->SetTitle( ("Efficiency - "+longString+Form(" - slice ctK=%1.2f ctL=%1.2f;#phi;Efficiency",centA,centB)).c_str() );

	// producing 1D slices of efficiency description
	ctL.setVal(centA);
	phi.setVal(centB*TMath::Pi());
	xframes.push_back( ctK.frame(Name( (Form("frameslice_ctk_%i_%i_",i,j)+shortString).c_str() )) );
	projectedFunc.plotOn(xframes.back(),LineColor(kRed),Name(Form("projectedFuncx_%i_%i",i,j))) ;

	ctK.setVal(centA);
	phi.setVal(centB*TMath::Pi());
	yframes.push_back( ctL.frame(Name( (Form("frameslice_ctl_%i_%i_",i,j)+shortString).c_str() )) );
	projectedFunc.plotOn(yframes.back(),LineColor(kRed),Name(Form("projectedFuncy_%i_%i",i,j))) ;

	ctK.setVal(centA);
	ctL.setVal(centB);
	zframes.push_back( phi.frame(Name( (Form("frameslice_phi_%i_%i_",i,j)+shortString).c_str() )) );
	projectedFunc.plotOn(zframes.back(),LineColor(kRed),Name(Form("projectedFuncz_%i_%i",i,j))) ;

	// plot in canvas
	cx1[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effHistsX.back()->Draw();
	cx1[confIndex]->cd(5*j+i+1)->Update(); 
	auto graphx = effHistsX.back()->GetPaintedGraph(); 
	graphx->SetMinimum(0);
	auto effValsX = graphx->GetY();
	for (int iBin=0; iBin<graphx->GetN(); ++iBin) if (maxEffX<effValsX[iBin]) maxEffX = effValsX[iBin];
	// if (maxEffX<graphx->GetYaxis()->GetXmax()) maxEffX = graphx->GetYaxis()->GetXmax();
	graphx->GetYaxis()->SetTitleOffset(1.7);
	cx1[confIndex]->cd(5*j+i+1)->Update();
	xframes.back()->Draw("same") ;

	if (i+j==0) {
	  leg->AddEntry(effHistsX.back(),"Binned efficiency" ,"lep");
	  leg->AddEntry(xframes.back()->findObject(Form("projectedFuncx_%i_%i",i,j)),"Projected spherical harmonics","l");
	}

	leg->Draw("same");

	cy1[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effHistsY.back()->Draw();
	cy1[confIndex]->cd(5*j+i+1)->Update(); 
	auto graphy = effHistsY.back()->GetPaintedGraph(); 
	graphy->SetMinimum(0);
	auto effValsY = graphy->GetY();
	for (int iBin=0; iBin<graphy->GetN(); ++iBin) if (maxEffY<effValsY[iBin]) maxEffY = effValsY[iBin];
	// if (maxEffY<graphy->GetYaxis()->GetXmax()) maxEffY = graphy->GetYaxis()->GetXmax();
	graphy->GetYaxis()->SetTitleOffset(1.7);
	cy1[confIndex]->cd(5*j+i+1)->Update();
	yframes.back()->Draw("same") ;
	leg->Draw("same");

	cz1[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effHistsZ.back()->Draw();
	cz1[confIndex]->cd(5*j+i+1)->Update(); 
	auto graphz = effHistsZ.back()->GetPaintedGraph(); 
	graphz->SetMinimum(0);
	auto effValsZ = graphz->GetY();
	for (int iBin=0; iBin<graphz->GetN(); ++iBin) if (maxEffZ<effValsZ[iBin]) maxEffZ = effValsZ[iBin];
	// if (maxEffZ<graphz->GetYaxis()->GetXmax()) maxEffZ = graphz->GetYaxis()->GetXmax();
	graphz->GetYaxis()->SetTitleOffset(1.7);
	cz1[confIndex]->cd(5*j+i+1)->Update();
	zframes.back()->Draw("same") ;
	leg->Draw("same");

      }    

    // set uniform y-axis ranges
    for (int i=0; i<effHistsX.size(); ++i) (effHistsX[i]->GetPaintedGraph())->SetMaximum(maxEffX*1.25);
    for (int i=0; i<effHistsY.size(); ++i) (effHistsY[i]->GetPaintedGraph())->SetMaximum(maxEffY*1.25);
    for (int i=0; i<effHistsZ.size(); ++i) (effHistsZ[i]->GetPaintedGraph())->SetMaximum(maxEffZ*1.25);
    for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {
	cx1[confIndex]->cd(5*j+i+1)->Update();
	cy1[confIndex]->cd(5*j+i+1)->Update();
	cz1[confIndex]->cd(5*j+i+1)->Update();
      }	

    cx1[confIndex]->SaveAs( ("effProj_"+shortString+Form("_%i_%i_%i_CTKslices_sh%io_dp%i.pdf",xbins,ybins,zbins,maxOrder,(int)(border*200))).c_str() );
    cy1[confIndex]->SaveAs( ("effProj_"+shortString+Form("_%i_%i_%i_CTLslices_sh%io_dp%i.pdf",xbins,ybins,zbins,maxOrder,(int)(border*200))).c_str() );
    cz1[confIndex]->SaveAs( ("effProj_"+shortString+Form("_%i_%i_%i_PHIslices_sh%io_dp%i.pdf",xbins,ybins,zbins,maxOrder,(int)(border*200))).c_str() );
  }
  
}
