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

// string plotName = "tentativeEff_b3rt_6o-sh";
string plotName = "kernelTestBin25_fixv1_b3rt_8o-sh";

int maxOrder = 8;
int q2Bin = 3;

bool effExport = false;

int xbins = 25;
int ybins = 25;
int zbins = 25;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

void fillHists(TH3F* denHist, TH3F* numHist, RooDataSet* data, RooDataSet* numData, int nev = 1000);

RooRealVar ctK ("ctK","cos(theta_k)",-1,1);
RooRealVar ctL ("ctL","cos(theta_L)",-1,1);
RooRealVar phi ("phi","phi",-TMath::Pi(),TMath::Pi());
RooArgSet vars (ctK, ctL, phi);

void projEff_spHarm_kernel()
{

  if (q2Bin >= nBins) return;

  // Load ntuples
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();

  t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN/2016MC_GEN_LMNR_double_sub*_p*.root/ntuple");
  t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/2016MC_RECO_p1p2p3_newtag_LMNR_addW_add4BDT_addvars_bestBDTv4.root/ntuple");

  // // Prepare denominator datasets
  // float genCosThetaK, genCosThetaL, genPhi, genDimuMass, genB0pT, genB0eta;
  // t_den->SetBranchAddress( "gen_cos_theta_k", &genCosThetaK );
  // t_den->SetBranchAddress( "gen_cos_theta_l", &genCosThetaL );
  // t_den->SetBranchAddress( "gen_phi"        , &genPhi       );
  // t_den->SetBranchAddress( "mumu_mass"      , &genDimuMass  );
  // t_den->SetBranchAddress( "b0_pt"          , &genB0pT      );
  // t_den->SetBranchAddress( "b0_eta"         , &genB0eta     );
  // RooDataSet* data = new RooDataSet( "data", "GEN distribution before GEN-filter", vars );
  // cout<<"Starting denominator dataset filling..."<<endl;
  // int counter=0;
  // for (int iCand=0; iCand<t_den->GetEntries(); ++iCand) {
  //   if ( iCand > 1.0*counter*t_den->GetEntries()/100 ) {
  //     cout<<counter<<"%"<<endl;
  //     counter += 10;
  //   }
  //   t_den->GetEntry(iCand);
  //   if ( genDimuMass*genDimuMass<binBorders[q2Bin] || genDimuMass*genDimuMass>binBorders[q2Bin+1] ) continue;
  //   // if ( fabs(genB0eta)>2.4 || genDimuMass*genDimuMass<binBorders[q2Bin] || genDimuMass*genDimuMass>binBorders[q2Bin+1] ) continue;
  //   // if ( genB0pT < 8 || fabs(genB0eta)>2.2 || genDimuMass*genDimuMass<binBorders[q2Bin] || genDimuMass*genDimuMass>binBorders[q2Bin+1] ) continue;
  //   ctK.setVal(genCosThetaK);
  //   ctL.setVal(genCosThetaL);
  //   phi.setVal(genPhi);
  //   data->add(vars);    
  // }
  // Prepare denominator datasets
  double genCosThetaK, genCosThetaL, genPhi, genDimuMass, genB0pT, genB0eta;
  t_den->SetBranchAddress( "gen_cos_theta_k" , &genCosThetaK );
  t_den->SetBranchAddress( "gen_cos_theta_l" , &genCosThetaL );
  t_den->SetBranchAddress( "gen_phi_kst_mumu", &genPhi       );
  t_den->SetBranchAddress( "genq2"           , &genDimuMass  );
  t_den->SetBranchAddress( "genbPt"          , &genB0pT      );
  t_den->SetBranchAddress( "genbEta"         , &genB0eta     );
  RooDataSet* data = new RooDataSet( "data", "GEN distribution before GEN-filter", vars );
  cout<<"Starting denominator dataset filling..."<<endl;
  int counter=0;
  for (int iCand=0; iCand<t_den->GetEntries(); ++iCand) {
    if ( iCand > 1.0*counter*t_den->GetEntries()/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    t_den->GetEntry(iCand);
    if ( genDimuMass*genDimuMass<binBorders[q2Bin] || genDimuMass*genDimuMass>binBorders[q2Bin+1] ) continue;
    // if ( fabs(genB0eta)>2.4 || genDimuMass*genDimuMass<binBorders[q2Bin] || genDimuMass*genDimuMass>binBorders[q2Bin+1] ) continue;
    // if ( genB0pT < 8 || fabs(genB0eta)>2.2 || genDimuMass*genDimuMass<binBorders[q2Bin] || genDimuMass*genDimuMass>binBorders[q2Bin+1] ) continue;
    ctK.setVal(genCosThetaK);
    ctL.setVal(genCosThetaL);
    phi.setVal(genPhi);
    data->add(vars);    
  }

  // Prepare numerator datasets
  double recoCosThetaK, recoCosThetaL, recoPhi;
  float recoDimuMass, recoB0pT, recoB0eta, genSignal, tagB0;
  t_num->SetBranchAddress( "cos_theta_k" , &recoCosThetaK );
  t_num->SetBranchAddress( "cos_theta_l" , &recoCosThetaL );
  t_num->SetBranchAddress( "phi_kst_mumu", &recoPhi       );
  t_num->SetBranchAddress( "mumuMass"    , &recoDimuMass  );
  t_num->SetBranchAddress( "bPt"         , &recoB0pT      );
  t_num->SetBranchAddress( "bEta"        , &recoB0eta     );
  t_num->SetBranchAddress( "genSignal"   , &genSignal     );
  t_num->SetBranchAddress( "tagB0"       , &tagB0     );
  RooDataSet* numData = new RooDataSet( "numData", "RECO distribution after selections", vars ); 
  cout<<"Starting numerator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<t_num->GetEntries(); ++iCand) {
    if ( iCand > 1.0*counter*t_num->GetEntries()/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    t_num->GetEntry(iCand);
    if ( recoDimuMass*recoDimuMass<binBorders[q2Bin] || recoDimuMass*recoDimuMass>binBorders[q2Bin+1] || genSignal!=tagB0+1) continue;
    // if ( fabs(recoB0eta)>2.4 || recoDimuMass*recoDimuMass<binBorders[q2Bin] || recoDimuMass*recoDimuMass>binBorders[q2Bin+1] || genSignal!=tagB0+1) continue;
    // if ( recoB0pT < 8 || fabs(recoB0eta)>2.2 || recoDimuMass*recoDimuMass<binBorders[q2Bin] || recoDimuMass*recoDimuMass>binBorders[q2Bin+1] || genSignal!=tagB0+1) continue;
    ctK.setVal(recoCosThetaK);
    ctL.setVal(recoCosThetaL);
    phi.setVal(recoPhi);
    numData->add(vars);    
  }
  // RooRealVar B0pT     ("bPt"      ,"B0pT"    ,0);
  // RooRealVar B0Eta    ("bEta"     ,"B0Eta"   ,0);
  // RooRealVar DimuMass ("mumuMass" ,"DimuMass",0);
  // RooRealVar recoTag  ("tagB0"    ,"recoTag" ,0);
  // RooRealVar genTag   ("genSignal","genTag"  ,1);
  // RooFormulaVar select("select","select",Form("bPt>8 && fabs(bEta)<2.2 && mumuMass*mumuMass>%f && mumuMass*mumuMass>%f && genSignal==tagB0+1",binBorders[q2Bin],binBorders[q2Bin+1]),
  // 		       RooArgList(B0pT,B0Eta,DimuMass,recoTag,genTag));

  cout<<"Done!"<<endl;

  double avgEff = numData->sumEntries() / data->sumEntries();
  cout<<"Average efficiency = "<<avgEff<<endl;

  // Plot numerator and denominator datasets
  RooPlot* xframe = ctK.frame(Title("Numerator and denominator cos(theta_K) distributions")) ;
  data->plotOn(xframe,Rescale(0.1),MarkerColor(kRed+1),Binning(30));
  numData->plotOn(xframe,Binning(30));
  RooPlot* yframe = ctL.frame(Title("Numerator and denominator cos(theta_L) distributions")) ;
  data->plotOn(yframe,Rescale(0.1),MarkerColor(kRed+1),Binning(30));
  numData->plotOn(yframe,Binning(30));
  RooPlot* zframe = phi.frame(Title("Numerator and denominator phi distributions")) ;
  data->plotOn(zframe,Rescale(0.1),MarkerColor(kRed+1),Binning(30));
  numData->plotOn(zframe,Binning(30));

  TCanvas* c = new TCanvas("NumDenCanvas","toy_numerator_denominator",1200,800) ;
  c->Divide(3,1);
  c->cd(1);
  gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->GetYaxis()->SetRangeUser(0,6000); xframe->Draw() ;
  c->cd(2);
  gPad->SetLeftMargin(0.15) ; yframe->GetYaxis()->SetTitleOffset(1.4) ; yframe->GetYaxis()->SetRangeUser(0,3400); yframe->Draw() ;
  c->cd(3);
  gPad->SetLeftMargin(0.15) ; zframe->GetYaxis()->SetTitleOffset(1.4) ; zframe->GetYaxis()->SetRangeUser(0,3400); zframe->Draw() ;

  //Declare and initialise the PDFs and all the needed objects
  // cout<<"Declaration"<<endl;
  vector < RooRealVar* > factors;
  vector < double > proj;
  vector < RooLegendre* > vectPdfLegCosK;
  vector < RooLegendre* > vectPdfLegCosL;
  vector < RooFormulaVar* > vectPdfPoly;
  vector < RooProduct* > vectPdf;

  RooArgList facList;
  RooArgList pdfList;

  RooPlot* xframeT = ctK.frame(Title("Numerator and denominator cos(theta_K) distributions")) ;
  RooPlot* yframeT = ctL.frame(Title("Numerator and denominator cos(theta_L) distributions")) ;
  RooPlot* zframeT = phi.frame(Title("Numerator and denominator phi distributions")) ;

  for (int xOrder=0; xOrder<=maxOrder; ++xOrder)
    for (int yOrder=0; yOrder<=maxOrder; ++yOrder)
      for (int zOrder=-1*TMath::Min(xOrder,yOrder); zOrder<=TMath::Min(xOrder,yOrder); ++zOrder) {

	factors  .push_back( new RooRealVar(Form( "l%i_k%i_m%i",xOrder,yOrder,zOrder),Form( "l%i_k%i_m%i",xOrder,yOrder,zOrder),0) );

	RooArgList prodList;

	if (zOrder>0) {
	  vectPdfPoly  .push_back( new RooFormulaVar( Form("pdfPoly%i_%i_%i"  ,xOrder,yOrder,zOrder), Form("pdfPoly%i_%i_%i"  ,xOrder,yOrder,zOrder),
						      Form("cos(%i*phi)",zOrder), phi ) );
	  prodList  .add( *vectPdfPoly  .back() );
	}
	if (zOrder<0) {
	  vectPdfPoly  .push_back( new RooFormulaVar( Form("pdfPoly%i_%i_%i"  ,xOrder,yOrder,zOrder), Form("pdfPoly%i_%i_%i"  ,xOrder,yOrder,zOrder),
						      Form("sin(%i*phi)",-1*zOrder), phi ) );
	  prodList  .add( *vectPdfPoly  .back() );
	}

	vectPdfLegCosK  .push_back( new RooLegendre ( Form("pdfLegctK%i_%i_%i"  ,xOrder,yOrder,zOrder),
						      Form("pdfLegctK%i_%i_%i"  ,xOrder,yOrder,zOrder), ctK, xOrder, abs(zOrder) ) );
	prodList  .add( *vectPdfLegCosK  .back() );

	vectPdfLegCosL  .push_back( new RooLegendre ( Form("pdfLegctL%i_%i_%i"  ,xOrder,yOrder,zOrder),
						      Form("pdfLegctL%i_%i_%i"  ,xOrder,yOrder,zOrder), ctL, yOrder, abs(zOrder) ) );
	prodList  .add( *vectPdfLegCosL  .back() );

	vectPdf  .push_back( new RooProduct ( Form("pdf%i_%i_%i"  ,xOrder,yOrder,zOrder), Form("pdf%i_%i_%i"  ,xOrder,yOrder,zOrder), prodList   ) );

	proj  .push_back(0);
	
	pdfList  .add( *vectPdf  .back() );
	facList  .add( *factors  .back() );

	if ( xOrder==yOrder ) {
	  vectPdfLegCosK.back()->plotOn(xframeT,LineColor(vectPdf.size()+10));
	  vectPdfLegCosL.back()->plotOn(yframeT,LineColor(vectPdf.size()+10));
	  if ( zOrder!=0 && abs(zOrder)==xOrder ) vectPdfPoly.back()->plotOn(zframeT,LineColor(vectPdf.size()+10));
	}
      }

  cout<<"Number of parameters used: "<<factors.size()<<endl;

  // Plot elements of the function basis

  TCanvas* cT = new TCanvas("NumDenCanvasT","Projections of function basis",1200,800) ;
  cT->Divide(3,1);
  cT->cd(1);
  gPad->SetLeftMargin(0.15) ; xframeT->GetYaxis()->SetTitleOffset(1.4) ; xframeT->Draw() ;
  cT->cd(2);
  gPad->SetLeftMargin(0.15) ; yframeT->GetYaxis()->SetTitleOffset(1.4) ; yframeT->Draw() ;
  cT->cd(3);
  gPad->SetLeftMargin(0.15) ; zframeT->GetYaxis()->SetTitleOffset(1.4) ; zframeT->Draw() ;

  // cout<<"Function declaration"<<endl;

  RooAddition projectedPdf   ("projectedPdf"  , "projectedPdf"  , pdfList  , facList  );

  //Create binned histos (only for projection method)
  // cout<<"Histo declaration"<<endl;
  TH3F* denHist = new TH3F("denHist","denHist",xbins,-1,1,ybins,-1,1,zbins,-TMath::Pi(),TMath::Pi());
  TH3F* numHist = new TH3F("numHist","numHist",xbins,-1,1,ybins,-1,1,zbins,-TMath::Pi(),TMath::Pi());

  TStopwatch t1;
  t1.Start();

  fillHists(denHist,numHist,data,numData,10000);
  denHist->Sumw2();
  numHist->Sumw2();

  t1.Stop();
  t1.Print();    

  // cout<<numHist->GetBinContent(3,3,3)<<endl;
  // cout<<denHist->GetBinContent(3,3,3)<<endl;
  // return;

  //Compute and set the coefficients
  // cout<<"Calculation"<<endl;
  // factors  [0]->setVal(avgEff);

  double fact;
  int iOrder=-1;

  int evNum = numData->sumEntries();
  int evDen =    data->sumEntries();

  TStopwatch t;
  t.Start();

  for (int xOrder=0; xOrder<=maxOrder; ++xOrder)
    for (int yOrder=0; yOrder<=maxOrder; ++yOrder)
      for (int zOrder=-1*TMath::Min(xOrder,yOrder); zOrder<=TMath::Min(xOrder,yOrder); ++zOrder) {
	
	++iOrder;

	// fill proj coefficients
	// cout<<"fill projected coefficients"<<endl;

	for (int xBin=1;xBin<=denHist->GetNbinsX();++xBin) for (int yBin=1;yBin<=denHist->GetNbinsY();++yBin) for (int zBin=1;zBin<=denHist->GetNbinsZ();++zBin) {

	      ctK.setVal( denHist->GetXaxis()->GetBinCenter(xBin) );
	      ctL.setVal( denHist->GetYaxis()->GetBinCenter(yBin) );
	      phi.setVal( denHist->GetZaxis()->GetBinCenter(zBin) );

 	      if ( denHist->GetBinContent(xBin,yBin,zBin)>0 )
		proj[iOrder] += numHist->GetBinContent(xBin,yBin,zBin) / denHist->GetBinContent(xBin,yBin,zBin) *
		  denHist->GetXaxis()->GetBinWidth(xBin) *
		  denHist->GetYaxis()->GetBinWidth(yBin) *
		  denHist->GetZaxis()->GetBinWidth(zBin) *
		  vectPdf[iOrder]->getVal( vars );

	    }
	// cout<<"-vector filled"<<endl;

	if (zOrder==0) proj[iOrder] = proj[iOrder]/2.0;
	
	// if (iOrder>0) {
	factors  [iOrder]->setVal( proj[iOrder]
				   *(2*xOrder+1)*TMath::Factorial(xOrder-abs(zOrder))/2/TMath::Factorial(xOrder+abs(zOrder))
				   *(2*yOrder+1)*TMath::Factorial(yOrder-abs(zOrder))/2/TMath::Factorial(yOrder+abs(zOrder))
				   /TMath::Pi() );
	// cout<<"Variables normalized"<<endl;
	// }

	cout<<xOrder<<" "<<yOrder<<" "<<zOrder<<"\t"<<iOrder<<"\t"<<factors[iOrder]->getValV()<<endl;
      }

  t.Stop();
  t.Print();    

  TCanvas* cx1 = new TCanvas("CanX_Proj_Func","Projected efficiency - cos(theta_k) slices",1500,1500) ;
  TCanvas* cy1 = new TCanvas("CanY_Proj_Func","Projected efficiency - cos(theta_l) slices",1500,1500) ;
  TCanvas* cz1 = new TCanvas("CanZ_Proj_Func","Projected efficiency - phi slices"         ,1500,1500) ;
  cx1->Divide(5,5);
  cy1->Divide(5,5);
  cz1->Divide(5,5);

  // double border = denHist->GetXaxis()->GetBinWidth(1)/2.1;
  double border = 0.1;

  vector <TEfficiency*> effHistsX; 
  vector <TEfficiency*> effHistsY;
  vector <TEfficiency*> effHistsZ;
  vector <RooPlot*> xframes;
  vector <RooPlot*> yframes;
  vector <RooPlot*> zframes;

  TLegend* leg = new TLegend (0.5,0.7,0.9,0.9);

  for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {
      
      // double centA = denHist->GetXaxis()->GetBinCenter(i+2);
      // double centB = denHist->GetXaxis()->GetBinCenter(j+2);
      // if (i==0) centA = denHist->GetXaxis()->GetBinCenter(2);
      // if (i==4) centA = denHist->GetXaxis()->GetBinCenter(8);
      // if (j==0) centB = denHist->GetXaxis()->GetBinCenter(2);
      // if (j==4) centB = denHist->GetXaxis()->GetBinCenter(8);
      double centA = -0.8 + 1.6*i/4;
      double centB = -0.8 + 1.6*j/4;
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
      effHistsX.back()->SetTitle(Form("Efficiency q^2 bin %i - slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Efficiency",q2Bin,centA,centB*TMath::Pi()) );
    
      effHistsY.push_back( new TEfficiency(*numProjY,*denProjY) );
      effHistsY.back()->SetName( Form("effHistY_%i_%i",i,j) );
      effHistsY.back()->SetTitle(Form("Efficiency q^2 bin %i - slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Efficiency",q2Bin,centA,centB*TMath::Pi()) );

      effHistsZ.push_back( new TEfficiency(*numProjZ,*denProjZ) );
      effHistsZ.back()->SetName( Form("effHistZ_%i_%i",i,j) );
      effHistsZ.back()->SetTitle(Form("Efficiency q^2 bin %i - slice ctK=%1.2f ctL=%1.2f;#phi;Efficiency",q2Bin,centA,centB) );

      ctL.setVal(centA);
      phi.setVal(centB*TMath::Pi());
      xframes.push_back( ctK.frame(Title( Form("Projected function - x projection %i %i",i,j))) );
      projectedPdf.  plotOn(xframes.back(),LineColor(kRed)  ,Name(Form("projectedPdfx_%i_%i"  ,i,j))) ;

      ctK.setVal(centA);
      phi.setVal(centB*TMath::Pi());
      yframes.push_back( ctL.frame(Title( Form("Projected function - y projection %i %i",i,j))) );
      projectedPdf.  plotOn(yframes.back(),LineColor(kRed)  ,Name(Form("projectedPdfy_%i_%i"  ,i,j))) ;

      ctK.setVal(centA);
      ctL.setVal(centB);
      zframes.push_back( phi.frame(Title( Form("Projected function - z projection %i %i",i,j))) );
      projectedPdf.  plotOn(zframes.back(),LineColor(kRed)  ,Name(Form("projectedPdfz_%i_%i"  ,i,j))) ;

      cx1->cd(5*j+i+1);
      effHistsX.back()->Draw();
      cx1->cd(5*j+i+1)->Update(); 
      auto graphx = effHistsX.back()->GetPaintedGraph(); 
      graphx->SetMinimum(0);
      // graphx->SetMaximum(0.2);
      graphx->GetYaxis()->SetTitleOffset(1.4);
      cx1->cd(5*j+i+1)->Update();
      xframes.back()->Draw("same") ;

      if (i+j==0) {
	leg->AddEntry(effHistsX.back()                 ,"Binned efficiency","lep");
	leg->AddEntry(Form("projectedPdfx_%i_%i"  ,i,j),"Standard projection","l");
      }

      leg->Draw("same");

      cy1->cd(5*j+i+1);
      effHistsY.back()->Draw();
      cy1->cd(5*j+i+1)->Update(); 
      auto graphy = effHistsY.back()->GetPaintedGraph(); 
      graphy->SetMinimum(0);
      // graphy->SetMaximum(0.2); 
      graphy->GetYaxis()->SetTitleOffset(1.4);
      cy1->cd(5*j+i+1)->Update();
      yframes.back()->Draw("same") ;
      leg->Draw("same");

      cz1->cd(5*j+i+1);
      effHistsZ.back()->Draw();
      cz1->cd(5*j+i+1)->Update(); 
      auto graphz = effHistsZ.back()->GetPaintedGraph(); 
      graphz->SetMinimum(0);
      // graphz->SetMaximum(0.2); 
      graphz->GetYaxis()->SetTitleOffset(1.4);
      cz1->cd(5*j+i+1)->Update();
      zframes.back()->Draw("same") ;
      leg->Draw("same");

    }    

  cx1->SaveAs( (plotName + "_x.pdf").c_str() );
  cy1->SaveAs( (plotName + "_y.pdf").c_str() );
  cz1->SaveAs( (plotName + "_z.pdf").c_str() );

  if (effExport) {

    RooWorkspace *w = new RooWorkspace("w","workspace");
    w->import(projectedPdf);
    w->writeToFile( (plotName+"_eff.root").c_str() );

  }
    
}

void fillHists(TH3F* denHist, TH3F* numHist, RooDataSet* data, RooDataSet* numData, int nev = 1000)
{
  int nBins = denHist->GetNbinsX() * denHist->GetNbinsY() * denHist->GetNbinsZ();
  double AvgRadSq = TMath::Power( 6.0*nev/TMath::Pi()/TMath::Pi()/data->sumEntries() , 2.0/3 );

  vector < TH1F* > distNum;
  vector < TH1F* > distDen;

  int iBin;

  for (iBin=0; iBin<nBins; ++iBin) {
    distNum.push_back( new TH1F(Form("distNum%i",iBin),Form("distNum%i",iBin),100,0,AvgRadSq*10) );
    distDen.push_back( new TH1F(Form("distDen%i",iBin),Form("distDen%i",iBin),100,0,AvgRadSq*10) );
  }

  double* xBinCenter = new double[denHist->GetNbinsX()];
  double* yBinCenter = new double[denHist->GetNbinsY()];
  double* zBinCenter = new double[denHist->GetNbinsZ()];
  denHist->GetXaxis()->GetCenter(xBinCenter);
  denHist->GetYaxis()->GetCenter(yBinCenter);
  denHist->GetZaxis()->GetCenter(zBinCenter);

  int iBinX, iBinY, iBinZ;
  double xVal, yVal, zVal;

  cout<<"Denominator preparation"<<endl;
  int counter=0;
  for (int iEv=0; iEv<data->sumEntries(); ++iEv) {
    if ( iEv > counter ) {
      cout<<counter*100/data->sumEntries()<<"%"<<endl;
      counter+=data->sumEntries()/10;
    }
    
    const RooArgSet *iPoint = data->get(iEv);
    xVal = iPoint->getRealValue("ctK");
    yVal = iPoint->getRealValue("ctL");
    zVal = iPoint->getRealValue("phi");

    for (iBinX=1; iBinX<=denHist->GetNbinsX(); ++iBinX)
      for(iBinY=1; iBinY<=denHist->GetNbinsY(); ++iBinY)
	for(iBinZ=1; iBinZ<=denHist->GetNbinsZ(); ++iBinZ) {
	  
	  iBin = (iBinX-1) + denHist->GetNbinsX() * ( (iBinY-1) + denHist->GetNbinsY() * (iBinZ-1) );
	  // iBin = denHist->GetBin(iBinX,iBinY,iBinZ);
	  if (iBin<0 || iBin>nBins-1) {cout<<"ERROR, bin index too high: "<<iBin<<" ("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<") vs. a max of "<<nBins<<endl; return;}
	  if ( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) > AvgRadSq*10 ) continue;

	  distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          // mirrored entries
          distDen[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]-2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distDen[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]+2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]-2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]+2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()-2,2) );
          distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()+2,2) );

	}

  }

  cout<<"Numerator preparation"<<endl;
  counter=0;
  for (int iEv=0; iEv<numData->sumEntries(); ++iEv) {
    if ( iEv > counter ) {
      cout<<counter*100/numData->sumEntries()<<"%"<<endl;
      counter+=numData->sumEntries()/10;
    }
    
    const RooArgSet *iPoint = numData->get(iEv);
    xVal = iPoint->getRealValue("ctK");
    yVal = iPoint->getRealValue("ctL");
    zVal = iPoint->getRealValue("phi");

    for (iBinX=1; iBinX<=denHist->GetNbinsX(); ++iBinX)
      for (iBinY=1; iBinY<=denHist->GetNbinsY(); ++iBinY)
	for (iBinZ=1; iBinZ<=denHist->GetNbinsZ(); ++iBinZ) {
	  
	  iBin = (iBinX-1) + denHist->GetNbinsX() * ( (iBinY-1) + denHist->GetNbinsY() * (iBinZ-1) );
	  // iBin = denHist->GetBin(iBinX,iBinY,iBinZ);
	  if (iBin<0 || iBin>nBins-1) {cout<<"ERROR, bin index too high: "<<iBin<<" ("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<") vs. a max of "<<nBins<<endl; return;}
	  if ( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) > AvgRadSq*10 ) continue;
	  
	  distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          // mirrored entries
          distNum[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]-2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distNum[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]+2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]-2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]+2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()-2,2) );
          distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()+2,2) );
	  
	}

  }

  for (iBinX=1; iBinX<=denHist->GetNbinsX(); ++iBinX)
    for (iBinY=1; iBinY<=denHist->GetNbinsY(); ++iBinY)
      for (iBinZ=1; iBinZ<=denHist->GetNbinsZ(); ++iBinZ) {

	iBin = denHist->GetBin(iBinX,iBinY,iBinZ);
	int binCount = (iBinX-1) + denHist->GetNbinsX() * ( (iBinY-1) + denHist->GetNbinsY() * (iBinZ-1) );

	int maxBin = 1;
	while ( maxBin < 100 && distDen[binCount]->Integral(1,maxBin) < nev ) ++maxBin;

	cout<<"("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<") range: 1->"<<maxBin<<
	  " num="<<distNum[binCount]->Integral(1,maxBin)<<" den="<<distDen[binCount]->Integral(1,maxBin)<<endl;

	if ( distDen[binCount]->Integral(1,maxBin) >= distNum[binCount]->Integral(1,maxBin) ) {
	  numHist->SetBinContent(iBin, distNum[binCount]->Integral(1,maxBin));
	  denHist->SetBinContent(iBin, distDen[binCount]->Integral(1,maxBin));
	} else {
	  cout<<"ERROR, numerator count exceeding denominator count: "<<
	    distNum[binCount]->Integral(1,maxBin)<<" vs. "<<distDen[binCount]->Integral(1,maxBin)<<" ("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<")"<<endl;
	  numHist->SetBinContent(iBin, distDen[binCount]->Integral(1,maxBin));
	  denHist->SetBinContent(iBin, distDen[binCount]->Integral(1,maxBin));
	}

	numHist->SetBinError(iBin, sqrt(numHist->GetBinContent(iBin)));
	denHist->SetBinError(iBin, sqrt(denHist->GetBinContent(iBin)));
	  
      }


}
