#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>

#include "PdfRT.h"
#include "PdfWT.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* c [2*nBins];

void fit_recoMC_projEffBin(int q2Bin, bool tagFlag, int maxOrder, int xbins, int ybins, int zbins)
{

  string shortString = Form(tagFlag?"b%ict_JpsiCh":"b%iwt_JpsiCh",q2Bin);
  string longString  = Form(tagFlag?"Jpsi q2 bin %i correct-tag":"Jpsi q2 bin %i wrong-tag",q2Bin);
  int confIndex = (tagFlag?q2Bin:q2Bin+nBins);

  // open file with efficiency and import efficiency coefficients and variables
  TFile* fin = new TFile( ( Form("effProjection_sh%io_",maxOrder)+shortString+Form("_%i_%i_%i.root",xbins,ybins,zbins)).c_str(), "READ" );
  if ( !fin || !fin->IsOpen() ) {
    cout<<Form("File not found: effProjection_sh%io_",maxOrder)<<shortString<<Form("_%i_%i_%i.root",xbins,ybins,zbins)<<endl;
    return;
  }
  RooWorkspace* ws = (RooWorkspace*)fin->Get("ws");
  if ( !ws || ws->IsZombie() ) {
    cout<<Form("Workspace not found in file: effProjection_sh%io_",maxOrder)<<shortString<<Form("_%i_%i_%i.root",xbins,ybins,zbins)<<endl;
    return;
  } 

  // import the efficiency coefficients
  RooArgList* EffCoeff = new RooArgList("EffCoeff");
  int k_ord, l_ord, m_ord;
  for (k_ord=0; k_ord<=maxOrder; ++k_ord)
    for (l_ord=0; l_ord<=maxOrder; ++l_ord)
      for (m_ord=-1*TMath::Min(k_ord,l_ord); m_ord<=TMath::Min(k_ord,l_ord); ++m_ord)
	EffCoeff->add(*ws->var(Form("l%i_k%i_m%i",k_ord,l_ord,m_ord)));

  RooRealVar* ctK = ws->var("ctK");
  RooRealVar* ctL = ws->var("ctL");
  RooRealVar* phi = ws->var("phi");
  RooArgSet vars (* ctK,* ctL,* phi);

  RooRealVar* Fl = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1 = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2 = new RooRealVar("P2","P_{2}",0,-1,1);
  RooRealVar* P3 = new RooRealVar("P3","P_{3}",0,-1,1);
  RooRealVar* P4p = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p = new RooRealVar("P6p","P'_{6}",0,-1,1);
  RooRealVar* P8p = new RooRealVar("P8p","P'_{8}",0,-1,1);

  RooAbsPdf* AnglesPDF = 0;
  if (tagFlag) AnglesPDF = new PdfRT("AnglesPDF","AnglesPDF",*ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*EffCoeff,maxOrder);
  else         AnglesPDF = new PdfWT("AnglesPDF","AnglesPDF",*ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*EffCoeff,maxOrder);

  // Load ntuples
  TChain* tChain = new TChain();
  tChain->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/2016MC_RECO_p1p2_newtag_JPsi_add4BDT_addvars_bestBDTv4.root/ntuple");
  int numEntries = tChain->GetEntries();

  double recoCosThetaK, recoCosThetaL, recoPhi;
  float recoDimuMass, recoB0pT, recoB0eta, genSignal, tagB0;
  tChain->SetBranchAddress( "cos_theta_k" , &recoCosThetaK );
  tChain->SetBranchAddress( "cos_theta_l" , &recoCosThetaL );
  tChain->SetBranchAddress( "phi_kst_mumu", &recoPhi       );
  tChain->SetBranchAddress( "mumuMass"    , &recoDimuMass  );
  tChain->SetBranchAddress( "bPt"         , &recoB0pT      );
  tChain->SetBranchAddress( "bEta"        , &recoB0eta     );
  tChain->SetBranchAddress( "genSignal"   , &genSignal     );
  tChain->SetBranchAddress( "tagB0"       , &tagB0         );

  RooDataSet* data = new RooDataSet( "data", "RECO distribution after selections", vars ); 

  int counter=0;
  for (int iCand=0; iCand<numEntries; ++iCand) {
    tChain->GetEntry(iCand);
    // selct q2 range and tag status
    if ( ( pow(recoDimuMass,2) > binBorders[q2Bin+1] ) ||
  	 ( pow(recoDimuMass,2) < binBorders[q2Bin]   ) || 
  	 ( ( tagFlag) && (genSignal == tagB0+3) ) ||
  	 ( (!tagFlag) && (genSignal != tagB0+3) ) ) continue;
    // status display
    if ( iCand > 1.0*counter*numEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // fill
    ctK->setVal(recoCosThetaK);
    ctL->setVal(recoCosThetaL);
    phi->setVal(recoPhi);
    data->add(vars);    
  }

  RooFitResult * fitResult = AnglesPDF->fitTo(*data,Minimizer("Minuit2","migrad"),Save(true),Timer(true)); 
  fitResult->Print("v");

  TFile f (("fitResult_recoMC_"+shortString+Form("_LegTest_%i_%i_%i_sh%io.root",xbins,ybins,zbins,maxOrder)).c_str(),"UPDATE") ;
  f.cd();
  fitResult->Write("fitResult");
  f.Close();

  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),2000,700);
  TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
  RooPlot* xframe = ctK->frame(Title((longString+" - cos(#theta_{K}) distribution").c_str()));
  RooPlot* yframe = ctL->frame(Title((longString+" - cos(#theta_{L}) distribution").c_str()));
  RooPlot* zframe = phi->frame(Title((longString+" - #phi distribution").c_str()));
  data->plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40),Name("plData"));
  data->plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40));
  data->plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40));
  AnglesPDF->plotOn(xframe,LineWidth(1),Name("plPDF"));
  AnglesPDF->plotOn(yframe,LineWidth(1));
  AnglesPDF->plotOn(zframe,LineWidth(1));
  xframe->GetYaxis()->SetTitleOffset(1.8);
  yframe->GetYaxis()->SetTitleOffset(1.8);
  zframe->GetYaxis()->SetTitleOffset(1.8);
  xframe->SetMaximum(xframe->GetMaximum()*1.15);
  yframe->SetMaximum(yframe->GetMaximum()*1.15);
  zframe->SetMaximum(zframe->GetMaximum()*1.15);
  xframe->SetMinimum(0);
  yframe->SetMinimum(0);
  zframe->SetMinimum(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(xframe->findObject("plData"),"Post-selection distribution" ,"lep");
  leg->AddEntry(xframe->findObject("plPDF" ),Form("Decay rate x efficiency (%ith order)",maxOrder) ,"l");

  c[confIndex]->Divide(3,1);
  c[confIndex]->cd(1);
  gPad->SetLeftMargin(0.19); 
  xframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(2);
  gPad->SetLeftMargin(0.19); 
  yframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(3);
  gPad->SetLeftMargin(0.19); 
  zframe->Draw();
  leg->Draw("same");

  c[confIndex]->SaveAs( ("fitResult_recoMC_"+shortString+Form("_LegTest_%i_%i_%i_sh%io.pdf",xbins,ybins,zbins,maxOrder)).c_str() );

}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin, [-1] for each bin recursively
  // tag format:    [1] correctly-tagged. [0] mistagged, [2] each tag, recursively

  if ( argc < 4) return 1;
  int q2Bin    = atoi(argv[1]);
  int tagFlag  = atoi(argv[2]);
  int maxOrder = atoi(argv[3]);

  if ( maxOrder < 0 ) return 1;

  int xbins = 25;
  int ybins = 0;
  int zbins = 0;

  if ( argc >= 5)  xbins = atoi(argv[4]);
  if ( argc >= 6)  ybins = atoi(argv[5]);
  if ( argc >= 7)  zbins = atoi(argv[6]);

  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return 1;

  if ( q2Bin > -1 && q2Bin < nBins ) {
    if (tagFlag < 2 && tagFlag > -1) {
      cout<<"Fitting post-selection distributions for q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
      fit_recoMC_projEffBin(q2Bin, (tagFlag==1), maxOrder, xbins, ybins, zbins);
    }
    if (tagFlag == 2) {
      cout<<"Fitting post-selection distributions for q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
      fit_recoMC_projEffBin(q2Bin, true,  maxOrder, xbins, ybins, zbins);
      fit_recoMC_projEffBin(q2Bin, false, maxOrder, xbins, ybins, zbins);
    }
  }
  if (q2Bin == -1) {
    cout<<"Fitting post-selection distributions for all q2 bins"<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if (tagFlag < 2 && tagFlag > -1) {
	cout<<endl<<"q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
	fit_recoMC_projEffBin(q2Bin, (tagFlag==1), maxOrder, xbins, ybins, zbins);
      }
      if (tagFlag == 2) {
	cout<<endl<<"q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
	fit_recoMC_projEffBin(q2Bin, true,  maxOrder, xbins, ybins, zbins);
	fit_recoMC_projEffBin(q2Bin, false, maxOrder, xbins, ybins, zbins);
      }
    }
  }

  return 0;

}

