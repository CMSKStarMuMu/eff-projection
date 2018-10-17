#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

void plotEffBin(int q2Bin, bool tagFlag, int maxOrder, int xbins, int ybins, int zbins);

TCanvas* c  [2*nBins];
TCanvas* c1 [2*nBins];
TCanvas* c3 [2*nBins];

void plotEff(int q2Bin, int tagFlag, int maxOrder, int xbins=25, int ybins = 0, int zbins = 0)
{
  // q2-bin format: [0-8] for one bin, [-1] for each bin recursively
  // tag format:    [1] correctly-tagged. [0] mistagged, [2] each tag, recursively

  if ( maxOrder < 0 ) return;

  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;

  if ( q2Bin > -1 && q2Bin < nBins ) {
    if (tagFlag < 2 && tagFlag > -1) {
      cout<<"Plotting efficiency for q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
      plotEffBin(q2Bin, (tagFlag==1), maxOrder, xbins, ybins, zbins);
    }
    if (tagFlag == 2) {
      cout<<"Plotting efficiency for q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
      plotEffBin(q2Bin, true,  maxOrder, xbins, ybins, zbins);
      plotEffBin(q2Bin, false, maxOrder, xbins, ybins, zbins);
    }
  }
  if (q2Bin == -1) {
    cout<<"Plotting efficiency for all q2 bins"<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if (tagFlag < 2 && tagFlag > -1) {
	cout<<endl<<"q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
	plotEffBin(q2Bin, (tagFlag==1), maxOrder, xbins, ybins, zbins);
      }
      if (tagFlag == 2) {
	cout<<endl<<"q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
	plotEffBin(q2Bin, true,  maxOrder, xbins, ybins, zbins);
	plotEffBin(q2Bin, false, maxOrder, xbins, ybins, zbins);
      }
    }
  }
  
}

void plotEffBin(int q2Bin, bool tagFlag, int maxOrder, int xbins, int ybins, int zbins)
{

  string shortString = Form(tagFlag?"b%ict_testLowStat":"b%iwt_testLowStat",q2Bin);
  string longString  = Form(tagFlag?"q2 bin %i correct-tag":"q2 bin %i wrong-tag",q2Bin);
  int confIndex = (tagFlag?q2Bin:q2Bin+nBins);

  // open file with efficiency and import efficiency function and variables
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
  RooAbsReal* eff = ws->function("projectedFunc");
  if ( !eff || eff->IsZombie() ) {
    cout<<Form("Efficiency not found in file: effProjection_sh%io_",maxOrder)<<shortString<<Form("_%i_%i_%i.root",xbins,ybins,zbins)<<endl;
    return;
  } 
  RooRealVar* ctK = ws->var("ctK");
  RooRealVar* ctL = ws->var("ctL");
  RooRealVar* phi = ws->var("phi");
  RooArgSet vars (* ctK,* ctL,* phi);
  // set a variable for weights, to correct GEN distributions in closure test
  RooRealVar wei ("wei","weight",0,1e4);

  // Plot efficiency projections
  gStyle->SetOptStat(0);
  auto h3_xyz = (TH3D*)eff->createHistogram("ctK,ctL,phi",50,50,50);
  auto h3_x  = h3_xyz->ProjectionX();
  auto h3_y  = h3_xyz->ProjectionY();
  auto h3_z  = h3_xyz->ProjectionZ();
  auto h3_xy = h3_xyz->Project3D("xy");
  auto h3_xz = h3_xyz->Project3D("xz");
  auto h3_yz = h3_xyz->Project3D("yz");

  // 1D projections
  c1[confIndex] = new TCanvas(Form("c1-%i",confIndex),("Efficiency 1D projections - "+longString).c_str(),1800,800);
  h3_x->SetTitle(("Efficiency cos(#theta_{K}) projection - "+longString).c_str());
  h3_y->SetTitle(("Efficiency cos(#theta_{L}) projection - "+longString).c_str());
  h3_z->SetTitle(("Efficiency #phi projection - "           +longString).c_str());
  h3_x->SetLineColor(2);
  h3_y->SetLineColor(2);
  h3_z->SetLineColor(2);
  h3_x->SetMinimum(0.0);
  h3_y->SetMinimum(0.0);
  h3_z->SetMinimum(0.0);
  c1[confIndex]->Divide(3,1);
  c1[confIndex]->cd(1);
  h3_x->Draw();
  c1[confIndex]->cd(2);
  h3_y->Draw();
  c1[confIndex]->cd(3);
  h3_z->Draw();
  c1[confIndex]->SaveAs( ("effProj_"+shortString+Form("_%i_%i_%i_1DProj_sh%io.pdf",xbins,ybins,zbins,maxOrder)).c_str() );

  // 2D projections
  c3[confIndex] = new TCanvas(Form("c3-%i",confIndex),("Efficiency 2D projections - "+longString).c_str(),1800,800);
  h3_xy->SetTitle(("Efficiency cos(#theta_{K})/cos(#theta_{L}) projection - "+longString).c_str());
  h3_xz->SetTitle(("Efficiency cos(#theta_{K})/#phi projection - "           +longString).c_str());
  h3_yz->SetTitle(("Efficiency cos(#theta_{L})/#phi projection - "           +longString).c_str());
  h3_xy->GetXaxis()->SetTitleOffset(1.4);
  h3_xy->GetYaxis()->SetTitleOffset(2);
  h3_xz->GetXaxis()->SetTitleOffset(1.4);
  h3_xz->GetYaxis()->SetTitleOffset(2);
  h3_yz->GetXaxis()->SetTitleOffset(1.4);
  h3_yz->GetYaxis()->SetTitleOffset(2);
  c3[confIndex]->Divide(3,1);
  c3[confIndex]->cd(1);
  h3_xy->Draw("SURF3");
  c3[confIndex]->cd(2);
  h3_xz->Draw("SURF3");
  c3[confIndex]->cd(3);
  h3_yz->Draw("SURF3");
  c3[confIndex]->SaveAs( ("effProj_"+shortString+Form("_%i_%i_%i_2DProj_sh%io.pdf",xbins,ybins,zbins,maxOrder)).c_str() );

  // Load ntuples
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
  t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN/gen_B0_miniaodWithoutGenCuts.root/tree");
  t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/2016MC_RECO_p1p2p3_newtag_LMNR_addW_add4BDT_addvars_bestBDTv4.root/ntuple");
  int denEntries = t_den->GetEntries();
  int numEntries = t_num->GetEntries();
  int counter;

  float genCosThetaK, genCosThetaL, genPhi, genDimuMass, genB0pT, genB0eta;
  double recoCosThetaK, recoCosThetaL, recoPhi;
  float recoDimuMass, recoB0pT, recoB0eta, genSignal, tagB0;
  t_den->SetBranchAddress( "gen_cos_theta_k" , &genCosThetaK );
  t_den->SetBranchAddress( "gen_cos_theta_l" , &genCosThetaL );
  t_den->SetBranchAddress( "gen_phi_kst_mumu", &genPhi       );
  t_den->SetBranchAddress( "genq2"           , &genDimuMass  );
  t_den->SetBranchAddress( "genbPt"          , &genB0pT      );
  t_den->SetBranchAddress( "genbEta"         , &genB0eta     );
  t_num->SetBranchAddress( "cos_theta_k" , &recoCosThetaK );
  t_num->SetBranchAddress( "cos_theta_l" , &recoCosThetaL );
  t_num->SetBranchAddress( "phi_kst_mumu", &recoPhi       );
  t_num->SetBranchAddress( "mumuMass"    , &recoDimuMass  );
  t_num->SetBranchAddress( "bPt"         , &recoB0pT      );
  t_num->SetBranchAddress( "bEta"        , &recoB0eta     );
  t_num->SetBranchAddress( "genSignal"   , &genSignal     );
  t_num->SetBranchAddress( "tagB0"       , &tagB0         );

  RooDataSet* data    = new RooDataSet( "data"   , "GEN distribution after GEN-filter" , RooArgSet(vars,wei) );
  RooDataSet* numData = new RooDataSet( "numData", "RECO distribution after selections", vars ); 

  // Prepare denominator datasets
  cout<<"Starting denominator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<denEntries; ++iCand) {
    t_den->GetEntry(iCand);
    // select q2 range
    if ( ( pow(genDimuMass,2) > binBorders[q2Bin+1] ) ||
	 ( pow(genDimuMass,2) < binBorders[q2Bin]   ) ) continue;
    // status display
    if ( iCand > 1.0*counter*denEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // fill
    ctK->setVal(genCosThetaK);
    ctL->setVal(genCosThetaL);
    phi->setVal(genPhi);
    // set event weight based on local efficiency value
    wei.setVal( eff->getVal( vars ) * 1e4 );
    data->add( RooArgSet(vars,wei) );    
  }
  // create the weighted dataset for GEN events
  RooDataSet* wdata = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,wei.GetName());

  // Prepare numerator datasets
  cout<<"Starting numerator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<numEntries; ++iCand) {
    t_num->GetEntry(iCand);
    // selct q2 range and tag status
    if ( ( pow(recoDimuMass,2) > binBorders[q2Bin+1] ) ||
	 ( pow(recoDimuMass,2) < binBorders[q2Bin]   ) || 
	 ( ( tagFlag) && (genSignal == tagB0+1) ) ||
	 ( (!tagFlag) && (genSignal != tagB0+1) ) ) continue;
    // status display
    if ( iCand > 1.0*counter*numEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // fill
    ctK->setVal(recoCosThetaK);
    ctL->setVal(recoCosThetaL);
    phi->setVal(recoPhi);
    numData->add(vars);    
  }
  cout<<"Dataset prepared"<<endl;

  // Plot projections for closure test (RECO vs. eff*GEN)
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Closure test - "+longString).c_str(),2000,700);
  TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
  RooPlot* xframe = ctK->frame(Title((longString+" cos(#theta_{K}) distributions").c_str()));
  RooPlot* yframe = ctL->frame(Title((longString+" cos(#theta_{L}) distributions").c_str()));
  RooPlot* zframe = phi->frame(Title((longString+" #phi distributions").c_str()));
  wdata->plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2),Name("plDenDist"));
  wdata->plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2));
  wdata->plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2));
  numData->plotOn(xframe,Binning(30),Name("plNumDist"));
  numData->plotOn(yframe,Binning(30));
  numData->plotOn(zframe,Binning(30));
  xframe->GetYaxis()->SetTitleOffset(1.6);
  yframe->GetYaxis()->SetTitleOffset(1.6);
  zframe->GetYaxis()->SetTitleOffset(1.6);
  xframe->SetMaximum(xframe->GetMaximum()*1.15);
  yframe->SetMaximum(yframe->GetMaximum()*1.15);
  zframe->SetMaximum(zframe->GetMaximum()*1.15);
  leg->SetTextSize(0.03);
  leg->AddEntry(xframe->findObject("plNumDist"),"Post-selection RECO distribution" ,"lep");
  leg->AddEntry(xframe->findObject("plDenDist"),"Efficiency-corrected GEN distribution" ,"lep");

  c[confIndex]->Divide(3,1);
  c[confIndex]->cd(1);
  gPad->SetLeftMargin(0.17); 
  xframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(2);
  gPad->SetLeftMargin(0.17); 
  yframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(3);
  gPad->SetLeftMargin(0.17); 
  zframe->Draw();
  leg->Draw("same");

  c[confIndex]->SaveAs( ("closure_"+shortString+Form("_%i_%i_%i_sh%io.pdf",xbins,ybins,zbins,maxOrder)).c_str() );

}

