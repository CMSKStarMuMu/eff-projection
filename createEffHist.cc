// #include "RooRealVar.h"
// #include "RooDataSet.h"
// #include "TCanvas.h"
// #include "TAxis.h"
// #include "TMath.h"
// #include "RooPlot.h"

using namespace RooFit;
using namespace std;

// want to produce and save plots of distributions and efficiency?
bool plot = true;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* c   [2*nBins];
TCanvas* cx1 [2*nBins];
TCanvas* cy1 [2*nBins];
TCanvas* cz1 [2*nBins];

void createEffHistBin(int q2Bin, bool tagFlag, int xbins, int ybins, int zbins);

RooRealVar ctK ("ctK","cos(#theta_{K})",-1,1);
RooRealVar ctL ("ctL","cos(#theta_{L})",-1,1);
RooRealVar phi ("phi","#phi",-TMath::Pi(),TMath::Pi());
RooArgSet vars (ctK, ctL, phi);
RooRealVar wei ("wei","weight",1e4);

void createEffHist(int q2Bin, int tagFlag=1, int xbins=25, int ybins = 0, int zbins = 0)
{
  // q2-bin format: [0-8] for one bin, [-1] for each bin recursively
  // tag format:    [1] correctly-tagged. [0] mistagged, [2] each tag, recursively

  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;

  if ( q2Bin > -1 && q2Bin < nBins ) {
    if (tagFlag < 2 && tagFlag > -1) {
      cout<<"Creating efficiency histo for q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
      createEffHistBin(q2Bin, (tagFlag==1), xbins, ybins, zbins);
    }
    if (tagFlag == 2) {
      cout<<"Creating efficiency histo for q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
      createEffHistBin(q2Bin, true,  xbins, ybins, zbins);
      createEffHistBin(q2Bin, false, xbins, ybins, zbins);
    }
  }
  if (q2Bin == -1) {
    cout<<"Creating efficiency histo for all q2 bins"<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if (tagFlag < 2 && tagFlag > -1) {
	cout<<endl<<"q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
	createEffHistBin(q2Bin, (tagFlag==1), xbins, ybins, zbins);
      }
      if (tagFlag == 2) {
	cout<<endl<<"q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
	createEffHistBin(q2Bin, true,  xbins, ybins, zbins);
	createEffHistBin(q2Bin, false, xbins, ybins, zbins);
      }
    }
  }
  
}

void createEffHistBin(int q2Bin, bool tagFlag, int xbins, int ybins, int zbins)
{

  string shortString = Form(tagFlag?"b%ict_testLowStat":"b%iwt_testLowStat",q2Bin);
  string longString  = Form(tagFlag?"q2 bin %i correct-tag":"q2 bin %i wrong-tag",q2Bin);
  int confIndex = (tagFlag?q2Bin:q2Bin+nBins);

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
  t_den->SetBranchAddress( "gen_cos_theta_k", &genCosThetaK );
  t_den->SetBranchAddress( "gen_cos_theta_l", &genCosThetaL );
  t_den->SetBranchAddress( "gen_phi"        , &genPhi       );
  t_den->SetBranchAddress( "mumu_mass"      , &genDimuMass  );
  t_den->SetBranchAddress( "b0_pt"          , &genB0pT      );
  t_den->SetBranchAddress( "b0_eta"         , &genB0eta     );
  t_num->SetBranchAddress( "cos_theta_k" , &recoCosThetaK );
  t_num->SetBranchAddress( "cos_theta_l" , &recoCosThetaL );
  t_num->SetBranchAddress( "phi_kst_mumu", &recoPhi       );
  t_num->SetBranchAddress( "mumuMass"    , &recoDimuMass  );
  t_num->SetBranchAddress( "bPt"         , &recoB0pT      );
  t_num->SetBranchAddress( "bEta"        , &recoB0eta     );
  t_num->SetBranchAddress( "genSignal"   , &genSignal     );
  t_num->SetBranchAddress( "tagB0"       , &tagB0         );

  RooDataSet* data    = new RooDataSet( "data"   , "GEN distribution" , RooArgSet(vars,wei) );
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
    ctK.setVal(genCosThetaK);
    ctL.setVal(genCosThetaL);
    phi.setVal(genPhi);
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
    ctK.setVal(recoCosThetaK);
    ctL.setVal(recoCosThetaL);
    phi.setVal(recoPhi);
    numData->add(vars);    
  }
  cout<<"Dataset prepared"<<endl;

  // compute and print average efficiency
  double avgEff = numData->sumEntries() / wdata->sumEntries();
  cout<<"Average efficiency = "<<avgEff<<endl;

  if (plot) {
    // Plot 1D distributions of numerator and denominator datasets
    double rescFac = 0.03;
    c[confIndex] = new TCanvas(("c_"+shortString).c_str(),("Num and Den 1D projections - "+longString).c_str(),2000,700);
    TLegend* leg = new TLegend(0.4,0.8,0.9,0.9);
    RooPlot* xframe = ctK.frame(Title((longString+" cos(#theta_{K}) distributions").c_str()));
    RooPlot* yframe = ctL.frame(Title((longString+" cos(#theta_{L}) distributions").c_str()));
    RooPlot* zframe = phi.frame(Title((longString+" #phi distributions").c_str()));
    wdata->plotOn(xframe,Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),Name("plDenDist"));
    wdata->plotOn(yframe,Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30));
    wdata->plotOn(zframe,Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30));
    numData->plotOn(xframe,Binning(30),Name("plNumDist"));
    numData->plotOn(yframe,Binning(30));
    numData->plotOn(zframe,Binning(30));
    xframe->GetYaxis()->SetTitleOffset(1.6);
    yframe->GetYaxis()->SetTitleOffset(1.6);
    zframe->GetYaxis()->SetTitleOffset(1.6);
    xframe->SetMaximum(xframe->GetMaximum()*rescFac*1.15);
    yframe->SetMaximum(yframe->GetMaximum()*rescFac*1.15);
    zframe->SetMaximum(zframe->GetMaximum()*rescFac*1.15);
    leg->AddEntry(xframe->findObject("plDenDist"),Form("Generator-level distributions (*%1.2f)",rescFac),"lep");
    leg->AddEntry(xframe->findObject("plNumDist"),"Post-selection RECO distributions","lep");
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

    c[confIndex]->SaveAs( ("dist_GEN_RECO_"+shortString+".pdf").c_str() );
  }

  // create numerator and denominator histograms
  TH3D* denHist = (TH3D*)wdata  ->createHistogram( ("denHist"+shortString).c_str(),
						   ctK,     Binning(xbins,-1,1) ,
						   YVar(ctL,Binning(ybins,-1,1)),
						   ZVar(phi,Binning(zbins,-TMath::Pi(),TMath::Pi())) );
  TH3D* numHist = (TH3D*)numData->createHistogram( ("numHist"+shortString).c_str(),
						   ctK,     Binning(xbins,-1,1) ,
						   YVar(ctL,Binning(ybins,-1,1)),
						   ZVar(phi,Binning(zbins,-TMath::Pi(),TMath::Pi())) );

  // save histograms
  TFile* fout = new TFile( ("effHist_"+shortString+Form("_%i_%i_%i.root",xbins,ybins,zbins)).c_str(), "RECREATE" );
  denHist->Write();
  numHist->Write();
  fout->Write();
  fout->Close();

  if (plot) {
    // plot 1D slices of the efficiency
    vector <TEfficiency*> effHistsX; 
    vector <TEfficiency*> effHistsY;
    vector <TEfficiency*> effHistsZ;
    cx1[confIndex] = new TCanvas(("cx1"+shortString).c_str(),("Binned efficiency - "+longString+" - cos(theta_k) slices").c_str(),1500,1500) ;
    cy1[confIndex] = new TCanvas(("cy1"+shortString).c_str(),("Binned efficiency - "+longString+" - cos(theta_l) slices").c_str(),1500,1500) ;
    cz1[confIndex] = new TCanvas(("cz1"+shortString).c_str(),("Binned efficiency - "+longString+" - phi slices"         ).c_str(),1500,1500) ;
    cx1[confIndex]->Divide(5,5);
    cy1[confIndex]->Divide(5,5);
    cz1[confIndex]->Divide(5,5);

    // width of the slices in the hidden variables ("border" is half of it)
    double border = 0.035;

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

	// plot in canvas
	cx1[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18); 
	effHistsX.back()->Draw();
	cx1[confIndex]->cd(5*j+i+1)->Update(); 
	auto graphx = effHistsX.back()->GetPaintedGraph(); 
	graphx->SetMinimum(0);
	// graphx->SetMaximum(0.2);
	graphx->GetYaxis()->SetTitleOffset(1.7);
	cx1[confIndex]->cd(5*j+i+1)->Update();

	cy1[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18); 
	effHistsY.back()->Draw();
	cy1[confIndex]->cd(5*j+i+1)->Update(); 
	auto graphy = effHistsY.back()->GetPaintedGraph(); 
	graphy->SetMinimum(0);
	// graphy->SetMaximum(0.2); 
	graphy->GetYaxis()->SetTitleOffset(1.7);
	cy1[confIndex]->cd(5*j+i+1)->Update();

	cz1[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18); 
	effHistsZ.back()->Draw();
	cz1[confIndex]->cd(5*j+i+1)->Update(); 
	auto graphz = effHistsZ.back()->GetPaintedGraph(); 
	graphz->SetMinimum(0);
	// graphz->SetMaximum(0.2); 
	graphz->GetYaxis()->SetTitleOffset(1.7);
	cz1[confIndex]->cd(5*j+i+1)->Update();

      }    

    cx1[confIndex]->SaveAs( ("effHist_"+shortString+Form("_%i_%i_%i_CTKslices_dp%i.pdf",xbins,ybins,zbins,(int)(border*200))).c_str() );
    cy1[confIndex]->SaveAs( ("effHist_"+shortString+Form("_%i_%i_%i_CTLslices_dp%i.pdf",xbins,ybins,zbins,(int)(border*200))).c_str() );
    cz1[confIndex]->SaveAs( ("effHist_"+shortString+Form("_%i_%i_%i_PHIslices_dp%i.pdf",xbins,ybins,zbins,(int)(border*200))).c_str() );
  }
    
}
