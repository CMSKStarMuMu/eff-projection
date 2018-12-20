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

#include "AngularRT.h"
#include "AngularWT.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* c [nBins];

void fit_genMCBin(int q2Bin)
{

  string shortString = Form("b%i_JpsiCh_LegTest",q2Bin);
  string longString  = Form("Jpsi q2 bin %i",q2Bin);
  int confIndex = q2Bin;

  RooRealVar* ctK = new RooRealVar("ctK","cos(#theta_{K})",-1,1);
  RooRealVar* ctL = new RooRealVar("ctL","cos(#theta_{L})",-1,1);
  RooRealVar* phi = new RooRealVar("phi","#phi",-TMath::Pi(),TMath::Pi());
  RooArgSet vars (* ctK,* ctL,* phi);

  RooRealVar* Fl = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1 = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2 = new RooRealVar("P2","P_{2}",0,-1,1);
  RooRealVar* P3 = new RooRealVar("P3","P_{3}",0,-1,1);
  RooRealVar* P4p = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p = new RooRealVar("P6p","P'_{6}",0,-1,1);
  RooRealVar* P8p = new RooRealVar("P8p","P'_{8}",0,-1,1);

  RooAbsPdf* _AnglesPDF = new AngularRT("_AnglesPDF","_AnglesPDF",*ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  // _AnglesPDF->Print("v");
  // P2->setVal(1);
  // ctL->setVal(1);
  // cout<<_AnglesPDF->getValV()<<endl;
  // return;

  // Load ntuples
  TChain* t_den = new TChain();
  t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN/2016MC_GEN_B0ToJpsiKStar_BFilter.root/ntuple");
  int denEntries = t_den->GetEntries() /1000; // DEBUG

  double genCosThetaK, genCosThetaL, genPhi, genDimuMass, genB0pT, genB0eta;
  t_den->SetBranchAddress( "cos_theta_k" , &genCosThetaK );
  t_den->SetBranchAddress( "cos_theta_l" , &genCosThetaL );
  t_den->SetBranchAddress( "phi_kst_mumu", &genPhi       );
  t_den->SetBranchAddress( "genq2"       , &genDimuMass  );
  t_den->SetBranchAddress( "genbPt"      , &genB0pT      );
  t_den->SetBranchAddress( "genbEta"     , &genB0eta     );

  RooDataSet* data = new RooDataSet( "data", "GEN distribution before GEN-filter" , vars );

  int counter=0;
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
    data->add( vars );
  }

  RooFitResult * fitResult = _AnglesPDF->fitTo(*data,Minimizer("Minuit2","migrad"),Save(true),Timer(true)); 
  // RooFitResult * fitResult = _AnglesPDF->fitTo(*data,Save(true),Timer(true),NumCPU(6));
  // RooFitResult * fitResult = _AnglesPDF->fitTo(*data,Extended(true),Save(true),Timer(true));

  fitResult->Print("v");
  // return; 			// DEBUG
  // TFile f (("fitResult_genMC_"+shortString+".root").c_str(),"UPDATE") ;
  // f.cd();
  // fitResult->Write("fitResult");
  // f.Close();

  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to GEN-level MC - "+longString).c_str(),2000,700);
  TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
  RooPlot* xframe = ctK->frame(Title((longString+" - cos(#theta_{K}) distribution").c_str()));
  RooPlot* yframe = ctL->frame(Title((longString+" - cos(#theta_{L}) distribution").c_str()));
  RooPlot* zframe = phi->frame(Title((longString+" - #phi distribution").c_str()));
  data->plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40),Name("plData"));
  data->plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40));
  data->plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40));
  _AnglesPDF->plotOn(xframe,LineWidth(1),Name("plPDF"));
  _AnglesPDF->plotOn(yframe,LineWidth(1));
  _AnglesPDF->plotOn(zframe,LineWidth(1));
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
  leg->AddEntry(xframe->findObject("plData"),"Generation-level distribution" ,"lep");
  leg->AddEntry(xframe->findObject("plPDF" ),"Angular decay rate" ,"l");

  c[confIndex]->Divide(3,1);
  c[confIndex]->cd(1);
  gPad->SetLeftMargin(0.15);
  xframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(2);
  gPad->SetLeftMargin(0.15);
  yframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(3);
  gPad->SetLeftMargin(0.15);
  zframe->Draw();
  leg->Draw("same");

  c[confIndex]->SaveAs( ("fitResult_genMC_"+shortString+".pdf").c_str() );

}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin, [-1] for each bin recursively

  if ( argc < 2) return 1;
  int q2Bin    = atoi(argv[1]);

  if ( q2Bin > -1 && q2Bin < nBins ) {
    cout<<"Fitting distributions for q2 bin "<<q2Bin<<endl;
    fit_genMCBin(q2Bin);
  }
  if (q2Bin == -1) {
    cout<<"Fitting distributions for all q2 bins"<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      cout<<endl<<"q2 bin "<<q2Bin<<endl;
      fit_genMCBin(q2Bin);
    }
  }

  return 0;

}

