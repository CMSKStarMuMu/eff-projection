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

using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* c [nBins];

void fit_genMCBin(int q2Bin, bool plot, bool save)
{

  string shortString = Form("b%i",q2Bin);
  string longString  = Form("q2 bin %i",q2Bin);
  int confIndex = q2Bin;

  // Load datasets
  TFile* fin = TFile::Open( ("effDataset_"+shortString+"wt.root").c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: effDataset_"+shortString+"wt.root"<<endl;
    return;
  }
  RooWorkspace* wsp = (RooWorkspace*)fin->Get(("ws_"+shortString+"wt").c_str());
  if ( !wsp || wsp->IsZombie() ) {
    cout<<"Workspace not found in file: effDataset_"+shortString+"wt.root"<<endl;
    return;
  }
  RooDataSet* data = (RooDataSet*)wsp->data(("data_"+shortString+"wt").c_str());
  RooRealVar* ctK = wsp->var("ctK");
  RooRealVar* ctL = wsp->var("ctL");
  RooRealVar* phi = wsp->var("phi");
  RooArgSet vars (* ctK,* ctL,* phi);

  // define angular parameters
  RooRealVar* Fl = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1 = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2 = new RooRealVar("P2","P_{2}",0,-1,1);
  RooRealVar* P3 = new RooRealVar("P3","P_{3}",0,-1,1);
  RooRealVar* P4p = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p = new RooRealVar("P6p","P'_{6}",0,-1,1);
  RooRealVar* P8p = new RooRealVar("P8p","P'_{8}",0,-1,1);

  RooAbsPdf* _AnglesPDF = new AngularRT("_AnglesPDF","_AnglesPDF",*ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  RooFitResult * fitResult = _AnglesPDF->fitTo(*data,Minimizer("Minuit2","migrad"),Save(true),Timer(true),NumCPU(6)); 
  fitResult->Print("v");

  if (save) {
    TFile f (("fitResult_genMC_"+shortString+".root").c_str(),"UPDATE") ;
    f.cd();
    fitResult->Write("fitResult");
    f.Close();
  }

  if (!plot) return;

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

  bool plot = true;
  if ( argc >= 3 && atoi(argv[2]) == 0 ) plot = false;

  bool save = true;
  if ( argc >= 4 && atoi(argv[3]) == 0 ) save = false;

  if ( q2Bin > -1 && q2Bin < nBins ) {
    cout<<"Fitting distributions for q2 bin "<<q2Bin<<endl;
    fit_genMCBin(q2Bin, plot, save);
  }
  if (q2Bin == -1) {
    cout<<"Fitting distributions for all q2 bins"<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      cout<<endl<<"q2 bin "<<q2Bin<<endl;
      fit_genMCBin(q2Bin, plot, save);
    }
  }

  return 0;

}

