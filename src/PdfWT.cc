/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "PdfWT.h" 

ClassImp(PdfWT) 

PdfWT::PdfWT(const char *name, const char *title, 
	     RooAbsReal& _ctK,
	     RooAbsReal& _ctL,
	     RooAbsReal& _phi,
	     RooAbsReal& _Fl,
	     RooAbsReal& _P1,
	     RooAbsReal& _P2,
	     RooAbsReal& _P3,
	     RooAbsReal& _P4p,
	     RooAbsReal& _P5p,
	     RooAbsReal& _P6p,
	     RooAbsReal& _P8p,
	     RooArgList& _EffCoeff,
	     int _degree) :
  RooAbsPdf(name,title), 
  ctK("ctK","ctK",this,_ctK),
  ctL("ctL","ctL",this,_ctL),
  phi("phi","phi",this,_phi),
  Fl("Fl","Fl",this,_Fl),
  P1("P1","P1",this,_P1),
  P2("P2","P2",this,_P2),
  P3("P3","P3",this,_P3),
  P4p("P4p","P4p",this,_P4p),
  P5p("P5p","P5p",this,_P5p),
  P6p("P6p","P6p",this,_P6p),
  P8p("P8p","P8p",this,_P8p),
  EffCoeff(_EffCoeff),
  degree(_degree)
{
}


PdfWT::PdfWT(const PdfWT& other, const char* name) :  
  RooAbsPdf(other,name), 
  ctK("ctK",this,other.ctK),
  ctL("ctL",this,other.ctL),
  phi("phi",this,other.phi),
  Fl("Fl",this,other.Fl),
  P1("P1",this,other.P1),
  P2("P2",this,other.P2),
  P3("P3",this,other.P3),
  P4p("P4p",this,other.P4p),
  P5p("P5p",this,other.P5p),
  P6p("P6p",this,other.P6p),
  P8p("P8p",this,other.P8p),
  EffCoeff(other.EffCoeff),
  degree(other.degree)
{
}



Double_t PdfWT::evaluate() const 
{ 
  // computing values of associated-Legendre polynomials at current coordinates
  // up to the maximum order of the efficiency 
  // (but not less than 2, needed for the decay rate)
  std::vector< double > LegCTL;
  std::vector< double > LegCTK;

  for (int l=0; l<=TMath::Max(2,degree); ++l) for (int m=0; m<=l; ++m) {
      LegCTL.push_back(ROOT::Math::assoc_legendre(l,m,ctL));
      LegCTK.push_back(ROOT::Math::assoc_legendre(l,m,ctK));
    }

  // computing the efficiency value at current coordinates
  double eff = 0;
  int i_ord = 0;
  int ctK_ord, ctL_ord, phi_ord;
  for (ctK_ord=0; ctK_ord<=degree; ++ctK_ord)
    for (ctL_ord=0; ctL_ord<=degree; ++ctL_ord)
      for (phi_ord=-1*TMath::Min(ctK_ord,ctL_ord); phi_ord<=TMath::Min(ctK_ord,ctL_ord); ++phi_ord) {

	if (phi_ord<0)      eff += ((RooAbsReal*)EffCoeff.at(i_ord))->getVal() * 
			      LegCTL[ctL_ord*(ctL_ord+1)/2.0-phi_ord] * 
			      LegCTK[ctK_ord*(ctK_ord+1)/2.0-phi_ord] * 
			      sin(-1.0*phi_ord*phi);

	else if (phi_ord>0) eff += ((RooAbsReal*)EffCoeff.at(i_ord))->getVal() * 
			      LegCTL[ctL_ord*(ctL_ord+1)/2.0+phi_ord] * 
			      LegCTK[ctK_ord*(ctK_ord+1)/2.0+phi_ord] * 
			      cos(phi_ord*phi);

	else                eff += ((RooAbsReal*)EffCoeff.at(i_ord))->getVal() * 
			      LegCTL[ctL_ord*(ctL_ord+1)/2.0]         * 
			      LegCTK[ctK_ord*(ctK_ord+1)/2.0];

	++i_ord;

      }

  // returning PDF value, as efficiency times the angular decay rate
  return (9./(32 * 3.14159265) * eff * (
					(4.0/9.0)             * LegCTL[0] * LegCTK[0] +
					(4.0*Fl/3.0-4.0/9.0)  * LegCTL[0] * LegCTK[3] +
					(2.0/9.0-2.0*Fl/3.0)  * LegCTL[3] * LegCTK[0] +
					(-2.0/9.0-2.0*Fl/3.0) * LegCTL[3] * LegCTK[3] +
					
					(-4.0/3.0)*P2*(1-Fl) * LegCTL[1] * LegCTK[0] +
					( 4.0/3.0)*P2*(1-Fl) * LegCTL[1] * LegCTK[3] +
					
					(1-Fl)/18.0*P1 * LegCTL[5] * LegCTK[5] * cos(2*phi) +
					(1-Fl)/ 9.0*P3 * LegCTL[5] * LegCTK[5] * sin(2*phi) +
					
					( 2.0/9.0)*sqrt(Fl-Fl*Fl)*P4p * LegCTL[4] * LegCTK[4] * cos(phi) +
					(-2.0/3.0)*sqrt(Fl-Fl*Fl)*P5p * LegCTL[2] * LegCTK[4] * cos(phi) +
					(-2.0/9.0)*sqrt(Fl-Fl*Fl)*P8p * LegCTL[4] * LegCTK[4] * sin(phi) +
					(-2.0/3.0)*sqrt(Fl-Fl*Fl)*P6p * LegCTL[2] * LegCTK[4] * sin(phi)
					
					)
  	  );

}



namespace {
  Bool_t fullRangeCosT(const RooRealProxy& x ,const char* range)
  {
    // set accepted integration range for cosTheta variables
    return range == 0 || strlen(range) == 0
      ? std::fabs(x.min() + 1.) < 1.e-5 && std::fabs(x.max() - 1.) < 1.e-5
      : std::fabs(x.min(range) + 1.) < 1.e-5 && std::fabs(x.max(range) - 1.) < 1.e-5;
  }
  Bool_t fullRangePhi(const RooRealProxy& x ,const char* range)
  {
    // set accepted integration range for phi variable
    return range == 0 || strlen(range) == 0
      ? std::fabs(x.min() + TMath::Pi()) < 1.e-3 && std::fabs(x.max() - TMath::Pi()) < 1.e-3
      : std::fabs(x.min(range) + TMath::Pi()) < 1.e-3 && std::fabs(x.max(range) - TMath::Pi()) < 1.e-3;
  }
}

Int_t PdfWT::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  // use analytical integral for the integration over all the three variables
  // note that only the case with efficiency order greater than one are accepted:
  // this can be easily extendend by adding a specific formula for other cases,
  // but so far there is no intention to use such efficiencies
  if (matchArgs(allVars,analVars,ctK,ctL,phi))
    if ( fullRangeCosT(ctK,rangeName) && fullRangeCosT(ctL,rangeName) && fullRangePhi(phi,rangeName) && degree>1 )
      return 1 ;

  // the lack of analytical integral for the subsets of angular variables does not slow down the fit
  // since only the complete integration is used there
  // if one wants to speed up also the PDF projection for plotting, the other analytical integrals can be computed
  // but it seems a huge effort to me
  return 0 ;

}

Double_t PdfWT::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code>0 && code<2) ;

  // use the analytical formula to compute the integral of the square of each associated Legendre polynomial
  std::vector< double > LegSqInt;

  for (int l=0; l<=2; ++l)
    for (int m=0; m<=l; ++m)
      LegSqInt.push_back(2.0/(2.0*l+1.0)*TMath::Factorial(l+m)/TMath::Factorial(l-m));

  // return the formula of the integrated PDF
  Double_t ret = 0;

  if (degree>1) 
    ret =  9./(32*3.14159265) * (
				 (4.0/9.0)             * ((RooAbsReal*)EffCoeff.at(0))->getVal()          * LegSqInt[0] * LegSqInt[0] * (2*TMath::Pi()) +
				 (4.0*Fl/3.0-4.0/9.0)  * ((RooAbsReal*)EffCoeff.at(4*degree+2))->getVal() * LegSqInt[0] * LegSqInt[3] * (2*TMath::Pi()) +
				 (2.0/9.0-2.0*Fl/3.0)  * ((RooAbsReal*)EffCoeff.at(2))->getVal()          * LegSqInt[3] * LegSqInt[0] * (2*TMath::Pi()) +
				 (-2.0/9.0-2.0*Fl/3.0) * ((RooAbsReal*)EffCoeff.at(4*degree+8))->getVal() * LegSqInt[3] * LegSqInt[3] * (2*TMath::Pi()) +
				 
				 (-4.0/3.0)*P2*(1-Fl) * ((RooAbsReal*)EffCoeff.at(1))->getVal()          * LegSqInt[1] * LegSqInt[0] * (2*TMath::Pi()) +
				 ( 4.0/3.0)*P2*(1-Fl) * ((RooAbsReal*)EffCoeff.at(4*degree+4))->getVal() * LegSqInt[1] * LegSqInt[3] * (2*TMath::Pi()) +
				 
				 (1-Fl)/18.0*P1 * ((RooAbsReal*)EffCoeff.at(4*degree+10))->getVal() * LegSqInt[5] * LegSqInt[5] * (TMath::Pi()) +
				 (1-Fl)/ 9.0*P3 * ((RooAbsReal*)EffCoeff.at(4*degree+ 6))->getVal() * LegSqInt[5] * LegSqInt[5] * (TMath::Pi()) +
				 
				 ( 2.0/9.0)*sqrt(Fl-Fl*Fl)*P4p * ((RooAbsReal*)EffCoeff.at(4*degree+9))->getVal() * LegSqInt[4] * LegSqInt[4] * (TMath::Pi()) +
				 (-2.0/3.0)*sqrt(Fl-Fl*Fl)*P5p * ((RooAbsReal*)EffCoeff.at(4*degree+5))->getVal() * LegSqInt[2] * LegSqInt[4] * (TMath::Pi()) +
				 (-2.0/9.0)*sqrt(Fl-Fl*Fl)*P8p * ((RooAbsReal*)EffCoeff.at(4*degree+7))->getVal() * LegSqInt[4] * LegSqInt[4] * (TMath::Pi()) +
				 (-2.0/3.0)*sqrt(Fl-Fl*Fl)*P6p * ((RooAbsReal*)EffCoeff.at(4*degree+3))->getVal() * LegSqInt[2] * LegSqInt[4] * (TMath::Pi())
				 
				 );

  return ret ;
 
}
