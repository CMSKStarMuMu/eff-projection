/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ANGULARWT
#define ANGULARWT

#include <math.h>
#include "Math/SpecFunc.h"
#include "TMath.h"

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooFit.h"
#include "Riostream.h"
 
class AngularWT : public RooAbsPdf {
 protected:

  RooRealProxy ctK ;
  RooRealProxy ctL ;
  RooRealProxy phi ;
  RooRealProxy Fl ;
  RooRealProxy P1 ;
  RooRealProxy P2 ;
  RooRealProxy P3 ;
  RooRealProxy P4p ;
  RooRealProxy P5p ;
  RooRealProxy P6p ;
  RooRealProxy P8p ;
  
  Double_t evaluate() const ;

 public:
  AngularWT() {} ; 
  AngularWT(const char *name, const char *title,
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
	    RooAbsReal& _P8p);
  AngularWT(const AngularWT& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new AngularWT(*this,newname); }
  inline virtual ~AngularWT() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  ClassDef(AngularWT,1) // PDF for angular decay rate description
    };
 
#endif
