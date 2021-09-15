#include <fstream>
#include "InverseBetaDecay.hh"
#include "Randomize.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//#include "G4VersionNumber.hh"

//InverseBetaDecaysTable	InverseBetaDecay::fInverseBetaDecaysTable;

InverseBetaDecay::InverseBetaDecay(){
	SetDecayProperties();
	BuildEmissionIntegral();
	
	//fInverseBetaDecaysTable.push_back(this);
}

InverseBetaDecay::~InverseBetaDecay()	{
	//InverseBetaDecaysTable::iterator iter = fInverseBetaDecaysTable.begin();
	//while((iter != fInverseBetaDecaysTable.end()) && (*iter != this)) iter++;
	//if(iter != fInverseBetaDecaysTable.end()) fInverseBetaDecaysTable.erase(iter);

	delete fEmissionIntegral;
}
/*
const	InverseBetaDecaysTable*	InverseBetaDecay::GetInverseBetaDecaysTable()	{
	return &fInverseBetaDecaysTable;
}

InverseBetaDecay*	InverseBetaDecay::GetInverseBetaDecay(G4String name, G4String type)	{
	size_t i=0;
	for(i=0; i<fInverseBetaDecaysTable.size(); i++)	{
		if((fInverseBetaDecaysTable[i]->GetNueSourceName() == name) &&
			(fInverseBetaDecaysTable[i]->GetNueSourceType() == type))
			return fInverseBetaDecaysTable[i];
	}

	InverseBetaDecay* theInverseBetaDecaySource = new InverseBetaDecay(name,type);
//	fInverseBetaDecaysTable.push_back(theInverseBetaDecaySource);
	return fInverseBetaDecaysTable[i];

}
*/
G4double	InverseBetaDecay::GetNeutrinoEnergy()	{
	return fNueEnergy;
}

G4double	InverseBetaDecay::DefinePositronEmissionEnergy()	{
G4double theEIIvalue;
G4double theEnergy = 0.0;
G4double theTotalPositronEnergy = 0.0;

G4double	acoeff,
			bcoeff,
			ccoeff;
G4double	x[2];

	theEIIvalue = G4UniformRand()*fEIImax;
	theEnergy = fEmissionIntegral->GetEnergy(theEIIvalue);
	theTotalPositronEnergy = theEnergy + electron_mass_c2;

	acoeff = 0.017 * sqrt(theTotalPositronEnergy*theTotalPositronEnergy-
				electron_mass_c2*electron_mass_c2)/theTotalPositronEnergy;
//	acoeff -= 1.2 * (theTotalPositronEnergy + fDelta)/proton_mass_c2;

	bcoeff = -1.0;
	
	ccoeff = bcoeff - acoeff + 2.0 * G4UniformRand();

	Quadratic(x,acoeff,bcoeff,ccoeff);
				
	fPositronScatteringAngle = x[1];

	fNueEnergy = (proton_mass_c2*(theTotalPositronEnergy+fDelta)+0.5*(fDelta*fDelta -
			electron_mass_c2*electron_mass_c2))/(proton_mass_c2-theTotalPositronEnergy+
			sqrt(theTotalPositronEnergy*theTotalPositronEnergy-electron_mass_c2*electron_mass_c2)*fPositronScatteringAngle);

	fNueEnergy *= MeV;
	fPositronEnergy	= theEnergy;	
	return theEnergy;
}

G4double	InverseBetaDecay::DefineNeutronEmissionEnergy()	{

G4double theEnergy = 0.0;

	theEnergy = fNueEnergy - fPositronEnergy - electron_mass_c2 + proton_mass_c2;
	theEnergy -= neutron_mass_c2;
				
	return theEnergy;
}

G4ThreeVector	InverseBetaDecay::DefinePositronEmissionDirection(G4ThreeVector Direction)
{
/* For the moment positron emission angular distribution is isotropic	*/

  G4ThreeVector theDirection;

  G4double	CosTheta,
		SinTheta,
		Phi,
		SinPhi,
		CosPhi;
  
  CosTheta = fPositronScatteringAngle;
    
  SinTheta = sqrt(1.0 - (CosTheta * CosTheta));

  Phi = 2*M_PI*G4UniformRand();
  fPositronScatteringAzimuth = Phi;
  SinPhi = sin(Phi);
  CosPhi = cos(Phi);

  theDirection.setX(SinTheta * CosPhi);
  theDirection.setY(SinTheta * SinPhi);
  theDirection.setZ(CosTheta);

  theDirection.rotateUz(Direction);

  return theDirection;
}

G4ThreeVector	InverseBetaDecay::DefineNeutronEmissionDirection(G4ThreeVector Direction)
{
/* For the moment neutron emission angular distribution is flat inside 0-max

	max angle = sqrt(2.0 * fNueEnergy * fDelta - (fDelta*fDelta-electron_mass_c2*electron_mass_c2))/fNueEnergy;
	*/
	
  G4ThreeVector 	theDirection;

  G4double	ThetaMax,
			CosTheta,
			SinTheta,
			Phi,
			SinPhi,
			CosPhi;

		
  ThetaMax = 	sqrt(2.0 * fNueEnergy * fDelta - (fDelta*fDelta-electron_mass_c2*electron_mass_c2))/fNueEnergy;
  	
  CosTheta = 1.0 - (1.0 - ThetaMax) * G4UniformRand();
  SinTheta = sqrt(1.0 - (CosTheta * CosTheta));

  //Phi = 2*M_PI*G4UniformRand();
  Phi = fPositronScatteringAzimuth + M_PI;
  SinPhi = sin(Phi);
  CosPhi = cos(Phi);

  theDirection.setX(SinTheta * CosPhi);
  theDirection.setY(SinTheta * SinPhi);
  theDirection.setZ(CosTheta);
  				
  theDirection.rotateUz(Direction);
//theEnergy = fNueEnergy - fPositronEnergy - electron_mass_c2 + proton_mass_c2 - neutron_mass_c2;  
//  G4cout << fNueEnergy << "	" << acos(CosTheta)*180.0/M_PI << G4endl; 
//  G4cout << theEnergy << "	" << acos(CosTheta)*180.0/M_PI << G4endl;    
  return theDirection;
}

G4ThreeVector	InverseBetaDecay::DefinePositronPolarization(G4ThreeVector Direction)
{
G4ThreeVector 	thePolarization;

  thePolarization = (-1)*Direction;

  thePolarization = thePolarization.unit();

  return thePolarization;
}

G4ThreeVector	InverseBetaDecay::DefineNeutronPolarization(G4ThreeVector Direction)
{
G4ThreeVector 	thePolarization;

  thePolarization = (0.0)*Direction;

//  thePolarization = thePolarization.unit();

  return thePolarization;
}



void	InverseBetaDecay::SetDecayProperties()
{
	fDelta = neutron_mass_c2 - proton_mass_c2;
}

void	InverseBetaDecay::BuildEmissionIntegral()	{
std::ifstream	ifs;
std::vector<G4double>	thePositronSpectra;
std::vector<G4double>	thePositronMomenta;
G4double	theNueEnergy,
		theNueValue,
		thePositronEnergy,
		thePositronValue;

G4double	theRecCorrection = 1.261*0.261/(1 + 3.0 * 1.261 * 1.261);
G4double	theWMCorrection  = 1.261*3.706/(1 + 3.0 * 1.261 * 1.261);
G4double	theCorrection;
G4double	acoeff = 1.00,
		bcoeff = proton_mass_c2 + fDelta,
		ccoeff;
G4double	x[2];
		
//	if(!getenv("NueSources"))
//		G4Exception("NueSources.cc","env NueSources is not set",RunMustBeAborted,"env NueSources should be set");
	
    G4String theEmissionSpectraFile = getenv("PWD");
	theEmissionSpectraFile += "/";
	//theEmissionSpectraFile += fNueSourceName;
	//theEmissionSpectraFile += "/";
	//theEmissionSpectraFile += fNueSourceType;
	theEmissionSpectraFile += "EmissionSpectra.dat";

  	ifs.open(theEmissionSpectraFile.c_str());
	if(!ifs.is_open()) {
		G4cout << "Couldn't open file: " << theEmissionSpectraFile << G4endl;
		return;
	}
	while (ifs >> theNueEnergy >> theNueValue) {
		ccoeff = proton_mass_c2*(fDelta - theNueEnergy) + (fDelta*fDelta - electron_mass_c2*electron_mass_c2)/2.0;
		Quadratic(x,acoeff,bcoeff,ccoeff);
		thePositronEnergy = x[1];
		if(thePositronEnergy - electron_mass_c2 >= 0.0)	{
			thePositronValue = theNueValue*
				sqrt(thePositronEnergy*thePositronEnergy-electron_mass_c2*electron_mass_c2)*
				thePositronEnergy;
			theCorrection = (thePositronEnergy + fDelta + (thePositronEnergy*thePositronEnergy-electron_mass_c2*electron_mass_c2)
					/thePositronEnergy)/proton_mass_c2;

			thePositronValue *= 1.0 + 2.0 * (theRecCorrection - theWMCorrection) * theCorrection;
			
			thePositronEnergy -= electron_mass_c2;	

			thePositronSpectra.push_back(thePositronValue);
			thePositronEnergy *= MeV;
			thePositronMomenta.push_back(thePositronEnergy);
		}		
	}
  	ifs.close();
	
  	G4int	theNum = thePositronSpectra.size();
	
	fEmissionIntegral = new G4PhysicsOrderedFreeVector();

	G4MaterialPropertyVector* theEmissionBetaVector =
			new G4MaterialPropertyVector(&thePositronMomenta[0],
						&thePositronSpectra[0], theNum);

	if (theEmissionBetaVector) {

		// Retrieve the first intensity point in vector
		// of (particle momentum, intensity) pairs
		

		G4double currentIN = (*theEmissionBetaVector)[0];

		
		if (currentIN >= 0.0) {

			// Create first (photon momentum, Emission
                        // Integral pair
                        

			G4double currentPM = theEmissionBetaVector->Energy(0);					 
					 
			G4double currentEII = 0.0;

			fEmissionIntegral->InsertValues(currentPM , currentEII);

			// Set previous values to current ones prior to loop

			G4double prevPM  = currentPM;
			G4double prevEII = currentEII;
                	G4double prevIN  = currentIN;

			// loop over all (photon momentum, intensity)
			// pairs stored for this material


            for (size_t i = 1;i < theEmissionBetaVector->GetVectorLength();i++)	{
                currentPM = theEmissionBetaVector->Energy(i);
                currentIN = (*theEmissionBetaVector)[i];				
				
				currentEII = 0.5 * (prevIN + currentIN);
				currentEII = prevEII + (currentPM - prevPM) * currentEII;

				fEmissionIntegral->InsertValues(currentPM, currentEII);

				prevPM  = currentPM;
				prevEII = currentEII;
				prevIN  = currentIN;

			}

		}
		delete theEmissionBetaVector;

		fEIImax = fEmissionIntegral->GetMaxValue();
	}

}

/* Quadratic equation solution. Real coefficients case.

   int Quadratic(double *x,double a,double b,double c);
   Parameters:
   x - solution array (size 2). On output:
       2 real roots -> then x is filled with them;
       2 complex-conjugate roots -> x[0] is real part, 
         x[1] is non-negative imaginary part. 
       other cases -> x[0] unique valid root if (-1) returned,
                     no valid roots otherwise.
   a, b, c - coefficients: ax^2 + bx + c = 0. 
   Returns: 2 - 2 real and distince roots;
            1 - 1 real distinct root (x[0]=x[1]);
            0 - 2 complex roots;
            -1 - one real root in case a==0;
            -2 - no roots in case a=0, b=0;
            -3 - infinite number of roots (a=b=c=0). 
*/

#include <math.h>   /* for sqrt() */

G4int InverseBetaDecay::Quadratic(G4double *x, G4double a, G4double b, G4double c) {
 G4double d;
 /* degenerated cases */
 if(a==0.) {
  if(b==0.) {
   if(c==0.) return(-3);
   return(-2);
  }
  x[0]=-c/b; return(-1);
 }
 /* the main case */
 d=b*b-4.*a*c;       /* the discriminant */
 /* one distinct root */
 if(d==0.) {
  x[0]=x[1]=-b/(2.*a); return(1);
 }
 /* conjugate complex roots */
 if(d<0.) {
  G4double t=0.5/a;
  x[0]=-b*t; x[1]=sqrt(-d)*t;
  return(0);
 }
 /* 2 real roots: avoid subtraction of 2 close numbers */
 if(b>=0.) d=(-0.5)*(b+sqrt(d));
 else d=(-0.5)*(b-sqrt(d));
 x[0]=d/a; x[1]=c/d;
 return(2);
}

G4int InverseBetaDecay::Cubic(double *x,double a,double b,double c) {
  double q,r,r2,q3;
  q=(a*a-3.*b)/9.; r=(a*(2.*a*a-9.*b)+27.*c)/54.;
  r2=r*r; q3=q*q*q;
  if(r2<q3) {
    double t=acos(r/sqrt(q3));
    a/=3.; q=-2.*sqrt(q);
    x[0]=q*cos(t/3.)-a;
    x[1]=q*cos((t+twopi)/3.)-a;
    x[2]=q*cos((t-twopi)/3.)-a;
    return(3);
  }
  else {
    double aa,bb;
    if(r<=0.) r=-r;
    aa=-pow(r+sqrt(r2-q3),1./3.); 
    if(aa!=0.) bb=q/aa;
    else bb=0.;
    a/=3.; q=aa+bb; r=aa-bb; 
    x[0]=q-a;
    x[1]=(-0.5)*q-a;
    x[2]=(sqrt(3.)*0.5)*fabs(r);
    if(x[2]==0.) return(2);
    return(1);
  }
}

