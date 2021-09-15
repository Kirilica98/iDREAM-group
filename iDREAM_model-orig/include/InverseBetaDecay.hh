#ifndef INVERSE_BETA_DECAY_HH
#define INVERSE_BETA_DECAY_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>
//#include "InverseBetaDecaysTable.hh"

#include "G4PhysicsOrderedFreeVector.hh"
#include "G4MaterialPropertyVector.hh"
 
class	InverseBetaDecay 
{

  public:

	InverseBetaDecay();
	virtual	~InverseBetaDecay();

	//G4String	GetNueSourceName()		const {return fNueSourceName;}
	//G4String	GetNueSourceType()		const {return fNueSourceType;}	

	//static	const	InverseBetaDecaysTable*		GetInverseBetaDecaysTable();
	static		InverseBetaDecay*		GetInverseBetaDecay(G4String name, G4String type);

	G4double	DefinePositronEmissionEnergy();
	G4double	DefineNeutronEmissionEnergy();
	G4double	GetNeutrinoEnergy();
	
	G4ThreeVector	DefinePositronEmissionDirection(G4ThreeVector Direction);
	G4ThreeVector	DefineNeutronEmissionDirection(G4ThreeVector Direction);
	
	G4ThreeVector	DefinePositronPolarization(G4ThreeVector Direction);
	G4ThreeVector	DefineNeutronPolarization(G4ThreeVector Direction);

   private:

	//G4String			fNueSourceName;
	//G4String			fNueSourceType;
		
	G4PhysicsOrderedFreeVector*	fEmissionIntegral;
	G4double			fEIImax;
	
	G4double			fNueEnergy;
	G4double			fPositronEnergy;	
	G4double			fDelta;
	G4double			fPositronScatteringAngle;
	G4double			fPositronScatteringAzimuth;
	G4double			fNeutronScatteringAngle;
	
	void				SetDecayProperties(void);
	void				BuildEmissionIntegral(void);
	
	//static	InverseBetaDecaysTable	fInverseBetaDecaysTable;
	
	G4int Quadratic(G4double *x, G4double a, G4double b, G4double c);
	G4int Cubic(G4double *x, G4double a, G4double b, G4double c);
};

#endif








