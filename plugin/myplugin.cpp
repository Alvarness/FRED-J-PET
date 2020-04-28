/*
 *  myplugin.cpp : Fred plugin
 *
 *  asch 2015 - 2017
 *
 */

#include "FredXfunc.h"
#include "Scorer.h"

#include <cstdlib>
#include <fstream>
using namespace std;


//========================================================================================
//=================      PLUGIN VARIABLES and DATA STRUCTURES     ========================
//========================================================================================

ofstream ftracks_phot;

//========================================================================================
//========================================================================================
//========================================================================================


extern "C" int UserHook_init(const char *vers){
	cout<<magentacolor;
	cout<<"Plugin user init:"<<endl<<endl;
// 	cout<<"\tFred "<<vers<<endl; // display Fred version

	ftracks_phot.open("out/photons.txt");

// 	cout<<"myFlag = "<<getBoolParam("myFlag",false)<<endl;
	// cout<<"myInt = "<<getIntParam("myInt",0)<<endl;
	// cout<<"myFloat = "<<getFloatParam("myFloat",0)<<endl;
	cout<<endl<<normalcolor;
	return 0; // OK
}

extern "C" int UserHook_close(){
	cout<<magentacolor;
	cout<<endl<<"Plugin user close:"<<endl;
	cout<<normalcolor;
	return 0; // return value is ignored
}

int lastTracked;

extern "C" void  UserHook_step_aft(Step *stp){
	vec3dRT xi;
	getPosition_A(stp,xi);
	double T = getKineticEnergy_A(stp);
	if(xi.z<0) return;
	int type = getType(stp);

	if(type!=GHOST_ID) lastTracked = type;
	if(type==PHOTON_ID)   ftracks_phot<<xi<<' '<<T<<' '<<endl;

	// cout<<getUID(stp)<<' '<<getParentUID(stp)<<' '<<getAncestorUID(stp)<<' '<<type<<' '<<getParticle_generation(stp)<<endl;
	if(type==PHOTON_ID) {
		setStepDeltaTime(stp,getStepLength(stp)/3e10);
		// cout<<"time now "<<getUID(stp)<<' '<<getTime_B(stp)<<endl;
	}
}

extern "C" void  UserHook_step_final(Step *stp){
	vec3dRT xi;
	getPosition_A(stp,xi);
	double T = getKineticEnergy_A(stp);

	if(lastTracked==PHOTON_ID) {
		ftracks_phot<<xi<<' '<<T<<' '<<endl;
		ftracks_phot<<endl;
	}
}

float massAttenuation_photon(Step *stp){
	// compute  the total mass attenuation coeff for the photon
	double mu_rho_tot_material = 0.;

	int I_mat = getImat_A(stp); //get material index at initial point A
	double T = getKineticEnergy_A(stp);

	string currReg = getRegion_A(stp);

	if(currReg.rfind("scint")==0){
		cout<<currReg<<endl;
		mu_rho_tot_material = 0.09443;
	}

	return mu_rho_tot_material;
}



void interactionEvent_Photon(Step *stp,const vec3dRT &v0,const vec3dRT &n1,const vec3dRT &n2){
	vec3dRT vnew;
	getVelocityVersor_A(stp,vnew);
	vec3dRT pol;
	getPolarizationDirection_A(stp,pol);
	// cout<<"v0 "<<vnew<<endl;
	//pushParticle(stp,ELECTRON_ID,getKineticEnergy_A(stp)/4,v0);
	cout<<"ciao "<<getUID(stp)<<' '<<getKineticEnergy_A(stp)<<' '<<pol<<' '<<getTime_A(stp)<<endl;
	extinguishRay(stp);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern "C" float  UserHook_Mass_Attenuation(Step *stp){
	// NB: returns the mass attenuation coefficient of interaction in the material
	switch(getType(stp)){
		case PHOTON_ID:
			return massAttenuation_photon(stp);
		break;
		default: break;
	}
	return 0.0f;
}

extern "C" void  UserHook_dint_event(Step *stp){
	// prepare particle FoR
	vec3dRT v0,n1,n2;
	getVelocityVersor_A(stp,v0);
	triad(v0,n1,n2);

	// vec3dRT epsilon,n3;
	// getPolarizationDirection_A(stp,epsilon);
	// n3 = cross(v0,epsilon);

	cout<<"ciao"<<endl;


	// dispatch event
	switch(getType(stp)){
		case PHOTON_ID:
			interactionEvent_Photon(stp,v0,n1,n2);
		break;
		default: break;
	}
}


extern "C" int UserHook_source(int argc,const char *argv[]){
	cout<<magentacolor<<"User Source called with "<<argc<<" arguments"<<endl;
	for(int i=0;i<argc;i++) cout<<"\t"<<i<<": "<<argv[i]<<endl;

	int maxPrimaryPerBatch = int(atof(argv[0]));
	int numAddedPrimary = int(atof(argv[1]));
	int totNumPrimary = int(atof(argv[2]));

	int nprimbatch = totNumPrimary - numAddedPrimary < maxPrimaryPerBatch ? totNumPrimary - numAddedPrimary : maxPrimaryPerBatch;

	cout<<"Maximum number of primary particle for this batch: "<<maxPrimaryPerBatch<<endl;

	cout<<"Adding primary particles from "<<numAddedPrimary<<" to "<<numAddedPrimary+nprimbatch-1<<endl;

	for(int iprim=0;iprim<nprimbatch;iprim++){
		vec3dRT O(0,0,0);
		float phi = pluginRandUnif(0,2*M_PI);
		vec3dRT v(sin(phi),cos(phi),0);
		vec3dRT polarization(0,0,1);
		double time = 0;
		float wei = 1;
		uint64 randState = pluginGetOneSeed();
		double T = 0.511;
		addPrimary_extended(PHOTON_ID,O,v,T,randState,wei,time,polarization);
	}
	cout<<normalcolor;
	return numAddedPrimary+nprimbatch<totNumPrimary;
}
