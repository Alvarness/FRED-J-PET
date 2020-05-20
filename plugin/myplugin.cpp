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

	//cout<<getUID(stp)<<' '<<getParentUID(stp)<<' '<<getAncestorUID(stp)<<' '<<type<<' '<<getGeneration(stp)<<endl;
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
		// cout<<currReg<<endl;
		mu_rho_tot_material = 0.09443;

	}

	return mu_rho_tot_material;
}



void interactionEvent_Photon(Step *stp,const vec3dRT &v0,const vec3dRT &n1,const vec3dRT &n2){
	vec3dRT vnew;
	getVelocityVersor_A(stp,vnew);
	vec3dRT pol;
	getPolarizationDirection_A(stp,pol);
	// vec3dRT Y(0,1,0);
	cout <<vnew << endl;
	cout <<pol << endl;
	cout << dot(vnew, pol) << endl;
	tuple<double, double, double> angles = compton_scattering(stp);

	// vec3dRT new_direction = pol*cos(get<1>(angles)) + Y*sin(get<1>(angles)) + vnew*cos(get<0>(angles));

	cout << get<0>(angles) << " " << get<1>(angles) << " " << get<2>(angles) << endl;

	// cout << acos(dot(vnew, new_direction)) << endl;

	// cout << acos(dot(new_direction - new_direction*dot(vnew, new_direction), pol)) << endl;

	vec3dRT rotated = rotate(vnew, get<0>(angles), cross(vnew,pol));
	// rotated = rotate(rotated, get<1>(angles), vnew);

	cout << acos(dot(rotated, vnew)) << endl;

	// setDirection_B(stp, rotated) ;
	// cout<<"v0 "<<vnew<<endl;
	//pushParticle(stp,ELECTRON_ID,getKineticEnergy_A(stp)/4,v0);
	// cout<<"ciao "<<getUID(stp)<<' '<<getKineticEnergy_A(stp)<<' '<<pol<<' '<<getTime_A(stp)<<endl;
	// exit(0);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vec3dRT rotate(const vec3dRT v, float angle, const vec3dRT axis){
	float sa = sin(angle), ca = cos(angle);
	float ll = norm(axis);
	float dx = axis[0]/ll, dy = axis[1]/ll, dz = axis[2]/ll;

	// ca+(1-ca)*dx*dx, (1-ca)*dx*dy-sa*dz, (1-ca)*dx*dz+sa*dy,
	// (1-ca)*dy*dx+sa*dz, ca+(1-ca)*dy*dy, (1-ca)*dy*dz-sa*dx,
	// (1-ca)*dz*dx-sa*dy, (1-ca)*dz*dy+sa*dx, ca+(1-ca)*dz*dz,

	return vec3dRT(	  (ca+(1-ca)*dx*dx)*v[0]	+ ((1-ca)*dx*dy-sa*dz)*v[1] 	+ ((1-ca)*dx*dz+sa*dy)*v[2],
				   ((1-ca)*dy*dx+sa*dz)*v[0] 	+ (ca+(1-ca)*dy*dy)*v[1] 		+ ((1-ca)*dy*dz-sa*dx)*v[2],
				   ((1-ca)*dz*dx-sa*dy)*v[0] 	+ ((1-ca)*dz*dy+sa*dx)*v[1] 	+ (ca+(1-ca)*dz*dz)*v[2]);

}

tuple<double, double, double> compton_scattering(Step *stp){

	double energy = getKineticEnergy_A(stp);
	energy = energy*1000.0;

    double Mec2 = 510.99895000;

    double E0 = energy / Mec2;

    double epsilon, onecos_t, sin2_t, greject;
    double epsilon0 = 1 / (1 + 2. * E0);
    double epsilon0sq = pow(epsilon0,2);

    double alpha1 = log(1 / epsilon0);
    double alpha2 = (1. - epsilon0sq) / 2;

    double r1, r2, r3;

    do{
        r1 = getRandUnif(stp);
        r2 = getRandUnif(stp);
        r3 = getRandUnif(stp);

        if(r1 < alpha1/(alpha1 + alpha2)){

            epsilon = exp(-alpha1 *r2);

        } else{

            epsilon = sqrt( epsilon0sq + (1 - epsilon0sq) * r2);
        }

        onecos_t = (1 - epsilon)/(epsilon);
        sin2_t = onecos_t*(2 - onecos_t);
        greject = 1 - (sin2_t*epsilon)/(1 + pow(epsilon, 2));

    } while( greject < r3);


    double compton_angle = acos(1 - onecos_t);
    double scattered_energy = epsilon*energy;




    double b = scattered_energy + 1.0/scattered_energy;
    double a = (sin(compton_angle),2);

    double X;
    double prob;
    double r4;

    do{
        r4 = getRandUnif(stp) * 2 * M_PI;
        X = getRandUnif(stp)*a*b;

        double b = epsilon + 1.0/epsilon;
    	double a = pow(sin(compton_angle),2);

        prob = 1 - pow(cos(r4),2)*2*a/b;
    } while(X > prob);

   double azimutlal_angle = r4;

   return tuple<double, double, double> {compton_angle, azimutlal_angle, scattered_energy/1000.0};

}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

extern "C" float  UserHook_Mass_Attenuation(Step *stp){
	// NB: returns the mass attenuation coefficient of interaction in the material
	switch(getType(stp)){
		case PHOTON_ID:
			return 1;//massAttenuation_photon(stp);
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

	// cout<<"ciao"<<endl;


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
		vec3dRT O(0, 0, 0);
		vec3dRT v(1, 1, 0);
		vec3dRT polarization, n3;
		triad(v, polarization, n3);
		double time = 0;
		float wei = 1;
		uint64 randState = pluginGetOneSeed();
		double T = 0.511;
		addPrimary_extended(PHOTON_ID,O,v,T,randState,wei,time,polarization);



		// vec3dRT O(0,0,0);
		// float theta = pluginRandUnif(0,2*M_PI);
		// float phi = acos(1 - pluginRandUnif(0,2));
		// vec3dRT v(sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi));
		// vec3dRT polarization(0,0,1);
		// double time = 0;
		// float wei = 1;
		// uint64 randState = pluginGetOneSeed();
		// double T = 0.511;
		// addPrimary_extended(PHOTON_ID,O,v,T,randState,wei,time,polarization);
		// addPrimary_extended(PHOTON_ID,O,-v,T,randState,wei,time,polarization);
	}
	cout<<normalcolor;
	return numAddedPrimary+nprimbatch<totNumPrimary;
}
