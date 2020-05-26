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
	ftracks_phot.close();
	return 0; // return value is ignored
}

int lastTracked;

extern "C" void  UserHook_step_aft(Step *stp){
	vec3dRT xi;
	getPosition_A(stp,xi);
	double T_in = getKineticEnergy_A(stp);
	double T_out = getKineticEnergy_B(stp);
	double time_A = getTime_A(stp);
	if(xi.z<0) return;
	int type = getType(stp);

	if(type!=GHOST_ID) lastTracked = type;
	// if(type==PHOTON_ID)   ftracks_phot << xi << ' ' << T_in - T_out << ' ' << time_A << " " \
	// 	<< getUID(stp) << ' ' << getParentUID(stp) << ' ' << getAncestorUID(stp) << ' ' << getGeneration(stp) << endl;
	if(type==PHOTON_ID) {
		setStepDeltaTime(stp,getStepLength(stp)/3e10);
		// cout<<"time now "<<getUID(stp)<<' '<<getTime_B(stp)<<endl;
	}
}

extern "C" void  UserHook_step_final(Step *stp){
	vec3dRT xi;
	getPosition_A(stp,xi);
	double T = getKineticEnergy_A(stp);

	// if(lastTracked==PHOTON_ID) {
	// 	ftracks_phot<<xi<<' '<<T<<' '<<endl;
	// 	ftracks_phot<<endl;
	// }
}

float massAttenuation_photon(Step *stp){
	// compute  the total mass attenuation coeff for the photon
	double mu_rho_tot_material = 0.;

	int I_mat = getImat_A(stp); //get material index at initial point A
	double T = getKineticEnergy_A(stp);

	string currReg = getRegion_A(stp);

	if(currReg.rfind("scint")==0){
		// cout<<currReg<<endl;
		// mu_rho_tot_material = 0.09443;
		mu_rho_tot_material = 4.5443;

	}

	return mu_rho_tot_material;
}



void interactionEvent_Photon(Step *stp){
	vec3dRT incoming;
	getVelocityVersor_A(stp,incoming);
	incoming = versor(incoming);
	vec3dRT pol;
	getPolarizationDirection_A(stp,pol);
	pol = versor(pol);
	vec3dRT n3 = versor(cross(incoming, pol));

	cout << dot(incoming, pol)<< " vectors -> " << incoming << " " << pol << endl;

	tuple<double, double, double, vec3dRT> angles = compton_scattering(stp);

	vec3dRT new_direction = versor(pol*sin(get<0>(angles))*cos(get<1>(angles)) + n3*sin(get<0>(angles))*sin(get<1>(angles)) + incoming*cos(get<0>(angles)));

	vec3dRT new_pol = versor(pol*get<3>(angles)[0] + n3*get<3>(angles)[1] + incoming*get<3>(angles)[2]);

	new_pol = versor(new_pol - new_direction*dot(new_pol, new_direction));

	setDirection_B(stp, new_direction);
	setKineticEnergy_B(stp, get<2>(angles));
	setPolarizationDirection_B(stp, new_pol);

	// pushParticle(stp,PHOTON_ID,getKineticEnergy_B(stp), new_direction);
	// extinguishRay(stp);

	cout << get<0>(angles) << " " << get<1>(angles) << " " << get<2>(angles) << endl;

	cout << acos(dot(incoming, new_direction)) << " ";

	cout << acos(dot(versor((new_direction - incoming*dot(new_direction, incoming))), pol)) << endl;

	cout << dot(new_pol, new_direction) << endl;

	cout << "/*--------------------*/" << endl;


	// cout<<"v0 "<<vnew<<endl;
	// cout<<"ciao "<<getUID(stp)<<' '<<getKineticEnergy_A(stp)<<' '<<pol<<' '<<getTime_A(stp)<<endl;

	// setDirection_B(stp, new_direction);
	// setKineticEnergy_B(stp, get<2>(angles));

	vec3dRT xi;
	getPosition_A(stp,xi);
	double T_in = getKineticEnergy_A(stp);
	double T_out = getKineticEnergy_B(stp);
	double time_A = getTime_A(stp);
	vec3dRT pol2;
	getPolarizationDirection_B(stp,pol2);
	ftracks_phot << pol << '\t' << pol2 << '\t' << new_pol << endl <<   xi << ' ' << T_in - T_out << ' ' << time_A << ' ' \
		<< getUID(stp) << ' ' << getParentUID(stp) << ' ' << getAncestorUID(stp) << ' ' << getGeneration(stp) << endl << endl;

}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tuple<double, double, double, vec3dRT> compton_scattering(Step *stp){

	/*----------  Generate Compton scattering angle  ----------*/

	double energy = getKineticEnergy_A(stp);
	energy = energy*1000.0;

    double Mec2 = 510.99895000;

    double E0 = energy / Mec2;

    double epsilon, onecos_t, SinSqrTh, greject;
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
        SinSqrTh = onecos_t*(2 - onecos_t);
        greject = 1 - (SinSqrTh*epsilon)/(1 + pow(epsilon, 2));

    } while( greject < r3);


    double compton_angle = acos(1 - onecos_t);
    double scattered_energy = epsilon*energy;


    /*----------  Generate azimuthal scattering angle   ----------*/


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

	/*----------  Generate polarization direction of scattered photon ----------*/


	double rand1 = getRandUnif(stp);
	double rand2 = getRandUnif(stp);

	double theta;

	double sinTheta = sin(compton_angle);
	double cosTheta = cos(compton_angle);

	double cosPhi = cos(azimutlal_angle);
	double sinPhi = sin(azimutlal_angle);

	double cosSqrPhi = cosPhi*cosPhi;
	double normalisation = sqrt(1. - cosSqrPhi*SinSqrTh);

	if (rand1 < (epsilon+1.0/epsilon-2) / (2.0*(epsilon+1.0/epsilon)-4.0*SinSqrTh*cosSqrPhi) ){

    	if (rand2<0.5) theta = M_PI/2.0;

    	else theta = 3.0*M_PI/2.0;

    } else {

     	if (rand2<0.5) theta = 0.;

    	else theta = M_PI;
    }

    double cosBeta = cos(theta);
 	double sinBeta = sqrt(1-cosBeta*cosBeta);

 	double xParallel = normalisation*cosBeta;
	double yParallel = -(SinSqrTh*cosPhi*sinPhi)*cosBeta/normalisation;
	double zParallel = -(cosTheta*sinTheta*cosPhi)*cosBeta/normalisation;
	double xPerpendicular = 0.;
	double yPerpendicular = (cosTheta)*sinBeta/normalisation;
	double zPerpendicular = -(sinTheta*sinPhi)*sinBeta/normalisation;

	double xTotal = (xParallel + xPerpendicular);
	double yTotal = (yParallel + yPerpendicular);
	double zTotal = (zParallel + zPerpendicular);

	vec3dRT newPolarizationDirection (xTotal, yTotal, zTotal);

   return tuple<double, double, double, vec3dRT> {compton_angle, azimutlal_angle - M_PI, scattered_energy/1000.0, newPolarizationDirection};

}


// void SystemOfRefChange(	vec3dRT& incoming,
// 						vec3dRT& scattered,
// 						vec3dRT& pol,
// 						vec3dRT& new_pol){
//   // direction0 is the original photon direction ---> z
//   // polarization0 is the original photon polarization ---> x
//   // need to specify y axis in the real reference frame ---> y
// 	vec3dRT Axis_Z0 = incoming;
// 	vec3dRT Axis_X0 = pol;
// 	vec3dRT Axis_Y0 = versor(cross(incoming, pol));

//   	double direction_x = scattered[0];
//   	double direction_y = scattered[1];
//   	double direction_z = scattered[2];

//   	scattered = versor(Axis_X0*direction_x + Axis_Y0*direction_y + Axis_Z0*direction_z);

//   	double polarization_x = new_pol[0];
//   	double polarization_y = new_pol[1];
//   	double polarization_z = new_pol[2];

//   	new_pol = versor(Axis_X0*polarization_x + Axis_Y0*polarization_y + Axis_Z0*polarization_z);

//   	// cout << dot(incoming, pol) << " "  << " " << dot(new_pol, scattered) << " zz" << endl;

// }

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
	vec3dRT v0;
	getVelocityVersor_A(stp,v0);
	// // triad(v0,n1,n2);

	vec3dRT epsilon;
	getPolarizationDirection_A(stp,epsilon);
	// n3 = cross(v0,epsilon);

	cout<<"dint "<< dot(v0, epsilon) << endl;


	// dispatch event
	switch(getType(stp)){
		case PHOTON_ID:
			interactionEvent_Photon(stp);
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
