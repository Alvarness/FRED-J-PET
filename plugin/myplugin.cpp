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

		mu_rho_tot_material =calculate_mass_attenuation_coefficient(T);

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


	tuple<double, double, double, vec3dRT> angles = compton_scattering(stp);

	vec3dRT new_direction = versor(pol*sin(get<0>(angles))*cos(get<1>(angles)) + n3*sin(get<0>(angles))*sin(get<1>(angles)) + incoming*cos(get<0>(angles)));

	vec3dRT new_pol = versor(pol*get<3>(angles)[0] + n3*get<3>(angles)[1] + incoming*get<3>(angles)[2]);

	new_pol = versor(new_pol - new_direction*dot(new_pol, new_direction));

	setDirection_B(stp, new_direction);
	setKineticEnergy_B(stp, get<2>(angles));
	setPolarizationDirection_B(stp, new_pol);


	extinguishRay(stp);
	pushParticle_extended(stp,PHOTON_ID,getKineticEnergy_B(stp), new_direction, new_pol);

	vec3dRT xi;
	getPosition_A(stp,xi);
	double T_in = getKineticEnergy_A(stp);
	double T_out = getKineticEnergy_B(stp);
	double time_A = getTime_A(stp);
	vec3dRT pol2;
	getPolarizationDirection_B(stp,pol2);
	ftracks_phot << pol << '\t' << pol2 << '\t' << new_pol << endl << incoming << endl <<   xi << ' ' << T_in - T_out << ' ' << time_A << ' ' \
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

double calculate_mass_attenuation_coefficient(double x) {
   const int fNp = 36, fKstep = 0;
   const double fDelta = -1, fXmin = 0.001, fXmax = 20;
   const double fX[36] = { 0.001, 0.0015, 0.002, 0.003, 0.004,
                        0.005, 0.006, 0.008, 0.01, 0.015,
                        0.02, 0.03, 0.04, 0.05, 0.06,
                        0.08, 0.1, 0.15, 0.2, 0.3,
                        0.4, 0.5, 0.6, 0.8, 1,
                        1.25, 1.5, 2, 3, 4,
                        5, 6, 8, 10, 15,
                        20 };
   const double fY[36] = { 2024, 640.9, 277, 82.7, 34.61,
                        17.53, 10.05, 4.22, 2.204, 0.7705,
                        0.4358, 0.2647, 0.2194, 0.1997, 0.1881,
                        0.1736, 0.1635, 0.1458, 0.1331, 0.1155,
                        0.1034, 0.09443, 0.08732, 0.07668, 0.06894,
                        0.06166, 0.05611, 0.0481, 0.03848, 0.03282,
                        0.02907, 0.02641, 0.0229, 0.02069, 0.0177,
                        0.01624 };
   const double fB[36] = { -4.40753e+06, -1.43594e+06, -330729, -93454.1, -22624.3,
                        -11558.6, -4821.13, -1580.93, -624.145, -98.7565,
                        -41.749, -4.96336, -3.31761, -1.26622, -1.00752,
                        -0.557452, -0.452675, -0.287149, -0.22273, -0.141322,
                        -0.102981, -0.0788551, -0.0639988, -0.0444968, -0.033714,
                        -0.0251511, -0.0196418, -0.0131071, -0.00705407, -0.00451666,
                        -0.0031093, -0.00227614, -0.00134953, -0.000905717, -0.000367641,
                        -0.000293717 };
   const double fC[36] = { 3.90479e+09, 2.0384e+09, 1.72012e+08, 6.52626e+07, 5.56726e+06,
                        5.49841e+06, 1.2391e+06, 381003, 97389.4, 7688.37,
                        3713.13, -34.5674, 199.143, 5.99566, 19.8744,
                        2.62891, 2.60996, 0.700555, 0.587823, 0.226253,
                        0.157164, 0.084091, 0.0644722, 0.0330379, 0.0208761,
                        0.0133756, 0.00866134, 0.00440818, 0.0016448, 0.000892613,
                        0.000514744, 0.000318412, 0.000144893, 7.70152e-05, 3.06e-05,
                        5 };
   const double fD[36] = { -1.24426e+12, -1.24426e+12, -3.55833e+10, -1.98984e+10, -2.29478e+07,
                        -1.41977e+09, -1.43016e+08, -4.7269e+07, -5.98007e+06, -265017,
                        -124923, 7790.34, -6438.24, 462.625, -287.425,
                        -0.315778, -12.7294, -0.751543, -1.20523, -0.230298,
                        -0.243577, -0.0653959, -0.0523905, -0.0202697, -0.0100006,
                        -0.00628574, -0.00283544, -0.000921125, -0.00025073, -0.000125957,
                        -6.5444e-05, -2.89197e-05, -1.1313e-05, -3.09435e-06, -3.09435e-06,
                        2.39746 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

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
		// vec3dRT O(0, 0, 0);
		// vec3dRT v(1, 1, 0);
		// vec3dRT polarization, n3;
		// triad(v, polarization, n3);
		// double time = 0;
		// float wei = 1;
		// uint64 randState = pluginGetOneSeed();
		// double T = 0.511;
		// addPrimary_extended(PHOTON_ID,O,v,T,randState,wei,time,polarization);



		vec3dRT O(0,0,0);
		float theta = pluginRandUnif(0,2*M_PI);
		float phi = acos(1 - pluginRandUnif(0,2));
		vec3dRT v(sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi));
		// vec3dRT polarization(0,0,1);
		vec3dRT p1,p2;
		triad(v, p1, p2);
		double time = 0;
		float wei = 1;
		uint64 randState = pluginGetOneSeed();
		double T = 0.511;
		addPrimary_extended(PHOTON_ID,O,v,T,randState,wei,time,p1);
		addPrimary_extended(PHOTON_ID,O,-v,T,randState,wei,time,p2);
	}
	cout<<normalcolor;
	return numAddedPrimary+nprimbatch<totNumPrimary;
}
