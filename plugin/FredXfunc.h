#ifndef FREDXFUNC_H
#define FREDXFUNC_H

/*
 *  FredXfunc.h
 *
 *  asch 2015-18
 *
 */

#include <cmath>
#include <iostream>
using namespace std;
#include "types.h"
#include "ParticleID.h"
namespace fred {
#include "vector3d.h"
}
#include "idx3d.h"
#include "Arr3d.h"

const char normalcolor[]  = "\e[0m", blackcolor[]   = "\e[0;30m",redcolor[]     = "\e[0;31m";
const char greencolor[]   = "\e[0;32m",yellowcolor[]  = "\e[0;33m",bluecolor[]    = "\e[0;34m";
const char magentacolor[] = "\e[0;35m",cyancolor[]    = "\e[0;36m",whitecolor[]   = "\e[0;37m";

const int SKIP=1,APPLY=0;
struct Step;

using namespace fred;

extern "C" {

// getter functions at the begin of the step, i.e. position A
void	getPosition_A(Step *stp, vec3dRT &pos);
void	getDirection_A(Step *stp, vec3dRT &vel);
double	getKineticEnergy_A(Step *stp);
double	getMomentum_A(Step *stp);
void    getRegion_A(Step *stp, char *name);
void	getPolarizationDirection_A(Step *stp, vec3dRT &polarization);
double	getTime_A(Step *stp);

float32	getNlambda(Step *stp);

// getter functions at the end of the step, i.e. position B
void	getPosition_B(Step *stp, vec3dRT &pos);
void	getDirection_B(Step *stp, vec3dRT &vel);
double	getKineticEnergy_B(Step *stp);
double	getMomentum_B(Step *stp);
void    getRegion_B(Step *stp, char *name);
void	getPolarizationDirection_B(Step *stp, vec3dRT &polarization);
double	getTime_B(Step *stp);
void	setPolarizationDirection_B(Step *stp, vec3dRT &polarization);

void 	setStepDeltaTime(Step *stp,double dt);


// setter functions at the end of the step, i.e. position B
void	setPosition_B(Step *stp, vec3dRT &pos);
void	setDirection_B(Step *stp,vec3dRT v) ;
void	setKineticEnergy_B(Step *stp,double T);
void	setMomentum_B(Step *stp,double p);
void	setNlambda_B(Step *stp,float nl);


int32  getType(Step *stp); // particle ID (using PDG 2006), e.g. PROTON = 2212, etc.

float  getParticle_m(Step *stp); // particle rest mass (MeV/c^2)
float  getParticle_Z(Step *stp); // electric charge number or atomic number, i.e. number of protons (P)
float  getParticle_A(Step *stp); // mass number, i.e. number of nucleons (P+N)

int32  getUID(Step *stp); // particle UID (unique identifier): this labels each single particle produced and tracked
int32  getParentUID(Step *stp); // parent UID, i.e. previous generation
int32  getAncestorUID(Step *stp); // ancestor UID, i.e. first generation
int32  getGeneration(Step *stp); // get particle generation: 1 = primary, 2 = secondary ,...


float	getStepLength(Step *stp);
float	getRangeStep(Step *stp);
float	getMassDensity(Step *stp);

int32	getPBindex(Step *stp);

float	getRandUnif(Step *stp);
float	getRandGauss(Step *stp);

short	getHU(Step *stp);


void extinguishRay(Step *stp); // kill present track immediately


int32 getPluginBuffSize();
void * getPluginBuff(Step *stp);

uint64 getPluginRandSeed();

int	pushParticle(Step *stp, int type, float T, vec3dRT v);
int	pushParticleAtArbitraryPosition(Step *stp, int type, float T, vec3dRT v, vec3dRT x);

int	addPrimary(int type,vec3dRT x, vec3dRT v,float T,uint64 randState, float wei);

int	addPrimary_extended(int type,vec3dRT x, vec3dRT v,float T,uint64 randState, float wei,double time, vec3dRT polarization);

float  pluginRandUnif(); // samples a uniform distribution
float  pluginRandGauss(); // samples a gaussian (normal) distribution
uint64 pluginGetOneSeed(); // get one 64-bit seed for random generators or random state of a particle


void	setLocallyDepositedEnergy(Step *stp,float Edep);
void	setLostEnergy(Step *stp,float Elost);

size_t	getVxlIdx(Step *stp);

bool	getBoolParam(const char *pname,bool defVal);
int		getIntParam(const char *pname,int defVal);
double	getFloatParam(const char *pname,double defVal);
char *	getStringParam(const char *pname,const char* defStr);
void getVec3dRT(const char *pname,vec3dRT &vec, vec3dRT default_vec);

// region functions
int  getRegion_Index(const char *regID); // obtain internal index of region named regID

void getRegion_Dims(int ireg, i3d &nn); // num of voxels in region
void getRegion_Spacing(int ireg, vec3dRT &hs); // voxel spacing
void getRegion_Origin(int ireg, vec3dRT &O); // region origin
void getRegion_Offset(int ireg, vec3dRT &x0); // region "lowest" corner

float getRegion_Mass(int ireg, int ivxl); // mass of the voxel
float getRegion_Density(int ireg, int ivxl); // density of the voxel
int   getRegion_Imat(int ireg, int ivxl); // material index of the voxel

// material functions
int getImat_A(Step *stp); // get material index at initial point A
int getImat_B(Step *stp); // get material index at final   point B

float getMat_Zmean(int imat); // mean Z for the material
float getMat_Amean(int imat); // mean A for the material
float getMat_RelStopPow(int imat);

// material nuclear composition: each element in the chemical formula corresponds to a different nucleus
int getMatNumElements(int imat); // num of elements in the material
float getMat_Z(int imat,int iel); // num of protons in element nucleus (P)
float getMat_A(int imat,int iel); // num of nucleons in element nucleus (N+P)
float getMat_m(int imat,int iel); // mass of element nucleus
float getMat_w(int imat,int iel); // mass weight of the element nucleus in the material
float getMat_x(int imat,int iel); // number fraction of the element nucleus in the material


float get_dEds(Step *stp, float T);

void lock_semaphore(int i);
void unlock_semaphore(int i);

void lock_ivxl_semaphore(size_t ivxl);
void unlock_ivxl_semaphore(size_t ivxl);


void getPhantomDims(int32 _nn[3]);
void getPhantomSpacing(float _v[3]);
void getPhantomOffset(float _v[3]);
void phantomGridScorer2map3d_32bit(double *data,const char *fname);

void saveGridScorer_f32(i3d &nn,vec3dRT &hs,vec3dRT &x0,double *data,const char *fname);

float getPhantom_Mass(int ivxl);
float getPhantom_Density(int ivxl);
short getPhantom_HU(int ivxl);

int set_NamedBool(const char *varname, bool *val);
int get_NamedBool(const char *varname, bool *val);

int get_NamedArr3d_float64(const char *varname, Arr3d<double> **arr);
int get_NamedArr3d_float32(const char *varname, Arr3d<float> **arr);

int getNumPB();
int getNumPrimary(int jpb);
int getNumTotPrimary();

// DEPRECATED FUNCTIONS that will be sooner or later removed!!!!!!
// try to avoid them in your codes: they are here just for backcompatibility... and for a short time!
void	getPosition(Step *stp, vec3dRT &pos); // << deprecated in favor of getPosition_B
void	getInitialPosition(Step *stp, vec3dRT &pos); // << deprecated in favor of getPosition_A
void	getFinalPosition(Step *stp, vec3dRT &pos); // << deprecated in favor of getPosition_B

void	getVelocityVersor_A(Step *stp, vec3dRT &vel); // << deprecated in favor of getDirection_A
void	getVelocityVersor_B(Step *stp, vec3dRT &vel); // << deprecated in favor of getDirection_B
void	setVelocityVersor_B(Step *stp,vec3dRT v) ;// << deprecated in favor of setDirection_B

double	getKineticEnergy(Step *stp); // << deprecated in favor of getKineticEnergy_B
void	setKineticEnergy(Step *stp,double T); // << deprecated in favor of setKineticEnergy_B

double	getMomentum(Step *stp); // << deprecated in favor of getMomentum_B
void	setMomentum(Step *stp,double p); // << deprecated in favor of setMomentum_B

float32	getNlambda(Step *stp); // << deprecated in favor of getNlambda_B
void	setNlambda(Step *stp,float nl); // << deprecated in favor of setNlambda_B

int32  getParticle_generation(Step *stp);

// END OF DEPRECATED FUNCTIONS LIST
#define getPosition getPosition_B
#define getInitialPosition getPosition_A
#define getFinalPosition getPosition_B
#define getVelocityVersor_A getDirection_A
#define getVelocityVersor_B getDirection_B
#define setVelocityVersor setDirection_B
#define getKineticEnergy getKineticEnergy_B
#define setKineticEnergy setKineticEnergy_B
#define getMomentum getMomentum_B
#define setMomentum setMomentum_B
// #define getNlambda getNlambda_B
// #define setNlambda setNlambda_B

}

inline float pluginRandUnif(float amin,float amax) // samples a uniform distribution in range [amin,amax)
{ return amin+(amax-amin)*pluginRandUnif();}
inline float pluginRandGauss(float mean,float stdev) // samples a gaussian distribution with mean and stdev
{ return mean+stdev*pluginRandGauss();}

inline vec3dRT getPosition_A(Step *stp){vec3dRT v; getPosition_A(stp,v); return v;}
inline vec3dRT getPosition_B(Step *stp){vec3dRT v; getPosition_B(stp,v); return v;}

inline string getRegion_A(Step *stp){char buff[256]; getRegion_A(stp,buff); return string(buff);}
inline string getRegion_B(Step *stp){char buff[256]; getRegion_B(stp,buff); return string(buff);}

#endif
