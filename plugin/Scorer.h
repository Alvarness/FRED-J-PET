/*
 *  Scorer.h
 *  Fred
 *
 *  Created by A.Schiavi on 17/10/2016.
 *  Copyright 2016 Università di Roma LA SAPIENZA. All rights reserved.
 *
 */

#pragma once
typedef int	int32;      // 32 bit signed integer
typedef unsigned long long    uint64;     // 64 bit unsigned

#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

class Scorer {
	public:
	vector<double> data;
	string name;
	vector<Scorer *> subScorers;

	Scorer(){name="noname";};

	void init(size_t N){ data.resize(N);}

	size_t N() {return data.size();}

	void score(size_t idx,double value){// threadsafe scoring!!!
		lock_ivxl_semaphore(idx);
		data[idx]+=value;
		unlock_ivxl_semaphore(idx);
	}
	
	void set(double value){for(size_t i=0;i<data.size();i++) data[i]=value;}
	void rescale(double factor){for(size_t i=0;i<data.size();i++) data[i]*=factor;}
	void add(double value){for(size_t i=0;i<data.size();i++) data[i]+=value;}
	void reset(){this->set(0);}
	
	double min(){double vmin=data[0];for(size_t i=1;i<data.size();i++) {if(data[i]<vmin) vmin=data[i];}return vmin;};
	double max(){double vmax=data[0];for(size_t i=1;i<data.size();i++) {if(data[i]>vmax) vmax=data[i];}return vmax;};
	double sum(){double sum=0;for(size_t i=0;i<data.size();i++) sum+=data[i];return sum;};
	double sumsq(){double sumsq=0;for(size_t i=0;i<data.size();i++) sumsq+=data[i]*data[i];return sumsq;};
};

class GridScorer : public Scorer {
	i3d nn; // dimensions (or shape)
	vec3dRT hs; // spacing
	vec3dRT x0; // lowest corner coords

	public:
	GridScorer(){}

	GridScorer(i3d NN,vec3dRT HS, vec3dRT X0, const char *Name){ 
		init(NN,HS,X0,Name);
	}

	void init(i3d NN,vec3dRT HS, vec3dRT X0, const char *Name){
		nn = NN; hs = HS; x0=X0;
		name = Name;
		Scorer::init(prod(NN));
		this->reset();
	}

	void save(){
		saveGridScorer_f32(nn,hs,x0,data.data(),name.c_str());
	}
};

class PhantomGridScorer : public Scorer {
	public:
	PhantomGridScorer(size_t N,const char *Name){ 
		name=Name;
		data.resize(N);
		this->reset();
	}
	void save(){
		phantomGridScorer2map3d_32bit(data.data(),name.c_str());
	}
};