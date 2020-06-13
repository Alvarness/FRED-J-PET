/*
 *  vector3d.h
 *  Fred
 *
 *  Created by A.Schiavi on 02/10/12.
 *  Copyright 2010 Università di Roma LA SAPIENZA. All rights reserved.
 *
 */

#ifndef VECTOR3D_H
#define VECTOR3D_H

template <typename T>
class vector3d {
private:	
public:
	union {
		T v[3];
		struct { T x,y,z;}; // CART
		struct { T Alpha,R,Z;}; // CYL
		struct { T phi,theta,rho;}; // SPHER
	};
	
	vector3d(){};
	
	vector3d(T a,T b,T c){v[0]=a;v[1]=b;v[2]=c;};
	
	vector3d(T a){v[0]=a;v[1]=a;v[2]=a;}
	
	vector3d(int n){v[0]=n;v[1]=n;v[2]=n;}
	
 	vector3d(T *a){v[0]=a[0];v[1]=a[1];v[2]=a[2];}
 	
	inline T & operator[](int j){return v[j];};
	inline const T & operator[](int j) const {return v[j];};
	
	inline void set(T x,T y,T z){v[0]=x;v[1]=y;v[2]=z;}
	
	
	inline vector3d& operator=(const vector3d &a){for (int j=0;j<3;j++) v[j]=a.v[j];return *this;}
	inline vector3d& operator+=(T a){v[0]+=a;v[1]+=a;v[2]+=a;return *this;};
	inline vector3d& operator-=(T a){v[0]-=a;v[1]-=a;v[2]-=a;return *this;};
	inline vector3d& operator*=(T a){v[0]*=a;v[1]*=a;v[2]*=a;return *this;};
	inline vector3d& operator/=(T a){v[0]/=a;v[1]/=a;v[2]/=a;return *this;};
	inline vector3d&  operator=(T a){v[0] =a;v[1] =a;v[2] =a;return *this;};


	inline vector3d& operator++(){v[0]++;v[1]++;v[2]++;return *this;};
	inline vector3d& operator--(){v[0]--;v[1]--;v[2]--;return *this;};
	inline vector3d& operator++(int ){v[0]++;v[1]++;v[2]++;return *this;};
	inline vector3d& operator--(int ){v[0]--;v[1]--;v[2]--;return *this;};
	
	
	inline vector3d operator-() const{return vector3d(-v[0],-v[1],-v[2]);}
	inline vector3d operator+(const vector3d &a) const{return vector3d(v[0]+a.v[0],v[1]+a.v[1],v[2]+a.v[2]);}
	inline vector3d operator-(const vector3d &a) const{return vector3d(v[0]-a.v[0],v[1]-a.v[1],v[2]-a.v[2]);}
	inline vector3d operator*(const vector3d &a) const{return vector3d(v[0]*a.v[0],v[1]*a.v[1],v[2]*a.v[2]);}
	inline vector3d operator/(const vector3d &a) const{return vector3d(v[0]/a.v[0],v[1]/a.v[1],v[2]/a.v[2]);}

	inline vector3d& operator+=(const vector3d &a){v[0]+=a.v[0];v[1]+=a.v[1];v[2]+=a.v[2];return *this;}
	inline vector3d& operator-=(const vector3d &a){v[0]-=a.v[0];v[1]-=a.v[1];v[2]-=a.v[2];return *this;}
	inline vector3d& operator*=(const vector3d &a){v[0]*=a.v[0];v[1]*=a.v[1];v[2]*=a.v[2];return *this;}
	inline vector3d& operator/=(const vector3d &a){v[0]/=a.v[0];v[1]/=a.v[1];v[2]/=a.v[2];return *this;}
    
	inline vector3d operator+(T a) const{vector3d x;for (int j=0;j<3;j++) x.v[j]=v[j]+a;return x;}
	inline vector3d operator-(T a) const{vector3d x;for (int j=0;j<3;j++) x.v[j]=v[j]-a;return x;}
	inline vector3d operator*(T a) const{vector3d x;for (int j=0;j<3;j++) x.v[j]=v[j]*a;return x;}
	inline vector3d operator/(T a) const{vector3d x;for (int j=0;j<3;j++) x.v[j]=v[j]/a;return x;}
	
};


template <typename T>
ostream& operator<<(ostream &os,const vector3d<T> &v){
	// space delimited output
	os<<v[0]<<' '<<v[1]<<' '<<v[2];
	return os;
}

// NB: functions hereafter are defined  inline for max performance
template <typename T>
inline vector3d<T> operator*(const T r,const vector3d<T> &v){return vector3d<T>(r*v[0],r*v[1],r*v[2]);}

template <typename T>  
inline T dot (const vector3d<T> &a,const vector3d<T> &b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}

template <typename T>
inline vector3d<T> cross(const vector3d<T> &a,const vector3d<T> &b){return vector3d<T>(a.v[1]*b.v[2]-a.v[2]*b.v[1],a.v[2]*b.v[0]-a.v[0]*b.v[2],a.v[0]*b.v[1]-a.v[1]*b.v[0]);}

template <typename T>
inline T norm(const vector3d<T> &a){return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);}

template <typename T>
inline vector3d<T> versor(const vector3d<T> &a){
    vector3d<T> x = a;
    T len = norm(a);
    if(len>0) x/=len;
    return x;
}

template <typename T>
void triad(const vector3d<T> &a, vector3d<T> &n1, vector3d<T> &n2){
  vector3d<T> v=versor(a);
  if(v[0]!=0 || v[2]!=0) {
    n1=vector3d<T>(v[2],0,-v[0]);
  }else if(v[1]!=0){
    n1=vector3d<T>(0,v[2],-v[1]);
  } else { //degenerate vector
    n1=0;
    n2=0;
    return;
  }
  n1=versor(n1);
  n2=cross(v,n1);
}	

typedef vector3d<real64> vec3d;
typedef vector3d<real32> vec3dRT;

#endif // VECTOR3D_H
