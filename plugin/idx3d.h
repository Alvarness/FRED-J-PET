#include <iostream>
#include <cmath>
#include "types.h"


using namespace std;

#ifndef IDX3D_H
#define IDX3D_H

namespace fred {

class i3d {
private:  
  static int32 outMode;
public:
    union {
		int32 v[3];
		struct { int32 x,y,z;}; // CART
		struct { int32 Alpha,R,Z;}; // CYL
		struct { int32 phi,theta,rho;}; // SPHER
	};

	i3d(){};

	i3d(int32 a,int32 b, int32 c){
 		v[0]=a;v[1]=b;v[2]=c;
 	};
	
	i3d(const i3d &a){
 		v[0]=a.v[0];v[1]=a.v[1];v[2]=a.v[2];
 	};
	
	i3d(int32 n){
 		v[0]=v[1]=v[2]=n;
 	}
	
 	i3d(int32 *a){
 		v[0]=a[0];v[1]=a[1];v[2]=a[2];
 	}
 	
	inline int32 & operator[](int j){return v[j];};
	inline const int32 & operator[](int j) const {return v[j];};
	
	inline void set(int32 x,int32 y,int32 z){v[0]=x;v[1]=y;v[2]=z;}
	
	
	inline i3d& operator=(const i3d &a){v[0]=a.v[0];v[1]=a.v[1];v[2]=a.v[2];return *this;}
    
    inline i3d& operator+=(int32 a){v[0]+=a;v[1]+=a;v[2]+=a;return *this;};
    inline i3d& operator-=(int32 a){v[0]-=a;v[1]-=a;v[2]-=a;return *this;};
    inline i3d& operator*=(int32 a){v[0]*=a;v[1]*=a;v[2]*=a;return *this;};
    inline i3d& operator/=(int32 a){v[0]/=a;v[1]/=a;v[2]/=a;return *this;};
    inline i3d&  operator=(int32 a){v[0] =a;v[1] =a;v[2] =a;return *this;};
    
    inline i3d& operator++(){v[0]++;v[1]++;v[2]++;return *this;};
    inline i3d& operator--(){v[0]--;v[1]--;v[2]--;return *this;};
    inline i3d& operator++(int ){v[0]++;v[1]++;v[2]++;return *this;};
    inline i3d& operator--(int ){v[0]--;v[1]--;v[2]--;return *this;};
    
    inline i3d& operator+=(const i3d &a){v[0]+=a.v[0];v[1]+=a.v[1];v[2]+=a.v[2];return *this;}
    inline i3d& operator-=(const i3d &a){v[0]-=a.v[0];v[1]-=a.v[1];v[2]-=a.v[2];return *this;}
    inline i3d& operator*=(const i3d &a){v[0]*=a.v[0];v[1]*=a.v[1];v[2]*=a.v[2];return *this;}
    inline i3d& operator/=(const i3d &a){v[0]/=a.v[0];v[1]/=a.v[1];v[2]/=a.v[2];return *this;}
    
    inline i3d operator-() const {return i3d(-v[0],-v[1],-v[2]);}
    
    friend i3d operator+(i3d lhs,       // passing first arg by value helps optimize chained a+b+c
                    const i3d& rhs) // alternatively, both parameters may be const references.
  	{
    	 return lhs += rhs; // reuse compound assignment and return the result by value
  	}
    friend i3d operator-(i3d lhs,const i3d& rhs){return lhs -= rhs;}
    friend i3d operator*(i3d lhs,const i3d& rhs){return lhs *= rhs;}
    friend i3d operator/(i3d lhs,const i3d& rhs){return lhs /= rhs;}

    friend i3d operator*(int32 lhs,const i3d& rhs){return rhs*lhs;}
    
	inline bool operator !=(i3d const& a) {return v[0]!=a.v[0]||v[1]!=a.v[1]||v[2]!=a.v[2];}
    inline bool operator ==(i3d const& a) {return !(*this!=a);}
    
    friend uint64 prod(i3d i) {return 1UL*i.x*i.y*i.z;}
    
    static void setOutMode(int kmode){outMode=kmode;}
    static int  getOutMode(){return outMode;}


  inline bool any_lt(const int32 a) const {return v[0]< a || v[1]< a || v[2]< a;}
  inline bool any_le(const int32 a) const {return v[0]<=a || v[1]<=a || v[2]<=a;}
  inline bool any_eq(const int32 a) const {return v[0]==a || v[1]==a || v[2]==a;}
  inline bool any_gt(const int32 a) const {return v[0]> a || v[1]> a || v[2]> a;}
  inline bool any_ge(const int32 a) const {return v[0]>=a || v[1]>=a || v[2]>=a;}
  inline bool all_lt(const int32 a) const {return v[0]< a && v[1]< a && v[2]< a;}
  inline bool all_le(const int32 a) const {return v[0]<=a && v[1]<=a && v[2]<=a;}
  inline bool all_eq(const int32 a) const {return v[0]==a && v[1]==a && v[2]==a;}
  inline bool all_gt(const int32 a) const {return v[0]> a && v[1]> a && v[2]> a;}
  inline bool all_ge(const int32 a) const {return v[0]>=a && v[1]>=a && v[2]>=a;}

	
};

ostream& operator<<(ostream &os,const i3d &v);

inline i3d min(const i3d &a, const i3d &b){
  i3d x;
  for (int j=0;j<3;j++) x[j]=a[j]<b[j]? a[j]: b[j];
  return x;
}

inline i3d max(const i3d &a, const i3d &b){
  i3d x;
  for (int j=0;j<3;j++) x[j]=a[j]>b[j]? a[j]: b[j];
  return x;
}

inline bool operator==(const i3d &a,const i3d &b){
    if(a.x != b.x || a.y != b.y || a.z != b.z ) return false;
	return true;
}

inline bool operator!=(const i3d &a,const i3d &b){
    if(a.x != b.x || a.y != b.y || a.z != b.z ) return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
class i3ds {
private:  
  static int32 outMode;
public:
    union {
		int16 v[3];
		struct { int16 x,y,z;}; // CART
		struct { int16 Alpha,R,Z;}; // CYL
		struct { int16 phi,theta,rho;}; // SPHER
	};

	i3ds(){};

	i3ds(int16 a,int16 b, int16 c){
 		v[0]=a;v[1]=b;v[2]=c;
 	};
	
	i3ds(const i3ds &a){
 		v[0]=a.v[0];v[1]=a.v[1];v[2]=a.v[2];
 	};
	
	i3ds(int16 n){
 		v[0]=v[1]=v[2]=n;
 	}
	
 	i3ds(int16 *a){
 		v[0]=a[0];v[1]=a[1];v[2]=a[2];
 	}
 	
	inline int16 & operator[](int j){return v[j];};
	inline const int16 & operator[](int j) const {return v[j];};
	
	inline void set(int16 x,int16 y,int16 z){v[0]=x;v[1]=y;v[2]=z;}
	
	
	inline i3ds& operator=(const i3ds &a){v[0]=a.v[0];v[1]=a.v[1];v[2]=a.v[2];return *this;}
    
    inline i3ds& operator+=(int16 a){v[0]+=a;v[1]+=a;v[2]+=a;return *this;};
    inline i3ds& operator-=(int16 a){v[0]-=a;v[1]-=a;v[2]-=a;return *this;};
    inline i3ds& operator*=(int16 a){v[0]*=a;v[1]*=a;v[2]*=a;return *this;};
    inline i3ds& operator/=(int16 a){v[0]/=a;v[1]/=a;v[2]/=a;return *this;};
    inline i3ds&  operator=(int16 a){v[0] =a;v[1] =a;v[2] =a;return *this;};
    
    inline i3ds& operator++(){v[0]++;v[1]++;v[2]++;return *this;};
    inline i3ds& operator--(){v[0]--;v[1]--;v[2]--;return *this;};
    inline i3ds& operator++(int ){v[0]++;v[1]++;v[2]++;return *this;};
    inline i3ds& operator--(int ){v[0]--;v[1]--;v[2]--;return *this;};
    
    inline i3ds& operator+=(const i3ds &a){v[0]+=a.v[0];v[1]+=a.v[1];v[2]+=a.v[2];return *this;}
    inline i3ds& operator-=(const i3ds &a){v[0]-=a.v[0];v[1]-=a.v[1];v[2]-=a.v[2];return *this;}
    inline i3ds& operator*=(const i3ds &a){v[0]*=a.v[0];v[1]*=a.v[1];v[2]*=a.v[2];return *this;}
    inline i3ds& operator/=(const i3ds &a){v[0]/=a.v[0];v[1]/=a.v[1];v[2]/=a.v[2];return *this;}
    
    inline i3ds operator-() const {return i3ds(-v[0],-v[1],-v[2]);}
    
    friend i3ds operator+(i3ds lhs,       // passing first arg by value helps optimize chained a+b+c
                    const i3ds& rhs) // alternatively, both parameters may be const references.
  	{
    	 return lhs += rhs; // reuse compound assignment and return the result by value
  	}
    friend i3ds operator-(i3ds lhs,const i3ds& rhs){return lhs -= rhs;}
    friend i3ds operator*(i3ds lhs,const i3ds& rhs){return lhs *= rhs;}
    friend i3ds operator/(i3ds lhs,const i3ds& rhs){return lhs /= rhs;}

    friend i3ds operator*(int16 lhs,const i3ds& rhs){return rhs*lhs;}
    
	inline bool operator !=(i3ds const& a) {return v[0]!=a.v[0]||v[1]!=a.v[1]||v[2]!=a.v[2];}
    inline bool operator ==(i3ds const& a) {return !(*this!=a);}
    
    friend uint64 prod(i3ds i) {return 1UL*i.x*i.y*i.z;}

  static void setOutMode(int kmode){outMode=kmode;}
  static int  getOutMode(){return outMode;}
	
};

ostream& operator<<(ostream &os,const i3ds &v);

inline i3ds min(const i3ds &a, const i3ds &b){
  i3ds x;
  for (int j=0;j<3;j++) x[j]=a[j]<b[j]? a[j]: b[j];
  return x;
}

inline i3ds max(const i3ds &a, const i3ds &b){
  i3ds x;
  for (int j=0;j<3;j++) x[j]=a[j]>b[j]? a[j]: b[j];
  return x;
}

inline bool operator==(const i3ds &a,const i3ds &b){
    if(a.x != b.x || a.y != b.y || a.z != b.z ) return false;
	return true;
}

inline bool operator!=(const i3ds &a,const i3ds &b){
    if(a.x != b.x || a.y != b.y || a.z != b.z ) return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
class ui3ds {
private:  
  static int32 outMode;
public:
    union {
		uint16 v[3];
		struct { uint16 x,y,z;}; // CART
		struct { uint16 Alpha,R,Z;}; // CYL
		struct { uint16 phi,theta,rho;}; // SPHER
	};

	ui3ds(){};

	ui3ds(uint16 a,uint16 b, uint16 c){
 		v[0]=a;v[1]=b;v[2]=c;
 	};
	
	ui3ds(const ui3ds &a){
 		v[0]=a.v[0];v[1]=a.v[1];v[2]=a.v[2];
 	};
	
	ui3ds(uint16 n){
 		v[0]=v[1]=v[2]=n;
 	}
	
 	ui3ds(uint16 *a){
 		v[0]=a[0];v[1]=a[1];v[2]=a[2];
 	}
 	
	inline uint16 & operator[](int j){return v[j];};
	inline const uint16 & operator[](int j) const {return v[j];};
	
	inline void set(uint16 x,uint16 y,uint16 z){v[0]=x;v[1]=y;v[2]=z;}
	
	
	inline ui3ds& operator=(const ui3ds &a){v[0]=a.v[0];v[1]=a.v[1];v[2]=a.v[2];return *this;}
    
    inline ui3ds& operator+=(uint16 a){v[0]+=a;v[1]+=a;v[2]+=a;return *this;};
    inline ui3ds& operator-=(uint16 a){v[0]-=a;v[1]-=a;v[2]-=a;return *this;};
    inline ui3ds& operator*=(uint16 a){v[0]*=a;v[1]*=a;v[2]*=a;return *this;};
    inline ui3ds& operator/=(uint16 a){v[0]/=a;v[1]/=a;v[2]/=a;return *this;};
    inline ui3ds&  operator=(uint16 a){v[0] =a;v[1] =a;v[2] =a;return *this;};
    
    inline ui3ds& operator++(){v[0]++;v[1]++;v[2]++;return *this;};
    inline ui3ds& operator--(){v[0]--;v[1]--;v[2]--;return *this;};
    inline ui3ds& operator++(int ){v[0]++;v[1]++;v[2]++;return *this;};
    inline ui3ds& operator--(int ){v[0]--;v[1]--;v[2]--;return *this;};
    
    inline ui3ds& operator+=(const ui3ds &a){v[0]+=a.v[0];v[1]+=a.v[1];v[2]+=a.v[2];return *this;}
    inline ui3ds& operator-=(const ui3ds &a){v[0]-=a.v[0];v[1]-=a.v[1];v[2]-=a.v[2];return *this;}
    inline ui3ds& operator*=(const ui3ds &a){v[0]*=a.v[0];v[1]*=a.v[1];v[2]*=a.v[2];return *this;}
    inline ui3ds& operator/=(const ui3ds &a){v[0]/=a.v[0];v[1]/=a.v[1];v[2]/=a.v[2];return *this;}
    
    inline ui3ds operator-() const {return ui3ds(-v[0],-v[1],-v[2]);}
    
    friend ui3ds operator+(ui3ds lhs,       // passing first arg by value helps optimize chained a+b+c
                    const ui3ds& rhs) // alternatively, both parameters may be const references.
  	{
    	 return lhs += rhs; // reuse compound assignment and return the result by value
  	}
    friend ui3ds operator-(ui3ds lhs,const ui3ds& rhs){return lhs -= rhs;}
    friend ui3ds operator*(ui3ds lhs,const ui3ds& rhs){return lhs *= rhs;}
    friend ui3ds operator/(ui3ds lhs,const ui3ds& rhs){return lhs /= rhs;}
    
	friend ui3ds operator*(uint16 lhs,const ui3ds& rhs){return rhs*lhs;}

	inline bool operator !=(ui3ds const& a) {return v[0]!=a.v[0]||v[1]!=a.v[1]||v[2]!=a.v[2];}
    inline bool operator ==(ui3ds const& a) {return !(*this!=a);}
    
    friend uint64 prod(ui3ds i) {return 1UL*i.x*i.y*i.z;}

    static void setOutMode(int kmode){outMode=kmode;}
    static int  getOutMode(){return outMode;}
	
};

ostream& operator<<(ostream &os,const ui3ds &v);

inline ui3ds min(const ui3ds &a, const ui3ds &b){
  ui3ds x;
  for (int j=0;j<3;j++) x[j]=a[j]<b[j]? a[j]: b[j];
  return x;
}

inline ui3ds max(const ui3ds &a, const ui3ds &b){
  ui3ds x;
  for (int j=0;j<3;j++) x[j]=a[j]>b[j]? a[j]: b[j];
  return x;
}

inline bool operator==(const ui3ds &a,const ui3ds &b){
    if(a.x != b.x || a.y != b.y || a.z != b.z ) return false;
	return true;
}

inline bool operator!=(const ui3ds &a,const ui3ds &b){
    if(a.x != b.x || a.y != b.y || a.z != b.z ) return true;
	return false;
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
class islice3d {
public:
	i3d imin,imax;
	islice3d(): imin(0),imax(0){}	
	islice3d(i3d a,i3d b): imin(a),imax(b){}	
	i3d dims() {return imax-imin+1;}
	uint64 size() {return prod(dims());}
};

} // namespace

#endif // IDX3D_H
