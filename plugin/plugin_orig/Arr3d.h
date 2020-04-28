#ifndef ARR3D_H
#define ARR3D_H

#include "idx3d.h"
#include "types.h"
#include <iostream>
#include <iomanip>
using namespace std;

namespace fred {

template <class T>
class Arr3d {
public:
	T *data;
	uint64 N;
	
	union { //size vector
		uint32 sizev[3];
		struct { uint32 nx,ny,nz;};
	};
	uint32 stridev[3]; // stride vector
	//basev[3]; // base vector

public:
	Arr3d();
	Arr3d(int n, int m, int k);
	Arr3d(i3d dims);

	
	inline T & operator[](const int i) {return data[i];} // access via voxel index
	inline const T & operator[](const int i) const {return data[i];}

	inline T& el(int i,int j,int k){return data[i+j*stridev[1]+k*stridev[2]];} // access via 3d idx ASSUMING COLUMN MAJOR
	inline T& el(size_t i){return data[i];}
	
	inline T& el(i3d idx){return data[idx.x+idx.y*stridev[1]+idx.z*stridev[2]];} // access via 3d idx using idx3d ASSUMING COLUMN MAJOR
	inline T& el(i3ds idx){return data[(int)idx.x+(int)idx.y*stridev[1]+(int)idx.z*stridev[2]];} 
	inline T& el(ui3ds idx){return data[(int)idx.x+(int)idx.y*stridev[1]+(int)idx.z*stridev[2]];} 
	
	inline T& operator()(int i,int j,int k){return data[i+j*stridev[1]+k*stridev[2]];}
	inline T  operator()(int i,int j,int k) const {return data[i+j*stridev[1]+k*stridev[2]];}

    inline T& operator()(i3d idx){return data[(int)idx.x+(int)idx.y*stridev[1]+(int)idx.z*stridev[2]];}
    inline T  operator()(i3d idx)const {return data[(int)idx.x+(int)idx.y*stridev[1]+(int)idx.z*stridev[2]];}

    inline T& operator()(i3ds idx){return data[(int)idx.x+(int)idx.y*stridev[1]+(int)idx.z*stridev[2]];}
    inline T  operator()(i3ds idx)const {return data[(int)idx.x+(int)idx.y*stridev[1]+(int)idx.z*stridev[2]];}

    inline T& operator()(ui3ds idx){return data[(int)idx.x+(int)idx.y*stridev[1]+(int)idx.z*stridev[2]];}
    inline T  operator()(ui3ds idx)const {return data[(int)idx.x+(int)idx.y*stridev[1]+(int)idx.z*stridev[2]];}


	inline uint32 vxlidx(int i,int j,int k) {return i+j*stridev[1]+k*stridev[2];}
	inline uint32 vxlidx(i3d idx) {return idx.x+idx.y*stridev[1]+idx.z*stridev[2];}
	
	inline i3d vxlcoord(int idx);
	

	inline size_t numel(){return N;}
	
	inline int dim1() const {return sizev[0];}
	
	inline int dim2() const {return sizev[1];}
	
	inline int dim3() const {return sizev[2];}
		
	i3d dims() const {return i3d(sizev[0],sizev[1],sizev[2]);};
	
	~Arr3d();

	Arr3d<T>&	operator=(T v);
	Arr3d<T>&	operator=(Arr3d<T>& B);
	
	Arr3d<T>&	operator+=(const T v);
	Arr3d<T>&	operator-=(const T v);
	Arr3d<T>&	operator*=(const T v);
	Arr3d<T>&	operator/=(const T v);

	void resize(i3d nn);
	void resize(int n, int m, int k);


	void reshape(int n, int m, int k);
	void realloc();
	
	void range(T &vmin,T &vmax);
	
	T sum();
	
	T max();
	T min();
	T max(i3d *loc);
	T min(i3d *loc);

	void clamp(T vmin, T vmax);

	void  indexes(int dim =-1);
	
	void info(ostream &os=std::cout);
};

template <class T>
Arr3d<T>::Arr3d() {sizev[0]=sizev[1]=sizev[2]=0; data=NULL; N=0;}

template <class T>
Arr3d<T>::Arr3d(int n, int m, int k)
{
	sizev[0]=sizev[1]=sizev[2]=0; data=NULL; N=0;
	resize(n,m,k);
}	

template <class T>
Arr3d<T>::Arr3d(i3d dims)
{
	sizev[0]=sizev[1]=sizev[2]=0; data=NULL; N=0;
	resize(dims[0],dims[1],dims[2]);
}	

template <class T>
void Arr3d<T>::reshape(int n, int m, int k)
{
	N = 1UL*n*m*k;
	sizev[0]=n;
	sizev[1]=m;
	sizev[2]=k;
	stridev[0]=1; // column major
	stridev[1]=n;
	stridev[2]=n*m; 
}

template <class T>
void Arr3d<T>::realloc()
{
	if (data != NULL) { //release memory
		delete[] (data);
	}
	data = new T[N];
}

template <class T>
void Arr3d<T>::resize(int n, int m, int k)
{
	reshape(n,m,k);
	realloc();	
}

template <class T>
void Arr3d<T>::resize(i3d nn)
{ 
	resize(nn[0],nn[1],nn[2]);
}


template <class T>
inline Arr3d<T>&  Arr3d<T>::operator=(const T v)
{
	T *p=data;
	for(size_t i=0; i<N; i++,p++) *p=v;
	return *this;
}

template <class T>
inline Arr3d<T>&  Arr3d<T>::operator+=(const T v)
{
	T *p=data;
	for(size_t i=0; i<N; i++,p++) *p+=v;
	return *this;
}

template <class T>
inline Arr3d<T>&  Arr3d<T>::operator-=(const T v)
{
	T *p=data;
	for(size_t i=0; i<N; i++,p++) *p-=v;
	return *this;
}

template <class T>
inline Arr3d<T>&  Arr3d<T>::operator*=(const T v)
{
	T *p=data;
	for(size_t i=0; i<N; i++,p++) *p*=v;
	return *this;
}

template <class T>
inline Arr3d<T>&  Arr3d<T>::operator/=(const T v)
{
	T *p=data;
	for(size_t i=0; i<N; i++,p++) *p/=v;
	return *this;
}


template <class T> 
Arr3d<T> & Arr3d<T>::operator=(Arr3d<T> &B)
{
	if (this != &B) {
	if(dims() != B.dims()) resize(B.dim1(),B.dim2(),B.dim3());

		memcpy((void *)data,(void *) B.data,1UL*B.numel()*sizeof(T));
	}
	return *this;
}

template <class T>
Arr3d<T>::~Arr3d()
{
	if (data != 0) {
		delete[] (data);
	}
}

template <class T>
void Arr3d<T>::range(T &vmin,T &vmax)
{
	T *p=data;
	vmin=vmax=p[0];
	for(size_t i=1; i<N; i++){
		if(p[i]<vmin) vmin=p[i];
		if(p[i]>vmax) vmax=p[i];
	}
}

template <class T>
T Arr3d<T>::sum()
{
	T *p=data;
	T s0=0;
	for(size_t i=0; i<N; i++) s0+=p[i];
	return s0;
}


template <class T>
T Arr3d<T>::max()
{
	T *p=data;
	T v=*p;
	for(;p<data+N; p++) if(v<*p) v=*p;
	return v;
}

template <class T>
T Arr3d<T>::min()
{
	T *p=data;
	T v=*p;
	for(;p<data+N; p++) if(v>*p) v=*p;
	return v;
}


template <class T>
T Arr3d<T>::max(i3d *loc)
{
	T *p=data;
	T v=*p;
	size_t ivxl=0;
	for(size_t i=0; i<N; i++,p++) if(v<*p) {v=*p;ivxl=i;}
	*loc = vxlcoord(ivxl);
	return v;
}

template <class T>
T Arr3d<T>::min(i3d *loc)
{
	T *p=data;
	T v=*p;
	size_t ivxl=0;
	for(size_t i=0; i<N; i++,p++) if(v>*p) {v=*p;ivxl=i;}
	*loc = vxlcoord(ivxl);
	return v;
}

template <class T>
void Arr3d<T>::clamp(T vmin, T vmax)
{
	size_t n=1UL*nx*ny*nz;
	T *p=data;
	//#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
	for(size_t i=0; i<n; i++) p[i] = (((p[i]) > (vmax)) ? (vmax) : (((p[i]) < (vmin)) ? (vmin) : (p[i])));
}



template <class T>
i3d Arr3d<T>::vxlcoord(int idx){
	i3d p;
	p[2] = idx/stridev[2];
	int rem = idx%stridev[2];
	p[1] = rem/stridev[1];
	p[0] = rem%stridev[1];
	return p;
}

template <class T> void	 Arr3d<T>::indexes(int dim){
	Arr3d<T> &A=*this;
	for(int i=0;i<nx;i++) {
		for(int j=0;j<ny;j++) {
			for(int k=0;k<nz;k++) {
				switch(dim){
					default:
					case -1:
						A(i,j,k)=100*(i+1)+10*(j+1)+k+1;
						break;
					case 0:
						A(i,j,k)=i+1;
						break;	
					case 1:
						A(i,j,k)=j+1;
						break;	
					case 2:
						A(i,j,k)=k+1;
						break;	
				}	
			}}}
}

template <class T> ostream& operator<<(ostream& os,Arr3d<T>& A){
	ios_base::fmtflags	oldFlags= os.flags();
	int prec = 3;
	long savedPrec =	os.precision();
	os.precision(prec);
	os<<"array dims=("<<A.nx<<','<<A.ny<<','<<A.nz<<")"<<endl;
	for(int k=0;k<A.nz;k++) {
		for(int i=0;i<A.nx;i++) {
			for(int j=0;j<A.ny;j++) {
					os << setw(prec) << A(i,j,k)<<' ';
			}
		os<<endl;	
		}
		os<<endl;	
	}
	os<<endl;
	os.precision(savedPrec);
	os.flags(oldFlags);
	return os;
}

template <class T> void	 Arr3d<T>:: info(ostream &os){
	os<<"\tdims: "<<this->dims()<<" => "<<this->N<<" elements"<<endl;
	os<<"\tdatatype: "<<typeid(T).name()<<endl;
	os<<"\tmemory occupation: ";
	size_t B = this->N*sizeof(T);
	os<<B<<" B ";
	if (B<1024L*1024L) { os<<" = "<<B/1024<<" KB";}
	else if (B<1024L*1024L*1024L) { os<<" = "<<B/1024/1024<<" MB";}
	else { os<<" = "<<B/1024/1024/1024<<" GB";}
	os<<endl;
}

} // namespace

#endif // ARR3D_H