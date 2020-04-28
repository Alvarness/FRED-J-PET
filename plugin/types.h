/*
 *  types.h
 *  Fred
 *
 *  Created by Angelo Schiavi on 18/01/08.
 *  Copyright 2008 A. Schiavi. All rights reserved.
 *
 */

#ifndef TYPES_H
#define TYPES_H

#include <typeinfo>
#include <cstring>

namespace fred {

typedef signed char           int8;       // 8 bit signed
typedef unsigned char         uint8;      // 8 bit unsigned
typedef short                 int16;      // 16 bit signed
typedef unsigned short        uint16;     // 16 bit unsigned
typedef int                   int32;      // 32 bit signed
typedef unsigned int          uint32;     // 32 bit unsigned
typedef long long             int64;      // 64 bit signed
typedef unsigned long long    uint64;     // 64 bit unsigned

typedef uint32 uint;
typedef uint64 ulong;

typedef float                 fp32;     // 32 bit floating point 
typedef double                fp64;     // 64 bit floating point

typedef float                 real32;     // 32 bit floating point 
typedef double                real64;     // 64 bit floating point

typedef float                 float32;     // 32 bit floating point 
typedef double                float64;     // 64 bit floating point

typedef float                 flt32;     // 32 bit floating point 
typedef double                flt64;     // 64 bit floating point


// typedef int32 rayType;



template<typename T>
void type2str(T v,char s[9]){
	if (typeid(v)==typeid(int8     )) {strncpy(s,"1int8   \0",8);return;}
	if (typeid(v)==typeid(uint8    )) {strncpy(s,"2uint8  \0",8);return;}
	if (typeid(v)==typeid(int16    )) {strncpy(s,"3int16  \0",8);return;}
	if (typeid(v)==typeid(uint16   )) {strncpy(s,"4uint16 \0",8);return;}
	if (typeid(v)==typeid(int32    )) {strncpy(s,"5int32  \0",8);return;}
	if (typeid(v)==typeid(uint32   )) {strncpy(s,"6uint32 \0",8);return;}
	if (typeid(v)==typeid(int64    )) {strncpy(s,"7int64  \0",8);return;}
	if (typeid(v)==typeid(uint64   )) {strncpy(s,"8uint64 \0",8);return;}
	if (typeid(v)==typeid(real32   )) {strncpy(s,"9real32 \0",8);return;}
	if (typeid(v)==typeid(real64   )) {strncpy(s,"Areal64 \0",8);return;}
	strncpy(s,"0unknown\0",8);
}

} // namespace

#endif // TYPES_H
