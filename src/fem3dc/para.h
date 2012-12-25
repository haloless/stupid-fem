#pragma once

#include <cstdio>

#include <mkl.h>
//#include <armadillo>

// blitz
#define BZ_HAVE_NAMESPACES
#include <blitz/array.h>
namespace bz = blitz;
typedef bz::Array<double, 1> FDArray1;
typedef bz::Array<double, 2> FDArray2;
typedef bz::Array<double, 3> FDArray3;
typedef bz::Array<MKL_INT, 1> FIArray1;
typedef bz::Array<MKL_INT, 2> FIArray2;

#define FORT_ARR(name, ...) name(__VA_ARGS__, bz::fortranArray)

template<typename T, int rank>
inline void zeroFortArray(bz::Array<T, rank> &a) {
	a = 0;
}


#define STRINGLIZE0(sth) #sth
#define STRINGLIZE(sth) STRINGLIZE0(sth)

/**
 * implicit integer (i-n)
 * implicit double precision (a-h,o-z)
 */

/**
c----------------------------------------------------------------------
c  IELM  : Max. number of an element 
c  INODE : Max. number of a node
c  IRD   : Max. number of integration point in a element
c  IGK   : Max. size of stiffness matrix  INODE*IDIM
c  IBN   : Max. size of half band width
c  IBOUND: Max. number of constrain (force and displacement)
c  INA   : Max. number of integration point  IELM*IRD
c  INMTE : Max. number of a node in a element --> 8
c  INMNM : Max. size of element stiffness matrix  INMTE*IDIM
c  IMNUM : Max. number of the kind of material 
c  IDM   : DIMENSION
c  ICOM  : number of stress component (2-Dim->4 ,3Dim->6)
c----------------------------------------------------------------------
ccc 3D-8node-2*2*2 integration
c      PARAMETER (IELM=1000,INODE=1000,IRD=8,IGK=3000,IBN=100,IBOUND=100,
c     +     INA=8000,INMTE=8,INMNM=24,IMNUM=10)
      PARAMETER (IELM=1000,INODE=3000,IRD=8,IGK=9000,IBN=300
     +     ,IBOUND=3000,INA=8000,INMTE=8,INMNM=24,IMNUM=10)
cccc 3D-20node-3*3*3 integration
c      PARAMETER (IELM=1000,INODE=3000,IRD=27,IGK=9000,IBN=300,
c     +     IBOUND=100,INA=27000,INMTE=20,INMNM=60,IMNUM=10)
c      PARAMETER (IELM=1000,INODE=3000,IRD=20,IGK=9000,IBN=300,
c     +     IBOUND=100,INA=27000,INMTE=20,INMNM=60,IMNUM=10)
      PARAMETER (IDIM=3,ICOM=6)
*/


const MKL_INT IDIM = 3;						// dimension
const MKL_INT ICOM = IDIM==3 ? 6 : 4;		// stress component

const MKL_INT IELM = 1000;					// max number of element
const MKL_INT INODE = 3000;					// max number of nodes
const MKL_INT IRD = 8;						// number of integration points in an element
const MKL_INT IGK = INODE * IDIM;			// stiffness matrix
const MKL_INT IBN = 300;					// max size of half band width
const MKL_INT IBOUND = 3000;				// max number of constraints (force/displacement)
const MKL_INT INA = IELM * IRD;				// max number of quadrature point
const MKL_INT INMTE = 8;					// number of node in an element
const MKL_INT INMNM = INMTE * IDIM;			// max size of element stiffness matrix
const MKL_INT IMNUM = 10;					// number of materials

#define M_MaxNumOfElems IELM
#define M_MaxNumOfNodes INODE
#define M_MaxNumOfQuadPointsInElem IRD
#define M_MaxSizeOfGlobalStiffMatrix IGK
#define M_MaxSizeOfHalfBandWidth IBN
#define M_MaxNumOfConstraints IBOUND
#define M_MaxNumOfQuadPoints INA
#define M_MaxNumOfNodesInElem INMTE
#define M_MaxSizeOfElementStiffMatrix INMNM
#define M_MaxNumOfMaterials IMNUM


#ifdef INTERN_SHARED_VARS
#define EXTERN
#else
#define EXTERN extern
#endif

/**
cc constant
c--------------------------------------------------------------
c FOR SUBROUTIE DATA
c  NGAUSS : Order of GAUSS integration
c  NPARA  : Analysis condition (0:plain stress,1:plain strain)
cccc  T      : Thickness
c  NTN    : total number of nodes
c  NTE    : total number of elements
c  NTF    : total number of freedoms
cccc  IE     : element number(variables)
c  NTFOR  : total number of freedoms with initial load
c  NTDIS  : total number of freedoms with initial displacement
c  NNA    : The total number of integration points
c  MNUM   : number of the kind of material
c  NB     : half band width
c--------------------------------------------------------------
	COMMON /CON1/ NGAUSS,NPARA,NTN,NTE,NTF,NNA,MNUM
	COMMON /CON2/ NTFOR,NTDIS,NB
	COMMON /CON3/ IE
*/
// CON1
EXTERN MKL_INT NGAUSS;
EXTERN MKL_INT NPARA;
EXTERN MKL_INT NTN;
EXTERN MKL_INT NTE;
EXTERN MKL_INT NTF;
EXTERN MKL_INT NNA;
EXTERN MKL_INT MNUM;
// CON2
EXTERN MKL_INT NTFOR;
EXTERN MKL_INT NTDIS;
EXTERN MKL_INT NB;
// CON3
EXTERN MKL_INT IE;

#define M_OrderOfGaussQuad NGAUSS
#define M_AnalysisCond NPARA
#define M_TotalNumOfNodes NTN
#define M_TotalNumOfElems NTE
#define M_TotalNumOfFreedoms NTF
#define M_TotalNumOfQuadPoints NNA
#define M_TotalNumOfMaterials MNUM

#define M_TotalNumOfFreedomsWithInitLoad NTFOR
#define M_TotalNumOfFreedomsWithInitDisp NTDIS
#define M_HalfBandWith NB

//#define M_NumOfElem IE

/**
c
cccc FOR EQUATION
c  NAEQ : NUMBER OF EQUATION
c     NEQB : ROOT FREEDOM
c     NTEQ : NUMBER OF COPY FREEDOM
c     NEQ  : COPY FREEDOM
	COMMON /EQ0/ NAEQ
	COMMON /EQ1/ NEQB(5),NTEQ(5)
	COMMON /EQ2/ NEQ(5,100
*/
// EQ0
EXTERN MKL_INT NAEQ;
// EQ1, EQ2
#ifdef INTERN_SHARED_VARS
//EXTERN FIArray1 NEQB(5, bz::fortranArray);
//EXTERN FIArray1 NTEQ(5, bz::fortranArray);
//EXTERN FIArray2 NEQ(5, 100, bz::fortranArray);
FIArray1 FORT_ARR(NEQB, 5);
FIArray1 FORT_ARR(NTEQ, 5);
FIArray2 FORT_ARR(NEQ, 5,100);
#else
//EXTERN MKL_INT NEQB[5];
//EXTERN MKL_INT NTEQ[5];
//EXTERN MKL_INT NEQ[5][100];
EXTERN FIArray1 NEQB;
EXTERN FIArray1 NTEQ;
EXTERN FIArray2 NEQ;
#endif

#define M_TotalNumOfEqns NAEQ
#define M_RootFreedoms NEQB
#define M_NumOfCopyFreedoms NTEQ
#define M_CopyFreedoms NEQ

// file pointers


//
#define ELEMTYPE_FINITE 0
#define ELEMTYPE_INFINITE 1

#define MATERIALCHAR_ISOTROPIC 1
#define MATERIALCHAR_ANISOTROPIC 2

