
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <algorithm>

#define INTERN_SHARED_VARS
#include "para.h"
#include "const.h"
#pragma comment(lib, "blitz.lib")

#define VTKHELPER_LINKER_PRAGMA
#include <vtkHelper.hpp>
#include <vtkHexahedron.h>
#include <vtkVoxel.h>
#include <vtkCellData.h>


static FILE *fpDatacheck = NULL;
static FILE *fpResult = NULL;

static void fatalError(const char *msg, int exitcode=1) {
	fprintf(stdout, "ERROR: %s\n", msg);
	exit(exitcode);
}

/**
 * set J, Jinv, detJ
 */
void JMat(const FDArray1 &x, const FDArray1 &y, const FDArray1 &z,
	const FDArray1 &dNdXi, const FDArray1 &dNdEta, const FDArray1 &dNdZeta,
	const FDArray1 &dNdXiInfty, const FDArray1 &dNdEtaInfty, const FDArray1 &dNdZetaInfty,
	FDArray2 &J, FDArray2 &Jinv, double &detJ,
	const FIArray2 &nodeTable, const FIArray1 &numOfNodes, const FIArray1 &elemType);

void data(const char* filename, FILE *fpLog,
	FIArray1 &numOfNodes, FDArray2 &elasticConstants, FIArray2 &nodesTable,
	FDArray1 &x, FDArray1 &y, FDArray1 &z,
	FIArray1 &freedomOfLoadPoint, FDArray1 &loadOfLoadPoint, 
	FIArray1 &freedomOfDispPoint, FDArray1 &dispOfDispPoint,
	FIArray1 &numOfQuadPoints, 
	FIArray1 &materialOfElem, FIArray1 &characterOfMaterial,
	FIArray1 &userElem, FIArray1 &userNode,
	FIArray1 &elemKind,
	int &nstep) 
{
	FILE *fpData = fopen(filename, "r");
	if(!fpData) {
		char msg[512];
		sprintf_s(msg, "Failed to open data file %s.", filename);
		fatalError(msg);
	}

	const int bufsize = 1024;
	char buf[bufsize];

#define DATACHECK(fmt, ...) fprintf_s(fpLog, fmt, __VA_ARGS__)
//#define DATASCAN(fmt, ...) fscanf_s(fpData, fmt, __VA_ARGS__)
#define DATASCAN(fmt, ...) while(fgets(buf,bufsize,fpData) && (buf[0]=='#')); \
	sscanf_s(buf, fmt, __VA_ARGS__)
	
	// line 1
	DATASCAN("%d", &nstep);			// increments
	DATASCAN("%d", &NGAUSS);		// order of quadrature
	DATASCAN("%d", &MNUM);			// number of materials
	//
	DATACHECK("Data Check: %s\n", filename);
	DATACHECK("Total # of increment %d\n", nstep);
	DATACHECK("NUMBER OF MATERIAL = %d\n", MNUM);
	DATACHECK("INTEGRATION ORDER = %d\n", NGAUSS);

	// line 4, materials 
	for(int iMatter=1; iMatter<=MNUM; iMatter++) {
		int materialChar = 1;
		DATASCAN("%d", &materialChar);

		characterOfMaterial(iMatter) = materialChar;
		if(materialChar == MATERIALCHAR_ISOTROPIC) {
			double young, poisson, yield;
			DATASCAN("%*s %lf %lf %lf", &young, &poisson, &yield);
			elasticConstants(1, iMatter) = young;
			elasticConstants(2, iMatter) = poisson;
			elasticConstants(3, iMatter) = yield;

			DATACHECK("Material %d\n", iMatter);
			DATACHECK("YOUNG=%lf, POISSON=%lf, Yield stress=%lf\n", young, poisson, yield);
		} else {
			// only isotropic material supported
			char msg[512];
			sprintf_s(msg, "Material %d is not ISOTROPIC, not supported.", iMatter);
			fatalError(msg);
		}
	}

	// line 6, meshes
	DATASCAN("%d %d", &NTN, &NTE);	// nodes and elements
	NTF = NTN * IDIM;
	//
	DATACHECK("Total Nodes = %d\n", NTN);
	DATACHECK("Total Elements = %d\n", NTE);
	DATACHECK("Total Freedoms = %d\n", NTF);

	// line 7, elements
	NNA = 0;						// number of quad. points
	for(int iElem=1; iElem<=NTE; iElem++) {
		int ii, iMatterType, iNodeNum, iElemKind;
		DATASCAN("%d %d %d %d", &ii, &iMatterType, &iNodeNum, &iElemKind);
		materialOfElem(iElem) = iMatterType;
		numOfNodes(iElem) = iNodeNum;
		elemKind(iElem) = iElemKind;
		// iElem: program element number; NEUS(IE): user element
		userElem(iElem) = ii;

		switch(iNodeNum) {
		case 8:			// 8nodes 8-men-tai?
			DATASCAN("%d %d %d %d %d %d %d %d", 
				&nodesTable(1,iElem),&nodesTable(2,iElem),&nodesTable(3,iElem),&nodesTable(4,iElem),
				&nodesTable(5,iElem),&nodesTable(6,iElem),&nodesTable(7,iElem),&nodesTable(8,iElem));
			// NGauss = 2 or 3
			numOfQuadPoints(iElem) = NGAUSS*NGAUSS*NGAUSS;
			//
			DATACHECK("ELEMENT %d: MATERIAL %d : 8 NODES 8-mentai\n", iElem, iMatterType);
			break;
		case 6:			// 6nodes sankaku-tyu
		case 15:		// 15nodes sankaku-tyu
		case 20:		// 20nodes 8-men-tai
		default:
			fatalError("Element not supported."); 
			break;
		}

		// 
		NNA += numOfQuadPoints(iElem);
	}
	//
	DATACHECK("TOTAL INT.POINTS = %d\n", NNA);
	DATACHECK("COORDINATE OF NODES\n");
	DATACHECK("PRO.NUM. , USR NUM. , COORDINATION OF NODE\n");

	// line 79, nodes
	for(int iNode=1; iNode<=NTN; iNode++) {
		int ii; 
		DATASCAN("%d %lf %lf %lf", &ii, &x(iNode), &y(iNode), &z(iNode));
		// iNode: program node number; NNUS(iNode): user node number
		userNode(iNode) = ii;

		DATACHECK("%d %d %le %le %le\n", iNode, ii, x(iNode), y(iNode), z(iNode));
	}

	// reconstruct nodes, not necessary?
	DATACHECK("*ELEMENT NODE IS RECONSTRUCTED\n");
	for(int iElem=1; iElem<=NTE; iElem++) {
		const int iNodesNum = numOfNodes(iElem);
		
		DATACHECK("USER    ");
		for(int i=1; i<=iNodesNum; i++) {
			DATACHECK("%4d ", nodesTable(i, iElem));
		}
		DATACHECK("\n");

		bool iNodesAllMatched = true;
		for(int j2=1; j2<=iNodesNum; j2++) {
			int nnj2 = nodesTable(j2, iElem);

			bool hasMatch = false;
			for(int i2=1; i2<=NTN; i2++) {
				int kk = userNode(i2);
				if(kk == nnj2) {
					nodesTable(j2, iElem) = i2;
					hasMatch = true;
					break;
				}
			}
			iNodesAllMatched = iNodesAllMatched && hasMatch;
		}

		if(iNodesAllMatched) {
			DATACHECK("PROGRAM ");
			for(int i=1; i<=iNodesNum; i++) {
				DATACHECK("%4d ", nodesTable(i, iElem));
			}
			DATACHECK("\n");
		} else {
			char msg[512]; sprintf_s(msg, "NO MATCH FOR element %d.", iElem);
			fatalError(msg);
		}
	}

	// line 179
	DATACHECK("*FORCE BOUNDARY\n");
	DATASCAN("%d", &NTFOR);			// number of freedoms with initial load
	if(NTFOR > 0) {
		for(int i=1; i<=NTFOR; i++) {
			int ifo, ixy, ifop=-1;
			double ff;
			DATASCAN("%d %d %lf", &ifo, &ixy, &ff); // node, direction, force
			//c  IFO(USER NODE NUMBER) --> IFOP(PROGRAM NODE NUMBER)
			for(int i2=1; i2<=NTN; i2++) {
				int kk = userNode(i2);
				if(kk == ifo) {
					ifop = i2;
					break;
				}
			}
			if(ifop == -1) { // not found node
				DATACHECK("INVALID LOAD BC!(NO NODE)\n");
				fatalError("Failed to set load BC.");
			}
			if(ixy > IDIM) {
				DATACHECK("INVALID LOAD BC!(IXY > IDIM)\n");
				fatalError("Failed to set load BC.");
			}

			int nFree = (ifop-1)*IDIM + ixy;
			freedomOfLoadPoint(i) = nFree;
			loadOfLoadPoint(i) = ff;
			//
			if(ixy == 1) DATACHECK("NODE %d -->FX= %lf %d\n", ifo, loadOfLoadPoint(i), nFree);
			if(ixy == 2) DATACHECK("NODE %d -->FY= %lf %d\n", ifo, loadOfLoadPoint(i), nFree);
			if(ixy == 3) DATACHECK("NODE %d -->FZ= %lf %d\n", ifo, loadOfLoadPoint(i), nFree);
		}
	}

	// line 182
	DATACHECK("*DISPLACEMENT BOUNDARY\n");
	DATASCAN("%d", &NTDIS);			// number of freedoms with displacement
	if(NTDIS > 0) {
		for(int i=1; i<=NTDIS; i++) {
			int idi, ixy, idip=-1;
			double uu;
			DATASCAN("%d %d %lf", &idi, &ixy, &uu); // node, direction, displacement
			if(fabs(uu) > 1e-9) {
				fatalError("Constrained displacement BC not supported.");
			}

			for(int i2=1; i2<=NTN; i2++) {
				int kk = userNode(i2);
				if(kk == idi) {
					idip = i2;
					break;
				}
			}
			if(idip == -1) {
				DATACHECK("INVALID BOUNDARY CONDITION!(NO NODE)");
				fatalError("Failed to set disp BC.");
			}
			if(ixy > IDIM) {
				DATACHECK("INVALID BOUNDARY CONDITION!(ixy > IDIM)");
				fatalError("Failed to set disp BC.");
			}

			int nFree = (idip-1)*IDIM + ixy;
			freedomOfDispPoint(i) = nFree;
			dispOfDispPoint(i) = uu;
			//
			if(ixy == 1) DATACHECK("NODE %d -->UX= %lf FREEDOM= %d\n", idi, dispOfDispPoint(i), nFree);
			if(ixy == 2) DATACHECK("NODE %d -->UY= %lf FREEDOM= %d\n", idi, dispOfDispPoint(i), nFree);
			if(ixy == 3) DATACHECK("NODE %d -->UZ= %lf FREEDOM= %d\n", idi, dispOfDispPoint(i), nFree);
		}
	}

	// line 213, equation
	DATASCAN("%d", &NAEQ);
	if(NAEQ != 0) {
		fatalError("EQUATION is not supported.");
	}

	fclose(fpData);
#undef DATACHECK
#undef DATASCAN
}

/**
 * detect the upper band width for compact band matrix storage
 */
void band(const FIArray2 &nodesTable, const FIArray1 &numberOfNodesInElem,
	FILE *fpLog) 
{
	if(NAEQ != 0) {
		fatalError("EQUATION not supported.");
		// TODO equations
	} else {
		int inb = 0;
		
		for(int iElem=1; iElem<=NTE; iElem++) {
			int iNodesNum = numberOfNodesInElem(iElem);
			int iMax = nodesTable(iNodesNum, iElem);
			int iMin = nodesTable(iNodesNum, iElem);
			
			for(int j=1; j<=iNodesNum; j++) {
				int nodeId = nodesTable(j, iElem);
				iMax = std::max(iMax, nodeId);
				iMin = std::min(iMin, nodeId);
			}

			int iDiff = iMax - iMin;
			if(iDiff > inb) inb = iDiff;
		}

		// set global band half width
		NB = IDIM * (inb+1);
	}

	if(fpLog) fprintf(fpLog, "BAND HALF WIDTH= %d\n", NB);
}


void check() {
	if(NTN > INODE) {
		fatalError("Number of node exceeds limitation!");
	}
	
	if(NTE > IELM) {
		fatalError("Number of element exceeds limitation!");
	}

	if(NB > IBN) {
		fatalError("Band width exceeds limitation!");
	}

	if(INODE*IDIM != IGK) {
		fatalError("INODE*IDIM =/= IGK");
	}

	if(IELM*IRD != INA) {
		fatalError("INA =/= IELM*IRD");
	}
}

void init(FDArray3 &stressAtElemForPoint, FDArray3 &strainAtElemForPoint) {
	//zeroFortArray(stressAtElemForPoint);
	//zeroFortArray(strainAtElemForPoint);
	for(int i=1; i<=ICOM; i++) {
		for(int j=1; j<=IELM; j++) {
			for(int k=1; k<=IRD; k++) {
				stressAtElemForPoint(i, j, k) = 0;
				strainAtElemForPoint(i, j, k) = 0;
			}
		}
	}
}

void clear(FDArray2 &elemStiffKMat, FDArray2 &globalStiffKMat, FDArray1 &forceVec) {
	//zeroFortArray(elemStiffKMat);
	//zeroFortArray(globalStiffKMat);
	//zeroFortArray(forceVec);

	for(int i=1; i<=INMNM; i++) {
		for(int j=1; j<=INMNM; j++) {
			elemStiffKMat(i, j) = 0;
		}
	}

	const int totalNumOfFreedoms = NTF;
	const int halfBandWidth = NB;
	for(int i=1; i<=totalNumOfFreedoms; i++) {
		forceVec(i) = 0;
		for(int j=1; j<=halfBandWidth; j++) {
			globalStiffKMat(i, j) = 0;
		}
	}
}

// 6face-8node-O2
static void newton_cotes_order2(int iElem, FDArray1 &localXi, FDArray1 &localEta, FDArray1 &localZeta,
	FDArray1 &weightH, const FIArray1 &numberOfNodesInElem) 
{
	const int numNodes = 8;

	// set weights
	for(int i=1; i<=numNodes; i++) {
		weightH(i) = 1.0;
	}

	// set positions of quad. points
	const double unit = 1.0;
	// lower 4 nodes
	localXi(1) = -unit; localEta(1) = -unit; localZeta(1) = -unit;
	localXi(2) =  unit; localEta(2) = -unit; localZeta(2) = -unit;
	localXi(3) = -unit; localEta(3) =  unit; localZeta(3) = -unit;
	localXi(4) =  unit; localEta(4) =  unit; localZeta(4) = -unit;
	// upper 4 nodes
	for(int i=1; i<=4; i++) {
		int j = i+4;
		localXi(j) = localXi(i); 
		localEta(j) = localEta(i);
		localZeta(j) = -localZeta(i);
	}
}
void newton_cotes(FDArray1 &localXi, FDArray1 &localEta, FDArray1 &localZeta, 
	FDArray1 &weightH, const FIArray1 &numberOfNodesInElem) 
{
	// TODO newton-cotes order1?
	const int iElem = IE;

	if(numberOfNodesInElem(iElem)==8 && NGAUSS==2) { // 6mentai-8node-2order
		newton_cotes_order2(iElem, localXi, localEta, localZeta, weightH, numberOfNodesInElem);
	} else {
		fatalError("Numer. quad. not supported.");
	}
}


// 3D D-matrix, set DMat
void DMat(const FDArray2 &elastoConstants, 
	FDArray2 &DMat, 
	const FIArray1 &materialTypeOfElem, const FIArray1 &materialCharacters) 
{
	const int iElem = IE; // current element
	const int iMaterial = materialTypeOfElem(iElem);
	const int iMatterType = materialCharacters(iMaterial); // isotropic?

	zeroFortArray(DMat);

	if(iMatterType == MATERIALCHAR_ISOTROPIC) { // isotropic: 1
		const double young = elastoConstants(1, iMaterial);
		const double poisson = elastoConstants(2, iMaterial);
		
		const double c = young / (1.0+poisson) / (1.0-2.0*poisson);
		const double ess = c * (1.0-poisson);
		const double c44 = c * 0.5 * (1.0-2.0*poisson);
		const double c12 = c * poisson;

		// diagonal comp.
		for(int iComp=1; iComp<=3; iComp++) {
			DMat(iComp, iComp) = ess;
			DMat(iComp+3, iComp+3) = c44;
		}
		// deviatoric comp.
		DMat(1,2) = c12; DMat(2,1) = c12;
		DMat(1,3) = c12; DMat(3,1) = c12;
		DMat(2,3) = c12; DMat(3,2) = c12;
	} else { // anisotropic: 2
		fatalError("anisotropic not supported.");
	}
}

// zero DMat first
inline void calcElasticPart(FDArray2 &DMat, 
	const double young, const double poisson) 
{
	const double c = young / (1.0+poisson) / (1.0-2.0*poisson);
	const double ess = c * (1.0-poisson);
	const double c44 = c * 0.5 * (1.0-2.0*poisson);
	const double c12 = c * poisson;

	// diagonal comp.
	for(int iComp=1; iComp<=3; iComp++) {
		DMat(iComp, iComp) = ess;
		DMat(iComp+3, iComp+3) = c44;
	}
	// deviatoric comp.
	DMat(1,2) = c12; DMat(2,1) = c12;
	DMat(1,3) = c12; DMat(3,1) = c12;
	DMat(2,3) = c12; DMat(3,2) = c12;
}
inline void getDeviatoricPart(const int jElem, const int kPoint,
	FDArray1 &devPart, const FDArray3 &stress) 
{
	for(int i=1; i<=ICOM; i++) {
		devPart(i) = stress(i,jElem,kPoint);
	}
	
	double ave = (devPart(1) + devPart(2) + devPart(3)) / 3;
	devPart(1) -= ave;
	devPart(2) -= ave;
	devPart(3) -= ave;
}
inline double calcVonMisesStress(const FDArray1 &devPart) {
	double s11=devPart(1), s22=devPart(2), s33=devPart(3);
	double s12=devPart(4), s13=devPart(5), s23=devPart(6);

	double mises = s11*s11 + s22*s22 + s33*s33 
		+ (s12*s12 + s13*s13 + s23*s23) * 2;
	mises = sqrt(1.5 * mises);

	return mises;
}

/**
 *
 */
int DpMat(const int jElem, const int kPoint,
	FDArray2 &DpMat, const FDArray3 &stress,
	FIArray2 &epState, FDArray2 &equivStress,
	const FDArray2 &elastoConstants, 
	const FIArray1 &materialTypeOfElem, const FIArray1 &materialCharacters) 
{
	int stat = 0; // elastic

	const int jMaterial = materialTypeOfElem(jElem); // 
	const int jMatterType = materialCharacters(jMaterial); // isotropic?

	FDArray1 FORT_ARR(devStress, ICOM);

	zeroFortArray(DpMat);
	zeroFortArray(devStress);

	if(jMatterType == MATERIALCHAR_ISOTROPIC) { // isotropic: 1
		const double young = elastoConstants(1, jMaterial);
		const double poisson = elastoConstants(2, jMaterial);
		const double yield = elastoConstants(3, jMaterial);
		
		// C^e
		calcElasticPart(DpMat, young, poisson);

		// deviatoric stress
		getDeviatoricPart(jElem, kPoint, devStress, stress);

		// von mises
		const double kMisesStress = calcVonMisesStress(devStress);

		if(kMisesStress-yield >= 0) { // yielded
			stat = 1; // plastic
			//printf("Elem%d Point%d yielded\n", jElem, kPoint);


			// calc C^e - C^p matrix
			const double shearMod = 0.5 * young / (1.0+poisson);
			const double cp = 3.0 * shearMod / (kMisesStress*kMisesStress);
			for(int i=1; i<=ICOM; i++) {
				for(int j=1; j<=ICOM; j++) {
					DpMat(i,j) -= devStress(i) * devStress(j) * cp;
				}
			}
		}

		//
		epState(jElem, kPoint) = stat;
		equivStress(jElem, kPoint) = kMisesStress;
	} else { // anisotropic: 2
		fatalError("anisotropic not supported.");
	}

	return stat;
}


/**
 * save D: ICOM x ICOM
 */
void DMMat(const FDArray2 &D, const int quadPointID, FDArray3 &DM, const FIArray1 &numOfNodes) {
	const int iElem = IE;

	for(int i=1; i<=ICOM; i++) {
		for(int j=1; j<=ICOM; j++) {
			DM(quadPointID, i, j) = D(i, j);
		}
	}
}

/**
 * shape function
 */
void integ(double xi, double eta, double zeta, 
	FDArray1 &shapeFuncN, FDArray1 &dNdXi, FDArray1 &dNdEta, FDArray1 &dNdZeta,
	FDArray1 &shapeFuncNInfty, FDArray1 &dNdXiInfty, FDArray1 &dNdEtaInfty, FDArray1 &dNdZetaInfty,
	const FIArray1 &numOfNodesInElem, const FIArray2 &nodesTable, const FIArray1 &elemKind) 
{
	const int iElem = IE;
	const int iNodesNum = numOfNodesInElem(iElem);
	const int iElemKind = elemKind(iElem);

	if(iNodesNum == 8) { // 6-mentai(liner element 8nodes)
		const double a1 = 1.0 + xi;
		const double b1 = 1.0 + eta;
		const double c1 = 1.0 + zeta;
		const double a2 = 1.0 - xi;
		const double b2 = 1.0 - eta;
		const double c2 = 1.0 - zeta;
		const double dd = 1.0 / 8;

		// N
		shapeFuncN(1) = a2 * b2 * c2 * dd;
		shapeFuncN(2) = a1 * b2 * c2 * dd;
		shapeFuncN(3) = a1 * b1 * c2 * dd;
		shapeFuncN(4) = a2 * b1 * c2 * dd;
		shapeFuncN(5) = a2 * b2 * c1 * dd;
		shapeFuncN(6) = a1 * b2 * c1 * dd;
		shapeFuncN(7) = a1 * b1 * c1 * dd;
		shapeFuncN(8) = a2 * b1 * c1 * dd;
		// dNdXi
		dNdXi(1) = -b2 * c2 * dd;
		dNdXi(2) =  b2 * c2 * dd;
		dNdXi(3) =  b1 * c2 * dd;
		dNdXi(4) = -b1 * c2 * dd;
		dNdXi(5) = -b2 * c1 * dd;
		dNdXi(6) =  b2 * c1 * dd;
		dNdXi(7) =  b1 * c1 * dd;
		dNdXi(8) = -b1 * c1 * dd;
		// dNdEta
		dNdEta(1) = -a2 * c2 * dd;
		dNdEta(2) = -a1 * c2 * dd;
		dNdEta(3) =  a1 * c2 * dd;
		dNdEta(4) =  a2 * c2 * dd;
		dNdEta(5) = -a2 * c1 * dd;
		dNdEta(6) = -a1 * c1 * dd;
		dNdEta(7) =  a1 * c1 * dd;
		dNdEta(8) =  a2 * c1 * dd;
		// dNdZeta
		dNdZeta(1) = -a2 * b2 * dd;
		dNdZeta(2) = -a1 * b2 * dd;
		dNdZeta(3) = -a1 * b1 * dd;
		dNdZeta(4) = -a2 * b1 * dd;
		dNdZeta(5) =  a2 * b2 * dd;
		dNdZeta(6) =  a1 * b2 * dd;
		dNdZeta(7) =  a1 * b1 * dd;
		dNdZeta(8) =  a2 * b1 * dd;

		if(iElemKind == ELEMTYPE_INFINITE) {
			// TODO infinite element
			fatalError(__FUNCTION__":infinite element not supported.");
		}
	} else { 
		// 6 sankaku-tyuu(liner element 6nodes)
		// 20
		// 15
		fatalError(__FUNCTION__":Element type not supported.");
	} 
}

/**
 * set J, Jinv, detJ
 */
void JMat(const FDArray1 &x, const FDArray1 &y, const FDArray1 &z,
	const FDArray1 &dNdXi, const FDArray1 &dNdEta, const FDArray1 &dNdZeta,
	const FDArray1 &dNdXiInfty, const FDArray1 &dNdEtaInfty, const FDArray1 &dNdZetaInfty,
	FDArray2 &J, FDArray2 &Jinv, double &detJ,
	const FIArray2 &nodeTable, const FIArray1 &numOfNodes, const FIArray1 &elemType)
{
	const int iElem = IE;
	const int iElemType = elemType(iElem);
	const int iNodesNum = numOfNodes(iElem);

	zeroFortArray(J);

	//const bool isInfElem = iElemType == ELEMTYPE_INFINITE;
	if(iElemType == ELEMTYPE_INFINITE) {
		fatalError("Infty elements not supported.");
	} else {
		//FDArray1 const &dXi = isInfElem ? dNdXiInfty : dNdXi;
		//FDArray1 const &dEta = isInfElem ? dNdEtaInfty : dNdEta;
		//FDArray1 const &dZeta = isInfElem ? dNdZetaInfty : dNdZeta;

		// Jacobian
		for(int jPoint=1; jPoint<=iNodesNum; jPoint++) { // 10
			const int pid = nodeTable(jPoint, iElem);
			double px = x(pid), py = y(pid), pz = z(pid);
			double d1 = dNdXi(jPoint), d2 = dNdEta(jPoint), d3 = dNdZeta(jPoint);
			J(1,1) += d1 * px; 
			J(2,1) += d2 * px;
			J(3,1) += d3 * px;
			J(1,2) += d1 * py;
			J(2,2) += d2 * py;
			J(3,2) += d3 * py;
			J(1,3) += d1 * pz;
			J(2,3) += d2 * pz;
			J(3,3) += d3 * pz;
		} // 10
	}

	// det[J]
	double j11 = J(1,1), j12 = J(1,2), j13 = J(1,3),
		j21 = J(2,1), j22 = J(2,2), j23 = J(2,3),
		j31 = J(3,1), j32 = J(3,2), j33 = J(3,3);

	detJ = j11 * (j22*j33 - j32*j23)
		+ j12 * (j31*j23 - j21*j33)
		+ j13 * (j21*j32 - j31*j22);
	double detJinv = 1.0 / detJ;
	if(fabs(detJ)<1e-9 || fabs(detJinv)<1e-9) {
		fatalError("Jacobian matrix singular.");
	}

	// invert J
	Jinv(1,1) = (j22*j33 - j32*j23) * detJinv;
	Jinv(2,2) = (j11*j33 - j13*j31) * detJinv;
	Jinv(3,3) = (j11*j22 - j12*j21) * detJinv;

	Jinv(2,1) = (j31*j23 - j21*j33) * detJinv;
	Jinv(1,2) = (j13*j32 - j12*j33) * detJinv;
	
	Jinv(3,1) = (j21*j32 - j31*j22) * detJinv;
	Jinv(1,3) = (j12*j23 - j13*j22) * detJinv;

	Jinv(3,2) = (j12*j31 - j32*j11) * detJinv;
	Jinv(2,3) = (j21*j13 - j23*j11) * detJinv;
	
}

// 
void BMat(const FDArray2 &Jinv, FDArray2 &B, const FDArray1 &shapeFuncN,
	const FDArray1 &dNdXi, const FDArray1 &dNdEta, const FDArray1 &dNdZeta,
	const FIArray1 &numOfNodes)
{
	const int iElem = IE;
	const int iNodesNum = numOfNodes(iElem);

	for(int i=1; i<=iNodesNum; i++) { // 10
		const int i1 = 3*i - 2;
		const int i2 = 3*i - 1;
		const int i3 = 3*i;
		
		const double dNd1 = dNdXi(i);
		const double dNd2 = dNdEta(i);
		const double dNd3 = dNdZeta(i);

		const double dNdx = Jinv(1,1)*dNd1 + Jinv(1,2)*dNd2 + Jinv(1,3)*dNd3;
		const double dNdy = Jinv(2,1)*dNd1 + Jinv(2,2)*dNd2 + Jinv(2,3)*dNd3;
		const double dNdz = Jinv(3,1)*dNd1 + Jinv(3,2)*dNd2 + Jinv(3,3)*dNd3;

		// ex,ey,ez
		B(1,i1) = dNdx;	B(1,i2) = 0;	B(1,i3) = 0;
		B(2,i1) = 0;	B(2,i2) = dNdy;	B(2,i3) = 0;
		B(3,i1) = 0;	B(3,i2) = 0;	B(3,i3) = dNdz;
		// exy, eyz, ezx
		B(4,i1) = dNdy;	B(4,i2) = dNdx;	B(4,i3) = 0;
		B(5,i1) = 0;	B(5,i2) = dNdz;	B(5,i3) = dNdy;
		B(6,i1) = dNdz;	B(6,i2) = 0;	B(6,i3) = dNdx;
	} // 10
}

/**
 * save B 
 */
void BMMat(const FDArray2 &B, const int quadPointID, FDArray3 &BM, const FIArray1 &numOfNodes) {
	const int iElem = IE;
	const int iNodesNum = numOfNodes(iElem);
	const int length = iNodesNum * IDIM; // B: ICOM x (nodes * IDIM)

	for(int i=1; i<=ICOM; i++) {
		for(int j=1; j<=length; j++) {
			BM(quadPointID, i, j) = B(i, j);
		}
	}
}

/**
 * C <- BT.D.B
 * Freedom: IDIM * nodes-per-elem
 * B: ICOM*Freedom; BT: Freedom*ICOM; D: ICOM*ICOM
 */
void BtDBMat(const FDArray2 &B, FDArray2 &C, 
	const FDArray2 &D, const FIArray1 &numOfNodes,
	FDArray2 &BT, FDArray2 &BTD)
{
	const int iElem = IE;
	const int iNodesNum = numOfNodes(iElem);
	const int iBMatColumns = iNodesNum * IDIM;

	// transpose B -> BT
	//BT = B.transpose(1, 0);
	for(int i=1; i<=iBMatColumns; i++) {
		for(int j=1; j<=ICOM; j++) {
			BT(i, j) = B(j, i);
		}
	}

	// BTD <- BT.D
	for(int i=1; i<=iBMatColumns; i++) { // 300
		for(int j=1; j<=ICOM; j++) { // 400
			BTD(i, j) = 0;
			for(int k=1; k<=ICOM; k++) { // 500
				BTD(i, j) += BT(i,k) * D(k,j);
			} // 500
		}
	}

	// C <- BTD.B
	for(int i=1; i<=iBMatColumns; i++) { // 600
		for(int j=1; j<=iBMatColumns; j++) { // 700
			C(i,j) = 0;
			for(int k=1; k<=ICOM; k++) {
				C(i,j) += BTD(i,k) * B(k,j);
			}
		}
	}
}

/**
 *
 */
void gauss2(const FDArray2 &BtDB, FDArray2 &elemK, 
	const double weightH, const FDArray1 &shapeFuncN, const double detJ,
	const FIArray1 &numOfNodes, const FIArray2 &nodesTable) 
{
	const int iElem = IE;
	const int iNodesNum = numOfNodes(iElem);
	const int iElemFreedom = IDIM * iNodesNum;

	for(int i=1; i<=iElemFreedom; i++) {
		for(int j=1; j<=iElemFreedom; j++) {
			elemK(i,j) += weightH * detJ * BtDB(i,j);
		}
	}
}

/**
 * Merge element-stiffness matrix with global-stiffness matrix.
 * Lower? triangle with band storage
 */
void merge(const FDArray2 &localK, const FIArray2 &nodesTable,
	FDArray2 &globalK, const FIArray1 &numOfNodes)
{
	const int iElem = IE;
	const int iNodesNum = numOfNodes(iElem);
	const int iFreedomNum = iNodesNum * IDIM;

	FIArray1 FORT_ARR(nf, INMNM); // what's this?
	zeroFortArray(nf);
	
	const int numOfEqns = NAEQ;
	if(numOfEqns != 0) {
		// TODO equation
		fatalError(__FUNCTION__":EQUATIONS not supported.");
	} else { // non-equation version
		for(int i=1; i<=iNodesNum; i++) { // 110
			int ix = IDIM*i - 2;
			int iy = IDIM*i - 1;
			int iz = IDIM*i;

			int nodeId = nodesTable(i, iElem);

			nf(ix) = IDIM*nodeId - 2;
			nf(iy) = IDIM*nodeId - 1;
			nf(iz) = IDIM*nodeId;
		} // 110
	}

	// modified for equation
	const int bandWidth = NB;
	// upper triangle with band storage
	for(int i=1; i<=iFreedomNum; i++) { // 20
		for(int j=1; j<=iFreedomNum; j++) { // 30
			int nfi = nf(i);
			int nfj = nf(j);

			if(nfi <= nfj) {
				int bandPos = nfj - nfi + 1;
				if(bandPos > bandWidth) {
					fatalError("BAND WIDTH ERROR (MERGE)");
				}
				// put in band matrix
				globalK(nfi, bandPos) += localK(i,j);
			} else {
				// ?
			}
		} // 30
	} // 20

}

void clear2(FDArray2 &elemK, const FIArray1 &numOfNodes) {
	zeroFortArray(elemK);
}

/**
 * BC
 */
void bound(const FIArray1 &freedomOfLoadPoint, const FDArray1 &loadOfLoadPoint,
	FDArray1 &forceVector,
	const FIArray1 &freedomOfDispPoint, const FDArray1 &dispOfDispPoint,
	FDArray2 &globalK)
{
	const int totalNumOfFreedomWithLoad = NTFOR;
	const int totalNumOfFreedomWithDisp = NTDIS;
	const int totalNumOfFreedom = NTF;
	const int bandWidth = NB;

	// load force
	for(int i=1; i<=totalNumOfFreedomWithLoad; i++) { // 10
		int iFreedom = freedomOfLoadPoint(i);
		forceVector(iFreedom) = loadOfLoadPoint(i);
	} // 10

	if(NAEQ != 0) {
		// TODO equation
		fatalError("EQUATION unsupported.");
	}

	// displacement
	for(int i=1; i<=totalNumOfFreedomWithDisp; i++) { // 20
		const int iFreedom = freedomOfDispPoint(i);
		const double iDisp = dispOfDispPoint(i);

		for(int j=1; j<=totalNumOfFreedom; j++) { // 30
			if(j <= iFreedom) {
				int absji = abs(iFreedom - j + 1);
				if(absji > bandWidth) {
					//fatalError(__FUNCTION__ STRINGLIZE(__LINE__));
					continue;
				}

				forceVector(j) -= globalK(j,iFreedom-j+1) * iDisp;
			} else {
				int absji = abs(j - iFreedom + 1);
				if(absji > bandWidth) {
					//fatalError(__FUNCTION__ STRINGLIZE(__LINE__));
					continue;
				}

				forceVector(j) -= globalK(iFreedom,j-iFreedom+1) * iDisp;
			}
		} // 30
	} // 20

	// alter solution, fixed displacement is set
	for(int i=1; i<=totalNumOfFreedomWithDisp; i++) { // 40
		int iFreedom = freedomOfDispPoint(i);
		forceVector(iFreedom) = dispOfDispPoint(i);
	} // 40

	// alter coef. matrix, fixed displacement is set to unity
	for(int i=1; i<=totalNumOfFreedomWithDisp; i++) { // 50
		const int iFreedom = freedomOfDispPoint(i);
		
		for(int j=1; j<=iFreedom; j++) { // 60
			int absji = abs(iFreedom - j + 1);
			if(absji > bandWidth) {
				//fatalError(__FUNCTION__ STRINGLIZE(__LINE__));
				continue;
			}

			globalK(j, iFreedom-j+1) = 0;
		} // 60

		for(int j=1; j<=bandWidth; j++) { // 70
			globalK(iFreedom, j) = 0;
		} // 70

		globalK(iFreedom, 1) = 1;
	} // 50
}


/************************************************************************
 *  SIMULTANEOUS LINEAR EQUATIONS WITH REAL SYMMETRIC POSITIVE DEFINITE *
 *      BAND MATRIX BY CHOLESKY METHOD.                                 *
 *  PARAMETERS                                                          *
 *    (1) A : 2-DIM. ARRAY CONTAINING THE MATRIX.                       *
 *    (2) N : ORDER OF THE MATRIX.                                      *
 *    (3) NUD : SIZE OF BAND'S HALF WIDTH.                              *
 *    (4) N1 : ROW SIZE OF THE ARRAY A IN THE 'DIMENSION' STATEMENT.    *
 *    (5) B : 1-DIM. ARRAY CONTAINING THE RIGHT HAND SIDE VECTOR.       *
 *    (6) EPS : PARAMETER TO CHECK SINGURARITY OFF THE MATRIX           *
 *              STANDARD VALUE = 1.0D-14                                *
 *    (7) DR : 1-DIM. WORKING ARRAY.                                    *
 *    (8) Z : 1-DIM. WORKING ARRAY.                                     *
 *    (9) IER : ERROR CODE.                                             *
 *  COPY RIGHT   T. OGUNI   JULY 30 1989   VERSION 1.0                  *
 ************************************************************************ 
 */
static int lisb(FDArray2 &A, int n, int nud, int n1,
	const FDArray1 &b, const double eps, 
	FDArray1 &dr, FDArray1 &z)
{
	int error = 0;
	
	int j;
	double xx, s;
	const int m = nud + 1;

	// check input
	if(n<=0 || nud<=0 || n1<m) {
		error = 2;
		printf("(SUBR. LISB) INVALID ARGUMENT. %d %d %d\n", n, nud, n1);
		goto final;
	}

	/**
	 * modified cholesky decomposition
	 */
	j = 1; // line 32
	if(fabs(A(m,1)) <= eps) {
		error = 1;
		printf("(SUBR. LISB) SINGULAR At STEP # %d\n", j);
		goto final;
	}
	dr(1) = 1.0 / A(m,1);
	xx = A(m-1,2);
	A(m-1,2) = A(m-1,2) * dr(1);
	s = A(m,2) - xx * A(m-1,2);

	//
	j = 2; // line 42
	if(fabs(s) <= eps) {
		error = 1;
		printf("(SUBR. LISB) SINGULAR At STEP # %d\n", j);
		goto final;
	}
	dr(2) = 1.0 / s;

	//
	//printf("(SUBR. LISB) band size %d\n", m);

	if(m < 3) { // line 49
		fatalError("(SUBR. LISB) should never reach here.");

		for(j=3; j<=n; j++) {
			xx = A(1,j);
			A(1,j) = xx * dr(j-1);
			s = A(2,j) - xx * A(1,j);

			if(fabs(s) <= eps) {
				error = 1;
				printf("(SUBR. LISB) SINGULAR At STEP # %d\n", j);
				goto final;
			}

			dr(j) = 1.0 / s;
		}
	} else { // line 61
		for(j=3; j<=n; j++) { // 30
			int k1 = 1;
			if(j >= m) k1 = j - m + 1;
			int mj = m - j;
			double sum = 0;

			for(int i=k1+1; i<=j-1; i++) { // 20
				sum = 0;
				for(int k=k1; k<=i-1; k++) { // 10
					sum += A(m-i+k,i) * A(mj+k,j);
				} // 10
				A(mj+i,j) = A(mj+i,j) -sum;
			} // 20

			sum = 0;
			for(int i=k1; i<=j-1; i++) { // 25
				xx = A(mj+i,j);
				double au = xx * dr(i);
				sum += xx * au;
				A(mj+i,j) = au;
			}

			double t = A(m,j) - sum;
			if(fabs(t) <= eps) {
				error = 1;
				printf("(SUBR. LISB) SINGULAR At STEP # %d\n", j);
				goto final;
			}

			dr(j) = 1.0 / t;
		} // 30
	}

final:
	return error;
}

/**
 * GK: IGK * IBN
 * GKC: IBN * IGK
 */
void solve_cholesky_decomp(
	const FDArray2 &globalK, FDArray2 &decompBuf, 
	const FDArray1 &forceVector, const FDArray1 &dispVector, 
	FDArray1 &dr, FDArray1 &zm) 
{
	const int totalNumOfFreedom = NTF;
	const int bandWidth = NB;

	for(int j=1; j<=totalNumOfFreedom; j++) { // 200
		for(int i=1; i<=bandWidth; i++) { // 100
			if(i+j <= bandWidth) {
				decompBuf(i,j) = 0;
			} else {
				int i1 = bandWidth + 1 - i;
				int i2 = j + 1 - i1;
				decompBuf(i,j) = globalK(i2,i1);
			}
		} // 100
	} // 200

	const int nud = bandWidth - 1;
	const double eps = 1.0e-12;

	//printf("%lf\n", decompBuf(nud+1,1));

	int status = lisb(decompBuf, totalNumOfFreedom, nud, bandWidth, 
		forceVector, eps, dr, zm);
	if(status != 0) {
		fatalError("Cholesky Decomposition failure.");
	}
}

/**
 * A: n1 x n
 */
static void sbsub(const FDArray2 &A, const int n, const int nud, const int n1,
	FDArray1 &b, const FDArray1 &dr, FDArray1 &z)
{
	const int m = nud + 1;

	if(m < 3) {
		// TODO m<3
		fatalError(__FUNCTION__" should never reach here.");
	} else {
		// forward
		z(1) = b(1);
		z(2) = b(2) - A(m-1,2) * z(1);
		
		for(int j=3; j<=n; j++) { // 80
			int i1 = j>m ? 1 : m-j+1;
			
			double sum = 0;
			for(int k=i1; k<=m-1; k++) { // 70
				sum += A(k,j) * z(j-m+k);
			}

			z(j) = b(j) - sum;
		} // 80

		for(int j=1; j<=n; j++) { // 90
			z(j) = z(j) * dr(j);
		}

		// backward
		b(n) = z(n);
		b(n-1) = z(n-1) - A(m-1,n) * z(n);

		for(int j=3; j<=n; j++) { // 110
			int j1 = n-j+1;
			int i1 = j<m ? m-j+1 : 1;

			double sum = 0;
			for(int k=i1; k<=m-1; k++) {
				sum += A(k,m-k+j1) * b(m-k+j1);
			}

			b(j1) = z(j1) - sum;
		} // 110
	}
}

// 
void solve_cholesky_disp(const FDArray2 &gkc, FDArray1 &f, 
	FDArray1 &u, const FDArray1 &dr, FDArray1 &zm)
{
	const int totalNumberOfFreedom = NTF;
	const int bandWidth = NB;
	const int nud = bandWidth - 1;

	sbsub(gkc, totalNumberOfFreedom, nud, bandWidth, f, dr, zm);

	for(int i=1; i<=totalNumberOfFreedom; i++) { // 300
		u(i) = f(i);
	}
}

// 
void strain(
	const FIArray2 &nodesTable, FDArray1 &displace, FDArray3 &strain, const FDArray3 &BMatrices, 
	const FIArray1 &numOfNodes, const FIArray1 &numOfQuadPoints, FDArray3 &strainIncr) 
{
	if(NAEQ != 0) {
		// TODO equation
		fatalError("EQUATION not supported.");
	}

	const int totalNumberOfElem = NTE;

	int nas = 0;
	FIArray1 FORT_ARR(nf, INMNM);
	zeroFortArray(nf);

	for(int iElem=1; iElem<=totalNumberOfElem; iElem++) { // 5
		const int iNodesNum = numOfNodes(iElem);
		const int iQuadPointsNum = numOfQuadPoints(iElem);

		for(int iNode=1; iNode<=iNodesNum; iNode++) {
			nf(3*iNode-2) = nodesTable(iNode,iElem) * IDIM - 2;
			nf(3*iNode-1) = nodesTable(iNode,iElem) * IDIM - 1;
			nf(3*iNode) = nodesTable(iNode,iElem) * IDIM;
		}

		for(int iPoint=1; iPoint<=iQuadPointsNum; iPoint++) { // 105
			nas += 1;

			for(int i=1; i<=ICOM; i++) { // 110
				strainIncr(i, iElem, iPoint) = 0;
			}

			for(int i=1; i<=ICOM; i++) {
				for(int j=1; j<=INMNM; j++) {
					int nfj = nf(j);
					strainIncr(i,iElem,iPoint) += BMatrices(nas,i,j) * displace(nfj);
				}
			}

			for(int i=1; i<=ICOM; i++) {
				strain(i,iElem,iPoint) += strainIncr(i,iElem,iPoint);
			}
		} // 105
	} // 5

}

//// old version, elasticity only
//void stress(FDArray2 &D, const FDArray3 &strainIncr, 
//	FDArray3 &stress, FDArray3 &stressIncr, const FIArray1 &numOfQuadPoints) 
//{
//	const int totalNumOfElem = NTE;
//
//	for(int iElem=1; iElem<=totalNumOfElem; iElem++) { // 100
//		const int iQuadPointNum = numOfQuadPoints(iElem);
//
//		for(int kk=1; kk<=iQuadPointNum; kk++) { // 105
//			for(int i=1; i<=ICOM; i++) {
//				stressIncr(i,iElem,kk) = 0;
//			}
//
//			for(int i=1; i<=ICOM; i++) {
//				for(int j=1; j<=ICOM; j++) {
//					stressIncr(i,iElem,kk) += D(i,j) * strainIncr(j,iElem,kk);
//				}
//			}
//
//			for(int i=1; i<=ICOM; i++) {
//				stress(i,iElem,kk) += stressIncr(i,iElem,kk);
//			}
//		} // 105
//	} // 100
//}

/*
 *
 */
void stress(FDArray3 &DM, const FDArray3 &strainIncr, 
	FDArray3 &stress, FDArray3 &stressIncr, const FIArray1 &numOfQuadPoints) 
{
	const int totalNumOfElem = NTE;

	int kQuadPoint = 0;

	for(int iElem=1; iElem<=totalNumOfElem; iElem++) { // 100
		const int iQuadPointNum = numOfQuadPoints(iElem);

		for(int kk=1; kk<=iQuadPointNum; kk++) { // 105
			// global id for current quad. point
			kQuadPoint += 1;

			for(int i=1; i<=ICOM; i++) {
				stressIncr(i,iElem,kk) = 0;
			}

			for(int i=1; i<=ICOM; i++) {
				for(int j=1; j<=ICOM; j++) {
					stressIncr(i,iElem,kk) += DM(kQuadPoint,i,j) * strainIncr(j,iElem,kk);
				}
			}

			for(int i=1; i<=ICOM; i++) {
				stress(i,iElem,kk) += stressIncr(i,iElem,kk);
			}
		} // 105
	} // 100
}

void udata(const FDArray1 &du, FDArray1 &u) {
	const MKL_INT numberOfFreedoms = NTF;
	for(int i=1; i<=numberOfFreedoms; i+=IDIM) {
		u(i) += du(i);
		u(i+1) += du(i+1);
		u(i+2) += du(i+2);
	}
}

void utotal(FDArray1 &finalDisp, const FDArray1 &u, const int nstep) {
	const MKL_INT numberOfFreedoms = NTF;
	for(int i=1; i<=numberOfFreedoms; i+=IDIM) {
		finalDisp(i) = u(i);
		finalDisp(i+1) = u(i+1);
		finalDisp(i+2) = u(i+2);
	}
}

/**
 *
 * c  NSGAV(i) :  summation number of global node i
 * c  SGN(J,I) :  j:component i:node number(global)
 */
void output_av(
	const FDArray1 &disp, const FDArray3 &pointStrain, const FDArray3 &pointStress,
	const FIArray1 &numOfQuadPoints, const FIArray1 &numOfNodes, 
	const FIArray1 &userElems, const FIArray1 &userNodes, 
	const FIArray2 &nodeTable,
	FDArray2 &nodeStress, FDArray2 &nodeStrain, FIArray1 &nsgav) 
{
	const int totalNumOfNode = NTN;
	const int totalNumOfElem = NTE;
	const int totalNumOfFree = NTF;
	
	//c  SGLOCAL(JJ,II) local stress -> II:node number(local),JJ:component
	FDArray2 FORT_ARR(localStress, ICOM,8), FORT_ARR(localStrain, ICOM,8);
	zeroFortArray(localStress);
	zeroFortArray(localStrain);
	
	FDArray2 FORT_ARR(qc88, 8,8), FORT_ARR(qc89, 8,8);
	{
		qc88(1,1) = 0.2549039E+01; qc88(1,2) = -0.6830133E+00; qc88(1,3) = 0.1830131E+00; qc88(1,4) = -0.6830133E+00; 
		qc88(1,5) = -0.6830133E+00; qc88(1,6) = 0.1830131E+00; qc88(1,7) = -0.4903840E-01; qc88(1,8) = 0.1830131E+00;

		qc88(2,2) = 0.2549039E+01; qc88(2,3) = -0.6830133E+00; qc88(2,4) = 0.1830131E+00; qc88(2,5) = 0.1830131E+00;
		qc88(2,6) = -0.6830133E+00; qc88(2,7) = 0.1830131E+00; qc88(2,8) = -0.4903840E-01;

		qc88(3,3) = 0.2549039E+01; qc88(3,4) = -0.6830133E+00; qc88(3,5) = -0.4903840E-01; qc88(3,6) = 0.1830131E+00;
		qc88(3,7) = -0.6830133E+00; qc88(3,8) = 0.1830131E+00; 

		qc88(4,4) = 0.2549039E+01; qc88(4,5) = 0.1830131E+00; qc88(4,6) = -0.4903840E-01; qc88(4,7) = 0.1830131E+00; 
		qc88(4,8) = -0.6830133E+00; 
		
		qc88(5,5) = 0.2549039E+01; qc88(5,6) = -0.6830133E+00; qc88(5,7) = 0.1830131E+00; qc88(5,8) = -0.6830133E+00; 

		qc88(6,6) = 0.2549039E+01; qc88(6,7) = -0.6830133E+00; qc88(6,8) = 0.1830131E+00;

		qc88(7,7) = 0.2549039E+01; qc88(7,8) = -0.6830133E+00; 
		
		qc88(8,8) = 0.2549039E+01;
	} {
		qc89(1,1) = 0.1503080E+01; qc89(1,2) = -0.1909162E+00; qc89(1,3) = 0.2424952E-01; qc89(1,4) = -0.1909162E+00; 
		qc89(1,5) = -0.1909162E+00; qc89(1,6) = 0.2424952E-01; qc89(1,7) = -0.3080088E-02; qc89(1,8) = 0.2424952E-01; 

		qc89(2,2) = 0.1503080E+01; qc89(2,3) = -0.1909162E+00; qc89(2,4) = 0.2424952E-01; qc89(2,5) = 0.2424952E-01; 
		qc89(2,6) = -0.1909162E+00; qc89(2,7) = 0.2424952E-01; qc89(2,8) = -0.3080088E-02; 

		qc89(3,3) = 0.1503080E+01; qc89(3,4) = -0.1909162E+00; qc89(3,5) = -0.3080088E-02; qc89(3,6) = 0.2424952E-01; 
		qc89(3,7) = -0.1909162E+00; qc89(3,8) = 0.2424952E-01; 

		qc89(4,4) = 0.1503080E+01; qc89(4,5) = 0.2424952E-01; qc89(4,6) = -0.3080088E-02; qc89(4,7) = 0.2424952E-01; 
		qc89(4,8) = -0.1909162E+00; 

		qc89(5,5) = 0.1503080E+01; qc89(5,6) = -0.1909162E+00; qc89(5,7) = 0.2424952E-01; qc89(5,8) = -0.1909162E+00; 

		qc89(6,6) = 0.1503080E+01; qc89(6,7) = -0.1909162E+00; qc89(6,8) = 0.2424952E-01; 

		qc89(7,7) = 0.1503080E+01; qc89(7,8) = -0.1909162E+00; 

		qc89(8,8) = 0.1503080E+01; 
	} {
		for(int i=1; i<=7; i++) {
			for(int j=i+1; j<=8; j++) {
				qc88(j,i) = qc88(i,j);
				qc89(j,i) = qc89(i,j);
			}
		}
	}

	for(int k=1; k<=totalNumOfNode; k++) {
		nsgav(k) = 0;
		for(int j=1; j<=ICOM; j++) {
			nodeStress(j,k) = 0;
			nodeStrain(j,k) = 0;
		}
	}

	for(int i=1; i<=8; i++) {
		for(int j=1; j<=ICOM; j++) {
			localStress(j,i) = 0;
			localStrain(j,i) = 0;
		}
	}

	for(int iElem=1; iElem<=totalNumOfElem; iElem++) { // 100
		const int iNodeNum = numOfNodes(iElem);
		const int iPointNum = numOfQuadPoints(iElem);

		if((iNodeNum==8 || iNodeNum==20) && iPointNum==8) {
			for(int ii=1; ii<=8; ii++) { // 120
				for(int jj=1; jj<=ICOM; jj++) { // 130
					localStress(jj,ii) += qc88(ii,1) * pointStress(jj,iElem,1);
					localStress(jj,ii) += qc88(ii,2) * pointStress(jj,iElem,2);
					localStress(jj,ii) += qc88(ii,3) * pointStress(jj,iElem,4);
					localStress(jj,ii) += qc88(ii,4) * pointStress(jj,iElem,3);
					localStress(jj,ii) += qc88(ii,5) * pointStress(jj,iElem,5);
					localStress(jj,ii) += qc88(ii,6) * pointStress(jj,iElem,6);
					localStress(jj,ii) += qc88(ii,7) * pointStress(jj,iElem,8);
					localStress(jj,ii) += qc88(ii,8) * pointStress(jj,iElem,7);
					//
					localStrain(jj,ii) += qc88(ii,1) * pointStrain(jj,iElem,1);
					localStrain(jj,ii) += qc88(ii,2) * pointStrain(jj,iElem,2);
					localStrain(jj,ii) += qc88(ii,3) * pointStrain(jj,iElem,4);
					localStrain(jj,ii) += qc88(ii,4) * pointStrain(jj,iElem,3);
					localStrain(jj,ii) += qc88(ii,5) * pointStrain(jj,iElem,5);
					localStrain(jj,ii) += qc88(ii,6) * pointStrain(jj,iElem,6);
					localStrain(jj,ii) += qc88(ii,7) * pointStrain(jj,iElem,8);
					localStrain(jj,ii) += qc88(ii,8) * pointStrain(jj,iElem,7);
				} // 130

				const int nodeId = nodeTable(ii, iElem);
				nsgav(nodeId) += 1;
				for(int jj=1; jj<=ICOM; jj++) { // 140
					nodeStress(jj,nodeId) += localStress(jj,ii);
					nodeStrain(jj,nodeId) += localStrain(jj,ii);
				} // 140
			} // 120
		} else {
			fatalError("Average not supported of this element.");
		}

		for(int i=1; i<=8; i++) {
			for(int j=1; j<=ICOM; j++) {
				localStress(j,i) = 0;
				localStrain(j,i) = 0;
			}
		}

	} // 100

	for(int i=1; i<=totalNumOfNode; i++) {
		//const int iUserElem = userElems(i);
		int iAveCount = nsgav(i);

		if(iAveCount != 0) {
			for(int j=1; j<=ICOM; j++) {
				nodeStress(j,i) = nodeStress(j,i) / iAveCount;
				nodeStrain(j,i) = nodeStrain(j,i) / iAveCount;
			}
		} else {
			char msg[512]; sprintf_s(msg, "Node %d is lonely node.", i);
			fatalError(msg);
		}
	}
}

// TODO average values
void output(FILE *fp,
	const FDArray1 &x, const FDArray1 &y, const FDArray1 &z, const FDArray1 &disp,
	const FDArray3 &pointStrain, const FDArray3 &pointStress, 
	const FIArray1 &numOfQuadPoints, const FIArray1 &userElems, const FIArray1 &userNodes,
	const FDArray2 &nodeStress, const FDArray2 &nodeStrain, const FIArray1 &nsgav)
{
	const int totalNumOfNode = NTN;
	const int totalNumOfElem = NTE;
	const int totalNumOfFree = NTF;

	int in = 0;

#define PRINT(fmt, ...) fprintf_s(fp, (fmt), __VA_ARGS__)
#define PRINTLN(fmt, ...) fprintf_s(fp, (fmt "\n"), __VA_ARGS__)
//#define DBL_FMT "%lE"

	PRINTLN("NODE,ELEMENT,INT.POINT");
	PRINTLN("%d %d %d", totalNumOfNode, totalNumOfElem, numOfQuadPoints(1));
	PRINTLN("DISPLACEMENT");
	PRINTLN("NODE,Ux,Uy,Uz");

	in = 0;
	for(int i=1; i<=totalNumOfFree; i+=IDIM) {
		in += 1;
		int iUserNode = userNodes(in);
		PRINTLN("%d %lE %lE %lE", iUserNode, disp(i), disp(i+1), disp(i+2));
	}

	PRINTLN("S T R E S S");
	PRINTLN("Element,Int.point,Sx,Sy,Sz,Sxy,Sxz,Syz");
	for(int i=1; i<=totalNumOfElem; i++) {
		const int iPointNum = numOfQuadPoints(i);
		const int iUser = userElems(i);

		for(int k=1; k<=iPointNum; k++) {
			PRINT("%d %d ", iUser, k);
			for(int j=1; j<=ICOM; j++) {
				PRINT("%lE ", pointStress(j,i,k));
			}
			PRINTLN("");
		}
	}

	PRINTLN("S T R A I N");
	PRINTLN("Element,Int.point,ex,ey,ez,exy,exz,eyz");
	for(int i=1; i<=totalNumOfElem; i++) {
		const int iPointNum = numOfQuadPoints(i);
		const int iUser = userElems(i);

		for(int k=1; k<=iPointNum; k++) {
			PRINT("%d %d ", iUser, k);
			for(int j=1; j<=ICOM; j++) {
				PRINT("%lE ", pointStrain(j,i,k));
			}
			PRINTLN("");
		}
	}

	PRINTLN("AVERAGED AT NODE (STRESS)");
	for(int i=1; i<=totalNumOfNode; i++) {
		const int iUser = userNodes(i);
		
		PRINTLN("%d %lE %lE %lE %lE %lE %lE", iUser, 
			nodeStress(1,i),nodeStress(2,i),nodeStress(3,i),
			nodeStress(4,i),nodeStress(5,i),nodeStress(6,i));
	}
	PRINTLN("");

	PRINTLN("AVERAGED AT NODE (STRAIN)");
	for(int i=1; i<=totalNumOfNode; i++) {
		const int iUser = userNodes(i);
		
		PRINTLN("%d %lE %lE %lE %lE %lE %lE", iUser, 
			nodeStrain(1,i),nodeStrain(2,i),nodeStrain(3,i),
			nodeStrain(4,i),nodeStrain(5,i),nodeStrain(6,i));
	}
	PRINTLN("");

	PRINTLN("NODE COORDINATE");
	for(int i=1; i<=totalNumOfNode; i++) {
		int iUserNode = userNodes(i);
		PRINTLN("%d %lE %lE %lE", iUserNode, x(i), y(i), z(i));
	}

	//in = 0;
	//for(int i=1; i<=totalNumOfFree; i+=IDIM) {
	//	in += 1;
	//	int iUserNode = userNodes(in);
	//	PRINTLN("%d %lE %lE %lE %lE %lE %lE", 
	//		iUserNode, disp(i), disp(i+1), disp(i+2), x(in), y(in), z(in));
	//}

#undef PRINT
#undef PRINTLN
}

void output_vtu(const char *filename, const int iStep,
	const FDArray1 &x, const FDArray1 &y, const FDArray1 &z, const FDArray1 &disp,
	const FDArray3 &strain, const FDArray3 &stress, 
	const FIArray1 &numOfNodes, const FIArray2 &nodesTable,
	const FIArray1 &numOfQuadPoints, /*const FIArray2 &quadPointsTable,*/
	const FIArray1 &userElems, const FIArray1 &userNodes,
	const FDArray2 &nodeStress, const FDArray2 &nodeStrain, const FIArray1 &nsgav, 
	const FIArray1 &freedomOfLoadPoint, const FDArray1 &loadOfLoadPoint,
	const FIArray2 &epStateOfQuadPoints
	)
{
	const int totalNumOfNode = NTN;
	const int totalNumOfElem = NTE;
	const int totalNumOfFree = NTF;

	vtkSmartPointer<vtkUnstructuredGrid> grid
		= vtkHelper_createGrid<vtkUnstructuredGrid>(totalNumOfNode, iStep);
	vtkSmartPointer<vtkPoints> points = grid->GetPoints();

	NEW_VTKOBJ(vtkCellArray, cells);
	for(int i=1; i<=totalNumOfElem; i++) {
		const int iNodeNum = numOfNodes(i);
		if(iNodeNum != 8) {
			fatalError("Element is not Hexahedron!");
		}

		NEW_VTKOBJ(vtkHexahedron, cell);
		cell->GetPointIds()->SetNumberOfIds(8);

		const int connectivity[] = {
			2, 3, 4, 1,
			6, 7, 8, 5,
		};

		cell->GetPointIds()->SetId(0, nodesTable(2,i)-1);
		cell->GetPointIds()->SetId(1, nodesTable(3,i)-1);
		cell->GetPointIds()->SetId(2, nodesTable(4,i)-1);
		cell->GetPointIds()->SetId(3, nodesTable(1,i)-1);

		cell->GetPointIds()->SetId(4, nodesTable(6,i)-1);
		cell->GetPointIds()->SetId(5, nodesTable(7,i)-1);
		cell->GetPointIds()->SetId(6, nodesTable(8,i)-1);
		cell->GetPointIds()->SetId(7, nodesTable(5,i)-1);

		cells->InsertNextCell(cell);

		if(false && i==1) { 
			std::cout << i << ": ";
			// 31-21-22-32, 11-1-2-12
			for(int j=1; j<=iNodeNum; j++) {
				int jNode = nodesTable(j,i);
				std::cout << jNode << ", ";
			}
			std::cout << std::endl;
		}
	}
	grid->SetCells(VTK_HEXAHEDRON, cells);
	
	vtkSmartPointer<vtkDoubleArray> dispArray
		= vtkHelper_declareField<vtkDoubleArray>(grid, "disp", 3, totalNumOfNode);
	vtkSmartPointer<vtkDoubleArray> nodeStrainArray
		= vtkHelper_declareField<vtkDoubleArray>(grid, "nStrain", 9, totalNumOfNode);
	vtkSmartPointer<vtkDoubleArray> nodeStressArray
		= vtkHelper_declareField<vtkDoubleArray>(grid, "nStress", 9, totalNumOfNode);

	vtkSmartPointer<vtkDoubleArray> loadArray
		= vtkHelper_declareField<vtkDoubleArray>(grid, "load", 3, totalNumOfNode);

	int in = 0;
	for(int i=1; i<=totalNumOfFree; i+=IDIM) {
		in += 1;
		int idx = in - 1;

		points->SetPoint(idx, x(in), y(in), z(in));

		dispArray->SetTuple3(idx, disp(i), disp(i+1), disp(i+2));
		{
			double a11 = nodeStrain(1,in), a22 = nodeStrain(2,in), a33 = nodeStrain(3,in);
			double a12 = nodeStrain(4,in), a13 = nodeStrain(5,in), a23 = nodeStrain(6,in);
			nodeStrainArray->SetTuple9(idx, a11,a12,a13, a12,a22,a23, a13,a23,a33);
		} {
			double a11 = nodeStress(1,in), a22 = nodeStress(2,in), a33 = nodeStress(3,in);
			double a12 = nodeStress(4,in), a13 = nodeStress(5,in), a23 = nodeStress(6,in);
			nodeStressArray->SetTuple9(idx, a11,a12,a13, a12,a22,a23, a13,a23,a33);
		} {
			// we zero the BC data first
			loadArray->SetTuple3(idx, 0, 0, 0);
		}
	}

	for(int iBC=1; iBC<=NTFOR; iBC++) {
		int iFree = freedomOfLoadPoint(iBC);
		int iNode = (iFree-1) / IDIM + 1;
		int iCoord = (iFree-1) % IDIM + 1;
		
		int idx = iNode - 1;
		double iLoad = loadOfLoadPoint(iBC);
		
		double val[IDIM];
		loadArray->GetTuple(idx, val);
		val[iCoord-1] = iLoad;
		loadArray->GetTuple3(idx)[iCoord-1] = iLoad;
		loadArray->SetTuple(idx, val);

		//std::cout << iNode << "," << iCoord << "," << iLoad << std::endl;
	}

	/*
	 * cell data
	 */
	NEW_VTKOBJ(vtkDoubleArray, epStateArray);
	epStateArray->SetName("e/p");
	epStateArray->SetNumberOfComponents(1);
	epStateArray->SetNumberOfTuples(totalNumOfElem);
	grid->GetCellData()->AddArray(epStateArray);

	for(int i=1; i<=totalNumOfElem; i++) {
		const int iNodeNum = numOfNodes(i);

		double iStat = 0;
		for(int k=1; k<=iNodeNum; k++) {
			iStat += epStateOfQuadPoints(i,k);
		}
		iStat /= iNodeNum;

		epStateArray->SetValue(i-1, iStat);
	}

	// save
	if(vtkHelper_saveGrid(grid, filename) == 0) {
		fatalError("Failed to write vtu file.");
	}
}


static void testLinearAlgebra() {
	//typedef ublas::matrix<double, ublas::column_major> FMat;
	//typedef ublas::vector<double> FVec;

	//{
	//	const int n = 3;
	//	FMat mat(n,n);
	//	for(int i=0; i<n*n; i++) {
	//		mat(i/n, i%n) = i;
	//	}

	//	std::cout << mat << std::endl;

	//	double *data = &mat(0,0);
	//	for(int i=0; i<n*n; i++) {
	//		std::cout << data[i] << " ";
	//	}
	//	std::cout << std::endl;
	//}
	{
		int n = 2;
		bz::Array<double, 2> mat(n, n, blitz::fortranArray);
		for(int i=0; i<n*n; i++) {
			mat(i/n+1,i%n+1) = i;
		}
		//mat = 1, 2, 3, 4;
		std::cout << mat << std::endl;

		//double *data = &mat(1, 1);
		double *data = mat.data();
		for(int i=0; i<n*n; i++) {
			std::cout << data[i] << " ";
		}
		std::cout << std::endl;

		FDArray2 matT(n, n, blitz::fortranArray);
		matT = mat.transpose(1, 0);
		std::cout << matT << std::endl;
		for(int i=0; i<n*n; i++) {
			std::cout << matT.data()[i] << " ";
		}
	}
}
int main(int argc, char **argv) {

	/**
	 * declare varaibles
	 */

	// Xi,Eta,Zeta
	FDArray1 FORT_ARR(dNdXi, INMTE), FORT_ARR(dNdEta, INMTE), FORT_ARR(dNdZeta, INMTE);
	// x,y,z
	FDArray1 FORT_ARR(x, INODE), FORT_ARR(y, INODE), FORT_ARR(z, INODE);

	//c  B(i,j)  : elemental stiffness matrix
	FDArray2 FORT_ARR(B, ICOM,INMNM);
	//c  AXI(i)  : x-coordinate of i-th integration point ( -1 < AXI  < 1 )
	//c  AETA(i) : y-coordinate of i-th integration point ( -1 < AETA < 1 )
	//c  AZETA(i): z-coordinate of i-th integration point ( -1 < AZETA < 1 )
	FDArray1 FORT_ARR(aXi, IRD), FORT_ARR(aEta, IRD), FORT_ARR(aZeta, IRD);
	//c  H(i)   : weight parameter of i-th integration point
	FDArray1 FORT_ARR(H, IRD);
	
	//c  D(i,j) : elasto-plastic matrix (1-x,2-y,3-z,4-xy,5-xz,6-yz)
	FDArray2 FORT_ARR(D, ICOM,ICOM);
	// DM(i,*,*): i-th point's elastic matrix 
	FDArray3 FORT_ARR(DM, INA,ICOM,ICOM);

	//c  C(i,j)  : BT(INMNM,6)*D(6,6)*B(6,INMNM) Bt*D*B
	//c  EK(i,j) : element stiffness matrix 
	//c  NN(i,j) : node number of i-th local node of element j 
	FDArray2 FORT_ARR(C, INMNM,INMNM), FORT_ARR(EK, INMNM,INMNM);
	FIArray2 FORT_ARR(NN, INMTE,IELM);

	//c  GK(i,j) : global stiffness matrix  i: matrix size  j:half band width
	//c  F(i)    : force vector of node i
	//c  IFOR(i) : freedom number of i-th initial load point
	FDArray2 FORT_ARR(GK, IGK,IBN);
	FDArray1 FORT_ARR(F, IGK);
	FIArray1 FORT_ARR(IFor, IBOUND);
	//c  FB(i)   : load value of i-th initial load point --> Load increment vector ! modified
	//c  fbnow(i) : present force vector  ! added
	FDArray1 FORT_ARR(FB, IBOUND), FORT_ARR(FBNow, IBOUND);

	//c  IDIS(i) : freedom number of i-th constrained disp. point
	FIArray1 FORT_ARR(IDis, IBOUND);
	//c  UB(i)   : disp. value of i-th constrained disp. point
	//c  U(i)    : displacement value of i-th freedom --> displcement increment vector ! modified
	//c  UNOW(i) : present displacement vector  ! added
	FDArray1 FORT_ARR(UB, IBOUND), FORT_ARR(U, IGK), FORT_ARR(UNow, IGK);
	//c  EP(i,j,k) : i-th component strain at j-th element
	//c  EPDASH(i,j,k) : strain increment ! added
	//c                k-th local integration point
	FDArray3 FORT_ARR(EP, ICOM,IELM,IRD), FORT_ARR(EPDash, ICOM,IELM,IRD);
	//c  BM(i,j,k)   : For strain (B matrix for i-th integration point)
	FDArray3 FORT_ARR(BM, INA,ICOM,INMNM);
	//c  SG(i,j,k) : i-th component stress at j-th element for k-th local integration point
	FDArray3 FORT_ARR(SG, ICOM,IELM,IRD), FORT_ARR(SGDASH, ICOM,IELM,IRD);
	//c  BT(i,j)   : transposed B(6,INMNM) matrix 
	//c  BTD(i,j)  : BT(INMNM,6)*D(6,6)
	FDArray2 FORT_ARR(BT, INMNM,ICOM), FORT_ARR(BTD, INMNM,ICOM);

	//c  ECON(i,j) : i-th elastic constants of j-th material 
	FDArray2 FORT_ARR(ECon, 9,IMNUM);
	//c  IECON(i)  : ELASTIC CHARACTER of i-th material, 1:isotropic; 2:anisotropic
	FIArray1 FORT_ARR(IECon, IMNUM);
	//c  NOM(i)   : material number of i-th element
	//c  NKIND(i) : kind of element 1--> infinite element
	FIArray1 FORT_ARR(NoM, IELM), FORT_ARR(NKind, IELM);
	//c  NMTE(i) : number of a node in i-th element
	//c  NGNG(i) : number of integration point in i-th element
	FIArray1 FORT_ARR(NMTE, IELM), FORT_ARR(NGNG, IELM);
	//c  NEUS(i) : i : PROGRAM NODE NUMBER , NEUS(I) : USER ELEMENT NUMBER 
	//c  NNUS(i) : i : PROGRAM NODE NUMBER , NNUS(I) : USER NODE NUMBER 
	FIArray1 FORT_ARR(NEUS, IELM), FORT_ARR(NNUS, INODE);

	//c  NSGAV(i) :  summation number of global node i
	//c  EPN(,)   :  ?
	//c  SGN(J,I) :  j:component i:node number(global)
	FIArray1 FORT_ARR(NSGAV, INODE);
	FDArray2 FORT_ARR(SGN, ICOM,INODE), FORT_ARR(EPN, ICOM,INODE);

	//ccc FOR CHOLESKY METHOD
	FDArray2 FORT_ARR(GKC, IBN,IGK);
	FDArray1 FORT_ARR(DR, IGK), FORT_ARR(ZM, IGK);

	//c  N    : shape function Ni
	FDArray1 FORT_ARR(N, INMTE);
	//c  JINV : Inverse of Jacobian matrix
	//c  JM   : Jacobian matrix
	FDArray2 FORT_ARR(JInv, IDIM,IDIM), FORT_ARR(JM, IDIM,IDIM);

	//ccc FOR INFINITE ELEMENT
	FDArray1 FORT_ARR(NIf, INMTE);
	FDArray1 FORT_ARR(DXIf, INMTE), FORT_ARR(DETIf, INMTE), FORT_ARR(DZEIf, INMTE);

	// for elasto-plastic behavior on each quad. point
	// 0: elastic, 1: plastic
	FIArray2 FORT_ARR(epState, IELM,IRD);
	FDArray2 FORT_ARR(misesStress, IELM,IRD);

	// 
	fpDatacheck = fopen("datacheck.dat", "w");
	fpResult = fopen("result.dat", "w");
	if(!fpDatacheck || !fpResult) {
		fatalError("IO failure.");
		return 1;
	}
	

	if(1) {
		int nStep = 1000;
		data("input.dat", fpDatacheck, 
			NMTE, ECon, NN, x, y, z,
			IFor, FB, IDis, UB, NGNG,
			NoM, IECon, NEUS, NNUS, NKind,
			nStep);

		// estimate (half) band with of the global matrix
		band(NN, NMTE, fpDatacheck);

		// intern some global variables
		const int TotalNumberOfElems = NTE;
		const int TotalNumberOfNodes = NTN;
		const int TotalNumberOfFrees = NTF;

		check();

		init(SG, EP);
		zeroFortArray(FBNow);
		zeroFortArray(epState);

		// calculate load increment
		for(int i=1; i<=NTFOR; i++) {
			FB(i) /= nStep;
		}

		// load incr loop
		for(int iStep=1; iStep<=nStep; iStep++) { // incr
			clear(EK, GK, F);

			int currentQuadPoint = 0; // through number of integration point

			for(int iElem=1; iElem<=TotalNumberOfElems; iElem++) { // 150
				// set the global element id-of-interest
				// TODO remove the global IE
				IE = iElem;

				// setup weights/positions of gauss/newton-cotes quadrature 
				newton_cotes(aXi, aEta, aZeta, H, NMTE);

				// for all quad. points
				const int numberOfQuadPointsInElem = NGNG(iElem);
				for(int kPoint=1; kPoint<=numberOfQuadPointsInElem; kPoint++) { // 100
					const double xi = aXi(kPoint);
					const double eta = aEta(kPoint);
					const double zeta = aZeta(kPoint);
					const double hNum = H(kPoint);

					double detJ = 0;

					// count current quad. point
					currentQuadPoint += 1;

					// elastic matrix
					//// TODO change this subroutine to dpmat.f
					//DMat(ECon, D, NoM, IECon);

					// elasto-plastic matrix
					DpMat(iElem, kPoint, D, SG, epState, misesStress, 
						ECon, NoM, IECon);

					// save D mat of quad. point
					DMMat(D, currentQuadPoint, DM, NMTE);
					
					// setup shape function N, dNdXi, dNdEta, dNdZeta
					integ(xi, eta, zeta, 
						N, dNdXi, dNdEta, dNdZeta,
						NIf, DXIf, DETIf, DZEIf,
						NMTE, NN, NKind);

					// Jacobian matrix and its inverse
					JMat(x, y, z, dNdXi, dNdEta, dNdZeta, DXIf, DETIf, DZEIf,
						JM, JInv, detJ, 
						NN, NMTE, NKind);

					// compute local B matrix
					BMat(JInv, B, N, dNdXi, dNdEta, dNdZeta, NMTE);
					// save B matrix
					BMMat(B, currentQuadPoint, BM, NMTE);

					// C <- BT.D.B
					BtDBMat(B, C, D, NMTE, BT, BTD);

					// add to element-stiffness matrix EK
					gauss2(C, EK, hNum, N, detJ, NMTE, NN);
				} // 100 end of loop quad. points in current element

				// element-stiffness -> global-stiffness
				merge(EK, NN, GK, NMTE);

				// clear element-stiffness 
				clear2(EK, NMTE);
			} // 150 end of loop elements

			//// dump global K
			//for(int i=1; i<=NTF; i++) {
			//	for(int j=1; j<=NB; j++) {
			//		fprintf_s(fpResult, "%lf\n", GK(i,j));
			//	}
			//}
			//goto cleanup;

			// boundary condition
			bound(IFor, FB, F, IDis, UB, GK);

			/**
			 * start solver
			 */

			// cholesky decomp.
			solve_cholesky_decomp(GK, GKC, F, U, DR, ZM);

			// solve displacement
			solve_cholesky_disp(GKC, F, U, DR, ZM);

			// strain B*a
			strain(NN, U, EP, BM, NMTE, NGNG, EPDash);

			// stress D*strain
			//// TODO D belongs to each element, not universal
			//stress(D, EPDash, SG, SGDASH, NGNG);
			stress(DM, EPDash, SG, SGDASH, NGNG);

			// update displacement and loading
			udata(U, UNow);
			FBNow += FB;


			if(iStep % 100 == 0) {
				printf("Increment Step %d in %d total\n", iStep, nStep);
			}

			if(1) {
				const int totalOutputCount = 100;
				int outputInterval = std::max(nStep / totalOutputCount, 1);
				if(iStep % outputInterval == 0) {
					// average output
					output_av(UNow, EP, SG, 
						NGNG, NMTE, NEUS, NNUS, NN, 
						SGN, EPN, NSGAV);

					char filename[128];
					sprintf_s(filename, "output/result%06d.vtu", iStep/outputInterval);
					output_vtu(filename, iStep, 
						x, y, z, UNow, EP, SG, 
						NMTE, NN, NGNG, NEUS, NNUS, 
						SGN, EPN, NSGAV, 
						IFor, FBNow,
						epState);
				}
			}
		} // end of load incr.

		// final displacement
		utotal(U, UNow, nStep);

		// average output
		output_av(U, EP, SG, 
			NGNG, NMTE, NEUS, NNUS, NN, 
			SGN, EPN, NSGAV);
		// output point values
		output(fpResult, x, y, z, U, EP, SG,
			NGNG, NEUS, NNUS, 
			SGN, EPN, NSGAV);
		output_vtu("result.vtu", nStep, x, y, z, U, EP, SG, 
			NMTE, NN, NGNG, NEUS, NNUS, 
			SGN, EPN, NSGAV, 
			IFor, FBNow,
			epState);

		// dump some instant values
		for(int kk=1; kk<=IRD; kk++) { // element 1
			printf("STR UPPER: %d %lf\n", kk, SG(1,1,kk));
		}
		for(int kk=1; kk<=IRD; kk++) { // element 28
			printf("STR BOTTOM: %d %lf\n", kk, SG(1,28,kk));
		}
		printf("DIS UPPER:%le %le\n", UNow(30), UNow(60));
		printf("DIS BOTTOM:%le %le\n", UNow(270), UNow(300));
	}

cleanup:
	fclose(fpDatacheck);
	fclose(fpResult);

	if(0) {
		//NEQB = 0;
		//std::cout << NEQB << std::endl;

		//NEQ = 0;
		//std::cout << NEQ << std::endl;

		testLinearAlgebra();

		const char *line = "1 205000 0.3 300";
		double young, poisson, yield;
		sscanf_s(line, "%*s %lf %lf %lf", &young, &poisson, &yield);
		std::cout << young <<" " << poisson << " " << yield << std::endl;

	}

	return 0;
}




