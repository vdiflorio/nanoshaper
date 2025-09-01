#include "BlobbySurface.h"


void BlobbySurface::init()
{
	B = DEFAULT_BLOBBYNESS;
	cutoff = DEFAULT_CUTOFF;
	surfType = MOLECULAR_SURFACE;
	providesAnalyticalNormals = false;
}

void BlobbySurface::init(ConfigFile* cf)
{
	Surface::init(cf);
	double blobby_B = cf->read<double>( "Blobbyness", -2.5 );
	// Set up inside value
	inside = 5;
	setBlobbyness(blobby_B);	
}
BlobbySurface::BlobbySurface()
{
	init();
}

BlobbySurface::BlobbySurface(DelPhiShared* ds):MeshSurface(ds)
{
	init();
}


BlobbySurface::BlobbySurface(ConfigFile* cf,DelPhiShared* ds):MeshSurface(ds)
{	
	init();
	Surface::init(cf);
	init(cf);	
}

BlobbySurface::~BlobbySurface()
{
	clear();
}

void BlobbySurface::clear()
{}

void BlobbySurface::setBlobbyness(double b)
{		
	if (b>=0)
	{
		cout << endl << WARN << "Blobbyness is always a stricly negative real: setting "<< DEFAULT_BLOBBYNESS;
		b = DEFAULT_BLOBBYNESS;
	}
	B=b;
}

double BlobbySurface::getBlobbyness()
{
	return B;
}


bool BlobbySurface::build()
{
	// build blobby in an atom-wise way in order to use a cutoff distance

	// numgrid points that are neighbour of an atom 
	int numgrid = (int)(cutoff/delphi->side+0.5);

	// disable cutoff
	//numgrid = delphi->nx;
	
	cout << endl << INFO << "Using cut-off " << cutoff << " num neighbour grid points " << numgrid;

	if (scalarField!=NULL)
		deleteMatrix3D<double>(last_nx,last_ny,scalarField);
	
	scalarField = allocateMatrix3D<double>(delphi->nx,delphi->ny,delphi->nz);

	if (scalarField==NULL)
	{
		cout << endl << ERR << "Cannot allocate scalar field!";
		exit(-1);
	}

	for (int i=0;i<delphi->nx;i++)
		for (int j=0;j<delphi->ny;j++)
			for (int k=0;k<delphi->nz;k++)
				scalarField[i][j][k]=0.0;

	//Atom** atoms = delphi->atoms;	
    std::vector<Atom>& atoms  = delphi->atoms;
    
	printf("\n");
	//int na = delphi->numAtoms;
    unsigned na = atoms.size();
    
	// fill surface class structures
	for (int i=0;i<na;i++)
	{
		printf("\r%sBlobby %.2f%%        ",INFO,((float)i+1)/na*100.0);
		//double* pos = atoms[i]->pos;
        double* pos = atoms[i].pos;
		//double r = atoms[i]->radius;
        double r = atoms[i].radius;
		double r2 = r*r;
		// get ref grid point
		int ix = (int)rintp((pos[0]-delphi->xmin)/delphi->side);
		int iy = (int)rintp((pos[1]-delphi->ymin)/delphi->side);
		int iz = (int)rintp((pos[2]-delphi->zmin)/delphi->side);

		int start_x = MAX(0,(ix-numgrid));
		int start_y = MAX(0,(iy-numgrid));
		int start_z = MAX(0,(iz-numgrid));

		int end_x = MIN((delphi->nx),(ix+numgrid));
		int end_y = MIN((delphi->ny),(iy+numgrid));
		int end_z = MIN((delphi->nz),(iz+numgrid));

		for (int ii=start_x;ii<end_x;ii++)
			for (int jj=start_y;jj<end_y;jj++)
				for (int kk=start_z;kk<end_z;kk++)
				{
					double dist2 = 0;
					double p[3];
					p[0] = delphi->x[ii]-delphi->hside;
					p[1] = delphi->y[jj]-delphi->hside;
					p[2] = delphi->z[kk]-delphi->hside;
					DIST2(dist2,p,pos);
					scalarField[ii][jj][kk] += exp(B*(dist2/r2-1)); 	     	
				}
	}	
	
	bool old = accurateTriangulation;
	bool old_saveMSMS = saveMSMS;
	bool old_vertexAtomsMapFlag = vertexAtomsMapFlag;	
	bool old_computeNormals = computeNormals;

	accurateTriangulation = true;	
	isAvailableScalarField = true;
	
	setSaveMSMS(false);
	setVertexAtomsMap(false);
	setComputeNormals(false);

	int octree_side_size = 2;
	bool found_label = 0;

	while (!found_label)
	{
		if (octree_side_size < delphi->nx ||
			octree_side_size < delphi->ny ||
			octree_side_size < delphi->nz)
		{
			octree_side_size *= 2;
		}
		else
		{
			found_label = 1;
		}
	}
	// allocate empty 
	allocIntersectionsMatrices(octree_side_size);

	/*
	if (computeNormals)
		allocNormalsMatrices();
	*/

	char buff[BUFLEN];
	//sprintf(buff,"%sblobby",conf.rootFile.c_str());
	sprintf(buff,"blobby");

	triangulateSurface(1.0,buff);

	cout << endl << INFO << "Smoothing blobby surface...";
	smoothSurface(buff);
	
	setSaveMSMS(old_saveMSMS);
	setVertexAtomsMap(old_vertexAtomsMapFlag);	
	setComputeNormals(old_computeNormals);

	accurateTriangulation = old;

	isAvailableScalarField = false;
	
	cout << "ok!";

	deleteMatrix3D<double>(delphi->nx,delphi->ny,scalarField);
	
	sprintf(buff,"%sblobby.off",conf.rootFile.c_str());

	// the triangulation tools of surface were used. Now they are cleaned

	////////////// reset surface status ////////////////////////////
	if (triList.size()>0)
	{
		vector<int*>::iterator it;
		for (it = triList.begin();it!=triList.end();it++)
			deleteVector<int>((*it));
	}

	if (vertList.size()>0) 
	{
		vector<double*>::iterator it;
		for (it = vertList.begin();it!=vertList.end();it++)
			deleteVector<double>((*it));
	}

	triList.clear();
	vertList.clear();
	totalSurfaceArea=0;
	totalVolume=0;	
	//////////////////////////////////////////////////////////

	// load its mesh
	load(buff);

	// the rest continue as a mesh surface
	MeshSurface::printSummary();
	return true;
}

void  BlobbySurface::printSummary()
{
	cout << endl << INFO << "Blobbyness value " << getBlobbyness();
	cout << endl << INFO << "Cut-off distance " << cutoff << " [A]";
}
