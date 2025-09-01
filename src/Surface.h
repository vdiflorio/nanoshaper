
//---------------------------------------------------------
/**    @file		Surface.h
*     @brief	Surface.h is the header for CLASS
*               Surface.cpp								*/
//---------------------------------------------------------

#ifndef Surface_h
#define Surface_h

#include "octree.h"
#include "globals.h"
#include "ConfigFile.h"
#include "SurfaceFactory.h"
#include <stack>  

#ifdef ENABLE_CGAL 
	//////////////////////// CGAL ///////////////////////////////////////////////////
	#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	#include <CGAL/Regular_triangulation_3.h>
	//#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
    //<CGAL/Regular_triangulation_euclidean_traits_3.h
    #include <CGAL/config.h>    
	//#include <CGAL/Regular_triangulation_filtered_traits_3.h>
	#include <CGAL/Triangulation_cell_base_with_info_3.h>
	#include <CGAL/Regular_triangulation_cell_base_3.h>
	#include <CGAL/Triangulation_vertex_base_with_info_3.h>
	#include <CGAL/Triangulation_data_structure_3.h>	
	////////////////////////////////////////////////////////////////////////////////

	// vor point container
	class VorPoint
	{
	public:
		double vor[3]; // Voronoi center
		bool visited;

		VorPoint()
		{
			visited = false;
		}
	};
#endif


/** grid coordinates + vector (usually a point or a normal vector)*/
class coordVec
{
public:

	coordVec(const int x,const int y,const int z,double* const v,const int d)
	{	
		ix = x;
		iy = y;
		iz = z;
		vec = v;
		dir = d;
	}

	/** shallow copy constructor*/
	coordVec(const coordVec& cv)
	{
		ix = cv.ix;
		iy = cv.iy;
		iz = cv.iz;
		vec = cv.vec;	
		dir = cv.dir;
	}

	int ix;
	int iy;
	int iz;
	int dir;
	double* vec;
};

class packet
{
public:
	packet()
	{
     first = nullptr;
     second = nullptr;
    }
	packet(const packet& p)
	{
		first = p.first;
		second = p.second;
	}
	vector<coordVec>* first;
	vector<coordVec>* second;
};


#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif

#include "tools.h"
#include "DelphiShared.h"

// ids for rays directions
#define X_DIR 0
#define Y_DIR 1
#define Z_DIR 2

// ids for panels
#define PANEL_YZ 0
#define PANEL_XY 1
#define PANEL_XZ 2

extern int num_cores;

#define DEFAULT_VOLUME 11.4
// tollerance on two intersections. If less than EPS_INT then the intersection is the same. def 1e-8
#define EPS_INT 1e-8
//#define MAX_ATOMS_MULTI_GRID 100
//#define DEBUG_SURFACE

#define DELTA 1e-10
#define RAND_DISPLACEMENT 1e-4

// boundary grid point codes
#define INTERNAL_BGP 0
#define EXTERNAL_BGP 1

// molecular surface
#define MOLECULAR_SURFACE 0
// analytical object
#define OBJECT_SURFACE 3
// unknown surface, e.g., mesh
#define GENERIC_SURFACE 2
// 'sum' of two different surface types
#define HYBRID_SURFACE 4

// file format
#define OFF 0 
#define OFF_A 1 
#define OFF_N 2
#define OFF_N_A 3 
#define MSMS_NO_A 4
#define MSMS 5

// now using inline templates
//#define GRID_MULTI_MAP(i,j,k,l,NX,NY,NZ) gridMultiMap[(l)+(MAX_ATOMS_MULTI_GRID)*((k)+(NZ)*((j)+(NY)*(i)))]

#ifdef ENABLE_BOOST_THREADS

#define THREAD_SAFE_SCOPE(SSS) \
	{ \
		boost::mutex::scoped_lock scopedLock(mutex); \
		SSS \
		cout.flush(); \
	}

#else

#define THREAD_SAFE_SCOPE(SSS) \
	{ \
		SSS \
		cout.flush(); \
	}
#endif




/** @brief Surface class is the general interface that a surface class should have
to be plugged inside DelPhi. Some functions implementations are mandatory
such as load, save, build, etc... Note that the surface is not necessarly a molecular surface.
Build function computes an internal representation of the surface; getSurf translates
that representation in the DelPhi compatible representation. 
Note that in order to put a new surface in DelPhi a surface must provide epsmap,
idebmap, computations of the surface area inside a grid cube, identification and projections
of boundary grid points and their surface normals; these computations must be done in getSurf while the surface construction
must be performed in build. \n \n

The global variables in the DelPhiShared object that must be filled are: \n
idebmap -> salt (in/out) map \n
epsmap  -> dielectric map \n
ibgp    -> indexes of the boundary grid points \n
scspos  -> coordinates of the projected boundary grid points \n
scsnor  -> coordinates of the outside pointing normal vector over the surface in correspondence of boundary grid points \n
scsarea -> value of the surface area in each cube of the grid where there is a boundary grid point (optional)\n

If a new surface is defined and it only provides the routines of ray intersection and projection
then most information can be still built; in this case one should not overload the getSurf method 
that in its base implementation is able to use projections/intersections to build almost all info needed by DelPhi to work; however this mode
of operating has restrictions in fact only two dielectrics can be used (in/out) and the Stern layer (idebmap)
is only in/out; this mode is used in the MeshSurface class where an arbitrary shape is loaded or by the SkinSurface module. Note 
that there is no need for the surface (or surface patches) to answer in/out queries; indeed in/out is worked
out by counting the number of times the ray intersects the surface starting from the outside.

To get a full implementation one has to derive the Surface class and re-implement the interface method
load,save,build and getSurf (or ray/projections routine) thus providing all the necessary info to DelPhi solver.


@author Sergio Decherchi 
@date 29/06/2013 
*/
class Surface
{
#ifdef ENABLE_CGAL 
private:
	//AV
	typedef CGAL::Exact_predicates_inexact_constructions_kernel _K;
 	// Regular T3


        typedef CGAL::Triangulation_cell_base_with_info_3< VorPoint*,_K> _Cell_with_info;
        typedef CGAL::Regular_triangulation_cell_base_3<_K, _Cell_with_info> _RT_Cell_with_info;
        
        // nopatch
        typedef CGAL::Triangulation_vertex_base_with_info_3<int,_K> _Vertex_with_info;
        typedef CGAL::Regular_triangulation_vertex_base_3<_K, _Vertex_with_info> _RT_Vertex_with_info;

        typedef CGAL::Triangulation_data_structure_3<
        _RT_Vertex_with_info,
        _RT_Cell_with_info,
        CGAL::Parallel_tag>                                           _tds;

        typedef CGAL::Regular_triangulation_3<_K, _tds>               _Rt;
        typedef _Rt::Weighted_point                                   _Weighted_point;
        typedef _Rt::Vertex_handle                                    _Vertex_handle;

        typedef _Rt::Bare_point                                       _Point3;
        typedef _Rt::Vertex_iterator                                  _Vertex_iterator;
        typedef _Rt::Finite_vertices_iterator                         _Finite_Vertex_Iterator;
        typedef _Rt::Finite_cells_iterator                            _Finite_Cells_Iterator;
        typedef _tds::Cell_circulator                                 _Cell_circulator;
        typedef _tds::Cell_handle                                     _Cell_handle;

        typedef _Rt::Finite_edges_iterator                            _Finite_Edges_Iterator;
        typedef _Rt::Finite_facets_iterator                           _Finite_Facets_Iterator;
        typedef _Rt::Facet                                            _Facet;

        typedef _K::FT                                                _Weight;


/*
        typedef CGAL::Triangulation_cell_base_with_info_3< VorPoint*,_K> _Cell_with_info;
        typedef CGAL::Regular_triangulation_cell_base_3<_K,_Cell_with_info> _RT_Cell_with_info;
        typedef CGAL::Triangulation_vertex_base_3<_K> _Vertex;
        //typedef CGAL::Triangulation_data_structure_3<_Vertex,_RT_Cell_with_info> _tds;

        typedef CGAL::Triangulation_data_structure_3<
            CGAL::Regular_triangulation_vertex_base_3<K>,
            CGAL::Regular_triangulation_cell_base_3<K>,
            CGAL::Parallel_tag>                                          _tds;
                                                                 

 	//typedef CGAL::Triangulation_data_structure_3<
	//CGAL::Regular_triangulation_vertex_base_3<_K>,
        //CGAL::Regular_triangulation_cell_base_3<_K>,
        //CGAL::Parallel_tag>                  	                     _tds;
        typedef CGAL::Regular_triangulation_3<_K, _tds>	             _Rt;
  	typedef _Rt::Bare_point                                      _Bare_point;
        typedef _Rt::Weighted_point                                  _Weighted_point;
        typedef _Rt::Vertex_handle                                   _Vertex_handle;
  
        typedef _Rt::Bare_point    	                             _Point3;
        typedef _Rt::Vertex_iterator                                 _Vertex_iterator;
        typedef _Rt::Finite_vertices_iterator                        _Finite_Vertex_Iterator;
        typedef _Rt::Finite_cells_iterator                           _Finite_Cells_Iterator;
        typedef _tds::Cell_circulator                                _Cell_circulator;
        typedef _tds::Cell_handle                                    _Cell_handle;

	typedef _K::FT                                               _Weight;
*/
	/*
	typedef CGAL::Exact_predicates_inexact_constructions_kernel  _K;
	typedef CGAL::Regular_triangulation_filtered_traits_3<_K>    _Traits;
        typedef CGAL::Regular_triangulation_euclidean_traits_3<_K>    _Traits;
	typedef CGAL::Triangulation_cell_base_with_info_3< VorPoint*,_Traits> _Cell_with_info;
	typedef CGAL::Regular_triangulation_cell_base_3<_Traits,_Cell_with_info> _RT_Cell_with_info;
	typedef CGAL::Triangulation_vertex_base_3<_Traits> _Vertex;
	typedef CGAL::Triangulation_data_structure_3<_Vertex,_RT_Cell_with_info> _tds;	
	
	typedef _Traits::RT                                          _Weight;
	typedef _Traits::Bare_point                                  _Point3;
	typedef _Traits::Weighted_point                              _Weighted_point;
	
	typedef CGAL::Regular_triangulation_3<_Traits,_tds>          _Rt;
	typedef _Rt::Vertex_iterator                                 _Vertex_iterator;
	typedef _Rt::Finite_vertices_iterator                        _Finite_Vertex_Iterator;
	typedef _Rt::Finite_cells_iterator			     _Finite_Cells_Iterator;
	typedef _Rt::Vertex_handle                                   _Vertex_handle;
	typedef _tds::Cell_circulator				     _Cell_circulator;
	typedef _tds::Cell_handle 				     _Cell_handle;
	*/

#endif

public:
	// current panel under analysis
	int panel;

protected:

	//////////////////////////////////// surf definition variables ///////////////////////
	/** says if the surface represents a molecule, an hybrid system or an object */
	int surfType;
	/** every class at Constructor time must say if it is Ray Casting based or not
	 If it is ray casting based than getSurf will derive the grid based on ray casting
	 if not RC based getSurf will assume that build method has already coloured the grid */
	bool isRCbased;
	/** every surface at constructor time must say if it will provide analytical normals during
	ray tracing or not*/
	bool providesAnalyticalNormals;
	////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////// flags/settings ////////////////////////////////////
	bool projBGP;
	double sternLayer;
	bool accurateTriangulation;
	bool doTriangulation;
	bool fillCavitiesFlag;
	bool computeNormals;
	/** if enabled loaded surface is checked for duplicated vertices*/
	bool checkDuplicatedVertices;
	/** if enabled the part of the cavities (or the entire cavities) that are not shaped as a water molecule are
	removed*/
	bool wellShaped;
	/** Probe radius. Currently used in the cavity shape filter*/
	double probe_radius;	
	bool useLoadBalancing;
	bool vertexAtomsMapFlag;
	bool saveMSMS;

	#ifdef LIMIT_NUM_ATOMS
	int maxNumAtoms;
	double domainShrinkage;
	#endif
	//////////////////////////////////////////////////////////////////////////////////////

	// delphi environment
	DelPhiShared* delphi;	
	double totalSurfaceArea;
	double totalVolume;
	
	// last row and column sizes of ind_2d matrix for ray casting acceleration
	int64_t last_rows_ind;
	int64_t last_cols_ind;
	/** how big is the random initial displacement of atoms*/
	double randDisplacement;
	// last nx,ny,nz dimensions seen by Surface class
	int64_t last_nx,last_ny,last_nz;
	/** 3d matrix of intersections along x rays*/
	Octree<int>* intersectionsMatrixAlongX;
	/** 3d matrix of intersections along y rays*/
	Octree<int>* intersectionsMatrixAlongY;
	/** 3d matrix of intersections along z rays*/
	Octree<int>* intersectionsMatrixAlongZ;
	/** 3d matrix of normals along x rays*/
	Octree<int>* normalsMatrixAlongX;
	/** 3d matrix of normals along y rays*/
	Octree<int>* normalsMatrixAlongY;
	/** 3d matrix of normals along z rays*/
	Octree<int>* normalsMatrixAlongZ;
	
	/** mark wich MC cubes contain triangles to allow fast reject in the second pass*/	
	bool* activeCubes;

	bool*** verticesInsidenessMap;
	//bool *verticesInsidenessMap;

	double*** scalarField;
	
	// when a scalar field is available scalarField is used instead of verticesInsidenessMap for vertex interpolation
	bool isAvailableScalarField;
	
	////////////////////////////// triangulation data structures ////////////////////////////////
	/** vector of vertex indices for the traingulation obtained by triangulateSurface function*/
	vector<int*> triList;
	/** vector of triangle vertices*/
	vector<double*> vertList;
	/** vector of normal to vertices*/
	vector<double*> normalsList;
	/////////////////////////////////////////////////////////////////////////////////////////////

	double delta_accurate_triangulation;
	
	/** type of bgp for each detected bgp **/
	int* bgp_type;
	
	/** grid multi-dielectric map*/
	int* gridMultiMap;
	double gxmin,gymin,gzmin,gside,gscale;
	int64_t ggrid;
	
	/** pointer to a vector which holds the information about the load for each grid slice*/
	int* gridLoad;
	int totalLoad;
	int* vertexAtomsMap;

	int MAX_ATOMS_MULTI_GRID;

	/** DelPhi code for inside.*/
	int inside;

	#ifdef ENABLE_BOOST_THREADS
		boost::mutex mutex;
	#endif

	// added keyword to support saving of files 
	string rootFile;

	/////////////////////////////////////////////////////////////////////////////////////////////////

	/** 3d cavity detection*/
	void floodFill(int ix,int iy,int iz,int idold,int idnew);
	
	/** scan line version*/
	void floodFill2(int ix,int iy,int iz,int idold,int idnew);
	
	/** parallel version with scaline */
	void floodFill4(int ix,int iy,int iz,int idold,int idnew);
	
	/** inner routine for scanline */
	void floodFill3(	pair<pair<int,int>,int > ind, 
						pair<int,int> z_limits, 
						pair<int,int> old_new, 
						pair< queue< pair<pair<int,int>,int> > *, queue< pair<pair<int,int>,int > > * > queues,
						queue< pair<pair<int,int>,int> > * in_queue);
	
	/** supports inner floodfill inside the 1.4 surface (given status1) to 
	understand if two cavities/pockets are status2-linked (check external communication
	between two points by looking at status2 map). Gives true if the cavities/pockets comunicate
	and uses as input two random indices of two cavities/pockets. It stops if a maximal number of
	moves have been reached. A 'move' means moving from one grid point to another one. In max moves is returned
	the number of moves done to get the target.*/
	bool innerFloodFill(int* start,int* target,int* status1,int* status2,int& maxMoves,bool debug);
	
	/** gives true if the point is outside vdw surface*/
	bool vdwAccessible(double* p,int& nearest);
	
	/**intersector routine, used to perform partial or full intersections. Usefully
	used together with boost threading routines. In order to get a 'robust' ray tracer a particular strategy is adopted during ray-tracing. 
	It can happen that, due to numerical imprecisions
	of the ray-intersection routine or due to problems of the surface (e.g. degenerate triangles, holes, non manifoldness in general) 
	the number of detected intersections is odd. In this case the ray-tracer randomly perturbs the direction of the ray in order to escape
	from the singularity. If this strategy fails for a given number of trails, the trace of the previous ray is copied to the current one. \n
	In some cases this strategy can also fill small holes in the surface.
    SD PB_NEW: added now the package to collect grid and normals intersections.
    */
	void intersector(double* volPanel,int nb,int start,int end,int jump,int* numIntersections,packet pack,packet gridPack=packet());
	
	/** projector routine, used to perform partial or full intersections. Usefully
	used together with boost threading routines*/
	void projector(int start,int end);
	
	/**return +1 if outside and -1 if inside for the vertex indexed by vertInd belonging to the grid cube
	given by i,j,k indexes*/
	char getInsidness(int i, int j, int k, int vertInd);
	
	/** marching cubes vertex interpolation*/
	void vertexInterp(double isolevel, double* p1, double* p2, double valp1, double valp2,double* p);
	
	/** build triangles for marching cubes cell*/
	//int getTriangles(double* vertexValues,double** vertexPos,double isolevel, int** triangles,int ix,int iy,int iz,
	//	int NX, int NY,int NZ,int* xedge,int* yedge,int* zedge,int* xedge_down,int* yedge_down);
	int getTriangles(double* vertexValues,double** vertexPos,double isolevel, int** triangles,int ix,int iy,int iz,
		int NX, int NY,int NZ);

	/** returns the vertices for a given section on z of the grid*/
	void getVertices(double isolevel,int start_z,int end_z,int jump,vector<coordVec>*,vector<double*>*);

	/** approximate normals based on the surrounding triangle planes normals*/
	void approximateNormals(vector<int>& appNormals,bool doOnlyList);
	
	/** triangulator thread*/	
	double triangulationKernel(double isolevel,bool revert,int start_z,int end_z,int jump,vector<int*>* localTriList,double* localArea);
	
	/** build a 3D grid for acellerating nearest atom queries*/
	void buildAtomsMap();
	
	/** deallocate memory of the 3D nearest atom query */
	void disposeAtomsMap();
	
	/** apply multidielectric correction after grid building. For each internal grid point detect the nearest
	atom and apply its dielectric constant.*/
	void applyMultidielectric();
	
	/** swap the state of a point in the epsmap from internal to the nearest atom dielectric*/
	void swap2multi(double gxmin,double gymin,double gzmin,double gside,unsigned int ggrid,int* gridMultiMap,int i,int j,int k,int l);
	
	/** build stern layer.*/
	void buildSternLayer();
	
	/** clean and alloc intersections*/
	void allocIntersectionsMatrices(int octree_side_size);
	
	/** clean and alloc normals*/
	void allocNormalsMatrices(int octree_side_size);
	
	/** based on provided info decides the best format to save the mesh*/
	int deduceFormat();

	/** gives true if the point is completely out: completely
	out means that the nearest grid point is out all its 1-neighbours are out*/
	bool isCompletelyOut(double* pos);
		
protected:
	
	Surface();
	Surface(ConfigFile* cf);

public:

	// In order to use the set of default surface methods, the following
	// interface methods must be provided.
	// If a fully custom solution is built the mandatory methods can be fake methods
	// and all the work can be made in getSurf();

	//////////////////////// INTERFACE MANDATORY METHODS //////////////////////////////////
	
	/** Build the surface internal representation. Returns if the
	building process was succesful or not*/
	virtual bool build() = 0;
	
	/** Save the surface in a specific class dependent format. Return
	if saving was succesful or not*/
	virtual bool save(char* fileName) = 0;
	
	/**Load the surface in a specific class dependent format. Return
	if loading was succesful or not*/
	virtual bool load(char* fileName) = 0;
	
	/** Print a summary of the surface type, status and other stuff*/
	virtual void printSummary() = 0;				
	
	/** Get a projection of a point on the surface. Return projection and normal*/	
	virtual bool getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3)=0;
	
	/** Get all the intersections of a ray that goes from P1 to P2 over the surface.
	The interesctions are returned with increasing distance order. 
	the double in the vector is the t parameter for the intersection of the parametric 
	line and the surface, the double pointer is the normal vector. The management of the memory
	of the normal is up to the derived class from surface*/
	virtual void getRayIntersection(double p1[3],double p2[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals)=0;

	/** function for the constructor without arguments*/
	virtual void init();

	/** functions for the constructor with config file argument*/
	virtual void init(ConfigFile* cf);

	/**function for the denstructor*/
	virtual void clear();
	/////////////////////////////////////////////////////////////

	//////////////////////// PROVIDED DEFAULT METHODS /////////////////////////////////
	/** Build DelPhi Surface by a specific method. If the defined surface
	provides at least projection and intersection routines then don't overload
	this function; it will already provide most of the information needed by DelPhi.
	The only restrictions are that only two dielectrics are allowed (in/out) and 
	there is no Stern layer. Overload this method if the surface does not provide
	intersect or project or all information is needed. This mode of 
	operation is used by the MeshSurface class which provides getRayIntersection and 
	getProjection primitives. If requested one can fill cavities by fill flag. In order to
	use parallel execution ray-tracing is performed using the itersectionFunctor.
    
    SD PB_NEW: added the intersectionsInfo input. If different from nullptr than the intersections and normals
    data is loaded into this vector
    */
	virtual bool getSurf(bool fill=false,double vol=0,vector<packet>* intersectionsInfo=nullptr);		
	
	/** Returns the volume computed during Surface::getSurf computations. This function
	can be overload to implement a custom volume computation method*/
	virtual double getVolume();
	
	/** Return the total surface area*/
	virtual double getArea();
	
	/** Compute the cavities of the surface by grid flooding. A list of cavities is returned 
	where each list is a list of triplet of indexes. It is up to the caller to free the memory
	of the int* pointers. The input decide wich is the first STATUS code to be checked.*/
	virtual int getCavities(int idStart=STATUS_POINT_TEMPORARY_OUT);
	
	/** Fill the previously detected cavities if their volume is bigger that the passed var. In the 
	baseline implementation the volume is approximated by the number of grid cubes volumes in that cavity.
	By default all cavities are filled*/
	virtual void fillCavities(double vol=0,bool silent=false);
	
	/** After the cavities are detected and marked this routine filter out the part or the entire cavities
	that are not able to fit the bounding box of a water molecule. That is we filter the bad shaped cavities*/
	virtual void filterCavities(bool modStatus=false);
	
	/** swap cavity status to temporary out*/
	virtual void cav2out();
	
	/** Triangulate the surface and save it in OFF format. In the Surface class the baseline method
	for triangulation is obtained by employing the marching cube method at each delphi grid cell.
	For each vertex in the grid cube its insideness is computed by voting of the insidness values of
	the incident cubes; thus the scalar field is given by the ensemble of the status map. The surface is saved in
	off format*/
	virtual double triangulateSurface(double iso=0.0,const char* fileName="triangulatedSurf",bool revert=false);	
	
	/** save mesh in a prescribed format, revert triangles (change plane sign) if requested*/
	virtual bool saveMesh(int format,bool revert,const char* fileName,vector<double*>& vertList,vector<int*>& triList,vector<double*>& normalsList);
	
	/** smooth a given mesh and overwrites the given file name.
	The input/output mesh is in .off format*/
	virtual void smoothSurface(const char* fn="triangulatedSurf",bool revert=false);

	/** this function is called before the ray tracing of the panel. It can be useful if a per panel pre-processing step is needed.
	By default this function does nothing*/
	virtual void preProcessPanel()
	{}
	
	/** this function is called after ray casting to perform any post processing, such as memory clean-up, after
	ray casting*/
	virtual void postRayCasting()
	{}
	
	/** this function is called before boundary grid projection to perform any pre-processing, such as memory setup-up.*/
	virtual bool preBoundaryProjection()
	{		
		return true;
	}

	/** translate a triangulation (a set of samples) into a set of balls (ray-casting sampling + power crust)*/
	void tri2Balls();
	
	/** save current status map in a temporary map*/
	void backupStatus();
	
	/** remove temporary status*/
	void removeBackupStatus();
	
	/** If between two any cavities there is communication in map st2 and this communication is all internal in map st1 these cavities are merged logically.*/
	int linkCavities(int* st1,int* st2);
	
	/** difference operator: said S1 the first surface (big probe), S2 (small probe) the second surface and S3 the
	output surface then the difference rule is the following rule:
	if S1 is not OUT and and S2 is OUT then that's the pocket
	*/	
	bool difference(Surface*);

	/** This is a Connolly regularized difference operator. Difference function is called
	and then Connolly filter is applied such that the noise can be filtered out.*/
	Surface& operator-=(Surface& surf2);

	/** for each cavity/pocket detect the atoms that compose the cavity*/
	void getCavitiesAtoms();

	/** given a triangulation coming from the current surface object says if the triangulations
	points are completely out or not according to the status map given in surf*/
	void triangulationPointsAreCompletelyOut(Surface* surf,vector<bool>& results);

	/** save selected points of the triangulation in file without normals, in filen with normals if they
	are available*/
	void saveSelectedPoints(vector<bool>& selection,char* file,char* filen);

	/**save a subset of the triangulation into triSubset file. The subset is given in a vertex wise way
	by the vector of boolean results. Revert must be true if the mesh comes from the
	triangulation of a pocket or a cavity. The value of the new area is returned.*/
	double saveTriSubSet(char* triSubset,vector<bool>& results,bool revert);

	/** gives MC cube index. If !=-1 then that cube contains triangles*/
	int classifyCube(double* vertexValues,double isolevel);

	///////////////////////////////////////////// getters/setters //////////////////////////////////////////////

	/** Set if surface has to project bgp or not. During cavity detection/surface visualization/triangulation
	bgp projection can be skipped*/
	void setProjBGP(bool flag)
	{
		projBGP = flag;
	}
	bool getProjBGP()
	{
		return projBGP;
	}
	double getSternLayer()
	{
		return sternLayer;
	}
	void setSternLayer(double l)
	{
		if (l<0)
		{
			cout << endl << WARN << "Cannot set a negative Stern Layer";
			return;
		}
		else
			sternLayer = l;
	}

	void setSaveMSMS(bool m)
	{
		saveMSMS = m;
	}

	bool getSaveMSMS(void)
	{
		return saveMSMS;
	}

	void setAccurateTriangulationFlag(bool flag)
	{
		accurateTriangulation = flag;
	}

	bool getAccurateTriangulationFlag(void)
	{
		return accurateTriangulation;
	}

	void setTriangulationFlag(bool flag)
	{
		doTriangulation = flag;
	}

	bool getTriangulationFlag(void)
	{
		return doTriangulation;
	}

	void setVertexAtomsMap(bool f)
	{
		vertexAtomsMapFlag = f;
	}

	bool getVertexAtomsMap()
	{
		return vertexAtomsMapFlag;
	}

	void setComputeNormals(bool cn)
	{
		computeNormals = cn;
	}

	bool getComputeNormals()
	{
		return computeNormals;
	}

	void setCheckDuplicatedVertices(bool cd)
	{
		checkDuplicatedVertices = cd;
	}
	bool getCheckDuplicatedVertices()
	{
		return checkDuplicatedVertices;
	}
	void setKeepWellShapedCavities(bool kwsc)
	{
		wellShaped = kwsc;
	}
	bool getKeepWellShapedCavities()
	{
		return wellShaped;
	}
	virtual void setProbeRadius(double probeRadius)
	{
		probe_radius = probeRadius;
	}
	virtual double getProbeRadius()
	{ 
		return probe_radius;
	}

	virtual double getRandDisplacement()
	{
		return randDisplacement;
	}

	virtual void setRandDisplacement(const double r)
	{		
		randDisplacement = r;
	}

	virtual void setLoadBalancing(bool doLoadBalancing)
	{
		useLoadBalancing = doLoadBalancing;
	}

	virtual bool getLoadBalancing()
	{
		return useLoadBalancing;
	}

	virtual int getNumTriangles()
	{
		return (int)triList.size();
	}

	virtual int getNumVertices()
	{		
		return (int)vertList.size();
	}

	virtual void setInsideCode(int i)
	{
		inside = i;
	}

	virtual int getInsideCode()
	{
		return inside;
	}

	#ifdef LIMIT_NUM_ATOMS
	void setMaxNumAtoms (int max_atoms)
	{
		maxNumAtoms = max_atoms;
	}

	void setDomainShrinkage (double domain_shrinkage)
	{
		domainShrinkage = domain_shrinkage;
	}
	#endif

	///////////////////////////////////////////////////////////////////////////////////////////////////////

	/** Destructor*/
	virtual ~Surface();
};

////////////////////////////////////////////////////////////////////////////////////////

#endif
