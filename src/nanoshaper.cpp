
#include "nanoshaper.h"
#include "DelphiShared.h"
#include "Surface.h"

using namespace NS;

#define SURF (static_cast<Surface*>(surf))
#define CONFIG (static_cast<ConfigFile*>(cf))
#define DS (static_cast<DelPhiShared*>(ds))

#ifdef ENABLE_CGAL 
	#ifdef CGAL_LINKED_WITH_TBB
		# include <tbb/global_control.h>
	#endif
#endif


NanoShaper::NanoShaper( const std::vector<Atom>& ats,
                        surface_type flag,
                        const double _parameter,
                        const double _stern_layer_thickness,
                        const unsigned nt,
                        const std::string* configFile)
{
    initConstructor(ats,flag,_parameter,_stern_layer_thickness,nt,configFile);
}


void NanoShaper::initConstructor( const std::vector<Atom>& ats,
                        surface_type flag,
                        const double _parameter,
                        const double _stern_layer_thickness,
                        const unsigned nt,
                        const std::string* configFile)
{
    cf = nullptr;
    ds = nullptr;
    surf = nullptr;
    collectGridRays = false;

    initConfig(configFile);
    // copy the atoms internally    
    atoms = ats;
    // unknown panel
    currentPanel = -1;
    // override with local info    
    CONFIG->add<int>( "Number_thread",nt);   

    // global variable for Boost and TBB
    num_cores = nt;

    #ifdef ENABLE_CGAL
	 	#ifdef CGAL_LINKED_WITH_TBB
			int TBB_num_thd = 1;
			if (num_cores<=0)
			{
				TBB_num_thd = (int)std::thread::hardware_concurrency();
				cout << endl << INFO << "Detected " << TBB_num_thd<< " logical cores";
				cout << endl << INFO << "TBB - Setting " << TBB_num_thd << " threads";
			}
			else
			{
				cout << endl << INFO << "TBB - User selected num threads " << num_cores;
				TBB_num_thd = num_cores;
			}
			auto mp = tbb::global_control::max_allowed_parallelism;
			tbb::global_control gc(mp, TBB_num_thd);			
		#endif
	#endif

    switch(flag) 
    {        
        case skin:
        {
            CONFIG->add<std::string>("Surface","skin");            
            CONFIG->add<double>( "Skin_Surface_Parameter",_parameter);                
            break;   
        }
        case ses:
        {
            CONFIG->add<std::string>("Surface","ses");
            CONFIG->add<double>( "Probe_Radius",_parameter);    
            break;  
        }
        default: 
            std::cout << ERR << "Unknown surface type. Existing\n"; std::exit(EXIT_FAILURE);
    }
}


/**it returns +1 for outside, -1 for inside, 0 for unknown (grid not computed)*/
int NanoShaper::getGridPointStatus(unsigned ix,unsigned iy,unsigned iz)
{      
   if (ix<0 || iy<0 || iz<0 || ix>=DS->nx || iy>=DS->ny || iz>=DS->nz || ds==nullptr || DS->status==nullptr)
   {    
    return 0;
   }
   // -1 -> point is inside
   // 2 -> point is temporary outside (may be really out or in cavity)
   // 3 -> point is really out	
   int val = read3DVector<int>(DS->status,ix,iy,iz,DS->nx,DS->ny,DS->nz);
   if (val<0)
    return -1;
   else
    return +1;   
}

/**it returns the grid size*/
bool NanoShaper::getGridSize(std::vector<unsigned>& p)
{    
   if (ds==nullptr || DS->status==nullptr)
    return false;
    
   p.clear(); 
   p.push_back(DS->nx);
   p.push_back(DS->ny);
   p.push_back(DS->nz);
   return true;
}

/**it returns false if the point cannot be accessed*/
bool NanoShaper::getGridPointCoordinates(unsigned ix,unsigned iy,unsigned iz,std::vector<double>& p)
{   
   if (ix<0 || iy<0 || iz<0 || ix>=DS->nx || iy>=DS->ny || iz>=DS->nz || ds==nullptr || DS->status==nullptr)
    return false;
   
   p.clear();
   p.push_back(DS->x[ix]);
   p.push_back(DS->y[iy]);
   p.push_back(DS->z[iz]);
   return true;
}

bool NanoShaper::setDirection(unsigned direction)
{
    if (direction>2)
        return false;
        
    if (surf==nullptr)
        return false;
                       
    currentPanel = panelResolver(direction);
    // set the panel
    SURF->panel = currentPanel;
    //every time we change direction we need to prepare
    //the right data structures.
    SURF->preProcessPanel();
    return true;
}


int NanoShaper::panelResolver(unsigned direction)
{
    //x axis
    if (direction == 0)
        return 0;
    //z axis
    else if (direction == 2)
        return 1;        
    //y axis
    else if (direction == 1)
        return 2;
    else
    {
        std::cout << ERR << "Please select a valid direction to cast the ray 0=x,1=y,z=2";
        return -1;
    }    
}

void NanoShaper::pointResolver(double *pb,double coord)
{
    if (currentPanel==0)
        pb[0]=coord;
    else if (currentPanel==1)
        pb[2]=coord;  
    else if (currentPanel==2)
        pb[1]=coord;  
}

short NanoShaper::castAxisOrientedRay(double* pa,double coord,std::vector<std::pair<double,double*> >& intersections,unsigned direction,bool computeNormals)
{
    int thd = 0;
    double pb[3];
    
    int desiredPanel = panelResolver(direction);
    
    if (desiredPanel<0)
    {
        std:cout << ERR << "Wrong requested ray direction";                
        return -1;
    }
        
    if (currentPanel!=desiredPanel)
    {
        std::cout << ERR << "Cannot cast a ray in a direction that is not the current one. Please use setDirection to set current direction before casting rays";                
        return -1;
    }

    // copy pa to pb
    pb[0]=pa[0];
    pb[1]=pa[1];
    pb[2]=pa[2];
    
    // update the correct coordinate based on the current panel
    pointResolver(pb,coord);
    
    double diff[3];
    diff[0]=pb[0]-pa[0];
    diff[1]=pb[1]-pa[1];
    diff[2]=pb[2]-pa[2];
        
    if (pa[direction]>=coord)
    {
        std::cout << ERR << "Please cast rays only in progressive directions";
        return -1;
    }        
    
    
    if (surf==nullptr)
    {
        std::cout << ERR << "Please compute a surface before casting rays";
        return -1;
    }
    
    std::vector<std::pair<double,double*> > intersecs;    
	SURF->getRayIntersection(pa,pb,intersecs,thd,computeNormals);

    if (intersecs.size()==0)
		return -1;               
						    
	bool closure = false;
	double lastValid = INFINITY;	
    std::vector<std::pair<double,double*> >::iterator itt = intersecs.begin();
    
    // performs a cleanup of the intersections
    // too near ones are removed
    for (;itt!=intersecs.end();itt++)
    {
        double a = (*itt).first;
        a = pa[direction]+diff[direction]*a;
        if (fabs(a-lastValid)>=EPS_INT)
        {
            std::pair<double,double*> pp;
            // put the interesection coordinate
            pp.first = a;
            pp.second = (*itt).second;
            // keep the intersection
            intersections.push_back(pp);
            lastValid = a;
        }        
        else
        {
            // discard the intersection and delete the normal if present
            if (computeNormals)
                delete (*itt).second;
        }
    }
    
    return +1;
}

bool NanoShaper::buildAnalyticalSurface(double randDisplacement)
{   
    if (ds!=nullptr)
        delete DS;
 
    if (surf!=nullptr)
        delete SURF;
                
    DelPhiShared* dsl = nullptr;

    if (!CONFIG->keyExists("PB_grid_mode") || !CONFIG->read<bool>("PB_grid_mode"))
    {
        dsl = new DelPhiShared(CONFIG->read<double>("Grid_scale"),
                          CONFIG->read<double>("Grid_perfil"),
                          atoms,
                          CONFIG->read<int>("Max_Num_Atoms"),
                               CONFIG->read<double>("Domain_Shrinkage"),
                               CONFIG->read<bool>("Build_epsilon_maps"),
                               CONFIG->read<bool>("Build_status_map"),
                               CONFIG->read<bool>("Multi_Dielectric"));
    }
    else 
    {

        double coord_min[3],coord_max[3];
        coord_min[0]=CONFIG->read<double>("xmin");
        coord_min[1]=CONFIG->read<double>("ymin");
        coord_min[2]=CONFIG->read<double>("zmin");

        coord_max[0]=CONFIG->read<double>("xmax");
        coord_max[1]=CONFIG->read<double>("ymax");
        coord_max[2]=CONFIG->read<double>("zmax");

        dsl = new DelPhiShared(CONFIG->read<double>("Grid_scale"),
                          coord_min,coord_max,  
                          atoms,
                          CONFIG->read<int>("Max_Num_Atoms"),
                               CONFIG->read<double>("Domain_Shrinkage"),
                               CONFIG->read<bool>("Build_epsilon_maps"),
                               CONFIG->read<bool>("Build_status_map"),
                               CONFIG->read<bool>("Multi_Dielectric"));
    }

    ds = dsl;    

    /*
    vector<double> p;
    getGridPointCoordinates(0,0,0,p);

    std::cout << "\n" << p[0];
    std::cout << "\n" << p[1];
    std::cout << "\n" << p[2];

    exit(-1);
    */

    
    // Get surface from the factory
    Surface* ss = surfaceFactory().create(CONFIG,DS);	
	surf = ss;
    if (randDisplacement>=0)
        ss->setRandDisplacement(randDisplacement);
    return ss->build();
}

bool NanoShaper::triangulate()
{   
    if (surf!=nullptr && ds!=nullptr)
    {
        double val = SURF->triangulateSurface();			
        if (val==0)
            return false;            
        if (CONFIG->read<bool>("Smooth_Mesh"))	
            SURF->smoothSurface();        
        return true;
    }
    return false;
}


bool NanoShaper::getNearestGridPoint(double p[3],std::vector<unsigned>& indices)
{
    if (ds==nullptr)
        return false;
    
    if (DS->status==nullptr)
        return false;
        
    indices.clear();
    
    int ix = (int)rintp((p[0]-DS->xmin)/DS->side);
    int iy = (int)rintp((p[1]-DS->ymin)/DS->side);
    int iz = (int)rintp((p[2]-DS->zmin)/DS->side);
    
    bool outside = false;
    
    if (ix<0)
    {
        ix = 0;        
        outside = true;
    }
    if (ix>=DS->nx)
    {
        ix = DS->nx-1;
        outside = true;
    }
        
    if (iy<0)
    {
        iy = 0;        
        outside = true;
    }
    
    if (iy>=DS->ny)
    {
        iy = DS->ny-1;
        outside = true;
    }
        
    if (iz<0)
    {
        iz = 0;        
        outside = true;
    }
    if (iz>=DS->nz)
    {
        iz = DS->nz-1;
        outside = true;
    }
        
    indices.push_back((unsigned)ix);
    indices.push_back((unsigned)iy);
    indices.push_back((unsigned)iz);
                
    return outside;
}	     	

/** read config var*/        
template<class T> T NanoShaper::readConfig(const std::string& keyword)
{
    return CONFIG->read<T>(keyword);
}
   
/** set config var*/
template<class T> void NanoShaper::setConfig(std::string keyword,const T val)
{
    CONFIG->add<T>(keyword,val);
}

template bool NanoShaper::readConfig<bool>(const std::string& keyword);
template int NanoShaper::readConfig<int>(const std::string& keyword);
template double NanoShaper::readConfig<double>(const std::string& keyword);
template float NanoShaper::readConfig<float>(const std::string& keyword);
template unsigned NanoShaper::readConfig<unsigned>(const std::string& keyword);
template short NanoShaper::readConfig<short>(const std::string& keyword);
template long NanoShaper::readConfig<long>(const std::string& keyword);
template std::string NanoShaper::readConfig<std::string>(const std::string& keyword);
    
template void NanoShaper::setConfig<bool>(std::string keyword,const bool val);
template void NanoShaper::setConfig<int>(std::string keyword,const int val);
template void NanoShaper::setConfig<double>(std::string keyword,const double val);
template void NanoShaper::setConfig<float>(std::string keyword,const float val);
template void NanoShaper::setConfig<unsigned>(std::string keyword,const unsigned val);
template void NanoShaper::setConfig<short>(std::string keyword,const short val);
template void NanoShaper::setConfig<long>(std::string keyword,const long val);
template void NanoShaper::setConfig<std::string>(std::string keyword,const std::string val);


bool NanoShaper::colourGrid()
{
    // conventional execution mode
    if (!collectGridRays) 
	    return SURF->getSurf(CONFIG->read<bool>("Cavity_Detection_Filling"),
                         CONFIG->read<double>("Conditional_Volume_Filling_Value"));		
    // new execution mode, keeps track of intersections and normals at grid points
    else
    {   
        vector<packet> intersectionsInfo;
        bool cvn = CONFIG->read<bool>("Compute_Vertex_Normals");
        // impose to compute normals in order to collect them
        SURF->setComputeNormals(true);

    	bool retVal = SURF->getSurf(CONFIG->read<bool>("Cavity_Detection_Filling"),
                         CONFIG->read<double>("Conditional_Volume_Filling_Value"),&intersectionsInfo);
        
        if (!retVal)
            return false;

        // clean up in case it is needed
        for (unsigned k=0;k<3;k++)        
            raysArrayMap[k].clear(); 

        // it initializes the map with empty crossing_t objects such that the guest code can verify 
        // that the rays have been cast and no intersections have been found.
        // later we fill where intersections have been found

        std::vector<unsigned> grid;
        getGridSize(grid);

        //panel 0, dir x
        for (int iy = 0; iy < grid[1]; ++iy) 
        {
            for (int iz = 0; iz < grid[2]; ++iz) 
            {
                std::vector<double> p;
                unsigned dir = 0;                
                getGridPointCoordinates(0, iy, iz, p);
                
                std::array<double,2> point;    

                // y,z panel, x direction      
                point[0] = p[1];
                point[1] = p[2];
                
                raysArrayMap[dir][point] = crossings_t();
                crossings_t& temp = raysArrayMap[dir][point];                                   
                temp.point[0]=point[0];
                temp.point[1]=point[1];   
                temp.dir = dir;  
                temp.init = 1;
            }
        }

        //panel 2, dir y
        for (int ix = 0; ix < grid[0]; ++ix) 
        {
            for (int iz = 0; iz < grid[2]; ++iz) 
            {
                unsigned dir = 1;                
                std::vector<double> p;
                getGridPointCoordinates(ix, 0, iz, p);
                
                std::array<double,2> point;    

                // y,z panel, x direction      
                point[0] = p[0];
                point[1] = p[2];
                
                raysArrayMap[dir][point] = crossings_t();
                crossings_t& temp = raysArrayMap[dir][point];                                   
                temp.point[0]=point[0];
                temp.point[1]=point[1];   
                temp.dir = dir;  
                temp.init = 1;
            }
        }

        //panel 1, dir z
        for (int ix = 0; ix < grid[0]; ++ix) 
        {
            for (int iy = 0; iy < grid[1]; ++iy) 
            {
                unsigned dir = 2;                
                std::vector<double> p;
                getGridPointCoordinates(ix, iy, 0, p);
                
                std::array<double,2> point;    

                // y,z panel, x direction      
                point[0] = p[0];
                point[1] = p[1];
                
                raysArrayMap[dir][point] = crossings_t();
                crossings_t& temp = raysArrayMap[dir][point];                                   
                temp.point[0]=point[0];
                temp.point[1]=point[1];   
                temp.dir = dir;  
                temp.init = 1;
            }
        }


        unsigned num_packs = intersectionsInfo.size();
        //printf("\npacks %d",num_packs);       
        auto it = intersectionsInfo.begin();

        // reorganize the data according to the required data structure
		for (unsigned i_pack=0;i_pack<num_packs;i_pack++,it++)
		{      
            packet& pack = (*it);
            auto it_intersections = pack.first->begin();
            auto it_normals = pack.second->begin();                                
            unsigned psize = pack.first->size();
		    cout << endl << INFO << "Unpacking rays packet " << i_pack << " of size " << psize;
            
            // unpacking and reshaping
            for (unsigned k=0;k<psize;k++,it_intersections++,it_normals++)
            {
                coordVec& cv_int = (*it_intersections);
                coordVec& cv_norm = (*it_normals);

                std::array<double,2> point;    
                unsigned dir=cv_int.dir;

                // y,z panel, x direction
                if (dir==0)
                {    
                   //point = mappa(array(cv_int.iy,cv_int.iz));
                   point[0] = cv_int.vec[1];
                   point[1] = cv_int.vec[2];
                }

                // x,z panel, y direction
                else if (dir==1)
                {                                                
                    point[0] = cv_int.vec[0];
                    point[1] = cv_int.vec[2];
                }

                // x,y panel, z direction
                else if (dir==2)
                {                    
                    point[0] = cv_int.vec[0];
                    point[1] = cv_int.vec[1];
                }

                //printf("\n%f %f %f , dir %d",cv_int.vec[0],cv_int.vec[1],cv_int.vec[2],dir);

                crossings_t& crossings = raysArrayMap[dir][point];               
    
                // store the coordinate of the intersection only, not the entire vector
                crossings.inters.push_back(cv_int.vec[dir]);      
                // store the normal              
                crossings.normals.push_back(cv_norm.vec[0]);
                crossings.normals.push_back(cv_norm.vec[1]);
                crossings.normals.push_back(cv_norm.vec[2]);                                        

                // TODO cleanup the intersection vector and normal, as by now is no more useful
            } // for intersections
        } // for packet
    

      
      // restore the previous configuration on normals computation
      SURF->setComputeNormals(cvn); 
      return true;
    } // else
}

void NanoShaper::clean()
{
    if (CONFIG!=nullptr)    
        delete CONFIG;    
        
    if (DS!=nullptr)
        delete DS;
 
    if (SURF!=nullptr)
        delete SURF;
}


void NanoShaper::initConfig(const std::string* confFile)
{
	NanoShaper::clean();
    ConfigFile* cfl;
    
    // apply default configuration values. These values are particularly
    // well suite to colour a coarse grid
	if (confFile==nullptr)
	{       
        cfl = new ConfigFile();
        cfl->add<std::string>("Surface","skin");
        cfl->add<std::string>("Output_Folder", "" );
        cfl->add<bool>("Build_epsilon_maps", false );
        cfl->add<double>("Probe_Radius",1.4);    
        cfl->add<double>("Skin_Surface_Parameter",0.45);                
        cfl->add<bool>("Save_eps_maps", false );            
        cfl->add<bool>("Save_Status_map", false );
        cfl->add<bool>("Build_status_map", true );
        cfl->add<bool>("Save_ideb_map", false );
        cfl->add<bool>("Cavity_Detection_Filling", false );
        cfl->add<bool>("Project_boundary_grid_points", false );
        cfl->add<bool>("Accurate_Triangulation", false );
        cfl->add<bool>("Triangulation", false );
        cfl->add<std::string>("Operative_Mode","normal");
        cfl->add<bool>("Debug_Internals", false );
        #ifdef LIMIT_NUM_ATOMS
        cfl->add<int>("Max_Num_Atoms", -1);
        cfl->add<double>("Domain_Shrinkage", 0.);
        #endif
        cfl->add<double>("Conditional_Volume_Filling_Value", 11.4 );
        cfl->add<int>("Num_Wat_Pocket",2);
        cfl->add<double>("Grid_scale", 1.0 );
        cfl->add<double>("Grid_perfil", 80.0 );
        cfl->add<std::string>("XYZR_FileName","null");
        cfl->add<bool>("Multi_Dielectric",false);
        cfl->add<bool>("Smooth_Mesh", true );	
        cfl->add<bool>("Tri2Balls", false );
        cfl->add<bool>("Save_eps_maps", false );
        cfl->add<bool>("Save_bgps", false );	
        cfl->add<bool>("Save_Cavities", false );
        cfl->add<std::string>("Sys_Name","mol_from_api");
        cfl->add<int>("Number_thread", 1 );	
        cfl->add<bool>("Print_Available_Surfaces", false );	
        cfl->add<int>("Seed", 1 );
        cfl->add<bool>("Pockets_And_Cavities",true);
        cfl->add<bool>("Link_Pockets",false);	
        cfl->add<double>("Pocket_Radius_Big",3.0);		
        cfl->add<double>("Pocket_Radius_Small",1.4);		
        cfl->add<double>("Pocket_Radius_Link",1.0);	
        cfl->add<double>("Memb_Height",25.0);
        cfl->add<double>("Memb_Shift",1.0);
        cfl->add<std::string>("Root_FileName","");
        cfl->add<bool>("Compute_Vertex_Normals",false);
        cfl->add<bool>("Save_Mesh_MSMS_Format",false);
        cfl->add<bool>("Load_Balancing",true);
        cfl->add<double>("Blobbyness",-2.5);        
        cfl->add<std::string>("Surface_File_Name","triangulatedSurf.off");        
        cfl->add<bool>("Keep_Water_Shaped_Cavities",false);        
        cfl->add<int>("Max_Probes_Self_Intersections",100);        
        cfl->add<int>("Self_Intersections_Grid_Coefficient",1.5);        
        cfl->add<bool>("Check_duplicated_vertices",true);        
        cfl->add<bool>("Save_PovRay",false);        
        cfl->add<int>("Max_mesh_auxiliary_grid_size", 100 );
        cfl->add<int>("Max_mesh_patches_per_auxiliary_grid_cell", 250 );
        cfl->add<int>("Max_mesh_auxiliary_grid_2d_size", 200 );
        cfl->add<int>("Max_mesh_patches_per_auxiliary_grid_2d_cell",  500 );
        cfl->add<int>("Max_ses_patches_auxiliary_grid_size", 100 );
        cfl->add<int>("Max_ses_patches_per_auxiliary_grid_cell", 400 );
        cfl->add<int>("Max_ses_patches_auxiliary_grid_2d_size", 4000 );
        cfl->add<int>("Max_ses_patches_per_auxiliary_grid_2d_cell", 1 );
        cfl->add<int>("Max_skin_patches_auxiliary_grid_size", 100 );
        cfl->add<int>("Max_skin_patches_per_auxiliary_grid_cell", 400 );
        cfl->add<int>("Max_skin_patches_auxiliary_grid_2d_size", 6400 );
        cfl->add<int>("Max_skin_patches_per_auxiliary_grid_2d_cell", 1 );

        // new grid mode keywords
        cfl->add<bool>("PB_grid_mode",false);
        cfl->add<double>("xmin",-1.0);
        cfl->add<double>("ymin",-1.0);
        cfl->add<double>("zmin",-1.0);
        cfl->add<double>("xmax",-1.0);
        cfl->add<double>("ymax",-1.0);
        cfl->add<double>("zmax",-1.0);
	}
	else
	{
        // get data from configuration file
        try
        {
            cfl = new ConfigFile(confFile->c_str());
        }
        catch(...)
        {
            cout << endl << ERR << "Cannot read " << confFile << endl;
            exit(-1);			
        }    
    }   
    cf = (void*)cfl;
}
