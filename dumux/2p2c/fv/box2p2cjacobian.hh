#ifndef DUNE_BOX2P2CJACOBIAN_HH
#define DUNE_BOX2P2CJACOBIAN_HH

#include<map>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<sstream>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dune/disc/functions/p1function.hh>
//#include"../boxjacobian_2p2c.hh"
#include"dumux/operators/boxjacobian.hh"
#include"dumux/2p2c/2p2cproblem.hh"
//#include "varswitch.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for box scheme for two-phase two-component flow equation
 * @author Peter Bastian
 */



namespace Dune
{
  /** @addtogroup DISC_Disc
   *
   * @{
   */
  /**
   * @brief compute local jacobian matrix for the boxfile for two-phase two-component flow equation
   *
   */


  //! Derived class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for 
	diffusion equation

	    div j = q; j = -K grad u; in Omega

		u = g on Gamma1; j*n = J on Gamma2.

	Uses conforming finite elements with the Lagrange shape functions.
	It should work for all dimensions and element types.
	All the numbering is with respect to the reference element and the
	Lagrange shape functions

	Template parameters are:

	- Grid  a DUNE grid type
	- RT    type used for return values 
  */
  template<class G, class RT, class BoxFunction = LeafP1Function<G, RT, 2> >
  class Box2P2CJacobian 
    : public BoxJacobian<Box2P2CJacobian<G,RT,BoxFunction>,G,RT,2,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef Box2P2CJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::MBlockType MBlockType;
    typedef FVElementGeometry<G> FVElementGeometry;

 	enum {pWIdx = 0, satNIdx = 1, numberOfComponents = 2};	// Solution vector index
	enum {wPhase = 0, nPhase = 1};									// Phase index
	enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};		// Phase state
	enum {water = 0, air = 1};										// Component index					
	
  public:
    // define the number of phases (m) and components (c) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=G::dimension};
    enum {m=2, c=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};

    typedef FieldMatrix<RT,dim,dim> FMatrix;
    typedef FieldVector<RT,dim> FVector;
    
    //! Constructor
    Box2P2CJacobian (TwoPTwoCProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_), 
      problem(params), 
      sNDat(this->vertexMapper.size()),
      sIPDat(this->vertexMapper.size())
    {
      this->analytic = false;
    }


    // Compute time dependent term (storage), loop over nodes / subcontrol volumes
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node)
    {
    	 GeometryType gt = e.geometry().type();
    	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
     	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
    	 
   	 int globalIdx = this->vertexMapper.template map<dim>(e, sfs[node].entity());
   	 //int globalCoord = this->fvGeom.subContVol[node].global;

   	 VBlockType result; 
   	 RT satN = vNDat[node].saturationN;
   	 RT satW = vNDat[node].saturationW;
   	    	    	                  
   	 // storage of component water
   	 result[water] = 
   		 sNDat[globalIdx].porosity*(vNDat[node].density[wPhase]*satW*vNDat[node].massfrac[water][wPhase]
   		                 +vNDat[node].density[nPhase]*satN*vNDat[node].massfrac[water][nPhase]);
   	 // storage of component air
   	 result[air] = 
   		 sNDat[globalIdx].porosity*(vNDat[node].density[nPhase]*satN*vNDat[node].massfrac[air][nPhase]
   	                    +vNDat[node].density[wPhase]*satW*vNDat[node].massfrac[air][wPhase]);   
   	 
   	 //std::cout << result << " " << node << std::endl;
   	 return result;
    };

    // loop over subcontrol volume faces
    virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
   	 int i = this->fvGeom.subContVolFace[face].i;
     	 int j = this->fvGeom.subContVolFace[face].j;
     	 
     	 // normal vector, value of the area of the scvf
    	 const FieldVector<RT,dim> normal(this->fvGeom.subContVolFace[face].normal);

    	 // get global coordinates of nodes i,j
    	 const FieldVector<DT,dim> global_i = this->fvGeom.subContVol[i].global;
		 const FieldVector<DT,dim> global_j = this->fvGeom.subContVol[j].global;

		 FieldMatrix<RT,m,dim> pGrad(0.), xGrad(0.); 
		 FieldVector<RT,dim> temp(0.); 
//		 FieldVector<RT,dim> Kij(0);
     	 VBlockType flux(0.);
		 
		 // effective permeability in edge direction 
     	 // RT Kij = sIPDat[global_j].K_eff[face]; 
     	 
     	 // calculate harmonic mean of permeabilities of nodes i and j
		 const FMatrix K = harmonicMeanK(global_i, global_j);
		 //const FMatrix K = harmonicMeanK(e, face);		  
     	 //K.umv(normal, Kij);  // Kij=K*n
		 
		 // calculate FE gradient (grad p for each phase)
		 for (int k = 0; k < this->fvGeom.nNodes; k++) // loop over adjacent nodes
		 {	 
			 // FEGradient at node k
			 const FieldVector<DT,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]);
			 FieldVector<RT,m> pressure(0.0), massfrac(0.0);

			 pressure[wPhase] = vNDat[k].pW;
			 pressure[nPhase] = vNDat[k].pN;
      	 
      	 // compute sum of pressure gradients for each phase
      	 for (int phase = 0; phase < m; phase++)
      	 {	      		 
      		 temp = feGrad;
      		 temp *= pressure[phase];
      		 pGrad[phase] += temp;

          	 // compute sum of concentration gradient
         	 temp = feGrad;
         	 temp *= vNDat[k].massfrac[air][phase];
         	 xGrad[phase] += temp;
      	 }

       }

       // deduce gravity*density of each phase
		 FieldMatrix<RT,m,dim> contribComp(0);
		 for (int phase=0; phase<m; phase++)
		 {
			 contribComp[phase] = problem.gravity();
			 contribComp[phase] *= vNDat[i].density[phase];  
			 pGrad[phase] -= contribComp[phase]; // grad p - rho*g
		 }
	 	 
		 VBlockType outward(0);  // Darcy velocity of each phase
		 FieldVector<RT,dim> v_tilde(0);

		 // calculate the advective flux using upwind: K*n(grad p -rho*g)
		 for (int phase=0; phase<m; phase++) 
		 	{
			 //pGrad[phase] *= Kij;
			 //outward[phase] = Kij * pGrad[phase];// * normal;
	     	 K.umv(pGrad[phase], v_tilde);  // v_tilde=K*gradP
	     	 outward[phase] = v_tilde*normal;
		 	}
		 
		 // evaluate upwind nodes
		 int up_w, dn_w, up_n, dn_n;
		 if (outward[wPhase] <= 0) {up_w = i; dn_w = j;}
		 else {up_w = j; dn_w = i;};
		 if (outward[nPhase] <= 0) {up_n = i; dn_n = j;}
		 else {up_n = j; dn_n = i;};


		 RT alpha = 1.0;  // Upwind parameter
 		 
		 // Water conservation
 		 flux[water] =   (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
 				               * vNDat[up_w].massfrac[water][wPhase] 
 				    + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase]
 				               * vNDat[dn_w].massfrac[water][wPhase])
 				               * outward[wPhase]; 		 		
 		 flux[water] +=  (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
 				               * vNDat[up_n].massfrac[water][nPhase] 
 				    + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase]
 				               * vNDat[dn_n].massfrac[water][nPhase])
 				               * outward[nPhase];
 		 // Air conservation
 		 flux[air]   =   (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
 				               * vNDat[up_n].massfrac[air][nPhase] 
 				    + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase]
 				               * vNDat[dn_n].massfrac[air][nPhase])
 				               * outward[nPhase];
 		 flux[air]  +=   (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
 				               * vNDat[up_w].massfrac[air][wPhase]
 				    + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase]
 				               * vNDat[dn_w].massfrac[air][wPhase])
 				               * outward[wPhase];
		 
		 // DIFFUSION
		 VBlockType normDiffGrad;

		 //FieldMatrix<DT,m,c> diffusionW;
		 RT diffusionAW, diffusionWW; // Diffusion in the water phase
		 RT diffusionWN, diffusionAN; // Diffusion in the non-wetting phase
		 VBlockType avgDensity, avgDpm;
		 avgDpm[wPhase]=1e-9; // needs to be changed !!!
		 avgDpm[nPhase]=1e-9;
		 
		 normDiffGrad[wPhase] = -(xGrad[wPhase]*normal);
		 normDiffGrad[nPhase] = -(xGrad[nPhase]*normal);

		 // calculate the arithmetic mean of densities
		 avgDensity[wPhase] = 0.5*(vNDat[i].density[wPhase] + vNDat[j].density[wPhase]);
		 avgDensity[nPhase] = 0.5*(vNDat[i].density[nPhase] + vNDat[j].density[nPhase]);

		 
		 diffusionAW = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad[wPhase];
		 diffusionWW = - diffusionAW;
		 diffusionWN = avgDpm[nPhase] * avgDensity[nPhase] * normDiffGrad[nPhase];
		 diffusionAN = - diffusionWN;
		 
//		 std::cout << "Diffusive Flux: " << diffusionWG << std::endl; 
		 

		 // add water diffusion to flux
		 flux[water] -= (diffusionWW + diffusionWN);
//		 std::cout << "Water Flux: " << flux[water] << std::endl; 
		 // air diffusion not implemeted
		 flux[air] -= (diffusionAW + diffusionAN);
//		 std::cout << "Air Flux: " << flux[air] << std::endl; 


		 return flux;
  };
    
   	// Integrate sources / sinks
   	virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
   	{
   		// ASSUME problem.q already contains \rho.q
   		return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
   	}


   	virtual void primaryVarSwitch (const int global, VBlockType* sol, const int local)
    {
        //int k = this->fvGeom.cellLocal;

        // neue Variable in SOL schreiben. Da Zeiger auf sol, sollte das funktionieren
        // Achtung: alte Loesung beibebehalten, falls Newton den Zeitschritt nicht schafft
   		// 1=gas, 2=water, 3=both
   		
        int state = sNDat[global].phaseState;
        bool switched = false;
    	RT pWSat = problem.multicomp().vaporPressure(vNDat[local].temperature);
    	RT henry = problem.multicomp().henry(vNDat[local].temperature);
        
        switch(state) 
        {
        case gasPhase :
        	RT pwg;
        	pwg = vNDat[local].massfrac[water][nPhase] * vNDat[local].pN;
  
        	if (pwg > 1.01*pWSat && switched == false)
            {
            	// appearance of water phase
            	std::cout << "Water appears at node " << global << std::endl;
            	sNDat[global].phaseState = bothPhases;
            	vNDat[local].saturationN = 1 - 1.E-6; // initialize
            	sol[global][satNIdx] = vNDat[local].saturationN;
            	switched = true;
            }
            break;

        case waterPhase :
        	RT pbub;
        	pbub = pWSat + vNDat[local].massfrac[air][wPhase]*henry;
        	if (pbub > vNDat[local].pN && switched == false)
            {
            	// appearance of gas phase
            	std::cout << "Gas appears at node " << global << std::endl;
            	sNDat[global].phaseState = bothPhases;
            	vNDat[local].saturationN = 1.E-6;  // initialize
            	sol[global][satNIdx] = vNDat[local].saturationN;
            	switched = true;
            }
            break;

        case bothPhases :
      	  	if (sol[global][satNIdx] < 0.0  && switched == false)
      	  	{
      		  	// disappearance of gas phase
      		  	std::cout << "Gas disappears " << global << std::endl;
      		  	sNDat[global].phaseState = waterPhase;
      		  	//vNDat[local].saturationN = ;  // Initialisierung
      		  	sol[global][satNIdx] = vNDat[local].massfrac[air][wPhase];
      		  	switched = true;
            }
      	  	else if ((1-sol[global][satNIdx]) < 0.0  && switched == false) //vNDat[i].saturationW 
      	  	{
      	  		// disappearance of water phase
      	  		std::cout << "Water disappears " << global << std::endl;
      	  		sNDat[global].phaseState = gasPhase;
      	  		//vNDat[local].massfrac[water][nPhase] = 1.E-6;  // Initialisierung
      	  		sol[global][satNIdx] = vNDat[local].massfrac[water][nPhase];
            	switched = true;
      	  	}
      	  	break;
        }
        
   	return;
    }
    
    // harmonic mean computed directly
    virtual FMatrix harmonicMeanK (const FVector global_i, const FVector global_j) const
    {
    	double eps = 1e-20;

    	FMatrix Ki, Kj;
   	 
    	Ki = this->problem.K(global_i); 
    	Kj = this->problem.K(global_j);

    	for (int kx=0; kx<dim; kx++){
    		for (int ky=0; ky<dim; ky++){
    			if (Ki[kx][ky] != Kj[kx][ky])
    			{
    				Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
    			}
    		}
    	}
   	 return Ki;
    }
    
//    // or from the staticNode Vector
//    virtual FMatrix harmonicMeanK (const Entity& e, int k) const
//    {
//	 FMatrix Ki, Kj;
//	 const RT eps = 1e-20;
//	 
//    GeometryType gt = e.geometry().type();
//    const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
//    sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
//
//     int i = this->fvGeom.subContVolFace[k].i;
//     int j = this->fvGeom.subContVolFace[k].j;
//
//     int global_i = this->vertexMapper.template map<dim>(e, sfs[i].entity());
//     int global_j = this->vertexMapper.template map<dim>(e, sfs[j].entity());
//    
//    
//     	Ki = sNDat[global_i].K;
//     	Kj = sNDat[global_j].K;
//       	 
//     	for (int kx=0; kx<dim; kx++){
//     		for (int ky=0; ky<dim; ky++){
//     			if (Ki[kx][ky] != Kj[kx][ky])
//     			{
//     				Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
//     			}
//     		}
//     	}
//     	return Ki;
//    }
    
    virtual void clearVisited ()
    {
    	for (int i = 0; i < this->vertexMapper.size(); i++)
   		sNDat[i].visited = false;
   	 
   	 return;
   	}
         
    
	  //*********************************************************
	  //*														*
	  //*	Calculation of Data at Elements (elData) 			*
	  //*						 								*
	  //*														*
	  //*********************************************************

    virtual void computeElementData (const Entity& e)
    {
//  	 // ASSUME element-wise constant parameters for the material law 
// 		 elData.parameters = problem.materialLawParameters
// 		 (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
//   	 
//		 // ASSUMING element-wise constant permeability, evaluate K at the cell center 
// 		 elData.K = problem.K(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);  
//
//		 // ASSUMING element-wise constant porosity 
// 		 elData.porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
   	 return;
    }

    
	  //*********************************************************
	  //*														*
	  //*	Calculation of Data at Nodes that has to be			*
	  //*	determined only once	(sNDat)					*
	  //*														*
	  //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData (const Entity& e, VBlockType* sol)
    {
   	 // size of the sNDat vector is determined in the constructor
   	 
   	 // local to global id mapping (do not ask vertex mapper repeatedly
   	 //int localToGlobal[LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize];

   	 // get access to shape functions for P1 elements
   	 GeometryType gt = e.geometry().type();
   	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
   	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

   	 // get local to global id map
   	 for (int k = 0; k < sfs.size(); k++) {
  		 const int globalIdx = this->vertexMapper.template map<dim>(e, sfs[k].entity());
  		  
  		 // if nodes are not already visited
  		 if (!sNDat[globalIdx].visited) 
  		  {
  			  // ASSUME parameters defined at the node/subcontrol volume for the material law 
//  			  sNDat[globalIdx].parameters = problem.materialLawParameters
//  			  (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

 			  // global coordinates
 			  FieldVector<DT,dim> global_i = this->fvGeom.subContVol[k].global;
  			  
  			  // ASSUME porosity defined at nodes
  			  sNDat[globalIdx].porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
  			  
  			  // evaluate primary variable switch
 			  primaryVarSwitch(globalIdx, sol, k);
  			  
  			  // mark elements that were already visited
  			  sNDat[globalIdx].visited = true;
  		  }
  	  }
  	  
	  return;
    }
    

	  //*********************************************************
	  //*														*
	  //*	Calculation of variable Data at Nodes				*
	  //*	(vNDat)											*
	  //*														*
	  //*********************************************************
   

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
    virtual void updateVariableData (const Entity& e, const VBlockType* sol)
    {
   	 vNDat.resize(this->fvGeom.nNodes);
   	 int size = vNDat.size();

   	 GeometryType gt = e.geometry().type();
   	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
   	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
     
     FieldVector<RT, 4> parameters = problem.materialLawParameters 
     (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

   	 for (int i = 0; i < size; i++) 
   	 {
   	   	 const int global = this->vertexMapper.template map<dim>(e, sfs[i].entity());

   	   	 if (sNDat[global].phaseState == bothPhases){
   			 vNDat[i].saturationN = sol[i][satNIdx];
   			 vNDat[i].saturationW = 1.0 - sol[i][satNIdx];
   	   		 // Solubilities of components in phases
   	   		 vNDat[i].massfrac[air][wPhase] = problem.multicomp().xAW(vNDat[i].pN, vNDat[i].temperature);
   	   		 vNDat[i].massfrac[water][wPhase] = 1.0 - vNDat[i].massfrac[air][wPhase];
   	   		 vNDat[i].massfrac[water][nPhase] = problem.multicomp().xWN(vNDat[i].pN, vNDat[i].temperature);
   	   		 vNDat[i].massfrac[air][nPhase] = 1.0 - vNDat[i].massfrac[water][nPhase];
   		 }
   		 else{
   	   		 vNDat[i].massfrac[air][wPhase] = sol[i][satNIdx];
   	   		 vNDat[i].massfrac[water][wPhase] = 1.0 - vNDat[i].massfrac[air][wPhase];
   	   		 vNDat[i].massfrac[water][nPhase] = problem.multicomp().xWN(vNDat[i].pN, vNDat[i].temperature);
   	   		 vNDat[i].massfrac[air][nPhase] = 1.0 - vNDat[i].massfrac[water][nPhase];
   			 
   		 }
   		 
   		 vNDat[i].pW = sol[i][pWIdx];
   		 vNDat[i].pC = problem.materialLaw().pC(vNDat[i].saturationW, parameters);
   		 vNDat[i].pN = sol[i][pWIdx] + vNDat[i].pC;
   		 vNDat[i].temperature = 283.15; // in [K]
   		 // Mobilities & densities
   		 vNDat[i].mobility[wPhase] = problem.materialLaw().mobW(vNDat[i].saturationW, parameters);
   		 vNDat[i].mobility[nPhase] = problem.materialLaw().mobN(sol[i][satNIdx], parameters);
   		 vNDat[i].density[wPhase] = problem.materialLaw().wettingPhase.density();
   		 vNDat[i].density[nPhase] = problem.materialLaw().nonwettingPhase.density();
   		 // Solubilities of components in phases
   		 vNDat[i].massfrac[air][wPhase] = problem.multicomp().xAW(vNDat[i].pN, vNDat[i].temperature);
   		 vNDat[i].massfrac[water][wPhase] = 1.0 - vNDat[i].massfrac[air][wPhase];
   		 vNDat[i].massfrac[water][nPhase] = problem.multicomp().xWN(vNDat[i].pN, vNDat[i].temperature);
   		 vNDat[i].massfrac[air][nPhase] = 1.0 - vNDat[i].massfrac[water][nPhase];
   		
         // CONSTANT solubility (for comparison with twophase)
//         vNDat[i].massfrac[air][wPhase] = 0.1; vNDat[i].massfrac[water][wPhase] = 0.9;
//         vNDat[i].massfrac[water][nPhase] = 0.1; vNDat[i].massfrac[air][nPhase] = 0.9;

         //std::cout << "water in gasphase: " << vNDat[i].massfrac[water][nPhase] << std::endl;
         //std::cout << "air in waterphase: " << vNDat[i].massfrac[air][wPhase] << std::endl;
   	 }   	 
   	 return;
    }
    

    
    struct StaticNodeData 
    {
   	 bool visited;
   	 int phaseState;
   	 RT cellVolume;
   	 RT porosity;
   	 FieldVector<RT, 4> parameters;
   	 FMatrix K;
    };
    
    struct VariableNodeData  
    {
   	 RT saturationN;
     RT saturationW;
     RT pW;
     RT pC;
     RT pN;
     RT temperature;
     VBlockType mobility;  //Vector with the number of phases
     VBlockType density;
     FieldMatrix<RT,c,m> massfrac;
    };

	 struct StaticIPData
	 {
		 bool visited;
		 FMatrix K;
	 };
    
    
    struct ElementData {
//   	 RT cellVolume;
//     	 RT porosity;
//   	 RT gravity;
//   	 FieldVector<RT, 4> parameters;
//   	 FieldMatrix<RT,dim,dim> K;
   	 } elData;
   	    	 
    
    // parameters given in constructor
    TwoPTwoCProblem<G,RT>& problem;
    CWaterAir multicomp;
    std::vector<StaticNodeData> sNDat;
    std::vector<StaticIPData> sIPDat;
    std::vector<VariableNodeData> vNDat;
  };  
  
}
#endif
