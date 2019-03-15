/**
 * @file Poisson 3D in hexahedra with shock problem
 */
#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "ProblemConfig.h"

#include "mixedpoisson.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"

#include "pzintel.h"


#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"


#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzmatmixedpoisson3d.h"
#include <tuple>
#include <memory>

TPZCompMesh *CreateFluxHDivMesh(const ProblemConfig &problem);
TPZCompMesh *CreatePressureMesh(const ProblemConfig &problem);
TPZCompMesh *CreateMixMesh(const ProblemConfig &problem, TPZVec<TPZCompMesh *> &meshvec);
void CloneMeshVec(TPZVec<TPZCompMesh *> &meshvec, TPZVec<TPZCompMesh *> &meshvec_clone);
/// Increase the approximation orders of the sides of the flux elements
void IncreaseSideOrders(TPZCompMesh *fluxmesh);
void UniformRefinement(int nDiv, TPZGeoMesh *gmesh);
TPZGeoMesh *LMesh();
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff/*, TPZFMatrix<STATE> &flux*/);
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void ForcingCte(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
TPZGeoMesh *GMeshSphere(int dim, bool triang, int ndiv);
void BCCondition(const TPZVec<REAL> &pt, TPZVec<STATE> &solp);

TPZCompMesh *CMeshH1(const ProblemConfig &problem);
void SolveSystem(TPZAnalysis &an, TPZCompMesh *fCmesh, TPZVec<TPZCompMesh *> &meshvector);
void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);
void ErrorL2(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);

bool gmeshreader=true;
bool isH1=false;
int matId=1;

bool datatest=false;

//int bc0 = -1;
int bc1 = -1;
//int bc2 = -3;
//int bc3 = -4;
//int bc4 = -5;
//int bc5 = -6;

#define DUALMat

//#define UnitSquare

int dirichlet = 0;
int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

// Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    
    TPZGeoMesh *gmesh = nullptr;
    
    ProblemConfig config;
    config.porder = 1;
    config.problemname = "LaplaceOnSphere";
    config.dimension=2;
    int ndiv=3;
    
    if(gmeshreader){
#ifdef UnitSquare
        std::string meshfilename = "../Mysquare.msh";
#else
        std::string meshfilename = "../Mysphere.msh";
#endif
         TPZGmshReader gmsh;
        gmsh.GetDimNamePhysical()[1]["boundary"] =2;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        
        gmsh.SetFormatVersion("4.0");
        
#ifdef MACOSX
       gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
#else
        gmesh = gmsh.GeometricGmshMesh("Mysphere.msh");
#endif
        gmesh->SetDimension(config.dimension);
        config.gmesh = gmesh;
    }
    else{
    
    gmesh= GMeshSphere(config.dimension,  true, ndiv);
        //gmesh->SetDimension(dim);
        config.gmesh = gmesh;
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        UniformRefinement( ndiv,gmesh);
    
    }
    
    

    {
        std::ofstream out("Sgmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out,true);
        ofstream arg1("GeoMesh.txt");
        gmesh->Print(arg1);
    }
//
    if(datatest)
    {
        TPZVec<STATE> xp(3);
        REAL x=1;
        REAL y=2;
        REAL z=1;
        
        xp[0]=x;
        xp[1]=y;
        xp[2]=z;
        REAL r = sqrt(x*x+y*y+z*z);
        REAL phi =acos(z/r); // (atan2(sqrt(x*x+y*y),z));//
        REAL theta = atan2(y,x);
        
        TPZVec<STATE>sol(1);
        TPZFMatrix<STATE> flux(3,1);
        SolExata(xp, sol, flux);
        std::cout<<"("<< xp[0]<<", "<<xp[1]<<", "<<xp[2]<<")"<<std::endl;
        std::cout<<"f(x)= "<<sol[0]<<std::endl;
        std::cout<<"flux \n"<<std::endl;
        std::cout<<"flux(0,0) "<<flux(0,0)<<std::endl;
        std::cout<<"flux(1,0) "<<flux(1,0)<<std::endl;
        std::cout<<"flux(2,0) "<<flux(2,0)<<std::endl;
        std::cout<<"--(r,theta,phi)-- "<<std::endl;
        std::cout<<"r "<<r<<std::endl;
        std::cout<<"theta "<<theta<<std::endl;
        std::cout<<"phi"<<phi<<std::endl;
        
    
    std::cout << "x= \n"<<std::endl;
    
    
        return 0;
    }

    if (isH1) {
        TPZCompMesh *cmeshH1 = CMeshH1(config);
        
        
        {
            
            ofstream arg1("CompMesh.txt");
            cmeshH1->Print(arg1);
        }
        
        
        TPZAnalysis an(cmeshH1, true);
#ifdef USING_MKL
        TPZSymetricSpStructMatrix strmat(cmeshH1);
        strmat.SetNumThreads(0);
        //        strmat.SetDecomposeType(ELDLt);
        an.SetStructuralMatrix(strmat);
#else
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
        strmat.SetNumThreads(0);
        //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
        //        strmat3.SetNumThreads(8);
#endif
        
        
        
        
        {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("ExactSolution");
        an.DefineGraphMesh(config.dimension, scalnames, vecnames, "ExactSol.vtk");
        // Post processing
        an.PostProcess(0,config.dimension);
            
        }
        
        TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
        direct->SetDirect(ELDLt);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;
        an.Assemble();
        an.Solve();//resolve o problema 
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Solution");
        an.DefineGraphMesh(config.dimension, scalnames, vecnames, "H1Solution.vtk");
        // Post processing
        an.PostProcess(0,config.dimension);
            ofstream file("Solutout");
            an.Solution().Print("solution", file);
        
        return EXIT_SUCCESS;
    }
    
    //REsolvendo problema misto
    TPZCompMesh *fluxmesh=CreateFluxHDivMesh(config);
    
    TPZCompMesh *pressuremesh=CreatePressureMesh(config);
    {
        
        ofstream file1("MeshPressure.txt");
        ofstream file2("MeshFlux.txt");
        pressuremesh->Print(file1);
         fluxmesh->Print(file2);
    }
    
    TPZVec<TPZCompMesh *> meshvec(2);
    meshvec[0]=fluxmesh;
    meshvec[1]=pressuremesh;
    TPZCompMesh *mphysics=CreateMixMesh(config, meshvec);
    
    {
        
        ofstream file3("MultPhysics.txt");
        pressuremesh->Print(file3);
    }
    
    TPZAnalysis an(mphysics,false);
    SolveSystem(an, mphysics,meshvec);
    //Calculo de erro
    
    
  std::ofstream Errorfile("Error.txt");
    ErrorHDiv(fluxmesh, Errorfile, config.porder, ndiv);
    ErrorL2(pressuremesh, Errorfile, config.porder, ndiv);
    
    

        
        
    return 0;
}

TPZCompMesh *CreatePressureMesh(const ProblemConfig &problem) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    
    
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId, cmesh->Dimension());
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    
    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    
    int64_t n_connects = cmesh->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }

    
    
    return cmesh;
}

TPZCompMesh *CreateFluxHDivMesh(const ProblemConfig &problem) {
    int dim = problem.gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);

   // problem.gmesh->ResetReference();
    
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    cmesh->SetDimModel(dim);
    
    
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExata,problem.porder);
    material->SetForcingFunctionExact(solexata);
    
    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(ForcingH1, problem.porder);
    material->SetForcingFunction(force1);
    
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 1.);
        int bctype=0;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        TPZDummyFunction<STATE> *forceBound = new TPZDummyFunction<STATE>(BCCondition,problem.porder);
        bc->TPZMaterial::SetForcingFunction(forceBound);
        cmesh->InsertMaterialObject(bc);

    }
    
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();
    
    
   
    cmesh->InitializeBlock();
    return cmesh;

}

TPZCompMesh *CreateMixMesh(const ProblemConfig &problem, TPZVec<TPZCompMesh *> &meshvector) {
   
    
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    
    int dim=cmesh-> Dimension();
   
#ifdef DUALMat
        TPZMatMixedPoisson3D *mix = new TPZMatMixedPoisson3D(matId, dim);
#else
        TPZMixedPoisson *mix = new TPZMixedPoisson(matId, dim);
#endif
        TPZAutoPointer<TPZFunction<STATE> > solexata;
        solexata = new TPZDummyFunction<STATE>(SolExata,problem.porder);
        mix->SetForcingFunctionExact(solexata);
    
        TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(ForcingH1, problem.porder);
        mix->SetForcingFunction(force1);
    

    
    
    
     //   mix->SetInternalFlux(1);

    TPZMaterial *mat(mix);
    
    cmesh->InsertMaterialObject(mat);

    
    
    REAL coefk=1.;
    mix->SetPermeability(coefk);
    
    TPZFMatrix<REAL> TensorK(dim,dim,0.);
    TPZFMatrix<REAL> InvK(dim,dim,0.);
    for(int i = 0; i<dim; i++) {
        TensorK(i,i)=coefk;
        InvK(i,i)=coefk;
    }
    mix->SetPermeabilityTensor(TensorK, InvK);
    
    
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        if (matid == -2) {
            bctype = 0;
            val2.Zero();
        }
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        TPZDummyFunction<STATE> *forceBound = new TPZDummyFunction<STATE>(BCCondition,problem.porder);
        bc->TPZMaterial::SetForcingFunction(forceBound);
        cmesh->InsertMaterialObject(bc);
    }
    
    
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();

    meshvector[0] = CreateFluxHDivMesh(problem);
    meshvector[1] = CreatePressureMesh(problem);
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh);
    cmesh->LoadReferences();
    bool keepmatrix = false;
//    TPZCompMeshTools::CreatedCondensedElements(cmesh, true, keepmatrix);


    return cmesh;
}

void UniformRefinement(int nDiv, TPZGeoMesh *gmesh) {
    
    TPZManVector<TPZGeoEl*> children;
    for(int division = 0; division < nDiv; division++) {
        
        int64_t nels = gmesh->NElements();
        
        for(int64_t elem = 0; elem < nels; elem++) {
            
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            
            if(!gel || gel->HasSubElement()) continue;
            if(gel->Dimension() == 0) continue;
            gel->Divide(children);
        }
    }
}


TPZCompMesh *CMeshH1(const ProblemConfig &problem)
{
    
   
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    int dim =problem.gmesh->Dimension();
    /// criar materiais
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    //    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExata,problem.porder);
    material->SetForcingFunctionExact(solexata);
    
    
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingH1,problem.porder);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    
    ///Inserir condicao de contorno
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        if (matid == -2) {
            bctype = 0;
            val2.Zero();
        }
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        TPZDummyFunction<STATE> *forceBound = new TPZDummyFunction<STATE>(BCCondition,problem.porder);
        bc->TPZMaterial::SetForcingFunction(forceBound);
        cmesh->InsertMaterialObject(bc);
    }
    
    
    cmesh->InsertMaterialObject(mat);
    cmesh->SetAllCreateFunctionsContinuous();
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
    
}

void ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff/*, TPZFMatrix<STATE> &flux*/)
{
    REAL x,y,z;
    // Defining data
    REAL a=M_PI/2.0;
   
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
#ifdef UnitSquare
    ff[0]=-2*(-1 + x)*x - 2*(-1 + y)*y;
#else
    REAL r = sqrt(x*x+y*y+z*z);
    REAL phi =atan2(y,x);
    REAL theta = acos(z/r) ;
  //  flux.Resize(3, 1);
//
//    flux(0,0)=0.;
//    flux(1,0)=-6.*pow(sin(theta), 5)*cos(theta)*sin(phi)*sin(phi);
//    flux(2,0)= -(pow(sin(theta), 6)*sin(2*phi));
    
    ff[0]=(-1.)*(pow(sin(theta),4)*(-6*pow(sin(theta),2) + pow(cos(phi),2)*(2 + 6*pow(sin(theta),2)) +(-2 + 36*pow(cos(theta),2))*pow(sin(phi),2)))/pow(r,2);
#endif
    

}

void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
     flux.Resize(3, 1);
    
    
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica
    
    
    REAL x,y,z;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];

#ifdef UnitSquare
    
    //
    solp[0] = x*y*(1-y)*(1-x);
    flux(0,0)= y*(-1 - 2*x*(-1 + y) + y);
    flux(1,0)= x*(-1 + x + 2*y - 2*x*y);
    flux(2,0)= -2*(-1 + x)*x - 2*(-1 + y)*y; /// Because TPZMatMixedPoisson3D is mess
    
#else
    
//    REAL r = sqrt(x*x+y*y+z*z);
//    REAL phi =atan2(y,x);
//    REAL theta = acos(z/r) ;
    
    STATE r = sqrt(x*x+y*y+z*z);
    STATE theta = acos(z/r);
    STATE phi = atan2(y,x);
    
    STATE dfdr       = 0;
    STATE dfdTheta   = ((6*cos(theta)*(1 - pow(cos(phi),2))*pow(sin(theta),5))/r);
    STATE dfdPhi     = ((2*cos(phi)*pow(sin(theta),5)*sin(phi))/r);

    
    STATE costheta = cos(theta);
    STATE cosphi = cos(phi);
    STATE sintheta = sin(theta);
    STATE sinphi = sin(phi);
    
    // Gradient computations
    
    STATE Radialunitx = sintheta*cosphi;
    STATE Radialunity = sintheta*sinphi;
    STATE Radialunitz = costheta;
    
    STATE Thetaunitx = cosphi*costheta;
    STATE Thetaunity = costheta*sinphi;
    STATE Thetaunitz = -sintheta;
    
    STATE Phiunitx = -sinphi;
    STATE Phiunity = cosphi;
    STATE Phiunitz = 0.0;
    
    //
    solp[0] = pow(sin(theta),6)*(1-cos(phi)*cos(phi));

    flux(0,0) = -1.0*(dfdr * Radialunitx + dfdTheta * Thetaunitx + dfdPhi * Phiunitx);
    flux(1,0) = -1.0*(dfdr * Radialunity + dfdTheta * Thetaunity + dfdPhi * Phiunity);
//    flux(2,0) = -1.0*(dfdr * Radialunitz + dfdTheta * Thetaunitz + dfdPhi * Phiunitz);
    flux(2,0) = (-1.)*(pow(sin(theta),4)*(-6*pow(sin(theta),2) + pow(cos(phi),2)*(2 + 6*pow(sin(theta),2)) +(-2 + 36*pow(cos(theta),2))*pow(sin(phi),2)))/pow(r,2);
    
#endif
    
    
}


TPZGeoMesh *GMeshSphere(int dim, bool triang, int ndiv)
{
    
    
    bool ftriangulo = triang;
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int64_t materialId = matId;
    //int64_t arc1 = bc1; // -1;
    
    int nodes = 6;//quantidade de nos da malha geometrica
    int elbc=4;
    int nel=8+elbc;
    geomesh->NodeVec().Resize(nodes);
    geomesh->SetDimension(dim);
    int r=1.;//rate of sphere
    TPZManVector<REAL,3> xc(3,0.);
    //center of sphere
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    //nos
    
    REAL coord[6][3] =
    {
        { xc[0]+r,0,0},
        {xc[0]-r,0,0},
        {0,xc[1]+r,0.},
        {0,xc[1]-r,0},
        {0,0,xc[2]+r},
        {0,0,xc[2]-r},
    };
    
    int elno[12][3]=
    {
        {0,2,4},
        {2,1,4},
        {1,3,4},
        {3,0,4},
        {0,2,5},
        {2,1,5},
        {1,3,5},
        {3,0,5},
        //el de contorno
        {0,2},
        {2,1},
        {1,3},
        {3,0}
    };
    
    
    
    
    
    geomesh->NodeVec().Resize(nodes);
    geomesh->SetDimension(dim);
    
    for(int no=0; no<nodes; no++)
    {
        TPZManVector<REAL,3> xco(3,0.);
        xco[0] = coord[no][0];
        xco[1] = coord[no][1];
        xco[2] = coord[no][2];
        geomesh->NodeVec()[no].Initialize(xco, *geomesh);
    }
    
    for(int el=0; el<nel-elbc; el++)
    {
        
        TPZVec<int64_t>elnodes(3);
        for (int i=0; i<3; i++) {
            elnodes[i] = elno[el][i];
        }
        //TPZGeoEl *gel = new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elnodes, matId,*geomesh);
        TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
        new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elnodes,matId,*geomesh);
        SphereRingT1->Geom().SetData(r,xc);
    }
   // geomesh->BuildConnectivity();
    
    int elementid=8;
    TPZVec<int64_t> topology(2);
    // Definition of Arc coordenates
        // Create Geometrical Arc #0
        topology[0] = 0;
        topology[1] = 2;
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, bc1 /** arc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
        }
        elementid++;
    
        // Create Geometrical Arc #1
        topology[0] = 2;
        topology[1] = 1;
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, bc1 /** arc4 */, *geomesh);
            // arc->Geom().Print(std::cout);
        }
        elementid++;
    
        // Create Geometrical Arc #2
        topology[0] = 1;
        topology[1] = 3;
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, bc1 /** arc1 */, *geomesh);
            //arc->Geom().Print(std::cout);
        }
        elementid++;
    
        // Create Geometrical Arc #3
        topology[0] = 3;
        topology[1] = 0;
        {
            TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, bc1 /** arc2 */, *geomesh);
            //arc->Geom().Print(std::cout);
        }
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    
    return geomesh;
    
    
}


void BCCondition(const TPZVec<REAL> &pt, TPZVec<STATE> &bc){
    
    TPZFMatrix<STATE> flux;
    SolExata(pt, bc, flux);
}

void SolveSystem(TPZAnalysis &an, TPZCompMesh *fCmesh, TPZVec<TPZCompMesh *> &meshvector)
{

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(fCmesh);
    strmat.SetNumThreads(2);
    //        strmat.SetDecomposeType(ELDLt);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
    strmat.SetNumThreads(2);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif
    
    int dim=fCmesh->Reference()->Dimension();
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();
    
    fCmesh->LoadSolution(an.Solution());
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, fCmesh);
    
//    meshvector[0]->Solution().Print("p = ",std::cout,EMathematicaInput);
    
    //Pos processamento
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    scalnames.Push("ExactPressure");
    vecnames.Push("ExactFlux");
    scalnames.Push("Divergence");
    scalnames.Push("ExactDiv");
    an.DefineGraphMesh(dim, scalnames, vecnames, "MixedSolution.vtk");
    // Post processing
    an.PostProcess(0,dim);

}


void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue;
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, NULL);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with HDiv space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
//    out << "L2 Norm for flux - "<< endl;
//    out<< "L2 Norm for divergence - Hdiv Norm for flux " << endl;
//    out <<  setw(16) << sqrt(globalerrors[1]) <<endl;
//    out<< setw(25)  << sqrt(globalerrors[2]) << setw(21)  << sqrt(globalerrors[3]) << endl;
    
        out << "L2 Norm for flux = "    << sqrt(globalerrors[1]) << endl;
        out << "L2 Norm for divergence = "    << sqrt(globalerrors[2])  <<endl;
        out << "Hdiv Norm for flux = "    << sqrt(globalerrors[3])  <<endl;
    
}

void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv)
{
    int64_t nel = l2mesh->NElements();
    //int dim = l2mesh->Dimension();
    TPZManVector<REAL,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, NULL);
        int nerr = elerror.size();
        globalerrors.resize(nerr);
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with L2 space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
    out << "L2 Norm for pressure = "    << sqrt(globalerrors[1]) << endl;
}
