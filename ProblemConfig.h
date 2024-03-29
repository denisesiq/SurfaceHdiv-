//
//  ProblemConfig.h
//  ErrorEstimation
//
//  Created by Philippe Devloo on 09/05/18.
//

#ifndef ProblemConfig_h
#define ProblemConfig_h

#include <set>
#include "TPZAnalyticSolution.h"

struct ProblemConfig
{
    TPZGeoMesh *gmesh = 0;
    int porder = 1;
    int hdivmais = 1;
    bool makepressurecontinuous = 0;
    std::string problemname;
    std::set<int> materialids;
    std::set<int> bcmaterialids;
    TLaplaceExample1 exact;
    TLaplaceExample1 forcingCte;
    int dimension=2;
    
    ProblemConfig(){};
    ProblemConfig(const ProblemConfig &cp) : gmesh(cp.gmesh), porder(cp.porder), hdivmais(cp.hdivmais),
    makepressurecontinuous(cp.makepressurecontinuous),
    problemname(cp.problemname),
    materialids(cp.materialids), bcmaterialids(cp.bcmaterialids),exact(cp.exact),forcingCte(cp.forcingCte),dimension(cp.dimension)
    {
    }
    
    ProblemConfig &operator=(const ProblemConfig &cp)
    {
        gmesh = cp.gmesh;
        porder = cp.porder;
        hdivmais = cp.hdivmais;
        makepressurecontinuous = cp.makepressurecontinuous;
        problemname = cp.problemname;
        materialids = cp.materialids;
        bcmaterialids = cp.bcmaterialids;
        exact = cp.exact;
        forcingCte=cp.forcingCte;
        dimension=cp.dimension;
        return *this;
    }
};

#endif /* ProblemConfig_h */
