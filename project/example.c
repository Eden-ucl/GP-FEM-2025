#include "../benchmark.h"
#include "fem.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

void elasticity_solve(const char *meshfile, const char *outfile, double E, double nu, double rho, double g) {
  // Read the mesh and the problem
  femGeo *theGeometry = geoMeshRead(meshfile);
  renumberNodes(theGeometry);
  femProblem *theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRAIN);

  // Boundary conditions are Dirichlet X and Y
  femElasticityAddBoundaryCondition(theProblem, "Base", DIRICHLET_X, 0.0);
  femElasticityAddBoundaryCondition(theProblem, "Base", DIRICHLET_Y, 0.0);    
  femElasticityAddBoundaryCondition(theProblem, "Symmetry", DIRICHLET_X, 0.0);
  femElasticityAddBoundaryCondition(theProblem, "Symmetry", DIRICHLET_Y, 0.0);    

  // Assemble and solve
  // femElasticityPrint(theProblem);
  //résolution sur le maillage renuméroté
  double *solNew = femElasticitySolve(theProblem);
  int nNodes     = theGeometry->theNodes->nNodes;

  //on re-crée un vecteur solOrig dans l’ordre initial
  double *solOrig = malloc(2 * nNodes * sizeof(double));
  for (int old = 0; old < nNodes; old++) {
    int newIdx = theGeometry->nodeInv[old];
    solOrig[2*old]   = solNew[2*newIdx];
    solOrig[2*old+1] = solNew[2*newIdx+1];
  }

  //on écrit solOrig, qui est dans l’ordre *original* des nœuds
  femSolutionWrite(nNodes, 2, solOrig, outfile);
  free(solOrig);


  // free the allocated ressources
  femElasticityFree(theProblem);
  femGeoFree(theGeometry);
}