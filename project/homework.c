#include "fem.h"
double *femSystemChebyschev(femFullSystem *theSystem, double l_max, double l_min);

double *femElasticitySolve(femProblem *theProblem) {

    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
    femGeo *theGeometry = theProblem->geometry;
    femNodes *theNodes = theGeometry->theNodes;
    femMesh *theMesh = theGeometry->theElements;

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int map[4], mapX[4], mapY[4];
    int iElem,iInteg,i,j;
    int nLocal = theMesh->nLocalNode;

    double a = theProblem->A;
    double b = theProblem->B;
    double c = theProblem->C;
    double rho = theProblem->rho;
    double g = theProblem->g;
    double **A = theSystem->A;
    double *B = theSystem->B;

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
        } 
            
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; 
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                        dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                        dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                        dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                        dphidx[i] * c * dphidx[j]) * jac * weight; 
                    }
            }
            for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; 
            }
        }
    }
    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (i = 0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
        double value = theProblem->conditions[theConstrainedNodes[i]]->value;
        femFullSystemConstrain(theSystem, i, value);
        }
    }
    double lambda_max, lambda_min;
    double *v = malloc((theSystem -> size)* sizeof(double));

    // Estimate largest eigenvalue
    power_iteration_system(theSystem, &lambda_max, v);

    // Estimate smallest eigenvalue
    inverse_power_system(theSystem, &lambda_min, v); 
    double *soluce = malloc(theSystem->size * sizeof(double));
    soluce = femSystemChebyschev(theSystem, lambda_max, lambda_min);
    free(v);
    return soluce;
}


double *femFullSystemEliminate(femFullSystem *mySystem) {
    double **A = mySystem->A;
    double *B = mySystem->B;
    int size = mySystem->size;
    double tol = 1e-16;
    int i, j, k;
    
    // Calculer la largeur de bande effective.
    // Pour chaque ligne i, on parcourt les colonnes j >= i et on retient le maximum de (j - i)
    // pour lesquelles A[i][j] n'est pas négligeable.
    int bw = 0;
    for (i = 0; i < size; i++) {
        for (j = i; j < size; j++) {
            if (fabs(A[i][j]) > 1e-12 && (j - i) > bw) {
                bw = j - i;
            }
        }
    }
    printf("Solveur bande : largeur de bande = %d\n", bw);

    // Elimination de Gauss restreinte à la bande.
    for (k = 0; k < size; k++) {
        if (fabs(A[k][k]) <= tol) {
            printf("Pivot index %d, valeur %e\n", k, A[k][k]);
            Error("Cannot eliminate with such a pivot");
        }
        // Pour la ligne k, seules les lignes i de k+1 à min(size, k + bw + 1) peuvent être non nulles.
        int i_max = (k + bw + 1 < size) ? (k + bw + 1) : size;
        for (i = k + 1; i < i_max; i++) {
            double factor = A[i][k] / A[k][k];
            // De même, dans la ligne i, seules les colonnes de k+1 à min(size, k + bw + 1) peuvent être non nulles.
            int j_max = (k + bw + 1 < size) ? (k + bw + 1) : size;
            for (j = k + 1; j < j_max; j++) {
                A[i][j] -= factor * A[k][j];
            }
            B[i] -= factor * B[k];
        }
    }

    // Rétro-substitution dans la bande.
    for (i = size - 1; i >= 0; i--) {
        double sum = 0.0;
        int j_max = (i + bw + 1 < size) ? (i + bw + 1) : size;
        for (j = i + 1; j < j_max; j++) {
            sum += A[i][j] * B[j];
        }
        B[i] = (B[i] - sum) / A[i][i];
    }
    
    return B;
}

void renumberNodes(femGeo *theGeometry) {
    int nNodes = theGeometry->theNodes->nNodes;
    /* --- 1. Construction de la liste d'adjacence --- */
    int *rcm_order = malloc(nNodes * sizeof(int));
    typedef struct {
        int count;
        int capacity;
        int *neighbors;
    } NodeNeighbors;
    NodeNeighbors *adj = malloc(nNodes * sizeof(NodeNeighbors));
    for (int i = 0; i < nNodes; i++) {
        adj[i].count = 0;
        adj[i].capacity = 4;
        adj[i].neighbors = malloc(adj[i].capacity * sizeof(int));
    }

    femMesh *mesh = (theGeometry->theElements) ? theGeometry->theElements : theGeometry->theEdges;
    if (!mesh) {
        fprintf(stderr, "Aucune connectivité disponible pour la renumérotation RCM.\n");
        free(rcm_order);
        free(adj);
        return;
    }
    int nElem = mesh->nElem;
    int nLocal = mesh->nLocalNode;
    int *elem = mesh->elem;
    for (int e = 0; e < nElem; e++) {
        for (int i = 0; i < nLocal; i++) {
            for (int j = i + 1; j < nLocal; j++) {
                int ni = elem[e * nLocal + i];
                int nj = elem[e * nLocal + j];
                /* Ajout réciproque des voisins */
                int exists = 0;
                for (int k = 0; k < adj[ni].count; k++) {
                    if (adj[ni].neighbors[k] == nj) { exists = 1; break; }
                }
                if (!exists) {
                    if (adj[ni].count == adj[ni].capacity) {
                        adj[ni].capacity *= 2;
                        adj[ni].neighbors = realloc(adj[ni].neighbors, adj[ni].capacity * sizeof(int));
                    }
                    adj[ni].neighbors[adj[ni].count++] = nj;
                }
                exists = 0;
                for (int k = 0; k < adj[nj].count; k++) {
                    if (adj[nj].neighbors[k] == ni) { exists = 1; break; }
                }
                if (!exists) {
                    if (adj[nj].count == adj[nj].capacity) {
                        adj[nj].capacity *= 2;
                        adj[nj].neighbors = realloc(adj[nj].neighbors, adj[nj].capacity * sizeof(int));
                    }
                    adj[nj].neighbors[adj[nj].count++] = ni;
                }
            }
        }
    }

    /* --- 2. Application de l'algorithme RCM --- */
    int *visited = calloc(nNodes, sizeof(int));
    int *queue = malloc(nNodes * sizeof(int));
    int queueHead, queueTail;
    int orderIndex = 0;
    for (int i = 0; i < nNodes; i++) {
        if (visited[i])
            continue;
        int start = i;
        for (int j = i; j < nNodes; j++) {
            if (!visited[j] && adj[j].count < adj[start].count)
                start = j;
        }
        queueHead = 0;
        queueTail = 0;
        visited[start] = 1;
        queue[queueTail++] = start;
        int *componentOrder = malloc(nNodes * sizeof(int));
        int compCount = 0;
        while (queueHead < queueTail) {
            int u = queue[queueHead++];
            componentOrder[compCount++] = u;
            int nUnvisited = 0;
            for (int k = 0; k < adj[u].count; k++) {
                int v = adj[u].neighbors[k];
                if (!visited[v])
                    nUnvisited++;
            }
            if (nUnvisited > 0) {
                int *unvisitedNeighbors = malloc(nUnvisited * sizeof(int));
                int idx = 0;
                for (int k = 0; k < adj[u].count; k++) {
                    int v = adj[u].neighbors[k];
                    if (!visited[v]) {
                        unvisitedNeighbors[idx++] = v;
                    }
                }
                /* Tri par degré croissant (simple bubble sort) */
                for (int a = 0; a < nUnvisited - 1; a++) {
                    for (int b = a + 1; b < nUnvisited; b++) {
                        if (adj[unvisitedNeighbors[a]].count > adj[unvisitedNeighbors[b]].count) {
                            int temp = unvisitedNeighbors[a];
                            unvisitedNeighbors[a] = unvisitedNeighbors[b];
                            unvisitedNeighbors[b] = temp;
                        }
                    }
                }
                for (int a = 0; a < nUnvisited; a++) {
                    int v = unvisitedNeighbors[a];
                    if (!visited[v]) {
                        visited[v] = 1;
                        queue[queueTail++] = v;
                    }
                }
                free(unvisitedNeighbors);
            }
        }
        /* Inversion de l'ordre */
        for (int a = 0; a < compCount / 2; a++) {
            int temp = componentOrder[a];
            componentOrder[a] = componentOrder[compCount - a - 1];
            componentOrder[compCount - a - 1] = temp;
        }
        for (int a = 0; a < compCount; a++) {
            rcm_order[orderIndex++] = componentOrder[a];
        }
        free(componentOrder);
    }
    free(queue);
    free(visited);

    /* --- 3. Mise à jour des coordonnées et connectivité --- */
    int *inv = malloc(nNodes * sizeof(int));
    for (int i = 0; i < nNodes; i++) {
        inv[rcm_order[i]] = i;
    }
    double *newX = malloc(nNodes * sizeof(double));
    double *newY = malloc(nNodes * sizeof(double));
    for (int i = 0; i < nNodes; i++) {
        newX[i] = theGeometry->theNodes->X[rcm_order[i]];
        newY[i] = theGeometry->theNodes->Y[rcm_order[i]];
    }
    free(theGeometry->theNodes->X);
    free(theGeometry->theNodes->Y);
    theGeometry->theNodes->X = newX;
    theGeometry->theNodes->Y = newY;
    if (theGeometry->theElements) {
        int nElem = theGeometry->theElements->nElem;
        int nLocal = theGeometry->theElements->nLocalNode;
        int *elem = theGeometry->theElements->elem;
        for (int i = 0; i < nElem * nLocal; i++) {
            int oldIndex = elem[i];
            elem[i] = inv[oldIndex];
        }
    }
    if (theGeometry->theEdges) {
        int nElem = theGeometry->theEdges->nElem;
        int nLocal = theGeometry->theEdges->nLocalNode;
        int *elem = theGeometry->theEdges->elem;
        for (int i = 0; i < nElem * nLocal; i++) {
            int oldIndex = elem[i];
            elem[i] = inv[oldIndex];
        }
    }
    for (int i = 0; i < nNodes; i++) {
        free(adj[i].neighbors);
    }
    free(adj);
    theGeometry->nodeOrder = rcm_order;
    theGeometry->nodeInv   = inv;
}


void femFullSystemConstrain(femFullSystem *mySystem, int myNode, double myValue) {
    double **A, *B;
    int i, size;

    A = mySystem->A;
    B = mySystem->B;
    size = mySystem->size;

    for (i = 0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0;
    }

    for (i = 0; i < size; i++)
        A[myNode][i] = 0;

    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}

#define MAX_ITER 1000

double *femSystemChebyschev(femFullSystem *theSystem, double l_max, double l_min) {
    double **A = theSystem->A;
    double *B = theSystem->B;
    int size = theSystem->size;
    
    double *x_k = calloc(size, sizeof(double));      // x^k (initialisé à zéro)
    double *x_k_prev = calloc(size, sizeof(double)); // x^{k-1}, initialisé à zéro
    double *r_k = malloc(size * sizeof(double));     // résidu courant
    
    double d = 0.5 * (l_max + l_min);
    double c = 0.5 * (l_max - l_min);
    double alpha_k = 1.0 / d;
    double alpha_k_prev;

    int k;
    for (k = 0; k < MAX_ITER; k++) {
        // Calcul du résidu r_k = B - A*x_k
        for (int i = 0; i < size; i++) {
            double Ax_i = 0.0;
            for (int j = 0; j < size; j++) {
                Ax_i += A[i][j] * x_k[j];
            }
            r_k[i] = B[i] - Ax_i;
        }

        // Critère d'arrêt (norme infinie du résidu)
        double res_norm = 0.0;
        for (int i = 0; i < size; i++)
            res_norm = fmax(res_norm, fabs(r_k[i]));

        if (res_norm < 1e-10) break;

        // Mise à jour de la solution
        double *x_k_next = malloc(size * sizeof(double));
        if (k == 0) {
            for (int i = 0; i < size; i++)
                x_k_next[i] = x_k[i] + alpha_k * r_k[i];
        } else {
            for (int i = 0; i < size; i++)
                x_k_next[i] = x_k[i] + alpha_k * r_k[i] + (alpha_k / alpha_k_prev - 1.0)*(x_k[i] - x_k_prev[i]);
        }

        // Mise à jour des alpha selon la récurrence de Chebyshev
        alpha_k_prev = alpha_k;
        alpha_k = 1.0 / (d - (c * c * alpha_k_prev) / 4.0);

        // Rotation des pointeurs pour prochaine itération
        free(x_k_prev);
        x_k_prev = x_k;
        x_k = x_k_next;
    }

    free(x_k_prev);
    free(r_k);
    return x_k; // solution finale (x_k)
}


#define MAX_IT_EIG 1000
#define TOL_EIG 1e-8
 
// dot product
double dot_vec(double *x, double *y, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += x[i] * y[i];
    return s;
}
 
/*
 * inverse_power_system
 *   Finds the eigenvalue of A closest to zero (λ of smallest magnitude),
 *   returning it in *mu (signed).  v[] is the eigenvector (length sys->size).
 */
int inverse_power_system(femFullSystem *sys, double *mu, double *v) {
    int n = sys->size;
    double mu_old = 0.0, norm;
    double *y = malloc(n * sizeof(double));
    if (!y) return 1;
 
    // 1) initialize v randomly and normalize
    for (int i = 0; i < n; i++) v[i] = (double)rand() / RAND_MAX;
    norm = sqrt(dot_vec(v, v, n));
    for (int i = 0; i < n; i++) v[i] /= norm;
 
    // 2) power‐inverse loop with shift=0
    for (int iter = 0; iter < MAX_IT_EIG; iter++) {
        // a) solve A y = v  via your FEM solver
        memcpy(sys->B, v, n * sizeof(double));
        double *sol = femFullSystemEliminate(sys);  // writes into sys->B
        // copy solution out
        for (int i = 0; i < n; i++) y[i] = sol[i];
 
        // b) normalize y → v
        norm = sqrt(dot_vec(y, y, n));
        if (norm < TOL_EIG) { free(y); return 1; }
        for (int i = 0; i < n; i++) v[i] = y[i] / norm;
 
        // c) Rayleigh quotient μ = vᵀ A v
        double Ray = 0.0;
        for (int i = 0; i < n; i++) {
            double row_dot = 0.0;
            for (int j = 0; j < n; j++)
                row_dot += sys->A[i][j] * v[j];
            Ray += v[i] * row_dot;
        }
        *mu = Ray;
 
        // d) convergence?
        if (iter > 0 && fabs(*mu - mu_old) < TOL_EIG * fabs(*mu))
            break;
        mu_old = *mu;
    }
 
    free(y);
    return 0;
}
 
/*
 * power_iteration_system
 *   Finds the eigenvalue of A of largest magnitude (λ_max),
 *   returning it in *mu.  v[] is the eigenvector.
 */
int power_iteration_system(femFullSystem *sys, double *mu, double *v) {
    int n = sys->size;
    double mu_old = 0.0, norm;
    double *y = malloc(n * sizeof(double));
    if (!y) return 1;
 
    // 1) initialize v randomly and normalize
    for (int i = 0; i < n; i++) v[i] = (double)rand() / RAND_MAX;
    norm = sqrt(dot_vec(v, v, n));
    for (int i = 0; i < n; i++) v[i] /= norm;
 
    // 2) standard power iteration
    for (int iter = 0; iter < MAX_IT_EIG; iter++) {
        // a) y = A v
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++)
                s += sys->A[i][j] * v[j];
            y[i] = s;
        }
 
        // b) normalize y → v
        norm = sqrt(dot_vec(y, y, n));
        if (norm < TOL_EIG) { free(y); return 1; }
        for (int i = 0; i < n; i++) v[i] = y[i] / norm;
 
        // c) Rayleigh quotient μ = vᵀ A v
        *mu = dot_vec(v, y, n);
 
        // d) convergence?
        if (iter > 0 && fabs(*mu - mu_old) < TOL_EIG * fabs(*mu))
            break;
        mu_old = *mu;
    }
 
    free(y);
    return 0;
}
