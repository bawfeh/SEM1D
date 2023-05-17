class GaussDistribution:

    # class variables

    tol = 1e-8
    epsilon = 2.220446049250313e-16 
    import numpy as np

    # Constructor
    def __init__(self,n = 0): # with numpy module
        self.n = int(n)        # polynomial degree, integer n >= 0
        self.ngll = self.n+1   # number of GLL points (min 2)
        self.mgl = self.n      # number of GL points required (min 1)
        if (n <= 0):
            self.ngll = 2
            self.mgl = 1
        self.glnodes = self.np.zeros((self.mgl,))
        self.glweights = self.np.zeros((self.mgl,))
        self.gllnodes = self.np.zeros((self.ngll,))
        self.gllweights = self.np.zeros((self.ngll,))
        self.GaussLobattoLegendre()
        self.N = 0
        self.M = 0
        self.ldiag = self.np.array([])
        self.ldeno = self.LagrangeBasisDenominator()# nodes = gll
        

    def LegendreDerivative(self, x, n = -1, Leg=False): 
        """ 1st derivative of n^th Legendre polynomial at x; 
        x = numpy array """

        if isinstance(x,float):
            if (n < 0): n = self.n
            L1 = 1.; Ld = 0.
            if n==1: 
                Ld = 1.; L1 = x
            elif n>1:
                dx = 1-x**2
                if abs(dx)<self.epsilon:
                    Ld = x**(n-1)*(.5*n*(n+1))
                else: # Bonnets recursion
                    L1,L0 = self.Legendre(x,n,True)
                    Ld = n/dx * (L0 - x*L1)
    
            if (Leg): return Ld,L1 
            return Ld

        m = len(x)
        if (n < 0): n = self.n
        L1 = self.np.ones((m,))
        Ld = self.np.zeros((m,))
        if n==1: 
            Ld = self.np.ones((m,)); L1 *= x
        elif n>1:
            L1,L0 = self.Legendre(x,n,True)
            for j in range(m):
                dx = 1-x[j]**2
                if abs(dx)<self.epsilon:
                    Ld[j] = x[j]**(n-1)*(.5*n*(n+1))
                else: # Bonnets recursion
                    Ld[j] = n/dx * (L0[j] - x[j]*L1[j])

        if (Leg): return Ld,L1 
        return Ld

    def Legendre(self, x, n = -1, prev=False): # inputs pointer to np array x
        """ n^th Legendre polynomial at x ( n = deg) """
        
        if (n < 0): n = self.n
        L0 = 1. if isinstance(x,float) else self.np.ones((len(x),))
        L1 = L0 if (n==0) else x*L0
        if (n > 1):
            # three-step [Bonnets] recurrence
            for i in self.np.linspace(1,n-1,n-1):
                L2 = ((2*i+1)/(i+1))*(x*L1) - (i/(i+1))*L0 # update
                L0 = L1; L1 = L2  # move pointers

        if prev: return L1, L0
        return L1

    def GaussLegendre(self):

        if (self.mgl <= 1):
            self.glweights[0] = 0.5
            return None

        np = self.mgl
        A = self.np.zeros((np,np))
        A[0][1] = 1.; k = 1
        for  i in self.np.linspace(2,np-1,np-2):
            A[k][k-1] = (i-1)/(2*i-1)
            A[k][k+1] = i/(2*i-1)
            k += 1
        A[-1][-2] = (np-1.)/(2*np-1)

        eig, _ = self.np.linalg.eig(A)
        eig = self.np.real(eig)
        idx = self.np.argsort(eig)
        self.glnodes = eig[idx]

        Ld = self.LegendreDerivative(self.glnodes, np)
        for i in range(self.mgl):
            self.glweights[i] = 2 / ((1-self.glnodes[i]**2) * (Ld[i]**2))

        return None

    def GaussLobattoLegendre(self):

        fac = (self.ngll-1) * self.ngll
        self.gllnodes[0] = -1.; self.gllnodes[-1] = 1.

        if (self.ngll >= 3):
            self.GaussLegendre()
            # suitable starting values for a Newton-iteration procedure
            gllnodes = (self.glnodes[0:-1]+self.glnodes[1:])/2

            for i in range(len(gllnodes)): #range(1,self.ngll-1):
                x = 0.
                while (abs(gllnodes[i]-x) > self.tol ):
                    x = gllnodes[i]
                    Ld, L = self.LegendreDerivative( x, Leg=True )
                    gllnodes[i] = x + ((1-x**2)*Ld) / (fac*L)
            self.gllnodes[1:-1] = gllnodes

        # This loop computes the associated weights
        Ln = self.Legendre (self.gllnodes)
        for i in range(self.ngll):
            self.gllweights[i] = 2 / (fac*Ln[i]**2)

        return None
                

    def LagrangeDerivativeMatrix(self, x=None):
        
        if (x is None): x = self.gllnodes
        N = len(x)
        D = self.np.zeros((N,N))

        if (N == 1): return D 
        
        Ln = self.Legendre( x )

        for i in range(N):
            for j in range(N):
                if (i != j): D[i][j] = Ln[i] / ( Ln[j] * (x[i] - x[j]))
        fac = float(N*(N-1)) / 4
        D[0][0] = -fac
        D[-1][-1] = fac

        return D
                

    def Lagrange2ndDerivativeMatrix(self, x=None):
        
        if (x is None): x = self.gllnodes
        N = len(x)
        D = self.np.zeros((N,N))

        if (N < 3): return D
        
        Ln = self.Legendre( x )
        dLn = self.LegendreDerivative( x )
#         dLn, Ln = self.LegendreDerivative( x, Leg=True)
        fac = N*(N-1)
        for j in range(1,N-1):
            for l in range(1,N-1):
                if (l == j):
                    D[l,l] = (2*x[l]*dLn[l]-fac*Ln[l]) / (3*Ln[l]*(1-x[l]**2))
                else:
                    D[j,l] = -2*Ln[j] / (Ln[l]*(x[j]-x[l])**2)
        for l in range(1,self.n):
            D[0,l] = ((-1)**(N-1) * fac*(1+x[l])-4) / (2*Ln[l]*(1+x[l])**2)#
            D[-1,l] = (fac*(1-x[l])-4) / (2*Ln[l]*(1-x[l])**2)
        l = -1
        D[l,0] = (fac*(1-x[0])-4) / (2*Ln[0]*(1-x[0])**2)
        D[0,l] = ((-1)**(N-1) * fac*(1+x[l])-4) / (2*Ln[l]*(1+x[l])**2)#
        for j in range(1,self.n):
            D[j,0] = -2*Ln[j] / (Ln[0]*(x[j]-x[0])**2)
            D[j,l] = -2*Ln[j] / (Ln[l]*(x[j]-x[l])**2)
        D[0,0] = fac * (fac-2.) / 24
        D[l,l] = D[0,0]

        return D

    def LagrangeBasisDenominator(self, glnodes = None):

        if (glnodes is None) : glnodes = self.gllnodes
        N = len(glnodes)
        X = self.np.tile(glnodes, (N,1))
        X = X.T - X
        self.np.fill_diagonal(X,1.)
        
        return self.np.prod(X,1)

    def LagrangeBasisMatrix(self, x, glnodes = None):
        
        if (glnodes is None) : glnodes = self.gllnodes
        N = len(glnodes); M = len(x)
        Z = self.np.ones((N,N,M))
        X = self.np.transpose(Z*x,(2,0,1)) - glnodes
        self.filldiagonal(self.np.shape(X))
        X[self.ldiag[:M]] = 1.
        
        if (glnodes is None):
            return self.np.prod(X,2) / self.ldeno 

        P = glnodes - self.np.tile(glnodes, (N,1)).T
        self.np.fill_diagonal(P,1.)
        return self.np.prod(X,2) / self.np.prod(P,0)

    def LagrangeBasisMatrix2(self, x, glnodes = None):
        
        if (glnodes is None) : 
            return self.LagrangeBasisMatrix(x)
        N = len(glnodes); M = len(x)
        P = glnodes - self.np.tile(glnodes, (N,1)).T
        Q, R = self.np.meshgrid(x, glnodes)
        Q = Q - R
        rgx = self.np.arange(1,N)
        for j in range(N):
            R[j] = self.np.prod(Q[rgx],0) / self.np.prod(P[rgx,j])
            rgx = self.np.mod(rgx+1,N)
        
        return R.T

    def filldiagonal(self, shape):
        if (shape[0]>self.M):
            X = self.np.eye(shape[1], dtype=bool).tolist()
            ldiag = self.ldiag.tolist()
            for i in range(self.M,shape[0]): ldiag.append(X)
            self.M = shape[0]
            self.ldiag = self.np.array(ldiag)
        return None

#===== Trials ================================

    def FirstPolynomialDerivativeMatrix(self, x=None):
        """
        PolynomialDerivativeMatrix(): Assemble the first Polynomial Derivative Matrix using matrix multiplication.
    
        Syntax:
            ``D = FirstPolynomialDerivativeMatrix(x)``
    
        Input:
            * ``x`` = (1d-array,float) set of points on which to evaluate the derivative matrix
    
        Output:
            * ``D`` = derivative matrix
    
        Notes:
            Algorithm (37) from :cite:`Kopriva2009`
        """
        if (x is None): x = self.gllnodes
        N = x.shape[0]
        w = self.BarycentricWeights(x)
        D = self.np.zeros((N,N))
        for i in range(0,N):
            for j in range(0,N):
                if (j != i):
                    D[i,j] = w[j]/w[i] * 1./(x[i] - x[j])
                    D[i,i] = D[i,i] - D[i,j]
        return D


    def PolynomialDerivativeMatrix(self, x=None, k=1):
        """
        PolynomialDerivativeMatrix(): Assemble the Polynomial ``k``-th Derivative Matrix using the matrix recursion. This algorithm is generic for every types of polynomials.
    
        Syntax:
            ``D = PolynomialDerivativeMatrix(x,k)``
    
        Input:
            * ``x`` = (1d-array,float) set of points on which to evaluate the derivative matrix
            * ``k`` = derivative order
    
        Output:
            * ``D`` = ``k``-th derivative matrix
    
        Notes:
            Algorithm (38) taken from :cite:`Kopriva2009`
        """
        if (x is None): x = self.gllnodes
        N = x.shape[0]
        w = self.BarycentricWeights(x)
        D = self.np.zeros((N,N,k))
        D[:,:,0] = self.FirstPolynomialDerivativeMatrix(x)
        if ( k == 1 ): return D[:,:,k-1]
        for m in range(2,k+1):
            for i in range(0,N):
                for j in range(0,N):
                    if ( j != i ):
                        D[i,j,m-1] = m/(x[i]-x[j]) * ( w[j]/w[i] * D[i,i,m-2] - D[i,j,m-2])
                        D[i,i,m-1] = D[i,i,m-1] - D[i,j,m-1]
        return D[:,:,k-1]


    def BarycentricWeights(self, x):
        """
        BarycentricWeights(): Returns a 1-d array of weights for Lagrange Interpolation
    
        Syntax:
            ``w = BarycentricWeights(x)``
    
        Input:
            * ``x`` = (1d-array,float) set of points
    
        Output:
            * ``w`` = (1d-array,float) set of barycentric weights
    
        Notes:
            Algorithm (30) from :cite:`Kopriva2009`
        """
        N = x.shape[0]
        w = self.np.zeros((N))
        for j in range(0,N):
            w[j] = 1.
        for j in range(1,N):
            for k in range(0,j):
                w[k] = w[k] * (x[k] - x[j])
                w[j] = w[j] * (x[j] - x[k])
        for j in range(0,N):
            w[j] = 1. / w[j]
        return w



    def LagrangeInterpolationMatrix(self, x, w, xi):
        """
        LagrangeInterpolationMatrix(): constructs the Lagrange Interpolation Matrix from points ``x`` to points ``xi``
    
        Syntax:
            ``T = LagrangeInterpolationMatrix(x, w, xi)``
    
        Input:
            * ``x`` = (1d-array,float) set of ``N`` original points
            * ``w`` = (1d-array,float) set of ``N`` barycentric weights
            * ``xi`` = (1d-array,float) set of ``M`` interpolating points
    
        Output:
            * ``T`` = (2d-array(``MxN``),float) Lagrange Interpolation Matrix
    
        Notes:
            Algorithm (32) from :cite:`Kopriva2009`
        """
        N = x.shape[0]
        M = xi.shape[0]
        T = self.np.zeros((M,N))
        for k in range(0,M):
            rowHasMatch = False
            for j in range(0,N):
                T[k,j] = 0.
                if self.np.isclose(xi[k],x[j]):
                    rowHasMatch = True
                    T[k,j] = 1.
            if (rowHasMatch == False):
                s = 0.
                for j in range(0,N):
                    t = w[j] / (xi[k] - x[j])
                    T[k,j] = t
                    s = s + t
                for j in range(0,N):
                    T[k,j] = T[k,j] / s
        return T


    def LagrangeInterpolate(self, x, f, xi):
        """
        LagrangeInterpolate(): Interpolate function values ``f`` from points ``x`` to points ``xi`` using Lagrange weights
    
        Syntax:
            ``fi = LagrangeInterpolate(x, f, xi)``
    
        Input:
            * ``x`` = (1d-array,float) set of ``N`` original points where ``f`` is evaluated
            * ``f`` = (1d-array/2d-array,float) set of ``N`` function values (if K functions are passed, the values are stored in a NxK matrix)
            * ``xi`` = (1d-array,float) set of ``M`` points where the function is interpolated
    
        Output:
            * ``fi`` = (1d-array,float) set of ``M`` function values (if K functions are passed, the values are stored in a MxK matrix)
    
        Notes:
            Modification of Algorithm (33) from :cite:`Kopriva2009`
        """
        # Obtain barycentric weights
        w = self.BarycentricWeights(x)
        # Generate the Interpolation matrix
        T = self.LagrangeInterpolationMatrix(x, w, xi)
        # Compute interpolated values
        fi = self.np.dot(T,f)
        return fi
# ============================================


class Mesh:

    # class variables
    tol = 1e-8
    import numpy as np

    # Constructor
#     def __init__(self,np, data, gnodes): # with numpy module
    def __init__(self, data, gll, constructMesh = True): # with Gauss distribution
        self.n = gll.mgl
        self.N = data.polyDeg+1
        self.data = data
        self.anyPeriodic = data.leftrightbc
        self.anyDirbc = (data.orientDirbc['left'] or data.orientDirbc['right']) or data.orientDirbc['all']
        self.anyNeumbc = data.orientNeumbc['left'] or data.orientNeumbc['right'] 
        self.dG0 = (data.polyDeg == 0)

        self.axesflag = None
        self.elementsflag = None
        self.periodicflag = None
        self.bcflag = None

        self.flags = {\
            'edges2' : True,\
            'axes' : True,\
            'elements' : True,\
            'elements2' : True,\
            'bc' : True,\
            'periodic' : True }

        self.edges2 = None
        self.Ned = None
        self.Ne = None
        self.mtol = None # mesh tol
        self.Nbed = None
        self.neighbours2 = None
        self.neighbours = None
        self.Bedges = None
        self.nodes = None
        self.elements = None   # global numbering element by element
        self.elementsL = None  # lobal numbering element by element
        self.edges = None
        self.Dnodes = []
        self.Nnodes = []
        self.Pnodes = []
        self.element_centres = None
        self.normal = None
        self.numNodes = None
        self.element_sizes = None
        self.xL = None

        if constructMesh : self.getAlldata(gll.gllnodes)

    def getAlldata(self,gnodes):
        self.flags['axes'] = self.createAxes(gnodes)
        self.flags['elements'] = self.fillNodeElementEdge()
        self.flags['bc'] = self.setboundaryConditions()
        self.processMesh()
#         self.mtol = max(self.tol / self.N, 1.0e-16)
        self.mtol = 2.220446049250313e-16
        self.xL = None # release this pointer
#         self.np = None # release this pointer
        self.data = None # release this pointer

        return None

    def createAxes(self,nodes):
        self.Ne = self.data.Ne;
        Ne = self.Ne
        xL = self.np.tile(nodes,(Ne,1))
        self.element_sizes = self.np.zeros((Ne,))
        nx = self.data.nx
        l = 0;

        for i in range(len(nx)):
            self.element_sizes[l:l+nx[i]] = self.data.dx[i]
            l += nx[i]

        boundaryvals = self.np.cumsum(self.element_sizes)
        boundaryvals[1:] = boundaryvals[0:-1]
        boundaryvals[0] = 0.
        boundaryvals += self.data.xb[0]

        for el in range(Ne):
            xL[el] = boundaryvals[el] + (self.element_sizes[el]/2) * (xL[el]+1.)
        self.xL = xL

        self.data.xb[1] = xL[-1][-1]
        print("Domain: x-axis (", self.np.around(self.data.xb[0]),", ", self.np.around(self.data.xb[1]),")")

        return False

    def edgeConnectivity(self):
        # assign unique global numbering to element edges (edge connectivity matrix)
        Ne = self.Ne
        edges2 = self.np.zeros((Ne,2),dtype=int) #  columns W,E
    
        # number all LEFT edges of elements in the order left < right
        l = 0;
        for el in range(Ne):
            edges2[el][0] = l; l += 1
    
        # number all RIGHT edges of elements in the order left < right
    
        for el in range(Ne-1): edges2[el][1] = edges2[el+1][0]

        edges2[-1][1] = l; l += 1
        Ned = l #  total number of edges
        
        neighbours2 = self.np.zeros((Ned,2),dtype=int)  #  W,E
        for  el in range(Ne):
            neighbours2[edges2[el][0]][1] = el;
            neighbours2[edges2[el][1]][0] = el;
        
        if (self.data.leftrightbc):
            neighbours2[edges2[0][0]][0] =  Ne-1 
            neighbours2[edges2[-1][1]][1] = 0    
        else:
            neighbours2[edges2[0][0]][0] =  -1 #  out of bounds
            neighbours2[edges2[-1][1]][1] = -1 #  out of bounds
        
        self.edges2 = edges2
        self.Ned = Ned
        self.neighbours2 = neighbours2
        
        return False

    def trackNeighbours(self):

        # keep record of the neighbouring elements to each given element
        # use -1 to denote imaginary (ghost) element neighbours

        Ne = self.Ne
        neighbours = self.np.zeros((Ne,2),dtype=int) # columns W,E
        Bedges = self.np.zeros((2,3),dtype=int)
    
        if (self.flags['edges2']): self.flags['edges2'] = self.edgeConnectivity()
    
        l = 0;
        for el in range(Ne):
    
            if (el > 0):  # not west boundary el==0
                neighbours[el][0] = el-1
            else:
                Bedges[l][0] = self.edges2[el][0]  # boundary edge numbering
                Bedges[l][1] = 0             # orientation (W,E)
                Bedges[l][2] = el; l += 1     # associated element
                if (self.data.leftrightbc): 
                    neighbours[el][0] = Ne-1
                else: 
                    neighbours[el][1] = -1 # no neighbour on the west of this element
            if (Ne > el+1):  # not east boundary el==Ne-1
                neighbours[el][1] = el+1
            else:
                Bedges[l][0] = self.edges2[el][1]
                Bedges[l][1] = 1
                Bedges[l][2] = el; l += 1
                if (self.data.leftrightbc): 
                    neighbours[el][1] = 0
                else: 
                    neighbours[el][1] = -1 
        self.Nbed = l
        self.neighbours = neighbours 
        self.Bedges = Bedges 
    
        return False

    def unique(self,x): # returns unique elements of a 1D sorted array
        lnext = x[1:].tolist()+[x[0]]
        return x[abs(self.np.array(lnext) - x) > self.tol]

    def fillNodeElementEdge(self):
        # count the nodes in each element and edge
        # collect the unique coordinates for all nodes in the mesh
        if (self.flags['edges2']): self.flags['edges2'] = self.edgeConnectivity()
        if (self.flags['elements2']): self.flags['elements2'] = self.trackNeighbours()
        ind = self.np.argsort(self.xL, axis=None)
        y = self.np.sort(self.xL, axis=None)
        self.nodes = self.unique(y)
        self.numNodes = len(self.nodes)
        N = len(ind)
        rank = self.np.arange(N)
        rank[0] = 0; l = 0
        for i in range(1,N):
            if abs(y[i]-y[i-1]) < self.tol: rank[i] = rank[i-1]
            else: l += 1; rank[i] = l
        lind = self.np.arange(N)
        lind[ind] = rank
        self.elements = self.np.reshape(lind,self.np.shape(self.xL))
        self.elementsL = self.np.reshape(self.np.arange(N),self.np.shape(self.xL))
        Ne = self.Ne
        edges = self.np.zeros((self.Ned,),dtype=int)
        for el in range(Ne):
            edges[self.edges2[el][0]] = self.elements[el][0]
            edges[self.edges2[el][-1]] = self.elements[el][-1]
        self.edges = edges

        return False

    def setboundaryConditions(self):
        # set boundary conditions locations
        if (self.data.leftrightbc): 
            self.Pnodes = self.edges[self.Bedges[:,0].tolist()]
            self.flags['periodic'] = False
            return False
    
        if (self.data.orientDirbc['left'] or self.data.orientDirbc['right']):
            if (self.data.orientDirbc['all']):
                self.Dnodes = self.edges[self.Bedges[:,0].tolist()]
                print("Dirichlet bc: ALL")
                return False
            else:
                for  l in range(self.Nbed):
                    if (self.Bedges[l][1]==0 and self.data.orientDirbc['left']):
                        self.Dnodes.append(self.edges[self.Bedges[l][0]])
                        print("Dirichlet bc: LEFT")
                    elif (self.Bedges[l][1]==1 and self.data.orientDirbc['right']):
                        self.Dnodes.append(self.edges[self.Bedges[l][0]])
                        print("Dirichlet bc: RIGHT")
                return False

#         if (len(self.data.Dirbc)>0): # for 2D or 3D domains
#             for edge in  self.data.Dirbc:
#                 if (self.BedgeLoc(edge)): self.Dnodes.append(self.edges[self.Bedges[edge][0]])
#             self.Dnodes = np.unique(self.Dnodes)
    
        if (self.data.orientNeumbc['left'] or self.data.orientNeumbc['right']):
            for  l in range(self.Nbed):
                if (self.Bedges[l][1]==0 and self.data.orientNeumbc['left']):
                    self.Nnodes.append(self.edges[self.Bedges[l][0]])
                    print("Neumann bc: LEFT")
                elif (self.Bedges[l][1]==1 and self.data.orientNeumbc['right']):
                    self.Nnodes.append(self.edges[self.Bedges[l][0]])
                    print("Neumann bc: RIGHT")
            return False
#         if (len(self.data.Neumbc)>0): # for 2D or 3D domains
#             for edge in self.data.Neumbc:
#                 if (self.BedgeLoc(edge)): self.Nnodes.append(self.edges[self.Bedges[edge][0]])
#             self.Nnodes = np.unique(self.Nnodes)
    
        return False

#     def BedgeLoc(self,l):
#         l = 0; rr = self.Nbed
#         while (ll<rr):
#             if (ll==l) :
#                 l = ll; return True
#             if (rr == l):
#                 l = rr; return True
#             ll += 1; rr -= 1
#         print("Edge ", l , " is not a valid boundary edge!")
#         return False

    def  isBedge(self,l):
        ll = 0; rr = self.Nbed-1
        while (ll<=rr) :
            if (self.Bedges[ll][0] == l): return True
            if (self.Bedges[rr][0] == l): return True
            ll += 1; rr -= 1
        return False

    def processMesh(self):
        # process each element of mesh,
        # to obtain element centres & volumes, edge/interface areas/length
        if ( self.flags['elements'] ) : self.flags['elements'] = self.fillNodeElementEdge();
        Ne = self.Ne
        self.element_centres = self.np.zeros((Ne,))
        v = self.np.zeros((2,),dtype=int)
        for el in range(Ne):
            v[0] = self.edges[self.edges2[el][0]]
            v[1] = self.edges[self.edges2[el][1]]
            self.element_centres[el] = (self.nodes[v[0]] + self.nodes[v[1]])/2;
        self.normal = self.np.array([-1.,1.])

        return None

    def grid(self):

        from matplotlib import pyplot as plt
        from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

        fig, ax = plt.subplots()
        x = self.np.unique(self.nodes)
        ax.plot(x, self.np.zeros(self.np.shape(x)), 'b.-')

        ax.set_xlim(min(x), max(x))
        ax.set_ylim(-0.5, 0.5)
        ax.xaxis.set_major_locator(MultipleLocator(.5))
        ax.xaxis.set_minor_locator(MultipleLocator(.25))
        plt.xlabel('x')
        plt.show()

    def plot(self, u):
        import matplotlib.pyplot as plt
        from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
        
        fig, ax = plt.subplots()
        
        h = ax.plot(self.nodes, u, linewidth=2)

        ax.set_xlim(min(self.nodes), max(self.nodes))
#         ax.set_ylim(min(u), max(u))
#         ax.yaxis.set_major_locator(MultipleLocator(.5))
#         ax.yaxis.set_minor_locator(MultipleLocator(.25))
        ax.xaxis.set_major_locator(MultipleLocator(.5))
        ax.xaxis.set_minor_locator(MultipleLocator(.25))
        ax.set_xlabel('x')
        ax.set_ylabel('u')
        
#         ax.grid(True)
        plt.show()

#     def inclusion(self, x, el):
#         xe = self.nodes[self.edges[self.edges2[el][0]]] # west = 0
#         boolb = abs(x-xe) < self.tol
#         boolW = ((x-xe) * self.normal[0] < self.tol ) | boolb
#         xe = self.nodes[self.edges[self.edges2[el][1]]] # east = 1
#         boolb = abs(x-xe) < self.tol
#         boolE = ((x-xe) * self.normal[1] < self.tol ) | boolb
#         return ( boolW & boolE )

    def inclusion(self, x, el):
        s = self.nodes[self.edges[self.edges2[el]]]
        if isinstance(x, float):
            return ( (x >= s[0]) and (x <= s[1]) )
        return ( (x >= s[0]) & (x <= s[1]) )

    def contains(self, x):
        """ find which element contains x """
        if isinstance(x, float):
            xe = self.nodes[self.edges[self.edges2[:,0]]] # west = 0
            return self.np.sum(xe <= x) - 1
        xe = self.nodes[self.edges[self.edges2[:,0]]] # west = 0
        return self.np.array([self.np.sum(xe <= s)-1 for s in x])

    def getNe(self):  
        return self.Ne
    def getN(self):
        return self.N
    def getnEdges(self) :
        return self.Ned
    def getnNodes(self) :
        return self.numNodes
    def getnBedges(self) :
        return self.Nbed
    def isValid(self,el):
        return (el>=0 and el<self.Ne)
    def centre(self, el) :
        return self.element_centres[el]
    def size(self,el) :
        return self.element_sizes[el]

# ==============================================================

class MeshDG (Mesh):  # derived class from Mesh class


    # Constructor
    def __init__(self, data, gll): # with numpy module
        super().__init__( data, gll, False)
        if (self.dG0): 
            self.n += 1; self.N += 1  # treat like n == 1

        self.getAlldata(gll.gllnodes)

    def edgeConnectivity(self):
        # assign unique global numbering to element edges (edge connectivity matrix)
        Ne = self.Ne
        edges2 = self.np.reshape(self.np.arange(2*Ne),(Ne,2))
        Ned = Ne*2
        
        neighbours2 = self.np.zeros((Ned,2),dtype=int)  #  columns W,E
        rangex = self.np.arange(Ne)
        neighbours2[edges2[:,0],1] = rangex
        neighbours2[edges2[:,1],0] = rangex
        neighbours2[edges2[:,0],0] = rangex-1
        neighbours2[edges2[:,1],1] = rangex+1
        
        if (self.data.leftrightbc):
            neighbours2[edges2[0][0]][0] =  Ne-1 
            neighbours2[edges2[-1][1]][1] = 0    
        else:
            neighbours2[edges2[0][0]][0] =  -1 #  out of bounds
            neighbours2[edges2[-1][1]][1] = -1 #  out of bounds
        
        self.edges2 = edges2
        self.Ned = Ned
        self.neighbours2 = neighbours2
        
        return False

    def fillNodeElementEdge(self):
        # count the nodes in each element and edge
        # collect the unique coordinates for all nodes in the mesh
        if (self.flags['edges2']): self.flags['edges2'] = self.edgeConnectivity()
        if (self.flags['elements2']): self.flags['elements2'] = self.trackNeighbours()
        self.nodes = self.np.reshape(self.xL, self.xL.size)
        self.numNodes = len(self.nodes)
        N = len(self.nodes)
        lind = self.np.arange(N)
        self.elements = self.np.reshape(lind,self.np.shape(self.xL))
        Ne = self.Ne
        edges = self.np.zeros((self.Ned,),dtype=int)
        edges[self.edges2[:,0]] = self.elements[:,0]
        edges[self.edges2[:,-1]] = self.elements[:,-1]
        self.edges = edges
        self.elementsL = self.elements

        return False
