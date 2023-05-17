class Assembly:

    tol = 1e-8
    epsilon = 2.220446049250313e-16 
        
    from scipy.sparse import coo_matrix
    from scipy.sparse import linalg as linalg
    

    def __init__(self, mesh, gll, Re):
        self.mesh = mesh
        self.gllweights = gll.gllweights
        self.nvglobal = mesh.getnNodes()
        self.nv = 0
        self.Re = Re
        for el in mesh.elements: self.nv += len(el)
        self.dfac = 1.
        self.mfac = 1.
        self.rfac = False  # rfac = True, for reaction-diffusion
        self.np = gll.np
        self.nDnodes = len(mesh.Dnodes)
        self.bc = []
        
        self.Vgl = self.np.zeros((self.nvglobal,))
        self.D = gll.LagrangeDerivativeMatrix()
        self.A = Re * (self.np.dot(self.D.T, (self.D.T * gll.gllweights).T))
        self.B = self.np.zeros(self.np.shape(self.D))
        di = self.np.diag_indices(len(self.B))
        self.B[di] += gll.gllweights
        self.Dmat = None  # diffusion matrix
        self.Amat = None  # stiffness matrix
        self.lu = None  # sparse lu factors of stiffness matrix
        self.r = None
        self.c = None
        self.bvec = None  # assembled lhs vector
        self.Mmat = None  # mass matrix
        self.Cmat = None  # convection matrix
        self.MmatDiag = None  # mass matrix as diagonal
        self.mvec = None  # assembled lhs vector
        self.bcontr = None # boundary contrition to lhs vector
        self.nDof = 0
        self.maskvec = [] # boolean to select all Dirichlet nodes from global vectors
        self.maskvecL = [] # boolean to select all Dirichlet nodes from local vectors
        self.uvec = None # boolean to select all Dof nodes from global vectors
        self.gvec = None # boolean to select all global nodes from local vectors
        self.MassLhs = True
        self.maskvector()
        


    def maskvector(self):
        if (self.nDof==0):
            self.nDof = self.nvglobal - len(self.mesh.Dnodes)
            self.maskvec  = self.np.zeros((self.nvglobal,), dtype=bool) # boolean vector
            for j in self.mesh.Dnodes: 
                self.maskvec[j] = True
            self.uvec = (self.maskvec == False)

            if (self.mesh.anyPeriodic):
                self.nDof -= 1
                self.uvec[self.mesh.Pnodes[1]] = False
            if (self.nDnodes>0): self.maskvecL = self.spread(self.maskvec)
            else: self.maskvec = []

            self.gvec = self.np.ones((self.nv,), dtype=bool) 
            gvec = self.np.zeros((self.nvglobal,), dtype=bool) 
            l = 0
            for el in self.mesh.elementsL:
                v = ~gvec[self.mesh.elements[l]] & self.gvec[el]
                self.gvec[el] = v
                gvec[self.mesh.elements[l]] = True; l += 1

        return None

    def zeros(self, n, u=None):
        if (u is None): return self.np.zeros((n,))
        elif (len(u)<n): return self.np.zeros((n,))
        u[:] = 0
        return u

    def setbc(self,bc):
        self.bc = bc
        return None

    def set_massfactor(self, mfac):
        self.mfac = mfac
        return None

    def set_difffactor(self, dfac):
        self.dfac = dfac
        return None

    def Scatter(self, ug, mesh):
        # local to global SEM map // scatter (copy) operator
        l = 0; u = self.zeros(self.np.prod(self.np.shape(mesh.elements)))
        for el in mesh.elements:
            u[mesh.elementsL[l]] = ug[el]; l += 1
        return u

    def scatter(self, ug):
        # local to global SEM map // scatter (copy) operator
        l = 0; u = self.zeros(self.nv)
        for el in self.mesh.elements:
            u[self.mesh.elementsL[l]] = ug[el]; l += 1
        return u

    def spread(self, ug): # ug is a global boolean array
        # local to global SEM map // scatter (copy) operator
        if (ug is None): return None
        l = 0; u = self.np.zeros((self.nv,), dtype=bool)
        for el in self.mesh.elements:
            u[self.mesh.elementsL[l]] = ug[el]; l += 1
        return u

    def Gather(self, u, mesh, add=True, mask=True):
        # local to global SEM map: gather (add) operator
        l = 0; ug = self.zeros(mesh.getnNodes())
        if add :
            for el in self.mesh.elements:
                ug[el] += u[mesh.elementsL[l]]; l += 1 
            self.periodic(ug)
        else :
            ug = u[self.gvec]
            self.periodic(ug,2.)
        if (mask): ug[self.maskvec] = 0.
        return ug


    def gather(self, u, add=True, mask=True):
        # local to global SEM map: gather (add) operator
        l = 0; ug = self.zeros(self.nvglobal)
        if add :
            for el in self.mesh.elements:
                ug[el] += u[self.mesh.elementsL[l]]; l += 1 
            self.periodic(ug)
        else :
            ug = u[self.gvec]
            self.periodic(ug,2.)
        if (mask): ug[self.maskvec] = 0.
        return ug

    def mask(self, u):  # remove Dirichlet Nodes from local vector
        # enforce Dirichlet boundary conditions (if any)
        if (self.mesh.anyDirbc):
            if (len(u)==self.nv): u[self.maskvecL] = 0.
            else : u[self.maskvec] = 0.
        return None

    def periodic(self, u, dfac = 1.): 
        # enforce periodic boundary conditions (if any), on lobal vector
        if (self.mesh.anyPeriodic):
            if (len(u)==self.nv): 
                self.Vgl = self.periodic(self.gather(u,False,False),dfac)
                return self.scatter(self.Vgl)
            else:
                l = self.mesh.Pnodes[0]; r = self.mesh.Pnodes[1]
                u[l] = (u[l] + u[r]) / dfac
                u[r] = u[l]
        return u

    def enforceBc(self, u): # on local vector
        u = self.periodic(u)
        self.mask(u)
        return u

    def dssum(self, u, mask = True):
        # Direct stiffness summation: scatter[gather]
        self.Vgl = self.gather(u,True,mask)
        return self.scatter(self.Vgl)

    def Mass(self, u, mask = True):
        l = 0
        for el in self.mesh.elementsL:
            dx = self.mesh.size(l) / 2; l += 1
            u[el] *= ((dx*self.mfac) * self.gllweights)
        if (mask): self.mask(u)
        return None

    def MassInv(self, u, assembleMatrix = False, mask = True):
        w = self.np.zeros((self.nvglobal,)); l = 0
        for el in self.mesh.elements:
            dx  = self.mesh.size(l) / 2; l += 1
            w[el] += (dx) * self.gllweights
        b = self.gather(u,False,mask)
        b[self.uvec] /= w[self.uvec]
        return self.scatter(b)

#     def MassInv(self, u, assembleMatrix = False, mask = True):
#         self.assembleMassMatrix(False)
#         b = self.gather(u,False,mask)
#         ut = self.linalg.spsolve(self.Mmat, b[self.uvec])
#         b[self.uvec] = ut
#         return self.scatter(b)

    def assembleMassMatrixDiag(self, assemblyMat):
        if (assemblyMat or (self.MmatDiag is None)):
            print("Assembling mass matrix...")
            w = self.np.zeros((self.nvglobal,)); l = 0
            for el in self.mesh.elements:
                dx  = self.mesh.size(l) / 2; l += 1
                w[el] += (dx) * self.gllweights
            w = self.periodic(w)
            self.MmatDiag = w[self.uvec]
        return None

#     def MassInv(self, u, assembleMatrix = False, mask = True):
#         self.assembleMassMatrixDiag(assembleMatrix)
#         b = self.gather(u, False, mask)
#         b[self.uvec] /= self.MmatDiag
#         return self.scatter(b)

    def Diffusion(self, u, mask = True, DSS = True):
        l = 0; Au = self.zeros(self.nv)
        for el in self.mesh.elementsL:
            dx  = self.Ajac( self.mesh.size(l) ); l += 1
            Au[el] = (self.dfac/dx) * self.A.dot(u[el])
        if DSS : return self.dssum(Au, mask)
        return Au

    def Diffusion1(self, u, a, mask = True, DSS = True): #function overloaded
#         u = a * u  # given diffusion coefficient a
        l = 0; Au = self.zeros(self.nv)
        for el in self.mesh.elementsL:
            dx  = self.Ajac( self.mesh.size(l) ); l += 1
            Au[el] = (self.dfac/dx) * self.np.dot(self.D.T,(self.gllweights * a[el]) * (self.D.dot(u[el])))
        if DSS : return self.dssum(Au, mask)
        return Au

    def Convection1(self, u, a, mask = True, DSS = True):
        # a = linear convecting velocity
        Cu = self.zeros(self.nv)
        for el in self.mesh.elementsL:
            Cu[el] = (self.gllweights * a[el]) * (self.D.dot(u[el]))
        if DSS : return self.dssum(Cu, mask)
        return Cu

    def Convection(self, u, mask = True, DSS = True):
        # a = linear convecting velocity
        Cu = self.zeros(self.nv)
        for el in self.mesh.elementsL:
            Cu[el] = (self.gllweights * u[el]) * (self.D.dot(u[el]))
        if DSS : return self.dssum(Cu, mask)
        return Cu

    def assembleMatrix(self, assemblyMat): # sparse matrix assembly
        # Assembles the stiffness matrix
        if (assemblyMat or (self.Amat is None)):

            print("Assembling stiffness matrix...")
            l = 0; I = []; J = []; V = []
            lnode = None
            if self.mesh.anyPeriodic: lnode = self.mesh.Pnodes[0]
            for el in self.mesh.elements:
                N = len(el)
                if self.mesh.anyPeriodic: el[el==self.mesh.Pnodes[1]] = lnode
                v = self.np.ones((N,))
                I += self.np.kron(el,v).tolist()
                J += self.np.kron(v,el).tolist()
                dx  = self.Ajac( self.mesh.size(l) ); l += 1
                A = (self.dfac/dx) * self.A
                if (self.rfac): A += (self.mfac*dx) * self.B
                V += A.reshape((N*N,)).tolist()

            Amat = self.coo_matrix((V,(I,J)),shape=(self.nvglobal,self.nvglobal)).tocsr()
            self.Amat = Amat[self.np.ix_(self.uvec,self.uvec)].tocsc()
                
            self.lu = self.linalg.splu(self.Amat)

            if (self.mesh.anyDirbc):
                u = self.np.zeros((self.nvglobal,))
                self.insertbc(u)
                self.bcontr = Amat.dot(u)

        if (self.mesh.anyDirbc): 
            return self.bcontr
        else: 
            return self.np.zeros((self.nvglobal,))


#     def assembleMassMatrix(self, assemblyMat): # sparse matrix assembly
#         # Assembles the mass matrix
#         if (assemblyMat or (self.Mmat is None)):
# 
#             l = 0; I = []; J = []; V = []
#             for el in self.mesh.elements:
#                 N = len(el)
#                 v = self.np.ones((N,))
#                 I += self.np.kron(el,v).tolist()
#                 J += self.np.kron(v,el).tolist()
#                 dx  = self.mesh.size(l) / 2; l += 1
#                 B = (dx) * self.B
#                 V += B.reshape((N*N,)).tolist()
#     
#             Mmat = self.coo_matrix((V,(I,J)),shape=(self.nvglobal,self.nvglobal)).tocsr()
#             self.Mmat = Mmat[self.np.ix_(self.uvec,self.uvec)]


    def assembleVector(self, f, assemblyMat):
        b = -self.assembleMatrix(assemblyMat); l = 0
        if self.MassLhs:
            for el in self.mesh.elements:
                dx  = self.mesh.size(l) / 2; l += 1
                b[el] += (dx * self.mfac * self.gllweights) * f[el]
        else: b += f

#         b = self.periodic(b)
        self.bvec = b[self.uvec]
        return None

    def Ajac(self, dx):
        return (dx/2)

    def insert(self, ut, ug):
        ug[self.uvec] = ut
        return None

    def insertbc(self,ug): # insert appropriate boundary conditions
        if (self.mesh.anyDirbc): 
            if (len(ug)==self.nvglobal):
                ug[self.maskvec] = self.bc
            elif (len(ug)==self.nv):
                ug[self.maskvecL] = self.bc
        if (self.mesh.anyPeriodic): 
            ug[self.mesh.Pnodes[1]] = ug[self.mesh.Pnodes[0]]
        return None

    def innerprod(self,u,*args):
        if (len(args)>0):
            if (len(u)==self.nv):
                ug = self.gather(u,False,False)
                vg = self.gather(args[0],False,False)
                return  ug.dot(vg)
            else: return  u.dot(args[0])
        else:
            if (len(u)==self.nv):
                ug = self.gather(u,False,False)
                return  ug.dot(ug)
            else: return  u.dot(u)


class AssemblyKdV(Assembly):

    tol = 1e-8
    epsilon = 2.220446049250313e-16 
        
    from scipy.sparse import coo_matrix
    from scipy.sparse import linalg as linalg
    
    def __init__(self, mesh, gll, Re):

        super().__init__(mesh, gll, Re)

        self.dfac = -1.

#         D = gll.Lagrange2ndDerivativeMatrix()
        D = gll.PolynomialDerivativeMatrix(k=2)
        self.A = Re * self.np.dot(self.D.T, (D.T * gll.gllweights).T)
#         self.A = Re * (self.np.dot(D.T, (self.D.T * gll.gllweights).T))

    def Ajac(self, dx):
        return ( (dx/2)**2 )



class AssemblyDG(Assembly):
    """ Methods required to solve hyperbolic conservation laws """

        
    def __init__(self, mesh, gll, Re):

        super().__init__(mesh, gll, Re)

        self.Dt = self.D.T
        self.Ne = mesh.getNe()
        self.Ned = mesh.getnEdges()
        self.initd = None
        self.flowtype = None
        self.limit = False
        self.order = 1
        self.ltype = None
        self.m0 = 0.
        self.M0 = 1.
        self.bc = None
        self.tol = 1e-8
        self.epsilon = 2.220446049250313e-16 
        self.a = None
        self.valida = False
        self.src = None
        self.ldeno = gll.ldeno
        self.gllnodes = gll.gllnodes
        self.upstream = None
        self.Ampp = None
#         self.LBM = gll.LagrangeBasisMatrix(self.gllweights)

    def MassInv(self, u, mask = True):
        l = 0
        for el in self.mesh.elements:
            dx = self.mesh.size(l) / 2; l += 1
            u[el] /= ((dx*self.mfac) * self.gllweights)
        if (mask): self.mask(u)
        return None
        
    def Derivative(self, u, t, Flux=False):
        """ Applies the derivative matrix """
        
        if self.mesh.dG0: return self.np.zeros((self.nv,))

        f = self.initd.compute4flux(u, False, t, self.a) if Flux\
                else u

        Du = self.np.zeros((self.nv,))
        for el in self.mesh.elements:
            Du[el] = self.np.dot(self.gllweights * f[el], self.D)
#             Du[el] = self.Dt.dot(self.gllweights * f[el])

        return Du

    def cellAverages(self, u):

        el = 0
        ubar = self.np.zeros((self.Ne,))
        for K in self.mesh.elements:
            ubar[el] = (self.gllweights.dot(u[K])) / 2
            el += 1

        return ubar

    def upStream(self, uc):
        self.upstream = self.np.zeros((self.Ned,),dtype=int)
        
        for el in range(self.Ne):
            r = self.mesh.edges2[el][1]
            neigh = self.mesh.neighbours2[r][1]
            node = self.mesh.edges[r]
            if uc[node] < 0 and neigh >= 0:
                self.upstream[r] = self.mesh.edges[self.mesh.edges2[neigh][0]]
            else: self.upstream[r] = node
            l = self.mesh.edges2[el][0]
            neigh = self.mesh.neighbours2[l][0]
            node = self.mesh.edges[l]
            if uc[node] > 0 and neigh >= 0:
                self.upstream[l] = self.upstream[self.mesh.edges2[neigh][1]]
            else: self.upstream[l] = node

        if self.mesh.anyPeriodic:
            self.upstream[self.mesh.Bedges[0][0]] = self.upstream[self.mesh.Bedges[1][0]]

        return None

    def LaxFriedrichs(self, u, t):
        """ Lax-Friedrichs approx. Riemann solver """

        F = self.np.zeros((self.Ned,))
        if self.mesh.anyDirbc: self.insertbc(u)
        f = self.initd.compute4flux(u, True, t, self.a)
        df = self.initd.compute4dflux(u, True, t, self.a)
        l = self.mesh.Bedges[0][0]
        r = self.mesh.Bedges[1][0]
        bfl = f[l]; bfr = f[r]
        bdfl = df[l]; bdfr = df[r]
        bul = u[self.mesh.edges[l]]
        bur = u[self.mesh.edges[r]]

        for el in range(self.Ne):
            r = self.mesh.edges2[el][1]
            neigh = self.mesh.neighbours2[r][1]
            offside = (neigh < 0)
            ur = bur; fr = bfr; dfr = bdfr
            if (neigh >= 0):
                ur = u[self.mesh.edges[self.mesh.edges2[neigh][0]]]
                fr = f[self.mesh.edges2[neigh][0]]
                dfr = df[self.mesh.edges2[neigh][0]]
            F[r] = ( f[r] + fr + max(abs(dfr),abs(df[r])) * (u[self.mesh.edges[r]] - ur) ) / 2
            l = self.mesh.edges2[el][0]
            neigh = self.mesh.neighbours2[l][0]
            if (neigh < 0):
                F[l] = ( bfl + f[l] + max(abs(bdfl),abs(df[r])) * (bul - u[self.mesh.edges[l]]) ) / 2
            else: F[l] = F[self.mesh.edges2[neigh][1]]

        if self.mesh.anyPeriodic:
            F[self.mesh.Bedges[0][0]] = F[self.mesh.Bedges[1][0]]

        return F

    def EnquistOsher(self, u, t):
        """ Enquist-Osher approx. Riemann solver """

        F = self.np.zeros((self.Ned,)); f0 = None
        if self.mesh.anyDirbc: self.insertbc(u)
        f = self.initd.compute4flux(u, True, t, self.a)

        if self.flowtype['nonlinear']:
            u0 = self.np.zeros(self.np.shape(u))
            f0 = self.initd.compute4flux(u0, True, t)
        else: 
            f0 = self.np.zeros(self.np.shape(f))
        F = f0

        dvec = self.src.optima()
        periodic = self.mesh.anyPeriodic
        bdr = self.mesh.Bedges[:,0]

        for el in range(self.Ne):
            r = self.mesh.edges2[el][1]
            neigh = self.mesh.neighbours2[r][1]
            ur = None; fr = None; f0k = None
            offside = (neigh < 0)
            if offside:
                ur = u[self.mesh.edges[bdr[0]]] if periodic \
                        else u[self.mesh.edges[bdr[1]]]
                fr = f[bdr[0]] if periodic else f[bdr[1]]
            else:
                nedge = self.mesh.edges2[neigh][0]
                ur = u[self.mesh.edges[nedge]]; fr = f[nedge]
            x = self.mesh.nodes[self.mesh.edges[r]]
            a = self.a[self.mesh.edges[r]] if self.valida else 0.
            if self.valida: dvec = self.src.optima(a)
            F[r] += self.EO(dvec, u[self.mesh.edges[r]], f[r], f0[r], x, t, a)
            if offside: f0k = f0[bdr[0]] if periodic else  f0[bdr[1]]
            else: f0k = f0[self.mesh.edges2[neigh][0]]
            F[r] += self.EO(dvec, ur, fr, f0k, x, t, a, False)

            l = self.mesh.edges2[el][0]
            neigh = self.mesh.neighbours2[l][0]
            x = self.mesh.nodes[self.mesh.edges[l]]
            if (neigh < 0):
                a = self.a[self.mesh.edges[l]] if self.valida else 0.
                if self.valida: dvec = self.src.optima(a)
                F[l] += self.EO(dvec, u[self.mesh.edges[bdr[0]]], f[bdr[0]], f0[bdr[0]], x, t, a)
                F[l] += self.EO(dvec, u[self.mesh.edges[l]], f[l], f0[l], x, t, a, False)
            else: 
                F[l] = F[self.mesh.edges2[neigh][1]]

        if periodic: F[bdr[0]] =  F[bdr[1]]

        return F

    def Upwind(self, u, t, f=None):
        """ Upwind approx. Riemann solver.
        To be used only in the case of linear xtic velocity! """

        F = self.np.zeros((self.Ned,))
        if self.mesh.anyDirbc: self.insertbc(u)
        uc = self.initd.compute4a(t)
        if (f is None): f = self.initd.compute4flux(u, True, t, self.a)
        F = f
        
        for el in range(self.Ne):
            r = self.mesh.edges2[el][1]
            neigh = self.mesh.neighbours2[r][1]
            node = self.mesh.edges[r]
            if uc[node] < 0 and neigh >= 0:
                F[r] = f[self.mesh.edges2[neigh][0]]
            l = self.mesh.edges2[el][0]
            neigh = self.mesh.neighbours2[l][0]
            node = self.mesh.edges[l]
            if uc[node] > 0 and neigh >= 0:
                F[l] = F[self.mesh.edges2[neigh][1]]

        if self.mesh.anyPeriodic:
            F[self.mesh.Bedges[0][0]] = F[self.mesh.Bedges[1][0]]

        return F

    def Godunov(self, u, t, fopt=0):
        """ Godunov approx. Riemann solver:
            We assume that the flux function is continuous! """

        F = self.np.zeros((self.Ned,))
        if self.mesh.anyDirbc: self.insertbc(u)
        if fopt==2:
            f = self.initd.compute4dflux(u, True, t) # f'(u)
        elif fopt==1:
            f = self.initd.compute4A(self.mesh.nodes[self.mesh.edges], t, u[self.mesh.edges])  # f(u)/u
        else:
            f = self.initd.compute4flux(u, True, t, self.a)

        l = self.mesh.Bedges[0][0]
        r = self.mesh.Bedges[1][0]
        bfl = f[l]; bfr = f[r]
        bul = u[self.mesh.edges[l]]
        bur = u[self.mesh.edges[r]]

        dvec = self.src.optima()

        for el in range(self.Ne):
            l = self.mesh.edges2[el][0]
            r = self.mesh.edges2[el][1]
            neigh = self.mesh.neighbours2[r][1]
            offside = (neigh < 0)
            ur = bur if (offside) else u[self.mesh.edges[self.mesh.edges2[neigh][0]]]
            fr = bfr if (offside) else f[self.mesh.edges2[neigh][0]]
            k = self.mesh.edges[r]
            if self.valida: dvec = self.src.optima(self.a[k])
            if (u[k] < ur):
                fr = min(f[r], fr)
                for opt in dvec:
                    if self.inside(opt[0], u[k], ur): fr = min(opt[1], fr)
            else:
                fr = max(f[r], fr)
                for opt in dvec:
                    if self.inside(opt[0], ur, u[k]): fr = max(opt[1], fr)
            neigh = self.mesh.neighbours2[l][0]
            offside = (neigh < 0)
            fl = bfl if (offside) else F[self.mesh.edges2[neigh][1]]
            if (offside):
                k = self.mesh.edges[l]
                if (bul < u[k]):
                    fl = min(fl, f[l])
                    for opt in dvec:
                        if self.inside(opt[0], bul, u[k]): fl = min(opt[1], fl)
                else:
                    fl = max(fl, f[l])
                    for opt in dvec:
                        if self.inside(opt[0], u[k], bul): fl = max(opt[1], fl)

            F[l] = fl; F[r] = fr

        if self.mesh.anyPeriodic:
            F[self.mesh.Bedges[0][0]] = F[self.mesh.Bedges[1][0]]

        return F


    def numericalFlux(self, f, u, dt, t, MassInv=True):
        """ Approximate numerical fluxes using appropriate Riemann solver """

        F = self.np.zeros((self.Ned,))
        if (self.order < 2):
            if self.flowtype['nonlinear']:
                F = self.Godunov(u, t)
            else: 
                F = self.Upwind(u, t)
        else:
            if self.flowtype['nonlinear']:
                F = self.EnquistOsher(u, t)
            else:
                F = self.LaxFriedrichs(u, t)

        # Add numerical flux and factor the mesh sizes """
        l = self.mesh.edges2[:,0]
        r = self.mesh.edges2[:,1]
        f[self.mesh.edges[l]] += F[l]
        f[self.mesh.edges[r]] -= F[r]
            
        if MassInv:
            self.MassInv(f, self.mesh.anyDirbc)
#             el = 0
#             for K in self.mesh.elements:
#                 dx = self.mesh.size(el) / 2; el += 1
#                 f[K] /= (self.gllweights * dx)
    
        return f

    def limiter(self, u):
        """ numerical flux limiter """
        if self.mesh.anyDirbc: self.insertbc(u)

        if (self.mesh.dG0): # re-average
            Q = self.cellAverages(u)
            u[self.mesh.edges[self.mesh.edges2[:,0]]] = Q
            u[self.mesh.edges[self.mesh.edges2[:,1]]] = Q
        elif self.flowtype['nonlinear']:
            u = self.minmod(u)
        else:
            u = self.minmod(u)

        return u

    def mpplim(self, u):
        """ mpp and pp flux limiter """
        if self.ltype is None: return u
        if self.mesh.anyDirbc: self.insertbc(u)

        l = 0
        if (self.mesh.dG0):
            Q = self.cellAverages(u)
            u[self.mesh.edges[self.mesh.edges2[:,0]]] = Q
            u[self.mesh.edges[self.mesh.edges2[:,1]]] = Q
        elif self.ltype['MPP']:
            Q = self.cellAverages(u)
            Mp = None; mp = None; deg = self.mesh.getN() - 1
#             deg2o3 = deg==2 or deg==3
#             if deg2o3: mp, Mp = self.MinMaxPoly2(u)
            for el in self.mesh.elements:
                m = min(u[el]); M = max(u[el])
#                 m = mp[l] if deg2o3 else min(u[el])
#                 M = Mp[l] if deg2o3 else max(u[el])
                if (m >= self.m0) and (M <= self.M0): l += 1; continue
                theta = 1.
                if abs(M - Q[l]) > self.epsilon:
                    theta = min(theta, abs( (self.M0 - Q[l]) / (M - Q[l]) ))
                if abs(m - Q[l]) > self.epsilon:
                    theta = min(theta, abs( (self.m0 - Q[l]) / (m - Q[l]) ))
                u[el] = theta * (u[el] - Q[l]) + Q[l]
#                 if min(u[el])<self.m0: print(min(u[el])); print(Q[l])
                l += 1
        elif self.ltype['PP']:
            Q = self.cellAverages(u)
            Mp = None; mp = None; deg = self.mesh.getN() - 1
            deg2o3 = deg==2 or deg==3
            if deg2o3: mp, Mp = self.MinMaxPoly(u)
            for el in self.mesh.elements:
#                 m = min(u[el])
                m = mp[l] if deg2o3 else min(u[el])
                theta = 1.
                if abs(m - Q[l]) > self.epsilon:
                    theta = min(theta, abs( Q[l] / (m - Q[l]) ) )
                u[el] = theta * (u[el] - Q[l]) + Q[l]
                l += 1
        return u

    def inside(self, x, a, b):
        return ( (x-a) * (b-x) > self.epsilon ) if (a < b) \
        else ( (x-b) * (a-x) > self.epsilon )

    def minmod(self, u):
        Q = self.cellAverages(u)
        Qbl = Q[self.mesh.Bedges[0][2]]
        Qbr = Q[self.mesh.Bedges[1][2]]
        el = 0; N = self.mesh.getN()
        for K in self.mesh.elements:
            Ql = Qbl if (self.mesh.neighbours[el][0] < 0) else Q[self.mesh.neighbours[el][0]]
            Qr = Qbr if (self.mesh.neighbours[el][1] < 0) else Q[self.mesh.neighbours[el][1]]
            l = self.mesh.edges2[el][0]; r = self.mesh.edges2[el][1]
            ul = Q[el] - u[self.mesh.edges[l]] 
            ur = u[self.mesh.edges[r]] -  Q[el]
            Dl = self.minmodf(Q[el]-Ql, Qr-Q[el])
            ul = self.minmodf(ul, Dl); ur = self.minmodf(ur, Dl)
#             limit = True
            limit = ( abs(ul - Q[el] + u[self.mesh.edges[l]]) > self.tol\
                    or abs(ur - u[self.mesh.edges[r]] + Q[el]) > self.tol )
            if (limit):
#                 print('minmod limiter')
                xl = None; xr = None; x = None; rhs = None
                if (N > 2):
                    xl = self.mesh.nodes[self.mesh.edges[l]]
                    xr = self.mesh.nodes[self.mesh.edges[r]]
                    rhs = self.np.array([Q[el], Q[el] - ul, Q[el] + ur])
                    x = self.mesh.nodes[K]
#                 if (N > 3):
#                     A = self.np.array([\
#                             [(xl**3 + xl*xl*xr + xl*xr*xr + xr**3)/4, (xl*xl + xl*xr + xr*xr)/3, (xl + xr)/2],\
#                             [xl**3, xl*xl, xl],\
#                             [xr**3, xr*xr, xr] ])
#                     coeffs = self.np.linalg.solve(A, rhs)
#                     u[K] = coeffs[0] * x*x*x + coeffs[1]*x*x + coeffs[2]*x
#                 elif (N > 2):
                    A = self.np.array([\
                            [(xl*xl + xl*xr + xr*xr)/3, (xl + xr)/2, 1],\
                            [xl*xl, xl, 1],\
                            [xr*xr, xr, 1] ])
                    coeffs = self.np.linalg.solve(A, rhs)
                    u[K] = coeffs[0] * x*x + coeffs[1]*x + coeffs[2]
                else:
                    u[self.mesh.edges[l]] = Q[el] - ul
                    u[self.mesh.edges[r]] = Q[el] + ur
            el += 1
        return u

    def minmodf(self, a, b):
        return ( a if (abs(a) < abs(b)) else b ) if (a*b > 0) else 0.

    def superbee(self, a, b):
        return max(0., max(min(2*a, b), min(a, 2*b)))

    def vanLeer(self, a, b):
        signb = 1. if (b > 0) else -1.
        return ( (2*a*b)/(a*signb + b) if (a > 0) else 0.)

    def EO(self, dvec, u1, f1, f0, x, t, a, positive=True):
        u0 = 0. 
        F = 0. if positive else f0
        dfluxsign = True
        for u in dvec:
            if self.inside(u[0], u0, u1):
                dfluxsign = self.src.dfluxsign(u0, u[0], x, t, a)
                if (not(dfluxsign or positive) or (dfluxsign and positive)):
                    F += (u[1] - f0)
                u0 = u[0]; f0 = u[1]
        dfluxsign = self.src.dfluxsign(u0, u1, x, t, a)
        if (not(dfluxsign or positive) or (dfluxsign and positive)):
            F += (f1 - f0)

        return F

    def computeJumps(self, u):

        jump = self.np.zeros(self.np.shape(u))
        l = self.mesh.edges[self.mesh.edges2]
        ln = l[self.mesh.neighbours[:,1],0]
        rn = l[self.mesh.neighbours[:,0],1]
        jump[l[:,0]] = u[rn] - u[l[:,0]]
        jump[l[:,1]] = u[ln] - u[l[:,1]]

        return ( jump )

    def TV(self, u):

        return sum(abs(u[1:] - u[:-1]))

    def slopeRatio(self, u):

        v = u[self.mesh.edges[self.mesh.edges2]]

        return (v[:,1]-v[:,0]) / self.mesh.element_sizes

    def MinMaxPoly(self, u):
        s = self.gllnodes
        deg = len(s)-1; l = 0
        k = self.ldeno
        m = self.np.zeros((self.mesh.getNe(),))
        M = self.np.zeros((self.mesh.getNe(),))
        if deg==2:
            def f(v,i, j=2, k=3):
                return v[self.np.mod(self.np.arange(i,i+j), k)].prod()
            def g(v,i):
                return v[self.np.mod(self.np.arange(i,i+2), 3)].sum()
            def P2(t, u):
                v = t - s
                return self.np.array([f(v,1), f(v,2), f(v,0)]).dot( u )
            A = self.np.array([[2., 2., 2.],\
                    [-g(s,1), -g(s,2), -g(s,0)]] )
            for el in self.mesh.elements:
                v = u[el]
                m[l] = min(v[0], v[-1])
                M[l] = max(v[0], v[-1])
                v = v / k
                coeffs = A.dot( v )
                if (abs(coeffs[0]) < self.epsilon): 
                    m[l] = min(m[l], P2(coeffs[1]/coeffs[0], v))
                    M[l] = max(M[l], P2(coeffs[1]/coeffs[0], v))
        elif deg==3:
            def f(v,i, j=3, k=4):
                return v[self.np.mod(self.np.arange(i,i+j), k)].prod()
            def g(v,i):
                return v[self.np.mod(self.np.arange(i,i+3), 4)].sum()
            def h(v,i):
                m = v[self.np.mod(self.np.arange(i,i+3), 4)]
                return self.np.array([f(m,1,2,3), f(m,2,2,3), f(m,0,2,3)]).sum()
            def P3(t, u):
                x = t - s
                return self.np.array([f(x,1), f(x,2), f(x,3), f(x,0)]).dot( u )
            if self.Ampp is None:
                self.Ampp = self.np.array([[3., 3., 3., 3.],\
                        [-2*g(s,1), -2*g(s,2), -2*g(s,3), -2*g(s,0) ],\
                        [h(s,1), h(s,2), h(s,3), h(s,0) ]])
            for el in self.mesh.elements:
                v = u[el]
                m[l] = min(v[0], v[-1])
                M[l] = max(v[0], v[-1])
                v = v / k
                coeffs = self.Ampp.dot( v )
                D = coeffs[1]**2 - 4 * coeffs[0] * coeffs[2]
                roots = self.np.roots( coeffs )
                if abs(coeffs[0]) > self.epsilon:
                    if (abs(D) < self.epsilon): # repeated roots
                        root = - coeffs[1] / (2 * coeffs[0])
                        m[l] = min(m[l], P3(roots[0], v))
                        M[l] = max(M[l], P3(roots[0], v))
                    elif ( D > self.epsilon ):
                        m[l] = min(m[l], P3(roots[0], v), P3(roots[1], v))
                        M[l] = max(M[l], P3(roots[0], v), P3(roots[1], v))
                elif abs(coeffs[1]) > self.epsilon: # degenerate case
                    m[l] = min(m[l], P3(coeffs[2]/coeffs[1], v))
                    M[l] = max(M[l], P3(coeffs[2]/coeffs[1], v))
        return m, M


    def MinMaxPoly2(self, u):
        s = self.gllnodes
        deg = len(s)-1; l = 0
        k = self.ldeno
        m = self.np.zeros((self.mesh.getNe(),))
        M = self.np.zeros((self.mesh.getNe(),))
        if deg==2:
            def f(v,i, j=2, k=3):
                return v[self.np.mod(self.np.arange(i,i+j), k)].prod()
            def g(v,i):
                return v[self.np.mod(self.np.arange(i,i+2), 3)].sum()
            def P2(t, u):
                v = t - s
                return self.np.array([f(v,1), f(v,2), f(v,0)]).dot( u )
            A = self.np.array([[2., 2., 2.],\
                    [-g(s,1), -g(s,2), -g(s,0)]] )
            for el in self.mesh.elements:
                v = u[el]
                m[l] = min(v[0], v[-1])
                M[l] = max(v[0], v[-1])
                v = v / k
                coeffs = A.dot( v )
                if (abs(coeffs[0]) < self.epsilon): 
                    m[l] = min(m[l], P2(coeffs[1]/coeffs[0], v))
                    M[l] = max(M[l], P2(coeffs[1]/coeffs[0], v))
        elif deg==3:
            def slice(i):
                return self.np.mod(self.np.arange(i,i+deg), deg+1)
            def f(v,i):
                return self.np.prod(v[slice(i)])
            def g(v,i):
                return -2 * self.np.sum(v[slice(i)])
            def h(v,i):
                return self.np.sum(v[slice(i)] * v[slice(i+1)])
            def P3(t, u):
                x = t - s
                return self.np.array([f(x,1), f(x,2), f(x,3), f(x,0)]).dot( u )
            if self.Ampp is None:
                self.Ampp = self.np.array([[3., 3., 3., 3.],\
                        [g(s,1), g(s,2), g(s,3), g(s,0) ],\
                        [h(s,1), h(s,2), h(s,3), h(s,0) ]])
            k = self.np.array([f(s,1), f(s,2), f(s,3), f(s,0)]) # test test
            for el in self.mesh.elements:
                v = u[el]
#                 m[l] = min(v[0], v[-1])
#                 M[l] = max(v[0], v[-1])
                m[l] = min(v); M[l] = max(v)
#                 if (M[l] <= self.M0 and m[l] <= self.m0): continue
                v = v / k
                coeffs = self.Ampp.dot( v )
                D = coeffs[1]**2 - 4 * coeffs[0] * coeffs[2]
                if D > 0:
                    roots = self.np.roots(coeffs)
                    root = roots[0]
                    if self.np.isreal(root) and abs(root)<1:
                        P3val = P3(root, v)
                        m[l] = min(m[l], P3val)
                        M[l] = max(M[l], P3val)
                    root = roots[1]
                    if self.np.isreal(root) and abs(root)<1:
                        P3val = P3(root, v)
                        m[l] = min(m[l], P3val)
                        M[l] = max(M[l], P3val)
#                 if abs(coeffs[0]) > self.epsilon:
#                     if (abs(D) < self.epsilon): # repeated roots
#                         root = - coeffs[1] / (2 * coeffs[0])
#                         if abs(root) < 1:
#                             P3val = P3(root, v)
#                             m[l] = min(m[l], P3val, self.m0)
#                             M[l] = max(M[l], P3val, self.M0)
#                     elif ( D > self.epsilon ):
#                         root = (-coeffs[1] - self.np.sqrt(D)) / (2 * coeffs[0]) 
#                         if abs(root) < 1:
#                             P3val = P3(root, v)
#                             m[l] = min(m[l], P3val, self.m0)
#                             M[l] = max(M[l], P3val, self.M0)
#                         root = (-coeffs[1] + self.np.sqrt(D)) / (2 * coeffs[0]) 
#                         if abs(root) < 1:
#                             P3val = P3(root, v)
#                             m[l] = min(m[l], P3val, self.m0)
#                             M[l] = max(M[l], P3val, self.M0)
#                 elif abs(coeffs[1]) > self.epsilon: # degenerate case
#                     root = - coeffs[2]/coeffs[1]
#                     if abs(root) < 1:
#                         P3val = P3(root, v)
#                         m[l] = min(m[l], P3val, self.m0)
#                         M[l] = max(M[l], P3val, self.M0)
        return m, M


