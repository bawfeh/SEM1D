class AssemblyKdV(Assembly):

    tol = 1e-8
    epsilon = 2.220446049250313e-16 
        
    def __init__(self, mesh, gll, Re):

        super().__init__(mesh, gll, Re)

        from scipy.sparse import csc_matrix, eye
        from scipy.linalg import toeplitz #, solve_toeplitz
        self.csc_matrix  = csc_matrix
        self.toeplitz  = toeplitz
#         self.solve_toeplitz  = solve_toeplitz
        self.eye  = eye
        self.dx = mesh.nodes[1:] - mesh.nodes[:-1] # for finite difference mesh sizes
        self.dx = self.np.append(mesh.nodes[-1]-mesh.nodes[-2], self.dx)
        self.MassLhs = False
        if (self.mesh.anyPeriodic):
            self.nDof += 1
            self.uvec[self.mesh.Pnodes[1]] = True

    def adjust(self, stencil):
        if self.mesh.anyPeriodic:
            return self.np.mod(stencil,self.nvglobal)
        else:
            lbool = (stencil < 0)
            rbool = (stencil >= self.nvglobal)
            stencil[lbool] = 0
            stencil[rbool] = self.nvglobal-1
        return stencil

    def Convection(self, u, mask = True, DSS = False): # Central differences
        # a = linear convecting velocity
        Cu = self.zeros(self.nvglobal)
        refstencil = self.np.array([-1,0,1]) # 3-point stencil
        w = self.gather(u,False,mask)
        w = w**2 / 2
        for j in range(self.nvglobal):
            stencil = self.adjust(refstencil+j)
            Cu[j] += ( (w[stencil[-1]] - w[stencil[1]]) / (2*self.dx[stencil[-1]])\
                    + (w[stencil[1]] - w[stencil[0]]) / (2*self.dx[stencil[1]]) )
#         if (mask and self.anyDirbc) : Cu[self.maskvec] = 0.
        return self.scatter(Cu)

    def Diffusion(self, u, mask = True, DSS = False): # Central differences
        Au = self.zeros(self.nvglobal)
        refstencil = self.np.array([-2,-1,0,1,2]) # 5-point stencil
        w = self.gather(u,False,mask)
        w = 1*u
        for j in range(self.nvglobal):
            stencil = self.adjust(refstencil+j)
            Au[j] += self.Re * ( (w[stencil[-1]] - w[stencil[3]]) / (self.dx[stencil[-1]] * self.dx[stencil[3]])\
                   - (w[stencil[3]] - w[stencil[2]]) / (self.dx[stencil[3]]**2)\
                   - (w[stencil[2]] - w[stencil[1]]) / (self.dx[stencil[2]] * self.dx[stencil[1]])\
                   + (w[stencil[1]] - w[stencil[0]]) / (self.dx[stencil[1]]**2) ) / (2*self.dx[stencil[2]])
        if (mask) : Au[self.maskvec] = 0.
        return self.scatter(Au)

    def assembleMatrix(self, assemblyMat): # sparse matrix assembly
        # Assembles the stiffness matrix
        if (assemblyMat or (self.Amat is None)):
            refstencil = self.np.array([-2,-1,0,1,2]) # 5-point stencil
#             v = self.np.zeros((5,))
# 
#             print("Assembling stiffness matrix...")
#             I = []; J = []; V = []
#             for j in range(self.nvglobal):
#                 stencil = self.adjust(refstencil+j)
#                 J += stencil.tolist()
#                 I += [j]*5
#                 v[0] =   -1. / (self.dx[stencil[1]] ** 2 )
#                 v[1] =  ( 1. / (self.dx[stencil[2]] * self.dx[stencil[1]]) +  1. / (self.dx[stencil[1]]**2) )
# #                 v[2] =  ( 1. / (self.dx[stencil[3]]**2) -  1. / (self.dx[stencil[2]] * self.dx[stencil[1]]) )
#                 v[3] =  - ( 1. / (self.dx[stencil[-1]] * self.dx[stencil[3]]) +  1. / (self.dx[stencil[3]]**2) )
#                 v[-1] =  1. / (self.dx[stencil[-1]] * self.dx[stencil[3]])
#                 v /= ( 2 * self.dx[stencil[2]] / self.Re )
#                 v[2] += self.mfac  # adding mass matrix
#                 V += v.tolist()
# 
#             Amat = self.coo_matrix((V,(I,J)),shape=(self.nvglobal,self.nvglobal)).tocsr()
# 
#             self.Amat = Amat[self.np.ix_(self.uvec,self.uvec)].tocsc()


            print("Assembling stiffness matrix...")
            c = self.np.zeros((self.nvglobal,))
            r = self.np.zeros((self.nvglobal,))
            stencil = self.adjust(refstencil)
#             r[-2] =   -1. / (self.dx[stencil[1]] ** 2 )
#             r[-1] =  ( 1. / (self.dx[stencil[2]] * self.dx[stencil[1]]) +  1. / (self.dx[stencil[1]]**2) )
#             r[0] =  ( 1. / (self.dx[stencil[3]]**2) -  1. / (self.dx[stencil[2]] * self.dx[stencil[1]]) )
#             r[1] =  - ( 1. / (self.dx[stencil[-1]] * self.dx[stencil[3]]) +  1. / (self.dx[stencil[3]]**2) )
#             r[2] =  1. / (self.dx[stencil[-1]] * self.dx[stencil[3]])
#             r /= ( 2 * self.dx[stencil[2]] / self.Re )
#             r[0] += self.mfac  # adding mass matrix
#             c[-2] =    1. / (self.dx[stencil[1]] ** 2 )
#             c[-1] =   - ( 1. / (self.dx[stencil[2]] * self.dx[stencil[1]]) +  1. / (self.dx[stencil[1]]**2) )
#             c[0] = - ( 1. / (self.dx[stencil[3]]**2) -  1. / (self.dx[stencil[2]] * self.dx[stencil[1]]) )
#             c[1] =  ( 1. / (self.dx[stencil[-1]] * self.dx[stencil[3]]) +  1. / (self.dx[stencil[3]]**2) )
#             c[2] =   -1. / (self.dx[stencil[-1]] * self.dx[stencil[3]])
#             c /= ( 2 * self.dx[stencil[2]] / self.Re )
#             c[0] += self.mfac  # adding mass matrix

            r[1] = -2
            r[2] = 1
            r[-1] = 2
            r[-2] = -1
            c[1] = 2
            c[2] = -1
            c[-1] = -2
            c[-2] = 1
            dx = self.dx[0]
            r /= (2*dx**3)
            c /= (2*dx**3)
            r *= self.Re
            c *= self.Re
            r[0] = self.mfac
            c[0] = self.mfac
            
#             self.r = r
#             self.c = c

            Amat = self.csc_matrix( self.toeplitz( c, r ) )

            self.Amat = Amat[self.np.ix_(self.uvec,self.uvec)]#.tocsc()
                
            self.lu = self.linalg.splu(self.Amat)

            if (self.mesh.anyDirbc):
                u = self.np.zeros((self.nvglobal,))
                self.insertbc(u)
                self.bcontr = Amat.dot(u)
                return self.bcontr

        return self.np.zeros((self.nvglobal,))

    def dssum(self, u, mask = True):
        # Direct stiffness summation: scatter[gather]
        return u

    def scatter(self, ug):
        # local to global SEM map // scatter (copy) operator
        return ug

    def gather(self, u, add=True, mask=True):
        # local to global SEM map: gather (add) operator
        return u

    def Mass(self, u, mask = True):
        u *= self.mfac
        if (mask): self.mask(u)
        return None

    def MassInv(self, u, assembleMatrix = False, mask = True):
#         if (mask): self.mask(u)
        return u

    def periodic(self, u, dfac = 1.): 
        # enforce periodic boundary conditions (if any), on lobal vector
        if (self.mesh.anyPeriodic):
            if (len(u)==self.nv): 
                self.Vgl = self.periodic(self.gather(u,False,False),dfac)
                return self.scatter(self.Vgl)
            else:
                l = self.mesh.Pnodes[0]; r = self.mesh.Pnodes[1]
                u[l] = u[r]
        return u


    def assembleConvection(self, assembleMat): # Central differences

        refstencil = self.np.array([-1,0,1]) # 3-point stencil

        if (assembleMat or (self.Cmat is None)):
            refstencil = self.np.array([-1,0,1]) # 3-point stencil
            v = self.np.zeros((3,))

            print("Assembling convection matrix...")
            I = []; J = []; V = []
            for j in range(self.nvglobal):
                stencil = self.adjust(refstencil+j)
                J += stencil.tolist()
                I += [j]*3
                v[0] =   -1 / (2*self.dx[stencil[1]]) 
                v[1] =  ( 1. / (2*self.dx[stencil[1]]) -  1. / (2*self.dx[stencil[-1]]) )
                v[-1] =  1. / (2*self.dx[stencil[-1]])
                V += v.tolist()

            self.Cmat = self.coo_matrix((V,(I,J)),shape=(self.nvglobal,self.nvglobal)).tocsr()

        return None



