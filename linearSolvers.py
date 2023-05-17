class LinearSolvers:

#     from scipy.sparse import linalg as linalg
    
    def __init__(self, semAss, tol = 1e-11, maxit = 1000):
        self.np = semAss.np
        self.semAss = semAss
        self.tol = min(tol,1e-11)
        self.rtol = 1e-3 * self.tol
        self.maxit = max(max(maxit,1000), self.semAss.nDof)
        self.u0 = self.np.array([]) # initial guess for iterative solvers
        self.residualError = None
        self.iterationNo = None
        self.validSol = True
        self.solverType = {"cholesky":False, "cg":False, "bicgstab":False, "gmres":True}
        self.directSolver = True
        self.initGuess = False; self.initGuessflag = False
        self.printSolverInfo("Default")
        self.matvec = False # flag to perform iterative solutions based on user-defined matrix-vector operations
        self.assembleMatrix = True
        self.displayInfo = True

    def printSolverInfo(self, dfault):
        llen = 60
        print('='*llen)
        print("%s parameters for linear solver:" % dfault)
        print('-'*llen)
        print("maximum allowed # of iterations = %i" % self.maxit)
        print("residual error tolerance (outer) = %2.2e" % self.tol)
        if (self.directSolver):
            print("Direct linear solver")
            if self.solverType["cholesky"]:
                print("Cholesky decomposition")
        else:
            print("Iterative linear solver:")
            if self.solverType["CG"]:
                print(" Conjugate Gradient")
            elif self.solverType["bigstab"]:
                print(" bi-conjugate gradient with stabilization")
            else: 
                print("other")
        print('-'*llen)
        return None

    def invalidVector(self, u):
        return  any(self.np.isnan(u)) or any(self.np.isinf(u))

    def setSolverType(self, solvertype, userLinearOperator=False):
        self.solverType = self.solverType.fromkeys(self.solverType,False)
        self.solverType[solvertype.lower()] = True
        self.directSolver = False
        self.matvec = userLinearOperator
        return None

    def useDirectSolver(self):
        self.directSolver = True
        self.matvec = False
        return None

    def setInitialGuess(self,u0):
        if self.directSolver: self.u0 = u0
        return None

    def solve(self, f):

        if self.displayInfo:
            print("Linear system ...")

        u = self.np.zeros(self.np.shape(f))
        ut = None
        flag = 0

        if self.directSolver:
            self.semAss.assembleVector(f, self.assembleMatrix)
#             ut = self.semAss.linalg.spsolve(self.semAss.Amat, self.semAss.bvec)
            ut = self.semAss.lu.solve(self.semAss.bvec)
            if self.displayInfo:
                print("... solved using a direct (matrix decomposition) method: Scipy's inbuilt superlu")

        elif self.solverType['cg']:
            if (self.matvec):
                if (self.semAss.nDnodes>0): self.BoundaryContributions(f)
                LinOp = self.semAss.linalg.LinearOperator( (self.semAss.nv,self.semAss.nv), matvec=self.A )
                u,flag = self.semAss.linalg.cg(LinOp, f, tol=self.rtol,maxiter=self.maxit)
            else: 
                self.semAss.assembleVector(f, self.assembleMatrix)
                ut,flag = self.semAss.linalg.cg(self.semAss.Amat, self.semAss.bvec,tol=self.rtol,maxiter=self.maxit)
                self.semAss.insert(ut, u)
                if self.displayInfo:
                    print("... solved using scipy's inbuilt iterative method , cg")
    
        elif self.solverType['bicgstab']:
            if (self.matvec):
                if (self.semAss.nDnodes>0): self.BoundaryContributions(f)
                LinOp = self.semAss.linalg.LinearOperator( (self.semAss.nv,self.semAss.nv), matvec=self.A )
                u,flag = self.semAss.linalg.bicgstab(LinOp, f, tol=self.rtol,maxiter=self.maxit)
            else: 
                self.semAss.assembleVector(f, self.assembleMatrix)
                ut,flag = self.semAss.linalg.bicgstab(self.semAss.Amat, self.semAss.bvec,tol=self.rtol,maxiter=self.maxit)
                if self.displayInfo:
                    print("... solved using scipy's inbuilt iterative method , bicgstab")
    
        elif self.solverType['gmres']:
            if (self.matvec):
                if (self.semAss.nDnodes>0): self.BoundaryContributions(f)
                LinOp = self.semAss.linalg.LinearOperator( (self.semAss.nvglobal,self.semAss.nvglobal), matvec=self.A )
                u,flag = self.semAss.linalg.gmres(LinOp, f, tol=self.rtol,maxiter=self.maxit)
            else: 
                self.semAss.assembleVector(f, self.assembleMatrix)
                ut,flag = self.semAss.linalg.gmres(self.semAss.Amat, self.semAss.bvec,tol=self.rtol,maxiter=self.maxit)
                if self.displayInfo:
                    print("... solved using scipy's inbuilt iterative method , gmres")
    
        else: print("Undefined solver type!")

        if (flag>0): print("@LinearSolver: Iterative tolerance not achieved!")
        elif (flag<0): 
            print("@LinearSolver: ilegal input or breakdown of iterative method")
            self.validSol = False

        if (ut is not None): self.semAss.insert(ut,u)

        self.semAss.insertbc(u)

        return u

    def A(self,u): # define linear operator
        if (self.matvec): return self.semAss.Diffusion(u)
        else: return self.semAss.Amat.dot(u)

    def BoundaryContributions(self, f):
        ug = self.np.zeros((self.semAss.nvglobal,))
        self.semAss.insertbc(ug)
        f -= self.A(self.semAss.scatter(ug))
        return None



class LinearSolverIMEX (LinearSolvers): # Derived class

    def __init__(self, semAss, tol = 1e-11, maxit = 1000):
        super().__init__(semAss, tol, maxit)

    def A(self,u): # define linear operator
        if (self.matvec):
            return self.semAss.Diffusion(u)
        else: return self.semAss.Amat.dot(u)

    def BoundaryContributions(self, f):
        ug = self.np.zeros((self.semAss.nvglobal,))
        self.semAss.insertbc(ug)
        f -= self.A(self.semAss.scatter(ug))
        return None
