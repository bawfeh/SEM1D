class timeStepping:


    def __init__(self, data, gll, mesh, bc, nsteps):

        from Assembly import Assembly
#         from Assembly import AssemblyKdV
        from linearSolvers import LinearSolvers
        import Initdata
        from Methods import DIRK_CF

        self.semAss = Assembly(mesh, gll, data.Re)
#         self.semAss = AssemblyKdV(mesh, gll, data.Re)
        self.semAss.rfac = True
        self.linearSolver = LinearSolvers(self.semAss, data.iterationTol, data.maxIters)
#         self.linearSolver = LinearSolvers(self.semAss,1e-6)
#         self.linearSolver.setSolverType('gmres')#, True)
        self.cf = DIRK_CF(gll.np,data.cfmethod)

        print('timestepping method: %s' % data.cfmethodName)

        self.bc = bc
        self.sources = Initdata.Sources(bc, mesh, data)
        self.mesh = mesh
        self.init = Initdata.Initdata(bc, mesh)
        self.gll = gll
        self.np = gll.np

        self.dt = data.dT / nsteps
        self.t0 = data.tb[0]
        self.tend = data.tb[-1]
        self.cdt = 1.
        self.t = 0.
        self.nsteps = nsteps
        self.nsteps2do = nsteps
        self.nstepsdone = 0
        self.Re = data.Re
        self.nv = self.semAss.nv
        self.epsilon = 2.220446049250313e-16 
        self.inf = 2.220446049250313e+16 
        self.stepBystepReport = True
        self.linearSolver.displayInfo = True

        self.ub = None; self.uc = None
        self.fs = None; self.Dub = None
        self.fs = self.semAss.scatter(self.sources.sources)
        self.sol = []
        self.u0 = self.init.compute4ex(self.t0)

        self.constantSteps = True
        self.varDirbc = False
        self.semAss.MassLhs = False
        self.exactsol = False
        self.storeAll = data.storeAll
        self.nonHomogeneousDirbc =  not(data.Dirbctype['homogeneous'] or data.Dirbctype['none']) and mesh.anyDirbc
        self.nonHomogeneousSource = not(data.sourcetype['homogeneous'])
        self.varSource = data.sourcetype['variable_s']
        self.varDirbc = data.Dirbctype['variable'] and mesh.anyDirbc
        if (data.Dirbctype['constant'] or data.Dirbctype['piecewise']):
            self.bc.set(self.mesh.Dnodes, self.t0)
            self.semAss.setbc(self.bc.DirichletBc())
        if (self.mesh.anyDirbc):
            self.bc.set(mesh.Dnodes, self.t0)
            self.semAss.setbc(bc.DirichletBc())

    def updateBoundaryData(self):
        if (self.nonHomogeneousDirbc):
            if (self.varDirbc):
                self.bc.set(self.mesh.Dnodes, self.t)
                self.semAss.setbc(self.bc.DirichletBc())
        return None

    def BoundaryFix(self, u):
        self.updateBoundaryData()
        self.semAss.enforceBc(u)#
#         if self.mesh.anyDirbc: u[self.semAss.maskvecL] = self.semAss.bc
#         else: self.semAss.enforceBc(u)#
    
    def updateSourceData(self):
        if (self.sources.var_source):
            self.sources.evaluate(self.t)
            self.fs = self.semAss.scatter(self.sources.sources)
        return None

    def setLsolverType(self, solvertype, userLinearOperator=False):
        self.linearSolver (solvertype, userLinearOperator)
        return None

    def validVector(self, u):
        if any(self.np.isnan(u)):
            return False
        elif any(self.np.isinf(u)):
            return False
        return  (self.np.linalg.norm(u)<self.inf)
#         return  (not(any(self.np.isnan(u))) or self.np.linalg.norm(u)>self.inf)

    def solveLinearDE(self, f):
        if self.linearSolver.matvec:
            return self.linearSolver.solve(f)
        else:
            ug = self.linearSolver.solve(self.semAss.gather(f,False,False))
            return self.semAss.scatter(ug)
        
#     def errorL1(self, u0, t, ue = None):
#         if (ue is None): ue = self.init.compute4ex(t)
#         err = 0.; l = 0
#         for el in self.mesh.elements:
#             du = u0[el] - ue[el]
#             dx  = self.mesh.size(l) / 2; l += 1
#             err += self.gll.gllweights.dot(abs(du)) * dx
#         return err
#         
#     def errorL2(self, u0, t, ue = None):
#         if (ue is None): ue = self.init.compute4ex(t)
#         err = 0.; l = 0
#         for el in self.mesh.elements:
#             du = u0[el] - ue[el]
#             dx  = self.mesh.size(l) / 2; l += 1
#             err += self.gll.gllweights.dot(du**2) * dx
#         return self.np.sqrt(err)
#         
#     def errorLinf(self, u0, t, ue = None):
#         if (ue is None): ue = self.init.compute4ex(t)
#         return self.np.linalg.norm(abs(u0 - ue), self.np.inf)
# 
        
    def errorL1(self, u0, t, ue = None, absVal = True):
        if (ue is None): ue = self.init.compute4ex(t)
        if isinstance(ue,float): ue = ue*self.np.ones(self.np.shape(u0))
        err = 0.; l = 0
        du = abs(u0 - ue) if (absVal) else (u0 - ue)
        for el in self.mesh.elements:
#             du = u0[el] - ue[el]
            dx  = self.mesh.size(l) / 2; l += 1
            err += self.gll.gllweights.dot(du[el]) * dx
#             err += self.gll.gllweights.dot(abs(du)) * dx if (absVal)\
#                     else self.gll.gllweights.dot(du) * dx
        return err
        
    def errorL2(self, u0, t, ue = None):
        if (ue is None): ue = self.init.compute4ex(t)
        if isinstance(ue,float): ue = ue*self.np.ones(self.np.shape(u0))
        err = 0.; l = 0
        for el in self.mesh.elements:
#             print(el)
            du = u0[el] - ue[el]
            dx  = self.mesh.size(l) / 2; l += 1
            err += self.gll.gllweights.dot(du**2) * dx
        return self.np.sqrt(err)
        
    def errorLinf(self, u0, t, ue = None):
        if (ue is None): ue = self.init.compute4ex(t)
        return self.np.linalg.norm(abs(u0 - ue), self.np.inf)
#     ============================================

class CF_Asch (timeStepping):
# working

    def __init__(self, data, gll, mesh, bc, nsteps):

        super().__init__(data, gll, mesh, bc, nsteps)

        from timeStepping import ExpSL

        self.expsl = ExpSL(mesh, gll, data.xb)
#         self.expsl = ExpSL(mesh, gll, data.xb, min(4,self.cf.order))
        self.massAssembly = True
#         self.semAss.assembleConvection(True)
#         self.expsl.C = self.semAss.Cmat

    def feval(self, u, Exp = True, mask = False):
        f = self.semAss.Diffusion(u, mask)
        f = self.semAss.MassInv(-f, self.massAssembly, mask)
        if (self.nonHomogeneousSource): f += self.fs
        self.BoundaryFix(f)
        if not(Exp): return f
        return self.expsl.evaluate(f, -self.uc, self.dt)

    def feval2(self, u, u0, mfac, Exp = True):
        if Exp:  u = self.expsl.evaluate(u, -self.uc, self.dt)
        return ( (u - u0) / mfac )

    def lhsVectors(self, f, expDone = False):
        self.updateSourceData()
        if not(expDone): f = self.expsl.evaluate(f, self.uc, self.dt)
        self.semAss.Mass(f)
        if (self.nonHomogeneousSource): 
            self.updateSourceData()
            f += self.fs
        if (self.nonHomogeneousDirbc): self.updateBoundaryData(self.t)
        return self.semAss.dssum(f)

    def coeff(self, a):
        return abs(a) > self.epsilon

    def evolve(self, u0 = None):
        tn  = self.t0
        u0 = self.semAss.scatter(self.u0)
        steps = [self.dt for _ in range(self.nsteps)]
        nstep = 1
        converged = True
        max_resize = 4; iresize = 0
        mfac_old = 0.; mfac = 0.
        lbeta = len(self.cf.beta)
        lalpha = 0

        if self.storeAll: self.sol.append(self.u0)

        nstages = len(self.cf.AI); 
#         F = [self.np.zeros((self.nv,)) for _ in range(nstages)]
#         U = [self.np.zeros((self.nv,)) for _ in range(nstages)]
#         f = None #ut = None
        instages = nstages
        if (self.cf.order>3): instages -= 1
        if (self.cf.alpha is not None): lapha = len(self.cf.alpha)

        self.linearSolver.displayInfo = self.stepBystepReport 
        AI = self.cf.AI
        AE = self.cf.AE
        beta = self.cf.beta
        alpha = self.cf.alpha
        c = self.cf.c
        b = self.cf.b

        for r in range(10):
#         while (len(steps)>0):

            dt = steps.pop()
            if self.stepBystepReport:
                print('='*33)
                print(" time step %i of %i;\tdt = %f" %(nstep, self.nsteps, dt))
                print('='*33)
            F = [self.np.zeros((self.nv,))] * nstages
            U = [self.np.zeros((self.nv,))] * nstages

            U[0] = u0
            if self.storeAll: u0 = 1.*u0 # renew location
            if self.coeff(AI[1][0]): 
                print("ESDIRK")
                F[0] = self.feval(u0, False); self.massAssembly = False
            for i in range(1,instages):
                self.t = tn + c[i] * dt
                mfac =  1./(dt*AI[i][i])
                if (abs(mfac-mfac_old) > self.epsilon):
                    self.linearSolver.assembleMatrix = True
                    self.semAss.set_massfactor( mfac )
                else:
                    self.linearSolver.assembleMatrix = False
                mfac_old = mfac
                self.dt = dt

                f = U[0] + dt * AI[i].dot(F)
                self.uc = AE[i].dot(U)
                f = self.lhsVectors(f)
                u0 = self.solveLinearDE(f)
                self.BoundaryFix(u0)

                U[i] = u0
                F[i] = self.feval(u0)
                if self.massAssembly: self.massAssembly = False
                converged = converged and self.linearSolver.validSol

            if (self.cf.order>3): # order 4
                self.t = tn + c[-1] * dt
                mfac =  1./(dt*AI[-1][-1])
                if (abs(mfac-mfac_old) > self.epsilon):
                    self.linearSolver.assembleMatrix = True
                    self.semAss.set_massfactor( mfac )
                else:
                    self.linearSolver.assembleMatrix = False
                mfac_old = mfac
                f = U[0] + dt * AI[-1,1:-1].dot(F[1:-1])
                for s in range(lalpha):
                    self.uc = self.cf.alpha[s].dot(U[:-1]) # linear combination
                    f = self.expsl.evaluate(f, self.uc, dt)
                f = self.lhsVectors(f, True)
                u0 = self.solveLinearDE(f)
                self.BoundaryFix(u0)
                U[-1] = u0
                f = self.feval(u0, False)
                for s in reversed(range(lalpha)):
                    self.uc = self.cf.alpha[s].dot(U[:-1]) # linear combination
                    f = self.expsl.evaluate(f, -self.uc, dt)
                converged = converged and self.linearSolver.validSol

            self.t = tn + dt

            if not(self.cf.stifflyAccurateDIRKCF and self.cf.order<=2):
                u0 = U[0] + dt * b.dot(F)  # update using weighted stage vectors
                self.BoundaryFix(u0)
                for s in range(lbeta):
                    self.uc = beta[s].dot(U) # linear combination
                    u0 = self.expsl.evaluate(u0, self.uc, dt)
                self.BoundaryFix(u0)
# 
            converged = (converged and self.validVector(u0))

            if converged:
                iresize = 0
                if self.storeAll: self.sol.append(self.semAss.gather(u0,False,False))
                else: self.sol = self.semAss.gather(u0,False,False)
                if self.exactsol: 
                    errL2 = self.errorL2(self.semAss.gather(u0,False,False), self.t)
                    if self.stepBystepReport:
                        print("Numerical ERROR (L2) at time (t = %f): %f" %(self.t, errL2))
                nstep += 1
            elif (iresize < max_resize):
                steps.append(dt/2)
                steps.append(dt/2)
                u0 = U[0]
                iresize += 1
                print("Step %i rejected!" %nstep)
            else: 
                print("invalid solution vector!")
                print(u0)
                if self.storeAll: 
                    self.sol.append(self.semAss.gather(u0,False,False))
                    while len(steps)>0: dt = steps.pop(); self.sol.append(self.semAss.gather(u0,False,False))
                else: self.sol = self.semAss.gather(u0,False,False)
                print("Cannot reduce time step further! dt = %f" % dt)
                print("invalid solution vector!"); return None

        return None
#     ============================================

    def evolve2(self, u0 = None):
        tn  = self.t0
        u0 = self.semAss.scatter(self.u0)
        steps = [self.dt for _ in range(self.nsteps)]
        nstep = 1
        converged = True
        max_resize = 4; iresize = 0
        mfac_old = 0.; mfac = 0.
        lbeta = len(self.cf.beta)

        if self.storeAll: self.sol.append(self.u0)

        nstages = len(self.cf.AI); 

        self.linearSolver.displayInfo = self.stepBystepReport 
        AI = self.cf.AI
        AE = self.cf.AE
        beta = self.cf.beta
        c = self.cf.c
        b = self.cf.b

#         for r in range(1):
        while (len(steps)>0):

            dt = steps.pop()
            if self.stepBystepReport:
                print('='*33)
                print(" time step %i of %i;\tdt = %f" %(nstep, self.nsteps, dt))
                print('='*33)
            if self.np.isnan(u0).any(): print("loop broken by NaN values at t = %f!" % t); break
            if self.np.isinf(u0).any(): print("loop broken by unbounded values at t = %f!" % t); break
            F = [self.np.zeros((self.nv,))] * nstages
            Y = [self.np.zeros((self.nv,))] * nstages
            mfac =  1./(dt*AI[1][1])
            self.semAss.set_massfactor( mfac )
            self.dt = dt

            self.t = tn + c[1] * dt
            Y[0] = u0
            F[0] = self.feval(Y[0], False, False); self.massAssembly = False
            rhs = u0 + dt * AI[1].dot(F)
            self.uc = AE[1].dot(Y)
            rhs = self.lhsVectors(rhs)
            Y[1] = self.solveLinearDE(rhs)
            self.BoundaryFix(Y[1])
            self.linearSolver.assembleMatrix = False
            if self.massAssembly: self.massAssembly = False
            converged = converged and self.linearSolver.validSol
            F[1] = self.feval(Y[1])

            self.t = tn + c[2] * dt
            rhs = u0 + dt * AI[2].dot(F)
            self.uc = AE[2].dot(Y)
            rhs = self.lhsVectors(rhs)
            Y[2] = self.solveLinearDE(rhs)
            self.BoundaryFix(Y[2])
            if self.massAssembly: self.massAssembly = False
            converged = converged and self.linearSolver.validSol
            F[2] = self.feval(Y[2])

            rhs = u0 + dt * b.dot(F)
            self.uc = beta[0].dot(Y)
            u0 = self.expsl.evaluate(u0, self.uc, dt)
            self.BoundaryFix(u0)
            self.uc = beta[1].dot(Y)
            u0 = self.expsl.evaluate(u0, self.uc, dt)
            self.BoundaryFix(u0)

            self.t = tn + dt

#             if not(self.cf.stifflyAccurateDIRKCF and self.cf.order<=2):
#                 u0 = U[0] + dt * b.dot(F)  # update using weighted stage vectors
#                 self.BoundaryFix(u0)
#                 for s in range(lbeta):
#                     self.uc = beta[s].dot(U) # linear combination
#                     u0 = self.expsl.evaluate(u0, self.uc, dt)
#                 self.BoundaryFix(u0)
# 
            if converged:
                iresize = 0
                if self.storeAll: self.sol.append(self.semAss.gather(u0,False,False))
                else: self.sol = self.semAss.gather(u0,False,False)
                if self.exactsol: 
                    errL2 = self.errorL2(self.semAss.gather(u0,False,False), self.t)
                    if self.stepBystepReport:
                        print("Numerical ERROR (L2) at time (t = %f): %f" %(self.t, errL2))
                nstep += 1
            elif (iresize < max_resize):
                steps.append(dt/2)
                steps.append(dt/2)
                u0 = U[0]
                iresize += 1
                print("Step %i rejected!" %nstep)
            else: 
                print("Cannot reduce time step further! dt = %f" %dt)
                print("invalid solution vector!"); return None
#             U = [self.np.zeros((self.nv,))] * nstages

        return None
#     ============================================

class IMEX_Asch (timeStepping): # derived class

    def __init__(self, data, gll, mesh, bc, nsteps):

        super().__init__(data, gll, mesh, bc, nsteps)

        self.massAssembly = True
#         self.semAss.assembleConvection(True)

    def feval(self, u, Convection=False, mask=False):
        f = self.np.zeros(self.np.shape(u))
        if (Convection):
            f = self.semAss.Convection(u, mask)#, not(self.cf.stifflyAccurateIMEX))
        else:
            f = self.semAss.Diffusion(u, mask)#, not(self.cf.stifflyAccurateIMEX))

        f = self.semAss.MassInv(f, self.massAssembly, not(mask))
        return -f

    def feval2(self, u, u0, mfac): # Diffusion only
        return ( (u - u0) / mfac )

    def lhsVectors(self, f):
        self.updateSourceData()
        self.semAss.Mass(f)
        if (self.nonHomogeneousSource): 
            f += self.fs
        if (self.nonHomogeneousDirbc): self.updateBoundaryData()
        return self.semAss.dssum(f)

    def coeff(self, a):
        return abs(a) > self.epsilon

    def evolve(self, u0 = None):
        tn  = self.t0
        u0 = self.semAss.scatter(self.u0)
        steps = [self.dt for _ in range(self.nsteps)]
        nstep = 1
        converged = True
        max_resize = 4; iresize = 0
        mfac_old = 0.; mfac = 0.

        if self.storeAll: self.sol.append(self.u0)

        nstages = len(self.cf.AI)
        self.linearSolver.displayInfo = self.stepBystepReport 
        AI = self.cf.AI
        AE = self.cf.AE
        c = self.cf.c
        b = self.cf.b
        bt = self.np.sum(self.cf.beta,0)

#         for r in range(1):
        while (len(steps)>0):
            F = [self.np.zeros((self.nv,))] * nstages
            Fc = [self.np.zeros((self.nv,))] * nstages
            U = [self.np.zeros((self.nv,))] * nstages

            dt = steps.pop()
            if self.stepBystepReport:
                print('='*33)
                print(" time step %i of %i;\tdt = %f" %(nstep, self.nsteps, dt))
                print('='*33)

            U[0] = u0
            if self.storeAll: u0 = 1.*u0 # renew location
            Fc[0] = self.feval(u0, Convection=True)
            if self.coeff(AI[1][0]): 
                print("ESDIRK")
                F[0] = self.feval(u0); self.massAssembly = False
            for i in range(1,nstages):
                self.t = tn + c[i] * dt
                mfac =  1./(dt*AI[i][i])
                if (abs(mfac-mfac_old) > self.epsilon):
                    self.linearSolver.assembleMatrix = True
                    self.semAss.set_massfactor( mfac )
                else:
                    self.linearSolver.assembleMatrix = False
                mfac_old = mfac
                self.dt = dt

                f = U[0] + dt * ( AI[i].dot(F) + AE[i].dot(Fc) )
                f = self.lhsVectors(f)
                u0 = self.solveLinearDE(f)
                self.BoundaryFix(u0)

                U[i] = u0
                F[i] = self.feval(u0)
                Fc[i] = self.feval(u0, Convection=True)
                if self.massAssembly: self.massAssembly = False
                converged = converged and self.linearSolver.validSol

            self.t = tn + dt

            if not(self.cf.stifflyAccurateIMEX):
                u0 = U[0] + dt * ( b.dot(F) + bt.dot(Fc) )  # update using weighted stage vectors
                self.BoundaryFix(u0)

            converged = (converged and self.validVector(u0))

            if converged:
#                 iresize = 0
                if self.storeAll: self.sol.append(self.semAss.gather(u0,False,False))
                else: self.sol = self.semAss.gather(u0,False,False)
                if self.exactsol: 
                    errL2 = self.errorL2(self.semAss.gather(u0,False,False), self.t)
                    if self.stepBystepReport:
                        print("Numerical ERROR (L2) at time (t = %f): %f" %(self.t, errL2))
                nstep += 1
            elif (iresize < max_resize):
                print("Step %i rejected!" %nstep)
                self.dt /= 2; self.nsteps *= 2; nstep = 1
                u0 = self.semAss.scatter(self.u0)
                steps = [self.dt for _ in range(self.nsteps)]
                self.sol = []
                if self.storeAll:
                    self.sol.append(self.u0)
#                 steps.append(dt/2)
#                 steps.append(dt/2)
#                 u0 = U[0]
                iresize += 1
                converged = True
            else: 
                print("invalid solution vector!")
                print(u0)
                if self.storeAll: 
                    self.sol.append(self.semAss.gather(u0,False,False))
                    while len(steps)>0: dt = steps.pop(); self.sol.append(self.semAss.gather(u0,False,False))
                else: self.sol = self.semAss.gather(u0,False,False)
                print("Cannot reduce time step further! dt = %f" % dt)
                print("invalid solution vector!"); return None

        return None

#     ============================================

class ExactSol (timeStepping):

    from scipy.integrate import solve_ivp

    def __init__(self, data, gll, mesh, bc):

        super().__init__(data, gll, mesh, bc, 1)

        self.MassAssembly(self.u0)

    def MassAssembly(self,u):
        f = self.semAss.scatter(u)
        f = self.semAss.MassInv(f,True)
        return None

    def feval(self, t, u):
        uL = self.semAss.scatter(u)
        f = self.semAss.Diffusion(uL, False)
        f += self.semAss.Convection(uL, False)
        f = self.semAss.MassInv(-f)
        self.t = t
        if (self.nonHomogeneousSource): 
            self.updateSourceData(); f += self.fs
#         self.BoundaryFix(f)
        return self.semAss.gather(f,False,False)

#     ============================================

class DIRK_Asch (timeStepping): # derived class

    def __init__(self, data, gll, mesh, bc, nsteps):

        super().__init__(data, gll, mesh, bc, nsteps)

    def feval2(self, u, f):
        fc = self.np.zeros(self.np.shape(f))
        f -= fc
        self.semAss.MassInv(f)
        if (self.nonHomogeneousSource): f = self.fs

    def evolve(self, u0):
        tn  = self.t0
        steps = [self.dt for _ in range(self.nsteps)]
        nstep = 1
        converged = False
        max_resize = 4; iresize = 0
        dt_old = 0.

        F = []; U = []
        nstages = len(self.cf.AI)

        self.linearSolver.printSolverInfo()

        while (len(steps)>0):

            f = None; ut = None

            dt = steps.pop()
            print('='*50)
            print(" TIME STEP %i OF %i;\tSTEPSIZE dt = %f" %(nstep, self.nsteps, dt))
            print('='*50)

            U.append(u0)

            for i in range(nstages):

                print("=== STAGE %i ===" % i+1)
                self.t = tn + dt * self.cf.c[i]

                if (i>0): f = u0 + dt * self.cf.AI[i,:i].dot(F[:i])
                else: f = u0

                self.cdt = self.cf.AI[i][i] * dt
                self.semAss.set_timefactor(1./self.cdt)
                f = self.lhsVectors(f)

                u0 = self.solveLinearDE(f, self.cdt)
                converged = self.linearSolver.validSol

                self.constantSteps = (abs(dt_old - dt) < self.epsilon)
                dt_old = self.cdt

                if (not(self.cf.stifflyAccurateDIRK) or i < nstages-1):
                    self.feval2(u0, f)
                    F.append(f)

                self.t = tn + dt
                tn  = self.t

                if (not(cf.stifflyAccurateDIRK)):
                    print("=== UPDATE STAGE ===")
                    u0 += dt * self.cf.b.dot(F)
#                     self.BoundaryFix(u0, self.t)

                converged = (converged and self.validVector(u0))
                errL2 = self.errorL2(u0, self.mesh, self.gll, self.init, self.t)
                print("Numerical ERROR (L2) at time (t = %f): %f" %(self.t, errL2))

                nstep += 1
                U = []; F = []
                self.sol = u0
                return None
#     ============================================

    def evolve2(self, u0 = None):
        tn  = self.t0
        u0 = self.semAss.scatter(self.u0)
        steps = [self.dt for _ in range(self.nsteps)]
        nstep = 1
        converged = True
        max_resize = 4; iresize = 0
        mfac_old = 0.; mfac = 0.
        lbeta = len(self.cf.beta)
        lalpha = 0

        if self.storeAll: self.sol.append(self.u0)

        nstages = len(self.cf.AI); 
        F = [self.np.zeros((self.nv,)) for _ in range(nstages)]
        U = [self.np.zeros((self.nv,)) for _ in range(nstages)]
        f = None #ut = None
        istart = self.cf.istart
        jstart = self.cf.jstart
        instages = nstages
        if (self.cf.order>3): instages -= 1
        if (self.cf.alpha is not None): lapha = len(self.cf.alpha)

        self.linearSolver.displayInfo = self.stepBystepReport 

#         for r in range(60):
        while (len(steps)>0):

            dt = steps.pop()
            if self.stepBystepReport:
                print('='*33)
                print(" time step %i of %i;\tdt = %f" %(nstep, self.nsteps, dt))
                print('='*33)

            U[0] = u0
            if (jstart==0): F[0] = self.feval(u0, False); #self.massAssembly = False
            for i in range(1,instages):
                self.t = tn + self.cf.c[i] * dt
                mfac =  1./(dt*self.cf.AI[i][i])
                if (abs(mfac-mfac_old) > self.epsilon):
                    self.linearSolver.assembleMatrix = True
                    self.semAss.set_massfactor( mfac )
                else:
                    self.linearSolver.assembleMatrix = False
                mfac_old = mfac
                self.dt = dt

                f0 = U[0] + dt * self.cf.AI[i,jstart:i].dot(F[jstart:i])
                self.uc = self.cf.AE[i,:i].dot(U[:i])
                f = self.lhsVectors(f0)
                u0 = self.solveLinearDE(f)
                self.BoundaryFix(u0)

                U[i] = u0
                F[i] = self.feval2(u0 , f0, mfac)
#                 if self.massAssembly: self.massAssembly = False
                converged = converged and self.linearSolver.validSol

            if (self.cf.order>3): # order 4
                self.t = tn + self.cf.c[-1] * dt
                mfac =  1./(dt*self.cf.AI[-1][-1])
                if (abs(mfac-mfac_old) > self.epsilon):
                    self.linearSolver.assembleMatrix = True
                    self.semAss.set_massfactor( mfac )
                else:
                    self.linearSolver.assembleMatrix = False
                mfac_old = mfac
                f0 = U[0] + dt * self.cf.AI[-1,1:-1].dot(F[1:-1])
                for s in range(lalpha):
                    self.uc = self.cf.alpha[s].dot(U[:-1]) # linear combination
                    f0 = self.expsl.evaluate(f0, self.uc, dt)
                f = self.lhsVectors(f0, True)
                u0 = self.solveLinearDE(f)
                self.BoundaryFix(u0)
                U[-1] = u0
                f = self.feval2(u0 , f0, mfac, False)
                for s in reversed(range(lalpha)):
                    self.uc = self.cf.alpha[s].dot(U[:-1]) # linear combination
                    f = self.expsl.evaluate(f, -self.uc, dt)
                converged = converged and self.linearSolver.validSol

            self.t = tn + dt

            if not(self.cf.stifflyAccurateDIRKCF and self.cf.order<=2):
                u0 = U[0] + dt * self.cf.b[jstart:].dot(F[jstart:])  # update using weighted stage vectors
                self.BoundaryFix(u0)
                for s in range(lbeta):
                    self.uc = self.cf.beta[s].dot(U) # linear combination
                    u0 = self.expsl.evaluate(u0, self.uc, dt)
                self.BoundaryFix(u0)

            if converged:
                iresize = 0
                if self.storeAll: self.sol.append(self.semAss.gather(u0,False,False))
                else: self.sol = self.semAss.gather(u0,False,False)
                if self.exactsol: 
                    errL2 = self.errorL2(self.semAss.gather(u0,False,False), self.t)
                    if self.stepBystepReport:
                        print("Numerical ERROR (L2) at time (t = %f): %f" %(self.t, errL2))
                nstep += 1
            elif (iresize < max_resize):
                steps.append(dt/2)
                steps.append(dt/2)
                u0 = U[0]
                iresize += 1
                print("Step %i rejected!" %nstep)
            else: 
                print("Cannot reduce time step further! dt = %f" %dt)
                print("invalid solution vector!"); return None

        return None


#     ============================================

class RKDG (timeStepping): # derived class

    def __init__(self, data, gll, mesh, bc, nsteps):

        super().__init__(data, gll, mesh, bc, nsteps)

        from Assembly import AssemblyDG
        from Methods import RKDG as RKDG1

        self.semAss = AssemblyDG(mesh, gll, data.Re)
        self.cf = RKDG1(gll.np, data.orderRK)
        self.semAss.initd = self.init
        self.semAss.flowtype = data.flowtype
#         self.semAss.limit = True
        self.semAss.order = data.orderRK
        self.semAss.ltype = data.DGlimiter 
        self.semAss.m0 = data.DGlimiterMinMax[0]
        self.semAss.M0 = data.DGlimiterMinMax[1]
        self.semAss.src = bc
        self.semAss.limit = data.DGlimiter['PP'] or \
                data.DGlimiter['MPP']
        if self.semAss.limit:
            self.semAss.ltype = data.DGlimiter
            self.semAss.m0 = data.DGlimiterMinMax[0]
            self.semAss.M0 = data.DGlimiterMinMax[1]

    def evolve(self):
        tn  = self.t0
        self.t = tn
        u0 = self.semAss.scatter(self.u0)
        self.semAss.bc = self.u0[self.mesh.Dnodes]
        steps = [self.dt]*self.nsteps
        nstep = 1
        converged = True
        max_resize = 4; iresize = 0
        # set controls
        nonlinearLimiter = self.semAss.flowtype['nonlinear'] or self.mesh.dG0
        temporary = self.semAss.flowtype['temporal']
        mpplim = self.semAss.limit
        anyDirbc = self.mesh.anyDirbc

        if self.mesh.dG0: u0 = self.semAss.limiter(u0)
        if self.storeAll: self.sol.append(self.u0)

        nstages = len(self.cf.alpha); 
        f = None #ut = None

        F = [self.np.zeros((self.nv,))] * nstages
        U = [self.np.zeros((self.nv,))] * nstages

        print('-'*60)
#         for r in range(305):
        while (len(steps)>0):

            f = None; ut = None

            dt = steps.pop()
            if self.stepBystepReport:
                print(" TIME STEP %i OF %i;\tdt = %f;\tt = %f" %(nstep, self.nsteps, dt, self.t+dt))

#             L2Stab = self.errorL2(u0, tn, 0.)

            for i in range(nstages):

                self.cdt = self.cf.c[i] * dt
                self.t = tn + self.cdt
                U[i] = u0
                f = self.semAss.Derivative(u0, self.t, True)
                F[i] =  self.semAss.numericalFlux(f, u0, dt, self.t)

                u0 = self.cf.alpha[i].dot(U)\
                     + dt * self.cf.beta[i].dot(F)

#                 u0 = self.semAss.limiter(u0)
                if anyDirbc: self.semAss.insertbc(u0)
                if nonlinearLimiter : u0 = self.semAss.limiter(u0)
                if mpplim: u0 = self.semAss.mpplim(u0)

            self.t = tn + dt
            tn  = self.t
            
#             L2Stab -= self.errorL2(u0, tn, 0.)
#             if L2Stab < -1e-2:
#                 print('L2-stability = %8.4E failed!' % L2Stab)

            converged = (converged and self.validVector(u0))

            if converged:
                if self.storeAll: 
                    self.sol.append(u0)
                else: self.sol = u0
                if self.exactsol: 
                    errL2 = self.errorL2(u0, self.t)
                    print("\tNumerical ERROR (L2) at time (t = %f): %f" %(self.t, errL2))
            else: 
                print("invalid solution vector!")#; self.sol.append(u0)
                print(u0)
                if self.storeAll: 
                    u0 = self.sol[-1]
                    self.nstepsdone = self.nsteps
                    self.sol.append(u0)
                    while len(steps)>0: dt = steps.pop(); self.sol.append(u0)
                return None


            self.nstepsdone = nstep
            nstep += 1
        if not(self.storeAll):self.sol = u0
        print('-'*60)

        return None
#     ============================================

class SLRKDG (RKDG): # derived class

    def __init__(self, data, gll, mesh, bc, nsteps):

        super().__init__(data, gll, mesh, bc, nsteps)

        from timeStepping import ExpSLDG

        self.expsl = ExpSLDG(mesh, gll, data.xb, self.semAss.initd, data.flowtype)# data.orderRK)
        self.semAss.tol = 1e-15
        self.semAss.limit = data.DGlimiter['PP'] or \
                data.DGlimiter['MPP']
        if self.semAss.limit:
            self.semAss.ltype = data.DGlimiter
            self.semAss.m0 = data.DGlimiterMinMax[0]
            self.semAss.M0 = data.DGlimiterMinMax[1]

    def continuousVelocity(self, f, u, t, fopt):
        if self.semAss.flowtype['nonlinear']:
            F = self.semAss.Godunov(u, t, fopt)
            l = self.mesh.edges2[:,0]
            r = self.mesh.edges2[:,1]
            f[self.mesh.edges[l]] = F[l]
            f[self.mesh.edges[r]] = F[r]
        return None

    def xticVelocity(self, u0, tn, psi = 0.0, theta = False):
            phi = 1.0
            g = None
            if self.semAss.flowtype['nonlinear']:
                if theta:
                    phi = max(abs(self.semAss.slopeRatio(u0)))
                    phi = (2 / self.np.pi) * self.np.arctan(phi)
                fzero = not(theta) and abs(phi) < self.epsilon
                gzero = not(theta) and abs(phi-1) < self.epsilon
                f = 0. if (fzero) else self.init.compute4a(tn, u0)  # f(u)/u
                if not(fzero): self.continuousVelocity(f, u0, tn, 1)
                g = 0.*u0 if (gzero) else self.init.compute4dflux(u0, False, tn) # f'(u)
                if not(gzero): self.continuousVelocity(g, u0, tn, 2)
                g = (1-phi) * g + phi * f
            else:
                g = self.init.compute4a(tn, u0)  # f(u)/u
            self.semAss.valida = (psi < 1)
            return (1-psi)*g

    def evolve(self):
        tn  = self.t0
        self.t = tn
        u0 = self.u0
        self.semAss.bc = self.u0[self.mesh.Dnodes]
        steps = [self.dt]*self.nsteps
        nstep = 1
        converged = True
        max_resize = 4; iresize = 0
        x = self.expsl.x
        # set controls
        nonlinearLimiter = self.semAss.flowtype['nonlinear'] or self.mesh.dG0
        temporary = self.semAss.flowtype['temporal']
        mpplim = self.semAss.limit
        anyDirbc = self.mesh.anyDirbc

        if self.mesh.dG0: u0 = self.semAss.limiter(u0)
        if self.storeAll: self.sol.append(self.u0)

        nstages = len(self.cf.alpha); 
        f = None #ut = None

        F = [self.np.zeros((self.nv,))] * nstages
        U = [self.np.zeros((self.nv,))] * nstages

        print('-'*60)
#         for r in range(1):
        while (len(steps)>0):

            f = None; ut = None

            dt = steps.pop()
            if self.stepBystepReport:
                print(" TIME STEP %i OF %i;\tdt = %f;\tt = %f" %(nstep, self.nsteps, dt, self.t+dt))

            uc = self.xticVelocity(u0, tn)
            self.semAss.a = uc 
#             self.expsl.t = tn + dt
            #===== Attempt 1: initial idea =====

            f = self.expsl.evaluate(u0, uc, dt, tn+dt, False)
            g = self.semAss.Derivative(f, tn)#, Flux=False, dx=True)
            l = self.mesh.edges[self.mesh.edges2]
            g[l[:,0]] += f[l[:,0]] # left edge
            g[l[:,1]] -= f[l[:,1]] # right edge

            self.semAss.MassInv(g, self.mesh.anyDirbc)

#             massConsv = self.errorL1(g, tn+dt, 0., False)
#             if abs(massConsv)>100*self.epsilon:
#                 print('mass conservation error = %8.4E failed!' % massConsv)

#             L2Stab = self.errorL2(u0, tn, 0.)

#             u0 = u0 + g
#             if nonlinearLimiter : u0 = self.semAss.limiter(u0)

#             if not(self.mesh.dG0):
#                 F = 0.5*(u0[l[1:,0]] + u0[l[:-1,1]])
#                 u0[l[1:,0]] = F
#                 u0[l[:-1,1]] = F
#                 Fl = 0.5*(u0[l[0,0]]+u0[l[-1,1]])
#                 u0[l[0,0]] = Fl
#                 u0[l[-1,1]] = Fl

            for i in range(nstages):

                self.cdt = self.cf.c[i] * dt
                self.t = tn + self.cdt
                if temporary : self.semAss.a = self.xticVelocity(u0, self.t)
                U[i] = u0
                f = self.semAss.Derivative(u0, self.t, True)
                f =  self.semAss.numericalFlux(f, u0, dt, self.t)

                F[i] = dt * f + g
                u0 = self.cf.alpha[i].dot(U) + self.cf.beta[i].dot(F)
                if anyDirbc: self.semAss.insertbc(u0)
                if nonlinearLimiter : u0 = self.semAss.limiter(u0)
                if mpplim : u0 = self.semAss.mpplim(u0)
#                 if self.mesh.anyDirbc: self.semAss.insertbc(u0)
                
#             #===== Attempt 2: Strang splitting =====
#             for i in range(nstages):
#                 self.cdt = self.cf.c[i] * dt/2
#                 self.t = tn + self.cdt
#                 if self.semAss.flowtype['temporal']: self.semAss.a = self.xticVelocity(u0, self.t)#, 1)
#                 U[i] = u0
#                 f = self.semAss.Derivative(u0, self.t, True)
#                 f =  self.semAss.numericalFlux(f, u0, dt/2, self.t)
#                 F[i] = dt/2 * f
#                 u0 = self.cf.alpha[i].dot(U) + self.cf.beta[i].dot(F)
#                 if self.semAss.flowtype['nonlinear']: u0 = self.semAss.limiter(u0)
# 
#             f = self.expsl.evaluate(u0, uc, dt)
#             g = self.semAss.Derivative(f, tn)
#             l = self.mesh.edges[self.mesh.edges2]
#             g[l[:,0]] += f[l[:,0]] # left edge
#             g[l[:,1]] -= f[l[:,1]] # right edge
#             self.semAss.MassInv(g)
#             u0 += g
#             if self.semAss.flowtype['nonlinear']: u0 = self.semAss.limiter(u0)
# 
#             for i in range(nstages):
#                 self.cdt = self.cf.c[i] * dt/2
#                 self.t = tn + self.cdt
#                 if self.semAss.flowtype['temporal']: self.semAss.a = self.xticVelocity(u0, self.t)#, 1)
#                 U[i] = u0
#                 f = self.semAss.Derivative(u0, self.t, True)
#                 f =  self.semAss.numericalFlux(f, u0, dt/2, self.t)
#                 F[i] = dt/2 * f
#                 u0 = self.cf.alpha[i].dot(U) + self.cf.beta[i].dot(F)
#                 if self.semAss.flowtype['nonlinear']: u0 = self.semAss.limiter(u0)

#             #===== Attempt 3: uc every stage =====
# 
#             for i in range(nstages):
# 
#                 self.cdt = self.cf.c[i] * dt
#                 self.t = tn + self.cdt
# 
#                 uc = self.xticVelocity(u0, self.t)#, 1)
#                 self.semAss.a = uc 
#                 self.expsl.t = self.t
#                 f = self.expsl.evaluate(u0, uc, dt)
#                 g = self.semAss.Derivative(f, self.t)
#                 l = self.mesh.edges[self.mesh.edges2]
#                 g[l[:,0]] += f[l[:,0]] # left edge
#                 g[l[:,1]] -= f[l[:,1]] # right edge
#                 self.semAss.MassInv(g)
# 
#                 U[i] = u0
#                 f = self.semAss.Derivative(u0, self.t, True)
#                 f =  self.semAss.numericalFlux(f, u0, dt, self.t)
# #                 print(self.np.linalg.norm(f))
#                 F[i] = dt * f + g
#                 u0 = self.cf.alpha[i].dot(U) + self.cf.beta[i].dot(F)
# 
#                 if self.semAss.flowtype['nonlinear']: u0 = self.semAss.limiter(u0)
#                 if self.mesh.anyDirbc: self.semAss.insertbc(u0)

            self.t = tn + dt
            tn  = self.t
            
#             L2Stab -= self.errorL2(u0, tn, 0.)
#             if L2Stab < -1e-2:
#                 print('L2-stability = %8.4E failed!' % L2Stab)

#             J = self.semAss.computeJumps(u0)
#             if (abs(J)>1e-15).any(): print('Jumps!')

            converged = (converged and self.validVector(u0))

            if converged:
                if self.storeAll: 
                    self.sol.append(u0)
                else: self.sol = u0
                if self.exactsol: 
                    errL2 = self.errorL2(u0, self.t)
                    print("\tNumerical ERROR (L2) at time (t = %f): %f" %(self.t, errL2))
            else: 
                print("invalid solution vector!")#; self.sol.append(u0)
                print(u0)
                if self.storeAll: 
                    u0 = self.sol[-1]
                    self.nstepsdone = self.nsteps
                    self.sol.append(u0)
                    while len(steps)>0: dt = steps.pop(); self.sol.append(u0)
                return None

            self.nstepsdone = nstep
            nstep += 1
        if not(self.storeAll):self.sol = u0
        print('-'*60)

        return None
#     ============================================

class ExpSL:


    def __init__(self, mesh, gll, xb, order = 4):
        self.N = mesh.getN()
        self.Ne = mesh.getNe()
        self.np = gll.np
        self.nv = mesh.elements.size
        self.mesh = mesh
        self.xb = xb # boundary points
        self.x = self.np.reshape(mesh.nodes[mesh.elements],(self.nv,))
        self.t = 0.
        self.gll = gll
        self.uc = None
        self.foot_loc = None
        self.order = order # order of RK method to solve xtics
        self.C = None
#         self.flowtype = {}
        self.constantFlow  = False
        self.nonlinearFlow  = True
        self.init = None
        self.nsteps = 1
        self.epsilon = 1e-12
        self.anyPeriodic = mesh.anyPeriodic
        self.Bel = None
        self.Bel = self.mesh.Bedges[:,-1] #  boundary elements

        from scipy.integrate import solve_ivp
        self.solve = solve_ivp

    def evaluate(self, u0, uc, dt, t = None, nonExact=True):

#         xi = self.x
#         self.uc = uc
#         unknownFlowtype = (self.flowtype is None)
#         nonlinear = (self.init is None) or (True if unknownFlowtype else self.flowtype['nonlinear'])
#         constantFlow = False if unknownFlowtype else self.flowtype['constant']
#         t = self.t if (t is None) else t
#         
#         if constantFlow:
#             xi = xi - dt*self.uc
#             self.uc = u0
#             return ( self.InterpSL(xi) if (nonExact) else self.Integrate(xi) )
# 
#         def feval(t, x):
#             return self.InterpSL(x) if (nonlinear) else self.init.compute4A(x, t)
# #         sol = self.solve(feval, [t,t-dt], xi, method='Radau',atol=1e-8, rtol=1e-8)
#         sol = self.solve(feval, [t,t-dt], xi, method='RK23',atol=1e-0, rtol=1e-0)
#         xi = sol.y[:,-1]
# 
#         self.uc = u0
#         return ( self.InterpSL(xi) if (nonExact) else self.Integrate(xi) )

        xi = self.x; x1 = None
        self.uc = uc
        dx = None; ds = None; F = None
        steps = [dt/self.nsteps] * self.nsteps; nsteps = 0
        nonlinear = self.nonlinearFlow
        if (t is None): t = self.t

        if (self.constantFlow):
            ds = steps.pop()
            dx = ds*self.uc
            xi = xi-dx
            nsteps += 1; t -= ds
#             while (len(steps)>0):
#                 ds = steps.pop()
#                 dx = ds*self.uc
#                 xi = xi-dx
#                 nsteps += 1; t -= ds
#                 if self.ResizeSteps(ds, steps, nsteps,self.np.linalg.norm(dx, ord=self.np.inf)): xi += dx
        elif (self.order==1): # Euler
            while (len(steps)>0):
                ds = steps.pop()
                F = self.InterpSL(xi) if (nonlinear) else self.init.compute4A(xi, t-ds)
                dx = ds*F
                xi = xi-dx;
                nsteps += 1; t -= ds
#                 if self.ResizeSteps(ds, steps, nsteps,self.np.linalg.norm(dx, ord=self.np.inf)): xi += dx; t += ds
        elif (self.order==2): # midpoint rule
            b = self.np.array([0.,1.])
            c = [0.,.5]
            F = self.np.zeros((2,self.nv))
            while (len(steps)>0):
                ds = steps.pop()
                F[0] = self.InterpSL(xi) if (nonlinear) else self.init.compute4A(xi, t-c[0]*ds)
                x1 = xi - (c[1]*ds)*F[0]
                F[1] = self.InterpSL(x1) if (nonlinear) else self.init.compute4A(x1, t-c[1]*ds)
                dx = ds * b.dot(F)
                xi = xi-dx;
                nsteps += 1; t -= ds
#                 if self.ResizeSteps(ds, steps, nsteps,self.np.linalg.norm(dx, ord=self.np.inf)): xi += dx; t += ds
        elif (self.order==3): # order 3 rule
            b = self.np.array([1./6,2./3,1./6])
            c = [.5,-1.,2.]
            F = self.np.zeros((3,self.nv))
            while (len(steps)>0):
                ds = steps.pop()
                F[0] = self.InterpSL(xi) if (nonlinear) else self.init.compute4A(xi, t-c[0]*ds)
                x1 = xi - (c[0]*ds)*F[0]
                F[1] = self.InterpSL(x1) if (nonlinear) else self.init.compute4A(x1, t-c[1]*ds)
                x1 = xi - (c[1]*ds)*F[0] - (c[2]*ds)*F[1]
                F[2] = self.InterpSL(x1) if (nonlinear) else self.init.compute4A(x1, t-c[2]*ds)
                dx = ds * b.dot(F)
                xi = xi-dx; 
                nsteps += 1; t -= ds
#                 if self.ResizeSteps(ds, steps, nsteps,self.np.linalg.norm(dx, ord=self.np.inf)): xi += dx; t += ds
        else:
            b = self.np.array([1./6, 1./3, 1./3, 1./6])
            c = [0.,.5,.5,1.]
            F = self.np.zeros((4,self.nv))
            while (len(steps)>0):
                ds = steps.pop()
                F[0] = self.InterpSL(xi) if (nonlinear) else self.init.compute4A(xi, t-c[0]*ds)
                x1 = xi - (c[1]*ds)*F[0]
                F[1] = self.InterpSL(x1) if (nonlinear) else self.init.compute4A(x1, t-c[1]*ds)
                x1 = xi - (c[2]*ds)*F[1]
                F[2] = self.InterpSL(x1) if (nonlinear) else self.init.compute4A(x1, t-c[2]*ds)
                x1 = xi - ds*F[2]
                F[3] = self.InterpSL(x1) if (nonlinear) else self.init.compute4A(x1, t-c[3]*ds)
                dx = ds * b.dot(F)
                xi = xi-dx; 
                nsteps += 1; t -= ds
#                 if self.ResizeSteps(ds, steps, nsteps,self.np.linalg.norm(dx, ord=self.np.inf)): xi += dx; t += ds
        self.uc = u0
        return ( self.InterpSL(xi) if (nonExact) else self.Integrate(xi) )

#     def ResizeSteps(self, ds, steps, nsteps, dx):
#             if (dx > 1.):
#                 print("too large step for solving xtics")
#                 steps.append(ds/2)
#                 steps.append(ds/2)
#                 nsteps -= 1
#                 return True
#             return False

    def Integrate(self, xi):
        return None

    def InterpSL(self, x, u = None):
        self.boundaryfix(x)
        Iu = self.np.zeros((self.nv,))
        if (u is not None): self.uc = u

#         sbool = abs(x[1:]-x[:-1]) < self.epsilon
#         x[1:][sbool] += 10*self.epsilon
#         x[:-1][sbool] -= 10*self.epsilon

        l = 0
        self.foot_loc = self.mesh.contains(x)

        if not(self.anyPeriodic):
            self.foot_loc[x<=self.xb[0]] = self.Bel[0]
            self.foot_loc[x>=self.xb[1]] = self.Bel[1]

        for el in self.mesh.elementsL:
            here = (self.foot_loc == l)
#             here = self.mesh.inclusion(x, l)
            if here.any():
                lx = self.gll.LagrangeBasisMatrix2(self.T(x[here],l))
                Iu[here] = lx.dot(self.uc[el])
            l += 1
        return Iu

    def T(self,x,el): # affine transformation into [-1,1]
        l = self.mesh.edges[self.mesh.edges2[el][0]]
        return (x-self.mesh.nodes[l])*(2/self.mesh.size(el))-1.

#     def boundaryfix(self, x):
#         xl = self.xb[0]; xr = self.xb[1]
#         if isinstance(x,float):
#             if self.anyPeriodic:
#                 Lx = xr - xl
#                 if x < xl or x > xr: x = self.np.mod(x-xl, Lx) + xl
#             else:
#                 if x < xl: x = xl
#                 if x > xr: x = xr
#             return x
#         if (self.anyPeriodic):
#             Lx = xr - xl
#             sbool = (x < xl) | (x > xr)
#             x[sbool] = self.np.mod(x[sbool]-xl, Lx) + xl
#         else: 
#             lbool = (x < xl); rbool = (x > xr)
#             if any(lbool): x[lbool] = xl
#             if any(rbool): x[rbool] = xr
#         return x

    def boundaryfix(self, x):
        xl = self.xb[0]; xr = self.xb[1]
        if isinstance(x,float):
            if self.anyPeriodic:
                Lx = xr - xl
                if x < xl : x = xr - self.np.mod(xl-x, Lx)
                if x > xr : x = xl + self.np.mod(x-xr, Lx)
            else:
                if x < xl: x = xl
                if x > xr: x = xr
            return x
        if (self.anyPeriodic):
            Lx = xr - xl; left = (x < xl); right = (x > xr)
            if any(left): x[left] = xr - self.np.mod(xl-x[left], Lx)
            if any(right): x[right] = xl + self.np.mod(x[right]-xr, Lx)
#         else: 
#             left = (x < xl); right = (x > xr)
#             if any(left): x[left] = xl
#             if any(right): x[right] = xr
        return x

#     ============================================

class ExpSLDG (ExpSL):

#     import time

    def __init__(self, mesh, gll, xb, init=None, flowtype=None, order = 4):

        super().__init__(mesh, gll, xb, order)

        if (flowtype is not None):
            self.constantFlow = flowtype['constant']
            self.nonlinearFlow = flowtype['nonlinear']

        self.init = init
        self.epsilon = 1e-15
        self.loc = self.np.zeros((self.nv,),dtype=int) # arrival locations
        self.GT = 0.
        el = 0
        for K in mesh.elements: self.loc[K] = el; el += 1
        self.Bel = self.mesh.Bedges[:,-1] #  boundary elements


    def Integrate(self, x): 

#         t = self.time.time()
#         nbool = abs(x[1:]-x[:-1]) < self.epsilon
#         x[1:][nbool] += self.epsilon
#         x[:-1][nbool] -= self.epsilon

        dx = self.x - x
        left = (dx > 0)
        tinyDx = abs(dx) < self.epsilon

        self.foot_loc = self.np.full((self.nv,), self.nv, dtype=int)
        out = ( x < self.xb[0] ) | ( x > self.xb[1] )
        self.foot_loc[~out] = self.mesh.contains(x[~out])

        Iu = self.np.zeros((self.nv,))
        if self.anyPeriodic:
            self.GT = self.Partition(self.xb[0], self.xb[1], self.Bel[0], self.Bel[-1])
            period = self.np.zeros((self.nv,))
            if out.any():
                L = self.xb[1] - self.xb[0]
                outside = (out & left)
                if outside.any(): period[outside] = 1. + (self.xb[0] - x[outside]) // L
                outside = (out & ~left)
                if outside.any(): period[outside] = 1. + (x[outside] - self.xb[1]) // L
                x = self.boundaryfix(x)
                self.foot_loc[out] = self.mesh.contains(x[out])

            for l in range(self.nv):
                if tinyDx[l]:  continue
                Iu[l] =  self.int_ab(x[l], l, left[l], period[l])
#             el = 0
#             for K in self.mesh.elements: # suitable for parallel implementation
#                 Iu[K] =  self.int_ab(x[K], self.mesh.nodes[K], self.foot_loc[K], el, left[K], period[K])
#                 el += 1
        else:
            outside = (out & left)
            self.foot_loc[outside] = self.Bel[0]
            outside = (out & ~left)
            self.foot_loc[outside] = self.Bel[-1]
            for l in range(self.nv):
                if tinyDx[l]:  continue
                Iu[l] =  self.int_ab(x[l], l, left[l])
#             el = 0
#             for K in self.mesh.elements: # suitable for parall implementation
#                 Iu[K] =  self.int_ab(x[K], self.mesh.nodes[K], self.foot_loc[K], el, left[K])
#                 el += 1
#         print("elapsed time:", self.time.time() - t)

        return Iu


    def int_ab(self, xt, l, left=True, period = 0):

        x = self.x[l]
        foot = self.foot_loc[l]
        el = self.loc[l]

        I = self.Partition(xt, x, foot, el) if (left)\
                else self.Partition(x, xt, el, foot)

        I += period * self.GT

#         I = 0.
#         if (period>0):
#             if (left):
#                 I += self.Partition(self.xb[0], x, self.Bel[0], el)
#                 I += self.Partition(xt, self.xb[-1], foot, self.Bel[-1])
#             else:
#                 I += self.Partition(self.xb[0], xt, self.Bel[0], foot)
#                 I += self.Partition(x, self.xb[-1], el, self.Bel[-1])
#             I += (period-1) * self.GT
#         else:
#             I = self.Partition(xt, x, foot, el) if (left)\
#                     else self.Partition(x, xt, el, foot)

        return (I if left else -I)


    def Partition(self, xt, x, foot, el):

        I = 0.
        ordered = (xt < x)
        if not(ordered): # swap
            tmp = xt; xt = x; x = tmp
            temp = foot; foot = el; el = temp
        xl = xt

        for er in range(foot,el+1):
            xr = self.mesh.nodes[self.mesh.edges[self.mesh.edges2[er,1]]]
            I += self.GLLsum(xl, min(xr, x), er)
            xl = xr

        return (I if ordered else -I)

    def GLLsum(self, xl, xr, el):

        dx = xr - xl

        if (abs(dx) < self.epsilon): return 0.

        K = self.mesh.elements[el]
        I = 0.
        
        if self.mesh.dG0:
            I = 2 * self.uc[K][0]
        elif K.size < 3:
            x0 = self.mesh.nodes[K[0]]
            I = 2*self.uc[K][0] + (self.uc[K][1]-self.uc[K][0])\
                        * (xl + xr - 2*x0) / self.mesh.size(el)
        else:
            x = xl + (dx/2)*(self.gll.gllnodes + 1)
            Mx = self.gll.LagrangeBasisMatrix(self.T(x, el)) # interpolation matrix
            I = self.gll.gllweights.dot( Mx.dot( self.uc[K] ) )
#             Iu = Mx.dot( self.uc[K] )
#             I = self.gll.gllweights.dot( Iu )

        return ( dx/2 * I )

    def checkXticsCross(self, x):
        qval = self.mesh.nodes[1:] - self.mesh.nodes[:-1]
        val = (x[1:] - x[:-1]) * qval
        mbool = (val < 0)

        return any(mbool)
#     ============================================


class ExpKrylov (ExpSL):


    def __init__(self, mesh, gll, xb, order = 4):

        super().__init__(mesh, gll, xb, order = 4)

        from scipy.linalg import expm
        from scipy.sparse import spdiags
        from scipy.linalg import norm
        self.expm = expm
        self.spdiags = spdiags
        self.norm = norm
        self.nv = mesh.getnNodes()
        self.x = mesh.nodes
        self.CC = None

    def evaluate(self, u0, uc, dt):
        # Computes exp(tC)u0 with C a matrix and u0 a vector using Krylov suspace projection

        self.CC = -0.5 * self.C.dot( self.spdiags(uc, 0, self.nv, self.nv, format='csr') )

        ee1 = self.np.zeros((self.nv,)); ee1[0] = 1
        V = self.np.zeros((self.nv,self.nv))

        normu0 = self.norm(u0)
        V[:,0] = u0 / normu0
        m = int(self.nv/4)
        H = self.np.zeros((m+1,m))
        renorm  = self.np.zeros((m,))
        exp = None

        for j in range(m):
            vj = self.CC.toarray().dot(V[:,j])
            tvj = self.np.zeros((self.nv,))

            for i in range(j+1):
                H[i][j] = V[:,i].dot(vj)
                tvj += H[i][j] * V[:,i]

            wj = vj - tvj
            H[j+1][j] = self.norm(wj)

            if H[j+1][j] <= 1e-20:
                print('krylov iterate')
                break

            j1 = j+1
            EtH = self.expm(dt*H[:j1,:j1])
            EtHe1 = EtH.dot(ee1[:j1])
            exp =  V[:,:j1].dot(normu0*EtHe1)
            renorm[j] = H[j+1][j] * abs(EtHe1[j]) * normu0

            if renorm[j] <= 1e-20: break

            if j < m: V[:,j+1] = wj / H[j+1][j]

        return exp
#     ============================================

# class ExpSLDG (ExpSL):
# 
# #     import time
# 
#     def __init__(self, mesh, gll, xb, init=None, flowtype=None, order = 4):
# 
#         super().__init__(mesh, gll, xb, order)
# 
#         if (flowtype is not None):
#             self.constantFlow = flowtype['constant']
#             self.nonlinearFlow = flowtype['nonlinear']
# 
#         self.init = init
#         self.epsilon = 1e-15
#         self.loc = self.np.zeros((self.nv,),dtype=int) # arrival locations
#         self.GT = 0.
#         el = 0
#         for K in mesh.elements: self.loc[K] = el; el += 1
#    
# 
#     def boundaryfix(self, x):
#         xl = self.xb[0]; xr = self.xb[1]
#         if isinstance(x,float):
#             if self.mesh.anyPeriodic:
#                 Lx = xr - xl
#                 if x < xl : x = xr - self.np.mod(xl-x, Lx)
#                 if x > xr : x = xl + self.np.mod(x-xr, Lx)
#             else:
#                 if x < xl: x = xl
#                 if x > xr: x = xr
#             return x
#         if (self.mesh.anyPeriodic):
#             Lx = xr - xl; left = (x < xl); right = (x > xr)
#             if any(left): x[left] = xr - self.np.mod(xl-x[left], Lx)
#             if any(right): x[right] = xl + self.np.mod(x[right]-xr, Lx)
#         else: 
#             left = (x < xl); right = (x > xr)
#             if any(left): x[left] = xl
#             if any(right): x[right] = xr
#         return x
# 
#     def Integrate(self, x): 
# 
# #         t = self.time.time()
#         nbool = abs(x[1:]-x[:-1]) < self.epsilon
#         x[1:][nbool] += self.epsilon
#         x[:-1][nbool] -= self.epsilon
# 
#         dx = self.x - x
#         left = (dx > 0)
#         tinyDx = abs(dx) < self.epsilon
# 
#         self.foot_loc = -self.np.ones((self.nv,), dtype=int)
#         Ne = self.mesh.getNe()
# #         for el in range(Ne):
# #             mbool = self.mesh.inclusion(x, el)
# #             self.foot_loc[mbool] = el
#         sbool = ( x >= self.xb[0] ) & ( x <= self.xb[1] )
#         self.foot_loc[sbool] = self.mesh.contains(x[sbool])
# #         print(self.foot_loc)
# 
#         Iu = self.np.zeros((self.nv,))
#         if self.mesh.anyPeriodic:
#             el = 0
#             L = self.xb[1] - self.xb[0]
#             self.GT = self.Partition(self.xb[0], self.xb[1], 0, Ne-1)
#             period = self.np.zeros((self.nv,))
#             sbool = ~sbool
#             if sbool.any():
#                 xl = self.boundaryfix(x[sbool])
#                 lbool = (sbool & left)
#                 if lbool.any(): period[lbool] = (xl[left[sbool]]-x[lbool]) // L
#                 lbool = (sbool & ~left)
#                 if lbool.any(): period[lbool] = (x[lbool]-xl[~left[sbool]]) // L
#                 self.foot_loc[sbool] = self.mesh.contains(xl)
#                 x[sbool] = xl
#             for l in range(self.nv):
#                 if tinyDx[l]:  continue
#                 Iu[l] =  self.int_ab(x[l], l, left[l], period[l])
# #             for K in self.mesh.elements: # suitable for parallel implementation
# #                 Iu[K] =  self.int_ab(x[K], self.mesh.nodes[K], self.foot_loc[K], el, left[K], period[K])
# #                 el += 1
#         else:
#             self.foot_loc[self.foot_loc < 0 & left] = self.mesh.Bedges[0][-1]
#             self.foot_loc[self.foot_loc < 0] = self.mesh.Bedges[1][-1]
#             el = 0
#             for l in range(self.nv):
#                 Iu[l] =  self.int_ab(x[l], l, left[l])
# #             for K in self.mesh.elements: # suitable for parall implementation
# #                 Iu[K] =  self.int_ab(x[K], self.mesh.nodes[K], self.foot_loc[K], el, left[K])
# #                 el += 1
# #         print("elapsed time:", self.time.time() - t)
# 
#         return Iu
# 
# 
#     def int_ab(self, xt, l, left=True, period = 0):
# 
#         x = self.x[l]
#         foot = self.foot_loc[l]
#         el = self.loc[l]
#         I = self.Partition(xt, x, foot, el) if (left)\
#                 else self.Partition(x, xt, el, foot)
# 
#         I += period * self.GT
# 
#         return (I if left else -I)
# 
# #     def int_ab(self, xt, x, foot, el, dx, left, GT = 0., period = None):
# #         # suitable for parallel implementation
# # 
# #         if isinstance(xt, float):
# #             return self.Partition(xt, x, foot, el, dx)
# # 
# #         sgn = self.np.ones(len(xt)); sgn[~left] = -1.
# #         lbool = xt < x
# # 
# #         if (period is None):
# #             I = self.np.zeros((len(xt),))
# #             for l in range(len(xt)):
# #                 if abs(dx[l]) < self.epsilon: continue
# #                 I[l] = self.Partition(xt[l], x[l], foot[l], el, lbool[l])
# #             return (sgn * I)
# # 
# #         I = period * GT
# #         sbool = (period > 0) & ~left
# #         for l in range(len(xt)):
# #             if abs(dx[l]) < self.epsilon: continue
# #             if sbool[l]:
# #                 I[l] -= self.Partition(xt[l], x[l], foot[l], el, lbool[l])
# #             else:
# #                 I[l] += self.Partition(xt[l], x[l], foot[l], el, lbool[l])
# # 
# #         return (sgn * I)
# 
#     def Partition(self, xt, x, foot, el):
# 
#         I = 0.; ordered = (xt < x)
#         if not(ordered): # swap
#             tmp = xt; xt = x; x = tmp
#             temp = foot; foot = el; el = temp
#         xl = xt
# # #         if (foot==el): return sgn * self.GLLsum(xl, x, el)
# # #         tag = False
# #         for er in range(foot,el+1):
# #             xr = self.mesh.nodes[self.mesh.edges[self.mesh.edges2[er][1]]]
# #             if (xr < x):
# #                 I += self.GLLsum(xl, xr, er)#, tag)
# #             else:
# #                 I += self.GLLsum(xl, x, er)
# #             xl = xr; #tag = True
# #         if (foot==el):
# #             I = self.GLLsum(xl, x, foot)
# #         else:
# #             for er in range(foot,el+1):
# #                 xr = self.mesh.nodes[self.mesh.edges[self.mesh.edges2[er][1]]]
# #                 if (xr < x):
# #                     I += self.GLLsum(xl, xr, er)
# #                 else:
# #                     I += self.GLLsum(xl, x, er)
# #                 xl = xr
# 
#         for er in range(foot,el+1):
#             xr = self.mesh.nodes[self.mesh.edges[self.mesh.edges2[er,1]]]
#             I += self.GLLsum(xl, min(xr, x), er)
#             xl = xr
# 
#         return (I if ordered else -I)
# 
#     def GLLsum(self, xl, xr, el):#, completeElement=False):
# 
#         dx = xr - xl
# 
#         if (abs(dx) < self.epsilon): return 0.
# 
#         K = self.mesh.elements[el]
#         I = 0.
#         
#         if self.mesh.dG0:
#             I = 2 * self.uc[K][0]
#         elif K.size < 3:
#             x0 = self.mesh.nodes[K[0]]
#             I = 2*self.uc[K][0] + (self.uc[K][1]-self.uc[K][0])\
#                         * (xl + xr - 2*x0) / self.mesh.size(el)
#         else:
#             x = xl + (dx/2)*(self.gll.gllnodes + 1)
#             lx = self.gll.LagrangeBasisMatrix(self.T(x, el))
#             Iu = lx.dot( self.uc[K] )
#             I = self.gll.gllweights.dot( Iu )
#         
# #         if (completeElement or self.mesh.dG0):
# #             I = self.gll.gllweights.dot( self.uc[K] )
# #         else:
# #             N = K.size
# #             x = xl + (dx/2)*(self.gll.gllnodes + 1)
# #             lx = self.gll.LagrangeBasisMatrix(self.T(x, el))
# #             Iu = lx.dot( self.uc[K] ) #* dxN
# #             I = self.gll.gllweights.dot( Iu )
# # #             lxj = 0; Iun = 0; xn = 0
# # #             for n in range(N):
# # #                 xn = xl + (dx/2) * (self.gll.gllnodes[n]+1)
# # #                 Iun = 0.
# # #                 for j in range(N):
# # #                     nodes = self.mesh.nodes[K]
# # #                     nodes -= xn
# # #                     nodes[j] = 1
# # #                     lxj = (self.np.prod(nodes) / self.gll.ldeno[j] )
# # #                     Iun += lxj * self.uc[K[j]]
# # #                 I += self.gll.gllweights[n] * Iun
# 
#         return ( dx/2 * I )
# 
#     def checkXticsCross(self, x):
#         qval = self.mesh.nodes[1:] - self.mesh.nodes[:-1]
#         val = (x[1:] - x[:-1]) * qval
#         mbool = (val < 0)
# 
#         return any(mbool)
# #     ============================================
# 
