class timeStepping:


    def __init__(self, data, gll, mesh, bc, nsteps):

        from Assembly import Assembly
        import Initdata

        self.semAss = Assembly(mesh, gll, data.Re)

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

        self.ub = None; self.uc = None
        self.sol = []
        self.u0 = self.init.compute4ex(self.t0)

        self.exactsol = False
        self.storeAll = data.storeAll
        if (data.Dirbctype['constant'] or data.Dirbctype['piecewise']):
            self.bc.set(self.mesh.Dnodes, self.t0)
            self.semAss.setbc(self.bc.DirichletBc())
        if (self.mesh.anyDirbc):
            self.bc.set(mesh.Dnodes, self.t0)
            self.semAss.setbc(bc.DirichletBc())

    def validVector(self, u):
        return  (not(any(self.np.isnan(u))) or self.np.linalg.norm(u)>self.inf)

        
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
            du = u0[el] - ue[el]
            dx  = self.mesh.size(l) / 2; l += 1
            err += self.gll.gllweights.dot(du**2) * dx
        return self.np.sqrt(err)
        
    def errorLinf(self, u0, t, ue = None):
        if (ue is None): ue = self.init.compute4ex(t)
        return self.np.linalg.norm(abs(u0 - ue), self.np.inf)
#     ============================================

class SLRKDG (timeStepping): # derived class

    def __init__(self, data, gll, mesh, bc, nsteps):

        super().__init__(data, gll, mesh, bc, nsteps)

        from Assembly import AssemblyDG

        self.semAss = AssemblyDG(mesh, gll, data.Re)
        self.semAss.initd = self.init
        self.semAss.flowtype = data.flowtype
        self.semAss.order = data.orderRK
        self.semAss.src = bc

#         from timeStepping1 import ExpSLDG

#         self.expsl = ExpSLDG(mesh, gll, data.xb, self.semAss.initd, data.flowtype)# data.orderRK)
        self.expsl = ExpSL(mesh, gll, data.xb, self.semAss.initd, data.flowtype)# data.orderRK)
        self.semAss.tol = 1e-15
        self.semAss.limit = data.DGlimiter['PP'] or \
                data.DGlimiter['MPP']
        if self.semAss.limit:
            self.semAss.ltype = data.DGlimiter
            self.semAss.m0 = data.DGlimiterMinMax[0]
            self.semAss.M0 = data.DGlimiterMinMax[1]


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

        if self.mesh.dG0: u0 = self.semAss.limiter(u0)
        if self.storeAll: self.sol.append(self.u0)

        print('-'*60)
#         for r in range(1):
        while (len(steps)>0):

            dt = steps.pop()
            if self.stepBystepReport:
                print(" TIME STEP %i OF %i;\tdt = %f;\tt = %f" %(nstep, self.nsteps, dt, self.t+dt))

            uc = self.init.compute4a(tn, u0)
            self.semAss.a = uc 
            #===== Attempt 1: initial idea =====

            f = self.expsl.evaluate(u0, uc, dt, tn+dt, False)
            g = self.semAss.Derivative(f, tn)
            l = self.mesh.edges[self.mesh.edges2]
#             print(f[l])
            g[l[:,0]] += f[l[:,0]] # left edge
            g[l[:,1]] -= f[l[:,1]] # right edge

            self.semAss.MassInv(g, self.mesh.anyDirbc)


            L2Stab = self.errorL2(u0, tn, 0.)

            u0 = u0 + g
#             u0 = self.semAss.mpplim(u0)

#             u0 = self.semAss.limiter(u0)
            if self.semAss.flowtype['nonlinear'] or self.mesh.dG0: u0 = self.semAss.limiter(u0)
            
            self.t = tn + dt
            tn  = self.t

            L2Stab -= self.errorL2(u0, tn, 0.)
            if L2Stab < 0.:
                print('L2-stability = %8.4E failed!' % L2Stab)

            converged = (converged and self.validVector(u0))

            if converged:
                if self.storeAll: 
                    self.sol.append(u0)
                if self.exactsol: 
                    errL2 = self.errorL2(u0, self.t)
                    print("\tNumerical ERROR (L2) at time (t = %f): %f" %(self.t, errL2))
            else: 
                print("invalid solution vector!"); self.sol.append(u0)
                print(u0)
                if self.storeAll: 
                    while len(steps)>0: dt = steps.pop(); self.sol.append(u0)
                return None

            self.nstepsdone = nstep
            nstep += 1
        if not(self.storeAll):self.sol = u0
        print('-'*60)

        return None
#     ============================================

class ExpSL:


#     def __init__(self, mesh, gll, xb, order = 4):
    def __init__(self, mesh, gll, xb, init, flowtype = None, order = 4):
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
        self.constantFlow  = False
        self.nonlinearFlow  = True
        self.init = None
        self.nsteps = 1
        self.epsilon = 1e-12

#         from scipy.integrate import solve_ivp
#         self.solve = solve_ivp

        #ExpSLDG
        self.init = init
#         self.epsilon = 1e-15
        self.loc = self.np.zeros((self.nv,),dtype=int) # arrival locations
        self.GT = 0.
        el = 0
        for K in mesh.elements: self.loc[K] = el; el += 1

        if (flowtype is not None):
            self.constantFlow = flowtype['constant']
            self.nonlinearFlow = flowtype['nonlinear']

    def evaluate(self, u0, uc, dt, t = None, nonExact=True):

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
        elif (self.order==1): # Euler
            while (len(steps)>0):
                ds = steps.pop()
                F = self.InterpSL(xi) if (nonlinear) else self.init.compute4A(xi, t-ds)
                dx = ds*F
                xi = xi-dx;
                nsteps += 1; t -= ds
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
        else:  # RK4
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
        self.uc = u0
#         xi = self.x + self.np.cos(t) - 1 # test test test
        if self.checkXticsCross(xi): print("Warning: Characteristic curves cross each other!")
        return ( self.InterpSL(xi) if (nonExact) else self.Integrate(xi) )


    def InterpSL(self, x, u = None):
        self.boundaryfix(x)
        Iu = self.np.zeros((self.nv,))
        if (u is not None): self.uc = u

        sbool = abs(x[1:]-x[:-1]) < self.epsilon
        x[1:][sbool] += 10*self.epsilon
        x[:-1][sbool] -= 10*self.epsilon

        l = 0
        self.foot_loc = self.mesh.contains(x)
        for el in self.mesh.elementsL:
            here = (self.foot_loc == l)
#             here = self.mesh.inclusion(x, l)
            if here.any():
                lx = self.gll.LagrangeBasisMatrix(self.T(x[here],l))
                Iu[here] = lx.dot(self.uc[el])
            l += 1
        return Iu

    def T(self,x,el): # affine transformation into [-1,1]
        l = self.mesh.edges[self.mesh.edges2[el][0]]
        return (x-self.mesh.nodes[l])*(2/self.mesh.size(el))-1.

    def boundaryfix(self, x):
        xl = self.xb[0]; xr = self.xb[1]
        if isinstance(x,float):
            if self.mesh.anyPeriodic:
                Lx = xr - xl
                if x < xl or x > xr: x = self.np.mod(x-xl, Lx) + xl
            else:
                if x < xl: x = xl
                if x > xr: x = xr
            return x
        if (self.mesh.anyPeriodic):
            Lx = xr - xl
            sbool = (x < xl) | (x > xr)
            x[sbool] = self.np.mod(x[sbool]-xl, Lx) + xl
        else: 
            lbool = (x < xl); rbool = (x > xr)
            if any(lbool): x[lbool] = xl
            if any(rbool): x[rbool] = xr
        return x
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

    def Integrate(self, x): 

#         t = self.time.time()
        nbool = abs(x[1:]-x[:-1]) < self.epsilon
        x[1:][nbool] += self.epsilon
        x[:-1][nbool] -= self.epsilon

        dx = self.x - x
        left = (dx > 0)
        tinyDx = abs(dx) < self.epsilon

        self.foot_loc = -self.np.ones((self.nv,), dtype=int)
        Ne = self.mesh.getNe()
#         for el in range(Ne):
#             mbool = self.mesh.inclusion(x, el)
#             self.foot_loc[mbool] = el
        sbool = ( x >= self.xb[0] ) & ( x <= self.xb[1] )
        self.foot_loc[sbool] = self.mesh.contains(x[sbool])
#         print(self.foot_loc)

        Iu = self.np.zeros((self.nv,))
        if self.mesh.anyPeriodic:
            el = 0
            L = self.xb[1] - self.xb[0]
            self.GT = self.Partition(self.xb[0], self.xb[1], 0, Ne-1)
            period = self.np.zeros((self.nv,))
            sbool = ~sbool
            if sbool.any():
                xl = self.boundaryfix(x[sbool])
                lbool = (sbool & left)
                if lbool.any(): period[lbool] = (xl[left[sbool]]-x[lbool]) // L
                lbool = (sbool & ~left)
                if lbool.any(): period[lbool] = (x[lbool]-xl[~left[sbool]]) // L
                self.foot_loc[sbool] = self.mesh.contains(xl)
                x[sbool] = xl
            for l in range(self.nv):
                if tinyDx[l]:  continue
                Iu[l] =  self.int_ab(x[l], l, left[l], period[l])
#             for K in self.mesh.elements: # suitable for parallel implementation
#                 Iu[K] =  self.int_ab(x[K], self.mesh.nodes[K], self.foot_loc[K], el, left[K], period[K])
#                 el += 1
        else:
            self.foot_loc[self.foot_loc < 0 & left] = self.mesh.Bedges[0][-1]
            self.foot_loc[self.foot_loc < 0] = self.mesh.Bedges[1][-1]
            el = 0
            for l in range(self.nv):
                Iu[l] =  self.int_ab(x[l], l, left[l])
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

        return (I if left else -I)


    def Partition(self, xt, x, foot, el):

        I = 0.; ordered = (xt < x)
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
        
        if K.size < 3:
            x0 = self.mesh.nodes[K[0]]
            I = 2*self.uc[K][0] + (self.uc[K][1]-self.uc[K][0])\
                        * (xl + xr - 2*x0) / self.mesh.size(el)
        else:
            x = xl + (dx/2)*(self.gll.gllnodes + 1)
            lx = self.gll.LagrangeBasisMatrix(self.T(x, el))
            Iu = lx.dot( self.uc[K] )
            I = self.gll.gllweights.dot( Iu )
        
        return ( dx/2 * I )

    def checkXticsCross(self, x):
        qval = self.mesh.nodes[1:] - self.mesh.nodes[:-1]
        val = (x[1:] - x[:-1]) * qval
        mbool = (val < 0)

        return any(mbool)
#     ============================================
