class polynomialEx:

    # class variables

    # Constructor
    def __init__(self,Re): # with numpy module
        self.fluxtype = {'convex':False,'concave':False,'linear':True}
        x = None
        t = None
        u0 = None
        self.Re = Re
        import numpy as np
        self.np = np
        self.epsilon = 2.220446049250313e-16
    
    def set(self, x, t, u0 = 0.):
        self.x = x
        self.t = t
        if isinstance(u0,float) and not(isinstance(x,float)): 
            self.u0 = u0*self.ones()
        else: self.u0 = u0

    def zeros(self):
        return ( self.np.zeros(self.np.shape(self.x)) )

    def ones(self):
        return ( self.np.ones(self.np.shape(self.x)) )

    def u(self):  # unknown variable u = u(x,t)
        return ( self.x * (self.x**2 - 3*self.t**2) )

    def a(self): # advection velocity
        return ( self.ones() ) # ones

    def u_x(self): # derivative of u wrt x 
        return ( 3 * (self.x**2 - self.t**2) )

    def u_t(self): # derivative of u wrt t
        return ( -6 * self.x*self.t  )

    def u_xx(self): # second derivative of u wrt x
        return ( 6 * self.x  )

    def f(self):  # load term for convection-diffusion
        return ( self.u_t() + self.u_x() - self.Re * self.u_xx() )
  
    def flux(self,alpha = 0.0): # conservative flux function
        return ( (self.a() - alpha) * self.u0 ) 
    
    def dflux(self,alpha = 0.0):
        return ( (self.a() - alpha) ) # derivative of flux wrt u

    def optima(self, alpha = 0.0): # optimal values of f(u) wrt u
        return []

#    =================================
class LinAdv0001 (polynomialEx): #  u_t + u_x = 0
    # u0 = \sin(\pi x), x \in [-1,1]

    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)

    def u(self):  # unknown variable u = u(x,t)
        L = 2.; a = -1.
        x = self.np.mod(self.x-self.a()*self.t - a, L) - a
        return ( self.np.sin(self.np.pi*(x)) )

    def a(self): # advection velocity
        return ( self.ones() ) # ones

    def u_x(self): # derivative of u wrt x 
        return  ( self.np.pi * self.np.cos(self.np.pi*(self.x-self.t)) )

    def u_t(self): # derivative of u wrt t:
        return  (  -self.np.pi * self.np.cos(self.np.pi*(self.x-self.t)) )

#    =================================

class LinAdv0001b (LinAdv0001): #  u_t + u_x = 0
    # u0 = \sin(x), x \in [0,2\pi]

    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)

    def u(self):  # unknown variable u = u(x,t)
        L = 2. * self.np.pi; a = 0.
        x = self.np.mod(self.x-self.a()*self.t - a, L) - a
        return ( self.np.sin(x) )

    def a(self): # advection velocity
        return ( self.ones() ) # ones

    def u_x(self): # derivative of u wrt x 
        return  (  self.np.cos(self.x-self.t) )

    def u_t(self): # derivative of u wrt t:
        return  (  - self.np.cos(self.x-self.t) )

#    =================================

class LinAdv0002 (polynomialEx):  #  u_t + (a(x)u)_x = 0

    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)

    def u(self):  # unknown variable u = u(x,t)
        x = self.x / 2; t = self.t
        if abs(t) < self.epsilon:
            return self.ones()
        return self.np.exp(-t) / ( self.np.cos(x)*self.np.cos(x) + self.np.exp(-2*t)*self.np.sin(x)*self.np.sin(x))

    def a(self): # advection velocity
        return ( self.np.sin(self.x) )

    def u_x(self): # derivative of u wrt x 
        return  ( 0.0 )

    def u_t(self): # derivative of u wrt t:
        return  ( 0.0 )

#    =================================

class LinAdv0003 (polynomialEx):  #  u_t + u_x = 0  (linear advection)

    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)
        self.fluxtype = self.fluxtype.fromkeys(self.fluxtype,False)
        self.fluxtype['convex'] = True

    def u(self):  # initial data, u0 = square wave function
        # Cockburn & Shu 2001
        u = self.zeros()
        x = self.np.mod(self.x - self.t, 2*self.np.pi)
        sbool = self.x - self.t < 0
        x[sbool] -= 2*self.np.pi
        if not(isinstance(x,float)):
            sbool = abs(x[1:]-x[:-1]) < 10*self.epsilon
            x[1:][sbool] += 10*self.epsilon
            x[:-1][sbool] -= 10*self.epsilon
        u[ abs(x-self.np.pi) < self.np.pi/2 ] = 1.
        return u

    def a(self): # advection velocity
        return ( self.ones() )

    def u_x(self): # derivative of u wrt x 
        return  ( 0.0 )

    def u_t(self): # derivative of u wrt t:
        return  ( 0.0 )

#    =================================

class LinAdv0004 (polynomialEx): #  u_t + (a(t)u)_x = 0

    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)

    def u(self):  # unknown variable u = u(x,t)

        return ( (3 + self.np.sin(self.np.pi*(self.x+1+self.np.cos(self.t))))/4  )

    def a(self): # advection velocity
        return ( self.np.sin(self.t) * self.ones() )

    def u_x(self): # derivative of u wrt x 
        return  ( 0.0 )

    def u_t(self): # derivative of u wrt t:
        return  ( 0.0 )

#    =================================

class LinAdv0005 (polynomialEx): #  u_t + (a(t)u)_x = 0

    # Constructor: Shu's linear test. JCP 126 (1996) 
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)
        self.l = -1.

    def u(self):  # unknown variable u = u(x,t)

        def G(x, bt, x0): return self.np.exp(-bt*(x-x0)**2)
        def F(x, alp, x0): 
            if isinstance(x, float): 
                return self.np.sqrt(max(1-alp**2*(x-x0)**2,0))
            u = 1-alp**2*(x-x0)**2
            u[u<0] = 0.
            return self.np.sqrt(u)
        x = self.np.mod(self.x-self.l - self.a()*self.t, 2.) + self.l
        if not(isinstance(x,float)):
            sbool = abs(x[1:]-x[:-1]) < 10*self.epsilon
            x[1:][sbool] += 10*self.epsilon
            x[:-1][sbool] -= 10*self.epsilon
        c = self.np.array([0.3,0.7,1.1,1.5]) + self.l; r = 0.1
        dlt = 0.005; alp = 10.; bt = self.np.log(2) / (36*dlt**2)
        if isinstance(x,float):
            if ( abs(x-c[0]) <= r ): return (G(x, bt, c[0]-dlt) + G(x, bt, c[0]+dlt) + 4*G(x, bt, c[0])) / 6
            elif ( abs(x-c[1]) <= r ): return 1.
            elif ( abs(x-c[2]) <= r ): return 1 - abs(10*(x-c[2]))
            elif ( abs(x-c[3]) <= r ): return (F(x, alp, c[3]-dlt) + F(x, alp, c[3]+dlt) + 4*F(x, alp, c[3])) / 6
            else: return 0.
        u = self.np.zeros(self.np.shape(x))
        sbool = ( abs(x-c[0]) <= r )
        u[sbool] = (G(x[sbool], bt, c[0]-dlt) + G(x[sbool], bt, c[0]+dlt) + 4*G(x[sbool], bt, c[0])) / 6
        u[ abs(x-c[1]) <= r ] = 1.
        sbool = ( abs(x-c[2]) <= r )
        u[sbool] = 1 - abs(10*(x[sbool]-c[2]))
        sbool = ( abs(x-c[3]) <= r )
        u[sbool] = (F(x[sbool], alp, c[3]-dlt) + F(x[sbool], alp, c[3]+dlt) + 4*F(x[sbool], alp, c[3])) / 6
        return u

    def a(self): # advection velocity
        return ( self.ones() ) # ones

class LinAdv0005b (LinAdv0005): #  u_t + u_x = 0

    # Constructor: Shu's linear test. JCP 126 (1996) 
    # See f.ex. Huang et al.2012, periodic bc
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)
        self.l = 0.
#    =================================

class NLinAdv0005 (polynomialEx): #  u_t + (u*u/2)_x = 0 (invicid Burgers)

    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)
        self.fluxtype = self.fluxtype.fromkeys(self.fluxtype,False)
        self.fluxtype['convex'] = True

    def u(self):  # Cockburn et al. 2001, x\in[0,1], t\in[0,0.4]
        return ( 0.25 + 0.5*self.np.sin(self.np.pi*(2*self.x-1)) )

    def a(self): # advection velocity
        return ( self.u0/2 )

    def u_x(self): # derivative of u wrt x 
        return  ( 0.0 )

    def u_t(self): # derivative of u wrt t:
        return  ( 0.0 )
    
    def dflux(self,alpha = 0.0):
        return ( (self.u0 - alpha) ) # derivative of flux wrt u


    def optima(self, alpha = 0.0): # optimal values of f(u) wrt u
        return [[alpha, alpha**2/2]] 

#    =================================

class NLinAdv0005b (NLinAdv0005): #  u_t + (u*u/2)_x = 0 (invicid Burgers)

    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)

    def u(self):  # initial data by Huang and Arbogast 2012
        # x \in [0,2]; t\in[0,T]; T = 0.25 or T = 3/(2\pi)
        return ( 0.5 + self.np.sin(self.np.pi*self.x) )

#    =================================

class NLinAdv0005c (NLinAdv0005): #  u_t + (u*u/2)_x = 0 (invicid Burgers)

    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)

    def u(self):  # Kometa et al. 2020
        return ( self.np.sin(self.np.pi*self.x) ) 

#    =================================

class NLinAdv0005d (NLinAdv0005): #  u_t + (u*u/2)_x = 0 (invicid Burgers)

    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)

    def u(self):  # Gottlieb et al. 2009, x\in[0,2], t\in[0,2]
        return ( 0.5 - 0.25*self.np.sin(self.np.pi*self.x) )

#    =================================

class NLinAdv0006 (polynomialEx): # u_t + (f(u))_x = 0 (Buckley-Leverett)
    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)
        self.fluxtype = self.fluxtype.fromkeys(self.fluxtype,False)

    def u(self):  #   u0, Huang & Arbogast 2012 x\in[0,1], t\in[0,0.4]
        x = self.x
#         if not(isinstance(x,float)):
#             sbool = abs(x[1:]-x[:-1]) < 10*self.epsilon
#             x[1:][sbool] += 10*self.epsilon
#             x[:-1][sbool] -= 10*self.epsilon
        if isinstance(x,float):
            if ( abs(x-0.025) <= 0.025 ): return 1 - 20*x
            elif ( abs(x-0.325) <= 0.075 ): return 0.5
            else: return 0.
        u = self.np.zeros(self.np.shape(x))
        sbool = ( abs(x-0.025) <= 0.025 )
        u[sbool] = 1 - 20*x[sbool]
        u[ abs(x-0.325) <= 0.075 ] = 0.5
        return u

    def a(self): # advection velocity
        return ( self.u0/(self.u0**2 + (1-self.u0)**2) )

    def u_x(self): # derivative of u wrt x 
        return  ( 0.0 )

    def u_t(self): # derivative of u wrt t:
        return  ( 0.0 )
    
    def dflux(self,alpha = 0.0):
        return ( (2*self.u0*(1-self.u0))/(self.u0**2 + (1-self.u0)**2)**2 - alpha) # derivative of flux wrt u

    def optima(self, alpha = 0.0): # optimal values of f(u) wrt u
        if abs(alpha + .125) < .125 - self.epsilon: #(alpha >= -.25 and alpha < 0.0):
            val = self.np.sqrt(4*alpha+1)
            n = [[-1.,-1.],[-1.,1.],[1.,-1.],[1.,1.]]
            for k in range(4):
                n[k][0] = (1 + n[k][0] * self.np.sqrt(-1 + n[k][1]*val/alpha - 1/alpha) ) / 2
                n[k][1] = n[k][0]**2 / (n[k][0]**2 + (1-n[k][0])**2) - alpha * n[k][0]
            return n
        elif abs(alpha - 1) < 1 - self.epsilon: #(alpha > 0 and alpha <= 2):
            val = self.np.sqrt(4*alpha+1)
            n = [[-1.,0.],[1.,0.]]
            for k in range(2):
                n[k][0] = (1 + n[k][0] * self.np.sqrt(-1 + val/alpha - 1/alpha) ) / 2
                n[k][1] = n[k][0]**2 / (n[k][0]**2 + (1-n[k][0])**2) - alpha * n[k][0]
            return n
        return [[0.,0.],[1.,1.]] if (alpha == 0.0) else []

#    =================================

class NLinAdv0006b (NLinAdv0006): # u_t + (f(u))_x = 0 (non-convex flux)
    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)

    def a(self): # advection velocity
        #   f(u0)/u0, Huang & Arbogast 2012 x\in[0,2], t\in[0,0.4]
        if isinstance(self.u0,float): 
            return 10. if (self.u0 <= 0.01) else (1 / self.np.sqrt(self.u0))
        vec = 10. * self.zeros()
        mbool = self.u0 > 0.01
        vec[mbool] = 1 / self.np.sqrt(self.u0[mbool])
        return vec

    def u_x(self): # derivative of u wrt x 
        return  ( 0.0 )

    def u_t(self): # derivative of u wrt t:
        return  ( 0.0 )
    
    def dflux(self,alpha = 0.0):
        if isinstance(self.u0,float): 
            vec = 10. if (self.u0 <= 0.01) else (0.5 / self.np.sqrt(self.u0))
            return (vec - alpha)
        vec = 10. * self.zeros()
        mbool = self.u0 > 0.01
        vec[mbool] = 0.5 / self.np.sqrt(self.u0[mbool])
        return (vec  - alpha) # derivative of flux wrt u

#     def optima(self): # optimal values of f(u) wrt u
#         return []
    def optima(self, alpha = 0.0): # optimal values of f(u) wrt u
        if (abs(alpha) < 5. - self.epsilon):
            if abs(alpha) < self.epsilon: return [[0.01,0.1]]
            val = 0.25 / (alpha**2)
            fval = self.np.sqrt(val) - alpha * val
            return [[val, fval]]
        elif (abs(alpha - 10.) < self.epsilon):
            return [[0.01,0.1]]
        return [] # breakpoints inclusive

#    =================================

class NLinAdv0007 (NLinAdv0005): #  u_t + (u*u/2)_x = 0 (invicid Burgers)

    # Constructor
    def __init__(self, Re = 0.): # with numpy module
        super().__init__(Re)
        self.fluxtype = self.fluxtype.fromkeys(self.fluxtype,False)
        self.fluxtype['convex'] = True

    def u(self):  # initial data for KdV by Pietra and Pohl 1999
        return ( 1.5 + self.np.cos(self.x) )

#    =================================
