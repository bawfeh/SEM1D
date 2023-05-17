class Initdata:

    # class variables

    # Constructor
    def __init__(self, BCnSource, mesh): 
        self.BCnSource = BCnSource
        self.mesh = mesh

    def compute4ex(self, t = 0.):

        self.BCnSource.set(self.mesh.nodes, t)

        if (self.mesh.dG0) : # re-average into piecewise constant data
            u = self.BCnSource.DirichletBc()
            e = self.mesh.edges[self.mesh.edges2]
            w = (u[e[:,0]] + u[e[:,1]])/2
            u[e[:,0]] = w
            u[e[:,1]] = w
            return u

        return self.BCnSource.DirichletBc()

    def compute4Ex(self, nodes, t = 0.):

        self.BCnSource.set(nodes, t)

        return self.BCnSource.DirichletBc()

    def compute4a(self, t = 0., u = 0.):
    
        self.BCnSource.set(self.mesh.nodes, t, u)

        return self.BCnSource.a()

    def compute4A(self, nodes, t, u = 0.):
    
        self.BCnSource.set(nodes, t, u)

        return self.BCnSource.a()

    def compute4f(self, t, u = 0.):
    
        self.BCnSource.set(self.mesh.nodes, t, u)

        return self.BCnSource.f()

    def compute4flux(self, u, onlyAtEdges, t = 0., a = None):
        
        if (onlyAtEdges):

            self.BCnSource.set(self.mesh.nodes[self.mesh.edges], t, u[self.mesh.edges])
            if (a is None): return self.BCnSource.flux()
            else: return self.BCnSource.flux(a[self.mesh.edges])
    
        self.BCnSource.set(self.mesh.nodes, t, u)

        if (a is None): return self.BCnSource.flux()
        return self.BCnSource.flux(a)
    
    def compute4dflux(self, u, onlyAtEdges, t = 0., a = None):
        
        if (onlyAtEdges):

            self.BCnSource.set(self.mesh.nodes[self.mesh.edges], t, u[self.mesh.edges])
            if (a is None): return self.BCnSource.dflux()
            else: return self.BCnSource.dflux(a[self.mesh.edges])
    
        self.BCnSource.set(self.mesh.nodes, t, u)

        if (a is None): return self.BCnSource.dflux()
        return self.BCnSource.dflux(a)

# ====================================================

class Source:

    # Constructor
    def __init__(self,example = 'polynomialEx', Re = 0.): 

        self.sol = None
        self.fluxtype = None
        import myExamples

        if self.compareStr(example,'polynomialEx'):
            self.sol = myExamples.polynomialEx(Re)

        elif self.compareStr(example,'LinAdv0001'):
            self.sol = myExamples.LinAdv0001(Re)

        elif self.compareStr(example,'LinAdv0001b'):
            self.sol = myExamples.LinAdv0001b(Re)

        elif self.compareStr(example,'LinAdv0002'):
            self.sol = myExamples.LinAdv0002(Re)

        elif self.compareStr(example,'LinAdv0003'):
            self.sol = myExamples.LinAdv0003(Re)

        elif self.compareStr(example,'LinAdv0004'):
            self.sol = myExamples.LinAdv0004(Re)

        elif self.compareStr(example,'LinAdv0005'):
            self.sol = myExamples.LinAdv0005(Re)

        elif self.compareStr(example,'LinAdv0005b'):
            self.sol = myExamples.LinAdv0005(Re)

        elif self.compareStr(example,'NLinAdv0005'):
            self.sol = myExamples.NLinAdv0005(Re)

        elif self.compareStr(example,'NLinAdv0005b'):
            self.sol = myExamples.NLinAdv0005b(Re)

        elif self.compareStr(example,'NLinAdv0005c'):
            self.sol = myExamples.NLinAdv0005c(Re)

        elif self.compareStr(example,'NLinAdv0005d'):
            self.sol = myExamples.NLinAdv0005d(Re)

        elif self.compareStr(example,'NLinAdv0006'):
            self.sol = myExamples.NLinAdv0006(Re)

        elif self.compareStr(example,'NLinAdv0006b'):
            self.sol = myExamples.NLinAdv0006b(Re)

        elif self.compareStr(example,'NLinAdv0007'):
            self.sol = myExamples.NLinAdv0007(Re)

        else:
            print("Test case example " + example + " is undefined!")

        if (self.sol is not None) : self.fluxtype = self.sol.fluxtype
          

    def compareStr(self,key1,key2):
        return (key1.upper()==key2.upper())

    def set(self, x, t = 0., u0 = 0.):
        self.sol.set(x, t, u0)

    def f(self):
        return self.sol.f()
    
    def a(self):
        return self.sol.a()

    def DirichletBc(self):
        return self.sol.u()    # u = uD

    def Gradux(self):
        return self.sol.u_x()  # Neuman BC: K.\Grad u . n = uN

    def flux(self, alpha = 0.):
        return self.sol.flux(alpha)  # Div(flux) = conservative convective term

    def dflux(self, alpha = 0.):   # Derivative of the flux
#         if (self.sol.fluxtype['linear']): 
#             return self.sol.a()
        return self.sol.dflux(alpha)

    def dfluxsign(self, u1, u2, x, t, alpha = 0.): # sign test for derivatives of C1-functions
        self.set(x, t, (u1+u2)/2)
        return ( self.dflux(alpha) > 0 )

    def optima(self, alpha = 0.0):
        return [] if (self.sol.fluxtype['linear']) else  self.sol.optima(alpha)
# ==========================================================

class Sources:

    # Constructor
    def __init__(self, source, mesh, data): 

        self.sources = mesh.np.zeros((mesh.getnNodes(),))
        self.var_source = False
        self.Re = data.Re
        self.mesh = mesh
        self.source = source

        if (data.sourcetype['piecewise']):
            if (self.validateSource()): 
                self.load(data.sourcetype, data.source_data, data.source_pos)
        elif (data.sourcetype['variable']):
            self.var_source = True
        elif (data.sourcetype['constant']):
            self.load(data.sourcetype, data.source_data)


    def load (self, sourcetype,sourcedata, sourcepos = []):
        
        if (sourcetype['constant']):
            if (len(sourcedata)>0): self.sources = sourcedata[0]
        elif (sourcetype['piecewise']):
            l = 0; N = self.mesh.getN()
            for val in sourcedata:
                elements = self.mesh.elements[sourcepos[l]]; l += 1
                for i in range(N): self.sources[elements[i]] = val
        return None

    def evaluate(self, t):
        self.source.set(mesh.nodes, t, 0.)
        self.sources = self.source.f()
        return None

    def validateSource(self):
        for el in self.sourcepos:
            if not(self.mesh.isValid(el)): 
                print("Element %i exceeds MAX ( %i )!"%(el,self.mesh.getNe()))
                return False
        return True

