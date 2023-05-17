class BoundaryConditions:
    ## Computes boundary data

    from numpy import zeros

    def __init__(self, BCnSource, mesh, data):
        
        self.numDnodes = len(mesh.Dnodes)
        self.numNnodes = len(mesh.Nnodes)
        self.dirbc = None

        if (mesh.anyDirbc and data.Dirbctype['variable']):
            self.dirbc = data.np.zeros((numDnodes,))
        if (mesh.anyNeumbc and data.Neumbctype['variable']):
            self.Neumbc = data.np.zeros((numNnodes,))

        self.zero = False
        if data.Dirbctype['homogeneous']: 
            self.zero = (mesh.anyNeumbc and data.Neumbctype['homogeneous'])
        if data.Neumbctype['homogeneous']: 
            self.zero = (mesh.anyDirbc and data.Dirbctype['homogeneous'])

        self.mesh = mesh
        self.data = data
        self.BCnSource = BCnSource

    def calculate(self, t):
        if (self.mesh.anyDirbc):
            self.BCnSource.set(self.mesh.nodes[self.mesh.Dnodes], t, 0., self.data.Re)
            self.dirbc =  self.BCnSource.u()
        if (self.mesh.anyNeumbc):
            self.BCnSource.set(self.mesh.nodes[self.mesh.Nnodes], t, 0., self.data.Re)
            self.Neumbc = -self.mesh.normal[[v for v in self.orientNeumbc.values()]] * self.BCnSource.u_x()
        return None

    def get(self, bc):
        if (self.mesh.anyDirbc):
            bc = self.zeros((self.mesh.getnNodes(),))
            if self.data.Dirbctype['constant']:
                bc[self.mesh.Dnodes] = self.data.Dirbc_data[0]
            elif self.data.Dirbctype['piecewise']:
                l = 0
                for val in self.data.Dirbc_data: 
                    bc[self.mesh.Dnodes[l]] = val; l += 1
            elif self.data.Dirbctype['variable']:
                l = 0
                for i in self.mesh.Dnodes: 
                    bc[i] = self.Dirbc[l]; l += 1
        if (self.mesh.anyNeumbc):
            bc = self.zeros((self.mesh.getnNodes(),))
            if self.data.Neumbctype['constant']:
                val = self.data.Neumbc_data[0]
                for i in self.mesh.Nnodes: bc[i] = val
            elif self.data.Neumbctype['piecewise']:
                l = 0
                for val in self.data.Neumbc_data: 
                    bc[self.mesh.Nnodes[l]] = val; l += 1
            elif self.data.Neumbctype['variable']:
                l = 0
                for i in self.mesh.Nnodes: 
                    bc[i] = self.Neumbc[l]; l += 1
        return None
