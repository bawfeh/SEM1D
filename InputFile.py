class InputFile:

    # class variables
    argv = []

    # Constructor
    def __init__(self, filename, printInfo = False,filenameOnly=False):
        self.filename = filename # Create an instance variable
        self.readfile()
        if printInfo : self.display_args(filenameOnly)

    # Instance methods

    def readfile(self):
        with open(self.filename) as fobj:
            argv = fobj.readline()
            alist = []
            while argv:
                alist = argv.split()
                if ('//' in alist): ind = alist.index('//'); del alist[ind:]
                if (not(argv.isspace()) and alist[0]!='#'): self.argv.append(alist)
                argv = fobj.readline()

    def display_args(self,filenameOnly=False):
        print("\nData read from your input file [%s]:\n" % self.filename)
        if not filenameOnly:
            for argv in self.argv:
                print(argv)

    def updateData(self,data):
        flag = True

        for argv in self.argv:
            flag1 = self.checkargs(argv,data)
            flag = flag and flag1

        if len(data.nx) != len(data.dx):
            print("The keys -Nx and -Dx must have the same number of arguments!")
            flag = False
        else:
            data.xb[-1] = data.xb[0] + \
                    sum([data.nx[j]*data.dx[j] for j in range(len(data.nx))])
            
        if data.sourcetype['piecewise']:
            if sum(data.nx) != len(data.source_data):
                print("Piece-wise source data does not match the number of elements!")
                flag = False
            if len(data.source_pos) != len(data.source_data):
                print("Piece-wise source data and position MUST be compatible!")
                flag = False
            
        if data.sourcetype['constant']:
            if len(data.source_data) > 1:
                print("Constant source data MUST be unique!")
                flag = False

        if data.Dirbctype['constant']:
            if len(data.Dirbc_data) != 1:
                print("Dirichlet data is inconsistent constant Dirichlet type")
                flag = False

        if data.Neumbctype['constant']:
            if len(data.Neumbc_data) != 1:
                print("Neumann data is inconsistent constant Neumman type")
                flag = False

        # needless to store every step if we do not need to export
#         data.storeAll = (data.storeAll and data.exportData)

        # needless to adapt stepping if dtmax=dtmin
        data.adaptiveStepping = (data.adaptiveStepping and (data.dtmax-data.dtmin)>1e-11)   
        if data.leftrightbc:
            data.orientNeumbc = data.orientNeumbc.fromkeys(data.orientNeumbc,False)
            data.orientDirbc = data.orientDirbc.fromkeys(data.orientDirbc,False)
        else:
            for key in data.orientNeumbc.keys():
                data.orientDirbc[key] = not(data.orientNeumbc[key])


        return flag

    def show_usage(self):
        print("\nUsage: ./<execfile> <option(s)>")
#         print("<execfile> is the name of the product file;")
#         print("Options:")
#         print("\t-h,--help\tGet help on how to run the program;")
#         print("\t-t t0 tN \tThe time interval [t0, tN];")
#         print("\t-dtmin val\t\tThe minimum time step to use;")
#         print("\t-dtmax val\t\tThe maximum time step to use;")
#         print("\t-save\t\tTo save the solution vector at each time step;")
#         print("\t-Re val\t\tThe Reynolds number;")
#         print("Use SINGLE SPACES to seperate a list of arguments/values!")

    def printErrorMsg(self,argv,str1='ONE',str2=''): 
        print(argv)
        print(str1)
        print("OptionArgumentError: %s option requires %s argument(s)! %s" % (argv[0], str1,str2));
        self.show_usage();

    def compareStr(self,key1,key2):
        return (key1.upper()==key2.upper())

    def checkargs(self,argv,data):

        if len(argv)==0:
            print("empy input vector!")
            return True

        argc = argv[0]
        if (argc.isspace() or len(argc)==0):
            print("Warning: Limit extra spaces between keywords/options!")
            return True

        if (self.compareStr(argc , '-h') or self.compareStr(argc ,'--help')):
            self.show_usage()

        elif self.compareStr(argc , 'origin'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.xb[0] = float(argv[1])

        elif self.compareStr(argc,'Nx'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.nx = [int(float(i)) for i in argv[1:]]
            data.Ne = sum(data.nx)

        elif self.compareStr(argc,'Dx'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.dx = [float(s) for s in argv[1:]]

        elif self.compareStr(argc,'-t'):
            if len(argv)<3:
                self.printErrorMsg(argv,'TWO')
                return False
            data.tb = [float(s) for s in argv[1:]]
            data.reset_dT()

        elif self.compareStr(argc,'-mpplimit'):
            if len(argv)<3:
                self.printErrorMsg(argv,'TWO')
                return False
            data.DGlimiterMinMax = [float(s) for s in argv[1:]]
            data.DGlimiter['PP'] = False 
            data.DGlimiter['MPP'] = True 

        elif self.compareStr(argc,'-pplimit'):
            data.DGlimiter['PP'] = True 
            data.DGlimiter['MPP'] = False 

        elif self.compareStr(argc,'-dtmax'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.dtmax = float(argv[1])
            if not(data.check_dtmax()):
                self.printErrorMsg(argv,'ONE', "dt < compuational time interval!")
                return False

        elif self.compareStr(argc,'-dtmin'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.dtmin = float(argv[1])
            if not(data.check_dtmin()):
                self.printErrorMsg(argv,'ONE', "dt < compuational time interval!")
                return False

        elif self.compareStr(argc,'-Re'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.Re = float(argv[1])

        elif self.compareStr(argc,'-Cr'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.Cr = float(argv[1])

        elif self.compareStr(argc,'-polyDeg'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.polyDeg = int(float(argv[1]))

        elif self.compareStr(argc,'-orderRK'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.orderRK = int(float(argv[1]))
            data.cfmethodName = 'RKDG{ord}'.format(ord=data.orderRK) 

        elif self.compareStr(argc,'-maxIters'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.maxIters = float(argv[1])

        elif self.compareStr(argc,'-errortype'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            if argv[1] in ['L1','L2','Linf']:
                data.error['L2'] = False
                data.error[argv[1]] = True
            else:
                print("Unknown error type %s %s" % (argc,argv[1]))
                print("Acceptable error types:'L1','L2','Linf'")
                return False

        elif self.compareStr(argc,'-outdir'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.outdir = argv[1]
            if not (data.outdir < data.outfile):
                data.outfile = data.outdir+data.outfile

        elif self.compareStr(argc,'-outfile'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.outfile = argv[1]
            if not (data.outdir < data.outfile):
                data.outfile = data.outdir+data.outfile
                data.exportData  = True

        elif self.compareStr(argc,'-storeall'):
            data.storeAll = True

        elif self.compareStr(argc,'-infile'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.infile = argv[1]

        elif self.compareStr(argc,'-info'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.info = ' '.join(argv[1:])

        elif self.compareStr(argc,'-iterationTol'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.iterationTol= float(argv[1])

        elif self.compareStr(argc,'-leftrightbc'):
            data.leftrightbc = True

        elif self.compareStr(argc,'-Dirbc'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            if argv[1].lower() in ['left','right','all']:
                data.orientDirbc[argv[1].lower()] = True
                data.leftrightbc = False
            else:
                print("Unknown boundary orientation %s %s" % (argc,argv[1]))
                print("Acceptable boundary orientations:'left','right','all'")
                return False

        elif self.compareStr(argc,'-Refbc'):
            # Reflective boundary conditions
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            if argv[1].lower() in ['left','right','all']:
                data.orientRefbc[argv[1].lower()] = True
                data.leftrightbc = False
            else:
                print("Unknown boundary orientation %s %s" % (argc,argv[1]))
                print("Acceptable boundary orientations:'left','right','all'")
                return False

        elif self.compareStr(argc,'-Neumbc'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            if argv[1] in ['left','right']:
                data.orientNeumbc = argv[1]
                data.leftrightbc = False
            else:
                print("Unknown boundary orientation %s %s" % (argc,argv[1]))
                print("Acceptable boundary orientations:'left','right'")
                return False

        elif self.compareStr(argc,'-Dirbctype'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            if argv[1].lower() in ['homogeneous','constant','piecewise','variable','none']:
                data.Dirbctype['homogeneous'] = False
                data.Dirbctype[argv[1].lower()] = True
            else:
                print("Unknown boundary type %s %s" % (argc,argv[1]))
                print("Acceptable boundary types:'homogeneous','constant','piecewise','variable','none'")
                return False
        elif self.compareStr(argc,'-Neumbctype'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            if argv[1].lower() in ['homogeneous','constant','piecewise','variable','none']:
                data.Neumbctype['homogeneous'] = False
                data.Neumbctype[argv[1].lower()] = True
            else:
                print("Unknown boundary type %s %s" % (argc,argv[1]))
                print("Acceptable boundary types:'homogeneous','constant','piecewise','variable','none'")
                return False

        elif self.compareStr(argc,'-sourcedata'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.source_data = [float(s) for s in argv[1:]]

        elif self.compareStr(argc,'-sourcepos'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.source_pos = [int(float(s)) for s in argv[1:]]

        elif self.compareStr(argc,'-Dirbcdata'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.Dirbc_data = [float(s) for s in argv[1:]]

        elif self.compareStr(argc,'-Neumbcdata'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            data.Neumbc_data = [float(s) for s in argv[1:]]

        elif self.compareStr(argc,'-sourcetype'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            if argv[1].lower() in ['homogeneous','constant','piecewise','variable','variable_s']:
                data.sourcetype['constant'] = False
                data.sourcetype[argv[1].lower()] = True
            else:
                print("Unknown source type %s %s" % (argc,argv[1]))
                print("Acceptable source types:'homogeneous','constant','piecewise','variable','variable_s'")
                return False

        elif (self.compareStr(argc,'-flowtype') or self.compareStr(argc,'-transtype')):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            if argv[1].lower() in ['constant','spatial','temporal','linear','nonlinear']:
                data.flowtype['nonlinear'] = False
                data.flowtype[argv[1].lower()] = True
            else:
                print("Unknown flow type %s %s" % (argc,argv[1]))
                print("Acceptable source types:'constant','spatial','temporal','linear','nonlinear'")
                return False

        elif self.compareStr(argc,'-CFmethod'):
            if len(argv)<2:
                self.printErrorMsg(argv)
                return False
            if argv[1].upper() in ['CF111','CF122','CF222','CF232','CF233','CF333','CF343','CF443','CK443','CF664']:
                data.cfmethod = data.cfmethod.fromkeys(data.cfmethod,False)
                data.cfmethodName = argv[1].upper()
                data.cfmethod[argv[1].upper()] = True
                if argv[1].upper() in ['CF122','CF222','CF232']:
                    data.orderRK = 2
                elif argv[1].upper() in ['CF233','CF333','CF343','CF443','CK443']:
                    data.orderRK = 3
                elif argv[1].upper() in ['CF664']:
                    data.orderRK = 4
            else:
                print("Unknown commutator-free method %s" % (argv[1].upper()))
                print("Acceptable method:'CF111','CF122','CF222','CF232','CF233','CF333','CF343','CF443','CK443','CF664'")
                return False

        else:
            print("Unknown parameter/key option: %s" % argc)
            self.show_usage()
            return False

        return True


# ===============================================================
class Parameters:

    # Constructor
    def __init__(self):
        self.storeAll = False         # to store just the final solution
        self.runtime_flag = False     # checks that runtime arguments are correct
        self.dt_max = 0.1
        self.dt_min = 1e-4
        
        self.leftrightbc = True # no periodic bc
        
        # homogeneous (0) , piece-wise constant (1), variable (2), none (3)
        self.Dirbctype = {'homogeneous':True,'constant':False,'piecewise':False,'variable':False,'none':False}   # Dirichlet boundary type
        self.Neumbctype = {'homogeneous':True,'constant':False,'piecewise':False,'variable':False,'none':False}  # Neumann boundary type

        # homogeneous (0) ,  constant (1), piece-wise constant(2), variable (3)   # Source type
        self.sourcetype = {'homogeneous':False,'constant':True,'piecewise':False,'variable':False,'variable_s':False}

        # constant (0), spatial (1), temporal (2), linear(3), nonlinear (4)
        self.flowtype = {'constant':False,'spatial':False,'temporal':False,'linear':False,'nonlinear':True}

        self.cfmethod = {'CF111':True,'CF122':False,'CF222':False,'CF232':False,'CF233':False,'CF333':False,'CF343':False,'CF443':False,'CF664':False}
        self.orderRK = 1         # order of RK method to use
        self.iterationTol = 1e-6 # tolerance for iterative methods
        self.maxIters = 10       # maximum no. of iterations for iterative methods
        self.cfmethodName = 'CF111'

        # type of errors to measure
        self.error = {'L1':False,'L2':True,'Linf':False}

        self.info = '' 
        self.infile = ''
        self.outdir = ''
        self.outputfile = ''

        self.exportData = False # to export into MATLAB
        self.printErrors = False

        self.DGlimiter = {'PP':False,'MPP':False}
        self.DGlimiterMinMax = [0.0,1.0]
        self.adaptiveStepping = False

        # number of elements of a given mesh-sizes (for non-uniform meshes)
        self.dx = [0.5]
        self.nx = [4]

        self.polyDeg = 1       # Polynomial degree
        self.xb = [-1.0,1.0]   # spartial domain endpoints [xl,xr]
        self.tb = [0.0,1.0]    # temporal domain endpoint [t0,tN]
        self.mplim = [0.0,1.0]    # MPP limit for RKDG limiter
        self.Cr = 1.0  # Courant number
        self.Ne = 4    # Number of elements
        self.Re = 1.0  # constant diffusion coeff.
        self.dT = 1    # length of time domain
        self.orderRK = 1 # order of RKDG method

        # number of elements with given piece-wise data
        self.source = []
        self.source_pos = []
        self.source_data = []
        self.Dirbc_data = []
        self.Neumbc_data = []

        # list of edges these bc conditions
        self.Dirbc = []
        self.Neumbc = []
        self.Refbc = []
        self.orientDirbc = {'left':False,'right':False,'all':False}  # left (0) or right (1)
        self.orientNeumbc = {'left':False,'right':False}  # left (0) or right (1)
        self.orientRefbc = {'left':False,'right':False,'all':False}  # left (0) or right (1)

    def reset_dT(self):
        self.dT = self.tb[1] - self.tb[0]
        self.dt_max = min(self.dT,self.dt_max)

    def check_dtmax(self):
        return (self.dtmax >= 0 and self.dtmax <= self.dT)

    def check_dtmin(self):
        return (self.dtmin >= 0 and self.dtmin <= self.dT)

    def set(self,flag):
        flag = True
