class Methods:

    def __init__(self):
        self.AI = None
        self.AE = None
        self.b = None
        self.c = None
        self.alpha = None
        self.beta = None
        self.order = None
        self.stifflyAccurateDIRK = False
        self.stifflyAccurateIMEX = False
        self.stifflyAccurateDIRKCF = False
        self.parms = []
        self.istart = 1
        self.jstart = 1

class DIRK_CF (Methods):

    def __init__(self, np, cfmethod, *argv):
        super().__init__()
        cf = None
        for arg in argv: self.parms.append(arg)
        if (cfmethod['CF111']):
            from Methods import CF111
            cf = CF111(np)
            self.stifflyAccurateDIRK = True
            self.stifflyAccurateIMEX = True
            self.stifflyAccurateDIRKCF = True
        elif (cfmethod['CF122']):
            from Methods import CF122
            cf = CF122(np)
        elif (cfmethod['CF222']):
            from Methods import CF222
            cf = CF222(np)
            self.stifflyAccurateDIRK = True
            self.stifflyAccurateIMEX = True
            self.stifflyAccurateDIRKCF = True
        elif (cfmethod['CF232']):
            from Methods import CF232
            cf = CF232(np)
        elif (cfmethod['CF233']):
            from Methods import CF233
            cf = CF233(np,self.parms)
        elif (cfmethod['CF333']):
            from Methods import CF333
            cf = CF333(np,self.parms)
            self.jstart = 0
        elif (cfmethod['CF343']):
            from Methods import CF343
            cf = CF343(np,self.parms)
            self.stifflyAccurateDIRK = True
        elif (cfmethod['CF443']):
            from Methods import CF443
            cf = CF443(np,self.parms)
            self.stifflyAccurateDIRK = True
            self.stifflyAccurateIMEX = True
        elif (cfmethod['CK443']):
            from Methods import CK443
            cf = CK443(np,self.parms)
            self.jstart = 0
            self.stifflyAccurateDIRK = True
            self.stifflyAccurateIMEX = True
        elif (cfmethod['CF664']):
            from Methods import CF664
            cf = CF664(np,self.parms)
            self.alpha = cf.alpha
            self.stifflyAccurateDIRK = True
            self.stifflyAccurateIMEX = True
            self.jstart = 0
        else:
            if not(cfmethod['CF111']): 
                print("Your desired method is not yet defined!")
            from Methods import CF111 # default option
            cf = CF111(np)
            self.stifflyAccurateDIRK = True
            self.stifflyAccurateIMEX = True
            self.stifflyAccurateDIRKCF = True

        self.AI = cf.AI; self.AE = cf.AE
        self.b = cf.b; self.c = cf.c
        self.beta = cf.beta; self.order = int(cf.order)

#================================================================================

class CF111:

    def __init__(self,np): # with numpy module
        # coefficients of DIRK-CF methods
        # from IMEX RK pairs by Ascher & co. 1997
        self.order = 1
        self.AI = np.array([[0.,0.],[0.,1.]]) # Implicit RK (DIRK)
        self.AE = np.array([[0.,0.],[1.,0.]]) # Explicit RK
        self.c = np.array([0.,1.])
        self.b = self.AI[-1]
        self.beta = self.AE[-1]

class CF122:

    def __init__(self,np): # with numpy module
        # coefficients of DIRK-CF methods
        # from IMEX RK pairs by Ascher & co. 1997
        self.order = 2
        self.AI = np.array([[0.,0.],[0.,.5]]) # Implicit RK (DIRK)
        self.AE = np.array([[0.,0.],[.5,0.]]) # Explicit RK
        self.c = np.array([0.,.5])
        self.b = np.array([0.,1.])
        self.beta = np.array([[0.,1.]])

class CF222:

    def __init__(self,np): # with numpy module
        # coefficients of DIRK-CF methods
        # from IMEX RK pairs by Ascher & co. 1997
        self.order = 2
        gma = (2.-np.sqrt(2))/2
        dta = 1. - 1./(2*gma)
        self.AI = np.array([[0.,0.,0.], [0.,gma,0.], [0.,1.-gma,gma]]) # Implicit RK (DIRK)
        self.AE = np.array([[0.,0.,0.], [gma,0.,0.], [dta,1.-dta,0.]]) # Explicit RK
        self.c = np.array([0.,gma,1.])
        self.b = self.AI[-1]
        self.beta = np.array([[dta,1.-dta,0.]])

class CF232:

    def __init__(self,np): # with numpy module
        # coefficients of DIRK-CF methods, using L-stable 2-stage order 2 DIRK
        # from IMEX RK pairs by Ascher & co. 1997
        self.order = 2
        gma = (2.-np.sqrt(2))/2
        dta = -2*np.sqrt(2)/3.
        self.AI = np.array([[0.,0.,0.], [0.,gma,0.], [0.,1.-gma,gma]]) # Implicit RK (DIRK)
        self.AE = np.array([[0.,0.,0.], [gma,0.,0.], [dta,1.-dta,0.]]) # Explicit RK
        self.c = np.array([0.,gma,1.])
        self.b = self.AI[-1]
        self.beta = np.array([[0.,1.-gma,gma]])

class CF233:

    def __init__(self, np, parms, x1 = 0.32192686, x2 = 0.26055649): # with numpy module
        # coefficients of DIRK-CF methods
        # from IMEX RK pairs by Ascher & co. 1997
        if(len(parms)>0): 
            if (len(parms)>1): x1 = parms[0]; x2 = parms[1]
            else: x1 = parms[0]
        self.order = 3
        gma = (3.+np.sqrt(3))/6
        self.AI = np.array([[0.,0.,0.], [0.,gma,0.], [0.,1.-2*gma,gma]]) # Implicit RK (DIRK)
        self.AE = np.array([[0.,0.,0.], [gma,0.,0.], [gma-1.,2*(1.-gma),0.]]) # Explicit RK
        self.c = np.array([0.,gma,1.-gma])
        self.b = np.array([0.,.5,.5])
        self.beta = np.zeros((2,3))
        self.beta[0,:-1] = [x1,x2]
        self.beta[0][-1] = (1./3 + (2*self.c-np.ones((3,))).dot(self.beta[0])) / (1-2*self.c[-1])
        self.beta[1] = self.b - self.beta[0]

class CF333:

    def __init__(self, np, parms, x1 = 0.422235, x2 = 0.35211646): # with numpy module
        # coefficients of DIRK-CF methods; A-stable 3-stage 3rd order DIRK - Griepentrog
        # from IMEX RK pairs by Ascher & co. 1997
        if(len(parms)>0): 
            if (len(parms)>1): x1 = parms[0]; x2 = parms[1]
            else: x1 = parms[0]
        self.order = 3
        bta = np.sqrt(3.)/3.
        self.AI = np.array([[0.,0.,0.],[-bta/2,(1+bta)/2,0.],[(3+5*bta)/2,-1-3*bta,(1+bta)/2]]) # Implicit RK (DIRK)
        self.AE = np.array([[0.,0.,0.],[.5,0.,0.],[-1.,2.,0.]]) # Explicit RK
        self.c = np.array([0.,.5,1.])
        self.b = np.array([1./6,2./3,1./6])
        self.beta = np.zeros((2,3))
        self.beta[0,:-1] = [x1,x2]
        self.beta[0][-1] = (1./3 + (2*self.c-np.ones((3,))).dot(self.beta[0])) / (1-2*self.c[-1])
        self.beta[1] = self.b - self.beta[0]

class CF343:

    def __init__(self, np, parms, x1 = 0.2580567, x2 = 0.85320449, x3 = -0.50889154): # with numpy module
        # coefficients of DIRK-CF methods; L-stable 3-stage 3rd order DIRK
        # from IMEX RK pairs by Ascher & co. 1997
        if(len(parms)>0): 
            if (len(parms)>2): x1 = parms[0]; x2 = parms[1]; x3 = parms[2]
            elif (len(parms)>1): x1 = parms[0]; x2 = parms[1]
            else: x1 = parms[0]
        self.order = 3
        gma = 0.4358665215
        dta = 0.5529291479
        b1 = -3*gma**2/2+4*gma-1./4
        a31 = (15./4-15*gma+21*gma**2/4)*dta - 7./2+13*gma-9*gma**2/2
        self.AI = np.array([[0.,0.,0.,0.],[0.,gma,0.,0.],[0.,(1-gma)/2,gma,0.],[0.,b1,1-gma-b1,gma]]) # Implicit RK (DIRK)
        self.AE = np.array([[0.,0.,0.,0.],[gma,0.,0.,0.],[a31,(1+gma)/2-a31,0.,0.],[1-2*dta,dta,dta,0.]]) # Explicit RK
        self.c = np.array([0.,gma,(1+gma)/2,1.])
        self.b = self.AI[-1]
        self.beta = np.zeros((2,4))
        self.beta[0,:-1] = [x1,x2,x3]
        self.beta[0][-1] = (1./3 + (2*self.c-np.ones((4,))).dot(self.beta[0])) / (1-2*self.c[-1])
        self.beta[1] = self.b - self.beta[0]

class CF443:

    def __init__(self, np, parms, x1 = 0.4336908, x2 = 1.11871198,x3 = 0.45944593, x4 = -1.10064624): # with numpy module
        # coefficients of DIRK-CF methods; L-stable 4-stage 4th order DIRK
        # from IMEX RK pairs by Ascher & co. 1997
        if(len(parms)>0): 
            if (len(parms)>3): x1 = parms[0]; x2 = parms[1]; x3 = parms[2]; x4 = parms[3]
            elif (len(parms)>2): x1 = parms[0]; x2 = parms[1]; x3 = parms[2]
            elif (len(parms)>1): x1 = parms[0]; x2 = parms[1]
            else: x1 = parms[0]
        self.order = 3
        self.AI = np.array([[0.,0.,0.,0.,0.],[0.,.5,0.,0.,0.],[0.,1./6,.5,0.,0.],[0.,-.5,.5,.5,0.],[0.,1.5,-1.5,.5,.5]]) # Implicit RK (DIRK)
        self.AE = np.array([[0.,0.,0.,0.,0.],[.5,0.,0.,0.,0.],[11./18,1./18,0.,0.,0.],[5./6,-5./6,.5,0.,0.],[.25,1.75,.75,-1.75,0.]]) # Explicit RK
        self.c = np.array([0.,.5,2./3,.5,1.])
        self.b = self.AI[-1]
        self.beta = np.zeros((2,5))
        self.beta[0,:-1] = [x1,x2,x3,x4]
#         b = self.AE[-1]
        self.beta[0][-1] = (1./3 + (2*self.c-np.ones((5,))).dot(self.beta[0])) / (1-2*self.c[-1])
        self.beta[1] = self.AE[-1] - self.beta[0]

class CK443:

    def __init__(self, np, parms, x1 = 1./5, x2 = 3./4, x3 = 2./3): # with numpy module
        # coefficients of DIRK-CF methods; L-stable 4-stage 4th order DIRK
        # from IMEX RK pairs by Ascher & co. 1997
        if(len(parms)>0): 
            if (len(parms)>2): x1 = parms[0]; x2 = parms[1]; x3 = parms[2]
            elif (len(parms)>1): x1 = parms[0]; x2 = parms[1]
            else: x1 = parms[0]
        self.order = 3
        self.AI = np.zeros((4,4)) # Implicit RK (SDIRK)
        self.AI[1] = [1767732205903./4055673282236,1767732205903./4055673282236,0.,0.]
        self.AI[2] = [2746238789719./10658868560708,-640167445237./6845629431997,1767732205903./4055673282236,0.]
        self.AI[3] = [1471266399579./7840856788654,-4482444167858./7529755066697,11266239266428./11593286722821,1767732205903./4055673282236]
        self.AE = np.zeros((4,4)) # Explicit RK
        self.AE[1] = [1767732205903./2027836641118,0.,0.,0.]
        self.AE[2] = [5535828885825./10492691773637,788022342437./10882634858940,0.,0.]
        self.AE[3] = [6485989280629./16251701735622,-4246266847089./9704473918619,10755448449292./10357097424841,0.]
        self.c = np.array([0,1767732205903./2027836641118,3./5,1.])
        self.b = self.AI[-1]
        self.beta = np.zeros((2,4))
        self.beta[0,:-1] = [x1,x2,x3]
        self.beta[0][-1] = (1./3 + (2*self.c-np.ones((4,))).dot(self.beta[0])) / (1-2*self.c[-1])
        self.beta[1] = self.AE[-1] - self.beta[0]

class CF664:

    def __init__(self, np, parms, x1 = -1./3, x2 = 1./2): # with numpy module
        # coefficients of DIRK-CF methods; L-stable 6-stage 4th order DIRK
        # from IMEX RK pairs by Kennedy and Carpenter
        if(len(parms)>0): 
            if (len(parms)>1): x1 = parms[0]; x2 = parms[1]
            else: x1 = parms[0]
        self.order = 4
        self.AI = np.zeros((7,7)) # Implicit RK (SDIRK)
        self.AI[2,1:3] = [1./4,1./4]
        self.AI[3,1:4] = [8611./62500,-1743./31250,1./4]
        self.AI[4,1:5] = [5012029./34652500,-654441./2922500,174375./388108,1./4]
        self.AI[5,1:6] = [15267082809./155376265600,-71443401./120774400,\
            730878875./902184768,2285395./8070912,1./4]
        self.AI[6,1:7] = [82889./524892,0.,15625./83664,69875./102672,-2260./8211,1./4]
        self.AE = np.zeros((7,7)) # Explicit RK
        self.AE[2][1] = 1./2
        self.AE[3,1:3] = [13861./62500,6889./62500]
        self.AE[4,1:4] = [-116923316275./2393684061468,-2731218467317./15368042101831,\
                9408046702089./11113171139209]
        self.AE[5,1:5] = [-451086348788./2902428689909,-2682348792572./7519795681897,\
                12662868775082./11960479115383,3355817975965./11060851509271]
        self.AE[6,1:6] = [647845179188./3216320057751,73281519250./8382639484533,\
      552539513391./3454668386233,3354512671639./8306763924573,4040./17871]
        self.c = np.array([0.,0., 1./2, 83./250, 31./50, 17./20,1.])
        self.b = np.array([0.,82889./524892,0.,15625./83664, 69875./102672, -2260./8211,1./4])
        self.beta = np.zeros((2,7))
        self.beta[0,5:] = [x1,x2]
        bb = np.array([x1,x2])
        b = np.array([.5-np.sum(bb),1./12-self.c[5:].dot(bb),-self.c[5:].dot(bb),-bb.dot(self.AE[5:].dot(self.c))])
        H = np.ones((4,4))
        H[1] = self.c[1:5]; H[2] = self.c[1:5]**2; H[3] = self.AE[1:5].dot(self.c)
        self.beta[0,1:5] = np.linalg.solve(H,b)
        self.beta[1] = self.b - self.beta[0]
        self.alpha = np.zeros((2,6))
        self.alpha[0,2:] = [1./6,-1./4,0.,5./12] # 4 free parameters
        Alp = self.AE[1:,2:]
        Alp[-1,1:] = self.alpha[0,2:]
        Gm = self.AE.dot(self.c)
        dta = Alp.dot(self.c[2:])
        dtah = Alp.dot(np.ones((5,)))
        self.alpha[0][1] = -1./(self.b[-1]*Gm[-1]) * (1./12 + self.b[1:-1].dot(Gm[1:-1]*self.AE[1:-1,1]) - self.b[1:].dot(self.c[1:]*dta + Gm[1:]*(self.c[1:]-dtah)))
        self.alpha[1,1:] = self.AE[-1,1:-1] - self.alpha[0,1:]


#================================================================================

class RKDG (Methods):

    def __init__(self, np, orderRK):
        super().__init__()
        obj = None
        if (orderRK == 1):
            from Methods import Euler
            obj = Euler(np)
        elif (orderRK == 2):
            from Methods import Midpoint
            obj = Midpoint(np)
        elif (orderRK == 3):
            from Methods import RKDG3
            obj = RKDG3(np)
        elif (orderRK == 4):
            from Methods import RKDG4
            obj = RKDG4(np)
        else:
            print("Your desired method is not defined!")
            from Methods import Euler # default 
            obj = Euler(np)
        self.alpha = obj.alpha; self.beta = obj.beta
        self.b = obj.b; self.c = obj.c; self.order = orderRK
#================================================================================

class Euler:

    def __init__(self, np):
        # Euler method for RKDG1
        self.alpha = np.array([[1.]])
        self.beta = np.array([[1.]])
        self.c = np.array([0.])
        self.b = np.array([1.])

class Midpoint:

    def __init__(self, np):
        # Midpoint method for RKDG2
        self.alpha = np.array([[1.,0.],[.5,.5]])
        self.beta = np.array([[1.,0.],[0.,.5]])
        self.c = np.array([0.,.5])
        self.b = np.array([0.,1.])

class RKDG3:

    def __init__(self, np):
        # method for SSPRK(3,3) or RKDG3
        self.alpha = np.array([[1., 0., 0.], [.75,.25,0.], [1./3, 0., 2./3]])
        self.beta = np.array([[1.,0.,0.], [0.,.25,0.], [0.,0.,2./3]])
        self.c = np.array([0.,.5,1.])
        self.b = np.array([0.,1.,1.])

# class RKDG3:
# 
#     def __init__(self, np):
#         # method for SSPRK(5,3) or RKDG3
#         self.alpha = np.array([[1.,0.,0.,0.,0.], [0.,1.,0.,0.,0.], [0.56656131914033,0.,0.43343868085967,0.,0.], [0.09299483444413,0.00002090369620,0.,0.90698426185967,0.], [0.00736132260920,0.20127980325145,0.00182955389682,0.,0.78952932024253] ])
#         self.beta = np.array([[0.37726891511710,0.,0.,0.,0.], [0.,0.37726891511710,0.,0.,0.], [0.,0.,0.16352294089771,0.,0.], [0.00071997378654,0.,0.,0.34217696850008,0.], [0.00277719819460,0.00001567934613,0.,0.,0.29786487010104] ])
#         self.c = np.zeros((5,))
#         self.b = np.zeros((5,))

class RKDG4:

    def __init__(self, np):
        # method for SSPRK(5,3) or RKDG4 cf. Gottlieb et al.2009
        self.alpha = np.array([[1., 0., 0., 0., 0.], [0.444370493651235, 0.555629506348765, 0., 0., 0.], [0.620101851488403, 0., 0.379898148511597, 0., 0.], [0.178079954393132, 0., 0., 0.821920045606868, 0.], [0., 0., 0.517231671970585, 0.096059710526147, 0.386708617503269] ] )
        self.beta = np.array([[0.391752226571890, 0., 0., 0., 0.], [0., 0.368410593050371, 0., 0., 0.], [0., 0., 0.251891774271694, 0., 0.], [0., 0., 0., 0.544974750228521, 0.], [0., 0., 0., 0.063692468666290, 0.226007483236906] ] )
        A = np.linalg.solve(np.eye(4) - self.alpha[1:,1:], self.beta[1:,1:])
        self.c = np.zeros((5,))
        self.c[1:] = np.sum(A,1)
        self.b = np.zeros((5,))
        self.b[1:] = self.beta[-1,1:] + self.alpha[-1,1:].dot(A)
