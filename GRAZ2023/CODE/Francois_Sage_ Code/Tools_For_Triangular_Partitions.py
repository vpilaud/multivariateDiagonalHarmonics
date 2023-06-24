from sage.graphs.graph_plot import GraphPlot

options = {
    'vertex_size': 100,
    'vertex_labels': True,
    'layout': None,
    'edge_style': 'solid',
    'edge_color': 'black',
    'edge_colors': None,
    'edge_labels': False,
    'iterations': False,
    'tree_orientation': 'up',
    'heights': None,
    'graph_border': False,
    'talk': False,
    'color_by_label': False,
    'partition': None,
    'dist': .075,
    'max_dist': 2,
    'loop_size': .075,
    'edge_labels_background': 'transparent'}


def ab_SturmFactors(n):
    return Set([w for w in Words([a,b],n) if w.is_sturmian_factor() and w[0]==a and w[-1]==b])

def ab_Factors(w):
    return [u for u in w.factor_iterator() if not u.is_empty() and  u[0]==a and u[-1]==b]

def AddableCells(mu):
    mu=Partition(list(mu))
    return [(j+1,i+1) for i,j in mu.outside_corners() 
            if Is_Triangular(mu.add_cell(i,j))]

def area_tau(alpha,tau):
    return tau.size()-alpha.size()


def A_comme_pente(tau,sl):
    tau=Partition(tau)
    return t_min(tau)<sl<t_max(tau)

def All_Integral_Triangular_Partitions(N):
    return Set([tau for k in range(N+1) for tau in Integral_Triangular_Partitions(k)])


@cached_function
def Bizley(k,n):
    d=gcd(k,n)
    a,b=k/d,n/d
    return add((a**(-mu.length()))/mu.aut()*mul(e[b*c](a*c*Z) for c in mu) for mu in Partitions(d))


def Cat_Mat_q(tau):
    tau=Partition(tau)
    return (matrix([[q_binomial(tau[i]+1,i-j+1)*q**(binomial(i-j+1,2))
                                       for j in range(tau.length())] 
                                      for i in range(tau.length())]))
def Cat(k,h=False):
    return Integer(1)/(k+1)*binomial(2*k,k,hold=h)

def Cat_Mat(tau):
    tau=Partition(tau)
    return (matrix([[hold_bin(tau[i]+1,i-j+1)*q**(binomial(i-j+1,2))
                                       for j in range(tau.length())] 
                                      for i in range(tau.length())]))

@cached_function
def Cat_rec(nu):
    nu=Partition(nu)
    if not Is_Triangular(nu):
        print(nu)
        return 0
    elif Partition(nu)==zero:
        return 1
    else:
        return add(Cat_rec(alpha)*Cat_rec(beta) 
                   for alpha,beta in Decoupes(nu))+Cat_rec(Interieur(nu))
    
def Cat_rec_sup(tau):
    tau=Partition(tau)
    rep=Set([])
    if not Is_Triangular(tau):
        return  rep
    elif Partition(tau)==zero:
        return rep
    else:
        for Dtau in Decoupes(tau):
            rep=rep.union(Set(Dtau))
            for nu in Dtau:
                rep=rep.union(Cat_rec_sup(nu))
        rep=rep.union(Set([Interieur(tau)]))
        rep=rep.union(Cat_rec_sup(Interieur(tau)))
        return rep
    
@cached_function
def Cat_Schur(tau):
    return InSchur(Cat_tau(tau))

def Cellules(mu):
    return [(b+1,a+1) for a in range(mu.length()) for b in range(mu[a])]

def Centre(tau):
    tau=Partition(tau)
    n=tau.size()
    if n==0:
        return [0,0]
    elif not Is_Dominant(tau):
        c=Centre(tau.conjugate())
        return [-c[0],c[1]]
    elif tau.length()==1 and tau.size()>2:
        return [(-n+2),n]
    elif tau.length()==1 and tau.size()==2:
        return [(-n+3/2),n]
    else:
        return [n*(2*t_tau(tau)-1),n]

def CentreLog(tau):
    tau=Partition(tau)
    n=tau.size()
    if n==0:
        return [0,-1/2]
    elif n==1:
        return [0,0]
    elif not Is_Dominant(tau):
        c=CentreLog(tau.conjugate())
        return [-c[0],c[1]]
    elif tau.length()==1 and n==2:
        return [-log(3/2,2),log(3/2,2)]
    elif tau.length()==1 and n>2:
        return [-log(n-1,2),log(n-1,2)]
    else:
        t=t_tau(tau)
        return [log(t/(1-t),2),2*log(n,2)+log(t*(1-t),2)]

def Coins(mu):
    mu=Partition(mu)
    return Set([tuple(vector((b,a))+vector((1,1))) 
            for a,b in mu.corners()])

def Cone(nu):
    return (t_min(nu),t_max(nu))

@cached_function
def Contruit_A(tau):
    tau=Partition(tau)
    if not Is_Dominant(tau):
        return Contruit_A(tau.conjugate())
    else:
        n=tau.length()+1
        A=C_tau(tau,q=q,t=t,Schur='true')
        if n<=3:
            return A
        else:
            F=qt_Frob(tau,N=n)
            for k in range(n//2):
                difference=InSchur(F.scalar(hook(n,k)))-A.skew_by(e[n-1-k])
                A=A+add(c*s(mu.Plus([1 for i in range(n-1-k)])) 
                        for mu,c in difference)
            return A
    

def CornerCut(tau):
    tau=Partition(tau)
    pp=(grid(ceil(r_tau(tau)),ceil(s_tau(tau)))
            +line_rs_tau(tau)+tau.diagram()
            +OutDiagram(OutDiag(tau),ceil(s_tau(tau)),ceil(r_tau(tau))))
    pp.axes(show=False)
    return pp.show(figsize=3)

def Coupures(n):
    return [0]+[((t_max(tau))) for tau in reversed(TriangularPartitions(n))] 

@cached_function
def C_tau(tau,q=q,t=t,Schur=False):
    tau=Partition(tau)
    if Schur==False:
        return add(q**area_tau(alpha,tau)*t**dinv_tau(alpha,tau) 
                   for alpha in Dyck_tau(tau))
    else:
        return InSchur(C_tau(tau,q=q,t=t))

@cached_function
def C_hook_formule(tau):
    tau=Partition(tau)
    return s[tau.size()]+ add(s[k,1] for k in range(tau.size()-tau.length(),tau.size()-1))

@cached_function
def Cp_tau(tau,q=q):
    return add(q**(mu.size()) for mu in Dyck_tau(tau)) 
    
@cached_function
def Cp_tau_det(tau,q=q):
    tau=Partition(tau)
    return det(matrix([[q_binomial(tau[i]+1,i-j+1,q=q)*q**(binomial(i-j+1,2))
                                       for j in range(tau.length())] 
                                      for i in range(tau.length())]))

def C_q(tau,q=q):
    return add(q**(mu.size()) for mu in Dyck_tau(tau)) 
    
def C_q_det(tau,q=q):
    tau=Partition(tau)
    return det(matrix([[q_binomial(tau[i]+1,i-j+1)*q**(binomial(i-j+1,2))
                                       for j in range(tau.length())] 
                                      for i in range(tau.length())]))

def C_tau_two_parts(tau):
    tau=Partition(tau)
    if tau.length()<=1:
        return s[tau.size()]
    elif not Is_Dominant(tau):
        return Alt_formule(tau.conjugate())
    elif tau.length()==2:
        n=tau.size()
        return add(s([n-2*i,i]) for i in range(tau[1]+1) if n>=3*i)
    else:
        show((tau,'Calcul direct'))   
        
def Decoupes(nu):
    return [(remove(shadow(nu,(1,j+1)),Diagonale(nu),j),shadow(nu,(i+1,1))) 
            for i,j in Diagonale(nu)]   

def delta(r,s):
    return Partition([floor((r*(s-k))/s) for k in range(1,floor(s+1))])

def Delta_p_Schur(f,G):
    return add(tensor([InSchur(c),s(mu)]) for mu,c in Delta(f(Z-1),G))

def DessinPoset(P):
    Edes=[Part_to_num(mu) for mu in P]
    Rdes=[(Part_to_num(mu),Part_to_num(nu)) for mu,nu in P.cover_relations()]
    return GraphPlot(Poset([Edes,Rdes],cover_relations=True).cover_relations_graph().to_undirected(), options).plot()

def DessinRessort(P,epaisseur=1):
    Edes=[Part_to_num(mu) for mu in P]
    Rdes=[(Part_to_num(mu),Part_to_num(nu)) for mu,nu in P.cover_relations()]
    return SpringDrawPoset(Poset([Edes,Rdes],cover_relations=True),epais=epaisseur)

def Diagonale(mu):
    mu=Partition(mu)
    u=RemovableCells(mu)[0]
    v=RemovableCells(mu)[-1]
    return [c for c in Segment(u,v) if c in Coins(mu)]

def dinv_tau(alpha,tau):
    return Dinv_tau(alpha,tau).cardinality()

def Dinv_tau(alpha,tau):
    return Set([c for c in Cellules(alpha) 
                if low_t(c,alpha) <= t_tau(tau) < top_t(c,alpha)])  # Avant +epsilon

def Dyck_tau(tau):
    tau=Partition(tau)
    return Set([Partition(list(mu)) for k in range(tau.size()+1) for mu in Partitions(k, outer=tau)])

def enleve(u,nu):
    nu=Partition(nu)
    i,j=u
    return nu.remove_cell(j-1,i-1)

@cached_function
def Est_Concave(tau):
    tau=Partition(tau)
    if tau.size()<=3:
        return True
    else:
        tau=Partition(tau)
        cell_tau=Set([tuple(c) for c in tau.cells()])
        rho=Partition([max(tau)+1 for i in range(tau.length()+1)])
        cell_rho=Set([tuple(c) for c in rho.cells()])
        P_tau=Polyhedron(vertices = cell_rho.difference(cell_tau))
        cell_P_tau=Set([tuple(c) for c in P_tau.integral_points()])
        return cell_tau==cell_rho.difference(cell_P_tau)

def Est_Convexe(mu):
    mu=Partition(mu)
    cell_mu=Set([tuple(c) for c in mu.cells()])
    P_mu=Polyhedron(vertices = cell_mu)
    cell_P_mu=Set([tuple(c) for c in P_mu.integral_points()])
    return cell_mu==cell_P_mu

def Est_Triangulaire(mu):
    return Est_Convexe(mu) and Est_Concave(mu)

@cached_function
def E_via_Caracterisation_Delta(n):
    Temp=Etend_Tensor(Delta_p_Schur(e[0],e[n]),n-1)
    for i in range(1,n):
        Temp=Temp+Etend_Tensor(Delta_p_Schur(e[i],e[n])-Temp.skew_by(e[n-1-i]).restrict_length(2),n-1-i)
    return Temp  

def Etend(F,k):
    return add(c*s(Partition([k]+list(mu.conjugate())).conjugate()) for mu,c in F)

def Etend_Tensor(F,k):
    return add(c*tensor([Etend(s(mu),k),s(nu)]) for (mu,nu),c in F)

@cached_function
def E_v(v,Formal=False):
    n=len(v)
    rep=s(add(Weight(mu,v)/T_mu(mu)*H(mu.conjugate())
                                  for mu in Partitions(n)))
    if Formal:
        return Formal_Symmetric(rep)
    else:
        return rep



def Farey_itere(n,L=[]):
    if n==0:
        return L
    else:
        return Farey_itere(n-1,Farey(L))

def FareyAdd(a,b):
    return (numerator(QQ(a))+numerator(QQ(b)))/(denominator(QQ(a))+denominator(QQ(b)))

def Farey(L):
    return [0]+[FareyAdd(L[i],L[i-1]) for i in range(1,len(L))]+[1]


def Formal_Symmetric(F):
    return add(tensor([InSchur(c),s(nu)]) for nu,c in F)

def Formule(j,N=var('n')):
    return (Cat(N)-add(Cat(k)*Cat(N-k-1) for k in range(j)))

def Formule_tau(tau):
    tau=Partition(tau)
    j=binomial(tau.length()+1,2)-tau.size()
    N=tau.length()+1
    return Formule(j,N)

def Formule_q(j,N=var('n')):
    return (q_Cat(N)-add(q**(k)*q_Cat(k)*q_Cat(N-k-1) for k in range(j)))


def grid(k,n):
    return (add(line([(0,i),(k,i)],color='lightgrey') for i in range(n))
            +add(line([(i,0),(i,n)],color='lightgrey') for i in range(k)))

def Hab(a,b,N=15,thk=1.5):
    if a==0 and b==0:
        return (parametric_plot((-1,t),
                           (t,-1,log(N,2)),
                           thickness=thk)+parametric_plot((t,-1),
                           (t,-1,log(N,2)),
                           thickness=thk))
    elif a==0 and not b==0:
        return parametric_plot((t,log(b,2)),
                           (t,-1,N),
                           thickness=thk/5)
    elif not a==0 and  b==0:
        return parametric_plot((log(a,2),t),
                           (t,-1,N),
                           thickness=thk/5)
    else:
        return parametric_plot((log(a*exp(t)+a,2),log(b*exp(-t)+b,2)),
                           (t,log(b/(N-b)),log((N-a)/a)),
                           thickness=(thk/(a+b)))

def Habpoint(a,b,N=15,thk=1,col='blue'):
    return parametric_plot((log(a*exp(t)+a),log(b*exp(-t)+b)),
                           (t,log(b/(N-b)),log((N-a)/a)),
                           thickness=(thk/(a+b)),color=col)

def HasseDraw(P,grosseur=300,vertex=False):
    return P.cover_relations_graph().plot(prog='neato',
                                          dim=3,
                                          vertex_labels=vertex,
                                          vertex_size=grosseur,
                                          vertex_color='yellow')

def hold_bin(a,b):
    if b<0 or b>a:
        return 0
    elif b==0 or a==b:
        return 1
    else:
        return binomial(a,b,hold=True)


def Inf_Aire_BaseHauteur(x,y):
    return (x*y+1+gcd(x-1,y-1))/2  

def InSchur(F):
    try:
        return s(Sym(S.from_polynomial(Rqt(F).numerator()))).restrict_partition_lengths(2,exact=false)
    except:
        return F

@cached_function
def Integral_Triangular_Partitions(n):
    if n<=4:
        return Set(Partitions(n))
    else:
        val=Set([delta(a,b) for a in range(1,2*n+1) for b in range(1,n+1) if delta(a,b).size()==n])
        return val.union(Set([v.conjugate() for v in val]))

@cached_function  
def Interieur(nu,k=1):
    if k<=1:
        rho=Partition(nu)
        for i,j in Diagonale(nu):
            rho=rho.remove_cell(j-1,i-1)
        return rho
    else:
        return Interieur(Interieur(nu),k-1)
    
def IntersectHab(u,v):
    (a,b)=u
    (c,d)=v
    return vector([log((b*c-a*d)/(b-d)),log((b*c-a*d)/(c-a))])

def IntervalAire_BaseHauteur(x,y):
    return (Inf_Aire_BaseHauteur(x,y),Sup_Aire_BaseHauteur(x,y))

def IntervalAire_BaseHauteur_qcount(x,y):
    return add(q**(mu.size()) for mu in Set(PartagesTriangles_BaseHauteur(x,y)))

def IntervalAire_BaseHauteur_Integral_tcount(x,y):
    return add(t**(mu.size()) for mu in Set(PartagesTriangles_BaseHauteur(x,y)) if Is_Integral(mu))

def Inv_Delta(f,G):
    return s(add(c*H(mu)/(f(B_mu(mu)*Un).scalar(Un)) for (mu,c) in H(G)))

def Is_Triangular(tau):
    tau=Partition(tau)
    if tau.size()==0:
        return True
    return t_min(tau)<t_max(tau)

def Is_Dominant(nu):
    return nu.dominates(nu.conjugate())

@cached_function
def Is_Integral(tau):
    if tau.size()<2:
        return True
    return tau in Integral_Triangular_Partitions(tau.size())

def Le_n(mu):
    mu=Partition(mu)
    if not Is_Concave(mu):
        return None
    else:
        l=mu.length()
        for n in range(l+1,2*l+1):
            if Is_Concave(mu.AjouteColonne(n)):
                return n
def line_rs_tau(tau):
    return line([(r_tau(tau),0),(0,s_tau(tau))])

def low_t(c,alpha):
    leg=alpha.leg_length(c[1]-1,c[0]-1)
    arm=alpha.arm_length(c[1]-1,c[0]-1)
    return leg/(arm+leg+1)

def Max(L):
    return max(L+[0])

def Milieu_Aire_BaseHauteur(x,y):
    return Integer(Inf_Aire_BaseHauteur(x,y)+Sup_Aire_BaseHauteur(x,y))//2

def Min(L):
    return min(L+[infinity])

def Montre_A(tau):
    tau=Partition(tau)
    return table([(Contruit_A(tau).skew_by(e[mu[0]-1]),a**(mu[0]-1)) 
                  for mu,c in qt_Frob(tau,N=tau.length()+1) if Is_Hook(mu)])

def NumTrianglesLargeurFixee(x):
    return 1/2+1/2*add((x-i+1)*euler_phi(i) for i in range(x+1))

def OrdreSchur(mu,nu):
    return min([0]+(C_tau(nu,Schur=True)-C_tau(mu,Schur=True)).coefficients())>=0


def OutDiagram(mu,k,n):
    mu=Partition(mu)
    return plot(add(carre(c,col='yellow') for c in RectangleCellules(k,n).difference(Set(Cellules(mu)))))

def OutDiag(tau):
    return Partition([tau[0]+2]+[c+1 for c in tau]+[1])


var('a','b')

def Park_tau(tau,n):
    return add(ParkFrob_mu(mu,n) for mu in Dyck_tau(tau))

def Parking_gen(tau):
    tau=Partition(tau)
    temp=add(q**area_tau(alpha,tau)*e(s[alpha.Plus(Partition([2*tau.size()+1]).conjugate())].skew_by(s[alpha])) 
             for alpha in Dyck_tau(tau))
    return add(tensor([pol_to_s(c),e[Bar(mu)]]) for mu,c in temp)

def ParkFrob_mu(mu,n):
    mu=Partition(mu)
    return s(mu.AjouteColonne(n)).skew_by(s(mu))

@cached_function
def part_to_word(mu):
    var('a','b')
    if mu.length()==0:
        return Word([])
    else:
        L=[0]+list(reversed(mu))
        DL=[L[i]-L[i-1] for i in range(1,len(L))]
        return mul(Word([a for i in range(DL[k])]+[b]) for k in range(len(DL)))
    
@cached_function
def PartagesTriangles_BaseHauteur(x,y):
    return Set([tau for k in range(Inf_Aire_BaseHauteur(x,y),Sup_Aire_BaseHauteur(x,y)+1)
            for tau in PartagesTriangulaires(k) if tau[0]==x and tau.length()==y])

@cached_function
def part_to_vect(mu,n=0):
    mu=Partition(mu)
    if mu.length()==0:
        return vector([0 for i in range(n)])
    else:
        if n==0: n=mu.length()+1
        L=list(mu)+[0]
        DL=[L[i]-L[i+1] for i in range(len(L)-1)]
        return vector([0]+DL+[0 for i in range(n-len(DL)-1)])
    
def Part_diff_vect(tau,n=None):
    tau=Partition(tau)
    w=list(tau)+[0]
    v=[w[0]]+[w[i]-w[i+1] for i in range(tau.length())]
    if n==None:
        n=Le_n(tau)
    return tuple(v+[0 for i in range(n-len(v))])

def Part_from_diff(v):
    return Partition([v[0]-add(v[1:i]) for i in range(1,v.__len__())])

def PermutationsTriangulaires(tau):
    tau=Partition(tau)
    return [Permutation(RSK_inverse(t1,t2)[1]) for t1 in TableauxTriangulaires(tau) 
            for t2 in TableauxTriangulaires(tau)]

def petit_carre(c,g=1,d=(0,0),col='yellow'):
    a,b=c
    dx=d[0]
    dy=d[1]
    return polygon([((a-1)*g+dx,(b-1)*g+dy),((a-1)*g+dx,b*g+dy),(a*g+dx,b*g+dy),(a*g+dx,(b-1)*g+dy)],
                   edgecolor="black",
                   color=col,axes=false,thickness=g)

def petit_diag(tau,g=1,d=(0,0)):
    tau=Partition(tau)
    delt=vector([tau.conjugate().length(),tau.length()])/2
    delt=vector([1,1])/2
    return add(petit_carre([c[0]-delt[0],c[1]-delt[1]],g=g,d=vector(d)) for c in Cellules(tau))

def pol_to_s(pol):
    return add(c[0]*s[c[1].degree()] for c in pol.numerator())

def Position(tau):
    tau=Partition(tau)
    if tau==Partition([0]):
        val=vector((.5,.5))
    elif tau.length()==1 or tau.conjugate().length()==1:
        L=AddableCells(tau)
        u=L[0]
        v=L[1]
        val=IntersectHab(u,v)-vector((.25,.25))/(tau.size())
    else:
        R=RemovableCells(tau)
        A=AddableCells(tau)
        if R.__len__()==2 and A.__len__()==2:
            val=(IntersectHab(A[0],A[1])+IntersectHab(R[0],R[1]))/2
        elif A.__len__()==2:
            val=IntersectHab(A[0],A[1])-vector((.25,.25))/(tau.size())
        else:
            val=IntersectHab(R[0],R[1])+vector((.4,.4))/(tau.size())
    return (RR(val[0]),RR(val[1]))



def P_tau(tau,n=None,x=1):
    if n==None:
        n=tau.length()+1
    return Sym_P_tau(tau,n).scalar(p[1]**n).substitute({q:x})

@cached_function    
def Qbarpoint(n,Thk=1.5):
    return add(Habpoint(a,b,N=2*n+1,thk=Thk) 
               for a in range(1,n+1) 
               for b in range(1,n+1))

def Qbar(n,Thk=1.5):
    return add(Hab(a,b,N=2*n+1,thk=Thk) 
               for a in range(1,n+1) 
               for b in range(1,n+1) if not (a==0 and b==0))

def qt_Frob(tau,N=n):
    tau=Partition(tau)
    return Xsi(alpha_tau(tau.conjugate()),n=N)

def q_Formule_tau(tau):
    tau=Partition(tau)
    j=binomial(tau.length()+1,2)-tau.size()
    N=tau.length()+1
    return Formule_q(j,N)/q**(j)

@cached_function
def q_Cat(n):
    return C_tau(delta(n,n),q=q)

def RectangleCellules(k,n):
    return Set([tuple((i,j)) for i in range(1,n+1) for j in range(1,k+1)])
               
def Relations(n):
    return add(line([Position(mu),Position(nu)],zorder=0,thickness=.5) 
               for k in range(1,n+1) for mu in TriangularPartitions(k-1) 
       for nu in TriangularPartitions(k) 
            if nu.contains(mu) 
            and mu.conjugate().length()<=6
         and nu.conjugate().length()<=6
         and mu.length()<=6 
         and nu.length()<=6)


def RemovableCells(mu):
    mu=Partition(mu)
    return [(j+1,i+1) for i,j in mu.corners() 
            if Is_Triangular(mu.remove_cell(i,j))]

def remove(alpha,cell_list,k):
    rho=Partition(alpha)
    for i,j in cell_list:
        if (i,j-k) in Cellules(rho):
            rho=rho.remove_cell(j-k-1,i-1)
    return rho

def r_tau(tau):
    tau=Partition(tau)
    a=t_tau(tau)
    b=1-a
    c=max([a*i+b*j for (i,j) in Cellules(tau)])
    return c/a

def Schur_Frob(tau,N):
    tau=Partition(tau)
    return add(tensor([InSchur(c),s(mu)]) for mu,c in Xsi(alpha_tau(tau),n=N))

def Segment(u,v):
    return Set([(i,j) for i in range(min(u[0],v[0]),max(u[0],v[0])+1) 
            for j in range(min(u[1],v[1]),max(u[1],v[1])+1) if SurLaDroite((i,j),u,v)])

def shadow(mu,c):
    mu=Partition(mu)
    i,j=c
    return Partition([mu[k]-(i-1) for k in range(j-1,mu.length()) if mu[k]>=(i-1)])

def SommeSturmFactors(x,y,n):
    if n==0:
        return 1
    else:
        return add(SturmFactorsNum(x,y, i, n-i) for i in range(n))

def SpringDrawPoset(P,epais=1):
    return P.cover_relations_graph().to_undirected().show(method="js",
                                                          vertex_size=3,edge_thickness=epais,link_distance=50)

@cached_function
def SturmFactorsNum(x,y,i,j) :
    if (x,y) == (b,b) : 
        return SturmFactorsNum(a,a,j,i)
    elif (x,y) == (b,a) : 
        return SturmFactorsNum(a,b,j,i)
    elif (x,y) == (a,a) : 
        if i==0 : 
            return 0    
        elif i==1 :
            if j==0 :   
                return 1     
            else :       
                return 0           
        elif (i >= j+1) :  
            return (SturmFactorsNum(a,a,i-j-1,j) 
                    + SturmFactorsNum(a,b,i-j-1,j) 
                    + SturmFactorsNum(a,b,j,i-j-1)
                    + SturmFactorsNum(b,b,i-j-1,j))
        else :          
            return SturmFactorsNum(a,a,i, j+1-i)       
    elif (x,y) == (a,b) : 
        if i==0 or j==0:       
            return 0 
        elif (i>=j) :    
            return SturmFactorsNum(a,b,i-j,j) + SturmFactorsNum(b,b,i-j, j)
        else :          
            return SturmFactorsNum(a,a,i,j-i) + SturmFactorsNum(a,b,i,j-i)   
        
def Super(tau):
    tau=Partition(tau)
    return add(c*a**(mu[Integer(0)]-Integer(1)) 
               for mu,c in qt_Frob(tau,N=tau.length()+Integer(1)) if Is_Hook(mu))

def Support_rec(tau):
    return sorted(Cat_rec_sup(tau))+[Partition(tau)]

def Sup_Aire_BaseHauteur(x,y):
    return ((x+1)*(y+1)-1-gcd(x,y))/2 


def Super_Schur(tau):
    tau=Partition(tau)
    return table([(InSchur(c),a**(mu[0]-1)) for mu,c in qt_Frob(tau,N=tau.length()+1) if Is_Hook(mu)])

def SurLaDroite(w,u,v):
    a,b=u
    c,d=v
    x,y=w
    return det(matrix([[x-a,c-a],
                       [y-b,d-b]]))==0

def Sym_P_tau(tau,n=None,x=q):
    tau=Partition(tau)
    if n==None:
        n=tau.length()+1
    return add(x**area_tau(Partition(mu),tau)*Sym_Park(mu,n) for mu in Dyck_tau(tau))

def Sym_Park(mu,n):
    mu=Partition(mu)
    return e(s(mu.Plus(Partition([n]).conjugate())).skew_by(s(mu)))


def s_tau(tau):
    tau=Partition(tau)
    a=t_tau(tau)
    b=1-a
    c=max([a*i+b*j for (i,j) in Cellules(tau)])
    return c/b



@cached_function
def TableauxTriangulaires(tau):
    tau=Partition(tau)
    if tau==zero:
        return Set([])
    elif tau==Partition([1]):
        return Set([StandardTableau([[1]])])
    else:
        rep=Set([])
        n=tau.size()
        for c in RemovableCells(tau):
            nu=tau.remove_cell(c[1]-1,c[0]-1)
            rep=Set([t.add_entry((c[1]-1,c[0]-1),n) 
                 for t in TableauxTriangulaires(nu)]).union(rep)
        return rep

@cached_function
def Triangles(n,condition=(lambda tau:True)):
    return [tau for k in range(n+1) 
            for tau in TriangularPartitions(k) if condition(tau)]

def TrianglesLargeurFixee(x):
    return Set([tau for y in range(1,x+1) for tau in PartagesTriangles_BaseHauteur(x,y)])

@cached_function
def TriangularPartitions(n):
    if n<4:
        return sorted(Partitions(n))
    else:
        return sorted(Set([mu.add_cell(c[0])
                       for mu in TriangularPartitions(n-1)
                       for c in mu.addable_cells()
                       if Is_Triangular(mu.add_cell(c[0]))]))
    
PartagesTriangulaires=TriangularPartitions
    
@cached_function
def Triangular_Partitions(n,condition=(lambda tau:True)):
    return  sorted([mu for mu in Partitions(n) if Is_Triangular(mu) and condition(mu)])

def Triangular_Partition_pente(sl,n):
    if n==0:
        return Partition([0])
    for nu in TriangularPartitions(n):
        if t_min(nu)<sl<t_max(nu):
            return nu

def top_t(c,alpha):
    leg=alpha.leg_length(c[1]-1,c[0]-1)
    arm=alpha.arm_length(c[1]-1,c[0]-1)
    return (leg+1)/(arm+leg+1)

def Tous_Partages(n):
    return list([mu for k in range(n+1) for mu in Partitions(k)])



def Tronque(F,k):
    return add(c*s(mu) for mu,c in F if mu.length()==2 and mu[1]==k)

def t_min(tau):
    tau=Partition(tau)
    if tau.is_empty():
        return 0
    else:
        return QQ(Max([low_t(c,tau) for c in Cellules(tau)]))

def t_max(tau):
    tau=Partition(tau)
    if tau.is_empty():
        return 1
    else:
        return QQ(Min([top_t(c,tau) for c in Cellules(tau)]))

def t_tau(tau):
    return QQ((t_min(tau)+t_max(tau))/2)

@cached_function
def word_to_part(w):
    return Partition(reversed([u.count(a) 
                               for u in w.prefixes_iterator() 
                               if len(u)>0 and u[-1]==b]))

def Yakob_Triangle(n,k=0):
    E=[nu for nu in PartagesTriangulaires(n) if Is_Dominant(nu) and nu.length()>=k]
    R=[(mu,nu) for mu in E for nu in E  
       if OrdreSchur(mu,nu) and not mu==nu]
    return Poset([E,R])

def Yakob_Triangle_Mots(n,k=0):
    E=[part_to_word(nu) for nu in PartagesTriangulaires(n) if Is_Dominant(nu) and nu.length()>=k]
    R=[(mu,nu) for mu in E for nu in E  
       if OrdreSchur(word_to_part(mu),word_to_part(nu)) and not mu==nu]
    return Poset([E,R])

@cached_function
def YoungTrianglesBaseHauteur(a,b):
    E=[part_to_word(nu) for nu in PartagesTriangles_BaseHauteur(a,b) if Is_Dominant(nu)]
    R=[(mu,nu) for mu in E for nu in E  
       if word_to_part(nu).contains(word_to_part(mu))]
    return Poset([E,R])

@cached_function
def Young_Triangle(n):
    E=[nu for k in range(n+1) for nu in TriangularPartitions(k)]
    R=[(mu,nu) for k in range(1,n+1) for mu in TriangularPartitions(k-1) 
       for nu in TriangularPartitions(k) if nu.contains(mu)]
    return Poset([E,R],cover_relations=True)


def Young_Triangle_Log(n,Test=lambda tau:False,Couleur='orange'):
    Sommets=add(point(CentreLog(tau),size=25/sqrt(max(tau.size(),1)),color='green',zorder=1) for k in range(n+1) 
                for tau in TriangularPartitions(k) if not Test(tau))
    Speciaux=add(point(CentreLog(tau),size=25/sqrt(max(tau.size(),1)),color=Couleur,zorder=1) for k in range(n+1) 
                for tau in TriangularPartitions(k) if Test(tau))
    Aretes=add(line([CentreLog(nu), CentreLog(tau)],thickness=.4,color='lightblue',zorder=0) 
               for k in range(1,n+1) for nu in TriangularPartitions(k-1) 
       for tau in TriangularPartitions(k) if tau.contains(nu))
    p=(Aretes+Sommets+Speciaux)
    p.axes(False)
    p.set_aspect_ratio(1)
    return p


def Young_Triangle_Quadran(n,Test=lambda tau:False,Couleur='red'):
    Sommets=add(point(Centre(tau),size=3,color='green',zorder=1) for k in range(n+1) 
                for tau in TriangularPartitions(k) if not Test(tau))
    Speciaux=add(point(Centre(tau),size=3,color=Couleur,zorder=1) for k in range(n+1) 
                for tau in TriangularPartitions(k) if Test(tau))
    Aretes=add(line([Centre(nu), Centre(tau)],thickness=.2,color='black',zorder=0) 
               for k in range(1,n+1) for nu in TriangularPartitions(k-1) 
       for tau in TriangularPartitions(k) if tau.contains(nu))
    p=(Aretes+Sommets+Speciaux)
    p.axes(False)
    p.set_aspect_ratio(1)
    return p



@cached_function
def Weight(mu,v):
    return add(Omega(tab)*poids(tab,v) for tab in StandardTableaux(mu).list())

def T_mu(mu):
    return mul(q**i*t**j for i,j in mu.cells())

def poids(tab,v):
    mu=tab.shape()
    return mul((q**i*t**j)**(v[tab.entry((i,j))-1]) for i,j in mu.cells())

def Omega(tab):
    mu=tab.shape()
    C_mu=mu.cells()
    return (mul((star(q**a*t**b-q**i*t**j)*star(q**a*t**b-q**(i+1)*t**(j+1))
                 /(star(q**a*t**b-q**(i+1)*t**j)*star(q**a*t**b-q**i*t**(j+1))))
                for i,j in C_mu for a,b in C_mu
                if tab.entry((i,j))<tab.entry((a,b)))
            *mul(q**a*t**b/(q**a*t**b-q**(i+1)*t**(j+1))
                 for i,j in C_mu for a,b in C_mu
                 if tab.entry((i,j))+1==tab.entry((a,b)))
            /mul(q**i*t**j-1
                 for i,j in C_mu if not (i,j)==(0,0)))

def star(x):
    if x==0: return 1
    else: return x




