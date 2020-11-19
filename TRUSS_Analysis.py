import numpy as np
import pandas as pd
nodes=int(input("Enter number of nodes in a given truss : "))
element=int(input("Enter number of elements in a given truss : "))
sn=[]
en=[]
xcoord=[]
ycoord=[]
cos=[]
sin=[]
le=[] 
M=[]
g=[]
F=np.zeros((2*nodes,1),dtype=float)
NWHS=[]
NWRS=[]
E1=[]
GSM = np.zeros((nodes*2, nodes*2), dtype=float)
def stiffmat(A,E,L,Mat):
    C=(A*E*10**3)/L
    T=np.multiply(C,Mat)
    return T
for i in range(nodes):
    x=float(input("Enter x coordinate of node %d : "%(i+1)))
    y=float(input("Enter y coordinate of node %d : "% (i+1)))
    xcoord.append(x)
    ycoord.append(y)
for j in range (element):
    snfe=int(input("Enter starting node number for element %d : "%(j+1)))
    enfe=int(input("Enter ending node number for element %d : "%(j+1)))
    sn.append(snfe)
    en.append(enfe)
    a=xcoord[snfe-1]
    b=ycoord[snfe-1]
    c=xcoord[enfe-1]
    d=ycoord[enfe-1]
    length=np.sqrt(np.square(c-a)+np.square(d-b))
    le.append(length)
    l=(c-a)/length
    cos.append(l)
    m=(d-b)/length
    sin.append(m)
    mat=np.array([[l**2,l*m,-l**2,-l*m],
                 [l*m,m**2,-l*m,-m**2],
                 [-l**2,-l*m,l**2,l*m],
                 [-l*m,-m**2,l*m,m**2]])
    M.append(mat)
print("\n")
print("Enter \'d\' if each element has different c/s area")
print("Enter \'i\' if each element has identical c/s area.")
type=str(input("Enter condition for c/s area of elements : "))
if type in ['D','d']:
    for k in range(element):
        A=float(input("\nEnter c/s area in mm^2 for element %d : "%(k+1)))
        E=float(input("Enter Young's Modulus in GPa for element %d : "%(k+1)))
        esm=stiffmat(A,E,le[k],M[k])
        print("\n k%d = \n"%(k+1))
        print(esm)
        g.append(esm)
        E1.append(E)
elif type in ['I','i']:
    A = float(input("Enter c/s area in mm^2 for elements : "))
    print("\n")
    print("Enter \'i\' if Young's Modulus is identical for all elements." )
    print("Enter \'d\' if Young's Modulus is different for all elements.")
    etype=str(input("Enter condition for Young's Modulus : "))
    if etype in ['D','d']:
        for k in range(element):
            E = float(input("Enter Young's Modulus in GPa for element %d : "%(k+1)))
            esm = stiffmat(A, E, le[k], M[k])
            print("\n k%d = \n"%(k+1))
            print(esm)
            print("\n")
            g.append(esm)
            E1.append(E)
    elif etype in ['I','i']:
        E = float(input("Enter Young's Modulus in GPa for elements  : "))
        for h in range(element):
            esm = stiffmat(A, E, le[h], M[h])
            print("\n k%d = \n"%(h+1))
            print(esm)
            g.append(esm)
            E1.append(E)
for l in range(element):
    gstmat=g[l]
    h = sn[l]*2-2
    j = en[l]*2-2
    h1 = (sn[l]*2+1)-2
    j1 = (en[l]*2+1)-2
    ele=[h,h1,j,j1]
    for m in range(4):
        for n in range(4):
            x=ele[m]
            y=ele[n]
            GSM[x,y]+=gstmat[m,n]
print("\n")
print("-----------------------Global stiffness Matrix--------------------------")
print("\n K= \n")
print(np.round(GSM*10**-3,3))
print("\n")
nnf=int(input("Enter total number of nodes subjected to external loads : "))
for p in range(nnf):
    nno=int(input("Enter node number subjected to loading : "))
    fh = float(input("Enter horizontal load in N acting at node %d : "%nno))
    fv = float(input("Enter vertical load in N acting at node %d : "%nno))
    F[2*nno-2,0]=fh 
    F[(2*nno+1)-2,0]=fv 
force = pd.DataFrame(data=F, index=range(2*nodes), columns=range(1))
gsd=pd.DataFrame(data=GSM,index=range(2*nodes),columns=range(2*nodes))
nnws=int(input("\nEnter number of nodes having supports: "))
for o in range(nnws):
    print("\n")
    nws=int(input("Enter node number having support : "))
    print("Enter h for hinge support.")
    print("Enter r for roller support.")
    st=input("Enter condition for support : ")
    if st in ['h','H']:
        NWHS.append(nws)
        K=gsd.drop([nws*2-2,(nws*2+1)-2],axis=0)
        K=K.drop([nws*2-2,(nws*2+1)-2],axis=1)
        gsd=K
        ff=force.drop(2*nws-2,axis=0)
        ff=ff.drop((2*nws+1)-2,axis=0)
        force=ff
    elif st in ['r','R']:
        NWRS.append(nws)
        K=gsd.drop((nws*2+1)-2,axis=0)
        K=K.drop((nws*2+1)-2,axis=1)
        gsd=K
        ff=force.drop((2*nws+1)-2, axis=0)
        force=ff
redmat=gsd.to_numpy()
redfmat=force.to_numpy()
X=np.linalg.solve(redmat,redfmat)
print("\n-------------------------Nodal Displacement (*10^-3 mm)---------------------------")
disp=np.ones((2*nodes))
for v in NWHS:
    disp[2*v-2]=0
    disp[(2*v+1)-2]=0
for w in NWRS:
    disp[(2*w+1)-2]=0
x=0
for z in range(2*nodes):
    if disp[z]==1:
        disp[z]=X[x]
        x+=1
disround=np.around(disp*10**3,6)
a=disround.tolist()
u=a[0:2*nodes:2]
v=a[1:2*nodes:2]
data={"Nodes":range(1,nodes+1),"dx":u,"dy":v}
u1=pd.DataFrame(data,columns=['Nodes','dx','dy'])
print(u1.to_string(index=False))
print("\n")
print("-------------------------Elemental Stresses (MPa)---------------------------")
D=np.array(disround)
S=[]
for y in range(element):
    C=E1[y]/le[y]
    h=sn[y]*2-2
    j=en[y]*2-2
    h1=(sn[y]*2+1)-2
    j1=(en[y]*2+1)-2
    B=np.array([-cos[y],-sin[y],cos[y],sin[y]])
    q=np.array([[D[h]],[D[h1]],[D[j]],[D[j1]]])
    multi=np.matmul(B,q)
    stress=C*multi[0]
    S.append(stress)
elestess={"Elements":range(1,element+1),"Stresses":S}
Eles=pd.DataFrame(elestess,columns=["Elements","Stresses"])
print(Eles.to_string(index=False))   
print("-------------------------Reaction at Supports (N)---------------------------")
hr=[]
vr=[]
for s in NWHS:
    hori=np.matmul(GSM[2*s-2],disp)-F[2*s-2]
    verti=np.matmul(GSM[(2*s+1)-2],disp)-F[(2*s+1)-2]
    hr.append(hori[0])
    vr.append(verti[0])
for t in NWRS:
    verti=np.matmul(GSM[(2*t+1)-2],disp)-F[(2*t+1)-2]
    hr.append(0)
    vr.append(verti[0])
L3=NWHS+NWRS
reactions={'Nodes':L3,'Rx':hr,'Ry':vr}
R=pd.DataFrame(reactions,columns=["Nodes","Rx","Ry"])
print(R.to_string(index=False))