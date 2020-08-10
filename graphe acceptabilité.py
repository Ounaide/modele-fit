from matplotlib.pyplot import *
from numpy import *
from easygui import choicebox
from fractions import Fraction
from scipy.optimize import curve_fit

Rxl=0
xmin=0

Rv=[-185.35,-261.95,-269.55,-422.05,-282.95,-257.425,-486.95,-354.2,-322.675,-417.45,-707.825,-750.6,-173.8,-216.6,-155.675,-179.975,-347.875,-379.975]
Rt=[460,514.6,347.45,432.35,403.825,417.275,956.725,504.125,568.1,639.375,607.275,730.725,301.625,318.775,353.6,343.55,366.375,388.375]
Rq=[27.6,28.925,21,22.925,25.6,27.375,59.15,32.425,33.025,45.675,77.45,80.15,17.825,20.725,17.075,20.05,29.4,23.25]
Acc=[0.833,0.833,0.5,1,0.833,0.667,0,0.167,0,0.167,0,0.167,0.833,1,1,0.667,0.833,0.833]
Rdqy=[0.0022225,0.0023625,0.0019175,0.002145,0.0021475,0.00212,0.00422,0.00244,0.002615,0.003155,0.0045125,0.004415,0.0016475,0.0017675,0.0015975,0.001715,0.0021225,0.0018675]
Rdqx=[0.0040125,0.00593,0.0043025,0.00561,0.0044675,0.005325,0.0072775,0.003885,0.0047075,0.00543,0.0064675,0.007365,0.0025625,0.0026075,0.0024775,0.00264,0.0028225,0.002905]


g=globals()

global Rqmax
Rqmax=max(Rq)

class Modeles():
    def __init__(self):      
        """Rq=array(Rq)
        Rt=array(Rt)
        Rv=array(Rv)"""
        pass
        
    def expo(self,x,k):
        global y0
        y0=exp(k-x)
        global fstr0
        fstr0="exp(k-x)"
        return y0
    
    def Uexponentiel(self,x,k):
        global y1
        y1=heaviside(x+1,0) + heaviside(x-k,0)*(exp((k-x)*10/Rqmax)-1)
        global fstr1
        fstr1="U(x+1,0) + U(x-k,0)*(exp((k-x)*10/max(x))-1"
        return y1
    
    def sigmoide(self,x,l,k):
        y2=(-1/(1+exp(-l*(x-k))))+1
        global fstr2
        fstr2="(-1/(1+exp(-λ*(Rx-Rx^))))+1"
        return y2
    
    def polynome(self,x,a,b,c):

        y3=a*x**2 + b*x + c
        global fstr3
        fstr3="(a*x²)+(b*x)+c"
        return y3

Md = Modeles()

def on_close(event):
    clf()
    fit()

def fit():
        
    sRx=choicebox('Choisissez la variable','zok bouk',['Rq','Rt','Rv','Rdqx','Rdqy'])
    Rx=g[sRx]
    pars,cov=curve_fit(Md.sigmoide,Rx,Acc)
    #r=round(corrcoef(Acc,yx(Rq,pars[0]))[0][1],2)
    
    r=round(corrcoef(Acc,Md.sigmoide(Rx,*pars))[0][1],2)
    fig=figure()
    fig.suptitle(f"Modèle :{fstr2}" +" ("+ Md.sigmoide.__name__+ ")"+ f"\n R²={r}")
    scatter(Rx,Acc,label='data') #points bruts
    Rx.sort()
    y=Md.sigmoide(Rx,*pars)
    #y=yx(Rq,pars[0]) #pour expo
    x2=linspace(min(Rx),max(Rx),100*len(Rx))
    y2=Md.sigmoide(x2,*pars)
    #y2=Md.sigmoidenentiel(x2,pars[0])

    scatter(Rx,y,marker='+',color='r',label=f'fit: λ=%.02f {sRx}^=%.02e'%tuple(pars)) #sigmoide
    #scatter(Rt,y,marker='+',color='r',label='fit: Rx^=%.02f'%pars[0]) #exponentielle
    #scatter(Rq,y,marker='+',color='r',label='fit: a=%.02f b=%.02f c=%.02f'%tuple(pars)) #polynome
    plot(x2,y2,color='r',linewidth=0.3)#courbe modèle
    
    xlabel(f'{sRx} (µm)')
    ylabel('Acc')
    legend()

    
    fig.canvas.mpl_connect('close_event',on_close)
    show()
    
    return pars
                       
                       
def Psolve(x1,y1,x2,y2,x3,y3):
    B=array([y1,y2,y3]).reshape(3,1)
    A=array([[x1**2,x1,1],[x2**2,x2,1],[x3**2,x3,1]])  
    resu=linalg.solve(A,B)
    a,b,c=resu.squeeze()
    a,b,c= Fraction(a).limit_denominator(),Fraction(b).limit_denominator(),Fraction(c).limit_denominator()
    f=figure(f"a={a} b={b} c={c}")
    f.suptitle("Modèle: ax²+bx+c")
    #x=linspace(0.9*min(x1,x2,x3),1.1*max(x1,x2,x3),1000)
    x=linspace(min(Rt),max(Rt),100*len(Rt))
    y=((a*x**2)+(b*x)+c)
    
    plot(x,y)
    plot(x1,y1,'ro')
    plot(x2,y2,'ro')
    plot(x3,y3,'ro')
    scatter(Rt,Acc)
    show()
    
    return a,b,c

def hv():
    Rx=choicebox('yes','aaa',['Rq','Rt'])
    if Rx in ['Rt','rt','RT','rT']:
        xmax=1000
        Rxl=500
    elif Rx in ['Rq','rq','RQ','rQ']:
        xmax=150
        Rxl=30
    nbr=xmax*10
    x=linspace(xmin,xmax,nbr)
    global Acc
    Acc=heaviside(x+1,0) + heaviside(x-Rxl,0)*(exp((Rxl-x)*10/xmax)-1)

    plot(x,Acc)
    axvline(x=Rxl,color='r',ymax=1)
    show()

def sg(l):
    k=30
    f=figure()
    f.suptitle(f'λ={l}')
    x=linspace(-10,10,100)
    y=(-1/(1+exp(-l*(x-k))))+1 #forme exacte
    #z=exp(-1.5)*exp(l*x) #forme approximée
    plot(x,y)
    #plot(x,z)
    show()
    


    