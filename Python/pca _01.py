from numpy import *
from matplotlib.pyplot import *

#Xp=loadtxt("Xp.csv",skiprows=1,delimiter=",")
X=loadtxt("Xs.csv",skiprows=1,delimiter=",")
Xd=loadtxt("Xd.csv",skiprows=1,delimiter=",") #zmerena data
Xw=loadtxt("Xw.csv",skiprows=1,delimiter=",") #smerodatne odchylky namerenych dat

#X=random.randn(1000,30)
#X[500:,:]=X[500:,:]+3


def fPCA2D(X):   # neig...number of principal components
    n=X.shape[1]
    N=X.shape[0]
    #X=(X-mean(X))/std(X)/N
    X=(X-Xd)**2/Xw/N
    covX=dot(X.T,X)/N
    d,V=linalg.eig(covX)
    d=abs(d)
    Vp=zeros((n,2), dtype='complex_')
    for i in range(2):
        Vp[:,i]=V[:,d==max(d)].squeeze()
        d[d==max(d)]=-1
    C=dot(X,Vp)
    return(C,d)

C,d=fPCA2D(X)

figure()
subplot(211)
plot(abs(d),'.');grid()
subplot(212)
plot(C[:,0],C[:,1],'.k');grid()
show()

