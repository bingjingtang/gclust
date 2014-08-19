def getAICc(nmf,W,H,Xorigin,Xmean):
    #Xorigin=np.matrix(np.loadtxt("/Users/tangbingjing/Downloads/Xorigin1.txt",dtype=int))
    #W=np.matrix(np.loadtxt("/Users/tangbingjing/Downloads/W.txt",dtype=float))
    #H=np.matrix(np.loadtxt("/Users/tangbingjing/Downloads/H.txt",dtype=float))    
    import numpy as np
    Xmean=W*H*np.matrix(np.diag(np.array(Xorigin.sum(axis=0))[0,:]))
    Xmean[Xmean<1e-12]=0
    Phat=W*H
    Phat[Phat<1e-12]=0
    nparams=(1.0*(W>0)).sum(axis=0)
    if W.shape[1]>1:
        nparams=nparams*np.diag(np.array(1.0/(Xorigin*H.T).sum(axis=0))[0,:])
    else:
        nparams=nparams/float((Xorigin*H.T).sum(axis=0))
    
    nparams=nparams.sum()
    
    Phat=np.hstack(Phat)
    Xorigin=np.hstack(Xorigin)
    zeroprob=np.argwhere(Phat<1e-12)
    zeroprob=zeroprob[:,0][:,1]
    zeroprob2=np.argwhere(Phat>1e-12)
    zeroprob2=zeroprob2[:,0][:,1]
    if any(tuple(np.array(Xorigin[0,zeroprob]>1e-12)[0,:])):
        retval=float('inf')
    elif zeroprob.size>0:
        Phat=Phat[0,zeroprob2]
    retval=2*(1)*(nparams)-2*np.sum(np.multiply(Phat,np.log(Phat)))
    nclust=W.shape[1]
    negloglikpart=-2*np.sum(np.multiply(Phat,np.log(Phat)))
    parampart=2*nparams*(1)
    AIC=retval
    return(nclust,negloglikpart,parampart,AIC)
    
    #zeroprob=np.where(Phat<1e-12)
    
    #if any(tuple(np.array(Xorigin[zeroprob]>1e-12)[0,:])):
        #retval=float('inf')
    #elif not(any(x.size==0 for x in zeroprob)):
        
        
        
