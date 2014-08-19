def gclust.rsvt(Xorigin,r=1,maxsvt=10,nmfout=False,maxit=10000):

    import numpy as np
    import nimfa
    Xorigin=np.matrix(np.loadtxt("/cis/home/btang/Xbtsensordata.txt",dtype=float))
    Xraw=Xorigin
    Xorigin=Xorigin*np.linalg.inv(np.diag(np.array(Xorigin.sum(axis=0))[0,:]))    

    u,s,v=np.linalg.svd(Xorigin)
    U=u[:,0:r]
    V=v[:,0:r]
    S=np.diag(s[0:r])
    if maxsvt<0:
       maxsvt=1 
    for itr in range(1,maxsvt):
        X=U*S*(V.T)
        X[X<0]=0
        u,s,v=np.linalg.svd(X) # try catch
        UU=u[:,0:r]
        VV=v[:,0:r]
        SS=np.diag(s[0:r]) 
        if np.linalg.norm(UU-U)+np.linalg.norm(VV-V)+np.linalg.norm(SS-S)<1e-12:
            break  
        U=UU
        V=VV
        S=SS   
    XX=X*np.linalg.inv(np.diag(np.array(X.sum(axis=0))[0,:]))   
    #a=np.array(X.sum(axis=0))[0,:]
    #a=np.array(Xcolsum)[0,:]
    #XX1=X*(np.linalg.inv(np.diag(a)))
    #XX＝X*(np.linalg.inv(np.diag((np.array(X.sum(axis=0)))[0,:])))
    #stash1=np.argwhere(XX.sum(axis=1)==0) 
    stash=np.argwhere(XX.sum(axis=1)!=0) 
    Nrow_XX=XX.shape[0]
    if stash[:,0][:,0].size<Nrow_XX:
        XX=XX[stash[:,0][:,0],:]
        #XX1=np.delete(XX,stash1[:,0][:,0],axis=0) 
    
    
    mynmf=nimfa.mf(XX,method='nmf',rank=r)
    fctr_res = nimfa.mf_run(mynmf)
    WW=np.zeros((X.shape[0],r))
    if stash[:,0][:,0].size<Nrow_XX:
        WW[stash[:,0][:,0],:]=fctr_res.fit.basis()
    else:  
        WW=fctr_res.fit.basis()
        
    HH=fctr_res.fit.coef()  
    if r>1:
        HH=np.diag(np.array(WW.sum(axis=0))[0,:])*HH
        WW=np.around(WW*np.linalg.inv(np.diag(np.array(WW.sum(axis=0))[0,:])),12)
    else:
        HH= WW.sum(axis=0)*HH
        WW=np.around(WW/float(WW.sum(axis=0)),12)  
    
    HH=np.around(HH*np.linalg.inv(np.diag(np.array(HH.sum(axis=0))[0,:])),12)
    Xmean=X*np.diag(np.array(Xraw.sum(axis=0))[0,:]) #why here? not used in getAICc 
    
    if nmfout:
        nmf=mynmf
    else:
        nmf= None
    
    W=WW
    H=HH
    Xorigin=Xraw
    return(nmf,W,H,Xorigin,Xmean)
    
     
