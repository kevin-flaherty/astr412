#An attempt to apply machine learning to turbulence results. Start with generating many spectra for different levels of turbulence, and using machine learning to pull out features of the spectra associated with turbulence.


def make_models():
    '''Make a series of models with different levels of turbulence. Write the flux as a function of wavelength into a file, which can then be read in later.'''

    from astropy.io import fits
    from helper import im_plot_spec
    import numpy as np
    #in order to generate models, you will need disk.py and raytrace.py from githud.com/kevin-flaherty/disk_model
    import disk
    import raytrace



    nmodels = 10000
    vturb = np.linspace(0,3.,nmodels)
    q = np.random.rand(nmodels)*(-.1+.8)-.8
    rc = np.random.rand(nmodels)*(2.5-2.0)+2.0
    tmid = np.random.rand(nmodels)*(20-10.)+10
    tatm = np.random.rand(nmodels)*(35.-20)+20
    incl = np.random.rand(nmodels)*(-25+45)-45 
    abund = np.random.rand(nmodels)*(-4+6)-6
    #vsys = np.random.rand(nmodels)*(2.92-2.9)+2.9
    
    
    datfile = '/Volumes/disks/kevin/DMTau_short/dmtau_co21sb'
    hdr = fits.getheader(datfile+'.vis.fits')
    freq=(np.arange(hdr['naxis4'])+1-hdr['crpix4'])*hdr['cdelt4']+hdr['crval4']
    obsv = (hdr['restfreq']-freq)/hdr['restfreq']*2.99e5
    chanstep = (np.abs(obsv[1]-obsv[0]))#-0.208337
    vsys=6.063
    nchans = int(2*np.ceil(np.abs(obsv-vsys).max()/chanstep)+1)
    chanmin = -(nchans/2.-.5)*chanstep
    offs = [-.075,-.017]#[.30,0.00]
    resolution = 0.05
    obs = [150,101,300,150] #rt grid nr,nphi,nz,zmax (250,101,250,150)
    

    file = open('vturb_spec3.txt','w')

    for i in range(nmodels):
        #(-.638,1.8,turb,24.73,25.71,1.)
        params = [q[i],
                  0.04, #Mdisk
                  1., #pp
                  1., #Rin
                  1000., #Rout
                  10**(rc[i]), #Rc
                  incl[i], #incl
                  .54, #Mstar
                  10**(abund[i]), #Xco
                  vturb[i], 
                  70., #Zq0
                  tmid[i], 
                  tatm[i], 
                  [.79,1e8],
                  [1., 800.],
                  0,
                  -1]
        disk_structure = disk.Disk(params,rtg=False)
        modfile='alma'

        disk_structure.set_obs(obs)
        disk_structure.set_rt_grid()
        disk_structure.set_line('co21')
        disk_structure.add_mol_ring(params[3],params[4],.79,3.,params[8],just_frozen=True)
        print i,params
        raytrace.total_model(disk=disk_structure,chanmin=chanmin,nchans=nchans,chanstep=chanstep,offs=offs,modfile=modfile,imres=resolution,obsv=obsv,vsys=vsys,freq0=230.538,Jnum=1,PA=154.8,distance=145.1)


        xaxis,spec=im_plot_spec('alma.fits',color='k',ls='-')

        np.append(vturb[i],spec).tofile(file,sep=',')
        file.write('\n')


    file.close()
    

def make_models2():
    '''Similar to make_models, but now use the image of the central velocity channel, rather than the spectra.'''

    from astropy.io import fits
    import numpy as np
#in order to generate models, you will need disk.py and raytrace.py from githud.com/kevin-flaherty/disk_model
    import disk
    import raytrace


    nmodels = 200
    vturb = np.linspace(0,3.,nmodels)
    q = np.random.rand(nmodels)*(-.1+.8)-.8
    rc = np.random.rand(nmodels)*(1.9-1.7)+1.7
    tmid = np.random.rand(nmodels)*(20-15.)+15
    tatm = np.random.rand(nmodels)*(35.-25)+25
    incl = np.random.rand(nmodels)*(-25+45)-45
    abund = np.random.rand(nmodels)*(-4+6)-6
    #vsys = np.random.rand(nmodels)*(2.92-2.9)+2.9
    
    datfile = '/Volumes/disks/kevin/v4046_short/v4046_co21sb'
    hdr = fits.getheader(datfile+'.vis.fits')
    freq=(np.arange(hdr['naxis4'])+1-hdr['crpix4'])*hdr['cdelt4']+hdr['crval4']
    obsv = (hdr['restfreq']-freq)/hdr['restfreq']*2.99e5
    chanstep = (np.abs(obsv[1]-obsv[0]))#-0.208337
    vsys=2.91
    obsv = [vsys-chanstep,vsys,vsys+chanstep]
    nchans = 3#int(2*np.ceil(np.abs(obsv-vsys).max()/chanstep)+1)
    chanmin = -(nchans/2.-.5)*chanstep
    offs = [.2721,.017]#[.30,0.00]
    resolution = 0.05
    obs = [200,101,250,150] #rt grid nr,nphi,nz,zmax (250,101,250,150)
    

    file = open('vturb_image2.txt','w')

    for i in range(nmodels):
        #(-.638,1.8,turb,24.73,25.71,1.)
        params = [q[i],
                  0.09, #Mdisk
                  1., #pp
                  1., #Rin
                  800., #Rout
                  10**(rc[i]), #Rc
                  incl[i], #incl
                  1.75, #Mstar
                  10**(abund[i]), #Xco
                  vturb[i], 
                  35., #Zq0
                  tmid[i], 
                  tatm[i], 
                  [.79,1e8],
                  [1., 800.],
                  -1]
        disk_structure = disk.Disk(params,rtg=False)
        modfile='alma'

        disk_structure.set_obs(obs)
        disk_structure.set_rt_grid()
        disk_structure.set_line('co21')
        disk_structure.add_mol_ring(params[3],params[4],.79,3.,params[8],just_frozen=True)
        print i,params
        raytrace.total_model(disk=disk_structure,chanmin=chanmin,nchans=nchans,chanstep=chanstep,offs=offs,modfile=modfile,imres=resolution,obsv=obsv,vsys=vsys,freq0=230.538,Jnum=1,PA=75.65,distance=73.)

        im=fits.open('alma.fits')
        data = im[0].data[1,125:385,175:325] #260x150
        np.append(vturb[i],data.flatten()).tofile(file,sep=',')
        file.write('\n')
        im.close()

    file.close()
        


def read_models():
    '''Read in the model file, splitting apart the 'y' array (containing the turbulence levels of each model) from the 'x' array (containing the flux in each spectra).'''

    import numpy as np

    data = np.recfromtxt('vturb_spec2.txt',delimiter=',')

    vturb = data[:,0]/3.438
    spec = data[:,1:]

    return vturb,spec
        

def pca_spec(n_components=6,degree=4):
    '''Read in the spectra and perform PCA to reduce the dimensionality of the data. The result that then can be regressed against the turbulence to find out how the spetrum varies with turbulence.'''
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LinearRegression
    from sklearn.neighbors import KernelDensity
    #from sklearn.grid_search import GridSearchCV
    #from sklearn.cross_validation import LeaveOneOut
    import matplotlib.pyplot as plt
    import numpy as np
    from random import sample

    vturb,spec = read_models()

    #Divide into training and validation sets
    nmodels = len(vturb)
    val_list = np.array(sample(range(nmodels),int(.1*nmodels)))
    vturb_train = np.zeros(int(.9*nmodels))
    spec_train = np.zeros((int(.9*nmodels),spec.shape[1]))
    vturb_val = np.zeros(int(.1*nmodels))
    spec_val = np.zeros((int(.1*nmodels),spec.shape[1]))
    val_count = 0
    train_count = 0
    for i in range(nmodels):
        if (i==val_list).sum()>0:
            #part of validation sample
            vturb_val[val_count] = vturb[i]
            spec_val[val_count,:] = spec[i,:]
            val_count+=1
        else:
            vturb_train[train_count] = vturb[i]
            spec_train[train_count,:] = spec[i,:]
            train_count+=1
            


    #Use PCA with 2 components. 
    #First plot the values of the two components, with turbulence highlighted
    #n_components = 8
    pca = PCA(n_components)
    projected = pca.fit_transform(spec_train)
    



    #Plot the components to see what they look like
    plt.figure()
    plt.subplot(221)
    for i in range(n_components):
        plt.plot(pca.components_[i],label='component {:0.0f}'.format(i),lw=n_components-i)
    #plt.legend(frameon=False,loc='best')

    print np.cumsum(pca.explained_variance_ratio_)


    #Use regression techniques to correlate the components with turbulence.
    #Start with Linear Regression
    model = LinearRegression(fit_intercept=True)
    model.fit(projected,vturb_train)

    #Try a polynomial transformation before the Linear Regression
    from sklearn.preprocessing import PolynomialFeatures
    poly = PolynomialFeatures(degree=degree,include_bias=False)
    projected2 = poly.fit_transform(projected)

    model2 = LinearRegression(fit_intercept=True)
    model2.fit(projected2,vturb_train)
    print 'R^2 of the fit: ',model2.score(projected2,vturb_train)
    
    
    #How well can we recover the input turbulence?
    vturb_predicted = model.predict(projected) #linear model
    vturb_predicted2 = model2.predict(projected2) #polynomial model

    plt.subplot(223)
    #plt.plot(vturb,vturb_predicted,'or',ms=1)
    plt.plot(vturb_train,vturb_predicted2,'ok',ms=1)
    plt.plot([0,np.max(np.append(vturb,vturb_predicted2))],[0,np.max(np.append(vturb,vturb_predicted2))],color='c',ls='--',lw=2)
    plt.xlabel('Actual v$_{turb}$')
    plt.ylabel('Predicted v$_{turb}$')


    #Run against the validation set and see the accuracy of the result
    #vturb_val_predicted = np.array(len(vturb_train))
#    for i in len(vturb_train):
    pspec = pca.transform(spec_val)
    pspec2 = poly.fit_transform(pspec)
    vturb_val_predicted = model2.predict(pspec2)
    plt.plot(vturb_val,vturb_val_predicted,'or',ms=1)
    error = (vturb_val-vturb_val_predicted)/vturb_val
    print 'Median fractional error in validation set: ',np.median(np.abs(error))

#n_components=2: .34 (degree=1),  .32 (degree=2), .33 (degree=3), .31 (degree=4), .32 (degree=5)
#n_components=3: .27 (degree=1), .25, .24, .22, .20
#n_components=4: .27, .24, .22, .23, .21
#n_components=5: .26, .25, .23, .20, .18, .19, .24 
#n_components=6: .27, .22, .18, .15*, .15, .19, .55
#n_components=7: .26, .22, .15, .14, .15, .37



    #Predict turbulence in our three disks from this decomposition
    disk_spec = np.array([0.07051709,  0.08370663,  0.09472561,  0.14514502,  0.19218614, 0.26022604,  0.3638115 ,  0.51294306,  0.77375685,  1.23985817, 1.86580414,  2.92085159,  4.45566949,  7.49341501,  9.57812521, 8.47671087,  7.11170747,  7.57635818,  9.61814015,  8.07201153, 5.09846799,  3.1666609 ,  1.96394516,  1.3104509 ,  0.82953736, 0.56494169,  0.38539208,  0.25932753,  0.18477959,  0.13578935, 0.10984786])
    vturb_actual = .275
    disk_spec = disk_spec[:,np.newaxis].T
    #disk_spec_pca = pca.transform(disk_spec)
    #disk_spec_pca2 = poly.fit_transform(disk_spec_pca)
    vturb_disk_predicted = []
    for i in range(100):
        noise = np.random.randn(len(disk_spec))*.01*disk_spec
        disk_spec_pca = pca.transform(disk_spec+noise)
    #    vturb_disk_predicted.append(model.predict(disk_spec_pca))
        disk_spec_pca2 = poly.fit_transform(disk_spec_pca)
        vturb_disk_predicted.append(model2.predict(disk_spec_pca2))
    vturb_disk_predicted = np.array(vturb_disk_predicted)
    #vturb_disk_predicted = .2#model2.predict(disk_spec_pca2)
    #20% higher/lower for systematic uncertainty
    disk_spec_pcaH = pca.transform(disk_spec*1.2)
    disk_spec_pcaL = pca.transform(disk_spec*0.8)
    disk_spec_pcaH2 = poly.fit_transform(disk_spec_pcaH)
    disk_spec_pcaL2 = poly.fit_transform(disk_spec_pcaL)
    vturb_disk_predictedH = model2.predict(disk_spec_pcaH2)
    vturb_disk_predictedL = model2.predict(disk_spec_pcaL2)
    #plt.subplot(236)
    #plt.plot([.01,],[vturb_disk_predicted,],'sb')
    plt.plot(np.zeros(100)+vturb_actual,vturb_disk_predicted,'sb',ms=2)
    plt.plot([vturb_actual,],[vturb_disk_predictedH,],'sc',ms=5)
    plt.plot([vturb_actual,],[vturb_disk_predictedL,],'sc',ms=5)
    #plt.axhline(.9*vturb_disk_predicted,color='k',ls=':')
    #plt.axhline(1.1*vturb_disk_predicted,color='k',ls=':')

    print 'Disk actual, predicted: ',vturb_actual,np.mean(vturb_disk_predicted),np.std(vturb_disk_predicted)
    print 'Median % error: ',np.median(np.abs(vturb_train-vturb_predicted2)/vturb_train)

    #plot the derived and original spectral to see if they match
    plt.subplot(222)
    plt.plot(disk_spec[0,:],color='k',lw=2)
    plt.plot((pca.inverse_transform(disk_spec_pca))[0,:],color='r',lw=2,ls='--')
    plt.legend(('Original','Reconstructed'),frameon=False,loc='best')



    #Derive the probability distribution of actual values of turbulence based on the predicted level of turbulence
    w=(vturb_predicted2>vturb_disk_predicted.min()) & (vturb_predicted2<vturb_disk_predicted.max())
    plt.subplot(224)
    #bandwidths = 10**(np.linspace(-2,0,100))
    #grid = GridSearchCV(KernelDensity(kernel='gaussian'),{'bandwidth':bandwidths},cv=LeaveOneOut(len(vturb_train[w])))
    #grid.fit(vturb_train[w][:,None])
    #print grid.best_params_

    kde = KernelDensity(bandwidth=.03,kernel='gaussian')
    kde.fit(vturb_train[w][:,None])

    v_d = np.linspace(0,1,100)
    logprob = kde.score_samples(v_d[:,None])

    plt.plot(v_d,np.exp(logprob),color='k')

    plt.axvline(vturb_actual,color='k',ls=':')
    plt.axvline(np.mean(vturb_disk_predicted),color='b',ls='--')
    plt.axvline(vturb_disk_predictedH,color='c',ls='--')
    plt.axvline(vturb_disk_predictedL,color='c',ls='--')
    plt.xlim(0,1)
    plt.ylabel('PDF')
    plt.xlabel('Actual v$_{turb}$')


    #Using three components with polynomial regression does OK, but not great. There is still a lot of scatter, although there is a general trend. It vastly overpredicts the turbulence in MWC 480...(PCA+LR_spec.png)


def pca_image():
    '''Similar to pca_spec, but using images instead of the spectra. This will fit principal components, and display the principal components, and then see if it can correlate the coefficients of those components with turbulence.'''

    from sklearn.decomposition import PCA
    from sklearn.linear_model import LinearRegression
    import matplotlib.pyplot as plt
    import numpy as np

    vturb,images = read_models()
    image_mean = np.mean(images,axis=0)
    
    #Use PCA with 4 components
    pca = PCA(4)
    projected = pca.fit_transform(images)

    #Plot the turbulence as a function of these components
    plt.figure()
    plt.subplot(241)
    plt.scatter(projected[:,0],projected[:,1],c=vturb,edgecolor='none',cmap=plt.cm.get_cmap('afmhot'))
    plt.xlabel('component 1')
    plt.ylabel('component 2')
    plt.subplot(242)
    plt.scatter(projected[:,2],projected[:,1],c=vturb,edgecolor='none',cmap=plt.cm.get_cmap('afmhot'))
    plt.xlabel('component 3')
    plt.ylabel('component 2')
    plt.subplot(243)
    plt.scatter(projected[:,2],projected[:,3],c=vturb,edgecolor='none',cmap=plt.cm.get_cmap('afmhot'))
    plt.xlabel('component 3')
    plt.ylabel('component 4')
    plt.colorbar(label='v$_{turb}$ (c$_s$)')

    #Display the components to see what they look like
    pca.fit(images)
    plt.subplot(2,8,9)
    plt.imshow(pca.components_[0].reshape(260,150),cmap=plt.cm.get_cmap('afmhot'))
    plt.title('1')
    plt.subplot(2,8,10)
    plt.imshow(pca.components_[1].reshape(260,150),cmap=plt.cm.get_cmap('afmhot'))
    plt.title('2')
    plt.subplot(2,8,11)
    plt.imshow(pca.components_[2].reshape(260,150),cmap=plt.cm.get_cmap('afmhot'))
    plt.title('3')
    plt.subplot(2,8,12)
    plt.imshow(pca.components_[3].reshape(260,150),cmap=plt.cm.get_cmap('afmhot'))
    plt.title('4')
    plt.subplot(247)
    plt.imshow(image_mean.reshape(260,150),cmap=plt.cm.get_cmap('afmhot'))
    plt.title('Mean')

    #Use regression techniques to correlate the components with turbulence
    #Start with Linear Regression
    model = LinearRegression(fit_intercept=True)
    model.fit(projected,vturb)

    #Try to recover the input values of turbulence
    vturb_predicted = model.predict(projected)
    
    plt.subplot(144)
    plt.plot(vturb,vturb/vturb_predicted,'ok')
    plt.axhline(1,color='k',lw=3,ls='--')
    #plt.plot([0,np.max(np.append(vturb,vturb_predicted))],[0,np.max(np.append(vturb,vturb_predicted))],color='k',ls='--',lw=2)
    plt.xlabel('Actual v$_{turb}$')
    plt.ylabel('Predicted v$_{turb}$')
    print np.std(vturb/vturb_predicted)
    

def pca_comp():
    '''Check to see how many components are needed to fit the spectra.'''

    from sklearn.decomposition import PCA
    import matplotlib.pyplot as plt
    import numpy as np

    vturb,spec = read_models()

    pca = PCA().fit(spec)
    plt.figure()
    plt.plot(np.cumsum(pca.explained_variance_ratio_))
    plt.xlabel('number of components')
    plt.ylabel('cumulative explained variance')



def isomap_spec(n_components=2,degree=2):
    '''Use IsoMap instead of PCA to do the dimensionality reduction among the spectra.'''

    from sklearn.manifold import Isomap
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.linear_model import LinearRegression
    from sklearn.neighbors import KernelDensity
    from random import sample

    vturb,spec = read_models()

    #Divide into training and validation sets
    nmodels = len(vturb)
    val_list = np.array(sample(range(nmodels),int(.1*nmodels)))
    vturb_train = np.zeros(int(.9*nmodels))
    spec_train = np.zeros((int(.9*nmodels),spec.shape[1]))
    vturb_val = np.zeros(int(.1*nmodels))
    spec_val = np.zeros((int(.1*nmodels),spec.shape[1]))
    val_count = 0
    train_count = 0
    for i in range(nmodels):
        if (i==val_list).sum()>0:
            #part of validation sample
            vturb_val[val_count] = vturb[i]
            spec_val[val_count,:] = spec[i,:]
            val_count+=1
        else:
            vturb_train[train_count] = vturb[i]
            spec_train[train_count,:] = spec[i,:]
            train_count+=1
        

    #Use MDS with 2 components
    iso = Isomap(n_components=n_components)
    proj = iso.fit_transform(spec_train)
    

    #Use polynomial regression to correlate the coefficients on the components with turbulence
    from sklearn.preprocessing import PolynomialFeatures
    poly = PolynomialFeatures(degree=degree,include_bias=False)
    proj2 = poly.fit_transform(proj)

    model2 = LinearRegression(fit_intercept=True)
    model2.fit(proj2,vturb_train)
    print 'R^2 of the fit: ',model2.score(proj2,vturb_train)
    
    vturb_predicted2 = model2.predict(proj2)

    plt.subplot(211)
    #plt.plot(vturb,vturb_predicted,'or',ms=1)
    plt.plot(vturb_train,vturb_predicted2,'ok',ms=1)
    plt.plot([0,np.max(np.append(vturb_train,vturb_predicted2))],[0,np.max(np.append(vturb_train,vturb_predicted2))],color='k',ls='--',lw=2)
    plt.xlabel('Actual v$_{turb}$')
    plt.ylabel('Predicted v$_{turb}$')

    #Run against the validation set and see the accuracy of the result
    #vturb_val_predicted = np.array(len(vturb_train))
#    for i in len(vturb_train):
    pspec = iso.transform(spec_val)
    pspec2 = poly.fit_transform(pspec)
    vturb_val_predicted = model2.predict(pspec2)
    plt.plot(vturb_val,vturb_val_predicted,'or',ms=1)
    error = (vturb_val-vturb_val_predicted)/vturb_val
    print 'Median fractional error in validation set: ',np.median(np.abs(error))

#n_comp = 2: .33, .32, .30, .37, .30
#n_comp = 3: .32, .32, .34, .32, .29
#n_comp = 4: .34, .28, .27, .25, .27
#n_comp = 5: .30, .25, .25, .31, .26
#n_comp = 6: .28, .26, .27, .25, .28
#n_comp = 7: .27, .24, .23*, .27, .36

    disk_spec = np.array([0.07051709,  0.08370663,  0.09472561,  0.14514502,  0.19218614, 0.26022604,  0.3638115 ,  0.51294306,  0.77375685,  1.23985817, 1.86580414,  2.92085159,  4.45566949,  7.49341501,  9.57812521, 8.47671087,  7.11170747,  7.57635818,  9.61814015,  8.07201153, 5.09846799,  3.1666609 ,  1.96394516,  1.3104509 ,  0.82953736, 0.56494169,  0.38539208,  0.25932753,  0.18477959,  0.13578935, 0.10984786])
    vturb_actual = .275
    disk_spec = disk_spec[:,np.newaxis].T
    #disk_spec_pca = pca.transform(disk_spec)
    #disk_spec_pca2 = poly.fit_transform(disk_spec_pca)
    vturb_disk_predicted = []
    for i in range(100):
        noise = np.random.randn(len(disk_spec))*.01*disk_spec
        disk_spec_pca = iso.transform(disk_spec+noise)
    #    vturb_disk_predicted.append(model.predict(disk_spec_pca))
        disk_spec_pca2 = poly.fit_transform(disk_spec_pca)
        vturb_disk_predicted.append(model2.predict(disk_spec_pca2))
    vturb_disk_predicted = np.array(vturb_disk_predicted)
    #vturb_disk_predicted = .2#model2.predict(disk_spec_pca2)
    #20% higher/lower for systematic uncertainty
    disk_spec_pcaH = iso.transform(disk_spec*1.2)
    disk_spec_pcaL = iso.transform(disk_spec*0.8)
    disk_spec_pcaH2 = poly.fit_transform(disk_spec_pcaH)
    disk_spec_pcaL2 = poly.fit_transform(disk_spec_pcaL)
    vturb_disk_predictedH = model2.predict(disk_spec_pcaH2)
    vturb_disk_predictedL = model2.predict(disk_spec_pcaL2)
    #plt.subplot(236)
    #plt.plot([.01,],[vturb_disk_predicted,],'sb')
    plt.plot(np.zeros(100)+vturb_actual,vturb_disk_predicted,'sb',ms=2)
    plt.plot([vturb_actual,],[vturb_disk_predictedH,],'sc',ms=5)
    plt.plot([vturb_actual,],[vturb_disk_predictedL,],'sc',ms=5)
    #plt.axhline(.9*vturb_disk_predicted,color='k',ls=':')
    #plt.axhline(1.1*vturb_disk_predicted,color='k',ls=':')

    print 'Disk actual, predicted: ',vturb_actual,np.mean(vturb_disk_predicted),np.std(vturb_disk_predicted)
    print 'Median % error: ',np.median(np.abs(vturb_train-vturb_predicted2)/vturb_train)

    #plot the derived and original spectral to see if they match
    #plt.subplot(222)
    #plt.plot(disk_spec[0,:],color='k',lw=2)
    #plt.plot((iso.inverse_transform(disk_spec_pca))[0,:],color='r',lw=2,ls='--')
    #plt.legend(('Original','Reconstructed'),frameon=False,loc='best')

    #This doesn't seem like an actual improvement over PCA+polynomial regression


    #Derive the probability distribution of actual values of turbulence based on the predicted level of turbulence
    w=(vturb_predicted2>vturb_disk_predicted.min()) & (vturb_predicted2<vturb_disk_predicted.max())
    plt.subplot(212)

    kde = KernelDensity(bandwidth=.03,kernel='gaussian')
    kde.fit(vturb_train[w][:,None])

    v_d = np.linspace(0,1,100)
    logprob = kde.score_samples(v_d[:,None])

    plt.plot(v_d,np.exp(logprob),color='k')
    plt.axvline(vturb_actual,color='k',ls=':')
    plt.axvline(np.mean(vturb_disk_predicted),color='b',ls='--')
    plt.axvline(vturb_disk_predictedH,color='c',ls='--')
    plt.axvline(vturb_disk_predictedL,color='c',ls='--')
    plt.xlim(0,1)
    plt.ylabel('PDF')
    plt.xlabel('Actual v$_{turb}$')


def pcann_spec(n_components=6):
    '''Read in the spectra and perform ICA (from astroML) to reduce the dimensionality of the data. The result that then can be regressed against the turbulence to find out how the spetrum varies with turbulence.'''
    from sklearn.decomposition import PCA
    from sklearn.neighbors import KernelDensity
    #from sklearn.grid_search import GridSearchCV
    #from sklearn.cross_validation import LeaveOneOut
    import matplotlib.pyplot as plt
    import numpy as np
    from random import sample

    vturb,spec = read_models()

    #Divide into training and validation sets
    nmodels = len(vturb)
    val_list = np.array(sample(range(nmodels),int(.1*nmodels)))
    vturb_train = np.zeros(int(.9*nmodels))
    spec_train = np.zeros((int(.9*nmodels),spec.shape[1]))
    vturb_val = np.zeros(int(.1*nmodels))
    spec_val = np.zeros((int(.1*nmodels),spec.shape[1]))
    val_count = 0
    train_count = 0
    for i in range(nmodels):
        if (i==val_list).sum()>0:
            #part of validation sample
            vturb_val[val_count] = vturb[i]
            spec_val[val_count,:] = spec[i,:]
            val_count+=1
        else:
            vturb_train[train_count] = vturb[i]
            spec_train[train_count,:] = spec[i,:]
            train_count+=1
            


    #Use PCA with 2 components. 
    #First plot the values of the two components, with turbulence highlighted
    #n_components = 8
    pca = PCA(n_components)
    projected = pca.fit_transform(spec_train)
    
    


    #Plot the components to see what they look like
    plt.figure()
    plt.subplot(221)
    for i in range(n_components):
        plt.plot(pca.components_[i],label='component {:0.0f}'.format(i),lw=n_components-i)
    #plt.legend(frameon=False,loc='best')

    print np.cumsum(pca.explained_variance_ratio_)


    #Standardize the test set
    #from sklearn.preprocessing import StandardScaler
    #scaler = StandardScaler()
    #scaler.fit(projected)
    #projected = scaler.transform(projected)

    #Try a neural network for predicting the turbulence
    from sklearn.neural_network import MLPRegressor
    #clf = MLPRegressor(solver='lbfgs',alpha=1e-5,hidden_layer_sizes=(5,2))
    clf = MLPRegressor(solver='lbfgs',alpha=1e-5,hidden_layer_sizes=(5,2),activation='logistic')

    clf.fit(projected,vturb_train)
    print 'R^2 of the fit: ',clf.score(projected,vturb_train)
    
    
    #How well can we recover the input turbulence?
    vturb_predicted2 = clf.predict(projected) #polynomial model

    plt.subplot(223)
    ##plt.plot(vturb,vturb_predicted,'or',ms=1)
    plt.plot(vturb_train,vturb_predicted2,'ok',ms=1)
    plt.plot([0,np.max(np.append(vturb,vturb_predicted2))],[0,np.max(np.append(vturb,vturb_predicted2))],color='c',ls='--',lw=2)
    plt.xlabel('Actual v$_{turb}$')
    plt.ylabel('Predicted v$_{turb}$')


    #Run against the validation set and see the accuracy of the result
    #vturb_val_predicted = np.array(len(vturb_train))
#    for i in len(vturb_train):
    pspec = pca.transform(spec_val)
    #pspec = scaler.transform(pspec) 
    vturb_val_predicted = clf.predict(pspec)
    plt.plot(vturb_val,vturb_val_predicted,'or',ms=1)
    error = (vturb_val-vturb_val_predicted)/vturb_val
    print 'Median fractional error in validation set: ',np.median(np.abs(error))

#layer_sizes = (5,2), alpha=1e-5, solver=lbfgs
#n_components=2: .32 (relu), .32 (logistic), .33 (tanh), .32 (identity)
#n_components=3: .24 (relu), .22 (logistic), .21 (tanh), .27 (identity)
#n_components=4: .22 (relu), .21 (logistic), .22 (tanh), .26 (identity)
#n_components=5: .20 (relu), .20 (logistic), .21 (tanh), .29 (identity)
#n_components=6: .24 (relu), .18 (logistic), .18 (tanh), .26 (identity)
#n_components=7: .22 (relu), .18 (logistic), .20 (tanh), .27 (identity)

#for logistic, the larger the number of PCA components, the more it prefers a weaker level of turbulence when fitting DM Tau 
#identity seems to recover the DM Tau level of turbulence most often, but with a wide posterior

#layer_sizes = (5,2), alpha=1e-5, activation=logistic
#n_comp = 2: .38 (sgd), .34 (adam)
#n_comp = 3: .46 (sgd), .29 (adam)
#n_comp = 4: .45 (sgd), .30 (adam)
#n_comp = 5: .37 (sgd), .27 (adam)
#n_comp = 6: .40 (sgd), .24 (adam)
#n_comp = 7: .41 (sgd), .32 (adam)

#stochastic gradient descent (sgd) and 'adam' (a stochastic gradient based optimizer) have no predictive power. 

#solver=lbfgs, activation=logistic, alpha=1e-5, one layer
#n_comp = 3: .24 (2), .22 (10), .24 (100)
#n_comp = 4: .22 (2), .22 (10), .21 (100)
#n_comp = 5: .24 (2), .22 (10), .20 (100)
#n_comp = 6: .20 (2), .19 (10), .20 (100)

#10-100 neurons seems to work okay, but not great

#solver=lbfgs, activation=logistic, alpha=1e-5, 2 layers
#n_comp = 3: .22 (10,5), .22 (50,5), .22 (50,20), .20 (100,20)
#n_comp = 4: .24 (10,5), .21 (50,5), .20 (50,20), .20 (100,20)
#n_comp = 5: .19 (10,5), .21 (50,5), .20 (50,20), .18 (100,20)
#n_comp = 6: .20 (10,5), .18 (50,5), .18 (50,20), .19 (100,20)

#solver=lbfgs, activation=logistic, alpha=1e-5, 3 layers
#n_comp = 3: .23 (10,5,2), .22 (100,50,20)
#n_comp = 4: .21 (10,5,2), .21 (100,50,20)
#n_comp = 5: .19 (10,5,2), .20 (100,50,20)
#n_comp = 6: .19 (10,5,2), .20 (100,50,20)

#changing the number of layers does not dramatically change the results





    #Predict turbulence in our three disks from this decomposition
    disk_spec = np.array([0.07051709,  0.08370663,  0.09472561,  0.14514502,  0.19218614, 0.26022604,  0.3638115 ,  0.51294306,  0.77375685,  1.23985817, 1.86580414,  2.92085159,  4.45566949,  7.49341501,  9.57812521, 8.47671087,  7.11170747,  7.57635818,  9.61814015,  8.07201153, 5.09846799,  3.1666609 ,  1.96394516,  1.3104509 ,  0.82953736, 0.56494169,  0.38539208,  0.25932753,  0.18477959,  0.13578935, 0.10984786])
    vturb_actual = .275
    disk_spec = disk_spec[:,np.newaxis].T
    disk_spec_pca = pca.transform(disk_spec)
    #disk_spec_pca = scaler.transform(disk_spec_pca)
    vturb_disk_predicted = []
    for i in range(100):
        noise = np.random.randn(len(disk_spec))*.01*disk_spec
        disk_spec_pca = pca.transform(disk_spec+noise)
        vturb_disk_predicted.append(clf.predict(disk_spec_pca))
    vturb_disk_predicted = np.array(vturb_disk_predicted)
    
    #20% higher/lower for systematic uncertainty
    disk_spec_pcaH = pca.transform(disk_spec*1.2)
    disk_spec_pcaL = pca.transform(disk_spec*0.8)
    vturb_disk_predictedH = clf.predict(disk_spec_pcaH)
    vturb_disk_predictedL = clf.predict(disk_spec_pcaL)
    plt.plot(np.zeros(100)+vturb_actual,vturb_disk_predicted,'sb',ms=2)
    plt.plot([vturb_actual,],[vturb_disk_predictedH,],'sc',ms=5)
    plt.plot([vturb_actual,],[vturb_disk_predictedL,],'sc',ms=5)
    
    print 'Disk actual, predicted: ',vturb_actual,np.mean(vturb_disk_predicted),np.std(vturb_disk_predicted)
    print 'Median % error: ',np.median(np.abs(vturb_train-vturb_predicted2)/vturb_train)

    #plot the derived and original spectral to see if they match
    plt.subplot(222)
    plt.plot(disk_spec[0,:],color='k',lw=2)
    plt.plot((pca.inverse_transform(disk_spec_pca))[0,:],color='r',lw=2,ls='--')
    plt.legend(('Original','Reconstructed'),frameon=False,loc='best')



    #Derive the probability distribution of actual values of turbulence based on the predicted level of turbulence
    w=(vturb_predicted2>vturb_disk_predicted.min()) & (vturb_predicted2<vturb_disk_predicted.max())
    plt.subplot(224)
    #bandwidths = 10**(np.linspace(-2,0,100))
    #grid = GridSearchCV(KernelDensity(kernel='gaussian'),{'bandwidth':bandwidths},cv=LeaveOneOut(len(vturb_train[w])))
    #grid.fit(vturb_train[w][:,None])
    #print grid.best_params_

    if w.sum()>0:
        kde = KernelDensity(bandwidth=.03,kernel='gaussian')
        kde.fit(vturb_train[w][:,None])

        v_d = np.linspace(0,1,100)
        logprob = kde.score_samples(v_d[:,None])

        plt.plot(v_d,np.exp(logprob),color='k')

    
    plt.axvline(vturb_actual,color='k',ls=':')
    plt.axvline(np.mean(vturb_disk_predicted),color='b',ls='--')
    plt.axvline(vturb_disk_predictedH,color='c',ls='--')
    plt.axvline(vturb_disk_predictedL,color='c',ls='--')
    plt.xlim(0,1)
    plt.ylabel('PDF')
    plt.xlabel('Actual v$_{turb}$')


    
def nn_spec():
    '''Run a neural network directly on the spectra, without any PCA.'''

    from sklearn.neighbors import KernelDensity
    #from sklearn.grid_search import GridSearchCV
    #from sklearn.cross_validation import LeaveOneOut
    import matplotlib.pyplot as plt
    import numpy as np
    from random import sample

    vturb,spec = read_models()

    #Divide into training and validation sets
    nmodels = len(vturb)
    val_list = np.array(sample(range(nmodels),int(.1*nmodels)))
    vturb_train = np.zeros(int(.9*nmodels))
    spec_train = np.zeros((int(.9*nmodels),spec.shape[1]))
    vturb_val = np.zeros(int(.1*nmodels))
    spec_val = np.zeros((int(.1*nmodels),spec.shape[1]))
    val_count = 0
    train_count = 0
    for i in range(nmodels):
        if (i==val_list).sum()>0:
            #part of validation sample
            vturb_val[val_count] = vturb[i]
            spec_val[val_count,:] = spec[i,:]
            val_count+=1
        else:
            vturb_train[train_count] = vturb[i]
            spec_train[train_count,:] = spec[i,:]
            train_count+=1


    #Standardize the test and training set
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    scaler.fit(spec_train)
    spec_train = scaler.transform(spec_train)
    spec_val = scaler.transform(spec_val)
      
    #Try a neural network for predicting the turbulence
    from sklearn.neural_network import MLPRegressor
    clf = MLPRegressor(solver='lbfgs',alpha=1e-4,hidden_layer_sizes=(100,50,20),activation='relu')

    clf.fit(spec_train,vturb_train)
    print 'R^2 of the fit: ',clf.score(spec_train,vturb_train)

    #How well can we recover the input turbulence?
    vturb_predicted2 = clf.predict(spec_train) #polynomial model

    plt.figure()
    plt.subplot(211)
    ##plt.plot(vturb,vturb_predicted,'or',ms=1)
    plt.plot(vturb_train,vturb_predicted2,'ok',ms=1)
    plt.plot([0,np.max(np.append(vturb,vturb_predicted2))],[0,np.max(np.append(vturb,vturb_predicted2))],color='c',ls='--',lw=2)
    plt.xlabel('Actual v$_{turb}$')
    plt.ylabel('Predicted v$_{turb}$')

    #Run against the validation set and see the accuracy of the result
    vturb_val_predicted = clf.predict(spec_val)
    plt.plot(vturb_val,vturb_val_predicted,'or',ms=1)
    error = (vturb_val-vturb_val_predicted)/vturb_val
    print 'Median fractional error in validation set: ',np.median(np.abs(error))

    #Predict turbulence in our three disks from this decomposition
    disk_spec = np.array([0.07051709,  0.08370663,  0.09472561,  0.14514502,  0.19218614, 0.26022604,  0.3638115 ,  0.51294306,  0.77375685,  1.23985817, 1.86580414,  2.92085159,  4.45566949,  7.49341501,  9.57812521, 8.47671087,  7.11170747,  7.57635818,  9.61814015,  8.07201153, 5.09846799,  3.1666609 ,  1.96394516,  1.3104509 ,  0.82953736, 0.56494169,  0.38539208,  0.25932753,  0.18477959,  0.13578935, 0.10984786])
    vturb_actual = .275
    disk_spec = disk_spec[:,np.newaxis].T
    disk_spec = scaler.transform(disk_spec)
    vturb_disk_predicted = []
    for i in range(100):
        noise = np.random.randn(len(disk_spec))*.01*disk_spec
        vturb_disk_predicted.append(clf.predict(disk_spec+noise))
    vturb_disk_predicted = np.array(vturb_disk_predicted)
 
    #20% higher/lower for systematic uncertainty
    vturb_disk_predictedH = clf.predict(disk_spec*1.2)
    vturb_disk_predictedL = clf.predict(disk_spec*0.8)
    plt.plot(np.zeros(100)+vturb_actual,vturb_disk_predicted,'sb',ms=2)
    plt.plot([vturb_actual,],[vturb_disk_predictedH,],'sc',ms=5)
    plt.plot([vturb_actual,],[vturb_disk_predictedL,],'sc',ms=5)
    
    print 'Disk actual, predicted: ',vturb_actual,np.mean(vturb_disk_predicted),np.std(vturb_disk_predicted)
    print 'Median % error: ',np.median(np.abs(vturb_train-vturb_predicted2)/vturb_train)

    #Derive the probability distribution of actual values of turbulence based on the predicted level of turbulence
    w=(vturb_predicted2>vturb_disk_predicted.min()) & (vturb_predicted2<vturb_disk_predicted.max())
    plt.subplot(212)
    #bandwidths = 10**(np.linspace(-2,0,100))
    #grid = GridSearchCV(KernelDensity(kernel='gaussian'),{'bandwidth':bandwidths},cv=LeaveOneOut(len(vturb_train[w])))
    #grid.fit(vturb_train[w][:,None])
    #print grid.best_params_

    if w.sum()>0:
        kde = KernelDensity(bandwidth=.03,kernel='gaussian')
        kde.fit(vturb_train[w][:,None])

        v_d = np.linspace(0,1,100)
        logprob = kde.score_samples(v_d[:,None])

        plt.plot(v_d,np.exp(logprob),color='k')

    
    plt.axvline(vturb_actual,color='k',ls=':')
    plt.axvline(np.mean(vturb_disk_predicted),color='b',ls='--')
    plt.axvline(vturb_disk_predictedH,color='c',ls='--')
    plt.axvline(vturb_disk_predictedL,color='c',ls='--')
    plt.xlim(0,1)
    plt.ylabel('PDF')
    plt.xlabel('Actual v$_{turb}$')

