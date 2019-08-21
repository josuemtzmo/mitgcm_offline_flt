import os
import re
import numpy as np
import pandas as pd
import pylab as plt

class ReadParticles():
    '''
    Class to read particle output from MITgcm.
    '''

    def __init__(self, filename):
        self.filename=filename

    def readBinary(self):
        '''
        Read particle binary file (modified mds.rdmds function)
        '''
        from operator import mul
        import functools

        dimList = [ 1,    1,    1,  1,    1,    1, 41,    1,   41 ];      

        gdims = tuple(dimList[-3::-3]) 
        i0s   = [ i-1 for i in dimList[-2::-3] ] 
        ies   = dimList[-1::-3] 
        byterate = 4
        tp = '>f{0}'.format(byterate)
        size = np.dtype(tp).itemsize

        ri0,rie,rj0,rje = 0,gdims[-1],0,gdims[-2]

        recshape = tuple( ie-i0 for i0,ie in zip(i0s,ies) )
        count = functools.reduce(mul, recshape)
        nrecords, = [ int(os.path.getsize(self.filename)/(count*byterate)) ]
        tileshape = (nrecords,) + recshape

        lev=()
        if not isinstance(lev,tuple):
            lev = (lev,)
        levs = tuple( aslist(l) for l in lev )
        levdims = tuple(len(l) for l in levs)

        itrs=[0]

        dims = levdims + gdims[len(levdims):-2] + (rje-rj0,rie-ri0)

        rdims = levdims + gdims[len(levdims):-2] + (rje-rj0,rie-ri0)

        reclist = range(nrecords)
        astype = tp

        arr = np.empty((len(itrs),len(reclist))+rdims, astype) 
        arr[...] = 0  

        ny,nx = arr.shape[-2:]
        sl = tuple( slice(i0,ie) for i0,ie in zip(i0s,ies) )

        for iit,it in enumerate(itrs):

            arrtile = arr[(iit,slice(None))+sl]
            f = open(self.filename)
            for irec,recnum in enumerate(reclist): 
                arrtile[irec] = np.fromfile(f, tp, count=count).reshape(recshape)
            f.close()
        self.arrtile=arrtile.squeeze()
        self.mitgcmtrack2pd()
        return self.df_tracks


    def mitgcmtrack2pd(self,numcols=13):
        cols_name=['x','y','z','i','j','k','p','u','v','t','s','c1','c2','c3'
        ,'c4','c5','c6','c7','c8','c9','c10','c11','c12','c13','c14','c15','c16']
        if numcols > len(cols_name):
            raise ValueError('''Most of the tracking data is contained in the 
                                first 13 collumns, if you need more outputs, 
                                modify the cols_name variable in this function.''')
        index = pd.MultiIndex.from_arrays(self.arrtile[:,0:2].T, names=('particle_id', 'time'))
        df_tracks = pd.DataFrame(np.array(self.arrtile[:,2:numcols+2],dtype=np.float32),
                    columns=cols_name[0:numcols],index=index)
        self.df_tracks=df_tracks.sort_index()
    
    def plot_tracks(self,latlon=True,**kargs):
        if latlon==True:
            import cartopy.crs as ccrs
            transform=ccrs.PlateCarree()
        else:
            transform=None
        
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1,projection=transform)

        for ii in np.sort(self.df_tracks.index.get_level_values(0).unique()):
            x=self.df_tracks.loc[ii]['x']
            y=self.df_tracks.loc[ii]['y']
            plt.plot(x,y,transform=transform)
            
        if latlon==True:
            ax.coastlines()

# if __name__ == "__main__":
#     p_file = ReadParticles('./float_profiles.001.001.data')
#     dic_file= p_file.readBinary()
#     p_file.plot_tracks()
