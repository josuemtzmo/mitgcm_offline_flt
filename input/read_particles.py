import os
import re
import numpy as np
import pandas as pd
import pylab as plt
import struct
import time

class ReadParticles():
    '''
    Class to read particle output from MITgcm.
    '''

    def __init__(self, filename,ffile=None):
        
        if ffile==None:
            ffile = filename.split('/')[-1][0:-1]
            
        if '*' in filename:
            folder = filename.split('*')[0]
            folder   = folder.split(folder.split('/')[-1])[0]
            
            files = sorted([os.path.join(folder, f) for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and ('.data' in f and ffile in f)])
            self.filename=files
            self.load_mftracks()
        else: 
            self.filename=filename
            self.load_tracks()

    def load_mftracks(self):
        printProgressBar(0, len(self.filename), prefix = 'Progress:', suffix = '', length = 10)
        count=0
        tic=time.time()
        for file in self.filename:
            if count==0:
                tmp_tracks=self.load_tracks(file)
            else:
                pd0=self.load_tracks(file)
                tmp_tracks=pd.concat([tmp_tracks,pd0])
            
            toc=time.time()    
            printProgressBar(count + 1, len(self.filename), prefix = 'Progress:', suffix = 'Ellapsed time: {0}'.format(round(toc-tic,2)), length = 10)
            count+=1
            
        self.df_tracks = tmp_tracks
        return self.df_tracks

    def load_tracks(self,file=None):
        '''
        Read particle binary file (modified mds.rdmds function)
        '''
        if file == None and type(self.filename)!=list:
            file=self.filename
        elif type(self.filename)==list and file == None:
            raise ValueError('Use the function load_mftracks to load multiple files.')

        dimList = [     1,    1,    1,     1,    1,    1,    16,    1,   16 ];     
        self.records=int(dimList[-1])
        
        byterate = 8
        tp = '>d{0}'

        recshape = self.records
        
        nrecords, = [ int(os.path.getsize(file)//(byterate)) ]
        tileshape = (nrecords//recshape, recshape)

        lev=()
        if not isinstance(lev,tuple):
            lev = (lev,)
        levs = tuple( aslist(l) for l in lev )
        levdims = tuple(len(l) for l in levs)

        reclist = range(nrecords)
        
        fmt = tp.format('d'*(nrecords-1))

        struct_size = struct.calcsize(fmt)

        with open(file, mode='rb') as file:
            fileContent = file.read()
            arrtile=np.array(struct.unpack(fmt, fileContent[:struct_size]))
        
        self.arrtile=arrtile.reshape(tileshape)
        self.bin2pd()
        return self.df_tracks

    def bin2pd(self,numcols=16):
        if self.records < numcols:
            numcols = self.records
        cols_name=['x','y','z','i','j','k','p','u','v','t','s','vort','ke','t2ave'
        ,'c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12','c13','c14']
        if numcols > len(cols_name):
            raise ValueError('''Most of the tracking data is contained in the 
                                first 13 collumns, if you need more outputs, 
                                modify the cols_name variable in this function.''')
        index = pd.MultiIndex.from_arrays(self.arrtile[:,0:2].T, names=('particle_id', 'time'))
        df_tracks = pd.DataFrame(np.array(self.arrtile[:,2:numcols],dtype=np.float64),
                    columns=cols_name[0:numcols-2],index=index)
        self.df_tracks=df_tracks.sort_index()
    
    def plot_tracks(self,latlon=True,**kargs):
        if latlon==True:
            import cartopy.crs as ccrs
            transform=ccrs.PlateCarree()
        else:
            transform=None
        
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1,projection=transform)

        index = np.sort(self.df_tracks.index.get_level_values(0).unique())
        
        printProgressBar(0, len(index), prefix = 'Progress:', suffix = '', length = 10)
        count=0
        for ii in index:
            x=self.df_tracks.loc[ii]['x']
            y=self.df_tracks.loc[ii]['y']
            plt.plot(x,y,transform=transform)
            printProgressBar(count+1, len(index), prefix = 'Progress:', suffix = 'Track: {0}'.format(ii), length = 10)
            count+=1
            
        if latlon==True:
            ax.coastlines()
    
    def plot_initp(self,latlon=True,**kargs):
        if latlon==True:
            import cartopy.crs as ccrs
            transform=ccrs.PlateCarree()
        else:
            transform=None
        
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1,projection=transform)

        for ii in np.sort(self.df_tracks.index.get_level_values(0).unique()):
            x=self.df_tracks.loc[ii]['x'].iloc[0]
            y=self.df_tracks.loc[ii]['y'].iloc[0]
            plt.plot(x,y,'og',transform=transform)
            
        if latlon==True:
            ax.coastlines()
            
    
    def savedf(self,path):
        self.df_tracks.to_csv(path)

    def LAVD(self):
        timediff=(self.df_tracks.loc[1].index[-1]-self.df_tracks.loc[1].index[0])/86400.0
        LAVD = (1/timediff) * self.df_tracks.sum(level='particle_id')
        return LAVD

    

def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()
        
# if __name__ == "__main__":
#     p_file = ReadParticles('./float_profiles.001.001.data')
#     dic_file= p_file.readBinary()
#     p_file.plot_tracks()