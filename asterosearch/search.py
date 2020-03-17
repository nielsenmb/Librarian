import astropy.units as u
import numpy as np
import pandas as pd
from astropy.table import Table 
from astropy.table import vstack as avstack
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.mast import Catalogs
Simbad.add_votable_fields('sptype', 'ids')
from tqdm import tqdm
import re, time

def kic_to_kplr(x):
    y = x.strip('KIC')
    y = 'kplr'+'0'*(9-len(y))+y
    return y

def kplr_to_kic(x):
    y = x.strip('kplr')
    return(str(int(y)))

def format_name(name):
    # Add naming exceptions here
    
    if any([x in name for x in ['KIC', 'kic', 'kplr', 'KPLR']]):
        name = name.replace('KIC','')
        name = re.sub(r"\s+", "", name, flags=re.UNICODE)
        name = 'KIC '+name

    if 'TIC' in name:
        name = name.replace('TIC','')
        name = re.sub(r"\s+", "", name, flags=re.UNICODE)
        name = 'TIC '+name
    
    if 'EPIC' in name:
        name = name.replace('EPIC','')
        name = re.sub(r"\s+", "", name, flags=re.UNICODE)
        name = 'EPIC '+name
        
    if 'GDR2' in name:
        for x in ['Gaia', 'GAIA', 'GDR2', 'Gaia DR2', 'GAIA DR2', 'DR2']:
            name = name.replace(x,'')
        name = re.sub(r"\s+", "", name, flags=re.UNICODE)
        name = 'Gaia DR2 '+name
        
    if 'GDR1' in name:
        for x in ['Gaia', 'GAIA', 'GDR2', 'Gaia DR2', 'GAIA DR2', 'DR2']:
            name = name.replace(x,'')
        name = re.sub(r"\s+", "", name, flags=re.UNICODE)
        name = 'Gaia DR1 '+name
    return name

def add_empty_row(tbl):
    tbl.add_row([0]*len(tbl.keys()), mask=[True]*len(tbl.keys()))

def add_to_table(job, tbl, identifier):
    if len(job) == 0:
        add_empty_row(tbl)
    else:
        for i, row in enumerate(job[0]):
            if row[identifier] in id:
                tbl = avstack([tbl, row])
                break

class search():
     
    def __init__(self, ID):
        
        self.cats_avail =  ['SPOCS', 'KIC', 'TIC', 'EPIC', 'Gaia DR1', 'Gaia DR2', 'HD', 
                            'HIP', 'HR', 'HIC', 'UBV', 'SAO', 'GEN#', 'TYC', '2MASS', 
                            'GJ', 'PPM', 'BD', 'AG', 'GSC', 'Ci', 'PLX', 'SKY#', 'WISEA', 
                            'WISE', 'PSO', 'ALLWISE']
               
        if isinstance(ID, str):
            ID = [ID]
        
        self.IDs = pd.DataFrame({'input': ID})
        
        for cat in self.cats_avail:
            self.IDs[cat] = np.empty(len(self.IDs),dtype=str)

            for i in self.IDs.index:
                self.IDs.input[i] = format_name(self.IDs.input[i])
           
                if cat in self.IDs.input[i]:                    
                    self.IDs.at[i, cat] = self.IDs.input[i]
         
    def __call__(self):
        self.query_simbad()
        self.query_KIC()
        #self.query_TIC() 
        self.query_GaiaDR2()
        
    def query_simbad(self, ID=None):
        time.sleep(1) # This is to keep you from being blacklisted by Simbad servers       

        if ID is None:
            ID = self.IDs.input
            
        job = Simbad.query_objects(ID)
        if job is None:
            print('Unable to find any of the targets with Simbad')      
            self.simbad = None
            return self.simbad
        
        tbl = job.copy()
        tbl.remove_rows(range(len(job)))
        
        for i, id in tqdm(enumerate(ID)):
            add_empty_row(tbl)
            for j, simIDs in enumerate(job['IDS']):
                simIDs = str(simIDs)
                if id in simIDs:
                    tbl[i] = job[j]
                    
                    simIDslist = simIDs.strip("b'").split('|')
                    if (len(simIDslist) == 0) or (len(simIDslist) > 100):
                        raise ValueError('Something went wrong trying to split the simbad IDS result')
                    else:
                        for k, simID in enumerate(simIDslist): # Assign names to ID dataframe
                            for c, cat in enumerate(self.cats_avail):
                                if (cat in simID) and (len(self.IDs.loc[i, cat])==0):
                                    self.IDs.loc[i, cat] = simID
        self.simbad = tbl
        return self.simbad        
                    
    def query_KIC(self, ID = None, radius = 10.0*u.arcsec):
        import warnings
        from astropy.utils.metadata import MergeConflictWarning
        warnings.filterwarnings("ignore", category = MergeConflictWarning)
        
        if ID is None:
            ID = self.IDs['KIC']
                
        tbl = Table(names = ('KIC',), dtype = (int,))
        for i, id in tqdm(enumerate(ID)):
            if not isinstance(id, str):
                add_empty_row(tbl)
            else:
                job = Vizier.query_object(object_name = id, catalog = 'V/133/kic', radius = radius)

                if len(job) > 0:
                    tbl = avstack([tbl, job[0]])      
                else:
                    add_empty_row(tbl)
                    
        self.KIC = tbl
        return self.KIC
    
    def query_GaiaDR2(self, ID=None):
        from astroquery.gaia import Gaia

        if ID is None:
            ID = self.IDs['Gaia DR2']

        tbl = Table(names = ('source_id',), dtype = (int,))
        for i, gid in tqdm(enumerate(ID)):
            if not isinstance(gid, str):
                add_empty_row(tbl)
            else:
                gid = int(gid.replace('Gaia DR2 ', ''))

                adql_query = "select * from gaiadr2.gaia_source where source_id=%i" % (gid)

                job = Gaia.launch_job(adql_query).get_results()

                idx = np.where(job['source_id'].quantity == gid)[0]
                
                if len(idx) > 0:
                    tbl = avstack([tbl, job[idx]])
                else:
                    add_empty_row(tbl)
           
        self.GDR2 = tbl
        return self.GDR2
    
    def query_EPIC(self, ID=None, radius = 10.0*u.arcsec):
        import warnings
        from astropy.utils.metadata import MergeConflictWarning
        warnings.filterwarnings("ignore", category = MergeConflictWarning)
        
        if ID is None:
            ID = self.IDs['EPIC']
                
        tbl = Table(names = ('ID',), dtype = (int,))
        
        for i, id in tqdm(enumerate(ID)):
            id = id.replace('EPIC ','')
            if not isinstance(id, str):
                add_empty_row(tbl)
            else:
                
                v = Vizier(column_filters={"ID":f"=={id}", 'OType':'STAR'})
                job = v.get_catalogs('IV/34/epic')
                
                ridx = job[0]['ID'].quantity == int(id)
                
                if len(job[0][ridx]) > 0:
                    tbl = avstack([tbl, job[0][ridx][0]])  
                else:
                    add_empty_row(tbl)
        
        self.EPIC = tbl
        
        # Fill in blank IDs where possible
        if not hasattr(self,'simbad'):
            self.query_simbad(ID)        
        for i in range(len(self.IDs)): 
            if len(self.IDs['2MASS'][i])==0:
                self.IDs['2MASS'][i] = '2MASS J'+self.EPIC['_2MASS'][i]
        self.query_simbad(self.IDs['2MASS'])
        
        return self.EPIC

    def query_TIC(self, ID=None, radius = 10.0*u.arcsec):
        
        if ID is None:
            ID = self.IDs['TIC']
        
        tbl = Table(names = ('ID',), dtype = (str,))
        
        for i, id in tqdm(enumerate(ID)):
            if not isinstance(id, str):
                add_empty_row(tbl)
            else:
                job = Catalogs.query_object(objectname=id, catalog='TIC', objType='STAR', radius = 10.0*u.arcsec)
                
                ridx = job['ID'] == str(id.replace('TIC ',''))
                if len(job[ridx][0]) > 0:
                    tbl = avstack([tbl, job[ridx][0]])
                else:
                    add_empty_row(tbl)
        
        self.TIC = tbl

        if not hasattr(self,'simbad'):
            self.query_simbad(ID)  
        for i in range(len(self.IDs)): 
            if len(self.IDs['2MASS'][i])==0 and (self.TIC['TWOMASS'][i] != 0):
                self.IDs['2MASS'][i] = '2MASS J'+self.TIC['TWOMASS'][i]
            
            if len(self.IDs['HIP'][i])==0 and (self.TIC['HIP'][i] != 0):
                self.IDs['HIP'][i] = 'HIP '+self.TIC['HIP'][i]
            
            if len(self.IDs['TYC'][i])==0 and (self.TIC['TYC'][i] != 0):
                self.IDs['TYC'][i] = 'TYC '+self.TIC['TYC'][i]
            
            if len(self.IDs['KIC'][i])==0 and (self.TIC['KIC'][i] != 0):
                self.IDs['KIC'][i] = 'KIC '+self.TIC['KIC'][i]
        
        return self.TIC
    
# TYC2: I/259/tyc2
# 2MASS: II/246/out
# WISE: II/311/wise
# ALLWISE: II/328/allwise


#twomassjob = Vizier.query_object(object_name = '2MASS J05065012+5236338', catalog = 'II/246/out', radius = 10.0*u.arcsec)
