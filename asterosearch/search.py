import astropy.units as u
import numpy as np
import pandas as pd
from astropy.table import Table 
from astropy.table import vstack as avstack
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
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
    if 'KIC' in name:
        name = name.replace('KIC','')
        name = re.sub(r"\s+", "", name, flags=re.UNICODE)
        name = 'KIC '+name
    
    if 'EPIC' in name:
        name = name.replace('EPIC','')
        name = re.sub(r"\s+", "", name, flags=re.UNICODE)
        name = 'EPIC '+name
        
    if 'GRD2' in name:
        name = name.replace('GDR2','Gaia DR2 ')
    if 'GRD1' in name:
        name = name.replace('GDR1','Gaia DR1 ')
    return name

def add_empty_row(tbl):
    if len(tbl) == 0:
        tbl.add_row([0], mask=[True])
    else:
        tbl.add_row([0]*len(tbl[-1]), mask = [True]*len(tbl[-1]))

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
                
        for i in self.IDs.index:
            self.IDs.input[i] = format_name(self.IDs.input[i])
            for cat in self.cats_avail:
                if cat in self.IDs.input[i]:                    
                    self.IDs.loc[i, cat] = self.IDs.input[i]

    def __call__(self):
        self.query_simbad()
        self.query_KIC()
        #self.query_TIC() 
        self.query_GaiaDR2()
        
    def query_simbad(self):
        time.sleep(1) # This is to keep you from being blacklisted by Simbad servers       
          
        self.simbad = Simbad.query_objects(self.IDs.input)
                   
        if not self.simbad:
            print('Unable to find any of the targets with Simbad')        
        else:            
            for i in self.IDs.index: # Loop through query targets
                inp = self.IDs.input[i]
            
                simID = [str(x) for x in self.simbad['IDS'] if inp in str(x)] # Find tgt simbad entry
                
                if len(simID) == 0: # Skip target if no entry found
                    continue
                else:
                    simID = simID[0].strip("b'").split('|')
                    if (len(simID) == 0) or (len(simID) > 100):
                        raise ValueError('Something went wrong trying to split the simbad IDS result')
 
                    for k, id in enumerate(simID): # Assign names to ID dataframe
                        for j, cat in enumerate(self.cats_avail):
                            if cat in id:
                                self.IDs.loc[i, cat] = simID[k]
    
        return self.simbad
    
    def query_KIC(self, ID = None, radius = 10.0*u.arcsec):
        import warnings
        from astropy.utils.metadata import MergeConflictWarning
        warnings.filterwarnings("ignore", category = MergeConflictWarning)
        
        if not ID:
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

        if not ID:
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
    
    def query_EPIC(self, ID=None, radius = 10.0*u.arcsec):
        import warnings
        from astropy.utils.metadata import MergeConflictWarning
        warnings.filterwarnings("ignore", category = MergeConflictWarning)
        
        if not ID:
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


#                elif len(job) == 1:
#                    tbl = avstack([tbl, job[0][ridx]])      
#                    
        self.EPIC = tbl
        return self.EPIC

# TYC2: I/259/tyc2
# 2MASS: II/246/out
# WISE: II/311/wise
# ALLWISE: II/328/allwise


#from astroquery.mast import Catalogs
#ticjob = Catalogs.query_object(objectname=self.TIC['ID'], catalog='TIC', objType='STAR', radius = radius)
