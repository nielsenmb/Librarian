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
import re, time, warnings
from astropy.utils.metadata import MergeConflictWarning

def kic_to_kplr(x):
    y = x.strip('KIC')
    y = 'kplr'+'0'*(9-len(y))+y
    return y

def kplr_to_kic(x):
    y = x.strip('kplr')
    return(str(int(y)))

def format_name(name):
    # Add naming exceptions here
    name = str(name)
    name = name.lower()
        
    variants = {'KIC': ['kic', 'kplr'],
                'Gaia DR3': ['gaia dr3', 'gdr3', 'dr3'],
                'Gaia DR2': ['gaia dr2', 'gdr2', 'dr2'],
                'Gaia DR1': ['gaia dr1', 'gdr1', 'dr1'], 
                'EPIC': ['epic', 'ktwo'],
                'TIC': ['tic']
               }
    fname = None
    for key in variants:   
        for x in variants[key]:
            if x in name:
                fname = name.replace(x,'')
                fname = re.sub(r"\s+", "", fname, flags=re.UNICODE)
                fname = key+' '+fname
                return fname
    
    if fname is None:
        warnings.warn(f'Unable to format target {name} for Simbad query.')

#def add_empty_row(tbl):
    #s.simbad.add_row(None)
    #s.simbad[-1]['MAIN_ID'] = 2
    #tbl.add_row([0]*len(tbl.keys()), mask=[True]*len(tbl.keys()))

def add_to_table(job, tbl, identifier):
    if len(job) == 0:
        add_empty_row(tbl)
    else:
        for i, row in enumerate(job[0]):
            if row[identifier] in id:
                tbl = avstack([tbl, row])
                break
            
def bytes_to_str(x):
    if not isinstance(x, (list, np.ndarray)):
        x = [x]
    return np.array(x).astype(str)

class search():
     
    def __init__(self, ID):
        
        self.cats_avail =  ['SPOCS', 'KIC', 'TIC', 'EPIC', 'Gaia DR1', 'Gaia DR2', 
                            'Gaia DR3', 'HD', 'HIP', 'HR', 'HIC', 'UBV', 'SAO', 
                            'GEN#', 'TYC', '2MASS',  'GJ', 'PPM', 'BD', 'AG', 
                            'GSC', 'Ci', 'PLX', 'SKY#', 'WISEA', 'WISE', 'PSO', 
                            'ALLWISE']
        
        if isinstance(ID, pd.Series):
            ID = ID.values

        elif not isinstance(ID, (list, tuple, np.ndarray, pd.Series)):       
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
        self.query_TIC() 
        self.query_GaiaDR3()


           
    def query_simbad(self, ID=None):
        time.sleep(1) # This is to keep you from being temporarily blacklisted by Simbad servers       

        if ID is None:
            ID = self.IDs.input.values
        
        job = Simbad.query_objects(ID)

        if job is None:
            print('Unable to find any of the targets with Simbad')      
            self.simbad = None
            return self.simbad
        
        # Create empty copy of query job to edit and parse.
        tbl = job.copy()
        tbl.remove_rows(np.arange(len(tbl)))
        
        # Loop through targets in the returned query
        for j, y in enumerate(ID):
            for i, jobIDS in enumerate(job.iterrows('IDS')):
                jobIDS = np.array(jobIDS).astype(str)[0]
            
                # Loop through requested targets and compare with returned query.
            
                if y in jobIDS:
                    tbl.add_row(job[i])
                    jobIDSlist = jobIDS.strip("b'").split('|')
                    if (len(jobIDSlist) == 0) or (len(jobIDSlist) > 100):
                        print(jobIDSlist)
                        raise ValueError('Something went wrong trying to split the simbad IDS result')
                    else:
                        for jobID in jobIDSlist: # Assign names to ID dataframe
                            for cat in self.cats_avail:
                                if (cat in jobID) and (len(self.IDs.loc[j, cat])==0):
                                    self.IDs.loc[j, cat] = jobID
                    break # ensure only 1 target is added. 
            
            if (len(tbl)-1 < j):
                tbl.add_row(None)
                tbl[-1]['MAIN_ID'] = y
                    
        self.simbad = tbl
        
        return self.simbad        
    
    
    def query_KIC(self, ID = None, radius = 10.0*u.arcsec):
    
        warnings.filterwarnings("ignore", category = MergeConflictWarning)
    
        if ID is None:
            ID = self.IDs['KIC']
    
        tbl = Table(names = ('KIC',), dtype = (int,))
        for i, id in tqdm(enumerate(ID)):
            
            if not isinstance(id, ):
                tbl.add_row(None)
                tbl[-1]['KIC'] = int(re.sub(r"\D", "", id))
                
            else:
                job = Vizier.query_object(object_name = id, catalog = 'V/133/kic', radius = radius)
                
                if len(job) > 0:
                    
                    id_int = int(id.replace('KIC',''))
                
                    idx = job[0]['KIC'] == id_int
                    
                    id_int = int(id.replace('KIC',''))
                
                    idx = job[0]['KIC'] == id_int
                    
                    tbl = avstack([tbl, job[0][idx]]) 
                
                else:
                    tbl.add_row(None)
                    tbl[-1]['KIC'] = int(re.sub(r"\D", "", id))
                    warnings.warn(f'Unable to find KIC entry for {id}')
    
        self.KIC = tbl
        return self.KIC
    
    def query_GaiaDR3(self, ID=None):
        from astroquery.gaia import Gaia
        
        key = 'Gaia DR3'

        if ID is None:
            ID = self.IDs[key]
    
        tbl = Table(names = ('source_id',), dtype = (int,))

        for i, gid in tqdm(enumerate(ID)):
            print(gid)

            if not isinstance(gid, str):
                tbl.add_row(None)

                tbl[-1][key] = int(re.sub(r"\D", "", gid))

            elif len(gid) == 0:
                tbl.add_row(None)

            else:
                gid = int(gid.replace(key+' ', ''))
    
                adql_query = "select * from gaiadr3.gaia_source where source_id=%i" % (gid)
                
                job = Gaia.launch_job(adql_query).get_results()
    
                idx = np.where(job['source_id'].quantity == gid)[0]
    
                if len(idx) > 0:
                    tbl = avstack([tbl, job[idx]])
                else:
                    tbl.add_row(None)
                    tbl[-1][key] = int(re.sub(r"\D", "", gid))
    
        self.GDR3 = tbl
        return self.GDR3
    
    def query_GaiaDR2(self, ID=None):
        from astroquery.gaia import Gaia
        
        key = 'Gaia DR2'

        if ID is None:
            ID = self.IDs[key]
    
        tbl = Table(names = ('source_id',), dtype = (int,))
        for i, gid in tqdm(enumerate(ID)):
            print(gid)
            if not isinstance(gid, str):
                tbl.add_row(None)
                tbl[-1][key] = int(re.sub(r"\D", "", gid))
            elif len(gid) == 0:
                tbl.add_row(None)
            else:
                gid = int(gid.replace(key+' ', ''))
    
                adql_query = "select * from gaiadr2.gaia_source where source_id=%i" % (gid)
                
                job = Gaia.launch_job(adql_query).get_results()
    
                idx = np.where(job['source_id'].quantity == gid)[0]
    
                if len(idx) > 0:
                    tbl = avstack([tbl, job[idx]])
                else:
                    tbl.add_row(None)
                    tbl[-1][key] = int(re.sub(r"\D", "", gid))
    
        self.GDR2 = tbl
        return self.GDR2

    def query_GaiaDR1(self, ID=None):
        from astroquery.gaia import Gaia
        
        key = 'Gaia DR1'
        
        if ID is None:
            ID = self.IDs[key]
    
        tbl = Table(names = ('source_id',), dtype = (int,))
        for i, gid in tqdm(enumerate(ID)):
            if not isinstance(gid, str):
                tbl.add_row(None)
                tbl[-1][key] = int(re.sub(r"\D", "", gid))
            elif len(gid) == 0:
                tbl.add_row(None)
            else:
                gid = int(gid.replace(key+' ', ''))
    
                adql_query = "select * from gaiadr1.gaia_source where source_id=%i" % (gid)
                
                job = Gaia.launch_job(adql_query).get_results()
    
                idx = np.where(job['source_id'].quantity == gid)[0]
    
                if len(idx) > 0:
                    tbl = avstack([tbl, job[idx]])
                else:
                    tbl.add_row(None)
                    tbl[-1][key] = int(re.sub(r"\D", "", gid))
    
        self.GDR1 = tbl
        return self.GDR1
    
    
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
        
        key = 'TIC'
        
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

        if not hasattr(self, 'simbad'):
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
