import astropy.units as u
import time
import numpy as np
import pandas as pd
from astropy.table import Table 
from astropy.table import vstack as avstack
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
Simbad.add_votable_fields('sptype', 'ids')

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
        name = name.replace('KIC','KIC ')
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

class search():
     
    def __init__(self, ID):
        
        self.cats_avail =  ['SPOCS', 'KIC', 'TIC', 'EPIC', 'Gaia DR1', 'Gaia DR2', 'HD', 
                            'HIP', 'HR', 'HIC', 'UBV', 'SAO', 'GEN#', 'TYC', '2MASS', 
                            'GJ', 'PPM', 'BD', 'AG', 'GSC', 'Ci', 'PLX', 'SKY#', 'WISEA', 
                            'WISE', 'PSO']
               
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
            
                simIDS = [str(x) for x in self.simbad['IDS'] if inp in str(x)] # Find tgt simbad entry
                
                if len(simIDS) == 0: # Skip target if no entry found
                    continue
                else:
                    simIDS = simIDS[0].strip("b'").split('|')
                    if (len(self.simbad['IDS']) == 0) or (len(self.simbad['IDS']) > 100):
                        raise ValueError('Something went wrong trying to split the simbad IDS result')
 
                    for k, id in enumerate(simIDS): # Assign names to ID dataframe
                        for j, cat in enumerate(self.cats_avail):
                            if cat in id:
                                self.IDs.loc[i, cat] = simIDS[k]
    
        return self.simbad
    
    def query_KIC(self, ID = None, radius = 10.0*u.arcsec):
        import warnings
        from astropy.utils.metadata import MergeConflictWarning
        warnings.filterwarnings("ignore", category = MergeConflictWarning)
        
        if not ID:
            ID = a.IDs['KIC']
                
        tbl = Table(names = ('KIC',), dtype = (int,))
        for i, id in enumerate(ID):
            if not isinstance(id, str):
                add_empty_row(tbl)
            else:
                job = Vizier.query_object(object_name = id, catalog = 'V/133/kic', radius = 10.0*u.arcsec)

                if len(job) > 0:
                    tbl = avstack([tbl, job[0]])      
                else:
                    add_empty_row(tbl)
                    
        self.KIC = tbl
        return self.KIC
    
    def query_GaiaDR2(self, ID=None):
        from astroquery.gaia import Gaia

        if not ID:
            ID = a.IDs['Gaia DR2']

        tbl = Table(names = ('source_id',), dtype = (int,))
        for i, gid in enumerate(ID):
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


#job = Vizier.query_object(object_name = self.KIC['ID'], catalog = 'IV/34/epic', radius = radius)
#from astroquery.mast import Catalogs
#ticjob = Catalogs.query_object(objectname=self.TIC['ID'], catalog='TIC', objType='STAR', radius = radius)
