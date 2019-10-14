import astropy.units as u
import time
import numpy as np

class search():
     
    def __init__(self, ID):
        
        if not isinstance(ID, str):
            TypeError('ID must be a string which is resolvable by Simbad')
            return None
        
        self.input_ID = ID
        self.KIC = {'ID': None, 'entry': None}
        self.EPIC = {'ID': None, 'entry': None}
        self.TIC = {'ID': None, 'entry': None}
        self.GDR2 = {'ID': None, 'entry': None}

    def __call__(self):
        self.query_simbad()
        self.query_KIC()
        #self.query_TIC() 
        self.query_GaiaDR2()
        
    def query_simbad(self):
        from astroquery.simbad import Simbad
        Simbad.add_votable_fields('sptype', 'ids')
        time.sleep(1) # This is to keep you from being blacklisted by Simbad servers       
        self.simbad = {'ID': None, 'entry': None}
        
        try:
            self.simbad['entry'] = Simbad.query_object(self.input_ID)
        except Exception as ex:
            message = "An exception of type {0} occurred when querying Simbad for {1}. Arguments:\n{2!r}".format(type(ex).__name__, self.input_ID, ex.args)
            print(message)
                   
        if not self.simbad['entry']:
            print(f'Unable to find {self.input_ID} with Simbad')        
        else:
            IDS = self.simbad['entry']['IDS']
            
            self.simbad['ID'] = str(IDS[0]).strip("b'").split('|') 
            if (len(self.simbad['ID']) == 0) or (len(self.simbad['ID']) > 40):
                raise ValueError('Something went wrong trying to split the simbad IDS result')
            
            prefixes = ['Gaia DR2', 'KIC', 'EPIC', 'TIC']
            #cats = ['KIC', 'TIC', 'EPIC', 'Gaia DR1', 'Gaia DR2', 'HD', 'HIP', 'HR', 'TYC']
            
            
            dicts = [self.GDR2, self.KIC, self.EPIC, self.TIC]
            for i, prefix in enumerate(prefixes):

                tempid = [id for id in self.simbad['ID'] if prefix in id]
                if len(tempid) == 1:
                    dicts[i]['ID'] = tempid[0]
                elif len(tempid) > 1:
                    print(f'Warning: for some reason Simbad returned more than one {prefix} entry for {self.input_ID}.')
                else:
                    pass

        return self.simbad['entry']
        
    def query_KIC(self, radius = 10.0*u.arcsec):
        from astroquery.vizier import Vizier
        time.sleep(1)
        kicjob = Vizier.query_object(object_name = self.KIC['ID'], catalog = 'V/133/kic', radius = radius)
        self.KIC['entry'] = kicjob
        return self.KIC['ID'], self.KIC['entry']
    
    def query_EPIC(self, radius = 10.0*u.arcsec):
        from astroquery.vizier import Vizier
        time.sleep(1)
        job = Vizier.query_object(object_name = self.KIC['ID'], catalog = 'IV/34/epic', radius = radius)
        self.EPIC['entry'] = job
        return self.EPIC['ID'], self.EPIC['entry']
    
    def query_TIC(self, radius = 10.0*u.arcsec):
        from astroquery.mast import Catalogs
        time.sleep(1)
        ticjob = Catalogs.query_object(objectname=self.TIC['ID'], catalog='TIC', objType='STAR', radius = radius)
        self.TIC['entry'] = ticjob
        return self.TIC['ID'], self.TIC['entry']
    
    def query_GaiaDR2(self):
        if not self.GDR2['ID']:
            
            if 'Gaia DR2' in self.input_ID:
                self.GDR2['ID'] = self.input_ID
            else:
                print(f'To query the Gaia DR2 catalog, {self.input_ID} must be a valid Gaia source id.')
                print('Try using the archives.query_simbad() method.')            
        
        if self.GDR2['ID']:
            from astroquery.gaia import Gaia
            time.sleep(1)

            gid = int(self.GDR2['ID'].replace('Gaia DR2 ', ''))
            adql_query = "select * from gaiadr2.gaia_source where source_id=%i" % (gid)

            job = Gaia.launch_job(adql_query).get_results()
            idx = np.where(job['source_id'].quantity == gid)[0][0]
            self.GDR2['entry'] = job[idx]

        return self.GDR2['ID'], self.GDR2['entry']