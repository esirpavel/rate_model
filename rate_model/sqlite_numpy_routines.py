"""
Created on Fri Mar 15 21:06:52 2019

@author: Pavel Esir
"""
import sqlite3
import hashlib
import json
import random
import string
from pathlib import Path
import numpy as np

def get_md5(str):
    return hashlib.md5(str.encode('utf-8')).hexdigest()

def generate_rand_str(N=20):
    rng = random.SystemRandom()
    return ''.join(rng.choice(list(string.ascii_lowercase +
        string.digits)) for _ in range(N))

class SqlWrapper():

    def __init__(self, root_folder, db_fname, db_tname):
       self.db_fname = Path(db_fname)
       self.db_tname = db_tname
       self.root_folder = Path(root_folder)
       self.create_table()

    def create_table(self):
        with sqlite3.connect(self.root_folder / self.db_fname) as conn:
            cur = conn.cursor()
            create_table_query = """
            CREATE TABLE IF NOT EXISTS {} (
                id integer PRIMARY KEY, 
                params TEXT, 
                stim_params TEXT, 
                noise_params TEXT,
                comments TEXT,
                md5_hash TEXT, 
                file_name TEXT,
                performed_at TDATETIME DEFAULT (DATETIME(CURRENT_TIMESTAMP, 'LOCALTIME')))
            """.format(self.db_tname)
            cur.execute(create_table_query)

    def save_results(self, x, u, hE, hI, params, stim_params, noise_params, comments=''):
        params_str = json.dumps(params)
        stim_params_str = json.dumps(stim_params)
        noise_params_str = json.dumps(noise_params)
        
        Path(self.root_folder / 'numpy_data').mkdir(exist_ok=True)
        file_name = 'numpy_data/' + generate_rand_str() + '.npz'
        
        with sqlite3.connect(self.root_folder / self.db_fname) as conn:
            cur = conn.cursor()
            insert_tbl_query = """
            INSERT INTO {} (
                params, 
                stim_params, 
                noise_params,
                comments,
                md5_hash,
                file_name)
            VALUES (?, ?, ?, ?, ?, ?)
            """.format(self.db_tname)
            
            args = (
                params_str, 
                stim_params_str, 
                noise_params_str,
                comments,
                get_md5(params_str + stim_params_str + noise_params_str +
                    comments),
                file_name
            )
            cur.execute(insert_tbl_query, args)
        
        np.savez_compressed(self.root_folder / file_name, 
                x=x, u=u, hE=hE, hI=hI)

    def get_results(self, params, stim_params, noise_params, comments=''):
        params_str = json.dumps(params)
        stim_params_str = json.dumps(stim_params)
        noise_params_str = json.dumps(noise_params)
        md5_hash = get_md5(params_str + stim_params_str + noise_params_str +
                comments)
        
        with sqlite3.connect(self.root_folder / self.db_fname) as conn:
            cur = conn.cursor()
            get_tbl_query = """ 
            SELECT file_name FROM {} WHERE md5_hash=(?)
            """.format(self.db_tname)
            
            cur.execute(get_tbl_query, (md5_hash, ))
            fetched_fnames = cur.fetchall()
        if len(fetched_fnames) > 1:
            raise RuntimeError("Multiple items with these parameters")
        elif len(fetched_fnames) == 0:
            raise RuntimeError("Simulations with specified parameters are not " 
                "found")
        
        with np.load(self.root_folder / fetched_fnames[0][0]) as data:
            x = data['x']
            u = data['u']
            hE = data['hE']
            hI = data['hI']
            return x, u, hE, hI

