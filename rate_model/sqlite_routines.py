"""
Created on Fri Mar 15 21:06:52 2019

@author: Pavel Esir
"""
import sqlite3
import hashlib
import json
import io
import numpy as np
from config import * 

def adapt_array(arr):
    """
    http://stackoverflow.com/a/31312102/190597 (SoulNibbler)
    """
    out = io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return sqlite3.Binary(out.read())

def convert_array(text):
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)

sqlite3.register_adapter(np.ndarray, adapt_array)
sqlite3.register_converter("array", convert_array)

def create_table(db_fname=db_fname, tname=db_tname):
    fname = test_results_db_path + db_fname
    with sqlite3.connect(fname, detect_types=sqlite3.PARSE_DECLTYPES) as conn:
        cur = conn.cursor()
        create_table_query = """
        CREATE TABLE {} (
            id integer PRIMARY KEY, 
            params TEXT, 
            stim_params TEXT, 
            md5_hash TEXT, 
            x array, 
            u array, 
            hE array, 
            hI array, 
            perfomed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP)
        """.format(tname)
        cur.execute(create_table_query)

def save_results(x, u, hE, hI, params, stim_params, 
                 db_fname=db_fname, tname=db_tname):
    params_str = json.dumps(params)
    stim_params_str = json.dumps(stim_params)
    md5_hash = get_params_md5(params, stim_params)
    
    fname = test_results_db_path + db_fname
    with sqlite3.connect(fname, detect_types=sqlite3.PARSE_DECLTYPES) as conn:
        cur = conn.cursor()
        insert_tbl_query = """
        INSERT INTO {} (
            params, 
            stim_params, 
            md5_hash, 
            x, 
            u, 
            hE, 
            hI)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """.format(tname)
        
        cur.execute(insert_tbl_query, 
                    (params_str, stim_params_str, md5_hash, x, u, hE, hI))

def get_results(params, stim_params, db_fname=db_fname, tname=db_tname):
    md5_hash = get_params_md5(params, stim_params)
    fname = test_results_db_path + db_fname
    with sqlite3.connect(fname, detect_types=sqlite3.PARSE_DECLTYPES) as conn:
        cur = conn.cursor()
        get_tbl_query = """ 
        SELECT x, u, hE, hI FROM {} WHERE md5_hash=(?)
        """.format(tname)
        
        cur.execute(get_tbl_query, (md5_hash, ))
        data = cur.fetchall()
    if len(data) > 1:
        raise RuntimeError("Number of rows is more than 1")
    else:
        x, u, hE, hI = data[0]
    return x, u, hE, hI

def get_params_md5(params, stim_params):
    all_params = {}
    all_params.update(params)
    all_params.update(stim_params)
    ser_string = json.dumps(all_params)
    return hashlib.md5(ser_string.encode('utf-8')).hexdigest()
