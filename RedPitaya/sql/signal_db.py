"""Functions compatible with the output of the RP for writing to an SQLite database
""" 

import sqlite3
import numpy as np
import time
import io

def adapt_array(arr):
    out = io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return sqlite3.Binary(out.read())

def convert_array(text):
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)

def create_db(dbname='./signaldb'):    
    """ Create a database"""    
    sqlite3.register_adapter(np.ndarray, adapt_array)
    sqlite3.register_converter("array", convert_array)
    db=sqlite3.connect(dbname,detect_types=sqlite3.PARSE_DECLTYPES)    
    cursor=db.cursor()    
    cursor.execute('''
        CREATE TABLE signal(time FLOAT PRIMARY KEY, macadd TEXT,
                            data ARRAY)
    ''')
    db.close()
    
def write_to_db(macaddress,arr,dbname='./signaldb',t=time.time()):
    """ Writes basic data to the signal table in a sqlite db."""
    db=sqlite3.connect(dbname,detect_types=sqlite3.PARSE_DECLTYPES)

    cursor=db.cursor()
    
    cursor.execute('''INSERT INTO signal(time, macadd, data) 
    VALUES(?,?,?)''',(t,macaddress,arr))
    
    db.close()