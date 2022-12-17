import pymysql
from datetime import datetime, timedelta
import math

# https://veritas.sao.arizona.edu/wiki/VOFFLINE_Database_Tables
# https://veritas.sao.arizona.edu/wiki/Ryan_Dickherber%27s_Wiki

def get_sec(time_str):
    """Get Seconds from time."""
    if time_str=='None': return 0
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + int(s)

# setup database connection
#dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)

# connect to database
crs=dbcnx.cursor()

# define query
# VOFFLINE
#query = 'SELECT * FROM tblRun_Analysis_Comments WHERE run_id=%s'
#query = 'SELECT run_id,time_cut_mask FROM tblRun_Analysis_Comments'
#query = 'SELECT run_id,data_category FROM tblRun_Analysis_Comments'
#query = 'SELECT run_id,usable_duration FROM tblRun_Analysis_Comments'
# VERITAS
#query = 'SELECT run_id,source_id,offsetRA,offsetDEC,offset_angle,weather FROM tblRun_Info WHERE run_id=%s'
query = 'SELECT run_id,source_id,offsetRA,offsetDEC,offset_angle,weather FROM tblRun_Info'

# get data to cursor
#crs.execute(query, '57156')
crs.execute(query)

# fetch from cursor
res = crs.fetchall()
for x in res:
    #print(x['run_id'],x['time_cut_mask'])
    #print(x['run_id'],x['data_category'])
    #print(x['run_id'],x['usable_duration'])
    #print(x['run_id'],get_sec(str(x['usable_duration'])))

    #print(x['run_id'],x['source_id'],x['weather'])
    source_name = x['source_id']
    if x['run_id']<46642: continue
    if x['weather']==None: continue
    if x['source_id']==None: continue
    if x['source_id']=='engineering': continue
    if x['source_id']=='other': continue
    if x['source_id']=='none': continue
    if 'C' in x['weather']: continue
    if 'D' in x['weather']: continue
    if 'F' in x['weather']: continue
    query2 = 'SELECT source_id,ra,decl FROM tblObserving_Sources WHERE source_id=%s'
    crs.execute(query2, source_name)
    res2 = crs.fetchall()
    print(x['run_id'],(x['offsetRA']+res2[0]['ra'])*180./math.pi,(x['offsetDEC']+res2[0]['decl'])*180./math.pi)

