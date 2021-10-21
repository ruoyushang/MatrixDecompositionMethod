import pymysql
from datetime import datetime, timedelta

# https://veritas.sao.arizona.edu/wiki/VOFFLINE_Database_Tables

def get_sec(time_str):
    """Get Seconds from time."""
    if time_str=='None': return 0
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + int(s)

# setup database connection
dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)

# connect to database
crs=dbcnx.cursor()

# define query
#query = 'SELECT * FROM tblRun_Analysis_Comments WHERE run_id=%s'
#query = 'SELECT run_id,time_cut_mask FROM tblRun_Analysis_Comments'
#query = 'SELECT run_id,data_category FROM tblRun_Analysis_Comments'
query = 'SELECT run_id,usable_duration FROM tblRun_Analysis_Comments'

# get data to cursor
#crs.execute(query, '88738')
crs.execute(query)

# fetch from cursor
res = crs.fetchall()
for x in res:
    #print(x['run_id'],x['time_cut_mask'])
    #print(x['run_id'],x['data_category'])
    print(x['run_id'],get_sec(str(x['usable_duration'])))

