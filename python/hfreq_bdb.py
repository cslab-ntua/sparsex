
""" Berkeley DB frequency hash """

from bsddb.db import DB, DB_HASH, DB_CREATE
from os import makedirs
from os.path import dirname, isdir
from struct import pack, unpack

def hfreq_bdb_init(db_file):
	""" initialize (open) a bdb frequency hash """
	db_dir = dirname(db_file)
	if not isdir(db_dir):
		makedirs(db_dir)
	db = DB()
	db.open(db_file, dbtype=DB_HASH, flags=DB_CREATE)
	return db


def hfreq_bdb_add(db, v, freq=1):
	try:
		db[v] = pack("L", unpack("L", db[v])[0] + long(freq))
	except KeyError:
		db[v] = pack("L", long(freq))

def hfreq_bdb_get(db, v):
	return unpack("L", db[v])[0]

def bdb_iteritems(db):
	c = db.cursor()
	while True:
		item = c.next()
		if item is None:
			break
		yield (item[0], unpack("L", item[1])[0])
