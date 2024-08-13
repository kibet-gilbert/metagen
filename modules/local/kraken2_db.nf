process KRAKEN2DB {
    // debug true
    
    input:
    path(kraken2_db)
    
    output:
    path("database/*.k2d"), emit: kraken2_db

    script:
    if ( kraken2_db.extension in ['gz', 'tgz'] )
	"""
	mkdir -p ./{database,db_tmp}
	tar -xf ${kraken2_db} -C db_tmp
	mv `find db_tmp/ -name "*.k2d"` database/
	"""
    else if ( kraken2_db.isDirectory() ) 
	"""
	mkdir -p database
	cp `find ${kraken2_db}/ -name "*.k2d"` database/
	"""
    else
	error "Path or Link to a kraken2 database not provided in ${params.kraken2_db}"
}

