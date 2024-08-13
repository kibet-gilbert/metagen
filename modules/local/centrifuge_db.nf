process CENTRIFUGEDB {
    // debug true
    
    input:
    tuple val(centrifugeDBName), path(centrifuge_db)
    
    output:
    tuple val(centrifugeDBName), path("database/*.cf"), emit: centrifuge_db

    script:
    if ( centrifuge_db.extension in ['gz', 'tgz'] )
	"""
	mkdir -p ./{database,db_tmp}
	tar -xf ${centrifuge_db} -C db_tmp
	mv `find db_tmp/ -name "*.cf"` database/
	"""
    else if ( centrifuge_db.isDirectory() ) 
	"""
	mkdir -p database
	cp `find ${centrifuge_db}/ -name "*.cf"` database/
	"""
    else
	error "Path or Link to a kraken2 database not provided in ${params.centrifuge_db}"
}
