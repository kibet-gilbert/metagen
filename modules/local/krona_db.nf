process KRONADB {
    // debug true

    input:
    path(taxonomyDB)

    output:
    path("./taxonomy/"), emit: kronaDB

    script:
    if ( taxonomyDB.extension in ['gz', 'tgz'] )
        """
        mkdir -p ./taxonomy
        tar -xf ${taxonomyDB} -C ./taxonomy/
        ktUpdateTaxonomy.sh --only-build ./taxonomy/
        """
    else if ( taxonomyDB.isDirectory() )
        """
        mkdir -p ./taxonomy
        cp -rf ${taxonomyDB}/* ./taxonomy/
        """
    else
        error "Path or Link to a krona database not provided in ${params.krona_db}"
}

