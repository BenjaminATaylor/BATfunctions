# A small edit to the taxizedb::db_download_ncbi() function
# We switch from ftp.ncbi.nih.gov, which sometimes fails to load, to ftp.ncbi.nlm.nih.gov
custom_db_download_ncbi = function (verbose = TRUE, overwrite = FALSE, path = NULL) 
{
  db_url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
  db_path_file <- file.path(tdb_cache$cache_path_get(), "taxdump.zip")
  db_path_dir <- file.path(tdb_cache$cache_path_get(), "taxdump")
  ncbi_names_file <- file.path(db_path_dir, "names.dmp")
  ncbi_nodes_file <- file.path(db_path_dir, "nodes.dmp")
  final_file <- file.path(tdb_cache$cache_path_get(), "NCBI.sql")
  assert(verbose, "logical")
  assert(overwrite, "logical")
  if (file.exists(final_file) && !overwrite) {
    mssg(verbose, "Database already exists, returning old file")
    return(final_file)
  }
  unlink(final_file, force = TRUE)
  tdb_cache$mkdir()
  if (is.null(path)) {
    mssg(verbose, "downloading...")
    curl::curl_download(db_url, db_path_file, quiet = TRUE)
  }
  else {
    file.copy(path, db_path_file)
  }
  mssg(verbose, "unzipping...")
  utils::unzip(db_path_file, files = c("names.dmp", "nodes.dmp"), 
               exdir = db_path_dir)
  mssg(verbose, "loading 'names.dmp'...")
  ncbi_names <- readr::read_tsv(ncbi_names_file, col_names = c("tax_id", 
                                                               "name_txt", "unique_name", "name_class"), col_type = "i_c_c_c_", 
                                quote = "")
  mssg(verbose, "loading 'nodes.dmp'...")
  ncbi_nodes <- readr::read_tsv(ncbi_nodes_file, col_names = c("tax_id", 
                                                               "parent_tax_id", "rank", "embl_code", "division_id", 
                                                               "inherited_div_flag", "genetic_code_id", "inherited_GC_flag", 
                                                               "mitochondrial_genetic_code_id", "inherited_MGC_flag", 
                                                               "GenBank_hidden_flag", "hidden_subtree_root_flag", "comments"), 
                                col_types = "i_i_c_c_i_i_i_i_i_i_i_i_c_", quote = "")
  mssg(verbose, "building hierarchy table...")
  hierarchs <- list()
  hierarchs[[1]] <- ncbi_nodes[, c("tax_id", "parent_tax_id")] %>% 
    magrittr::set_names(c("tax_id", "ancestor"))
  hierarchs[[1]]$level <- 1
  child2parent <- ncbi_nodes$parent_tax_id
  names(child2parent) <- ncbi_nodes$tax_id
  while (TRUE) {
    top <- tail(hierarchs, 1)[[1]]
    incomplete <- top$ancestor != 1L
    top <- top[incomplete, ]
    if (nrow(top) == 0) {
      break
    }
    hierarchs[[length(hierarchs) + 1]] <- tibble::tibble(tax_id = top$tax_id, 
                                                         ancestor = child2parent[as.character(top$ancestor)], 
                                                         level = rep(length(hierarchs) + 1, nrow(top)))
  }
  hierarchy <- do.call(rbind, hierarchs)
  hierarchy$level <- as.integer(hierarchy$level)
  mssg(verbose, "building SQLite database...")
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = final_file)
  RSQLite::dbExecute(conn = db, "\n    CREATE TABLE names (\n      tax_id INTEGER,\n      name_txt TEXT COLLATE NOCASE,\n      unique_name TEXT,\n      name_class TEXT\n    )\n    ")
  RSQLite::dbWriteTable(conn = db, name = "names", value = as.data.frame(ncbi_names), 
                        append = TRUE)
  RSQLite::dbWriteTable(conn = db, name = "nodes", value = as.data.frame(ncbi_nodes), 
  )
  RSQLite::dbWriteTable(conn = db, name = "hierarchy", value = as.data.frame(hierarchy), 
  )
  RSQLite::dbExecute(db, "CREATE INDEX tax_id_index_names ON names (tax_id)")
  RSQLite::dbExecute(db, "CREATE INDEX name_txt_index_names ON names (name_txt COLLATE NOCASE)")
  RSQLite::dbExecute(db, "CREATE INDEX tax_id_index_nodes ON nodes (tax_id)")
  RSQLite::dbExecute(db, "CREATE INDEX tax_id_index_hierarchy ON hierarchy (tax_id)")
  RSQLite::dbExecute(db, "CREATE INDEX tax_id_ancestor_hierarchy ON hierarchy (ancestor)")
  RSQLite::dbDisconnect(db)
  mssg(verbose, "cleaning up...")
  unlink(db_path_file)
  unlink(db_path_dir, recursive = TRUE)
  return(final_file)
}
