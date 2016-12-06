## Load an EnsDb package matching Ensembl version 86
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

## Retrieve all lincRNAs encoded on chromosome Y.
## Create the filter objects
sf <- SeqnameFilter("Y")
gbf <- GenebiotypeFilter("lincRNA")

## Retrieve the data.
genes(edb, filter = list(sf, gbf))

library(AnnotationHub)
ah <- AnnotationHub()

## Query for available Ensembl gtf files for release 83.
query(ah, pattern=c("ensembl", "release-83", "gtf"))

## Select one; in this case: Anolis carolinensis (lizard)
edbSql83 <- ensDbFromAH(ah=ah["AH50353"])

## BUT: lacks NCBI Entrezgene IDs and protein annotation.

## Load the database and use it.
db <- EnsDb(edbSql83)
genes(db, filter=SeqnameFilter("2"))

## Optionally make a package.
makeEnsembldbPackage(ensdb=edbSql83, version="1.0.0", author="J Rainer",
                     maintainer="Johannes Rainer <johannes.rainer@eurac.edu>")

## Load an EnsDb from an MySQL database server.
library(RMySQL)
dbc <- dbConnect(MySQL(), host = "localhost", user = "anonuser", pass = "")

## list all available EnsDb databases.
listEnsDbs(dbc)

## Connect to a database.
dbc <- dbConnect(MySQL(), host = "localhost", user = "anonuser", pass = "",
                 dbname = "ensdb_dmelanogaster_v86")
edb_2 <- EnsDb(dbc)
edb_2

## To insert an EnsDb to a MySQL: useMySQL
db_my <- useMySQL(edb, host = "localhost", user = "anonuser", pass = "")

## Get all genes with a C2H2 Zinc finger domain and
## return all of their Uniprot IDs
pfam <- ProtdomidFilter("PF13912")

genes(edb, filter = pfam, return.type = "DataFrame",
    columns = c("gene_name", "uniprot_id"))

## Fetch all proteins for ZBTB16 return the result
## as a AAStringSet:
prts <- proteins(edb, filter = GenenameFilter("ZBTB16"),
                 columns = c("tx_id", "tx_biotype"),
                 return.type = "AAStringSet")
prts

## Additional columns are available as mcols:
mcols(prts)

## load Pbase - we need the "ensembldb" branch.
library(Pbase)

## Fetch proteins including protein domains for ZBTB16
prts <- Proteins(edb, filter = GenenameFilter("ZBTB16"))

## Amino acid sequence:
aa(prts)

## Peptide features:
pranges(prts)

## Map all protein domains from each protein/tx to the genome
gen_map <- mapToGenome(prts, edb)

## Plot the results for the first protein (transcript)
txid <- gen_map[[1]]$tx_id
## Get the gene region track for the first transcript
tx <- getGeneRegionTrackForGviz(edb, filter = TxidFilter(txid))

## Add a protein ID column
map_1 <- gen_map[[1]]
map_1$id <- names(map_1)

## Plot using Gviz
library(Gviz)
plotTracks(list(GenomeAxisTrack(),
                GeneRegionTrack(tx, name = "tx"),
                AnnotationTrack(map_1, groupAnnotation = "id",
                                just.group = "above",
                                name = "Protein domains")),
           transcriptAnnotation = "transcript")

## Get all data for the gene SKA2
Res <- select(edb, keys="SKA2", keytype="GENENAME")
head(Res, n=3)

## Or: pass filters with keys parameter to have more control:
## For the gene SKA2: get all exons except exons 1 and 2
## for all tx targeted for nonsense mediated decay.
select(edb, keys=list(GenenameFilter("SKA2"),
  		    TxbiotypeFilter("nonsense_mediated_decay"),
  		    ExonrankFilter(1:2, condition="!=")))

## Get chromosome names
head(seqlevels(edb))
## Different from UCSC style: chr1...

## Get genes on chromosome Y, UCSC style.
genes(edb, filter=SeqnameFilter("chrY"))

## Solution: change the chromosome naming style:
seqlevelsStyle(edb) <- "UCSC"

## Get chromosome names
head(seqlevels(edb))

genes(edb, filter=SeqnameFilter("chrY"))


## Use case:
## Get mRNA sequences for SKA2 using BSgenome.
library(BSgenome.Hsapiens.UCSC.hg38)  ## <- UCSC based

## Get exons by transcript
ska2tx <- exonsBy(edb, by="tx", filter=GenenameFilter("SKA2"))

## Use GenomicFeatures' extractTranscriptSeqs
head(extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, ska2tx))


## Alternative (preferred) way:
seqlevelsStyle(edb) <- "Ensembl"
## Using AnnotationHub:
## Get the genomic fasta file matching the package's genome version:
faf <- getGenomeFaFile(edb)
extractTranscriptSeqs(faf, exonsBy(edb, by="tx",
                                   filter=GenenameFilter("SKA2")))
