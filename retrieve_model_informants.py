"""Contains code to automatically download informants for given models
"""

from io import StringIO
from math import floor
from pathlib import Path
import shutil
import sys
from time import sleep
from BCBio import GFF
from Bio import SeqIO, SearchIO, Entrez
from Bio.SearchIO import QueryResult, Hit
from Bio.Blast import NCBIWWW
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from find_model_diffs import collapse_record_cds


MINIMUM_COVERAGE = 0.9
NUM_MODELS = 8


def ncbi_protein_accession_to_genome(protein_accession: str, context_len=400):
    """Retrieves the genome region associated with a given protein accession
    additional DNA upstream and downstream is selected based on context_len

    returns a list of seqrecords for the GBK associated with this portion of the genome
    """

    # get the protein details. we need to know where in the genome the protein starts and stops
    protein_handle = Entrez.efetch(
        db="protein",
        id=protein_accession,
        rettype="gb",
        retmode="text"
    )
    protein_seq_rec: list[SeqRecord] = list(SeqIO.parse(protein_handle, 'genbank'))
    protein_cds_feature = None
    for feature in protein_seq_rec[0].features:
        feature: SeqFeature
        if feature.type == "CDS":
            protein_cds_feature = feature

    if protein_cds_feature is None:
        return # could not find the cds feature. shouldn't happen

    coded_by = protein_cds_feature.qualifiers["coded_by"][0]
    # remove any spaces. there's a very fun issue where one of the start/stop coordinates
    # is split by a newline, which somehow turns into a space. since spaces aren't important,
    # we can just remove all of them
    coded_by = coded_by.replace(" ", "")

    # we need to trim the 'join()' and 'complement()' if they exist
    start_trim = 0
    end_trim = len(coded_by)
    if 'join' in coded_by:
        start_trim += 5
        end_trim -= 1


    # reverse strand? we may need to correct for this in getting start and stop if so
    if "complement" in coded_by:
        start_trim += 11
        end_trim -= 1

    join_string = coded_by[start_trim:end_trim]


    # gbks can contain other characters in their ranges. these shouldn't appear if all goes well,
    # but just to be sure we will remove them anyway
    join_string = join_string.replace(">", "")
    join_string = join_string.replace("<", "")
    join_parts = join_string.split(",")


    nucleotide_accession = None
    accession_join_parts = []
    for join_part in join_parts:
        range_accession, range = join_part.split(":")

        if nucleotide_accession is None:
            nucleotide_accession = range_accession

        if range_accession != nucleotide_accession:
            continue

        split_range = range.split("..")
        accession_join_parts.append(split_range)

    # get start and stop
    start = int(accession_join_parts[0][0])
    stop = int(accession_join_parts[-1][1])

    # in some cases we got a link to
    # link_handle = Entrez.elink(dbfrom="protein", db="nucleotide", id=protein_accession)
    # link_record = Entrez.read(link_handle, validate=False)
    # nucleotide_accession = link_record[0]["LinkSetDb"][0]["Link"][0]["Id"]


    # get the genome bit
    nucleotide_handle = Entrez.efetch(
        db="nucleotide",
        id=nucleotide_accession,
        rettype="gb",
        retmode="text",
        seq_start=str(start - context_len),
        seq_stop=str(stop + context_len)
    )
    nucleotide_seq_rec: list[SeqRecord] = list(SeqIO.parse(nucleotide_handle, 'genbank'))

    # make sure we don't overload the ncbi server. we are limited to 3 requests per second, and
    # in this function we make 2
    sleep(1)

    return nucleotide_seq_rec


def download_accessions(accessions, gbk_path_base):
    """executes ncbi_protein_accession_to_genome for each accession supplied, and skips any
    accession for which a gbk already exists
    """
    for idx, accession in enumerate(accessions):
        accession_gbk_path: Path = gbk_path_base / Path(accession + '.gbk')
        # skip if it already exists in the output directory
        if accession_gbk_path.exists():
            print("e", end="")
            continue

        # copy if it already exists in the cache
        gbk_cache_path: Path = Path("gbk_cache") / gbk_path_base.name / Path(accession + '.gbk')
        if gbk_cache_path.exists():
            shutil.copy(gbk_cache_path, accession_gbk_path)
            print("c", end="")
            continue

        nucleotide_seq_rec = ncbi_protein_accession_to_genome(accession)

        with open(accession_gbk_path, 'w', encoding='utf-8') as gbk_handle:
            SeqIO.write(nucleotide_seq_rec, accession_gbk_path, 'genbank')
        print("d", end="")

def find_best_accession(ids: list[str]):
    """Returns the first hit accession that is a non-XM_ accession

    if it cannot find one, returns the first accession"""
    for id in ids:
        parts = id.split('|')
        if parts[0] != 'ref':
            return parts[1]

    return ids[0].split('|')[1]


def select_blast_xml_accessions(xml_path, coverage_threshold=0.9, ident_threshold=None):
    """Selects a number of accessions from a blast result by:
    1. filtering the results to have >90% coverage
    2. taking the top hit and bottom hit
    3. taking every 1 in len(hits) / 10 hits

    returns a list of sequence ids which can be used to download the genome context of a protein
    using ncbi_protein_accession_to_genome
    """
    accessions = []

    with open(xml_path, encoding='utf-8') as xml_handle:
        parsed_results: list[QueryResult] = list(SearchIO.parse(xml_handle, 'blast-xml'))

    hits: list[Hit] = parsed_results[0].hits


    # 1. filter
    if coverage_threshold is not None:
        hits = list(filter(lambda hit: 1 - hit.hsps[0].gap_num/hit.hsps[0].aln_span > coverage_threshold, hits))
    if ident_threshold is not None:
        hits = list(filter(lambda hit: hit.hsps[0].ident_num/hit.hsps[0].aln_span > ident_threshold, hits))

    # 2 top and bottom
    hits = list(reversed(sorted(hits, key=lambda hit: hit.hsps[0].ident_num/hit.hsps[0].aln_span)))
    top_accession = find_best_accession(hits[0].id_all)
    bot_accession = find_best_accession(hits[-1].id_all)
    accessions.append(top_accession)
    accessions.append(bot_accession)

    # 3. take a spread of hits
    remaining_len = len(hits) - 2
    # if we now have less than or equal to 10, take all the hits
    if remaining_len <= 10:
        interval = 1
    else:
        # otherwise pick 10 hits
        interval = floor(remaining_len / 10)
    for i in range(1, len(hits) - 1, interval):
        accessions.append(find_best_accession(hits[i].id_all))

    return accessions


def extract_source_gbk(gbk_seq_recs: list[SeqRecord], target_feature_id: str, gbk_base_path: Path):
    """Retrieves the GBK entry associated with a particular feature id and writes it to a gbk
    under gbk_base_path/feature_id/feature_id.gbk

    this is inefficient, but only needs to be done once for all xmls
    """

    for seq_rec in gbk_seq_recs:
        for feature in seq_rec.features:
            # actual id
            gbk_feature_id = feature.qualifiers["ID"][0]

            # skip if not target
            if gbk_feature_id != target_feature_id:
                continue

            # folder for all gbk files related to this feature
            target_gbk_dir = gbk_base_path / Path(target_feature_id)
            target_gbk_dir.mkdir(exist_ok=True, parents=True)

            # where the actual gbk file ends up
            target_gbk_path = target_gbk_dir / Path(target_feature_id + '.gbk')

            if target_gbk_path.exists():
                return

            # generate new seq record
            new_seq_rec = SeqRecord(
                seq_rec.seq,
                id=seq_rec.id,
                name=seq_rec.name,
                description=seq_rec.description
            )
            new_seq_rec.features = [feature]
            new_seq_rec.annotations["molecule_type"] = "DNA"

            # write gbk
            gbk_handle = open(target_gbk_path, mode='w', encoding='utf-8')
            SeqIO.write(new_seq_rec, gbk_handle, 'genbank')
            gbk_handle.close()
            return

    return


def extract_xmls_source_gbk(gbk_seq_recs: list[SeqRecord], xml_dir_path: Path, gbk_base_path: Path):
    """Loops through all xml files in the xml dir path and executes source gbk extraction for each
    filename

    this should result in a number of folders under gbk_base_path, each of which contains a gbk file
    for the xml itself. Later this folder will be used to add informant gbks as well"""
    xml_files = xml_dir_path.glob("*.xml")
    for xml_file in xml_files:
        extract_source_gbk(gbk_seq_recs, xml_file.stem, gbk_base_path)


def retrieve_xmls_informants(
    xml_dir_path: Path,
    gbk_path_base: Path,
    coverage_threshold=0.9,
    ident_threshold=None
):
    """Loops through XML files in the xml folder, selects relevant accessions for each, and
    downloads the nucleotide sequence for each

    xml_dir_path: the directory in which the xml files are located
    gbk_path_base: the directory to output gbk files to. subdirectories will be created (or have
        already been created in this folder for each xml file
    email: email address to use for requests to the entrez databases
    """
    xml_files = list(xml_dir_path.glob("*.xml"))
    for idx, xml_file in enumerate(xml_files):
        print(f"XML: {idx+1}/{len(xml_files)}")
        xml_accessions = select_blast_xml_accessions(
            xml_file,
            coverage_threshold,
            ident_threshold
        )
        xml_gbk_path_base = gbk_path_base / Path(xml_file.stem)
        xml_gbk_path_base.mkdir(parents=True, exist_ok=True)
        print("|0%" + " " * (len(xml_accessions) - 8) + "100%|")
        download_accessions(xml_accessions, xml_gbk_path_base)
        print("")
    return


def blastp_model_translation(model: SeqFeature) -> StringIO:
    """runs a blastp with a model translation as a query

    this function assumes the the model parameter is a seq feature that
    contains a translation qualifier (these can be gotten from gbk files
    which include the translation, like those from the GEMCAPY gff converter
    or antismash for example)
    """
    translation = model.qualifiers["translation"][0]

    blast_output = NCBIWWW.qblast(
        program = "blastp",
        database = "nr",
        descriptions = 100,
        alignments = 100,
        sequence=str(translation)
    )

    return blast_output


def blastp_gbk_seq_recs(xml_base_path, gbk_seq_recs, include_ids=None):
    """Performs a blastp for all entries found in a gbk, creating an XML file containing hits in
    the blastp_out directory. if include_ids is passed to this function, will skip any feature ID
    that is not present in include_ids. otherwise it does everything

    this function also skips a blast query if an xml file for it already exists
    """
    count = 1
    max_seq_recs = sum([len(seq_rec.features) for seq_rec in gbk_seq_recs])
    # go through records
    for seq_rec in gbk_seq_recs:
        seq_rec: SeqRecord

        # cds level
        for feature in seq_rec.features:
            feature: SeqFeature
            if "ID" not in feature.qualifiers:
                continue
            feature_id = feature.qualifiers["ID"][0]

            xml_path = xml_base_path / Path(feature_id + ".xml")
            if xml_path.exists():
                if include_ids is not None:
                    print(f"({count}/{len(include_ids)}) XML exists: {feature_id}")
                else:
                    print(f"({count}/{len(max_seq_recs)}) XML exists: {feature_id}")
                count += 1
                continue

            if include_ids is not None and feature_id not in include_ids:
                continue

            if include_ids is not None:
                print(f"({count}/{len(include_ids)}) Running blast: {feature_id}")
            else:
                print(f"({count}/{len(max_seq_recs)}) Running blast: {feature_id}")

            blastp_results = blastp_model_translation(feature)

            write_xml(blastp_results, xml_path)

            count += 1


def get_file_models(gff_file: Path):
    gff_iter = GFF.parse(gff_file, None, {"gff_type": ["CDS"]})

    models = []
    for gff_seq_rec in gff_iter:
        new_record: SeqRecord = collapse_record_cds(gff_seq_rec)
        record_features: list[SeqFeature] = new_record.features
        models.extend(record_features)

    return models


def read_include_list(include_file: Path):
    include_ids = set()
    with open(include_file, 'r', encoding='utf-8') as include_handle:
        for line in include_handle:
            include_ids.add(line.rstrip())
    return include_ids


def write_xml(blast_results: StringIO, xml_path: Path):
    xml_path.parent.mkdir(parents=True, exist_ok=True)
    with open(xml_path, 'w', encoding='utf-8') as xml_handle:
        xml_handle.write(blast_results.read())


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"usage: {__file__} <models.gbk> <gbk output path> <your@email.com> [include models.txt] [coverage threshold] [identity threshold]")
        exit()

    gbk_file = Path(sys.argv[1])
    if not gbk_file.exists():
        print(f"{gbk_file} does not exist")
        exit()

    gbk_base_path = Path(sys.argv[2])
    if gbk_base_path.exists() and gbk_base_path.is_file():
        print(f"{gbk_file} exists and is a file!")
        exit()

    gbk_base_path.mkdir(parents=True, exist_ok=True)

    Entrez.email = sys.argv[3]

    include_ids = None
    if len(sys.argv) > 4:
        include_file = Path(sys.argv[4])
        if not include_file.exists():
            print(f"{include_file} does not exist")
            exit()
        include_ids = read_include_list(include_file)

    coverage_treshold = 0.9
    if len(sys.argv) > 5:
        coverage_treshold = float(sys.argv[5])

    ident_threshold = None
    if len(sys.argv) > 6:
        ident_threshold = float(sys.argv[6])

    gbk_seq_recs = list(SeqIO.parse(gbk_file, format="genbank"))

    # perform blasts on whatever is relevant
    xml_base_path = Path('blastp_out')
    # blastp_gbk_seq_recs(xml_base_path, gbk_seq_recs, include_ids)

    extract_xmls_source_gbk(gbk_seq_recs, xml_base_path, gbk_base_path)
    retrieve_xmls_informants(xml_base_path, gbk_base_path, coverage_treshold, ident_threshold)
