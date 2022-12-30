"""Contains code to automatically download informants for given models
"""

from io import StringIO
from math import floor
from pathlib import Path
import sys
from BCBio import GFF
from Bio import SeqIO, SearchIO, Entrez
from Bio.SearchIO import QueryResult, Hit
from Bio.Blast import NCBIWWW
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from find_model_diffs import collapse_record_cds


MINIMUM_COVERAGE = 0.9
NUM_MODELS = 8


def ncbi_protein_accession_to_genome(accession: str, context_len=100):
    """Retrieves the genome region associated with a given protein accession
    additional DNA upstream and downstream is selected based on context_len"""

    # get the protein details. we need to know where in the genome the protein starts and stops
    protein_handle = Entrez.efetch(
        db="protein",
        id=accession,
        rettype="gb",
        retmode="text"
    )
    protein_seq_rec = SeqIO.parse(protein_handle, 'gbk')

    # get the link to the genome
    link_handle = Entrez.elink(dbfrom="protein", db="nucleotide", id=accession)
    link_record = Entrez.read(link_handle, validate=False)

    nucleotide_accession = link_record[0]["LinkSetDb"][0]["Link"][0]["Id"]

    # get the genome bit
    nucleotide_handle = Entrez.efetch(
        db="nucleotide",
        id=nucleotide_accession,
        rettype="gb",
        retmode="text"
    )

    return

def download_accessions(accessions, gbk_path_base):
    for accession in accessions:
        ncbi_protein_accession_to_genome(accession)


def select_blast_xml_accessions(xml_path):
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

    hsps: list[Hit] = parsed_results[0].hsps


    # 1. filter
    hsps = list(filter(lambda hsp: 1 - hsp.gap_num/hsp.aln_span > 0.9, hsps))

    # 2 top and bottom
    hsps = list(reversed(sorted(hsps, key=lambda hsp: hsp.ident_num/hsp.aln_span)))
    top_accession = hsps[0].hit_id.split("|")[1]
    bot_accession = hsps[-1].hit_id.split("|")[1]
    accessions.append(top_accession)
    accessions.append(bot_accession)

    # 3. take a spread of hits
    remaining_len = len(hsps) - 2
    interval = floor(remaining_len / 10)
    for i in range(1, len(hsps) - 1, interval):
        accessions.append(hsps[i].hit_id.split("|")[1])

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
            new_seq_rec = SeqRecord(seq_rec.seq, id=seq_rec.id, name=seq_rec.name, description=seq_rec.description)
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


def process_xmls(xml_dir_path: Path, gbk_path_base: Path):
    """Loops through XML files in the xml folder, selects relevant accessions for each, and
    downloads the nucleotide sequence for each

    xml_dir_path: the directory in which the xml files are located
    gbk_path: the directory to output gbk files to. subdirectories will be created in this folder
    for each xml file
    email: email address to use for requests to the entrez databases
    """
    xml_files = xml_dir_path.glob("*.xml")
    for xml_file in xml_files:
        xml_accessions = select_blast_xml_accessions(xml_file)
        download_accessions(xml_accessions, gbk_path_base)
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
    if len(sys.argv) < 3:
        print(f"usage: {__file__} <models.gbk> <include models.txt>")
        exit()

    gbk_file = Path(sys.argv[1])
    if not gbk_file.exists():
        print(f"{gbk_file} does not exist")
        exit()

    Entrez.email = sys.argv[2]

    include_ids = None
    if len(sys.argv) > 3:
        include_file = Path(sys.argv[3])
        if not include_file.exists():
            print(f"{include_file} does not exist")
            exit()
        include_ids = read_include_list(include_file)

    gbk_seq_recs = list(SeqIO.parse(gbk_file, format="genbank"))

    # perform blasts on whatever is relevant
    xml_base_path = Path('blastp_out')
    # blastp_gbk_seq_recs(xml_base_path, gbk_seq_recs, include_ids)

    gbk_base_path = Path('gbk_out')
    extract_xmls_source_gbk(gbk_seq_recs, xml_base_path, gbk_base_path)
    # process_xmls(xml_base_path, gbk_path_base)

