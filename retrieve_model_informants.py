"""Contains code to automatically download informants for given models
"""

from pathlib import Path
import sys
from BCBio import GFF
from Bio import SeqIO, SearchIO
from Bio.Blast import NCBIWWW
from Bio.SeqFeature import (CompoundLocation, ExactPosition, FeatureLocation,
                            SeqFeature)
from Bio.SeqRecord import SeqRecord

from find_model_diffs import collapse_record_cds


MINIMUM_COVERAGE = 0.9
NUM_MODELS = 8



def blastp_model_translation(model: SeqFeature):
    """runs a blastp with a model translation as a query

    this function assumes the the model parameter is a seq feature that
    contains a translation qualifier (these can be gotten from gbk files
    which include the translation, like those from the GEMCAPY gff converter
    or antismash for example)
    """
    translation = model.qualifiers["translation"][0]

    blast_output = NCBIWWW.qblast(
        program = "blastp",
        database = "nr_clustered",
        descriptions = 100,
        alignments = 0,
        sequence=str(translation)
    )

    return blast_output

def get_file_models(gff_file: Path):
    gff_iter = GFF.parse(gff_file, None, {"gff_type": ["CDS"]})

    models = []
    for gff_seq_rec in gff_iter:
        new_record: SeqRecord = collapse_record_cds(gff_seq_rec)
        record_features: list[SeqFeature] = new_record.features
        models.extend(record_features)

    return models


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"usage: {__file__} models.gbk")
        exit()

    gbk_file = Path(sys.argv[1])

    gbk_seq_recs = list(SeqIO.parse(gbk_file, format="genbank"))

    # maybe multiple scaffolds
    for seq_rec in gbk_seq_recs:
        seq_rec: SeqRecord

        # cds level
        for feature in seq_rec.features:
            blastp_results = blastp_model_translation(feature)
            exit()
