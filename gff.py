"""Contains code to convert gff files to gbk files

Authors:
Arjan Draisma
Jorge Navarro
Jérôme Collemare
"""

import os
from pathlib import Path
import sys

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation
from Bio.SeqRecord import SeqRecord


def get_exons(feature: SeqFeature):
    """Retrieves the single feature location or a compound location representing the exons in a CDS
    """
    if len(feature.sub_features) == 0:
        return feature.location

    location_parts = []
    sub_features = feature.sub_features
    for subfeature in sub_features if feature.strand == 1 else reversed(sub_features):
        subfeature: SeqFeature
        location_parts.append(subfeature.location)
    return CompoundLocation(location_parts)

def add_translations(record: SeqRecord):
    """Adds translation qualifiers to each feature of a record"""
    for feature in record.features:
        translation = feature.translate(record.seq, cds=False)
        feature.qualifiers["translation"] = translation

def remove_stop_codons(record: SeqRecord):
    """Removes stop codons and updates feature locations if one was detected"""
    for feature in record.features:
        if feature.qualifiers["translation"][-1:] != "*":
            continue
        feature: SeqFeature
        feature.qualifiers["translation"] = feature.qualifiers["translation"][:-1]
        if feature.strand == 1:
            if len(feature.location.parts) == 1:
                feature.location = FeatureLocation(feature.location.start, feature.location.end - 3, 1)
                continue
            new_locations = []
            # add old
            new_locations.extend(feature.location.parts[:-1])
            # add fixed
            new_locations.append(FeatureLocation(feature.location.parts[-1].start, feature.location.parts[-1].end - 3, 1))
        else:
            if len(feature.location.parts) == 1:
                feature.location = FeatureLocation(feature.location.start + 3, feature.location.end, -1)
                continue
            new_locations = []
            # add fixed
            new_locations.append(FeatureLocation(feature.location.parts[0].start + 3, feature.location.parts[0].end, -1))
            # add old
            new_locations.extend(feature.location.parts[1:])
        feature.location = CompoundLocation(new_locations)




def collapse_record_cds(gff_seq_rec):
    """Collapses several CDS features in a SeqRecord into a single CDS with a compound location"""
    new_record = SeqRecord(
        gff_seq_rec.seq,
        gff_seq_rec.id,
        gff_seq_rec.name,
        gff_seq_rec.dbxrefs,
        ""
    )
    new_record.annotations["molecule_type"] = "DNA"

    # gene level
    for cds_id, feature in enumerate(gff_seq_rec.features):
        feature: SeqFeature
        gene_id = gff_seq_rec.id + "+CDS" + str(cds_id + 1)


        # exon level
        exon_locations = get_exons(feature)

        new_feature = SeqFeature(
            location = exon_locations,
            type = "CDS",
            strand = feature.strand,
            id = gene_id,
            qualifiers={
                "ID": gene_id
            }
        )

        new_record.features.append(new_feature)
    return new_record


def sort_locations(record: SeqRecord):
    """Sorts the location parts in each model by the start location of each feature location part"""
    for feature in record.features:
        if feature.strand == -1 and len(feature.location.parts) > 1:
            location_parts = sorted(feature.location.parts, key=lambda parts: parts.start)
            feature.location = CompoundLocation(location_parts)

def convert_gff(gff_file: Path, fasta_file: Path):
    """Converts a GFF file with an associated FASTA file to Genbank format"""
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    gff_iter = GFF.parse(gff_file, fasta_input, {"gff_type": ["CDS"]})

    output_path = gff_file.parent

    if not output_path.exists():
        os.mkdir(output_path)

    # create new gbk file
    gbk_path = gff_file.parent / Path(gff_file.stem + ".gbk")

    # chr level
    with open(gbk_path, mode="w", encoding="utf-8") as out_file_handle:
        for gff_seq_rec in gff_iter:
            new_record = collapse_record_cds(gff_seq_rec)

            add_translations(new_record)

            sort_locations(new_record)

            remove_stop_codons(new_record)

            SeqIO.write(new_record, out_file_handle, "genbank")


if __name__ == "__main__":
    gff_path = Path(sys.argv[1])
    if not gff_path.exists():
        print(f"{gff_path} does not exist")
        exit()

    fasta_path = Path(sys.argv[2])
    if not fasta_path.exists():
        print(f"{fasta_path} does not exist")
        exit()


    convert_gff(gff_path, fasta_path)
