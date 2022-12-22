
"""Contains code to compare two GFF files and determine which models contain differences

Authors:
Arjan Draisma
Jorge Navarro
Jérôme Collemare
"""

from pathlib import Path
import sys

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation, ExactPosition
from Bio.SeqRecord import SeqRecord


def get_exons(feature: SeqFeature) -> CompoundLocation|FeatureLocation:
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

def fix_start(models: list[SeqFeature], offset=-3):
    for model in models:
        if len(model.location.parts) == 1:
            model.location = FeatureLocation(int(model.location.start) + offset, model.location.end, model.location.strand)
            continue
        new_locations = []
        if model.strand == 1:
            fixed_location = FeatureLocation(int(model.location.parts[0].start + offset), model.location.parts[0].end, 1)
            new_locations.append(fixed_location)
            new_locations.extend(model.location.parts[1:])

        else:
            fixed_location = FeatureLocation(int(model.location.parts[-1].start + offset), model.location.parts[-1].end, -1)
            new_locations.extend(model.location.parts[:-1])
            new_locations.append(fixed_location)
        model.location = CompoundLocation(new_locations)


def collapse_record_cds(gff_seq_rec):
    """Collapses several CDS features in a SeqRecord into a single CDS with a compound location"""
    new_record = SeqRecord(gff_seq_rec.seq, gff_seq_rec.id, gff_seq_rec.name, gff_seq_rec.description, gff_seq_rec.dbxrefs)
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
            id = gene_id
        )
        new_feature.qualifiers = {
            "translation": new_feature.translate(gff_seq_rec.seq)
        }

        new_record.features.append(new_feature)
    return new_record

def exons_match(model_a: SeqFeature, model_b: SeqFeature):
    if len(model_a.location.parts) != len(model_b.location.parts):
        return False

    for exon_idx in range(len(model_a.location.parts)):
        a_start = int(model_a.location.parts[exon_idx].start)
        a_end = int(model_a.location.parts[exon_idx].end)
        b_start = int(model_b.location.parts[exon_idx].start)
        b_end = int(model_b.location.parts[exon_idx].end)
        if a_start != b_start or a_end != b_end:
            return False

    return True



def calculate_overlap(a_start, a_end, b_start, b_end):
    if a_start > b_end:
        return 0
    if b_start > a_end:
        return 0

    if a_start >= b_start:
        start = a_start
    else:
        start = b_start

    if b_end >= a_end:
        end = a_end
    else:
        end = b_end

    return end - start


def calculate_exon_len(model: SeqFeature):
    sum = 0
    for exon in model.location.parts:
        start = exon.start
        end = exon.end
        sum += end - start
    return sum


def calculate_exon_overlap(model_a: SeqFeature, model_b: SeqFeature):
    sum = 0
    for exon_a in model_a.location.parts:
        for exon_b in model_b.location.parts:
            a_start = exon_a.start
            a_end = exon_a.end
            b_start = exon_b.start
            b_end = exon_b.end

            sum += calculate_overlap(a_start, a_end, b_start, b_end)
    return sum



def get_match_candidates(models_a: list[SeqFeature], models_b: list[SeqFeature]):
    """Returns two lists:
        - a list of tuples, each of which is a set of models that correspond exactly
        - a list of tuples, each of which is a set of models that correspond to some degree, but not exactly

        the format of the data in the tuples are as follows:
        (model_a_id, model_b_id, sense, len_a, len_b, overlap, percent_overlap_a, percent_overlap_b, num_exons_a, num_exons_b, exon_overlap, exon_len_a, exon_len_b, perc_exon_overlap_a, perc_exon_overlap_b)
    """
    inexact = []
    exact = []

    # very inefficient. but we only need to do it once
    start_at = 0
    first_match = -1
    for idx, model_a in enumerate(models_a):
        print(f"{idx + 1}/{len(models_a)}")
        for b_idx, model_b in enumerate(models_b[start_at:]):
            if model_a.strand != model_b.strand:
                continue


            if model_a.strand == 1:
                a_start = int(model_a.location.parts[0].start)
                a_end = int(model_a.location.parts[-1].end)
                b_start = int(model_b.location.parts[0].start)
                b_end = int(model_b.location.parts[-1].end)
            else:
                a_start = int(model_a.location.parts[-1].start)
                a_end = int(model_a.location.parts[0].end)
                b_start = int(model_b.location.parts[-1].start)
                b_end = int(model_b.location.parts[0].end)

            if a_end < b_start:
                break

            # no overlap at all
            total_overlap = calculate_overlap(a_start, a_end, b_start, b_end)
            if total_overlap == 0:
                continue

            if first_match == -1:
                first_match = b_idx

            if exons_match(model_a, model_b):
                a_len = a_end - a_start
                b_len = b_end - b_start
                sense = "+" if model_a.strand == 1 else "-"
                perc_tot_overlap_a = total_overlap / a_len
                perc_tot_overlap_b = total_overlap / b_len
                exons_a = len(model_a.location.parts)
                exons_b = len(model_b.location.parts)
                exon_len_a = calculate_exon_len(model_a)
                exon_len_b = calculate_exon_len(model_b)
                exon_overlap = calculate_exon_overlap(model_a, model_b)
                perc_exon_overlap_a = exon_overlap / exon_len_a
                perc_exon_overlap_b = exon_overlap / exon_len_b
                exact.append((model_a.id, model_b.id, sense, a_len, b_len, total_overlap, perc_tot_overlap_a, perc_tot_overlap_b, exons_a, exons_b, exon_len_a, exon_len_b, exon_overlap, perc_exon_overlap_a, perc_exon_overlap_b))
                continue


            # some sort of overlap
            a_len = a_end - a_start
            b_len = b_end - b_start
            sense = "+" if model_a.strand == 1 else "-"
            perc_tot_overlap_a = total_overlap / a_len
            perc_tot_overlap_b = total_overlap / b_len
            exons_a = len(model_a.location.parts)
            exons_b = len(model_b.location.parts)
            exon_len_a = calculate_exon_len(model_a)
            exon_len_b = calculate_exon_len(model_b)
            exon_overlap = calculate_exon_overlap(model_a, model_b)
            perc_exon_overlap_a = exon_overlap / exon_len_a
            perc_exon_overlap_b = exon_overlap / exon_len_b
            inexact.append((model_a.id, model_b.id, sense, a_len, b_len, total_overlap, perc_tot_overlap_a, perc_tot_overlap_b, exons_a, exons_b, exon_len_a, exon_len_b, exon_overlap, perc_exon_overlap_a, perc_exon_overlap_b))

        if first_match != -1:
            start_at = first_match
            first_match = -1

    return exact, inexact

def find_diffs(gff_file_a, gff_file_b):
    """Returns a list of tuples, each of which contains a combination of models that correspond, but are different somehow"""

    gff_iter_a = GFF.parse(gff_file_a, None, {"gff_type": ["CDS"]})
    gff_iter_b = GFF.parse(gff_file_b, None, {"gff_type": ["CDS"]})

    a_models = []
    for gff_seq_rec in gff_iter_a:
        new_record: SeqRecord = collapse_record_cds(gff_seq_rec)
        record_features: list[SeqFeature] = new_record.features
        a_models.extend(record_features)

    b_models = []
    for gff_seq_rec in gff_iter_b:
        new_record: SeqRecord = collapse_record_cds(gff_seq_rec)
        record_features: list[SeqFeature] = new_record.features
        b_models.extend(record_features)

    fix_start(b_models)

    exact_matches, inexact_matches = get_match_candidates(a_models, b_models)

    with open("exact_matches.csv", "w", encoding="utf-8") as exact_matches_csv:
        exact_matches_csv.write("model_a_id,model_b_id,sense,len_a,len_b,overlap,percent_overlap_a,percent_overlap_b,num_exons_a,num_exons_b,exon_overlap,exon_len_a,exon_len_b,perc_exon_overlap_a,perc_exon_overlap_b\n")
        for exact_match in exact_matches:
            fields = map(str, list(exact_match))
            exact_matches_csv.write(",".join(fields) + "\n")

    with open("inexact_matches.csv", "w", encoding="utf-8") as inexact_matches_csv:
        inexact_matches_csv.write("model_a_id,model_b_id,sense,len_a,len_b,overlap,percent_overlap_a,percent_overlap_b,num_exons_a,num_exons_b,exon_overlap,exon_len_a,exon_len_b,perc_exon_overlap_a,perc_exon_overlap_b\n")
        for inexact_match in inexact_matches:
            fields = map(str, list(inexact_match))
            inexact_matches_csv.write(",".join(fields) + "\n")







if __name__ == "__main__":
    gff_a_path = Path(sys.argv[1])
    gff_b_path = Path(sys.argv[2])

    find_diffs(gff_a_path, gff_b_path)
