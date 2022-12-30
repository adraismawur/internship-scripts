
"""Contains code to compare two GFF files and determine which models contain differences

Authors:
Arjan Draisma
Jorge Navarro
Jérôme Collemare
"""

from math import ceil
from pathlib import Path
from pyexpat import model
import sys

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation, ExactPosition
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

def sort_locations(models: list[SeqFeature]):
    """Sorts the location parts in each model by the start location of the model part"""
    for model in models:
        if model.strand == -1 and len(model.location.parts) > 1:
            location_parts = sorted(model.location.parts, key=lambda parts: parts.start)
            model.location = CompoundLocation(location_parts)

def fix_ends(models: list[SeqFeature]):
    for model in models:
        if len(model.location.parts) == 1:
            if model.location.end - model.location.start < 3:
                continue
            if model.strand == 1:
                model.location = FeatureLocation(int(model.location.start), model.location.end - 3, 1)
            else:
                model.location = FeatureLocation(int(model.location.start) + 3, model.location.end, -1)
            continue
        new_locations = []
        if model.strand == 1:
            if model.location.parts[-1].end - model.location.parts[-1].start < 3:
                continue
            new_locations.extend(model.location.parts[:-1])
            fixed_location = FeatureLocation(int(model.location.parts[-1].start), model.location.parts[-1].end - 3, 1)
            new_locations.append(fixed_location)

        else:
            if model.location.parts[0].end - model.location.parts[0].start < 3:
                continue
            fixed_location = FeatureLocation(int(model.location.parts[0].start + 3), model.location.parts[0].end, -1)
            new_locations.append(fixed_location)
            new_locations.extend(model.location.parts[1:])
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

    no_matches_a = set([model.id for model in models_a])
    no_matches_b = set([model.id for model in models_b])

    # very inefficient. but we only need to do it once
    start_at = 0
    first_match = -1
    for idx, model_a in enumerate(models_a):
        print(f"{idx + 1}/{len(models_a)}")
        for b_idx, model_b in enumerate(models_b[start_at:]):
            if model_a.strand != model_b.strand:
                continue

            if model_a.id[0:7] != model_b.id[0:7]:
                continue


            # if model_a.strand == 1:
            a_start = int(model_a.location.parts[0].start)
            a_end = int(model_a.location.parts[-1].end)
            b_start = int(model_b.location.parts[0].start)
            b_end = int(model_b.location.parts[-1].end)
            # else:
            #     a_start = int(model_a.location.parts[-1].start)
            #     a_end = int(model_a.location.parts[0].end)
            #     b_start = int(model_b.location.parts[-1].start)
            #     b_end = int(model_b.location.parts[0].end)

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

                if model_a.id in no_matches_a:
                    no_matches_a.remove(model_a.id)

                if model_b.id in no_matches_b:
                    no_matches_b.remove(model_b.id)

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

            if model_a.id in no_matches_a:
                no_matches_a.remove(model_a.id)

            if model_b.id in no_matches_b:
                no_matches_b.remove(model_b.id)

        if first_match != -1:
            start_at = first_match
            first_match = -1

    no_matches_a = sorted(list(no_matches_a))
    no_matches_b = sorted(list(no_matches_b))

    return exact, inexact, no_matches_a, no_matches_b

def write_matches(matches, output_path):
    with open(output_path, "w", encoding="utf-8") as matches_csv:
        matches_csv.write("model_a_id,model_b_id,sense,len_a,len_b,overlap,percent_overlap_a,percent_overlap_b,num_exons_a,num_exons_b,exon_overlap,exon_len_a,exon_len_b,perc_exon_overlap_a,perc_exon_overlap_b\n")
        for match in matches:
            fields = map(str, list(match))
            matches_csv.write(",".join(fields) + "\n")

def read_matches(input_path: Path):
    if not input_path.exists:
        return None

    matches = []

    with open(input_path):
        for line in input_path:
            line: str
            matches.append(tuple(line.rstrip().split(",")))
    return matches


def write_list(list_obj, output_path):
    with open(output_path, "w", encoding="utf-8") as list_txt:
        for item in list_obj:
            list_txt.write(f"{item}\n")



def find_diffs(gff_file_a, gff_file_b):
    """Returns a list of tuples, each of which contains a combination of models that correspond, but are different somehow"""

    exact_match_path = Path("exact_matches.csv")
    inexact_match_path = Path("inexact_matches.csv")
    no_matches_a_path = Path("no_matches_a.txt")
    no_matches_b_path = Path("no_matches_b.txt")


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

    print(f"{len(a_models)} models in {gff_file_a.stem}")
    print(f"{len(b_models)} models in {gff_file_b.stem}")

    # sort the model locations. they start out in reverse order which is annoying to work with
    sort_locations(a_models)
    sort_locations(b_models)

    # this script was made to compare two models, one of which included stop codons in its cds regions
    # this function removes those codons
    fix_ends(a_models)

    exact_matches, inexact_matches, no_matches_a, no_matches_b = get_match_candidates(a_models, b_models)

    write_matches(exact_matches, exact_match_path)
    write_matches(inexact_matches, inexact_match_path)

    write_list(no_matches_a, no_matches_a_path)
    write_list(no_matches_b, no_matches_b_path)


HELP = (
    f"usage: {__file__} annotations_a.gff annotations_b.gff"
)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Not enough arguments.")
        print(HELP)
        exit()

    gff_a_path = Path(sys.argv[1])
    if not gff_a_path.exists():
        print(f"{gff_a_path} does not exist")
        exit()

    gff_b_path = Path(sys.argv[2])
    if not gff_b_path.exists():
        print(f"{gff_b_path} does not exist")
        exit()

    find_diffs(gff_a_path, gff_b_path)

    # generate an include list. should be deterministic
    with open("inexact_matches.csv", encoding='utf-8') as inexact_matches_handle:
        # skip header
        inexact_matches_handle.readline()
        inexact_matches = inexact_matches_handle.readlines()

    inexact_matches = [row.split(",") for row in inexact_matches]

    # filter out any % identity < 50%
    inexact_matches = filter(lambda row: float(row[6]) > 0.5 and float(row[7]) > 0.5, inexact_matches)

    # filter out any % exon identity < 50%
    inexact_matches = filter(lambda row: float(row[13]) > 0.5 and float(row[14]) > 0.5, inexact_matches)

    # to list
    inexact_matches = list(inexact_matches)

    include_ids_lines = []
    # we want 100. start at 0, divide the length by 100, and increase by that number until we hit len
    for i in range(0, len(inexact_matches), ceil(len(inexact_matches) / 100)):
        include_ids_lines.append(inexact_matches[i][1] + '\n')

    with open('include_ids.txt', 'w', encoding='utf-8') as include_ids_handle:
        include_ids_handle.writelines(include_ids_lines)
