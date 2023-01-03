"""Functions to compare results from gemcapy with assumed ground truth model(s)"""



from pathlib import Path
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, CompoundLocation
from Bio.SeqRecord import SeqRecord


def get_truth_model_map(inexact_matches_path: Path):
    with open(inexact_matches_path, encoding='utf-8') as inexact_matches_handle:
        lines = inexact_matches_handle.readlines()

    matches = [line.split(",")[0:2] for line in lines]

    return {match[1]: match[0] for match in matches}


def run_succeeded(output_folder: Path):
    last_log = list(sorted(output_folder.glob("*.log")))[-1]

    with open(last_log, encoding='utf-8') as log_handle:
        for line in log_handle:
            if "Differences in model " + output_folder.name in line:
                return True
    return False

def get_model(feature_id: str, gbk_path: Path):
    """Returns the first feature in a gbk file that matches the feature_id"""
    with open(gbk_path, encoding='utf-8') as gbk_handle:
        gbk_seq_recs = SeqIO.parse(gbk_handle, 'genbank')
        for gbk_seq_rec in gbk_seq_recs:
            gbk_seq_rec: SeqRecord
            gbk_seq: Seq = gbk_seq_rec.seq

            for feature in gbk_seq_rec.features:
                feature: SeqFeature
                gbk_feature_id = feature.qualifiers["ID"][0]
                if gbk_feature_id == feature_id:
                    return feature, gbk_seq


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(f"usage: {__file__} <internship-scripts output dir> <inexact_matches.csv> <ground truth gbk>")
        exit()

    output_dir = Path(sys.argv[1])
    if not output_dir.exists():
        print(f"{output_dir} does not exist")
        exit()

    inexact_matches_path = Path(sys.argv[2])
    if not inexact_matches_path.exists():
        print(f"{inexact_matches_path} does not exist")
        exit()

    truth_gbk_path = Path(sys.argv[3])
    if not truth_gbk_path.exists():
        print(f"{truth_gbk_path} does not exist")
        exit()

    old_gbk_path = Path(sys.argv[4])
    if not old_gbk_path.exists():
        print(f"{old_gbk_path} does not exist")
        exit()


    output_folders = sorted(list(output_dir.glob("*")))
    truth_model_map = get_truth_model_map(inexact_matches_path)


    print("id,total_len,exon_len,start,stop,n_exons,sequence\n")


    for output_folder in output_folders:
        predicted_model_name = output_folder.name
        truth_model_name = truth_model_map[predicted_model_name]

        if not run_succeeded(output_folder):
            print(f"{predicted_model_name} ended with an error. skipping...")
            continue


        # "truth" gbk
        truth_feature, truth_nucleotide_seq = get_model(truth_model_name, truth_gbk_path)
        truth_feature: SeqFeature
        truth_nucleotide_seq: Seq
        pad = (3 - len(truth_nucleotide_seq) % 3) * "N"
        truth_nucleotide_seq = truth_nucleotide_seq + pad

        truth_annotated_aa_sequence = truth_feature.qualifiers["translation"][0]
        truth_translated_aa_sequence = truth_feature.translate(truth_nucleotide_seq, cds=False)

        print(f"{predicted_model_name}|TRUE,{0},{0},{0},{0},{0},{0},{0}")

        # prediction gbk
        output_gbk_path = output_folder / Path(f"output/{predicted_model_name}.gbk")
        with open(output_gbk_path, encoding='utf-8') as output_gbk_handle:
            predicted_seq_recs = SeqIO.parse(output_gbk_handle, 'genbank')
            # we know there will only be one
            predicted_seq_rec: SeqRecord = list(predicted_seq_recs)[0]

        pad = (3 - len(predicted_seq_rec.seq) % 3) * "N"
        predicted_nucleotide_seq = predicted_seq_rec.seq + pad
        predicted_feature: SeqFeature = predicted_seq_rec.features[0]

        if len(predicted_feature.location.parts) > 1 and predicted_feature.strand == -1:
            predicted_feature.location = CompoundLocation(list(reversed(predicted_feature.location.parts)))

        predicted_annotated_aa_sequence = predicted_feature.qualifiers["translation"][0]
        predicted_translated_aa_sequence = predicted_feature.translate(predicted_nucleotide_seq, cds=False)

        print(f"{predicted_model_name}|PRED,{0},{0},{0},{0},{0},{0},{0}")


        # old gbk
        old_feature, old_nucleotide_seq = get_model(predicted_model_name, old_gbk_path)
        old_feature: SeqFeature
        old_nucleotide_seq: Seq
        pad = (3 - len(old_nucleotide_seq) % 3) * "N"
        old_nucleotide_seq = old_nucleotide_seq + pad

        old_annotated_aa_sequence = old_feature.qualifiers["translation"][0]
        old_translated_aa_sequence = old_feature.translate(old_nucleotide_seq, cds=False)


        print(f"{predicted_model_name}|ORIG,{0},{0},{0},{0},{0},{0},{0}")
