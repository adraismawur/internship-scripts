"""Functions to compare results from gemcapy with assumed ground truth model(s)"""



from pathlib import Path
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, CompoundLocation
from Bio.SeqRecord import SeqRecord


def get_inexact_matches_map(inexact_matches_path: Path):
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

def fix_feature_locations(feature):
    first_start = feature.location.parts[0].start
    last_start = feature.location.parts[-1].start
    is_ascending = last_start > first_start
    if len(feature.location.parts) > 1:
        sense_wrong_order = feature.location.strand == 1 and not is_ascending
        antisense_wrong_order = feature.location.strand == -1 and is_ascending
        if antisense_wrong_order or sense_wrong_order:
            reversed_order = list(reversed(feature.location.parts))
            feature.location = CompoundLocation(reversed_order)


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


def unpack_gbk_to_map(feature_ids: set, gbk_path: Path):
    """returns two maps, one with feature ids as keys and SeqFeatures as items, and another with feature ids as keys and translated aa Seq objects as items"""
    feature_map = {}
    seq_map = {}
    with open(gbk_path, encoding='utf-8') as gbk_handle:
        gbk_seq_recs = SeqIO.parse(gbk_handle, 'genbank')
        for gbk_seq_rec in gbk_seq_recs:
            gbk_seq_rec: SeqRecord
            gbk_seq: Seq = gbk_seq_rec.seq

            for feature in gbk_seq_rec.features:
                feature: SeqFeature
                gbk_feature_id = feature.qualifiers["ID"][0]
                if gbk_feature_id in feature_ids:
                    fix_feature_locations(feature)
                    feature_map[gbk_feature_id] = feature
                    seq_map[gbk_feature_id] = feature.translate(gbk_seq, cds=False)

    return feature_map, seq_map


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print(f"usage: {__file__} <internship-scripts output dir> <inexact_matches.csv> <ground truth gbk> <augustus result gbk> <result output dir>")
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

    result_base_path = Path(sys.argv[5])
    result_base_path.mkdir(parents=True, exist_ok=True)


    output_folders = sorted(list(output_dir.glob("*")))

    # map of curated model ID to truth model id
    truth_model_map = get_inexact_matches_map(inexact_matches_path)

    # get the sets of ids from items and keys
    truth_id_set = set(truth_model_map.values())
    old_id_set = set(truth_model_map.keys())

    truth_features, truth_aa_seqs = unpack_gbk_to_map(truth_id_set, truth_gbk_path)
    old_features, old_aa_seqs = unpack_gbk_to_map(old_id_set, old_gbk_path)

    summary_lines = []
    results_handle = open(result_base_path / Path('results.csv'), mode='w', encoding='utf-8')

    header = "id,total_len,exon_len,start,stop,n_exons,sequence\n"
    print(header, end="")
    results_handle.write(header)


    # stats
    total_curations = 0 # total number of curations which we can analyze
    total_old_equal_truth = 0 # total number of old models that already match the truth
    perfect_curations = 0
    total_exons = 0
    curated_matching_exons = 0
    old_matching_exons = 0
    fixed_exons = 0
    broken_exons = 0
    unchanged_exons = 0

    for output_folder in output_folders:
        curated_model_id = output_folder.name
        truth_model_id = truth_model_map[curated_model_id]

        if not run_succeeded(output_folder):
            print(f"{curated_model_id} ended with an error. skipping...")
            continue


        # "truth" gbk
        truth_feature = truth_features[truth_model_id]
        truth_aa_seq = truth_aa_seqs[truth_model_id]
        truth_feature: SeqFeature

        truth_row = f"{curated_model_id}|TRUE,{0},{0},{0},{0},{truth_aa_seq}\n"
        print(truth_row, end="")
        results_handle.write(truth_row)

        # curated gbk
        curated_gbk_path = output_folder / Path(f"output/{curated_model_id}.gbk")
        with open(curated_gbk_path, encoding='utf-8') as output_gbk_handle:
            predicted_seq_recs = SeqIO.parse(output_gbk_handle, 'genbank')
            # we know there will only be one
            predicted_seq_rec: SeqRecord = list(predicted_seq_recs)[0]

        # retrieve nucleotide seq and features
        pad = (3 - len(predicted_seq_rec.seq) % 3) * "N"
        curated_dna_seq = predicted_seq_rec.seq + pad
        curated_feature: SeqFeature = predicted_seq_rec.features[0]

        # ensure location is in the right order
        first_start = curated_feature.location.parts[0].start
        last_start = curated_feature.location.parts[-1].start
        is_ascending = last_start > first_start
        if len(curated_feature.location.parts) > 1:
            sense_wrong_order = curated_feature.location.strand == 1 and not is_ascending
            antisense_wrong_order = curated_feature.location.strand == -1 and is_ascending
            if antisense_wrong_order or sense_wrong_order:
                reversed_order = list(reversed(curated_feature.location.parts))
                curated_feature.location = CompoundLocation(reversed_order)

        # translate the feature
        curated_aa_seq = curated_feature.translate(curated_dna_seq, cds=False)

        curated_row = f"{curated_model_id}|CRTD,{0},{0},{0},{0},{curated_aa_seq}\n"
        print(curated_row, end="")
        results_handle.write(curated_row)


        # old gbk
        # same id as the curated model, since the old gbk is where those come from
        old_feature = old_features[curated_model_id]
        old_feature: SeqFeature
        old_aa_seq = old_aa_seqs[curated_model_id]

        old_row = f"{curated_model_id}|ORIG,{0},{0},{0},{0},{old_aa_seq}\n"
        print(old_row, end="")
        results_handle.write(old_row)


        # summaries
        summary_lines.append(f"{curated_model_id} (original {truth_model_id}):\n")

        summary_lines.append("Translated sequence lengths:\n")
        summary_lines.append(f"GROUND TRUTH: {len(truth_aa_seq)}\n")
        summary_lines.append(f"CURATION:     {len(curated_aa_seq)}\n")
        summary_lines.append(f"ORIGINAL:     {len(old_aa_seq)}\n")

        # figure out sequence equality
        o_t_equal = old_aa_seq == truth_aa_seq
        if o_t_equal:
            total_old_equal_truth += 1
        o_p_equal = old_aa_seq == curated_aa_seq
        p_t_equal = curated_aa_seq == truth_aa_seq
        if p_t_equal:
            perfect_curations += 1

        summary_lines.append("Equality:\n")
        summary_lines.append(f"ORIGINAL TO TRUTH (expect False): {o_t_equal}\n")
        summary_lines.append(f"ORIGINAL TO PRED (expect False): {o_p_equal}\n")
        summary_lines.append(f"PREDICTED TO TRUTH (expect True): {p_t_equal}\n")

        total_curations += 1

        # exons
        truth_exons = set([(int(p.start), int(p.end)) for p in truth_feature.location.parts])
        curated_exons = set([(int(p.start), int(p.end)) for p in curated_feature.location.parts])
        old_exons = set([(int(p.start), int(p.end)) for p in old_feature.location.parts])

        # total number of expected exons
        total_exons += len(truth_exons)

        # of which old matched truth
        old_matching_exons += len(truth_exons & old_exons)

        # of which curated match truth
        curated_matching_exons += len(truth_exons & curated_exons)

        # of which were incorrect before, but are correct now
        fixed_exons += len((truth_exons - old_exons) & curated_exons)

        # of which were correct before, but are not now
        broken_exons += len((truth_exons & old_exons) - curated_exons)

        # of which are unchanged between truth, curation and old model
        unchanged_exons += len(truth_exons & old_exons & curated_exons)


        summary_lines.append("\n")

    with open(result_base_path / Path('summary.txt'), mode='w', encoding='utf-8') as summary_handle:
        summary_handle.write("Summary stats:\n")
        summary_handle.write(f"Total curations: {total_curations}\n")
        summary_handle.write(f"Of which original==truth: {total_old_equal_truth}\n")
        summary_handle.write(f"Perfect curations: {perfect_curations} ({perfect_curations/total_curations})\n")
        summary_handle.write(f"Total exons in truth: {total_exons}\n")
        summary_handle.write(f"Old exons matching truth: {old_matching_exons} ({old_matching_exons/total_exons})\n")
        summary_handle.write(f"Curated exons matching truth: {curated_matching_exons} ({curated_matching_exons/total_exons})\n")
        summary_handle.write(f"Exons that did not match truth, but do match after curation: {fixed_exons} ({fixed_exons/total_exons})\n")
        summary_handle.write(f"Exons that matched truth, but do not match after curation: {broken_exons} ({broken_exons/total_exons})\n")
        summary_handle.write(f"Exons which did not change at all after curation: {unchanged_exons} ({unchanged_exons/total_exons})\n")
        summary_handle.write("\n\n")
    results_handle.close()
