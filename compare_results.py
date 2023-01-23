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

    results_handle = open(result_base_path / Path('results.csv'), mode='w', encoding='utf-8')

    header = "id,strand,start,stop,n_exons,sequence\n"
    print(header, end="")
    results_handle.write(header)


    # stats
    total_curations = 0 # total number of curations which we can analyze
    total_old_equal_truth = 0 # total number of old models that already match the truth
    perfect_curations = 0
    # exons
    total_introns = 0
    curated_matching_introns = 0
    curated_nonmatching_introns = 0
    old_matching_introns = 0
    improvable_introns = 0
    improved_introns = 0
    regressable_introns = 0
    regressed_introns = 0
    unchanged_introns = 0
    # starts
    total_starts = 0
    regressable_starts = 0
    curated_matching_starts = 0
    curated_nonmatching_starts = 0
    unchanged_starts = 0
    regressable_starts =  0
    regressed_starts = 0
    improvable_starts = 0
    improved_starts = 0
    # stops
    total_stops = 0
    regressable_stops = 0
    curated_matching_stops = 0
    curated_nonmatching_stops = 0
    unchanged_stops = 0
    regressable_stops = 0
    regressed_stops = 0
    improvable_stops = 0
    improved_stops = 0

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



        # old gbk
        # same id as the curated model, since the old gbk is where those come from
        old_feature = old_features[curated_model_id]
        old_feature: SeqFeature
        old_aa_seq = old_aa_seqs[curated_model_id]



        # figure out sequence equality
        o_t_equal = old_aa_seq == truth_aa_seq
        if o_t_equal:
            total_old_equal_truth += 1
        o_p_equal = old_aa_seq == curated_aa_seq
        p_t_equal = curated_aa_seq == truth_aa_seq
        if p_t_equal:
            perfect_curations += 1


        total_curations += 1

        # intron sets
        truth_introns = set()
        curated_introns = set()
        old_introns = set()

        if len(truth_feature.location.parts) > 1:
            # get introns from locations
            # go through all but the last feature locations. these correspond to exons.
            # first sort by start
            truth_feature_locations = list(sorted(
                truth_feature.location.parts,
                key = lambda part: part.start
            ))
            for idx, exon in enumerate(truth_feature_locations[:-1]):
                # get intron buy taking end of this exon
                intron_start = exon.end
                # and start of next exon
                intron_end = truth_feature_locations[idx+1].end
                truth_introns.add((intron_start, intron_end))

        # old
        if len(old_feature.location.parts) > 1:
            old_feature_locations = list(sorted(
                old_feature.location.parts,
                key = lambda part: part.start
            ))
            for idx, exon in enumerate(old_feature_locations[:-1]):
                # get intron buy taking end of this exon
                intron_start = exon.end
                # and start of next exon
                intron_end = old_feature_locations[idx+1].end
                old_introns.add((intron_start, intron_end))

        # curated
        if len(curated_feature.location.parts) > 1:
            curated_feature_locations = list(sorted(
                curated_feature.location.parts,
                key = lambda part: part.start
            ))
            for idx, exon in enumerate(curated_feature_locations[:-1]):
                # get intron buy taking end of this exon
                intron_start = exon.end
                # and start of next exon
                intron_end = curated_feature_locations[idx+1].end
                curated_introns.add((intron_start, intron_end))

        # total number of expected exons
        total_introns += len(truth_introns)

        # of which old matched truth
        old_matching_introns += len(truth_introns & old_introns)

        # of which curated match truth
        curated_matching_introns += len(truth_introns & curated_introns)

        curated_nonmatching_introns += len(curated_introns - truth_introns)

        # of which can be improved
        improvable_introns += len(truth_introns - old_introns)

        # of which can regress
        regressable_introns += len(truth_introns & old_introns)

        # of which were incorrect before, but are correct now
        improved_introns += len((truth_introns - old_introns) & curated_introns)

        # of which were correct before, but are not now
        regressed_introns += len((truth_introns & old_introns) - curated_introns)

        # of which are unchanged between curation and old model
        unchanged_introns += len(old_introns & curated_introns)


        # find start and stop

        if truth_feature.strand == 1:
            truth_start = truth_feature.location.start
            curated_start = curated_feature.location.start
            old_start = old_feature.location.start

            truth_stop = truth_feature.location.end
            curated_stop = curated_feature.location.end
            old_stop = old_feature.location.end
        else:
            truth_start = truth_feature.location.end
            curated_start = curated_feature.location.end
            old_start = old_feature.location.end

            truth_stop = truth_feature.location.start
            curated_stop = curated_feature.location.start
            old_stop = old_feature.location.start

        # starts
        total_starts += 1

        if truth_start == curated_start:
            curated_matching_starts += 1
        else:
            curated_nonmatching_starts += 1

        if old_start == curated_start:
            unchanged_starts += 1

        if truth_start == old_start:
            regressable_starts += 1

        if truth_start == old_start and truth_start != curated_start:
            regressed_starts += 1

        if truth_start != old_start:
            improvable_starts += 1

        if truth_start != old_start and truth_start == curated_start:
            improved_starts += 1


        # stops
        total_stops += 1

        if truth_stop == curated_stop:
            curated_matching_stops += 1
        else:
            curated_nonmatching_stops += 1

        if old_stop == curated_stop:
            unchanged_stops += 1

        if truth_stop == old_stop:
            regressable_stops += 1

        if truth_stop == old_stop and truth_stop != curated_stop:
            regressed_stops += 1

        if truth_stop != old_stop:
            improvable_stops += 1

        if truth_stop != old_stop and truth_stop == curated_stop:
            improved_stops += 1


        truth_row = f"{curated_model_id}|TRUE,{truth_feature.strand},{truth_start},{truth_stop},{len(truth_introns)},{truth_aa_seq}\n"
        print(truth_row, end="")
        results_handle.write(truth_row)

        curated_row = f"{curated_model_id}|CRTD,{curated_feature.strand},{curated_start},{curated_stop},{len(curated_introns)},{curated_aa_seq}\n"
        print(curated_row, end="")
        results_handle.write(curated_row)

        old_row = f"{curated_model_id}|ORIG,{old_feature.strand},{old_start},{old_stop},{len(old_introns)},{old_aa_seq}\n"
        print(old_row, end="")
        results_handle.write(old_row)

    with open(result_base_path / Path('summary.txt'), mode='w', encoding='utf-8') as summary_handle:
        summary_handle.write("Summary stats:\n")
        summary_handle.write(f"Total curations: {total_curations}\n")
        summary_handle.write(f"Of which original==truth: {total_old_equal_truth}\n")
        summary_handle.write(f"Perfect curations: {perfect_curations} ({perfect_curations/total_curations})\n")
        summary_handle.write("\n\n")
        summary_handle.write(f"Total introns in truth: {total_introns}\n")
        summary_handle.write(f"Old introns matching truth: {old_matching_introns} ({old_matching_introns/total_introns})\n")
        summary_handle.write(f"Curated introns matching truth: {curated_matching_introns} ({curated_matching_introns/total_introns})\n")
        summary_handle.write(f"Introns that did not match truth, but do match after curation: {improved_introns} ({improved_introns/total_introns})\n")
        summary_handle.write(f"Introns that matched truth, but do not match after curation: {regressed_introns} ({regressed_introns/total_introns})\n")
        summary_handle.write(f"Introns which did not change at all after curation: {unchanged_introns} ({unchanged_introns/total_introns})\n")
        summary_handle.write("\n\n")
        summary_handle.write(f"Total starts in truth: {total_starts}\n")
        summary_handle.write(f"Old starts matching truth: {regressable_starts} ({regressable_starts/total_starts})\n")
        summary_handle.write(f"Curated starts matching truth: {curated_matching_starts} ({curated_matching_starts/total_starts})\n")
        summary_handle.write(f"Starts that did not match truth, but do match after curation: {improved_starts} ({improved_starts/total_starts})\n")
        summary_handle.write(f"Starts that matched truth, but do not match after curation: {regressed_starts} ({regressed_starts/total_starts})\n")
        summary_handle.write(f"Starts which did not change at all after curation: {unchanged_starts} ({unchanged_starts/total_starts})\n")
        summary_handle.write("\n\n")
        summary_handle.write(f"Total stops in truth: {total_stops}\n")
        summary_handle.write(f"Old stops matching truth: {regressable_stops} ({regressable_stops/total_stops})\n")
        summary_handle.write(f"Curated stops matching truth: {curated_matching_stops} ({curated_matching_stops/total_stops})\n")
        summary_handle.write(f"Stops that did not match truth, but do match after curation: {improved_stops} ({improved_stops/total_stops})\n")
        summary_handle.write(f"Stops that matched truth, but do not match after curation: {regressed_stops} ({regressed_stops/total_stops})\n")
        summary_handle.write(f"Stops which did not change at all after curation: {unchanged_stops} ({unchanged_stops/total_stops})\n")
        summary_handle.write("\n")
        summary_handle.write("")
        # line for excel summaries
        excel_fields = map(str, [
            total_curations,
            perfect_curations,
            total_introns,
            len(old_introns),
            improvable_introns,
            regressable_introns,
            len(curated_introns),
            improved_introns,
            regressed_introns,
            unchanged_introns,
            improved_introns / improvable_introns,
            regressed_introns / regressable_introns,
            total_starts,
            improvable_starts,
            regressable_starts,
            curated_matching_starts,
            curated_nonmatching_starts,
            improved_starts,
            regressed_starts,
            unchanged_starts,
            improved_starts / improvable_starts,
            regressed_starts / regressable_starts,
            total_stops,
            improvable_stops,
            regressable_stops,
            curated_matching_stops,
            curated_nonmatching_stops,
            improved_stops,
            regressed_stops,
            unchanged_stops,
            improved_stops / improvable_stops,
            regressed_stops / regressable_stops
        ])
        summary_handle.write(",".join(excel_fields))
        summary_handle.write("\n")
    results_handle.close()
