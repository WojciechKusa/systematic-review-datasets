import argparse
import os
from subprocess import Popen, PIPE

import pandas
from Bio import Entrez
from pandas import *
from tqdm import tqdm


def read_date(date_file):
    date_dict = {}
    with open(date_file) as f:
        for line in f:
            topic, start_day, end_day = line.split()
            new_start_day = start_day[0:4] + "/" + start_day[4:6] + "/" + start_day[6:8]
            new_end_day = end_day[0:4] + "/" + end_day[4:6] + "/" + end_day[6:8]
            date_dict[topic] = (new_start_day, new_end_day)
    return date_dict


def write_ranking(topic, final_result, output):
    for rank, r in enumerate(final_result):
        output.write(f"{topic} 0 {r} {rank + 1} {1/(rank+1)} rank\n")


def temporal_submission(date, query, email):
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=10000,
        email=email,
        mindate=date[0].replace("-", "/"),
        maxdate=date[1].replace("-", "/"),
    )
    record = Entrez.read(handle)
    count = int(record["Count"])
    id_list = record["IdList"]
    return count, id_list


def divide_dates(mindate, maxdate, query, email):
    original_chunks = [[mindate, maxdate]]
    id_lists = []
    # final_chunks = []
    while len(original_chunks) > 0:
        current_date_range_count, current_id_list = temporal_submission(
            original_chunks[0], query, email
        )
        if current_date_range_count > 10000:
            times = 3
            # math.ceil(current_date_range_count/10000)
            ts1 = pandas.Timestamp(original_chunks[0][0])
            ts2 = pandas.Timestamp(original_chunks[0][1])
            ading_date = (ts2 - ts1) / times
            previous_mean = original_chunks[0][0]
            for index in range(1, times):
                chunk_mean = ts1 + (index * ading_date)
                chunk_mean = str(chunk_mean)[:10]
                chunk_now = [previous_mean, str(chunk_mean)]
                previous_mean = str(chunk_mean)
                original_chunks.append(chunk_now)
            chunk_last = [previous_mean, original_chunks[0][1]]
            original_chunks.append(chunk_last)
            original_chunks.pop(0)
        else:
            id_lists.extend(current_id_list)
            original_chunks.pop(0)
            print(original_chunks)
    return set(id_lists)


def submit_result(query, email, date_info):
    min_date = date_info[0].replace("/", "-")
    max_date = date_info[1].replace("/", "-")
    results = divide_dates(min_date, max_date, query, email)
    # results = []
    # for date_chunk in date_chunks:
    #     handle = Entrez.esearch(db="pubmed", term=query, retmax=10000, email=email, mindate=date_chunk[0].replace('-', '/'),
    #                             maxdate=date_chunk[1].replace('-', '/'))
    #     record = Entrez.read(handle)
    #     result = record["IdList"]
    #     results += result
    return set(results)


# def submit_result(query, email, date_info):
#     min_date = date_info[0]
#     max_date = date_info[1]
#     query = query.replace('"', '')
#     handle = Entrez.esearch(db="pubmed", term=query, retmax=10000, email=email, mindate=min_date,
#                             maxdate=max_date)
#
#     record = Entrez.read(handle)
#     results = []
#     count = int(record["Count"])
#     if count<=10000:
#         result = record["IdList"]
#         results += result
#     else:
#         query_string = f'esearch -db pubmed -query "{query}" | efilter -mindate {min_date} -maxdate {max_date} | efetch -format uid'
#         print(query_string)
#         p = subprocess.Popen(query_string ,stdout=subprocess.PIPE, shell=True)
#         (output, err) = p.communicate()
#         result = output.decode("utf-8").strip().split('\n')
#         results += result
#     return results


def read_run(run_file):
    run_results = {}
    with open(run_file) as f:
        for l in f:
            try:
                qid, _, docid, rank, score, _ = l.strip().split()
            except ValueError:
                raise ValueError("Wrong run format.")
            if qid not in run_results:
                run_results[qid] = {}
            run_results[qid][docid] = float(score)
    return run_results


def read_qrel(qrel_file):
    qrels = {}
    with open(qrel_file, "r") as f:
        for l in f:
            try:
                qid, _, docid, rel = l.strip().split()
            except ValueError:
                raise ValueError("Wrong qrel format.")
            if (qid not in qrels) and int(rel) != 0:
                qrels[qid] = {}
            if int(rel) != 0:
                qrels[qid][docid] = int(rel)
    return qrels


def eval(trec_eval, qrel_path, output_path):
    result_dict = {}
    command = (
        trec_eval + " -m set_recall -m set_P -m set_F " + qrel_path + " " + output_path
    )
    print(command)
    results = Popen(
        command, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True
    ).stdout.readlines()
    for result in results:
        items = result.split()
        if (len(items) == 3) and (items[1] == "all"):

            result_dict[items[0]] = float(items[-1])
    return result_dict


def main(args):
    date_dict = read_date(args.date_file)

    original_file = os.path.join(args.out_folder, "original.trec")
    original_set = set()

    f_original = open(original_file, "r+")
    for line in f_original.readlines():
        qid = line.split(" ", 1)[0]
        if qid not in f_original:
            original_set.add(qid)

    print(len(original_set))

    # q1_file = os.path.join(args.out_folder, 'q1_result.trec')
    # q2_file = os.path.join(args.out_folder, 'q2_result.trec')
    # q1_set = set()
    # q2_set = set()
    # f1 = open(q1_file, 'r+')
    # f2 = open(q2_file, 'r+')
    # #
    # for line in f1.readlines():
    #     qid = line.split(' ', 1)[0]
    #     if qid not in q1_set:
    #         q1_set.add(qid)
    #
    # for line in f2.readlines():
    #     qid = line.split(' ', 1)[0]
    #     if qid not in q2_set:
    #         q2_set.add(qid)

    # print(len(q1_set), len(q2_set))
    csv_data = read_csv(args.input_file)
    ids = csv_data["topicid"].tolist()
    ql_dict = read_qrel(args.input_qrel)

    originals = csv_data["original_query"].to_list()
    for id, original in tqdm(zip(ids, originals)):
        if original is None:
            continue
        if id not in ql_dict:
            continue
        if id in original_set:
            continue
        print(id)
        date = ("1946/01/01", "2018/12/31")
        if id in date_dict:
            date = date_dict[id]
        if (id == "CD010772") or (id == "CD010771") or (id == "CD008054"):
            continue
        original_result = submit_result(original, args.email, date)
        write_ranking(id, original_result, f_original)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True)
    parser.add_argument("--input_qrel", type=str, required=True)
    parser.add_argument("--date_file", type=str, default="data/combined_pubdates")
    parser.add_argument("--email", type=str, default="tester@gmail.com")
    parser.add_argument("--out_folder", type=str, default="output/")
    args = parser.parse_args()
    main(args)
