#! /usr/bin/env python

import math
from scipy.stats import nbinom
import polars as pl
import numpy as np


def calc_loglikelihood(PR1, PR2, D_count, marker_list, max_copy_number=2):
    # mprob_k = sum of prob_i * prob_j where i + j = k
    mprob_exprs = []
    for k in range(2 * max_copy_number + 1):
        terms = [pl.col(f"prob_{i}") * pl.col(f"prob_{k-i}_PR2")
                 for i in range(max_copy_number + 1) if 0 <= k - i <= max_copy_number]
        mprob_exprs.append(sum(terms).alias(f"mprob_{k}"))

    PR12 = PR1.filter(pl.col("Marker").is_in(marker_list)) \
        .join(PR2, on="Marker", suffix="_PR2") \
        .with_columns(mprob_exprs)

    tprob_expr = sum([pl.col(f"mprob_{k}") * pl.col(f"dprob_{k}")
                      for k in range(2 * max_copy_number + 1)]).alias("tprob")

    PR12_count = PR12.join(D_count, on="Marker").with_columns([tprob_expr])

    tL = PR12_count["tprob"].log().sum()

    return(tL)


def calc_loglikelihood_single(PR1, D_count, marker_list):

    PR1 = PR1.filter(pl.col("Marker").is_in(marker_list)) \
        .with_columns([
            (pl.col("prob_0")).alias("mprob_0"),
            (pl.col("prob_1")).alias("mprob_1"),
            (pl.col("prob_2")).alias("mprob_2")
        ])

    PR1_count = (
        PR1.join(D_count, on="Marker") \
        .with_columns([
            (pl.col("mprob_0") * pl.col("dprob_0") +
             pl.col("mprob_1") * pl.col("dprob_1") +
             pl.col("mprob_2") * pl.col("dprob_2")).alias("tprob_012")
        ])
    )

    tL = PR1_count["tprob_012"].log().sum()

    return(tL)


def calc_posterior_prob(PR1, PR2, D_count, max_copy_number=2):
    # Generate (N+1)² mprob states: mprob_ij = prob_i * prob_j
    mprob_exprs = [(pl.col(f"prob_{i}") * pl.col(f"prob_{j}_PR2")).alias(f"mprob_{i}{j}")
                   for i in range(max_copy_number + 1) for j in range(max_copy_number + 1)]

    PR12 = PR1.join(PR2, on="Marker", suffix="_PR2").with_columns(mprob_exprs)

    # tprob_ij = mprob_ij * dprob_{i+j}
    tprob_exprs = [(pl.col(f"mprob_{i}{j}") * pl.col(f"dprob_{i+j}")).alias(f"tprob_{i}{j}")
                   for i in range(max_copy_number + 1) for j in range(max_copy_number + 1)]

    # Sum all tprob values
    tprob_sum_expr = sum([pl.col(f"tprob_{i}{j}") for i in range(max_copy_number + 1)
                          for j in range(max_copy_number + 1)]).alias("tprob_sum")

    # Normalize to get posterior probabilities
    prob_exprs = [(pl.col(f"tprob_{i}{j}") / pl.col("tprob_sum")).alias(f"Prob_{i}{j}")
                  for i in range(max_copy_number + 1) for j in range(max_copy_number + 1)]

    PR12_count = PR12.join(D_count, on="Marker") \
        .with_columns(tprob_exprs) \
        .with_columns([tprob_sum_expr]) \
        .with_columns(prob_exprs)

    return(PR12_count)


def calc_posterior_prob_single(PR1, D_count):

    PR1 = PR1 \
        .with_columns([
            (pl.col("prob_0")).alias("mprob_0"),
            (pl.col("prob_1")).alias("mprob_1"),
            (pl.col("prob_2")).alias("mprob_2")
        ])

    PR1_count = PR1.join(D_count, on="Marker") \
        .with_columns([
            (pl.col("mprob_0") * pl.col("dprob_0")).alias("tprob_0"),
            (pl.col("mprob_1") * pl.col("dprob_1")).alias("tprob_1"), 
            (pl.col("mprob_2") * pl.col("dprob_2")).alias("tprob_2")
        ]) \
        .with_columns([
            (pl.col("tprob_0") + pl.col("tprob_1") + pl.col("tprob_2")).alias("tprob_sum")]) \
        .with_columns([
            (pl.col("tprob_0") / pl.col("tprob_sum")).alias("Prob_0"),
            (pl.col("tprob_1") / pl.col("tprob_sum")).alias("Prob_1"),
            (pl.col("tprob_2") / pl.col("tprob_sum")).alias("Prob_2")
        ]) 

    return(PR1_count)


def generate_samples(probability_matrix, num_samples):
    """
    Generate samples based on a given probability matrix.

    Parameters:
        probability_matrix (numpy.ndarray):
            A 2D array where each row contains probabilities for categories.
        num_samples (int):
            Number of samples to generate per row.

    Returns:
        numpy.ndarray:
            A 2D array where each column contains the generated samples for one row.
    """
    # Define category labels based on the number of columns in the probability matrix
    categories = np.arange(probability_matrix.shape[1])

    # Preallocate memory for the samples
    num_rows = probability_matrix.shape[0]
    samples = np.empty((num_samples, num_rows), dtype=int)

    # Generate samples for each row in the probability matrix
    for i, probs in enumerate(probability_matrix):
        # Compute cumulative probabilities
        cumulative_probs = np.cumsum(probs)

        # Generate uniform random numbers
        random_uniform = np.random.uniform(0, 1, num_samples)

        # Use vectorized operations to assign categories based on random uniform values
        indices = np.searchsorted(cumulative_probs, random_uniform, side="right")
        samples[:, i] = categories[indices]

    return samples  # No need to transpose as preallocation ensures correct shape



def estimage_cosine_dist(probability_matrix, hap1_vec, hap2_vec, num_samples = 1000):

    target_columns = ["Prob_00", "Prob_01", "Prob_10", "Prob_11", "Prob_02", "Prob_20", "Prob_12", "Prob_21", "Prob_22"]
    column_indeces = [probability_matrix.columns.index(col) for col in target_columns]

    probability_matrix_np = probability_matrix.to_numpy()
    # Ensure the probabilities are normalized
    probability_matrix_np = probability_matrix_np / probability_matrix_np.sum(axis=1, keepdims=True)

    # Generate samples
    pD_samples = generate_samples(probability_matrix_np, num_samples)

    hap1_samples = (pD_samples == column_indeces[2]).astype(int) + \
               (pD_samples == column_indeces[3]).astype(int) + \
               (pD_samples == column_indeces[6]).astype(int) + \
               2 * (pD_samples == column_indeces[5]).astype(int) + \
               2 * (pD_samples == column_indeces[7]).astype(int) + \
               2 * (pD_samples == column_indeces[8]).astype(int)

    hap2_samples = (pD_samples == column_indeces[1]).astype(int) + \
               (pD_samples == column_indeces[3]).astype(int) + \
               (pD_samples == column_indeces[7]).astype(int) + \
               2 * (pD_samples == column_indeces[4]).astype(int) + \
               2 * (pD_samples == column_indeces[6]).astype(int) + \
               2 * (pD_samples == column_indeces[8]).astype(int)


    # Compute cosine distances
    cos_dist1 = np.mean(1 - (hap1_samples @ hap1_vec) / (np.sqrt(np.sum(hap1_samples**2, axis=1)) * np.sqrt(np.sum(hap1_vec**2))))
    cos_dist2 = np.mean(1 - (hap2_samples @ hap2_vec) / (np.sqrt(np.sum(hap2_samples**2, axis=1)) * np.sqrt(np.sum(hap2_vec**2))))

    # Print cosine distances
    return(round(cos_dist1, 8), round(cos_dist2, 8))


def estimage_cosine_dist_single(probability_matrix, hap1_vec, num_samples = 1000):

    target_columns = ["Prob_0", "Prob_1", "Prob_2"]
    column_indeces = [probability_matrix.columns.index(col) for col in target_columns]

    probability_matrix_np = probability_matrix.to_numpy()
    # Ensure the probabilities are normalized
    probability_matrix_np = probability_matrix_np / probability_matrix_np.sum(axis=1, keepdims=True)

    # Generate samples
    pD_samples = generate_samples(probability_matrix_np, num_samples)

    hap1_samples = (pD_samples == column_indeces[1]).astype(int) + \
               2 * (pD_samples == column_indeces[2]).astype(int)

    # Compute cosine distances
    cos_dist1 = np.mean(1 - (hap1_samples @ hap1_vec) / (np.sqrt(np.sum(hap1_samples**2, axis=1)) * np.sqrt(np.sum(hap1_vec**2))))

    # Print cosine distances
    return(round(cos_dist1, 8))


def match_cluster_haplotype(kmer_count_file, output_prefix, kmer_info_file, cluster_kmer_count_file, depth,
    cluster_haplotype_file, max_copy_number=None, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0_unit = 0.8 / 30, nbinom_mu_unit = 0.4):

    cluster_marker_count_df = pl.read_csv(cluster_kmer_count_file, infer_schema_length = None, separator = '\t')

    # Auto-detect max_copy_number from Count_* columns if not specified
    if max_copy_number is None:
        count_cols = [col for col in cluster_marker_count_df.columns if col.startswith("Count_")]
        max_copy_number = len(count_cols) - 1

    max_depth_thres = math.ceil(depth * 3)
    # max_depth_thres = 100

    # prob_0: background distribution
    nbinom_mu_0 = nbinom_mu_0_unit * depth
    nbinom_p_0 = nbinom_size_0 / (nbinom_size_0 + nbinom_mu_0)
    prob_list = [[nbinom.pmf(x, nbinom_size_0, nbinom_p_0) for x in range(max_depth_thres + 1)]]
    # prob_1 through prob_{2*max_copy_number}: signal distributions
    for i in range(1, 2 * max_copy_number + 1):
        nbinom_p = nbinom_size / (nbinom_size + i * nbinom_mu_unit * depth)
        prob_list.append([nbinom.pmf(x, nbinom_size, nbinom_p) for x in range(max_depth_thres + 1)])

    count = pl.read_csv(kmer_count_file, separator = '\t', new_columns = ["Marker", "Count"]) \
            .filter(pl.col("Count") <= max_depth_thres)

    D_count_dict = {"Marker": count["Marker"]}
    for i, prob in enumerate(prob_list):
        D_count_dict[f"dprob_{i}"] = [prob[c] for c in count["Count"]]
    D_count = pl.DataFrame(D_count_dict)


    marker_list =  pl.Series(cluster_marker_count_df \
            .filter((pl.col("Marker_num") >= 2) & (pl.col("Hap_minus_marker_num") >= 2)) \
            .with_columns([pl.col("Rel_pos_std").cast(pl.Float64)]) \
            .select("Marker").unique()
            # .filter((pl.col("Rel_pos_mean") > -0.1) & (pl.col("Rel_pos_mean") < 1.1) & (pl.col("Rel_pos_std") < 0.5)) \
            # .select("Marker").unique()
        )

    cluster_num = cluster_marker_count_df["Cluster"].max()

    count_cols = [pl.col(f"Count_{i}") for i in range(max_copy_number + 1)]
    total = sum(count_cols) + (max_copy_number + 1) * pseudo_count
    prob_exprs = [((pl.col(f"Count_{i}") + pseudo_count) / total).alias(f"prob_{i}")
                  for i in range(max_copy_number + 1)]
    cluster_marker_count_df_2 = cluster_marker_count_df.with_columns(prob_exprs)

    # cluster_marker_count_df_2.write_csv("cluster_marker_count_df_2.tsv", separator = '\t')

    cl1_list = []
    cl2_list = []
    LL_list = []

    LL_max = -float("Inf")
    PR1_max = None
    PR2_max = None

    for i in range(1, cluster_num + 1):
        for j in range(i, cluster_num + 1):
     
            PR1 = cluster_marker_count_df_2.filter(pl.col("Cluster") == i)
            PR2 = cluster_marker_count_df_2.filter(pl.col("Cluster") == j)

            tL = calc_loglikelihood(PR1, PR2, D_count, marker_list)
            cl1_list.append(i)
            cl2_list.append(j)
            LL_list.append(tL)
            # print(f"({i}, {j}, {tL})")

            if tL > LL_max:
                PR1_max = PR1
                PR2_max = PR2
                LL_max = tL


    D_LL = pl.DataFrame({"Cluster1": cl1_list, "Cluster2": cl2_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)

    D_LL.write_csv(output_prefix + ".cluster.hap_pair.txt", separator = '\t')

    # for debug
    # PR1_max.write_csv("PR1_max.tsv", separator = '\t')
    # PR2_max.write_csv("PR2_max.tsv", separator = '\t')

    prob_cols = [f"Prob_{i}{j}" for i in range(max_copy_number + 1)
                 for j in range(max_copy_number + 1)]
    PR12_count = calc_posterior_prob(PR1_max, PR2_max, D_count, max_copy_number) \
                    .select(["Marker", "Marker_num", "Hap_minus_marker_num", "Rel_pos_mean", "Rel_pos_std"] + prob_cols)


    PR12_count.write_csv(output_prefix + ".cluster.marker_prob.txt", separator = '\t')


    if cluster_haplotype_file is None: return


    hap_marker_count_df_2_tmp = pl.read_csv(kmer_info_file, infer_schema_length = None, separator = '\t') \
        .select("Marker", "Haplotype") \
        .group_by("Marker", "Haplotype") \
        .agg(pl.len().alias("Count"))

    unique_markers = hap_marker_count_df_2_tmp.select("Marker").unique()
    unique_haplotypes = hap_marker_count_df_2_tmp.select("Haplotype").unique()

    # Create all combinations of Markers and Haplotypes
    all_combinations = (
        unique_markers.join(unique_haplotypes, how="cross")
    )

    # Count_0 through Count_N: exact match
    count_exprs = [(pl.col("Count") == i).cast(pl.Int64).alias(f"Count_{i}")
                   for i in range(max_copy_number + 1)]

    hap_count_cols = [pl.col(f"Count_{i}") for i in range(max_copy_number + 1)]
    hap_total = sum(hap_count_cols)
    hap_prob_exprs = [(pl.col(f"Count_{i}") / hap_total).alias(f"prob_{i}")
                      for i in range(max_copy_number + 1)]

    hap_marker_count_df_2 = (
        all_combinations.join(hap_marker_count_df_2_tmp, on=["Marker", "Haplotype"], how="left")
        .select(["Marker", "Haplotype", pl.col("Count").fill_null(0)])
        .with_columns(count_exprs)
        .drop("Count")
        .with_columns(hap_prob_exprs)
    )


    target_cluster_1 = D_LL[0, "Cluster1"]
    target_cluster_2 = D_LL[0, "Cluster2"]

    prob_c_cols = ["Marker"] + [pl.col(f"prob_{i}").alias(f"prob_{i}_c") for i in range(max_copy_number + 1)]
    PR1_c = cluster_marker_count_df_2 \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select(prob_c_cols)

    PR2_c = cluster_marker_count_df_2 \
        .filter(pl.col("Cluster") == target_cluster_2) \
        .select(prob_c_cols)  

    # marker_list.to_frame().write_csv("marker_list.tsv", separator = '\t')
    # PR2_c.write_csv("PR1_c.tsv", separator = '\t')
    # hap_marker_count_df_2.write_csv("hap_marker_count_df_2.tsv", separator = '\t')

    cluster_haplotype_list = pl.read_csv(cluster_haplotype_file, separator = '\t')

    target_cluster_hap_1 = cluster_haplotype_list \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select("Haplotype") \
        .to_series() \
        .to_list()

    target_cluster_hap_2 = cluster_haplotype_list \
        .filter(pl.col("Cluster") == target_cluster_2) \
        .select("Haplotype") \
        .to_series() \
        .to_list()


    cl1_list = []
    cl2_list = []
    LL_list = []

    LL_max = -float("Inf")
    PR1_max = None
    PR2_max = None


    for i in range(len(target_cluster_hap_1)):
        for j in range(len(target_cluster_hap_2)):

            if target_cluster_hap_1[i] not in [None, "NA", "None", "NONE"]:

                PR1 = hap_marker_count_df_2 \
                    .filter(pl.col("Haplotype") == target_cluster_hap_1[i]) \
                    .select(["Marker"] + [pl.col(f"prob_{k}").alias(f"prob_{k}_h") for k in range(max_copy_number + 1)]) \
                    .join(PR1_c, on="Marker", how="inner") \
                    .with_columns([((1 - cluster_ratio) * pl.col(f"prob_{k}_h") + cluster_ratio * pl.col(f"prob_{k}_c")).alias(f"prob_{k}")
                                   for k in range(max_copy_number + 1)])

            else:
                PR1 = PR1_c \
                    .select(["Marker"] + [pl.col(f"prob_{k}_c").alias(f"prob_{k}") for k in range(max_copy_number + 1)])

            if target_cluster_hap_2[j] not in [None, "NA", "None", "NONE"]:

                PR2 = hap_marker_count_df_2 \
                    .filter(pl.col("Haplotype") == target_cluster_hap_2[j]) \
                    .select(["Marker"] + [pl.col(f"prob_{k}").alias(f"prob_{k}_h") for k in range(max_copy_number + 1)]) \
                    .join(PR2_c, on="Marker", how="inner") \
                    .with_columns([((1 - cluster_ratio) * pl.col(f"prob_{k}_h") + cluster_ratio * pl.col(f"prob_{k}_c")).alias(f"prob_{k}")
                                   for k in range(max_copy_number + 1)])

            else:
                PR2 = PR2_c \
                    .select(["Marker"] + [pl.col(f"prob_{k}_c").alias(f"prob_{k}") for k in range(max_copy_number + 1)])

     
            tL = calc_loglikelihood(PR1, PR2, D_count, marker_list)
            cl1_list.append(target_cluster_hap_1[i])
            cl2_list.append(target_cluster_hap_2[j])
            LL_list.append(tL)
            # print(f"({target_cluster_hap_1[i]}, {target_cluster_hap_2[j]}, {tL})")

            if tL > LL_max:
                PR1_max = PR1
                PR2_max = PR2
                LL_max = tL


    D_LL2 = pl.DataFrame({"Haplotype1": cl1_list, "Haplotype2": cl2_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)

    D_LL2.write_csv(output_prefix + ".haplotype.hap_pair.txt", separator = '\t')


    hap_info = pl.read_csv(kmer_info_file, separator = '\t', infer_schema_length = None) \
        .filter(pl.col("Haplotype").is_in([D_LL2[0, "Haplotype1"], D_LL2[0, "Haplotype2"]])) \
        .group_by(["Marker", "Haplotype", "Marker_num"]) \
        .agg([
            pl.col("Marker_pos").mean().alias("Mean_marker_pos"),
            pl.col("Marker_pos").count().alias("Marker_count")
        ])

    hap_info1 = (
        hap_info.filter(pl.col("Haplotype") == D_LL2[0, "Haplotype1"]) 
        .rename({"Haplotype": "Haplotype1", "Mean_marker_pos": "Mean_marker_pos1", "Marker_count": "Marker_count1"})
    )

    hap_info2 = (
        hap_info.filter(pl.col("Haplotype") == D_LL2[0, "Haplotype2"]) 
        .rename({"Haplotype": "Haplotype2", "Mean_marker_pos": "Mean_marker_pos2", "Marker_count": "Marker_count2"})
    )



    prob_cols = [f"Prob_{i}{j}" for i in range(max_copy_number + 1)
                 for j in range(max_copy_number + 1)]
    PR12_count = calc_posterior_prob(PR1_max, PR2_max, D_count, max_copy_number) \
        .join(hap_info1, on = "Marker", how = "left") \
        .join(hap_info2, on = "Marker", how = "left") \
        .select(["Marker", pl.col("Haplotype1").fill_null("NA"), pl.col("Mean_marker_pos1").fill_null("NA"), pl.col("Marker_count1").fill_null("NA"),
                 pl.col("Haplotype2").fill_null("NA"), pl.col("Mean_marker_pos2").fill_null("NA"), pl.col("Marker_count2").fill_null("NA")] + prob_cols)

    # PR1_max.join(PR2_max, on="Marker", suffix = "_PR2").write_csv("PR12_max.tsv", separator = '\t')
    # PR1_max.write_csv("PR1_max_hap.tsv", separator = '\t')
    # PR2_max.write_csv("PR2_max_hap.tsv", separator = '\t')
    PR12_count.write_csv(output_prefix + ".haplotype.marker_prob.txt", separator = '\t')

    prob_mat = PR12_count.select(prob_cols)

    # Create hap1_vec and hap2_vec examples
    hap1_vec = PR12_count.select(["Marker_count1"]).to_numpy().flatten()
    hap1_vec = np.where(hap1_vec == 'NA', 0, hap1_vec).astype(float)
    hap2_vec = PR12_count.select(["Marker_count2"]).to_numpy().flatten()
    hap2_vec = np.where(hap2_vec == 'NA', 0, hap2_vec).astype(float)
    
    ##########
    # this is experimental
    # cdist1, cdist2 = estimage_cosine_dist(prob_mat, hap1_vec, hap2_vec)
    # 
    # with open(output_prefix + ".haplotype.cosine_dist.txt", 'w') as hout:
    #     print("Cosine_dist1\tCosine_dist2", file = hout)
    #     print(f'{cdist1}\t{cdist2}', file = hout)
    ##########


def match_cluster_haplotype_single(kmer_count_file, output_prefix, kmer_info_file, cluster_kmer_count_file, depth,
    cluster_haplotype_file, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0_unit = 0.8 / 30, nbinom_mu_unit = 0.4):

    max_depth_thres = math.ceil(depth * 1.5)
    # max_depth_thres = 200

    nbinom_mu_0 = nbinom_mu_0_unit * depth
    nbinom_p_0 = nbinom_size_0 / (nbinom_size_0 + nbinom_mu_0)
    prob_0 = [nbinom.pmf(x, nbinom_size_0, nbinom_p_0) for x in range(max_depth_thres + 1)]
    nbinom_p_1 = nbinom_size / (nbinom_size + 1.0 * nbinom_mu_unit * depth)
    prob_1 = [nbinom.pmf(x, nbinom_size, nbinom_p_1) for x in range(max_depth_thres + 1)]
    nbinom_p_2 = nbinom_size / (nbinom_size + 2.0 * nbinom_mu_unit * depth)
    prob_2 = [nbinom.pmf(x, nbinom_size, nbinom_p_2) for x in range(max_depth_thres + 1)]


    count = pl.read_csv(kmer_count_file, separator = '\t', new_columns = ["Marker", "Count"]) \
            .filter(pl.col("Count") <= max_depth_thres)

    D_count = pl.DataFrame({
        "Marker": count["Marker"],
        "dprob_0": [prob_0[c] for c in count["Count"]],
        "dprob_1": [prob_1[c] for c in count["Count"]],
        "dprob_2": [prob_2[c] for c in count["Count"]],
    })

    cluster_marker_count_df = pl.read_csv(cluster_kmer_count_file, infer_schema_length = None, separator = '\t')


    marker_list =  pl.Series(cluster_marker_count_df \
            .filter((pl.col("Marker_num") >= 2) & (pl.col("Hap_minus_marker_num") >= 2)) \
            .with_columns([pl.col("Rel_pos_std").cast(pl.Float64)]) \
            .filter((pl.col("Rel_pos_mean") > -0.1) & (pl.col("Rel_pos_mean") < 1.1) & (pl.col("Rel_pos_std") < 0.5)) \
            .select("Marker").unique()
        )

    cluster_num = cluster_marker_count_df["Cluster"].max()

    cluster_marker_count_df_2 = cluster_marker_count_df.with_columns([
        ((pl.col("Count_0") + pseudo_count) /
         (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2") + 3 * pseudo_count)).alias("prob_0"),
        ((pl.col("Count_1") + pseudo_count) /
         (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2") + 3 * pseudo_count)).alias("prob_1"),
        ((pl.col("Count_2") + pseudo_count) /
         (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2") + 3 * pseudo_count)).alias("prob_2")
    ])

    # cluster_marker_count_df_2.write_csv("cluster_marker_count_df_2.tsv", separator = '\t')

    cl1_list = []
    LL_list = []

    LL_max = -float("Inf")
    PR1_max = None


    for i in range(1, cluster_num + 1):

        PR1 = cluster_marker_count_df_2.filter(pl.col("Cluster") == i)

        tL = calc_loglikelihood_single(PR1, D_count, marker_list)
        cl1_list.append(i)
        LL_list.append(tL)

        if tL > LL_max:
            PR1_max = PR1
            LL_max = tL


    D_LL = pl.DataFrame({"Cluster1": cl1_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)
    D_LL.write_csv(output_prefix + ".cluster.hap_pair.txt", separator = '\t')

    PR1_count = calc_posterior_prob_single(PR1_max, D_count) \
                    .select(["Marker", "Marker_num", "Hap_minus_marker_num", "Rel_pos_mean", "Rel_pos_std",
                    "Prob_0", "Prob_1", "Prob_2"])

    PR1_count.write_csv(output_prefix + ".cluster.marker_prob.txt", separator = '\t')


    if cluster_haplotype_file is None: return


    hap_marker_count_df_2_tmp = pl.read_csv(kmer_info_file, infer_schema_length = None, separator = '\t') \
        .select("Marker", "Haplotype") \
        .group_by("Marker", "Haplotype") \
        .agg(pl.len().alias("Count"))

    unique_markers = hap_marker_count_df_2_tmp.select("Marker").unique()
    unique_haplotypes = hap_marker_count_df_2_tmp.select("Haplotype").unique()

    # Create all combinations of Markers and Haplotypes
    all_combinations = (
        unique_markers.join(unique_haplotypes, how="cross")
    )

    hap_marker_count_df_2 = (
        all_combinations.join(hap_marker_count_df_2_tmp, on=["Marker", "Haplotype"], how="left")
        .select(["Marker", "Haplotype", pl.col("Count").fill_null(0)])
        .with_columns([
        # Creating binary columns based on `Count` values
            (pl.col("Count") == 0).cast(pl.Int64).alias("Count_0"),
            (pl.col("Count") == 1).cast(pl.Int64).alias("Count_1"),
            (pl.col("Count") == 2).cast(pl.Int64).alias("Count_2"),
        ])
        .drop("Count")
        .with_columns([
            (pl.col("Count_0") / (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2"))).alias("prob_0"),
            (pl.col("Count_1") / (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2"))).alias("prob_1"),
            (pl.col("Count_2") / (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2"))).alias("prob_2"),
        ])

    )


    target_cluster_1 = D_LL[0, "Cluster1"]

    PR1_c = cluster_marker_count_df_2 \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select([
            "Marker",
            pl.col("prob_0").alias("prob_0_c"),
            pl.col("prob_1").alias("prob_1_c"),
            pl.col("prob_2").alias("prob_2_c")
        ])

    cluster_haplotype_list = pl.read_csv(cluster_haplotype_file, separator = '\t')

    target_cluster_hap_1 = cluster_haplotype_list \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select("Haplotype") \
        .to_series() \
        .to_list()

    cl1_list = []
    LL_list = []

    LL_max = -float("Inf")
    PR1_max = None


    for i in range(len(target_cluster_hap_1)):

        if target_cluster_hap_1[i] not in [None, "NA", "None", "NONE"]:

            PR1 = hap_marker_count_df_2 \
                .filter(pl.col("Haplotype") == target_cluster_hap_1[i]) \
                .select([
                    "Marker",
                    pl.col("prob_0").alias("prob_0_h"),
                    pl.col("prob_1").alias("prob_1_h"),
                    pl.col("prob_2").alias("prob_2_h")
                ]) \
                .join(PR1_c, on="Marker", how="inner") \
                .with_columns([
                    ((1 - cluster_ratio) * pl.col("prob_0_h") + cluster_ratio * pl.col("prob_0_c")).alias("prob_0"),
                    ((1 - cluster_ratio) * pl.col("prob_1_h") + cluster_ratio * pl.col("prob_1_c")).alias("prob_1"),
                    ((1 - cluster_ratio) * pl.col("prob_2_h") + cluster_ratio * pl.col("prob_2_c")).alias("prob_2")
                ])

        else:
            PR1 = PR1_c \
                .select([
                    "Marker",
                    pl.col("prob_0_c").alias("prob_0"),
                    pl.col("prob_1_c").alias("prob_1"),
                    pl.col("prob_2_c").alias("prob_2")
                ])

        tL = calc_loglikelihood_single(PR1, D_count, marker_list)
        cl1_list.append(target_cluster_hap_1[i])
        LL_list.append(tL)

        if tL > LL_max:
            PR1_max = PR1
            LL_max = tL


    D_LL2 = pl.DataFrame({"Haplotype1": cl1_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)

    D_LL2.write_csv(output_prefix + ".haplotype.hap_pair.txt", separator = '\t')



    hap_info = pl.read_csv(kmer_info_file, separator = '\t', infer_schema_length = None) \
        .filter(pl.col("Haplotype").is_in([D_LL2[0, "Haplotype1"]])) \
        .group_by(["Marker", "Haplotype", "Marker_num"]) \
        .agg([
            pl.col("Marker_pos").mean().alias("Mean_marker_pos"),
            pl.col("Marker_pos").count().alias("Marker_count")
        ])

    hap_info1 = (
        hap_info.filter(pl.col("Haplotype") == D_LL2[0, "Haplotype1"])
        .rename({"Haplotype": "Haplotype1", "Mean_marker_pos": "Mean_marker_pos1", "Marker_count": "Marker_count1"})
    )

    PR1_count = calc_posterior_prob_single(PR1_max, D_count) \
        .join(hap_info1, on = "Marker", how = "left") \
        .select(["Marker", pl.col("Haplotype1").fill_null("NA"), pl.col("Mean_marker_pos1").fill_null("NA"), pl.col("Marker_count1").fill_null("NA"),
                 "Prob_0", "Prob_1", "Prob_2"])

    # for debug
    # PR1_max.write_csv("PR1_max.tsv", separator = '\t')
    # D_count.write_csv("D_count.tsv", separator = '\t')

    # PR1_max.join(PR2_max, on="Marker", suffix = "_PR2").write_csv("PR12_max.tsv", separator = '\t')
    # PR1_max.write_csv("PR1_max_hap.tsv", separator = '\t')
    # PR2_max.write_csv("PR2_max_hap.tsv", separator = '\t')
    PR1_count.write_csv(output_prefix + ".haplotype.marker_prob.txt", separator = '\t')


    prob_mat = PR1_count.select(["Prob_0", "Prob_1", "Prob_2"])

    hap1_vec = PR1_count.select(["Marker_count1"]).to_numpy().flatten()
    hap1_vec = np.where(hap1_vec == 'NA', 0, hap1_vec).astype(float)

    ##########
    # this is experimental
    # cdist1 = estimage_cosine_dist_single(prob_mat, hap1_vec)
    # 
    # with open(output_prefix + ".haplotype.cosine_dist.txt", 'w') as hout:
    #     print("Cosine_dist", file = hout)
    #     print(f'{cdist1}', file = hout)
    ##########



