import copy


def bonferroni_correction(enrichment_results):
    """Performs Bonferroni correction on the p-values in the enrichment_results dictionary to adjust for multiple hypothesis testing."""

    adjusted_results = copy.deepcopy(enrichment_results)

    # Extract raw p-values
    class_p_list = [(cls, res["p_value"]) for cls, res in enrichment_results.items()]

    m = len(class_p_list)  # number of tests
    correction_map = {}

    for cls, raw_p in class_p_list:
        if raw_p is None:
            corrected_p = None
        else:
            corrected_p = min(raw_p * m, 1.0)  # Bonferroni correction

        adjusted_results[cls]["p_value_corrected"] = corrected_p
        correction_map[cls] = corrected_p

    return adjusted_results, correction_map


def benjamini_hochberg_fdr_correction(enrichment_results):
    """Performs Benjamini-Hochberg FDR (false discovery rate) correction on the p-values in the enrichment_results dictionary to adjust for multiple hypothesis testing."""
    # Extract p-values into a sorted list, filtering out None values
    p_values_with_class = [
        (cls, info["p_value"])
        for cls, info in enrichment_results.items()
        if info["p_value"] is not None
    ]
    p_values_with_class.sort(key=lambda x: x[1])  # Sort by p-value ascending

    m = len(p_values_with_class)  # number of non-None tests
    if m == 0:
        # No valid p-values to correct
        for cls in enrichment_results:
            enrichment_results[cls]["p_value_corrected"] = None
        return enrichment_results

    adj = [None] * m  # list for adjusted p-values

    # Apply Benjamini-Hochberg procedure
    running_min = 1.0  # p-value cannot be higher than 1
    for i in range(m, 0, -1):
        raw_p = p_values_with_class[i - 1][1]
        bh_value = (raw_p * m) / i
        running_min = min(running_min, bh_value)
        adj[i - 1] = running_min

    # Map adjusted p-values back into the results dictionary
    for (cls, _), corrected_p in zip(p_values_with_class, adj):
        enrichment_results[cls]["p_value_corrected"] = corrected_p

    # Mark classes with None p-values as having None corrected p-values
    for cls, info in enrichment_results.items():
        if info["p_value"] is None:
            info["p_value_corrected"] = None

    return enrichment_results
