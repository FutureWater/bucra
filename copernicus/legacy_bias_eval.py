rice_threshold = metrics.ThresholdMetric(
    name="Rice Threshold",
    variable="tas",
    threshold_value=[293, 308],
    threshold_type="outside",
)

marginal_bias_data = marginal.calculate_marginal_bias(
    metrics=[rice_threshold],
    statistics=["mean", 0.05, 0.95],
    percentage_or_absolute="absolute",
    obs=reanalysis_validation_np,
    QM=QM_val[0],
    raw=projections_validation_np[0],
)

tas_marginal_bias_plot = marginal.plot_marginal_bias(
    variable="tas",
    bias_df=marginal_bias_data,
    metrics_title="Absolute bias [days / year]",
    statistics_title="Absolute bias [K]",
)

tas_bias_map_mean = marginal.plot_bias_spatial(
    variable="tas", metric="Rice Threshold", bias_df=marginal_bias_data
)
# %%
i = 5
print(rice_threshold.calculate_instances_of_threshold_exceedance(QM_val[i]).sum())
print(
    rice_threshold.calculate_instances_of_threshold_exceedance(
        projections_validation_np[i]
    ).sum()
)
print(
    f"Difference (corrected vs uncorrected): {rice_threshold.calculate_instances_of_threshold_exceedance(projections_validation_np[i]).sum()-rice_threshold.calculate_instances_of_threshold_exceedance(QM_val[i]).sum()}"
)
print(
    rice_threshold.calculate_instances_of_threshold_exceedance(
        reanalysis_validation_np
    ).sum()
)
print(
    f"Difference (corrected vs reality): {rice_threshold.calculate_instances_of_threshold_exceedance(QM_val[i]).sum()-rice_threshold.calculate_instances_of_threshold_exceedance(reanalysis_validation_np).sum()}"
)
