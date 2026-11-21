"""Post hoc compare the prespecified models to progression status alone."""
import pandas as pd
from sklearn.model_selection import StratifiedKFold

from analyze import OUT, SEED, cohort_data, fit_predict, horizon_data, log_model, metrics
from prepare_data import load_all

frames, _ = cohort_data(load_all())
development = frames["TCGA"]
fit = fit_predict(log_model(), development, frames, ["early"], 12)
y, weights, _, _ = horizon_data(development, 12)
split = StratifiedKFold(5, shuffle=True, random_state=SEED)
for train, test in split.split(development, y + 2 * (weights == 0)):
    fit["TCGA"][test] = fit_predict(log_model(), development.iloc[train],
                                   {"test": development.iloc[test]}, ["early"], 12)["test"]
rows = [{"cohort": name, **metrics(data, fit[name], 12)} for name, data in frames.items()]
pd.DataFrame(rows).to_csv(OUT / "progression_only_benchmark.csv", index=False)
