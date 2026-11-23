"""Compare confoundvis::sens_dml() with the Python package doubleml.

Run after 03_validation.R, from this folder:  python 04_doubleml_check.py
Requires: pip install doubleml==0.11.4 scikit-learn pandas numpy
The R fits' cross-fitted predictions are passed to doubleml as external
predictions, so both packages work from the same nuisance estimates.
"""
import warnings

import doubleml as dml
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, LogisticRegression

warnings.filterwarnings("ignore")
path = "output/validation/"
data = pd.read_csv(path + "k401_data.csv")
r_res = pd.read_csv(path + "dml_r_results.csv").set_index("model")
xc = ["age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown"]
obj = dml.DoubleMLData(data, y_col="net_tfa", d_cols="e401", x_cols=xc)

specs = {
    "plr": (dml.DoubleMLPLR, dict(ml_l=LinearRegression(), ml_m=LinearRegression()),
            ["ml_l", "ml_m"], 3),
    "irm": (dml.DoubleMLIRM,
            dict(ml_g=LinearRegression(), ml_m=LogisticRegression(),
                 ps_processor_config=dml.utils.PSProcessorConfig(clipping_threshold=1e-12)),
            ["ml_g0", "ml_g1", "ml_m"], 1),
}
worst = {"bounds": 0.0, "rv": 0.0}
for name, (cls, kw, learners, n_rep) in specs.items():
    pred = pd.read_csv(path + f"{name}_predictions.csv")
    ext = {"e401": {k: pred[[f"{k}_{i + 1}" for i in range(n_rep)]].to_numpy()
                    for k in learners}}
    m = cls(obj, n_folds=5, n_rep=n_rep, **kw)
    m.fit(external_predictions=ext)
    m.sensitivity_analysis(cf_y=0.03, cf_d=0.03, rho=1.0, level=0.95)
    sp = m.sensitivity_params
    py = {"coef": np.median(m.all_coef), "rv": sp["rv"][0], "rva": sp["rva"][0],
          "theta_lower": sp["theta"]["lower"][0], "se_lower": sp["se"]["lower"][0],
          "ci_lower": sp["ci"]["lower"][0], "ci_upper": sp["ci"]["upper"][0]}
    print(f"{name.upper()} (n_rep = {n_rep})")
    for k, v in py.items():
        rv = r_res.loc[name, k]
        rel = abs(v - rv) / abs(v)
        group = "rv" if k in ("rv", "rva") else "bounds"
        worst[group] = max(worst[group], rel)
        print(f"  {k:12s} python {v:14.6f}  R {rv:14.6f}  rel. diff {rel:.1e}")
print(f"\nLargest relative difference: bounds {worst['bounds']:.1e}; "
      f"robustness values {worst['rv']:.1e}")
