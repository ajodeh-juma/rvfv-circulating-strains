import argparse
import json
import re
import numpy as np


def calc_weighted_quantiles(values, weights):
    """Calculates weighted median and 95% credible interval (2.5% and 97.5% percentiles)."""
    sort_order = np.argsort(values)
    sorted_vals = values[sort_order]
    cum_weights = np.cumsum(weights[sort_order])

    lo = sorted_vals[np.searchsorted(cum_weights, 0.025)]
    med = sorted_vals[np.searchsorted(cum_weights, 0.50)]
    hi = sorted_vals[np.searchsorted(cum_weights, 0.975)]

    return med, lo, hi


def parse_fubar_log(log_path):
    """Parses tabulated site-level results from a HyPhy FUBAR log file and calculates omega (dN/dS)."""
    site_log_data = {}
    with open(log_path, "r") as f:
        in_table = False
        for line in f:
            if "Codon" in line and "Partition" in line and "alpha" in line:
                in_table = True
                continue
            if in_table:
                if line.strip().startswith("|:") or line.strip().startswith("|---"):
                    continue
                if not line.strip().startswith("|"):
                    in_table = False
                    continue
                cols = [c.strip() for c in line.split("|")[1:-1]]
                if len(cols) >= 5:
                    try:
                        codon = int(cols[0])
                        partition = cols[1]
                        alpha = float(cols[2])
                        beta = float(cols[3])
                        pos_post_str = cols[4]

                        # Calculate point-estimate dN/dS (omega)
                        omega = beta / alpha if alpha > 0 else 0.0

                        # Extract posterior probability from 'Pos. posterior = X.XXXX'
                        match = re.search(r"=\s*([0-9.]+)", pos_post_str)
                        pos_prob = float(match.group(1)) if match else None

                        site_log_data[codon] = {
                            "partition": partition,
                            "alpha": alpha,
                            "beta": beta,
                            "omega": omega,
                            "pos_prob": pos_prob,
                            "raw_text": pos_post_str,
                        }
                    except ValueError:
                        pass
    return site_log_data


def main():
    parser = argparse.ArgumentParser(
        description="Extract alpha, beta, omega (dN/dS), and posterior probabilities from HyPhy FUBAR JSON and log files."
    )
    parser.add_argument("--json", required=True, help="Path to the FUBAR JSON file")
    parser.add_argument("--log", required=False, help="Path to the FUBAR log file (optional)")
    parser.add_argument("--site", type=int, default=278, help="1-indexed codon site number (default: 278)")

    args = parser.parse_args()

    # --- 1. Extract Grid & Posterior Distributions from JSON ---
    with open(args.json, "r") as f:
        json_data = json.load(f)

    grid = json_data["grid"]
    posterior_raw = json_data.get("posterior") or json_data.get("posteriors")

    if isinstance(posterior_raw, dict):
        if "0" in posterior_raw:
            post_partition = posterior_raw["0"]
        elif 0 in posterior_raw:
            post_partition = posterior_raw[0]
        else:
            post_partition = posterior_raw
    else:
        post_partition = posterior_raw

    site_idx = args.site - 1
    weights = None
    matched_key = None

    if isinstance(post_partition, dict):
        possible_keys = [str(site_idx), str(args.site), site_idx, args.site]
        for k in possible_keys:
            if k in post_partition:
                weights = post_partition[k]
                matched_key = k
                break

        if weights is None:
            raise KeyError(
                f"Codon site {args.site} (index '{site_idx}') not found in JSON partition dictionary."
            )
    elif isinstance(post_partition, list):
        if 0 <= site_idx < len(post_partition):
            weights = post_partition[site_idx]
            matched_key = site_idx
        else:
            raise IndexError(
                f"Codon site {args.site} is out of bounds for sequence length {len(post_partition)}."
            )

    grid_arr = np.array(grid)
    alpha_grid = grid_arr[:, 0]
    beta_grid = grid_arr[:, 1]
    omega_grid = np.where(alpha_grid > 0, beta_grid / alpha_grid, 0.0)

    weights_flat = np.array(weights).flatten()
    weights_flat = weights_flat / weights_flat.sum()

    alpha_med, alpha_lo, alpha_hi = calc_weighted_quantiles(alpha_grid, weights_flat)
    beta_med, beta_lo, beta_hi = calc_weighted_quantiles(beta_grid, weights_flat)
    omega_med, omega_lo, omega_hi = calc_weighted_quantiles(omega_grid, weights_flat)

    print(f"\n=========================================")
    print(f"  FUBAR Results for Codon Site {args.site}")
    print(f"=========================================")
    print("\n[JSON Grid Posterior Distributions]")
    print(f"  α (synonymous rate):     {alpha_med:.4f}  [95% CI: {alpha_lo:.4f} - {alpha_hi:.4f}]")
    print(f"  β (non-synonymous rate):  {beta_med:.4f}  [95% CI: {beta_lo:.4f} - {beta_hi:.4f}]")
    print(f"  ω (dN/dS):              {omega_med:.4f}  [95% CI: {omega_lo:.4f} - {omega_hi:.4f}]")

    # --- 2. Extract Point Estimates & Selection Metrics from Log File ---
    if args.log:
        log_data = parse_fubar_log(args.log)
        print("\n[Log File Tabulated Summary]")
        if args.site in log_data:
            site_info = log_data[args.site]
            print(f"  Partition:                            {site_info['partition']}")
            print(f"  α (point estimate):                  {site_info['alpha']:.4f}")
            print(f"  β (point estimate):                  {site_info['beta']:.4f}")
            print(f"  ω (dN/dS point estimate):            {site_info['omega']:.4f}")
            if site_info["pos_prob"] is not None:
                print(f"  Posterior Prob (Positive Selection): {site_info['pos_prob']:.4f}")
            else:
                print(f"  Selection Result:                    {site_info['raw_text']}")
        else:
            print(f"  Codon site {args.site} is not listed in the log file table.")
            print("  (Note: HyPhy log tables only report sites exceeding threshold Prob(β > α) ≥ 0.9)")


if __name__ == "__main__":
    main()