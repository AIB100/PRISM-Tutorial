#!/usr/bin/env python3
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


def _load_distribution_builder():
    repo_root = Path(__file__).resolve().parents[4] / 'PRISM'
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from prism.fep.gromacs.mdp_templates import _generate_distribution
    return _generate_distribution


def main() -> None:
    build_distribution = _load_distribution_builder()
    n_points = 32

    distributions = [
        ('Linear', build_distribution(n_points, 'linear'), '#2E86C1'),
        ('Nonlinear (default)', build_distribution(n_points, 'nonlinear'), '#C0392B'),
        ('Quadratic (p=2)', build_distribution(n_points, 'quadratic', quadratic_exponent=2.0), '#27AE60'),
        ('Quadratic (p=4)', build_distribution(n_points, 'quadratic', quadratic_exponent=4.0), '#8E44AD'),
    ]

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 14

    fig, ax = plt.subplots(1, 1, figsize=(8.5, 6.5))
    x = np.arange(n_points)

    markers = ['o', 's', '^', 'D']
    for (title, values, color), marker in zip(distributions, markers):
        ax.plot(x, values, marker=marker, linewidth=2.8, markersize=4.8, color=color, label=title)

    ax.set_title('PRISM lambda point distributions', fontweight='bold', fontsize=28, pad=14)
    ax.set_xlabel('Point index', fontsize=20)
    ax.set_ylabel('Lambda value', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(-0.5, len(x) - 0.5)
    ax.set_xticks(np.linspace(0, len(x) - 1, 5, dtype=int))
    ax.grid(True, alpha=0.25, linestyle='--')
    ax.legend(loc='lower right', frameon=False, fontsize=16)

    plt.tight_layout(rect=[0, 0.02, 1, 0.98])

    out = Path(__file__).with_name('lambda_distributions.png')
    fig.savefig(out, dpi=220, bbox_inches='tight')
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
