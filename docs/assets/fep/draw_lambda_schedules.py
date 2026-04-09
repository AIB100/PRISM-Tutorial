#!/usr/bin/env python3
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


def _load_schedule_builder():
    repo_root = Path(__file__).resolve().parents[4] / 'PRISM'
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from prism.fep.gromacs.mdp_templates import _generate_lambda_schedule

    return _generate_lambda_schedule


def main() -> None:
    build_schedule = _load_schedule_builder()

    decoupled = build_schedule(
        strategy='decoupled', distribution='nonlinear', total_windows=32, coul_windows=12, vdw_windows=20
    )
    coupled = build_schedule(strategy='coupled', distribution='nonlinear', total_windows=32)
    custom = build_schedule(
        strategy='custom',
        custom_coul=[0.0, 0.25, 0.70, 1.0],
        custom_vdw=[0.0, 0.00, 0.40, 1.0],
    )

    schedules = [
        ('Decoupled (default)', decoupled, ''),
        ('Coupled', coupled, ''),
        ('Custom (example)', custom, ''),
    ]

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 14

    fig, axes = plt.subplots(3, 1, figsize=(7.5, 14))

    colors = {
        'coul': '#C0392B',
        'vdw': '#2E86C1',
        'bonded': '#7F8C8D',
        'mass': '#27AE60',
    }

    for ax, (title, schedule, subtitle) in zip(axes, schedules):
        coul = [float(x) for x in schedule['coul_lambdas'].split()]
        vdw = [float(x) for x in schedule['vdw_lambdas'].split()]
        bonded = [float(x) for x in schedule['bonded_lambdas'].split()]
        mass = [float(x) for x in schedule['mass_lambdas'].split()]
        x = np.arange(len(coul))

        ax.plot(x, coul, marker='o', linewidth=2.5, markersize=4, color=colors['coul'], label='coul-lambdas')
        ax.plot(x, vdw, marker='s', linewidth=2.5, markersize=4, color=colors['vdw'], label='vdw-lambdas')
        ax.plot(x, bonded, linestyle='--', linewidth=1.8, color=colors['bonded'], label='bonded-lambdas')
        ax.plot(x, mass, linestyle=':', linewidth=2.0, color=colors['mass'], label='mass-lambdas')

        ax.set_title(title, fontweight='bold', fontsize=24)
        ax.set_xlabel('Window index', fontsize=20)
        ax.set_ylabel('Lambda value', fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlim(-0.5, len(x) - 0.5)
        if len(x) <= 10:
            ax.set_xticks(x)
        else:
            ax.set_xticks(np.linspace(0, len(x) - 1, 5, dtype=int))
        ax.grid(True, alpha=0.25, linestyle='--')

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.03), fontsize=21)
    fig.suptitle('PRISM FEP lambda schedule modes', fontsize=28, fontweight='bold', y=0.988)

    plt.tight_layout(rect=[0, 0.095, 1, 0.982])

    out = Path(__file__).with_name('lambda_schedule_modes.png')
    fig.savefig(out, dpi=220, bbox_inches='tight')
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
