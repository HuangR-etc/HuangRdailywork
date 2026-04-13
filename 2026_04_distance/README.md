# Phylogenetic Dispersed Subset Analysis Framework

A comprehensive R framework for constructing phylogenetically dispersed subsets on phylogenetic trees and performing four-part analysis as described in the design document.

## Overview

This framework implements methods for selecting phylogenetically dispersed subsets from phylogenetic trees using a lexicographic multi-objective optimization approach. The main algorithm maximizes three metrics in order of priority:
1. **MinPD** (Minimum Pairwise Distance)
2. **MeanPD** (Mean Pairwise Distance)
3. **MeanNND** (Mean Nearest Neighbor Distance)

## Project Structure

```
projects/2026_04_distance/
├── README.md                    # This file
├── main.R                       # Main analysis program
├── test_implementation.R        # Test script
├── config/
│   └── analysis_config.R        # Configuration settings
├── R/
│   ├── tree_generators.R        # Tree generation functions
│   ├── distance_metrics.R       # Distance metric calculations
│   ├── objective_compare.R      # Lexicographic comparison
│   ├── subset_greedy.R          # Greedy construction algorithm
│   ├── subset_exchange.R        # Exchange refinement algorithm
│   ├── subset_random.R          # Random subset sampling
│   ├── subset_exhaustive.R      # Exhaustive search (for small trees)
│   ├── trait_simulation.R       # Trait simulation functions
│   ├── signal_metrics.R         # Phylogenetic signal calculations
│   ├── result1_analysis.R       # Result 1: Comparison with random
│   ├── result2_analysis.R       # Result 2: Design component analysis
│   ├── result3_analysis.R       # Result 3: Trait-level consequences
│   ├── result4_analysis.R       # Result 4: Heuristic vs exact optimum
│   ├── plotting.R               # Visualization functions
│   └── utils_io.R               # Utility and I/O functions
└── outputs/                     # Generated outputs (created at runtime)
    ├── tables/                  # CSV tables with results
    ├── figures/                 # PDF figures
    ├── rds/                     # RDS objects for intermediate results
    ├── logs/                    # Analysis logs
    └── reports/                 # Summary reports
```

## Four-Part Analysis

### Result 1: Method vs Random Sampling
Tests whether the main method selects more dispersed subsets than random sampling from the same candidate pool.

**Key outputs:**
- Null distributions of three metrics
- Empirical p-values and z-scores
- Visualization of observed vs null distributions

### Result 1b: Random-Tree Replicated Analysis
Extends Result 1 by testing the method on 100 random phylogenetic trees to assess robustness across different topologies.

**Key features:**
- Generates 100 random trees with 256 tips each
- Performs Result 1 analysis on each tree independently
- Uses tree-specific null distributions (not mixed across trees)
- Summarizes results using standardized metrics (z-scores, percentiles)
- Compares with original Result 1 fixed tree results

**Key outputs:**
- Per-tree summary with observed metrics and null comparisons
- Cross-tree summary statistics
- Distribution plots of z-scores and percentiles across trees
- Proportion of trees where observed exceeds 95th percentile
- Comparison with fixed tree results from original Result 1

**Directory structure:**
```
outputs/result1_random_trees/
├── per_tree/
│   ├── tree_001/
│   │   ├── tree.rds
│   │   ├── observed_subset.csv
│   │   ├── null_metrics.csv
│   │   └── tree_summary.csv
│   ├── tree_002/
│   └── ...
├── summary/
│   ├── per_tree_summary.csv
│   ├── cross_tree_summary.csv
│   └── random_tree_config.json
└── figures/
    ├── example_tree_null_dist.pdf
    ├── zscore_boxplot.pdf
    ├── percentile_boxplot.pdf
    ├── above95_barplot.pdf
    └── compare_with_existing_result1.pdf
```

**Usage:**
```r
# Run the full random-tree analysis
source("R/result1_random_trees.R")
results <- run_result1_random_trees()

# Run with custom parameters
results <- run_result1_random_trees(
  n_trees = 100,
  n_tips = 256,
  subset_size = 20,
  n_null_reps = 1000,
  global_seed = 20260409
)

# Run test version (smaller parameters)
test_result1_random_trees()
```

### Result 2: Design Component Analysis
Evaluates the necessity of different design components:
- Exchange refinement vs greedy-only
- Multi-criterion vs single-objective optimization
- Comparison of five algorithms

**Key outputs:**
- Performance metrics for all algorithms
- Pairwise comparisons
- Improvement from exchange refinement

### Result 3: Trait-level Consequences
Investigates whether phylogenetic dispersion reduces trait phylogenetic signal.

**Key outputs:**
- Blomberg's K and Pagel's λ for different subset types
- Comparison: clustered > random > dispersed
- Statistical tests of pattern consistency

### Result 4: Heuristic vs Exact Optimum
Validates the heuristic algorithm by comparing it with exact optimum on small trees.

**Key outputs:**
- Exact optimum vs heuristic solution
- Gap analysis and relative differences
- Rank and percentile of heuristic solution

## Installation and Requirements

### Required R Packages
```r
install.packages(c(
  "ape",        # Phylogenetic trees
  "phytools",   # Phylogenetic tools
  "picante",    # Phylogenetic community analysis
  "geiger",     # Comparative methods
  "ggplot2",    # Plotting
  "gridExtra",  # Grid graphics
  "viridis",    # Color scales
  "dplyr"       # Data manipulation
))
```

The framework will automatically check for and optionally install missing packages.

## Usage

### Running the Full Analysis
```bash
cd /home/huangr/projects/2026_04_distance
Rscript main.R
```

### Running Individual Results
```bash
# Run only Result 1
Rscript main.R result1

# Run only Result 2
Rscript main.R result2

# Run only Result 3
Rscript main.R result3

# Run only Result 4
Rscript main.R result4
```

### Running Test Analysis
```bash
# Run test with smaller parameters
Rscript main.R test
```

### Interactive Testing
```r
# In R console
setwd("/home/huangr/projects/2026_04_distance")
source("test_implementation.R")
```

## Configuration

The default configuration is defined in `config/analysis_config.R`:

```r
get_default_config <- function() {
  list(
    seed = 123,
    tree_reps = 1,
    large_n = 256,
    small_n = 32,
    subset_large = 20,
    subset_small = 5,
    null_reps_large = 1000,
    null_reps_small = 1000,
    exhaustive_small = TRUE,
    trait_reps = 500,
    random_subset_reps_for_trait = 100,
    ou_alpha = c(0.2, 1, 5),
    run_result1 = TRUE,
    run_result2 = TRUE,
    run_result3 = TRUE,
    run_result4 = TRUE
  )
}
```

### Key Parameters
- **large_n**: Number of tips in large trees (default: 256)
- **small_n**: Number of tips in small tree (default: 32)
- **subset_large**: Subset size for large trees (default: 20)
- **subset_small**: Subset size for small tree (default: 5)
- **null_reps_large**: Random subsets for large trees (default: 1000)
- **null_reps_small**: Random subsets for small tree (default: 1000 or exhaustive)
- **trait_reps**: Trait simulation replicates (default: 500)
- **ou_alpha**: OU process alpha values (default: 0.2, 1, 5)

## Algorithm Details

### Main Algorithm (Lexicographic Multi-objective)
1. **Greedy Construction**: Start with random species, add species that maximizes the three metrics in lexicographic order
2. **Exchange Refinement**: Iteratively swap selected/unselected species to improve metrics
3. **Convergence**: Stop when no single exchange improves the solution

### Metrics
1. **MinPD**: Minimum pairwise patristic distance in subset
2. **MeanPD**: Mean pairwise patristic distance in subset
3. **MeanNND**: Mean nearest neighbor distance in subset

### Comparison Algorithms
- **A_main**: Full algorithm (MinPD > MeanPD > MeanNND, greedy + exchange)
- **B_greedy_only**: Greedy only (MinPD > MeanPD > MeanNND)
- **C_meanpd_only**: Maximize MeanPD only (greedy + exchange)
- **D_minpd_only**: Maximize MinPD only (greedy + exchange)
- **E_meannnd_only**: Maximize MeanNND only (greedy + exchange)

## Output Files

### Tables (`outputs/tables/`)
- `result1_summary.csv`: Summary of Result 1 analysis
- `result2_method_metrics.csv`: Algorithm performance metrics
- `result3_signal_metrics.csv`: Trait signal metrics
- `result4_exact_vs_heuristic.csv`: Heuristic vs exact comparison

### Figures (`outputs/figures/`)
- Null distribution plots for each metric
- Algorithm comparison bar plots
- Trait signal box plots
- Heuristic vs exact comparison plots

### RDS Objects (`outputs/rds/`)
- Intermediate results for reproducibility
- Distance matrices and tree objects
- Complete analysis results

### Reports (`outputs/reports/`)
- `analysis_summary.txt`: Text summary of all results

## Design Principles

1. **Modularity**: Each analysis component is in a separate, reusable module
2. **Reproducibility**: Random seeds, logging, and saved intermediate results
3. **Flexibility**: Configurable parameters and separate result execution
4. **Performance**: Efficient distance matrix calculations and subset evaluation
5. **Validation**: Comprehensive testing and comparison with exact solutions

## Testing

The framework includes comprehensive testing:

1. **Unit Tests**: Each module has built-in test functions
2. **Integration Test**: `test_implementation.R` verifies all modules work together
3. **Pipeline Test**: `Rscript main.R test` runs the full pipeline with smaller parameters

## Performance Considerations

- Distance matrices are pre-computed and cached
- Subset evaluation uses efficient matrix operations
- Exhaustive search is only used for small trees (n ≤ 32, k ≤ 5)
- Parallelization can be added for trait simulations if needed

## Extending the Framework

### Adding New Algorithms
1. Create a new function in the appropriate module
2. Add it to the algorithm comparison in `result2_analysis.R`
3. Update the plotting functions if needed

### Adding New Metrics
1. Add calculation function to `distance_metrics.R`
2. Update the lexicographic comparison in `objective_compare.R`
3. Update analysis and plotting modules

### Adding New Tree Types
1. Add generation function to `tree_generators.R`
2. Update the main analysis to include the new tree type

## Citation

If you use this framework in your research, please cite the original methodology paper (to be provided).

## License

This project is for academic research use.

## Contact

For questions or issues, please contact the maintainer.
