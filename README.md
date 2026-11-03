# crossing-ratio

Code for the paper "Do the Answer Patterns Cross? An Alternative Method for Comparing Group Means on Likert Items."

The crossing ratio shows how much a difference in mean Likert scores between two groups depends on the assumption of equal intervals between answer categories. It is infinite when the answer patterns of the two groups do not cross.

## Files

- `crossing_ratio.R`: the function `crossing_ratio(y, group)` for one item and two groups
- `sim2_crossing.R`: the simulation (2,160 conditions, plus summed scales)
- `analysis2_crossing.R`: rates by condition and hypothesis tests
- `audit2_survey.R`: analysis of the General Social Survey and the psychTools bfi and spi data
- `PROTOCOL.md`: study protocol and hypotheses
- `results/`: rates for every condition, hypothesis results, and survey results

## Use

Copy a script into the R console. Output folders are created automatically.

```r
source("crossing_ratio.R")
y <- c(rep(3, 2), rep(4, 8), 1, rep(2, 5), 3, rep(5, 13))
group <- c(rep("men", 10), rep("women", 20))
crossing_ratio(y, group, seed = 1)
```

The function returns the crossing ratio, an approximate one-sided 95% bootstrap lower bound, and a provisional reading:

- lower bound of at least 2: the direction holds for interval ratios below 2
- crossing ratio below 2 with at least 300 respondents per group: potentially sensitive to unequal intervals
- otherwise: undecided

Required packages: `ordinal`, `psychTools`, and `gssr` (for the General Social Survey data).

## License

MIT
