figure-cycle-thresholds.pdf: figure-cycle-thresholds.R cycle.thresholds.RData
	R --no-save < $<
cycle.thresholds.RData: cycle.thresholds.R
	R --no-save < $<

