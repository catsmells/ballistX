package drag
type Function string
const (
	G1 Function = "G1"
	G7 Function = "G7"
)
type point struct {
	mach float64
	cd   float64
}
func Cd(fn Function, mach float64) float64 {
	tbl := g1Table
	if fn == G7 {
		tbl = g7Table
	}
	if mach <= tbl[0].mach {
		return tbl[0].cd
	}
	last := tbl[len(tbl)-1]
	if mach >= last.mach {
		return last.cd
	}
	lo, hi := 0, len(tbl)-1
	for lo+1 < hi {
		mid := (lo + hi) / 2
		if tbl[mid].mach <= mach {
			lo = mid
		} else {
			hi = mid
		}
	}
	a, b := tbl[lo], tbl[hi]
	if b.mach == a.mach {
		return a.cd
	}
	frac := (mach - a.mach) / (b.mach - a.mach)
	return a.cd + frac*(b.cd-a.cd)
}
