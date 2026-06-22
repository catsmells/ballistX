package ballistics
import (
	"testing"
	"ballistx/internal/atmosphere"
	"ballistx/internal/drag"
)
func TestSanity308Win(t *testing.T) {
	atm := atmosphere.Conditions{TemperatureF: 59, PressureInHg: 29.92, HumidityPct: 0}
	in := Inputs{
		BC:                0.462,
		DragFunction:      drag.G1,
		MuzzleVelocityFPS: 2650,
		WeightGrains:      168,
		SightHeightIn:     1.5,
		ZeroRangeYards:    100,
		DensityRatio:      atm.DensityRatio(),
		SpeedOfSoundFPS:   atmosphere.SpeedOfSoundFPS(59),
		MaxRangeYards:     1200,
	}
	pts, _ := in.Solve(1000, 100)
	got := map[int]Point{}
	for _, p := range pts {
		got[int(p.RangeYards+0.5)] = p
	}
	p500, ok := got[500]
	if !ok {
		t.Fatalf("missing 500yd sample; points=%d", len(pts))
	}
	if p500.DropInches > -40 || p500.DropInches < -75 {
		t.Errorf("500yd drop out of expected band: got %.1fin, want roughly -50 to -60in", p500.DropInches)
	}
	p1000, ok := got[1000]
	if !ok {
		t.Fatalf("missing 1000yd sample")
	}
	if p1000.DropInches > -330 || p1000.DropInches < -480 {
		t.Errorf("1000yd drop out of expected band: got %.1fin, want roughly -400 to -430in", p1000.DropInches)
	}
	t.Logf("100yd v=%.0ffps, 500yd drop=%.1fin v=%.0ffps, 1000yd drop=%.1fin v=%.0ffps",
		got[100].VelocityFPS, p500.DropInches, p500.VelocityFPS, p1000.DropInches, p1000.VelocityFPS)
}
