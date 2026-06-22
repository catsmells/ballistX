package ballistics
import (
	"math"
	"ballistx/internal/drag"
)
const dragConstant = 0.00020863
const gravityFPS2 = 32.174
type Inputs struct {
	BC                float64
	DragFunction      drag.Function
	MuzzleVelocityFPS float64
	WeightGrains      float64
	SightHeightIn  float64
	SightAngleDeg  float64
	ZeroRangeYards float64
	DensityRatio    float64
	SpeedOfSoundFPS float64
	HeadwindFPS  float64
	CrosswindFPS float64
	Barrel        BarrelSpec
	MaxRangeYards float64
}
type Point struct {
	RangeYards      float64
	TimeOfFlightSec float64
	VelocityFPS     float64
	Mach            float64
	DropInches      float64
	WindageInches   float64
	DropMOA         float64
	DropMIL         float64
	WindageMOA      float64
	WindageMIL      float64
	EnergyFtLbs     float64
}
type state struct {
	x, y, z    float64
	vx, vy, vz float64
	t          float64
}
func (in Inputs) step(s *state, dt float64) {
	relVX := s.vx + in.HeadwindFPS
	relVY := s.vy
	relVZ := s.vz - in.CrosswindFPS
	relSpeed := math.Sqrt(relVX*relVX + relVY*relVY + relVZ*relVZ)
	if relSpeed <= 0 {
		return
	}
	mach := relSpeed / in.SpeedOfSoundFPS
	cd := drag.Cd(in.DragFunction, mach)
	k := dragConstant * in.DensityRatio * cd / in.BC
	decel := k * relSpeed
	s.vx -= decel * relVX * dt
	s.vy -= decel*relVY*dt + gravityFPS2*dt
	s.vz -= decel * relVZ * dt
	s.x += s.vx * dt
	s.y += s.vy * dt
	s.z += s.vz * dt
	s.t += dt
}
func (in Inputs) simulate(elevationRad float64, maxRangeFt, stepFt float64) []Point {
	s := &state{
		y:  -in.SightHeightIn / 12.0,
		vx: in.MuzzleVelocityFPS * math.Cos(elevationRad),
		vy: in.MuzzleVelocityFPS * math.Sin(elevationRad),
	}
	dt := 0.0005
	var points []Point
	nextSample := 0.0
	stab := in.Barrel.MillerStability(in.WeightGrains, in.MuzzleVelocityFPS, in.DensityRatio)
	for steps := 0; steps < 2_000_000; steps++ {
		if s.x >= nextSample {
			rangeYards := s.x / 3.0
			velocity := math.Sqrt(s.vx*s.vx + s.vy*s.vy + s.vz*s.vz)
			mach := velocity / in.SpeedOfSoundFPS
			windageIn := s.z * 12.0
			if in.Barrel.TwistIn > 0 {
				windageIn += SpinDriftInches(stab, s.t, in.Barrel.TwistDirection)
			}
			dropIn := s.y * 12.0
			cantRad := in.SightAngleDeg * math.Pi / 180.0
			adjDrop := dropIn*math.Cos(cantRad) - windageIn*math.Sin(cantRad)
			adjWindage := dropIn*math.Sin(cantRad) + windageIn*math.Cos(cantRad)
			energy := in.WeightGrains * velocity * velocity / 450400.0
			points = append(points, Point{
				RangeYards:      rangeYards,
				TimeOfFlightSec: s.t,
				VelocityFPS:     velocity,
				Mach:            mach,
				DropInches:      adjDrop,
				WindageInches:   adjWindage,
				DropMOA:         -inchesToMOA(adjDrop, rangeYards),
				DropMIL:         -inchesToMIL(adjDrop, rangeYards),
				WindageMOA:      -inchesToMOA(adjWindage, rangeYards),
				WindageMIL:      -inchesToMIL(adjWindage, rangeYards),
				EnergyFtLbs:     energy,
			})
			nextSample += stepFt
		}
		if s.x >= maxRangeFt || s.vx <= 10 {
			break
		}
		in.step(s, dt)
	}
	return points
}
func inchesToMOA(inches, rangeYards float64) float64 {
	if rangeYards <= 0 {
		return 0
	}
	return inches / (1.0471975512 * rangeYards / 100.0)
}
func inchesToMIL(inches, rangeYards float64) float64 {
	if rangeYards <= 0 {
		return 0
	}
	return inches / (3.6 * rangeYards / 100.0)
}
func (in Inputs) heightAtRange(elevationRad, rangeFt float64) float64 {
	s := &state{
		y:  -in.SightHeightIn / 12.0,
		vx: in.MuzzleVelocityFPS * math.Cos(elevationRad),
		vy: in.MuzzleVelocityFPS * math.Sin(elevationRad),
	}
	dt := 0.0005
	for steps := 0; steps < 2_000_000; steps++ {
		if s.x >= rangeFt || s.vx <= 10 {
			break
		}
		in.step(s, dt)
	}
	return s.y * 12.0
}
func (in Inputs) SolveZeroAngle() float64 {
	zeroFt := in.ZeroRangeYards * 3.0
	if zeroFt <= 0 {
		return 0
	}
	lo, hi := -0.05, 0.15
	for i := 0; i < 60; i++ {
		mid := (lo + hi) / 2
		h := in.heightAtRange(mid, zeroFt)
		if h < 0 {
			lo = mid
		} else {
			hi = mid
		}
	}
	return (lo + hi) / 2
}
func (in Inputs) Solve(requestedMaxRangeYards, stepYards float64) (points []Point, cappedAt float64) {
	maxRange := requestedMaxRangeYards
	capped := 0.0
	if in.MaxRangeYards > 0 && in.MaxRangeYards < maxRange {
		capped = in.MaxRangeYards
		maxRange = in.MaxRangeYards
	}
	elevation := in.SolveZeroAngle()
	pts := in.simulate(elevation, maxRange*3.0, stepYards*3.0)
	return pts, capped
}
