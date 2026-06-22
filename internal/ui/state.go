package ui
import (
	"ballistx/internal/ammo"
	"ballistx/internal/atmosphere"
	"ballistx/internal/ballistics"
	"ballistx/internal/drag"
	"ballistx/internal/units"
)
type State struct {
	Catalog *ammo.Catalog
	Preset  ammo.Preset
	UnitSys units.System
	AltitudeFt         float64
	TemperatureF       float64
	HumidityPct        float64
	PressureInHg       float64
	UseStationPressure bool
	WindSpeedFPS     float64
	WindDirectionDeg float64
	SightHeightIn  float64
	SightAngleDeg  float64
	SightType      string
	ZeroRangeYards float64
	BarrelLengthIn float64
	TwistIn        float64
	TwistDirection ballistics.TwistDirection
	BCAdjustmentPct float64
	MaxRangeYards  float64
	RangeStepYards float64
}
func DefaultState(cat *ammo.Catalog) *State {
	preset := ammo.Preset{}
	if len(cat.Presets) > 0 {
		preset = cat.Presets[0]
	}
	return &State{
		Catalog:            cat,
		Preset:             preset,
		UnitSys:            units.Imperial,
		AltitudeFt:         0,
		TemperatureF:       59,
		HumidityPct:        50,
		PressureInHg:       29.92,
		UseStationPressure: false,
		WindSpeedFPS:       0,
		WindDirectionDeg:   0,
		SightHeightIn:      1.5,
		SightAngleDeg:      0,
		SightType:          "Scope",
		ZeroRangeYards:     100,
		BarrelLengthIn:     preset.RefBarrelLengthIn,
		TwistIn:            12,
		TwistDirection:     ballistics.RightHand,
		BCAdjustmentPct:    preset.BCAdjustmentPct,
		MaxRangeYards:      800,
		RangeStepYards:     50,
	}
}
func (s *State) Solve() ([]ballistics.Point, float64) {
	atm := atmosphere.Conditions{
		AltitudeFt:         s.AltitudeFt,
		TemperatureF:       s.TemperatureF,
		HumidityPct:        s.HumidityPct,
		PressureInHg:       s.PressureInHg,
		UseStationPressure: s.UseStationPressure,
	}
	densityRatio := atm.DensityRatio()
	sos := atmosphere.SpeedOfSoundFPS(s.TemperatureF)
	barrel := ballistics.BarrelSpec{
		LengthIn:           s.BarrelLengthIn,
		RefLengthIn:        s.Preset.RefBarrelLengthIn,
		VelocityPerInchFPS: s.Preset.VelocityPerInchFPS,
		TwistIn:            s.TwistIn,
		TwistDirection:     s.TwistDirection,
		BulletDiameterIn:   s.Preset.BulletDiameterIn,
		BulletLengthIn:     s.Preset.BulletLengthIn,
	}
	muzzleVelocity := barrel.EstimatedMuzzleVelocity(s.Preset.MuzzleVelocityFPS)
	headwind, crosswind := units.WindComponents(s.WindSpeedFPS, s.WindDirectionDeg)
	bc := s.Preset.BC * (1.0 + s.BCAdjustmentPct/100.0)
	in := ballistics.Inputs{
		BC:                bc,
		DragFunction:      drag.Function(s.Preset.DragFunction),
		MuzzleVelocityFPS: muzzleVelocity,
		WeightGrains:      s.Preset.BulletWeightGrains,
		SightHeightIn:     s.SightHeightIn,
		SightAngleDeg:     s.SightAngleDeg,
		ZeroRangeYards:    s.ZeroRangeYards,
		DensityRatio:      densityRatio,
		SpeedOfSoundFPS:   sos,
		HeadwindFPS:       headwind,
		CrosswindFPS:      crosswind,
		Barrel:            barrel,
		MaxRangeYards:     s.Preset.MaxRangeYards,
	}
	return in.Solve(s.MaxRangeYards, s.RangeStepYards)
}
