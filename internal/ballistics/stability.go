package ballistics
import "math"
type TwistDirection int
const (
	RightHand TwistDirection = iota
	LeftHand
)
type BarrelSpec struct {
	LengthIn           float64
	RefLengthIn        float64
	VelocityPerInchFPS float64
	TwistIn            float64
	TwistDirection     TwistDirection
	BulletDiameterIn   float64
	BulletLengthIn     float64
}
func (b BarrelSpec) EstimatedMuzzleVelocity(refVelocityFPS float64) float64 {
	if b.RefLengthIn <= 0 || b.LengthIn <= 0 {
		return refVelocityFPS
	}
	delta := b.LengthIn - b.RefLengthIn
	return refVelocityFPS + delta*b.VelocityPerInchFPS
}
func (b BarrelSpec) MillerStability(weightGrains, muzzleVelocityFPS, densityRatio float64) float64 {
	if b.TwistIn <= 0 || b.BulletDiameterIn <= 0 || b.BulletLengthIn <= 0 || weightGrains <= 0 {
		return 0
	}
	t := b.TwistIn / b.BulletDiameterIn
	l := b.BulletLengthIn / b.BulletDiameterIn
	sg := millerStabilityRaw(weightGrains, t, b.BulletDiameterIn, l)
	return sg * stabilityVelocityFactor(muzzleVelocityFPS) * stabilityDensityFactor(densityRatio)
}
func millerStabilityRaw(weightGrains, twistCalibers, diameterIn, lengthCalibers float64) float64 {
	denom := twistCalibers * twistCalibers * diameterIn * diameterIn * diameterIn * lengthCalibers * (1 + lengthCalibers*lengthCalibers)
	if denom == 0 {
		return 0
	}
	return (30 * weightGrains) / denom
}
func stabilityVelocityFactor(velocityFPS float64) float64 {
	if velocityFPS <= 0 {
		return 1
	}
	return math.Pow(velocityFPS/2800.0, 1.0/3.0)
}
func stabilityDensityFactor(densityRatio float64) float64 {
	if densityRatio <= 0 {
		return 1
	}
	return 1.0 / densityRatio
}
func SpinDriftInches(stabilityFactor, timeOfFlightSec float64, dir TwistDirection) float64 {
	if stabilityFactor <= 0 || timeOfFlightSec <= 0 {
		return 0
	}
	drift := 1.25 * (stabilityFactor + 1.2) * math.Pow(timeOfFlightSec, 1.83)
	if dir == LeftHand {
		return -drift
	}
	return drift
}
