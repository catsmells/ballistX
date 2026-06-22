package units
import "math"
type System int
const (
	Imperial System = iota
	Metric
)
func (s System) String() string {
	if s == Metric {
		return "Metric"
	}
	return "Imperial"
}
const metersPerFoot = 0.3048
const metersPerYard = 0.9144
func MetersToFeet(m float64) float64 { return m / metersPerFoot }
func FeetToMeters(f float64) float64 { return f * metersPerFoot }
func MetersToYards(m float64) float64 { return m / metersPerYard }
func YardsToMeters(y float64) float64 { return y * metersPerYard }
func YardsToFeet(y float64) float64 { return y * 3.0 }
func FeetToYards(f float64) float64 { return f / 3.0 }
func InchesToFeet(i float64) float64 { return i / 12.0 }
func FeetToInches(f float64) float64 { return f * 12.0 }
func InchesToMeters(i float64) float64 { return i * 0.0254 }
func MetersToInches(m float64) float64 { return m / 0.0254 }
func InchesToMM(i float64) float64  { return i * 25.4 }
func MMToInches(mm float64) float64 { return mm / 25.4 }
func InchesToCentimeters(i float64) float64  { return i * 2.54 }
func CentimetersToInches(cm float64) float64 { return cm / 2.54 }
func MPSToFPS(mps float64) float64 { return mps * 3.280839895 }
func FPSToMPS(fps float64) float64 { return fps / 3.280839895 }
func MPHToFPS(mph float64) float64 { return mph * 1.46666667 }
func FPSToMPH(fps float64) float64 { return fps / 1.46666667 }
func MPHToMPS(mph float64) float64 { return mph * 0.44704 }
func MPSToMPH(mps float64) float64 { return mps / 0.44704 }
const grainsPerGram = 15.4323584
const grainsPerPound = 7000.0
func GrainsToGrams(gr float64) float64 { return gr / grainsPerGram }
func GramsToGrains(g float64) float64  { return g * grainsPerGram }
func GrainsToPounds(gr float64) float64 { return gr / grainsPerPound }
func PoundsToGrains(lb float64) float64 { return lb * grainsPerPound }
const inHgPerHPa = 0.0295299830714
func HPaToInHg(hpa float64) float64  { return hpa * inHgPerHPa }
func InHgToHPa(inHg float64) float64 { return inHg / inHgPerHPa }
func CelsiusToFahrenheit(c float64) float64 { return c*9.0/5.0 + 32.0 }
func FahrenheitToCelsius(f float64) float64 { return (f - 32.0) * 5.0 / 9.0 }
func FahrenheitToRankine(f float64) float64 { return f + 459.67 }
func CelsiusToKelvin(c float64) float64     { return c + 273.15 }
const (
	InchesPerMOAAt100Yards = 1.0471975512
	InchesPerMILAt100Yards = 3.6
)
func InchesToMOA(inches, rangeYards float64) float64 {
	if rangeYards <= 0 {
		return 0
	}
	return inches / (InchesPerMOAAt100Yards * rangeYards / 100.0)
}
func InchesToMIL(inches, rangeYards float64) float64 {
	if rangeYards <= 0 {
		return 0
	}
	return inches / (InchesPerMILAt100Yards * rangeYards / 100.0)
}
func CentimetersToMOAAtMeters(cm, rangeMeters float64) float64 {
	if rangeMeters <= 0 {
		return 0
	}
	inchesPerMOAAtRange := InchesPerMOAAt100Yards * (rangeMeters * 1.0936132983) / 100.0
	cmPerMOAAtRange := inchesPerMOAAtRange * 2.54
	return cm / cmPerMOAAtRange
}
func CentimetersToMILAtMeters(cm, rangeMeters float64) float64 {
	if rangeMeters <= 0 {
		return 0
	}
	cmPerMILAtRange := 10.0 * (rangeMeters / 100.0)
	return cm / cmPerMILAtRange
}
func WindComponents(speed, directionDeg float64) (headwind, crosswind float64) {
	rad := directionDeg * (math.Pi / 180.0)
	headwind = speed * math.Cos(rad)
	crosswind = speed * math.Sin(rad)
	return
}
