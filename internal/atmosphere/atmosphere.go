package atmosphere
import "math"
const (
	StdTemperatureF   = 59.0
	StdPressureInHg   = 29.92
	StdHumidityPct    = 0.0
	StdTempRankine    = StdTemperatureF + 459.67
	dryAirGasConstant = 287.058
	vaporGasConstant  = 461.495
)
type Conditions struct {
	AltitudeFt         float64
	TemperatureF       float64
	HumidityPct        float64
	PressureInHg       float64
	UseStationPressure bool
}
func standardPressureAtAltitude(altitudeFt float64) float64 {
	return StdPressureInHg * math.Pow(1.0-6.8755856e-6*altitudeFt, 5.2558797)
}
func saturationVaporPressureInHg(tempF float64) float64 {
	tc := (tempF - 32.0) * 5.0 / 9.0
	hPa := 6.1121 * math.Exp((18.678-tc/234.5)*(tc/(257.14+tc)))
	return hPa * 0.0295299830714
}
func (c Conditions) DensityRatio() float64 {
	pressure := c.PressureInHg
	if !c.UseStationPressure || pressure <= 0 {
		pressure = standardPressureAtAltitude(c.AltitudeFt)
	}
	tempR := c.TemperatureF + 459.67
	if tempR <= 0 {
		tempR = StdTempRankine
	}
	vaporPressure := saturationVaporPressureInHg(c.TemperatureF) * (c.HumidityPct / 100.0)
	humidityFactor := 1.0 - dryAirGasConstant/vaporGasConstant
	effectivePressure := pressure - humidityFactor*vaporPressure
	return (effectivePressure / StdPressureInHg) * (StdTempRankine / tempR)
}
func SpeedOfSoundFPS(tempF float64) float64 {
	tempR := tempF + 459.67
	if tempR <= 0 {
		tempR = StdTempRankine
	}
	return 49.0223 * math.Sqrt(tempR)
}
