package ui
import (
	"strconv"
	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/widget"
	"ballistx/internal/ballistics"
	"ballistx/internal/units"
)
const windDirHelp = "0°=12 o'clock (headwind) • 90°=3 o'clock • 180°=6 o'clock (tailwind) • 270°=9 o'clock"
type conditionsPanel struct {
	state *State
	win   fyne.Window
	presetSelect *widget.Select
	presetNotes  *widget.Label
	bcAdjust     *widget.Entry
	unitsRadio *widget.RadioGroup
	altitude    *widget.Entry
	temperature *widget.Entry
	humidity    *widget.Entry
	pressure    *widget.Entry
	useStation  *widget.Check
	windSpeed *widget.Entry
	windDir   *widget.Entry
	sightType   *widget.Select
	sightHeight *widget.Entry
	sightAngle  *widget.Entry
	zeroRange   *widget.Entry
	barrelLength *widget.Entry
	twist        *widget.Entry
	twistDir     *widget.RadioGroup
	maxRange  *widget.Entry
	rangeStep *widget.Entry
	content fyne.CanvasObject
	unitForms                                     []*widget.Form
	altitudeItem, temperatureItem, pressureItem   *widget.FormItem
	windSpeedItem, sightHeightItem, zeroRangeItem *widget.FormItem
	barrelLengthItem, maxRangeItem, rangeStepItem *widget.FormItem
	onCalculate func()
}
func newConditionsPanel(win fyne.Window, state *State, onCalculate func()) *conditionsPanel {
	p := &conditionsPanel{state: state, win: win, onCalculate: onCalculate}
	p.build()
	return p
}
func sectionCard(title string, content fyne.CanvasObject) fyne.CanvasObject {
	return widget.NewCard(title, "", content)
}
func numEntry(value float64) *widget.Entry {
	e := widget.NewEntry()
	e.SetText(formatNum(value))
	return e
}
func formatNum(v float64) string {
	rounded := strconv.FormatFloat(v, 'f', 4, 64)
	v2, _ := strconv.ParseFloat(rounded, 64)
	return strconv.FormatFloat(v2, 'f', -1, 64)
}
func parseNum(e *widget.Entry, fallback float64) float64 {
	v, err := strconv.ParseFloat(e.Text, 64)
	if err != nil {
		return fallback
	}
	return v
}
func (p *conditionsPanel) build() {
	s := p.state
	names := make([]string, len(s.Catalog.Presets))
	for i, pr := range s.Catalog.Presets {
		names[i] = pr.DisplayName()
	}
	p.presetSelect = widget.NewSelect(names, func(selected string) {
		for _, pr := range s.Catalog.Presets {
			if pr.DisplayName() == selected {
				s.Preset = pr
				p.bcAdjust.SetText(formatNum(pr.BCAdjustmentPct))
				p.presetNotes.SetText(pr.Notes)
				if s.BarrelLengthIn == 0 {
				p.barrelLength.SetText(formatNum(pr.RefBarrelLengthIn))
				}
				break
			}
		}
	})
	p.presetNotes = widget.NewLabel(s.Preset.Notes)
	p.presetNotes.Wrapping = fyne.TextWrapWord
	p.bcAdjust = numEntry(s.Preset.BCAdjustmentPct)
	manageBtn := widget.NewButton("Manage Ammo Presets…", func() {
		showPresetEditor(p.win, s.Catalog, func() {
			names := make([]string, len(s.Catalog.Presets))
			for i, pr := range s.Catalog.Presets {
				names[i] = pr.DisplayName()
			}
			p.presetSelect.Options = names
			p.presetSelect.Refresh()
		})
	})
	ammoBox := container.NewVBox(
		widget.NewForm(
			widget.NewFormItem("Ammunition", p.presetSelect),
			widget.NewFormItem("BC Adjustment (%)", p.bcAdjust),
		),
		p.presetNotes,
		manageBtn,
	)
	p.unitsRadio = widget.NewRadioGroup([]string{"Imperial", "Metric"}, nil)
	p.unitsRadio.Horizontal = true
	p.unitsRadio.SetSelected("Imperial")
	p.altitude = numEntry(0)
	p.temperature = numEntry(59)
	p.humidity = numEntry(50)
	p.pressure = numEntry(29.92)
	p.useStation = widget.NewCheck("Use entered station pressure (instead of estimating from altitude)", nil)
	p.altitudeItem = widget.NewFormItem("Altitude (ft)", p.altitude)
	p.temperatureItem = widget.NewFormItem("Temperature (°F)", p.temperature)
	p.pressureItem = widget.NewFormItem("Barometric Pressure (inHg)", p.pressure)
	atmForm := widget.NewForm(
		p.altitudeItem,
		p.temperatureItem,
		widget.NewFormItem("Humidity (%)", p.humidity),
		p.pressureItem,
	)
	atmBox := container.NewVBox(atmForm, p.useStation)
	p.windSpeed = numEntry(0)
	p.windDir = numEntry(0)
	windHelp := widget.NewLabel(windDirHelp)
	windHelp.Wrapping = fyne.TextWrapWord
	p.windSpeedItem = widget.NewFormItem("Wind Speed (mph)", p.windSpeed)
	windForm := widget.NewForm(
		p.windSpeedItem,
		widget.NewFormItem("Wind Direction (°)", p.windDir),
	)
	windBox := container.NewVBox(windForm, windHelp)
	p.sightHeight = numEntry(1.5)
	p.sightType = widget.NewSelect([]string{"Scope", "Iron Sights", "Red Dot", "Holographic"}, func(t string) {
		switch t {
		case "Iron Sights":
			p.sightHeight.SetText(formatNum(p.fromInches(0.6)))
		case "Red Dot", "Holographic":
			p.sightHeight.SetText(formatNum(p.fromInches(2.0)))
		default:
			p.sightHeight.SetText(formatNum(p.fromInches(1.5)))
		}
	})
	p.sightType.SetSelected("Scope")
	p.sightAngle = numEntry(0)
	p.zeroRange = numEntry(100)
	p.sightHeightItem = widget.NewFormItem("Sight Height over Bore (in)", p.sightHeight)
	p.zeroRangeItem = widget.NewFormItem("Zero Range (yd)", p.zeroRange)
	sightForm := widget.NewForm(
		widget.NewFormItem("Sight Type", p.sightType),
		p.sightHeightItem,
		widget.NewFormItem("Sight Cant / Angle (°)", p.sightAngle),
		p.zeroRangeItem,
	)
	p.barrelLength = numEntry(s.Preset.RefBarrelLengthIn)
	p.twist = numEntry(12)
	p.twistDir = widget.NewRadioGroup([]string{"Right-Hand", "Left-Hand"}, nil)
	p.twistDir.Horizontal = true
	p.twistDir.SetSelected("Right-Hand")
	p.barrelLengthItem = widget.NewFormItem("Barrel Length (in)", p.barrelLength)
	barrelForm := widget.NewForm(
		p.barrelLengthItem,
		widget.NewFormItem("Twist Rate (in/turn)", p.twist),
	)
	barrelBox := container.NewVBox(barrelForm, p.twistDir)
	p.maxRange = numEntry(800)
	p.rangeStep = numEntry(50)
	p.maxRangeItem = widget.NewFormItem("Max Range (yd)", p.maxRange)
	p.rangeStepItem = widget.NewFormItem("Table Step (yd)", p.rangeStep)
	rangeForm := widget.NewForm(
		p.maxRangeItem,
		p.rangeStepItem,
	)
	p.unitForms = []*widget.Form{atmForm, windForm, sightForm, barrelForm, rangeForm}
	p.unitsRadio.OnChanged = func(string) { p.applyUnitSystem() }
	calcBtn := widget.NewButton("Calculate Trajectory", func() {
		p.readInto(s)
		p.onCalculate()
	})
	calcBtn.Importance = widget.HighImportance
	if len(names) > 0 {
		p.presetSelect.SetSelected(s.Preset.DisplayName())
	}
	p.content = container.NewVBox(
		sectionCard("Ammunition & Material", ammoBox),
		sectionCard("Units", p.unitsRadio),
		sectionCard("Atmosphere", atmBox),
		sectionCard("Wind", windBox),
		sectionCard("Sight", sightForm),
		sectionCard("Barrel", barrelBox),
		sectionCard("Trajectory Range", rangeForm),
		calcBtn,
	)
}
func (p *conditionsPanel) readInto(s *State) {
	metric := p.unitsRadio.Selected == "Metric"
	s.UnitSys = units.Imperial
	if metric {
		s.UnitSys = units.Metric
	}
	if metric {
		s.AltitudeFt = units.MetersToFeet(parseNum(p.altitude, 0))
		s.TemperatureF = units.CelsiusToFahrenheit(parseNum(p.temperature, 15))
		s.PressureInHg = units.HPaToInHg(parseNum(p.pressure, 1013.25))
		s.WindSpeedFPS = units.MPSToFPS(parseNum(p.windSpeed, 0))
		s.SightHeightIn = units.CentimetersToInches(parseNum(p.sightHeight, 3.8))
		s.ZeroRangeYards = units.MetersToYards(parseNum(p.zeroRange, 100))
		s.BarrelLengthIn = units.CentimetersToInches(parseNum(p.barrelLength, 60))
		s.MaxRangeYards = units.MetersToYards(parseNum(p.maxRange, 800))
		s.RangeStepYards = units.MetersToYards(parseNum(p.rangeStep, 50))
	} else {
		s.AltitudeFt = parseNum(p.altitude, 0)
		s.TemperatureF = parseNum(p.temperature, 59)
		s.PressureInHg = parseNum(p.pressure, 29.92)
		s.WindSpeedFPS = units.MPHToFPS(parseNum(p.windSpeed, 0))
		s.SightHeightIn = parseNum(p.sightHeight, 1.5)
		s.ZeroRangeYards = parseNum(p.zeroRange, 100)
		s.BarrelLengthIn = parseNum(p.barrelLength, 24)
		s.MaxRangeYards = parseNum(p.maxRange, 800)
		s.RangeStepYards = parseNum(p.rangeStep, 50)
	}
	s.HumidityPct = parseNum(p.humidity, 50)
	s.UseStationPressure = p.useStation.Checked
	s.WindDirectionDeg = parseNum(p.windDir, 0)
	s.SightAngleDeg = parseNum(p.sightAngle, 0)
	s.SightType = p.sightType.Selected
	s.TwistIn = parseNum(p.twist, 12)
	s.TwistDirection = ballistics.RightHand
	if p.twistDir.Selected == "Left-Hand" {
		s.TwistDirection = ballistics.LeftHand
	}
	s.BCAdjustmentPct = parseNum(p.bcAdjust, 0)
}
func (p *conditionsPanel) fromInches(in float64) float64 {
	if p.unitsRadio.Selected == "Metric" {
		return units.InchesToCentimeters(in)
	}
	return in
}
func (p *conditionsPanel) applyUnitSystem() {
	metric := p.unitsRadio.Selected == "Metric"

	convert := func(e *widget.Entry, toMetric func(float64) float64, toImperial func(float64) float64, fallback float64) {
		v := parseNum(e, fallback)
		if metric {
			e.SetText(formatNum(toMetric(v)))
		} else {
			e.SetText(formatNum(toImperial(v)))
		}
	}
	if metric {
		p.altitudeItem.Text = "Altitude (m)"
		p.temperatureItem.Text = "Temperature (°C)"
		p.pressureItem.Text = "Barometric Pressure (hPa)"
		p.windSpeedItem.Text = "Wind Speed (m/s)"
		p.sightHeightItem.Text = "Sight Height over Bore (cm)"
		p.zeroRangeItem.Text = "Zero Range (m)"
		p.barrelLengthItem.Text = "Barrel Length (cm)"
		p.maxRangeItem.Text = "Max Range (m)"
		p.rangeStepItem.Text = "Table Step (m)"
	} else {
		p.altitudeItem.Text = "Altitude (ft)"
		p.temperatureItem.Text = "Temperature (°F)"
		p.pressureItem.Text = "Barometric Pressure (inHg)"
		p.windSpeedItem.Text = "Wind Speed (mph)"
		p.sightHeightItem.Text = "Sight Height over Bore (in)"
		p.zeroRangeItem.Text = "Zero Range (yd)"
		p.barrelLengthItem.Text = "Barrel Length (in)"
		p.maxRangeItem.Text = "Max Range (yd)"
		p.rangeStepItem.Text = "Table Step (yd)"
	}
	convert(p.altitude, units.FeetToMeters, units.MetersToFeet, 0)
	convert(p.temperature, units.FahrenheitToCelsius, units.CelsiusToFahrenheit, 59)
	convert(p.pressure, units.InHgToHPa, units.HPaToInHg, 29.92)
	convert(p.windSpeed, units.MPHToMPS, units.MPSToMPH, 0)
	convert(p.sightHeight, units.InchesToCentimeters, units.CentimetersToInches, 1.5)
	convert(p.zeroRange, units.YardsToMeters, units.MetersToYards, 100)
	convert(p.barrelLength, units.InchesToCentimeters, units.CentimetersToInches, 24)
	convert(p.maxRange, units.YardsToMeters, units.MetersToYards, 800)
	convert(p.rangeStep, units.YardsToMeters, units.MetersToYards, 50)
	for _, f := range p.unitForms {
		f.Refresh()
	}
}
