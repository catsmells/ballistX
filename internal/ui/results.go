package ui
import (
	"fmt"
	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/widget"
	"ballistx/internal/ballistics"
	"ballistx/internal/units"
)
var tableHeaders = []string{"Range", "Time (s)", "Velocity", "Energy", "Drop", "Drop MOA", "Drop MIL", "Wind Drift", "Wind MOA", "Wind MIL"}
type resultsView struct {
	state *State
	table *widget.Table
	chart *trajectoryChart
	note  *widget.Label
	rows [][]string
}
func newResultsView(state *State) *resultsView {
	r := &resultsView{state: state, rows: [][]string{tableHeaders}}
	r.table = widget.NewTable(
		func() (int, int) { return len(r.rows), len(tableHeaders) },
		func() fyne.CanvasObject {
			l := widget.NewLabel("")
			l.Alignment = fyne.TextAlignCenter
			return l
		},
		func(id widget.TableCellID, obj fyne.CanvasObject) {
			label := obj.(*widget.Label)
			label.SetText(r.rows[id.Row][id.Col])
			if id.Row == 0 {
				label.TextStyle = fyne.TextStyle{Bold: true}
			} else {
				label.TextStyle = fyne.TextStyle{}
			}
		},
	)
	for i := range tableHeaders {
		r.table.SetColumnWidth(i, 96)
	}
	r.table.SetRowHeight(0, 32)
	r.chart = newTrajectoryChart()
	r.note = widget.NewLabel("")
	r.note.Wrapping = fyne.TextWrapWord
	return r
}
func (r *resultsView) canvasObject() fyne.CanvasObject {
	tableCard := widget.NewCard("Trajectory Table", "", container.NewStack(r.table))
	chartCard := widget.NewCard("Drop & Wind Drift", "", r.chart)
	split := container.NewVSplit(chartCard, tableCard)
	split.Offset = 0.4
	return container.NewBorder(r.note, nil, nil, nil, split)
}
func (r *resultsView) Update(points []ballistics.Point, cappedAt float64) {
	metric := r.state.UnitSys == units.Metric
	rows := [][]string{tableHeaders}
	if metric {
		rows[0] = []string{"Range (m)", "Time (s)", "Velocity (m/s)", "Energy (J)", "Drop (cm)", "Drop MOA", "Drop MIL", "Drift (cm)", "Wind MOA", "Wind MIL"}
	} else {
		rows[0] = []string{"Range (yd)", "Time (s)", "Velocity (fps)", "Energy (ft·lb)", "Drop (in)", "Drop MOA", "Drop MIL", "Drift (in)", "Wind MOA", "Wind MIL"}
	}
	for _, pt := range points {
		var rangeVal, vel, energy, drop, drift float64
		if metric {
			rangeVal = units.YardsToMeters(pt.RangeYards)
			vel = units.FPSToMPS(pt.VelocityFPS)
			energy = pt.EnergyFtLbs * 1.355818
			drop = units.InchesToMM(pt.DropInches) / 10.0
			drift = units.InchesToMM(pt.WindageInches) / 10.0
		} else {
			rangeVal = pt.RangeYards
			vel = pt.VelocityFPS
			energy = pt.EnergyFtLbs
			drop = pt.DropInches
			drift = pt.WindageInches
		}
		rows = append(rows, []string{
			fmt.Sprintf("%.0f", rangeVal),
			fmt.Sprintf("%.3f", pt.TimeOfFlightSec),
			fmt.Sprintf("%.0f", vel),
			fmt.Sprintf("%.0f", energy),
			fmt.Sprintf("%.1f", drop),
			fmt.Sprintf("%.2f", pt.DropMOA),
			fmt.Sprintf("%.2f", pt.DropMIL),
			fmt.Sprintf("%.1f", drift),
			fmt.Sprintf("%.2f", pt.WindageMOA),
			fmt.Sprintf("%.2f", pt.WindageMIL),
		})
	}
	r.rows = rows
	r.table.Refresh()
	r.chart.SetData(points, metric)
	if cappedAt > 0 {
		unit := "yd"
		val := cappedAt
		if metric {
			unit = "m"
			val = units.YardsToMeters(cappedAt)
		}
		r.note.SetText(fmt.Sprintf("Note: trajectory capped at %.0f%s, the reasonable maximum range for this ammunition preset.", val, unit))
	} else {
		r.note.SetText("")
	}
}
