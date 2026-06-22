package ui
import (
	"fmt"
	"math"
	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/canvas"
	"fyne.io/fyne/v2/widget"
	"ballistx/internal/ballistics"
	"ballistx/internal/units"
)
type trajectoryChart struct {
	widget.BaseWidget
	points []ballistics.Point
	metric bool
}
func newTrajectoryChart() *trajectoryChart {
	c := &trajectoryChart{}
	c.ExtendBaseWidget(c)
	return c
}
func (c *trajectoryChart) SetData(points []ballistics.Point, metric bool) {
	c.points = points
	c.metric = metric
	c.Refresh()
}
func (c *trajectoryChart) CreateRenderer() fyne.WidgetRenderer {
	bg := canvas.NewRectangle(colInputBg)
	zeroLine := canvas.NewLine(colSeparator)
	return &chartRenderer{chart: c, bg: bg, zeroLine: zeroLine}
}
type chartRenderer struct {
	chart    *trajectoryChart
	bg       *canvas.Rectangle
	zeroLine *canvas.Line
	lines    []fyne.CanvasObject
	labels   []fyne.CanvasObject
}
func (r *chartRenderer) Layout(size fyne.Size) {
	r.bg.Resize(size)
	r.draw(size)
}
func (r *chartRenderer) MinSize() fyne.Size {
	return fyne.NewSize(300, 200)
}
func (r *chartRenderer) Refresh() {
	r.draw(r.chart.Size())
	canvas.Refresh(r.chart)
}
func (r *chartRenderer) Objects() []fyne.CanvasObject {
	objs := []fyne.CanvasObject{r.bg, r.zeroLine}
	objs = append(objs, r.lines...)
	objs = append(objs, r.labels...)
	return objs
}
func (r *chartRenderer) Destroy() {}
func (r *chartRenderer) draw(size fyne.Size) {
	r.lines = nil
	r.labels = nil
	r.bg.Resize(size)
	pts := r.chart.points
	if len(pts) < 2 || size.Width <= 0 || size.Height <= 0 {
		r.zeroLine.Hide()
		return
	}
	const padL, padR, padT, padB = 48, 16, 16, 28
	plotW := float32(size.Width) - padL - padR
	plotH := float32(size.Height) - padT - padB
	if plotW <= 0 || plotH <= 0 {
		return
	}
	minDrop, maxDrop := math.Inf(1), math.Inf(-1)
	maxRange := 0.0
	for _, p := range pts {
		d := p.DropInches
		if r.chart.metric {
			d = units.InchesToMM(d) / 10.0
		}
		if d < minDrop {
			minDrop = d
		}
		if d > maxDrop {
			maxDrop = d
		}
		rangeVal := p.RangeYards
		if r.chart.metric {
			rangeVal = units.YardsToMeters(rangeVal)
		}
		if rangeVal > maxRange {
			maxRange = rangeVal
		}
	}
	if maxDrop == minDrop {
		maxDrop += 1
		minDrop -= 1
	}
	if minDrop > 0 {
		minDrop = 0
	}
	if maxDrop < 0 {
		maxDrop = 0
	}
	xAt := func(rangeVal float64) float32 {
		return padL + float32(rangeVal/maxRange)*plotW
	}
	yAt := func(drop float64) float32 {
		frac := (drop - minDrop) / (maxDrop - minDrop)
		return padT + plotH - float32(frac)*plotH
	}
	zy := yAt(0)
	r.zeroLine.Position1 = fyne.NewPos(padL, zy)
	r.zeroLine.Position2 = fyne.NewPos(padL+plotW, zy)
	r.zeroLine.StrokeWidth = 1
	r.zeroLine.Show()
	dropUnit := "in"
	rangeUnit := "yd"
	if r.chart.metric {
		dropUnit = "cm"
		rangeUnit = "m"
	}
	topLabel := canvas.NewText(fmt.Sprintf("%.0f%s", maxDrop, dropUnit), colForeground)
	topLabel.TextSize = 11
	topLabel.Move(fyne.NewPos(2, padT-6))
	botLabel := canvas.NewText(fmt.Sprintf("%.0f%s", minDrop, dropUnit), colForeground)
	botLabel.TextSize = 11
	botLabel.Move(fyne.NewPos(2, padT+plotH-8))
	rangeLabel := canvas.NewText(fmt.Sprintf("0 — %.0f%s", maxRange, rangeUnit), colForeground)
	rangeLabel.TextSize = 11
	rangeLabel.Move(fyne.NewPos(padL, padT+plotH+6))
	r.labels = []fyne.CanvasObject{topLabel, botLabel, rangeLabel}
	prevX, prevY := xAt(0), yAt(0)
	for i, p := range pts {
		rangeVal := p.RangeYards
		drop := p.DropInches
		if r.chart.metric {
			rangeVal = units.YardsToMeters(rangeVal)
			drop = units.InchesToMM(drop) / 10.0
		}
		x, y := xAt(rangeVal), yAt(drop)
		if i > 0 {
			seg := canvas.NewLine(colPrimary)
			seg.StrokeWidth = 2.5
			seg.Position1 = fyne.NewPos(prevX, prevY)
			seg.Position2 = fyne.NewPos(x, y)
			r.lines = append(r.lines, seg)
		}
		prevX, prevY = x, y
	}
	prevX, prevY = xAt(0), yAt(0)
	for i, p := range pts {
		rangeVal := p.RangeYards
		drift := p.WindageInches
		if r.chart.metric {
			rangeVal = units.YardsToMeters(rangeVal)
			drift = units.InchesToMM(drift) / 10.0
		}
		x, y := xAt(rangeVal), yAt(drift)
		if i > 0 {
			seg := canvas.NewLine(colPlum)
			seg.StrokeWidth = 1.5
			seg.Position1 = fyne.NewPos(prevX, prevY)
			seg.Position2 = fyne.NewPos(x, y)
			r.lines = append(r.lines, seg)
		}
		prevX, prevY = x, y
	}
}
