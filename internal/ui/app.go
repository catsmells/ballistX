package ui
import (
	"log"
	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/app"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/widget"
	"ballistx/internal/ammo"
)
func Run() {
	cat, err := ammo.LoadOrCreateDefault(ammo.DefaultPath())
	if err != nil {
		log.Printf("ballistx: failed to load ammo catalog, using built-in defaults: %v", err)
		cat = &ammo.Catalog{Path: ammo.DefaultPath(), Presets: ammo.DefaultPresets()}
	}
	a := app.New()
	a.Settings().SetTheme(&ballistXTheme{})
	win := a.NewWindow("ballistX")
	state := DefaultState(cat)
	results := newResultsView(state)
	panel := newConditionsPanel(win, state, func() {
		points, cappedAt := state.Solve()
		results.Update(points, cappedAt)
	})
	title := widget.NewLabel("ballistX")
	title.TextStyle = fyne.TextStyle{Bold: true}
	subtitle := widget.NewLabel("Modular ballistics calculator for competition shooting.")
	header := container.NewVBox(title, subtitle, widget.NewSeparator())
	leftScroll := container.NewVScroll(panel.content)
	leftScroll.SetMinSize(fyne.NewSize(360, 0))
	split := container.NewHSplit(leftScroll, results.canvasObject())
	split.Offset = 0.32
	win.SetContent(container.NewBorder(header, nil, nil, nil, split))
	win.Resize(fyne.NewSize(1180, 760))
	points, cappedAt := state.Solve()
	results.Update(points, cappedAt)
	win.ShowAndRun()
}
