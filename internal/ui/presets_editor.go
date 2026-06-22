package ui
import (
	"strconv"
	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/dialog"
	"fyne.io/fyne/v2/widget"
	"ballistx/internal/ammo"
)
func showPresetEditor(win fyne.Window, cat *ammo.Catalog, onChange func()) {
	var d dialog.Dialog
	list := widget.NewList(
		func() int { return len(cat.Presets) },
		func() fyne.CanvasObject { return widget.NewLabel("template") },
		func(i widget.ListItemID, obj fyne.CanvasObject) {
			obj.(*widget.Label).SetText(cat.Presets[i].DisplayName())
		},
	)
	selectedIdx := -1
	list.OnSelected = func(id widget.ListItemID) { selectedIdx = id }
	refresh := func() {
		list.Refresh()
		onChange()
	}
	addBtn := widget.NewButton("Add New…", func() {
		showPresetForm(win, ammo.Preset{DragFunction: "G1", Construction: ammo.Jacketed}, func(p ammo.Preset) {
			cat.Add(p)
			cat.Save()
			refresh()
		})
	})
	editBtn := widget.NewButton("Edit Selected…", func() {
		if selectedIdx < 0 || selectedIdx >= len(cat.Presets) {
			return
		}
		showPresetForm(win, cat.Presets[selectedIdx], func(p ammo.Preset) {
			cat.Update(p)
			cat.Save()
			refresh()
		})
	})
	deleteBtn := widget.NewButton("Delete Selected", func() {
		if selectedIdx < 0 || selectedIdx >= len(cat.Presets) {
			return
		}
		target := cat.Presets[selectedIdx]
		dialog.NewConfirm("Delete preset?", "Remove "+target.DisplayName()+" from the catalog?", func(ok bool) {
			if !ok {
				return
			}
			cat.Remove(target.ID)
			cat.Save()
			selectedIdx = -1
			refresh()
		}, win).Show()
	})
	deleteBtn.Importance = widget.DangerImportance
	buttons := container.NewHBox(addBtn, editBtn, deleteBtn)
	content := container.NewBorder(nil, buttons, nil, nil, container.NewVScroll(list))
	content.Resize(fyne.NewSize(420, 420))
	d = dialog.NewCustom("Manage Ammo Presets", "Close", content, win)
	d.Resize(fyne.NewSize(460, 480))
	d.Show()
}
func showPresetForm(win fyne.Window, p ammo.Preset, onSave func(ammo.Preset)) {
	brand := widget.NewEntry()
	brand.SetText(p.Brand)
	line := widget.NewEntry()
	line.SetText(p.Line)
	cartType := widget.NewEntry()
	cartType.SetText(p.Type)
	weight := widget.NewEntry()
	weight.SetText(formatNum(p.BulletWeightGrains))
	bc := widget.NewEntry()
	bc.SetText(formatNum(p.BC))
	dragFn := widget.NewSelect([]string{"G1", "G7"}, nil)
	if p.DragFunction == "" {
		p.DragFunction = "G1"
	}
	dragFn.SetSelected(p.DragFunction)
	velocity := widget.NewEntry()
	velocity.SetText(formatNum(p.MuzzleVelocityFPS))
	constructionNames := make([]string, len(ammo.AllConstructions))
	for i, c := range ammo.AllConstructions {
		constructionNames[i] = string(c)
	}
	construction := widget.NewSelect(constructionNames, nil)
	if p.Construction == "" {
		p.Construction = ammo.Jacketed
	}
	construction.SetSelected(string(p.Construction))
	refLen := widget.NewEntry()
	refLen.SetText(formatNum(p.RefBarrelLengthIn))
	velPerInch := widget.NewEntry()
	velPerInch.SetText(formatNum(p.VelocityPerInchFPS))
	diameter := widget.NewEntry()
	diameter.SetText(formatNum(p.BulletDiameterIn))
	bulletLen := widget.NewEntry()
	bulletLen.SetText(formatNum(p.BulletLengthIn))
	maxRange := widget.NewEntry()
	maxRange.SetText(formatNum(p.MaxRangeYards))
	notes := widget.NewMultiLineEntry()
	notes.SetText(p.Notes)
	notes.Wrapping = fyne.TextWrapWord
	form := widget.NewForm(
		widget.NewFormItem("Brand", brand),
		widget.NewFormItem("Line", line),
		widget.NewFormItem("Cartridge / Type", cartType),
		widget.NewFormItem("Bullet Weight (gr)", weight),
		widget.NewFormItem("Ballistic Coefficient", bc),
		widget.NewFormItem("Drag Function", dragFn),
		widget.NewFormItem("Muzzle Velocity (fps)", velocity),
		widget.NewFormItem("Construction / Material", construction),
		widget.NewFormItem("Test Barrel Length (in)", refLen),
		widget.NewFormItem("Velocity per Inch (fps/in)", velPerInch),
		widget.NewFormItem("Bullet Diameter (in)", diameter),
		widget.NewFormItem("Bullet Length (in)", bulletLen),
		widget.NewFormItem("Max Reasonable Range (yd)", maxRange),
		widget.NewFormItem("Notes", notes),
	)
	scroll := container.NewVScroll(form)
	scroll.SetMinSize(fyne.NewSize(420, 460))
	cd := dialog.NewCustomConfirm("Ammo Preset", "Save", "Cancel", scroll, func(save bool) {
		if !save {
			return
		}
		p.Brand = brand.Text
		p.Line = line.Text
		p.Type = cartType.Text
		p.BulletWeightGrains = atof(weight.Text)
		p.BC = atof(bc.Text)
		p.DragFunction = dragFn.Selected
		p.MuzzleVelocityFPS = atof(velocity.Text)
		p.Construction = ammo.Construction(construction.Selected)
		p.RefBarrelLengthIn = atof(refLen.Text)
		p.VelocityPerInchFPS = atof(velPerInch.Text)
		p.BulletDiameterIn = atof(diameter.Text)
		p.BulletLengthIn = atof(bulletLen.Text)
		p.MaxRangeYards = atof(maxRange.Text)
		p.Notes = notes.Text
		onSave(p)
	}, win)
	cd.Resize(fyne.NewSize(460, 540))
	cd.Show()
}
func atof(s string) float64 {
	v, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0
	}
	return v
}
