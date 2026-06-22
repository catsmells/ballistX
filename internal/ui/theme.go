package ui
import (
	"image/color"
	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/theme"
)
type ballistXTheme struct{}
var _ fyne.Theme = (*ballistXTheme)(nil)
var (
	colBackground  = color.NRGBA{R: 0xF8, G: 0xF6, B: 0xF1, A: 0xFF}
	colInputBg     = color.NRGBA{R: 0xFF, G: 0xFF, B: 0xFE, A: 0xFF}
	colHeaderBg    = color.NRGBA{R: 0xEF, G: 0xEA, B: 0xE1, A: 0xFF}
	colSeparator   = color.NRGBA{R: 0xDE, G: 0xD6, B: 0xC8, A: 0xFF}
	colInputBorder = color.NRGBA{R: 0xD3, G: 0xCA, B: 0xB9, A: 0xFF}
	colForeground  = color.NRGBA{R: 0x32, G: 0x3A, B: 0x45, A: 0xFF}
	colDisabled    = color.NRGBA{R: 0xA3, G: 0x9E, B: 0x93, A: 0xFF}
	colPlaceholder = color.NRGBA{R: 0xAE, G: 0xA6, B: 0x99, A: 0xFF}
	colPrimary      = color.NRGBA{R: 0x4F, G: 0x6E, B: 0x8C, A: 0xFF}
	colPrimaryHover = color.NRGBA{R: 0x5E, G: 0x7E, B: 0x9C, A: 0xFF}
	colSage         = color.NRGBA{R: 0x6F, G: 0x91, B: 0x76, A: 0xFF}
	colPlum         = color.NRGBA{R: 0x9C, G: 0x62, B: 0x76, A: 0xFF}
	colHyperlink    = color.NRGBA{R: 0x4F, G: 0x6E, B: 0x8C, A: 0xFF}
	colSelection = color.NRGBA{R: 0xCF, G: 0xDC, B: 0xE4, A: 0xB0}
	colHover     = color.NRGBA{R: 0xEA, G: 0xE3, B: 0xD8, A: 0xFF}
	colPressed   = color.NRGBA{R: 0xE0, G: 0xD8, B: 0xCA, A: 0xFF}
)
func (ballistXTheme) Color(name fyne.ThemeColorName, _ fyne.ThemeVariant) color.Color {
	switch name {
	case theme.ColorNameBackground:
		return colBackground
	case theme.ColorNameMenuBackground, theme.ColorNameOverlayBackground:
		return colBackground
	case theme.ColorNameHeaderBackground:
		return colHeaderBg
	case theme.ColorNameInputBackground:
		return colInputBg
	case theme.ColorNameInputBorder:
		return colInputBorder
	case theme.ColorNameSeparator:
		return colSeparator
	case theme.ColorNameForeground:
		return colForeground
	case theme.ColorNameDisabled, theme.ColorNameDisabledButton:
		return colDisabled
	case theme.ColorNamePlaceHolder:
		return colPlaceholder
	case theme.ColorNamePrimary:
		return colPrimary
	case theme.ColorNameHyperlink:
		return colHyperlink
	case theme.ColorNameButton:
		return colInputBg
	case theme.ColorNameHover:
		return colHover
	case theme.ColorNamePressed:
		return colPressed
	case theme.ColorNameFocus:
		return colPrimaryHover
	case theme.ColorNameSelection:
		return colSelection
	case theme.ColorNameSuccess:
		return colSage
	case theme.ColorNameError, theme.ColorNameWarning:
		return colPlum
	case theme.ColorNameForegroundOnPrimary, theme.ColorNameForegroundOnError,
		theme.ColorNameForegroundOnSuccess, theme.ColorNameForegroundOnWarning:
		return colBackground
	case theme.ColorNameScrollBar:
		return colInputBorder
	case theme.ColorNameScrollBarBackground:
		return colBackground
	case theme.ColorNameShadow:
		return color.NRGBA{A: 0x22}
	}
	return theme.DefaultTheme().Color(name, theme.VariantLight)
}
func (ballistXTheme) Font(style fyne.TextStyle) fyne.Resource {
	return theme.DefaultTheme().Font(style)
}
func (ballistXTheme) Icon(name fyne.ThemeIconName) fyne.Resource {
	return theme.DefaultTheme().Icon(name)
}
func (ballistXTheme) Size(name fyne.ThemeSizeName) float32 {
	switch name {
	case theme.SizeNamePadding:
		return 6
	case theme.SizeNameInputRadius, theme.SizeNameSelectionRadius:
		return 8
	}
	return theme.DefaultTheme().Size(name)
}
