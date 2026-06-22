package ammo
import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"sync/atomic"
	"time"
)
type Construction string
const (
	Jacketed     Construction = "Jacketed Lead-Core"
	Monolithic   Construction = "Monolithic Copper"
	PolymerTip   Construction = "Polymer Tip"
	OpenTipMatch Construction = "Open-Tip Match (OTM)"
	SoftPoint    Construction = "Soft Point"
	FMJ          Construction = "Full Metal Jacket"
)
var AllConstructions = []Construction{Jacketed, Monolithic, PolymerTip, OpenTipMatch, SoftPoint, FMJ}
type Preset struct {
	ID    string `json:"id"`
	Brand string `json:"brand"`
	Line  string `json:"line"`
	Type  string `json:"type"`
	BulletWeightGrains float64      `json:"bulletWeightGrains"`
	BC                 float64      `json:"bc"`
	DragFunction       string       `json:"dragFunction"`
	MuzzleVelocityFPS  float64      `json:"muzzleVelocityFps"`
	Construction       Construction `json:"construction"`
	BCAdjustmentPct    float64      `json:"bcAdjustmentPct"`
	RefBarrelLengthIn  float64 `json:"refBarrelLengthIn"`
	VelocityPerInchFPS float64 `json:"velocityPerInchFps"`
	BulletDiameterIn float64 `json:"bulletDiameterIn"`
	BulletLengthIn   float64 `json:"bulletLengthIn"`
	MaxRangeYards float64 `json:"maxRangeYards"`
	Notes         string  `json:"notes"`
}
func (p Preset) DisplayName() string {
	return fmt.Sprintf("%s %s %s (%.0fgr)", p.Brand, p.Line, p.Type, p.BulletWeightGrains)
}
func (p Preset) EffectiveBC() float64 {
	return p.BC * (1.0 + p.BCAdjustmentPct/100.0)
}
type Catalog struct {
	Path    string   `json:"-"`
	Presets []Preset `json:"presets"`
}
var idCounter atomic.Int64
func NewID() string {
	idCounter.Add(1)
	return fmt.Sprintf("p-%d-%d", time.Now().UnixNano(), idCounter.Load())
}
func DefaultPath() string {
	dir, err := os.UserConfigDir()
	if err != nil || dir == "" {
		dir = "."
	}
	return filepath.Join(dir, "ballistx", "catalog.json")
}
func LoadOrCreateDefault(path string) (*Catalog, error) {
	if _, err := os.Stat(path); os.IsNotExist(err) {
		cat := &Catalog{Path: path, Presets: DefaultPresets()}
		if err := cat.Save(); err != nil {
			return nil, err
		}
		return cat, nil
	}
	return Load(path)
}
func Load(path string) (*Catalog, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var cat Catalog
	if err := json.Unmarshal(data, &cat); err != nil {
		return nil, fmt.Errorf("parse catalog %s: %w", path, err)
	}
	cat.Path = path
	return &cat, nil
}
func (c *Catalog) Save() error {
	if err := os.MkdirAll(filepath.Dir(c.Path), 0o755); err != nil {
		return err
	}
	data, err := json.MarshalIndent(c, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(c.Path, data, 0o644)
}
func (c *Catalog) Add(p Preset) Preset {
	if p.ID == "" {
		p.ID = NewID()
	}
	c.Presets = append(c.Presets, p)
	c.sort()
	return p
}
func (c *Catalog) Update(p Preset) bool {
	for i := range c.Presets {
		if c.Presets[i].ID == p.ID {
			c.Presets[i] = p
			c.sort()
			return true
		}
	}
	return false
}
func (c *Catalog) Remove(id string) bool {
	for i := range c.Presets {
		if c.Presets[i].ID == id {
			c.Presets = append(c.Presets[:i], c.Presets[i+1:]...)
			return true
		}
	}
	return false
}
func (c *Catalog) sort() {
	sort.Slice(c.Presets, func(i, j int) bool {
		a, b := c.Presets[i], c.Presets[j]
		if a.Brand != b.Brand {
			return a.Brand < b.Brand
		}
		if a.Line != b.Line {
			return a.Line < b.Line
		}
		return a.Type < b.Type
	})
}
