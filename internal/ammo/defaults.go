package ammo
func DefaultPresets() []Preset {
	return []Preset{
		{
			ID: "ppu-8x57js-196sp", Brand: "PPU", Line: "Standard", Type: "8x57 JS",
			BulletWeightGrains: 196, BC: 0.331, DragFunction: "G1",
			MuzzleVelocityFPS: 2510, Construction: SoftPoint,
			RefBarrelLengthIn: 24, VelocityPerInchFPS: 20,
			BulletDiameterIn: 0.323, BulletLengthIn: 1.25,
			MaxRangeYards: 875,
			Notes:         "Prvi Partizan 196gr soft-point. Max range limited to ~800m as a sane upper bound for this bullet/velocity combination.",
		},
		{
			ID: "fed-gmm-308-168", Brand: "Federal", Line: "Gold Medal Match", Type: ".308 Winchester",
			BulletWeightGrains: 168, BC: 0.462, DragFunction: "G1",
			MuzzleVelocityFPS: 2650, Construction: OpenTipMatch,
			RefBarrelLengthIn: 24, VelocityPerInchFPS: 18,
			BulletDiameterIn: 0.308, BulletLengthIn: 1.21,
			MaxRangeYards: 1000,
			Notes:         "Sierra 168gr MatchKing load, the classic NRA High Power/F-Class reference round.",
		},
		{
			ID: "horn-match-223-77", Brand: "Hornady", Line: "Match", Type: ".223 Remington",
			BulletWeightGrains: 77, BC: 0.235, DragFunction: "G7",
			MuzzleVelocityFPS: 2750, Construction: OpenTipMatch,
			RefBarrelLengthIn: 20, VelocityPerInchFPS: 25,
			BulletDiameterIn: 0.224, BulletLengthIn: 1.10,
			MaxRangeYards: 700,
			Notes:         "77gr ELD/OTM-class loading, common service-rifle and PRS sub-caliber stage round.",
		},
		{
			ID: "horn-eldm-65cm-140", Brand: "Hornady", Line: "Precision Hunter / ELD Match", Type: "6.5 Creedmoor",
			BulletWeightGrains: 140, BC: 0.326, DragFunction: "G7",
			MuzzleVelocityFPS: 2700, Construction: PolymerTip,
			RefBarrelLengthIn: 24, VelocityPerInchFPS: 15,
			BulletDiameterIn: 0.264, BulletLengthIn: 1.33,
			MaxRangeYards: 1300,
			Notes:         "140gr ELD-class bullet, the dominant modern PRS/NRL match cartridge.",
		},
		{
			ID: "lapua-65x55-156", Brand: "Lapua", Line: "Scenar-L", Type: "6.5x55 Swedish",
			BulletWeightGrains: 156, BC: 0.324, DragFunction: "G7",
			MuzzleVelocityFPS: 2625, Construction: OpenTipMatch,
			RefBarrelLengthIn: 24, VelocityPerInchFPS: 15,
			BulletDiameterIn: 0.264, BulletLengthIn: 1.42,
			MaxRangeYards: 1200,
			Notes:         "156gr Scenar-L, popular long-range/F-Class loading for the 6.5x55.",
		},
	}
}
